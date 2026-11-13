/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Raven-MODFLOW 6 groundwater coupling: GeoTIFF reader (see GWGeoTIFF.h)

  References: TIFF 6.0 specification (Adobe, 1992) with the TIFF Technical Note 2 DEFLATE and the floating-point
  predictor of Adobe's Photoshop TIFF Technical Note 3; BigTIFF (libtiff design note); OGC GeoTIFF 1.1 standard;
  RFC 1950/1951 (zlib and DEFLATE). Every offset and count read from the file is checked against the file size and
  every decoder against its input and output buffers, so a truncated or corrupt file gives an error, never a read
  outside the buffers.
----------------------------------------------------------------*/
#include "GWGeoTIFF.h"
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <set>
#include <stdint.h>
#include "GWEPSGTable.h"

gw_geotiff::gw_geotiff()
{
  ncols=nrows=0; xll=yll=0.0; dx=dy=1.0; hasNodata=false; nodata=0.0;
  modelType=0; rasterType=1; epsg=0; geographic=false;
}

//=================================================================
// decoders
//=================================================================
namespace
{
  //---------------------------------------------------------------
  // DEFLATE (RFC 1951): canonical Huffman codes decoded bit by bit
  //---------------------------------------------------------------
  struct InBits
  {
    const unsigned char *p; size_t n,pos; uint32_t buf; int cnt; bool eof;
    InBits(const unsigned char *d,size_t len):p(d),n(len),pos(0),buf(0),cnt(0),eof(false){}
    int Bits(const int need)  //need <= 16
    {
      uint32_t v=buf;
      while (cnt<need){
        if (pos>=n){eof=true; return 0;}
        v|=(uint32_t)p[pos++]<<cnt; cnt+=8;
      }
      buf=v>>need; cnt-=need;
      return (int)(v&((1u<<need)-1u));
    }
  };
  struct Huffman
  {
    short count[16];   ///< number of codes of each length
    short symbol[320]; ///< symbols ordered by code
  };
  //builds the decoding table from code lengths; returns 0 for a complete code, >0 incomplete, <0 over-subscribed
  int Build(Huffman &h,const short *length,const int n)
  {
    for (int l=0;l<16;l++){h.count[l]=0;}
    for (int s=0;s<n;s++){h.count[length[s]]++;}
    if (h.count[0]==n){return 0;}
    int left=1;
    for (int l=1;l<16;l++){left<<=1; left-=h.count[l]; if (left<0){return left;}}
    short offs[16]; offs[1]=0;
    for (int l=1;l<15;l++){offs[l+1]=(short)(offs[l]+h.count[l]);}
    for (int s=0;s<n;s++){if (length[s]!=0){h.symbol[offs[length[s]]++]=(short)s;}}
    return left;
  }
  int Decode(InBits &b,const Huffman &h)  //-1 end of input, -2 invalid code
  {
    int code=0,first=0,index=0;
    for (int l=1;l<16;l++){
      code|=b.Bits(1); if (b.eof){return -1;}
      int c=h.count[l];
      if (code-c<first){return h.symbol[index+(code-first)];}
      index+=c; first+=c; first<<=1; code<<=1;
    }
    return -2;
  }
  const short LBASE[29]={3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258};
  const short LEXT [29]={0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0};
  const short DBASE[30]={1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
  const short DEXT [30]={0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13};

  //one compressed block (fixed or dynamic codes); returns 0 at the end-of-block, 1 when the output is full, <0 on error
  int Codes(InBits &b,vector<unsigned char> &out,const size_t limit,const Huffman &lc,const Huffman &dc,string &err)
  {
    while (true){
      int s=Decode(b,lc);
      if (s<0){err=(s==-1)?"the DEFLATE data end in the middle of a block":"invalid DEFLATE code"; return -1;}
      if (s<256){
        if (out.size()>=limit){return 1;}
        out.push_back((unsigned char)s);
      }
      else if (s==256){return 0;}
      else {
        s-=257;
        if (s>=29){err="invalid DEFLATE length code"; return -1;}
        int len=LBASE[s]+b.Bits(LEXT[s]);
        int ds=Decode(b,dc);
        if (ds<0){err=(ds==-1)?"the DEFLATE data end in the middle of a block":"invalid DEFLATE distance code"; return -1;}
        if (ds>=30){err="invalid DEFLATE distance code"; return -1;}
        size_t dist=(size_t)(DBASE[ds]+b.Bits(DEXT[ds]));
        if (b.eof){err="the DEFLATE data end in the middle of a block"; return -1;}
        if (dist>out.size()){err="DEFLATE distance too far back"; return -1;}
        for (int k=0;k<len;k++){
          if (out.size()>=limit){return 1;}
          out.push_back(out[out.size()-dist]);
        }
      }
    }
  }
  //raw DEFLATE stream; returns 0 after the last block, 1 when the output is full, <0 on error (pos: bytes used)
  int Inflate(InBits &b,vector<unsigned char> &out,const size_t limit,string &err)
  {
    static const short ORDER[19]={16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15};
    int last=0;
    do {
      last=b.Bits(1);
      int type=b.Bits(2);
      if (b.eof){err="the DEFLATE data end before the last block"; return -1;}
      if (type==0){            //stored block
        b.buf=0; b.cnt=0;      //(to the next byte)
        if (b.pos+4>b.n){err="the DEFLATE data end in a stored block header"; return -1;}
        unsigned len =b.p[b.pos]|(b.p[b.pos+1]<<8);
        unsigned nlen=b.p[b.pos+2]|(b.p[b.pos+3]<<8);
        b.pos+=4;
        if (len!=(~nlen&0xffffu)){err="corrupt DEFLATE stored block (length check)"; return -1;}
        if (b.pos+len>b.n){err="the DEFLATE data end in a stored block"; return -1;}
        for (unsigned k=0;k<len;k++){
          if (out.size()>=limit){return 1;}
          out.push_back(b.p[b.pos++]);
        }
      }
      else if (type==1){       //fixed codes
        Huffman lc,dc; short L[288];
        for (int s=0;s<144;s++){L[s]=8;} for (int s=144;s<256;s++){L[s]=9;} for (int s=256;s<280;s++){L[s]=7;} for (int s=280;s<288;s++){L[s]=8;}
        Build(lc,L,288);
        for (int s=0;s<30;s++){L[s]=5;}
        Build(dc,L,30);
        int r=Codes(b,out,limit,lc,dc,err); if (r!=0){return r;}
      }
      else if (type==2){       //dynamic codes
        int nlen=b.Bits(5)+257,ndist=b.Bits(5)+1,ncode=b.Bits(4)+4;
        if (b.eof){err="the DEFLATE data end in a block header"; return -1;}
        if ((nlen>286) || (ndist>30)){err="corrupt DEFLATE block header (too many codes)"; return -1;}
        short L[320]; Huffman lc,dc;
        for (int i=0;i<19;i++){L[ORDER[i]]=(short)((i<ncode)?b.Bits(3):0);}
        if (b.eof){err="the DEFLATE data end in a block header"; return -1;}
        if (Build(lc,L,19)!=0){err="corrupt DEFLATE block header (code-length code)"; return -1;}
        int i=0;
        while (i<nlen+ndist){
          int s=Decode(b,lc);
          if (s<0){err="corrupt DEFLATE block header (code lengths)"; return -1;}
          if (s<16){L[i++]=(short)s; continue;}
          int rep=0; short val=0;
          if (s==16){if (i==0){err="corrupt DEFLATE block header (repeat with no previous length)"; return -1;} val=L[i-1]; rep=3+b.Bits(2);}
          else if (s==17){rep=3+b.Bits(3);}
          else {rep=11+b.Bits(7);}
          if (b.eof){err="the DEFLATE data end in a block header"; return -1;}
          if (i+rep>nlen+ndist){err="corrupt DEFLATE block header (too many lengths)"; return -1;}
          while (rep--){L[i++]=val;}
        }
        if (L[256]==0){err="corrupt DEFLATE block header (no end-of-block code)"; return -1;}
        int e=Build(lc,L,nlen);
        if ((e<0) || ((e>0) && (nlen-lc.count[0]!=1))){err="corrupt DEFLATE block header (literal/length code)"; return -1;}
        e=Build(dc,L+nlen,ndist);
        if ((e<0) || ((e>0) && (ndist-dc.count[0]!=1))){err="corrupt DEFLATE block header (distance code)"; return -1;}
        int r=Codes(b,out,limit,lc,dc,err); if (r!=0){return r;}
      }
      else {err="invalid DEFLATE block type"; return -1;}
    } while (!last);
    return 0;
  }
}

bool GWGeoTIFF::InflateZlib(const unsigned char *in,const size_t n,vector<unsigned char> &out,const size_t limit,bool &full,string &err)
{
  out.clear(); full=false;
  if (n<2){err="DEFLATE data too short"; return false;}
  unsigned cmf=in[0],flg=in[1];
  if (((cmf&15)!=8) || ((cmf>>4)>7) || (((cmf<<8)|flg)%31!=0)){err="not a zlib (DEFLATE) stream (bad header)"; return false;}
  if (flg&0x20){err="zlib stream with a preset dictionary"; return false;}
  out.reserve(limit);
  InBits b(in+2,n-2);
  int r=Inflate(b,out,limit,err);
  if (r<0){return false;}
  if (r==1){full=true; return true;}
  //Adler-32 of the output (big-endian, after the last block)
  if (b.pos+4>b.n){err="the DEFLATE data end before the checksum"; return false;}
  uint32_t want=((uint32_t)b.p[b.pos]<<24)|((uint32_t)b.p[b.pos+1]<<16)|((uint32_t)b.p[b.pos+2]<<8)|(uint32_t)b.p[b.pos+3];
  uint32_t s1=1,s2=0;
  for (size_t i=0;i<out.size();i++){s1=(s1+out[i])%65521u; s2=(s2+s1)%65521u;}
  if (((s2<<16)|s1)!=want){err="DEFLATE checksum mismatch (corrupt data)"; return false;}
  return true;
}

//TIFF LZW: codes of 9 to 12 bits, most significant bit first, code 256 clears the table, 257 ends the data, and the code
//width grows one code early (at 511, 1023 and 2047 entries)
bool GWGeoTIFF::LZWDecode(const unsigned char *in,const size_t n,vector<unsigned char> &out,const size_t limit,bool &full,string &err)
{
  out.clear(); full=false; out.reserve(limit);
  if ((n>=2) && (in[0]==0) && (in[1]&1)){err="old-style LZW data (written by libtiff before version 5) are not supported"; return false;}
  vector<int> prefix(4096,-1),len(4096,1); vector<unsigned char> suffix(4096),first(4096);
  for (int i=0;i<256;i++){suffix[i]=first[i]=(unsigned char)i;}
  int next=258,width=9,old=-1;
  uint32_t buf=0; int cnt=0; size_t pos=0;
  unsigned char tmp[4096];
  while (true){
    while (cnt<width){
      if (pos>=n){return true;}        //(no end code: the data end here; a short output is found by the caller)
      buf=(buf<<8)|in[pos++]; cnt+=8;
    }
    int code=(int)((buf>>(cnt-width))&((1u<<width)-1u)); cnt-=width;
    if (code==257){return true;}
    if (code==256){next=258; width=9; old=-1; continue;}
    int c=code;
    if (old==-1){
      if (code>255){err="corrupt LZW data (first code after a clear is not a byte)"; return false;}
    }
    else if (code<next){
      if (next<4096){prefix[next]=old; suffix[next]=first[code]; first[next]=first[old]; len[next]=len[old]+1; next++;}
    }
    else if (code==next){
      if (next<4096){prefix[next]=old; suffix[next]=first[old]; first[next]=first[old]; len[next]=len[old]+1; next++;}
      else {err="corrupt LZW data (code beyond a full table)"; return false;}
    }
    else {err="corrupt LZW data (code not yet defined)"; return false;}
    if ((next>=(1<<width)-1) && (width<12)){width++;}
    int L=len[c]; if ((L<1) || (L>4096)){err="corrupt LZW data"; return false;}
    int k=L,q=c;
    while ((k>0) && (q>=0)){tmp[--k]=suffix[q]; q=prefix[q];}
    if (k!=0){err="corrupt LZW data"; return false;}
    for (int i=0;i<L;i++){
      if (out.size()>=limit){full=true; return true;}
      out.push_back(tmp[i]);
    }
    old=code;
  }
}

bool GWGeoTIFF::PackBits(const unsigned char *in,const size_t n,vector<unsigned char> &out,const size_t limit,bool &full,string &err)
{
  out.clear(); full=false; out.reserve(limit);
  size_t i=0;
  while (i<n){
    if (out.size()>=limit){full=true; return true;}
    int c=(signed char)in[i++];
    if (c>=0){
      size_t cnt=(size_t)c+1;
      if (i+cnt>n){err="PackBits data end in the middle of a literal run"; return false;}
      for (size_t k=0;k<cnt;k++){if (out.size()>=limit){full=true; return true;} out.push_back(in[i+k]);}
      i+=cnt;
    }
    else if (c!=-128){
      if (i>=n){err="PackBits data end in the middle of a repeat run"; return false;}
      unsigned char v=in[i++];
      for (int k=0;k<1-c;k++){if (out.size()>=limit){full=true; return true;} out.push_back(v);}
    }
  }
  return true;
}

//=================================================================
// TIFF structure
//=================================================================
namespace
{
  struct TEntry { uint16_t type; uint64_t count; uint64_t off; }; ///< off: file offset of the values
  int TypeSize(const int t)
  {
    switch (t){
      case 1: case 2: case 6: case 7: return 1;
      case 3: case 8: return 2;
      case 4: case 9: case 11: case 13: return 4;
      case 5: case 10: case 12: case 16: case 17: case 18: return 8;
      default: return 0;
    }
  }
  string CompressionName(const uint64_t c)
  {
    switch (c){
      case 1: return "none"; case 2: return "CCITT RLE"; case 3: return "CCITT Group 3 fax"; case 4: return "CCITT Group 4 fax";
      case 5: return "LZW"; case 6: return "old-style JPEG"; case 7: return "JPEG"; case 8: return "DEFLATE"; case 9: return "JBIG";
      case 32773: return "PackBits"; case 32946: return "DEFLATE"; case 34661: return "JBIG"; case 34676: return "SGILog";
      case 34677: return "SGILog24"; case 34712: return "JPEG 2000"; case 34887: return "LERC"; case 34925: return "LZMA";
      case 50000: return "ZSTD"; case 50001: return "WebP"; case 50002: return "JPEG XL";
      default: return "unknown";
    }
  }
  class TIFF
  {
  public:
    const unsigned char *d; uint64_t n; bool be,big; map<int,TEntry> tags; uint64_t first,next; ///< (first and next image directory)
    TIFF(const unsigned char *data,uint64_t len):d(data),n(len),be(false),big(false),first(0),next(0){}
    bool In(const uint64_t off,const uint64_t len) const {return (off<=n) && (len<=n-off);}
    uint64_t U(const uint64_t off,const int bytes) const //(caller checks In)
    {
      uint64_t v=0;
      if (be){for (int i=0;i<bytes;i++){v=(v<<8)|d[off+i];}}
      else   {for (int i=bytes-1;i>=0;i--){v=(v<<8)|d[off+i];}}
      return v;
    }
    bool Header(string &err)
    {
      if (n<8){err="the file is too short for a TIFF header"; return false;}
      if      ((d[0]=='I') && (d[1]=='I')){be=false;}
      else if ((d[0]=='M') && (d[1]=='M')){be=true;}
      else {err="not a TIFF file (no II or MM byte-order mark)"; return false;}
      uint64_t magic=U(2,2);
      if (magic==42){big=false;}
      else if (magic==43){
        big=true;
        if (n<16){err="the file is too short for a BigTIFF header"; return false;}
        if ((U(4,2)!=8) || (U(6,2)!=0)){err="unsupported BigTIFF header (offset size not 8)"; return false;}
      }
      else {err="not a TIFF file (version number not 42 or 43)"; return false;}
      return true;
    }
    bool FirstIFD(string &err)
    {
      first=big?U(8,8):U(4,4);
      return ReadIFD(first,err);
    }
    //the entries of the image directory at `ifd` (into tags, cleared) and the offset of the next one (next; 0: none)
    bool ReadIFD(const uint64_t ifd,string &err)
    {
      tags.clear(); next=0;
      int csz=big?8:2,esz=big?20:12,vsz=big?8:4;
      if ((ifd<(uint64_t)(big?16:8)) || !In(ifd,csz)){err="an image directory lies outside the file"; return false;}
      uint64_t ne=U(ifd,csz);
      if ((ne==0) || (ne>4096) || !In(ifd+csz,ne*esz)){err="corrupt image directory (entry count "+to_string((unsigned long long)ne)+")"; return false;}
      if (In(ifd+csz+ne*esz,vsz)){next=U(ifd+csz+ne*esz,vsz);}
      for (uint64_t i=0;i<ne;i++){
        uint64_t e=ifd+csz+i*esz;
        int tag=(int)U(e,2),type=(int)U(e+2,2);
        uint64_t count=U(e+4,big?8:4);
        int ts=TypeSize(type);
        if (ts==0){continue;}                            //(unknown type: the entry is skipped)
        if (count>n/(uint64_t)ts){
          err="corrupt image directory (tag "+to_string(tag)+" has "+to_string((unsigned long long)count)+" values, more than the file holds)"; return false;
        }
        uint64_t size=count*ts,off;
        if (size<=(uint64_t)vsz){off=e+4+(big?8:4);}
        else {
          off=U(e+4+(big?8:4),vsz);
          if (!In(off,size)){err="corrupt image directory (the values of tag "+to_string(tag)+" lie outside the file)"; return false;}
        }
        if (tags.find(tag)==tags.end()){TEntry t; t.type=(uint16_t)type; t.count=count; t.off=off; tags[tag]=t;}
      }
      return true;
    }
    const TEntry *Tag(const int tag) const {map<int,TEntry>::const_iterator it=tags.find(tag); return (it==tags.end())?NULL:&it->second;}
    //value i of an entry as an unsigned integer (integer types only; negative values refused)
    bool UInt(const TEntry &t,const uint64_t i,uint64_t &v) const
    {
      if (i>=t.count){return false;}
      int ts=TypeSize(t.type); uint64_t o=t.off+i*ts;
      switch (t.type){
        case 1: case 3: case 4: case 13: case 16: case 18: v=U(o,ts); return true;
        case 6: case 8: case 9: case 17: {
          uint64_t u=U(o,ts); int64_t s=(ts==8)?(int64_t)u:((int64_t)(u<<(64-8*ts))>>(64-8*ts));
          if (s<0){return false;} v=(uint64_t)s; return true;
        }
        default: return false;
      }
    }
    bool Double(const TEntry &t,const uint64_t i,double &v) const
    {
      if (i>=t.count){return false;}
      int ts=TypeSize(t.type); uint64_t o=t.off+i*ts;
      switch (t.type){
        case 1: case 3: case 4: case 13: case 16: case 18: v=(double)U(o,ts); return true;
        case 6: case 8: case 9: case 17: {
          uint64_t u=U(o,ts); int64_t s=(ts==8)?(int64_t)u:((int64_t)(u<<(64-8*ts))>>(64-8*ts)); v=(double)s; return true;
        }
        case 5: {uint64_t a=U(o,4),b=U(o+4,4); if (b==0){return false;} v=(double)a/(double)b; return true;}
        case 10:{int32_t a=(int32_t)(uint32_t)U(o,4),b=(int32_t)(uint32_t)U(o+4,4); if (b==0){return false;} v=(double)a/(double)b; return true;}
        case 11:{uint32_t u=(uint32_t)U(o,4); float f; memcpy(&f,&u,4); v=f; return true;}
        case 12:{uint64_t u=U(o,8); memcpy(&v,&u,8); return true;}
        default: return false;
      }
    }
    //a single integer tag with a default
    bool One(const int tag,const uint64_t def,uint64_t &v,string &err) const
    {
      const TEntry *t=Tag(tag); if (t==NULL){v=def; return true;}
      if ((t->count<1) || !UInt(*t,0,v)){err="invalid value of TIFF tag "+to_string(tag); return false;}
      return true;
    }
    bool List(const int tag,vector<uint64_t> &v,const uint64_t need,string &err) const
    {
      const TEntry *t=Tag(tag);
      if (t==NULL){err="TIFF tag "+to_string(tag)+" is missing"; return false;}
      if (t->count<need){err="TIFF tag "+to_string(tag)+" has "+to_string((unsigned long long)t->count)+" values, "+to_string((unsigned long long)need)+" expected (corrupt file)"; return false;}
      v.resize((size_t)need);
      for (uint64_t i=0;i<need;i++){if (!UInt(*t,i,v[(size_t)i])){err="invalid value of TIFF tag "+to_string(tag); return false;}}
      return true;
    }
    bool Doubles(const int tag,vector<double> &v) const
    {
      const TEntry *t=Tag(tag); v.clear(); if (t==NULL){return false;}
      if (t->count>4096){return false;}
      v.resize((size_t)t->count);
      for (uint64_t i=0;i<t->count;i++){if (!Double(*t,i,v[(size_t)i])){v.clear(); return false;}}
      return true;
    }
    string Ascii(const int tag,const uint64_t maxlen=256) const
    {
      const TEntry *t=Tag(tag); if ((t==NULL) || (t->type!=2) || (t->count>maxlen)){return "";}
      string s((const char*)d+t->off,(size_t)t->count);
      size_t z=s.find('\0'); if (z!=string::npos){s.erase(z);}
      return s;
    }
  };

  //value of one sample (bytes in big-endian order when be)
  inline double SampleValue(const unsigned char *p,const int bytes,const int fmt,const bool be)
  {
    uint64_t u=0;
    if (be){for (int i=0;i<bytes;i++){u=(u<<8)|p[i];}}
    else   {for (int i=bytes-1;i>=0;i--){u=(u<<8)|p[i];}}
    if (fmt==3){
      if (bytes==4){uint32_t w=(uint32_t)u; float f; memcpy(&f,&w,4); return (double)f;}
      double x; memcpy(&x,&u,8); return x;
    }
    if (fmt==2){
      int64_t s=(bytes==8)?(int64_t)u:((int64_t)(u<<(64-8*bytes))>>(64-8*bytes));
      return (double)s;
    }
    return (double)u;
  }
  //undoes the horizontal predictor (2): integer differences of the samples of each row, in the file's byte order
  void UndoHorizontal(unsigned char *buf,const size_t rows,const size_t w,const int bytes,const bool be)
  {
    for (size_t r=0;r<rows;r++){
      unsigned char *row=buf+r*w*bytes;
      if (bytes==1){for (size_t i=1;i<w;i++){row[i]=(unsigned char)(row[i]+row[i-1]);} continue;}
      uint64_t prev=0;
      for (size_t i=0;i<w;i++){
        unsigned char *p=row+i*bytes; uint64_t u=0;
        if (be){for (int k=0;k<bytes;k++){u=(u<<8)|p[k];}} else {for (int k=bytes-1;k>=0;k--){u=(u<<8)|p[k];}}
        if (i>0){u+=prev;}
        if (bytes<8){u&=((uint64_t)1<<(8*bytes))-1;}
        prev=u;
        if (be){for (int k=bytes-1;k>=0;k--){p[k]=(unsigned char)(u&0xff); u>>=8;}} else {for (int k=0;k<bytes;k++){p[k]=(unsigned char)(u&0xff); u>>=8;}}
      }
    }
  }
  //undoes the floating-point predictor (3): byte differences along the row, then the byte planes (most significant first)
  //put back together; the samples come out big-endian whatever the file's byte order
  void UndoFloat(unsigned char *buf,const size_t rows,const size_t w,const int bytes,vector<unsigned char> &tmp)
  {
    size_t rb=w*bytes; tmp.resize(rb);
    for (size_t r=0;r<rows;r++){
      unsigned char *row=buf+r*rb;
      for (size_t i=1;i<rb;i++){row[i]=(unsigned char)(row[i]+row[i-1]);}
      memcpy(&tmp[0],row,rb);
      for (size_t i=0;i<w;i++){for (int k=0;k<bytes;k++){row[i*bytes+k]=tmp[k*w+i];}}
    }
  }
  string Fmt(const double v){ostringstream s; s.precision(12); s<<v; return s.str();}
}

namespace
{
  //strips or tiles of one image directory
  struct Blocks
  {
    bool tiled; uint64_t bw,bh,nacross,ndown,n; vector<uint64_t> offs,cnts;
    uint64_t Rows(const uint64_t b,const uint64_t h) const {return tiled?bh:std::min(bh,h-(b/nacross)*bh);} //(a tile is always whole; the last strip may be shorter)
  };
  bool Layout(const TIFF &T,const uint64_t w,const uint64_t h,Blocks &B,string &e)
  {
    B.tiled=(T.Tag(322)!=NULL) || (T.Tag(324)!=NULL);
    if (B.tiled){
      if (!T.One(322,0,B.bw,e) || !T.One(323,0,B.bh,e)){return false;}
      if ((B.bw==0) || (B.bh==0)){e="tile width or length is zero or missing (corrupt file)"; return false;}
      if ((double)B.bw*(double)B.bh>GWGeoTIFF::MAX_PIXELS){e="tiles of "+to_string((unsigned long long)B.bw)+" x "+to_string((unsigned long long)B.bh)+" pixels are larger than Raven reads (corrupt file)"; return false;}
      B.nacross=(w+B.bw-1)/B.bw; B.ndown=(h+B.bh-1)/B.bh;
    }
    else {
      uint64_t rps; if (!T.One(278,h,rps,e)){return false;}
      if (rps==0){e="RowsPerStrip is zero (corrupt file)"; return false;}
      if (rps>h){rps=h;}
      B.bw=w; B.bh=rps; B.nacross=1; B.ndown=(h+rps-1)/rps;
    }
    B.n=B.nacross*B.ndown;
    if ((double)B.n*(double)B.bw*(double)B.bh>4.0*GWGeoTIFF::MAX_PIXELS){  //(bounds the decoding work of a crafted file)
      e="its "+to_string((unsigned long long)B.n)+(B.tiled?" tiles of ":" strips of ")+to_string((unsigned long long)B.bw)+" x "+to_string((unsigned long long)B.bh)+
        " pixels cover more than "+to_string((unsigned long long)(4.0*GWGeoTIFF::MAX_PIXELS/1e6))+" million pixels (corrupt file)"; return false;
    }
    return T.List(B.tiled?324:273,B.offs,B.n,e) && T.List(B.tiled?325:279,B.cnts,B.n,e);
  }
  //the sizes of the blocks against the bytes they must give (the largest compression ratios of the codecs: PackBits 64,
  //DEFLATE 1032, LZW about 1300), so a corrupt size stops before anything is allocated
  bool Plausible(const Blocks &B,const uint64_t w,const uint64_t h,const uint64_t comp,const double rowBytes,string &e)
  {
    double have=0,want=0;
    const double ratio=(comp==1)?1.0:(comp==32773)?64.0:(comp==5)?2048.0:1100.0;
    for (uint64_t b=0;b<B.n;b++){
      if (B.cnts[(size_t)b]==0){continue;}
      have+=(double)B.cnts[(size_t)b]; want+=(double)B.Rows(b,h)*rowBytes;
    }
    if (have*ratio<want){e="the file holds "+to_string((unsigned long long)have)+" bytes of image data, too few for "+to_string((unsigned long long)w)+" x "+
                         to_string((unsigned long long)h)+" pixels (truncated or corrupt file)"; return false;}
    return true;
  }
  //block b decompressed into buf (need bytes) and its predictor undone (rows of bw samples of nb bytes)
  bool DecodeBlock(const TIFF &T,const Blocks &B,const uint64_t b,const uint64_t comp,const uint64_t pred,const size_t need,
                   const uint64_t rows,const int nb,vector<unsigned char> &buf,vector<unsigned char> &tmp,string &e)
  {
    const string what=string(B.tiled?"tile ":"strip ")+to_string((unsigned long long)b);
    if (!T.In(B.offs[(size_t)b],B.cnts[(size_t)b])){e="the data of "+what+" lie outside the file (truncated or corrupt file)"; return false;}
    const unsigned char *src=T.d+B.offs[(size_t)b]; size_t srcn=(size_t)B.cnts[(size_t)b];
    bool full=false,ok=true; string de;
    if (comp==1){
      if (srcn<need){e=what+" is shorter than its pixels (truncated or corrupt file)"; return false;}
      buf.assign(src,src+need);
    }
    else if (comp==5)    {ok=GWGeoTIFF::LZWDecode  (src,srcn,buf,need,full,de);}
    else if (comp==32773){ok=GWGeoTIFF::PackBits   (src,srcn,buf,need,full,de);}
    else                 {ok=GWGeoTIFF::InflateZlib(src,srcn,buf,need,full,de);}
    if (!ok){e="decompression error in "+what+" ("+de+")"; return false;}
    if (buf.size()<need){e=what+" decompresses to "+to_string((unsigned long long)buf.size())+" bytes, "+to_string((unsigned long long)need)+" expected (truncated or corrupt file)"; return false;}
    if      (pred==2){UndoHorizontal(&buf[0],(size_t)rows,(size_t)B.bw,nb,T.be);}
    else if (pred==3){UndoFloat     (&buf[0],(size_t)rows,(size_t)B.bw,nb,tmp);}
    return true;
  }
  string UnitName(const uint64_t u)
  {
    switch (u){
      case 9001: return "metre"; case 9002: return "foot"; case 9003: return "US survey foot"; case 9005: return "Clarke's foot";
      case 9014: return "fathom"; case 9030: return "nautical mile"; case 9036: return "kilometre"; case 9101: return "radian";
      case 9102: return "degree"; case 9103: return "arc-minute"; case 9104: return "arc-second"; case 9105: return "grad";
      case 9122: return "degree"; case 32767: return "user-defined"; default: return "code "+to_string((unsigned long long)u);
    }
  }
  //the scale, offset and unit of band 1 in GDAL's metadata (tag 42112, XML; GDAL writes them as
  //<Item name="SCALE" sample="0" role="scale">0.1</Item>): upper-case role (SCALE, OFFSET, UNITTYPE) -> text (entities
  //decoded, blanks trimmed). As GDAL reads them, only items with such a role count: an item without one is ordinary
  //band metadata
  void GDALBandItems(const string &x,map<string,string> &items)
  {
    size_t p=0;
    while ((p=x.find("<Item",p))!=string::npos){
      size_t gt=x.find('>',p); if (gt==string::npos){break;}
      const string head=x.substr(p,gt-p);
      if (!head.empty() && (head[head.size()-1]=='/')){p=gt+1; continue;}     //(an empty item)
      size_t end=x.find("</Item>",gt); if (end==string::npos){break;}
      string body=x.substr(gt+1,end-gt-1); p=end+7;
      string role,sample;
      for (int k=0;k<2;k++){
        const string a=(k==0)?" role=\"":" sample=\""; size_t q=head.find(a); if (q==string::npos){continue;}
        q+=a.size(); size_t r=head.find('"',q); if (r!=string::npos){((k==0)?role:sample)=head.substr(q,r-q);}
      }
      if (sample!="0"){continue;}                                             //(dataset items and other bands)
      for (size_t i=0;i<role.size();i++){role[i]=(char)toupper((unsigned char)role[i]);}
      if ((role!="SCALE") && (role!="OFFSET") && (role!="UNITTYPE")){continue;} //(other band metadata)
      const char *ent[5][2]={{"&lt;","<"},{"&gt;",">"},{"&quot;","\""},{"&apos;","'"},{"&amp;","&"}};
      for (int e=0;e<5;e++){size_t q=0; while ((q=body.find(ent[e][0],q))!=string::npos){body.replace(q,strlen(ent[e][0]),ent[e][1]); q+=strlen(ent[e][1]);}}
      size_t a=body.find_first_not_of(" \t\r\n"),b=body.find_last_not_of(" \t\r\n");
      body=(a==string::npos)?"":body.substr(a,b-a+1);
      items[role]=body;
    }
  }
  const char *Lookup(const gw_epsg_name *tab,const int nt,const int code)
  {
    int lo=0,hi=nt-1;
    while (lo<=hi){int mid=(lo+hi)/2; if (tab[mid].code==code){return tab[mid].name;} if (tab[mid].code<code){lo=mid+1;} else {hi=mid-1;}}
    return NULL;
  }
  //an internal mask (a GDAL per-dataset mask: an image directory with NewSubfileType 4 and the raster's size, 1 or 8 bits
  //per pixel, 0 = no data): looked for in the chain of image directories after the raster (at most MAX_IFDS); its no-data
  //pixels become NaN. A chain that is damaged, loops or goes on beyond the bound ends the search: the raster is used without
  //a mask, and `chain` says why (the caller warns). Returns false only for a mask Raven cannot read (masked: number of
  //pixels set to no data; found: a mask)
  const int MAX_IFDS=4096;
  bool ApplyMask(const unsigned char *data,const size_t n,const TIFF &T0,const uint64_t w,const uint64_t h,vector<double> &v,
                 uint64_t &masked,bool &found,string &chain,string &e)
  {
    masked=0; found=false; chain="";
    TIFF M(data,n); M.be=T0.be; M.big=T0.big;
    std::set<uint64_t> seen; seen.insert(T0.first);
    uint64_t off=T0.next;
    for (int k=0;off!=0;k++){
      if (k>=MAX_IFDS){chain="the search for an internal mask stopped after "+to_string(MAX_IFDS)+" image directories"; return true;}
      if (seen.count(off)){chain="its chain of image directories loops"; return true;}
      seen.insert(off);
      string de;
      if (!M.ReadIFD(off,de)){chain="its chain of image directories is damaged after the raster ("+de+")"; return true;}
      off=M.next;
      uint64_t st,mw,mh;
      if (!M.One(254,0,st,de) || !M.One(256,0,mw,de) || !M.One(257,0,mh,de)){chain="its chain of image directories is damaged after the raster ("+de+")"; return true;}
      if (!(st&4) || (st&1) || (mw!=w) || (mh!=h)){continue;}         //(not the raster's mask: overviews, their masks, other pages)
      found=true;
      const string cannot="the file has an internal mask (no-data pixels) that Raven cannot read (";
      const string fix="): save the raster with a NODATA value instead (e.g. gdalwarp -dstnodata -9999 in.tif out.tif)";
      uint64_t bits,spp,comp,pred,fill;
      if (!M.One(258,1,bits,de) || !M.One(277,1,spp,de) || !M.One(259,1,comp,de) || !M.One(317,1,pred,de) || !M.One(266,1,fill,de)){e=cannot+de+fix; return false;}
      if ((comp==1) || (comp==32773)){pred=1;}
      if (((bits!=1) && (bits!=8)) || (spp!=1)){e=cannot+to_string((unsigned long long)bits)+" bits, "+to_string((unsigned long long)spp)+" samples per pixel"+fix; return false;}
      if ((comp!=1) && (comp!=5) && (comp!=8) && (comp!=32946) && (comp!=32773)){e=cannot+"compression "+CompressionName(comp)+fix; return false;}
      if ((pred!=1) && !((pred==2) && (bits==8))){e=cannot+"predictor "+to_string((unsigned long long)pred)+fix; return false;}
      if (fill!=1){e=cannot+"FillOrder "+to_string((unsigned long long)fill)+fix; return false;}
      Blocks B;
      if (!Layout(M,w,h,B,de)){e=cannot+de+fix; return false;}
      if ((bits==1) && B.tiled && (B.bw%8!=0)){e=cannot+"1-bit tiles of a width that is not a multiple of 8"+fix; return false;}
      const uint64_t rb=(bits==1)?(B.bw+7)/8:B.bw;                 //(bytes per row of a block; 1-bit rows start on a byte)
      if (!Plausible(B,w,h,comp,(double)rb,de)){e=cannot+de+fix; return false;}
      vector<unsigned char> buf,tmp;
      for (uint64_t b=0;b<B.n;b++){
        uint64_t r0=(b/B.nacross)*B.bh,c0=(b%B.nacross)*B.bw,rows=B.Rows(b,h);
        if (B.cnts[(size_t)b]==0){e=cannot+"an empty block"+fix; return false;}
        if (!DecodeBlock(M,B,b,comp,pred,(size_t)(rows*rb),rows,1,buf,tmp,de)){e=cannot+de+fix; return false;}
        for (uint64_t r=0;r<rows && r0+r<h;r++){
          for (uint64_t c=0;c<B.bw && c0+c<w;c++){
            bool on=(bits==1)?((buf[(size_t)(r*rb+c/8)]>>(7-c%8))&1):(buf[(size_t)(r*rb+c)]!=0);
            size_t i=(size_t)((r0+r)*w+c0+c);
            if (!on && (v[i]==v[i])){v[i]=NAN; masked++;}
          }
        }
      }
      return true;
    }
    return true;
  }
}

bool GWGeoTIFF::IsTIFFName(const string &filename)
{
  size_t p=filename.find_last_of('.'); if (p==string::npos){return false;}
  string e=filename.substr(p+1); for (size_t i=0;i<e.size();i++){e[i]=(char)tolower((unsigned char)e[i]);}
  size_t sl=filename.find_last_of("/\\"); if ((sl!=string::npos) && (sl>p)){return false;}
  return (e=="tif") || (e=="tiff") || (e=="gtiff");
}

bool GWGeoTIFF::Read(const string &filename,gw_geotiff &R,string &err)
{
  ifstream IN(filename.c_str(),ios::binary);
  if (IN.fail()){err="cannot open raster file "+filename; return false;}
  IN.seekg(0,ios::end); streamoff sz=IN.tellg(); IN.seekg(0,ios::beg);
  if (sz<0){err="cannot read raster file "+filename; return false;}
  if ((unsigned long long)sz>MAX_FILE_BYTES){
    err="GeoTIFF file "+filename+" is "+to_string((unsigned long long)sz/(1024ULL*1024ULL))+" MB, more than the "+to_string(MAX_FILE_BYTES/(1024ULL*1024ULL))+
        " MB Raven reads: resample it to a coarser resolution (a few raster cells per model cell are enough; each cell averages 25 samples)";
    return false;
  }
  vector<unsigned char> data((size_t)sz);
  if ((sz>0) && !IN.read((char*)&data[0],sz)){err="cannot read raster file "+filename; return false;}
  return ReadMemory(data.empty()?NULL:&data[0],data.size(),filename,R,err);
}

bool GWGeoTIFF::ReadMemory(const unsigned char *data,const size_t n,const string &name,gw_geotiff &R,string &err)
{
  R=gw_geotiff();
  string e;
  const string F="GeoTIFF file "+name+": ";
  TIFF T(data,(data==NULL)?0:n);
  if (!T.Header(e) || !T.FirstIFD(e)){err=F+e; return false;}
  //-- image layout
  uint64_t w,h,spp,bps,fmt,comp,pred,photo,planar,subfile,fill;
  if (!T.One(256,0,w,e) || !T.One(257,0,h,e) || !T.One(277,1,spp,e) || !T.One(339,1,fmt,e) || !T.One(259,1,comp,e) ||
      !T.One(317,1,pred,e) || !T.One(262,1,photo,e) || !T.One(284,1,planar,e) || !T.One(254,0,subfile,e) || !T.One(266,1,fill,e)){err=F+e; return false;}
  { const TEntry *t=T.Tag(258); bps=1;
    if (t!=NULL){
      if (!T.UInt(*t,0,bps)){err=F+"invalid BitsPerSample"; return false;}
      for (uint64_t i=1;i<t->count && i<spp;i++){uint64_t b; if (!T.UInt(*t,i,b) || (b!=bps)){err=F+"samples of different bit depths"; return false;}}
    }
  }
  if ((w==0) || (h==0)){err=F+"the image width or height is zero or missing (corrupt file)"; return false;}
  if ((double)w*(double)h>MAX_PIXELS){
    err=F+"the raster has "+to_string((unsigned long long)w)+" x "+to_string((unsigned long long)h)+" pixels, more than the "+Fmt(MAX_PIXELS/1e6)+
      " million Raven reads: resample it to a coarser resolution (a few raster cells per model cell are enough; each cell averages 25 samples)";
    return false;
  }
  if (subfile&5){err=F+"the first image of the file is a reduced-resolution copy or a mask, not the raster itself"; return false;}
  if (spp!=1){err=F+"the file has "+to_string((unsigned long long)spp)+" samples (bands) per pixel; Raven reads single-band rasters: save the elevation band alone (e.g. gdal_translate -b 1)"; return false;}
  if ((photo!=0) && (photo!=1)){
    const char *pn[]={"WhiteIsZero","BlackIsZero","RGB colour","palette colour","transparency mask","CMYK colour","YCbCr colour"};
    err=F+"the photometric interpretation is "+string((photo<7)?pn[photo]:"not a grey-scale value")+" ("+to_string((unsigned long long)photo)+"); Raven reads single-band rasters of values (grey scale)"; return false;
  }
  if (fill!=1){err=F+"FillOrder "+to_string((unsigned long long)fill)+" (bits of each byte in reverse order) is not supported: re-save the file (e.g. gdal_translate)"; return false;}
  if ((fmt<1) || (fmt>3) || ((bps!=8) && (bps!=16) && (bps!=32) && (bps!=64)) || ((fmt==3) && (bps<32))){
    const char *fn[]={"","unsigned integer","signed integer","floating-point","undefined","complex integer","complex floating-point"};
    err=F+"samples of "+to_string((unsigned long long)bps)+" bits, "+string((fmt>=1 && fmt<=6)?fn[fmt]:"unknown format")+
      ", are not supported: Raven reads 8, 16, 32 and 64 bit integers and 32 and 64 bit floating-point values"; return false;
  }
  if ((comp!=1) && (comp!=5) && (comp!=8) && (comp!=32946) && (comp!=32773)){
    err=F+"compression "+CompressionName(comp)+" ("+to_string((unsigned long long)comp)+") is not supported: Raven reads uncompressed, LZW, DEFLATE and PackBits files "
      "(re-save it, e.g. gdal_translate -co COMPRESS=DEFLATE)"; return false;
  }
  if ((comp==1) || (comp==32773)){pred=1;}           //(a predictor belongs to the LZW and DEFLATE codecs; libtiff ignores it otherwise)
  if ((pred!=1) && (pred!=2) && (pred!=3)){err=F+"predictor "+to_string((unsigned long long)pred)+" is not supported (1 none, 2 horizontal, 3 floating point)"; return false;}
  if ((pred==3) && (fmt!=3)){err=F+"the floating-point predictor (3) is given for integer samples (corrupt file)"; return false;}
  const int nb=(int)(bps/8);
  //-- strips or tiles
  Blocks B;
  if (!Layout(T,w,h,B,e)){err=F+e; return false;}
  //-- georeferencing
  vector<double> scale,tie,mt;
  bool hasScale=T.Doubles(33550,scale)&&(scale.size()>=2),hasTie=T.Doubles(33922,tie)&&(tie.size()>=6),hasMT=T.Doubles(34264,mt)&&(mt.size()>=16);
  double sx,sy,xul,yul;
  if (hasScale && hasTie){
    sx=scale[0]; sy=scale[1];
    if (!std::isfinite(sx) || !std::isfinite(sy) || (sx<=0) || (sy==0)){err=F+"invalid ModelPixelScale ("+Fmt(sx)+", "+Fmt(sy)+")"; return false;}
    if (sy<0){err=F+"the raster is not north-up (rows run south to north, negative ModelPixelScale y): re-save it north-up (e.g. gdalwarp)"; return false;}
    xul=tie[3]-tie[0]*sx; yul=tie[4]+tie[1]*sy;
  }
  else if (hasMT){
    double tol=1e-10*std::max(fabs(mt[0]),fabs(mt[5]));
    if (!(fabs(mt[1])<=tol) || !(fabs(mt[4])<=tol)){err=F+"its ModelTransformation rotates or shears the raster; Raven reads north-up rasters only (re-save it north-up, e.g. gdalwarp)"; return false;}
    if (!(mt[0]>0) || !(mt[5]<0)){err=F+"the raster is not north-up (ModelTransformation scale "+Fmt(mt[0])+", "+Fmt(mt[5])+"): re-save it north-up (e.g. gdalwarp)"; return false;}
    sx=mt[0]; sy=-mt[5]; xul=mt[3]; yul=mt[7];
  }
  else if (hasTie && (tie.size()>6)){err=F+"the raster is georeferenced by control points (tie points without ModelPixelScale); Raven reads north-up rasters with a pixel size: warp it first (e.g. gdalwarp)"; return false;}
  else {err=F+"no georeferencing (ModelPixelScale and ModelTiepoint, or ModelTransformation): not a GeoTIFF"; return false;}
  if (!std::isfinite(xul) || !std::isfinite(yul) || !std::isfinite(sx) || !std::isfinite(sy)){err=F+"invalid georeferencing values"; return false;}
  //-- GeoKeys: model type, raster type (PixelIsArea/PixelIsPoint), EPSG code of the CRS
  int projCS=0,geogCS=0; bool hasProj=false,hasGeog=false;
  uint64_t projLin=0,geogLin=0,geogAng=0,vertU=0,geogPM=0,vertCS=0,geogDatum=0; //(units, prime meridian, datum and vertical keys; 0: not given)
  double pmLong=0.0; bool hasPMLong=false;   //(GeogPrimeMeridianLongGeoKey, a double held in GeoDoubleParams, tag 34736)
  { const TEntry *t=T.Tag(34735);
    if (t!=NULL){
      vector<uint64_t> k; uint64_t nk=t->count;
      if (nk>4096){nk=4096;}
      if (!T.List(34735,k,nk,e) || (k.size()<4)){err=F+"corrupt GeoKeyDirectory"; return false;}
      uint64_t nkeys=std::min(k[3],(uint64_t)(k.size()-4)/4);
      for (uint64_t i=0;i<nkeys;i++){
        uint64_t id=k[4+4*i],loc=k[5+4*i],cnt=k[6+4*i],val=k[7+4*i];
        if ((loc==34736) && (id==2061)){                //(the prime meridian of a user-defined geographic system, as GDAL writes it)
          vector<double> dbl;
          if (!T.Doubles(34736,dbl) || (cnt<1) || (val>=dbl.size())){err=F+"corrupt GeoKeyDirectory (GeogPrimeMeridianLongGeoKey outside GeoDoubleParams)"; return false;}
          pmLong=dbl[(size_t)val]; hasPMLong=true; continue;
        }
        if (loc!=0){continue;}                          //(other values held in other tags are not needed)
        if      (id==1024){R.modelType=(int)val;}
        else if (id==1025){R.rasterType=(int)val;}
        else if (id==2048){geogCS=(int)val; hasGeog=true;}
        else if (id==3072){projCS=(int)val; hasProj=true;}
        else if (id==2052){geogLin=val;}
        else if (id==2054){geogAng=val;}
        else if (id==3076){projLin=val;}
        else if (id==4099){vertU=val;}
        else if (id==2051){geogPM=val;}
        else if (id==2050){geogDatum=val;}
        else if (id==4096){vertCS=val;}
      }
    }
  }
  if (R.modelType==3){err=F+"geocentric coordinates (GTModelTypeGeoKey 3) are not supported: give the raster in longitude/latitude or a projection"; return false;}
  if ((R.rasterType!=1) && (R.rasterType!=2)){R.rasterType=1;}
  if      (R.modelType==1){if (hasProj && (projCS>0) && (projCS!=32767)){R.epsg=projCS;}}
  else if (R.modelType==2){if (hasGeog && (geogCS>0) && (geogCS!=32767)){R.epsg=geogCS; R.geographic=true;}}
  else if (R.modelType==0){
    if      (hasProj && (projCS>0) && (projCS!=32767)){R.epsg=projCS;}
    else if (!hasProj && hasGeog && (geogCS>0) && (geogCS!=32767)){R.epsg=geogCS; R.geographic=true;}
  }
  //units: Raven reads coordinates in metres (projections) or degrees, and elevations in metres
  const bool geo=(R.modelType==2) || ((R.modelType==0) && !hasProj && hasGeog);
  if ((projLin!=0) && (projLin!=9001)){err=F+"the linear units of its projection (ProjLinearUnitsGeoKey) are "+UnitName(projLin)+"; Raven reads projected coordinates in metres: reproject the raster to a projection in metres (e.g. gdalwarp -t_srs EPSG:...)"; return false;}
  if ((geogLin!=0) && (geogLin!=9001)){err=F+"the linear units of its geographic coordinate system (GeogLinearUnitsGeoKey) are "+UnitName(geogLin)+"; Raven expects metres: reproject the raster (e.g. gdalwarp -t_srs EPSG:4326)"; return false;}
  if (geo && (geogAng!=0) && (geogAng!=9102) && (geogAng!=9122)){err=F+"its angular units (GeogAngularUnitsGeoKey) are "+UnitName(geogAng)+"; Raven reads longitude/latitude in degrees: reproject the raster (e.g. gdalwarp -t_srs EPSG:4326)"; return false;}
  if (R.geographic){
    const char *pm=Lookup(GW_EPSG_GEOG_PRIME_MERIDIAN,GW_EPSG_GEOG_PRIME_MERIDIAN_N,R.epsg),*au=Lookup(GW_EPSG_GEOG_NOT_DEGREES,GW_EPSG_GEOG_NOT_DEGREES_N,R.epsg);
    if (pm!=NULL){err=F+"it is in EPSG:"+to_string(R.epsg)+", whose longitudes count from the prime meridian of "+string(pm)+" from Greenwich; Raven reads longitudes from Greenwich: reproject the raster (e.g. gdalwarp -t_srs EPSG:4326)"; return false;}
    if (au!=NULL){err=F+"it is in EPSG:"+to_string(R.epsg)+", whose angles are in "+string(au)+"; Raven reads longitude/latitude in degrees: reproject the raster (e.g. gdalwarp -t_srs EPSG:4326)"; return false;}
  }
  if (geo && hasPMLong && (pmLong!=0.0)){err=F+"its prime meridian lies at longitude "+Fmt(pmLong)+" from Greenwich (GeogPrimeMeridianLongGeoKey); Raven reads longitudes from Greenwich: reproject the raster (e.g. gdalwarp -t_srs EPSG:4326)"; return false;}
  if (geo && (geogPM!=0) && (geogPM!=8901) && !((geogPM==32767) && hasPMLong)){err=F+"its prime meridian (GeogPrimeMeridianGeoKey "+to_string((unsigned long long)geogPM)+") is not Greenwich; Raven reads longitudes from Greenwich: reproject the raster (e.g. gdalwarp -t_srs EPSG:4326)"; return false;}
  if (geo && (geogDatum!=0) && (geogDatum!=32767)){
    const char *dp=Lookup(GW_EPSG_DATUM_PRIME_MERIDIAN,GW_EPSG_DATUM_PRIME_MERIDIAN_N,(int)geogDatum);
    if (dp!=NULL){err=F+"its geodetic datum EPSG:"+to_string((unsigned long long)geogDatum)+" (GeogGeodeticDatumGeoKey) counts longitudes from the prime meridian of "+string(dp)+" from Greenwich; Raven reads longitudes from Greenwich: reproject the raster (e.g. gdalwarp -t_srs EPSG:4326)"; return false;}
  }
  if ((R.epsg>0) && !R.geographic && (projLin==0)){   //(a projected code that implies its unit; GeoTIFF 1.1 lets a writer leave ProjLinearUnitsGeoKey out)
    const char *pu=Lookup(GW_EPSG_PROJECTED_NOT_METRE,GW_EPSG_PROJECTED_NOT_METRE_N,R.epsg);
    if (pu!=NULL){err=F+"it is in EPSG:"+to_string(R.epsg)+", a projection whose coordinates are in "+string(pu)+"; Raven reads projected coordinates in metres: reproject the raster to a projection in metres (e.g. gdalwarp -t_srs EPSG:...)"; return false;}
  }
  if ((vertCS!=0) && (vertCS!=32767)){
    const char *fu=Lookup(GW_EPSG_VERTICAL_FEET,GW_EPSG_VERTICAL_FEET_N,(int)vertCS);
    if (fu!=NULL){err=F+"its elevations are in "+string(fu)+" (vertical coordinate system EPSG:"+to_string((unsigned long long)vertCS)+"); Raven reads elevations in metres: convert the values to metres and save the raster with a vertical system in metres or none (e.g. gdal_calc.py)"; return false;}
  }
  const bool vMetre=std::binary_search(GW_EPSG_VERTICAL_METRE,GW_EPSG_VERTICAL_METRE+GW_EPSG_VERTICAL_METRE_N,(int)vertCS) ||
                    std::binary_search(GW_EPSG_GEOGRAPHIC,GW_EPSG_GEOGRAPHIC+GW_EPSG_GEOGRAPHIC_N,(int)vertCS); //(a geographic 3D system: ellipsoidal heights in metres, as GDAL writes EPSG:4979)
  if ((vertCS!=0) && (vertU!=9001) && ((vertCS==32767) || !vMetre)){
    R.notes+=string(R.notes.empty()?"":"\n")+"raster "+name+" has a vertical coordinate system ("+((vertCS==32767)?string("user-defined"):"EPSG:"+to_string((unsigned long long)vertCS))+
             ") whose unit Raven does not know; its elevations are read as metres";
  }
  if ((vertU!=0) && (vertU!=9001)){err=F+"its elevations are in "+UnitName(vertU)+" (VerticalUnitsGeoKey); Raven reads elevations in metres: convert the values to metres and save the raster without that key (e.g. gdal_calc.py)"; return false;}
  //-- GDAL's band scale, offset and unit (GDAL_METADATA): the values are read unscaled, as GDAL's -unscale and GIS
  //   programs show them; a unit of feet stops, another unit than metres gives a warning
  double bscale=1.0,boffset=0.0;
  { const uint64_t MAXMETA=16ULL<<20;
    const TEntry *tm=T.Tag(42112);
    if ((tm!=NULL) && (tm->type==2) && (tm->count>MAXMETA)){
      err=F+"its GDAL_METADATA tag holds "+to_string((unsigned long long)(tm->count>>20))+" MB, more than the 16 MB Raven reads, so its band scale, offset and unit cannot be checked: "
          "save the raster without the extra metadata (e.g. gdal_translate -unscale -co PROFILE=GeoTIFF in.tif out.tif)"; return false;
    }
    const string meta=T.Ascii(42112,MAXMETA);
    map<string,string> items; if (meta!=""){GDALBandItems(meta,items);}
    for (int k=0;k<2;k++){
      const char *nm=(k==0)?"SCALE":"OFFSET"; map<string,string>::const_iterator it=items.find(nm); if (it==items.end()){continue;}
      char *endp; const double v=strtod(it->second.c_str(),&endp);
      if (it->second.empty() || (*endp!='\0') || !std::isfinite(v)){err=F+"its band "+string(nm)+" in GDAL_METADATA ("+it->second+") is not a number"; return false;}
      if (k==0){bscale=v;} else {boffset=v;}
    }
    if (bscale==0.0){err=F+"its band SCALE in GDAL_METADATA is 0 (every value would be the OFFSET): save the raster with its elevations (e.g. gdal_translate -unscale)"; return false;}
    map<string,string>::const_iterator iu=items.find("UNITTYPE");
    if (iu!=items.end()){
      string u=iu->second; for (size_t i=0;i<u.size();i++){u[i]=(char)tolower((unsigned char)u[i]);}
      const bool metre=(u=="") || (u=="m") || (u=="metre") || (u=="metres") || (u=="meter") || (u=="meters");
      const char *len[]={"ft","foot","feet","cm","centimet","mm","millimet","dm","decimet","km","kilomet","inch","yard","yd","mile"};
      bool other=(u=="in");
      for (int i=0;i<15;i++){if (u.find(len[i])!=string::npos){other=true;}}
      if (!metre && other){
        err=F+"its elevations are in "+iu->second+" (the band unit, UNITTYPE in GDAL_METADATA); Raven reads elevations in metres: convert the values to metres (e.g. gdal_calc.py) and set the unit to metres (gdal_edit.py -units m)"; return false;
      }
      if (!metre){
        R.notes+=string(R.notes.empty()?"":"\n")+"raster "+name+" has the band unit '"+iu->second+"' (UNITTYPE in GDAL_METADATA), which Raven does not know; its elevations are read as metres";
      }
    }
  }
  if (R.rasterType==2){xul-=0.5*sx; yul+=0.5*sy;}    //PixelIsPoint: the tie point is the centre of the pixel
  //-- NODATA
  string nd=T.Ascii(42113);
  { size_t a=nd.find_first_not_of(" \t\r\n"),b=nd.find_last_not_of(" \t\r\n");
    nd=(a==string::npos)?"":nd.substr(a,b-a+1); }
  if (nd!=""){
    string l=nd; for (size_t i=0;i<l.size();i++){l[i]=(char)tolower((unsigned char)l[i]);}
    if ((l=="nan") || (l=="-nan") || (l=="+nan")){R.hasNodata=true; R.nodata=NAN;}
    else {
      char *endp; double v=strtod(nd.c_str(),&endp);
      if ((endp==nd.c_str()) || (*endp!='\0')){err=F+"its GDAL_NODATA tag ("+nd+") is not a number"; return false;}
      if ((fmt==3) && (bps==32) && std::isfinite(v)){  //(32-bit float pixels: the value as a float holds it; just beyond the float range: the largest float)
        const double FM=3.4028234663852886e38;
        if (fabs(v)<=FM){v=(double)(float)v;} else if (fabs(v)<FM*(1.0+1e-6)){v=(v<0)?-FM:FM;}
      }
      R.hasNodata=true; R.nodata=v;
    }
  }
  //-- pixels
  if (!Plausible(B,w,h,comp,(double)B.bw*nb,e)){err=F+e; return false;}
  R.ncols=(int)w; R.nrows=(int)h; R.dx=sx; R.dy=sy; R.xll=xul; R.yll=yul-(double)h*sy;
  R.v.assign((size_t)(w*h),0.0);
  const bool outBE=(pred==3)?true:T.be;          //(after the floating-point predictor the samples are big-endian)
  vector<unsigned char> buf,tmp;
  for (uint64_t b=0;b<B.n;b++){
    uint64_t r0=(b/B.nacross)*B.bh,c0=(b%B.nacross)*B.bw,rows=B.Rows(b,h);
    if (B.cnts[(size_t)b]==0){                        //(an empty block of a sparse file: NODATA, as GDAL reads it)
      if (!R.hasNodata){err=F+string(B.tiled?"tile ":"strip ")+to_string((unsigned long long)b)+" is empty and the file has no NODATA value (GDAL_NODATA) to fill it with"; return false;}
      for (uint64_t r=0;r<rows && r0+r<h;r++){for (uint64_t c=0;c<B.bw && c0+c<w;c++){R.v[(size_t)((r0+r)*w+c0+c)]=R.nodata;}}
      continue;
    }
    if (!DecodeBlock(T,B,b,comp,pred,(size_t)(rows*B.bw*nb),rows,nb,buf,tmp,e)){err=F+e; return false;}
    for (uint64_t r=0;r<rows && r0+r<h;r++){
      for (uint64_t c=0;c<B.bw && c0+c<w;c++){
        double v=SampleValue(&buf[(size_t)((r*B.bw+c)*nb)],nb,(int)fmt,outBE);
        if (std::isinf(v)){v=NAN;}
        R.v[(size_t)((r0+r)*w+c0+c)]=v;
      }
    }
  }
  //-- internal mask
  uint64_t nmask=0; bool hasMask=false; string chain;
  if (!ApplyMask(data,(size_t)T.n,T,w,h,R.v,nmask,hasMask,chain,e)){err=F+e; return false;}
  if (chain!=""){R.notes+=string(R.notes.empty()?"":"\n")+"raster "+name+": "+chain+"; the raster is read without an internal mask";}
  //-- no NODATA value: values that usually mark missing data are reported (they stay values, as GDAL reads them)
  if (!R.hasNodata){
    vector<double> fills; fills.push_back(-9999.0);
    if (fmt==3){fills.push_back(-3.4028234663852886e38); fills.push_back(3.4028234663852886e38);}
    if ((fmt==3) && (bps==64)){fills.push_back(-1.7976931348623157e308); fills.push_back(1.7976931348623157e308);}
    if ((fmt==2) && (bps<=32)){fills.push_back(-ldexp(1.0,(int)bps-1));}
    if ((fmt==2) && (bps==64)){fills.push_back(-9223372036854775808.0);}
    ostringstream s; int nf=0;
    for (size_t k=0;k<fills.size();k++){
      size_t c=0; for (size_t i=0;i<R.v.size();i++){if (R.v[i]==fills[k]){c++;}}
      if (c>0){s<<(nf?", ":"")<<c<<((c==1)?" pixel of ":" pixels of ")<<Fmt(fills[k]); nf++;}
    }
    if (nf>0){
      R.notes+=string(R.notes.empty()?"":"\n")+"raster "+name+" has no NODATA value (GDAL_NODATA) but holds "+s.str()+", values that often mark missing data; they are read as elevations. "
                 "If they are no data, declare it (e.g. gdal_translate -a_nodata -9999 in.tif out.tif)";
    }
  }
  //-- band scale and offset applied (no-data pixels become NaN, so that no scaled value is taken for NODATA)
  const bool scaled=(bscale!=1.0) || (boffset!=0.0);
  if (scaled){
    for (size_t i=0;i<R.v.size();i++){
      double &v=R.v[i];
      if ((v!=v) || (R.hasNodata && (fabs(v-R.nodata)<=1e-6))){v=NAN; continue;}  //(as CGWRaster::Valid tells no data)
      v=v*bscale+boffset;
      if (!std::isfinite(v)){v=NAN;}                                            //(beyond the range of a double: no data, as infinite values are)
    }
    R.hasNodata=false; R.nodata=0.0;
  }
  (void)planar;                                   //(irrelevant with one sample per pixel)
  //-- description
  ostringstream s;
  s<<(T.big?"BigTIFF":"TIFF")<<", "<<(T.be?"big-endian":"little-endian")<<", ";
  if (B.tiled){s<<"tiled "<<B.bw<<"x"<<B.bh;} else {s<<"strips of "<<B.bh<<" rows";}
  s<<", "<<((comp==1)?string("uncompressed"):CompressionName(comp))<<", predictor "<<pred<<", "<<bps<<"-bit "<<((fmt==1)?"unsigned":(fmt==2)?"signed":"float")
   <<", "<<((R.rasterType==2)?"PixelIsPoint":"PixelIsArea");
  if (hasMask){s<<", internal mask ("<<nmask<<" pixels no data)";}
  if (scaled){s<<", scale "<<Fmt(bscale)<<" and offset "<<Fmt(boffset)<<" applied (GDAL_METADATA)";}
  R.format=s.str();
  return true;
}
