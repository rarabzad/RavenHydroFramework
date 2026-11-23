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
#include <iomanip>
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
  //a text for a message: bytes that are not UTF-8, and control characters, written as \xNN
  string Shown(const string &s){
    string o;
    for (size_t i=0;i<s.size();){
      const unsigned char c=(unsigned char)s[i]; int n=0;
      if (c<0x80){n=0;} else if ((c&0xE0)==0xC0){n=1;} else if ((c&0xF0)==0xE0){n=2;} else if ((c&0xF8)==0xF0){n=3;} else {n=-1;}
      bool ok=(n>=0) && (i+n<s.size()) && !((n==1) && (c<0xC2));
      for (int k=1;ok && (k<=n);k++){if (((unsigned char)s[i+k]&0xC0)!=0x80){ok=false;}}
      if (ok && (n==0) && ((c<0x20) || (c==0x7F))){ok=false;}
      if (!ok){char b[8]; snprintf(b,sizeof(b),"\\x%02X",c); o+=b; i++; continue;}
      o.append(s,i,n+1); i+=n+1;
    }
    return o;
  }
  //-- GDAL's metadata XML (the GDAL_METADATA tag and the side file <raster>.aux.xml) is read only in the form GDAL's own
  //   writer lays it out (CPLSerializeXMLTree, cpl_minixml.cpp, with the escapes of CPLEscapeString, cpl_string.cpp): each
  //   element on its own line, indented by two blanks a level; attributes as ` name="value"`, in double quotes, one blank
  //   before each; an element with only text as <name>text</name>, one without content as <name />, one with child elements
  //   as <name>, its children, then </name> at its own indentation. In text GDAL writes &lt; &gt; &amp; for < > & (in an
  //   attribute also &quot; for "), &#xFEFF; for a byte-order mark, and leaves out other control characters; no CDATA or
  //   other entities. In the metadata of an xml: domain, which GDAL keeps as parsed XML, it also writes comments
  //   (<!--...-->), literals (<!DOCTYPE ...>) and processing instructions (<?name attributes?>), each on its own line at
  //   the indentation of an element, and mixed content: a text where it stands, with no line end of its own, the first child
  //   on a new line, and the closing tag at its indentation after the last text (CPLSerializeXMLNode); Raven reads nothing
  //   in them. A file in any other form may be read differently by GDAL (its parser takes
  //   other quoting, entities and layouts, and ignores a file that does not start with its element), so the caller stops
  //   when such a file can hold band values
  struct gx_node {string name; vector<pair<string,string> > attr; bool text; string value; vector<gx_node> kids; gx_node():text(false){}};
  bool GXNameChar(const char c){return isalnum((unsigned char)c) || (c=='_') || (c=='-') || (c=='.') || (c==':') || (c=='#');}
  //the text of an element or attribute value as GDAL escapes it: decoded (false: not GDAL's escaping)
  bool GXUnescape(const string &s,const bool attr,string &out)
  {
    const char *ent[5]={"&lt;","&gt;","&amp;","&quot;","&#xFEFF;"},*chr[5]={"<",">","&","\"","\xEF\xBB\xBF"};
    out.clear();
    for (size_t i=0;i<s.size();){
      const unsigned char c=(unsigned char)s[i];
      if (c=='&'){
        int k=0; for (k=0;k<5;k++){if (s.compare(i,strlen(ent[k]),ent[k])==0){break;}}
        if ((k==5) || ((k==3) && !attr)){return false;}
        out+=chr[k]; i+=strlen(ent[k]); continue;
      }
      if ((c=='<') || (c=='>') || (attr && (c=='"')) || ((c<0x20) && (c!='\t') && (c!='\n') && (c!='\r'))){return false;}
      if ((c==0xEF) && (s.compare(i,3,"\xEF\xBB\xBF")==0)){return false;}
      out+=s[i++];
    }
    return true;
  }
  //a comment or a literal at an indentation, and its line end, as GDAL writes them (above): passed over
  bool GXSkip(const string &x,size_t &p,const size_t ind)
  {
    if ((x.compare(p,ind,string(ind,' '))!=0) || (x.compare(p+ind,2,"<!")!=0)){return false;}
    size_t q=p+ind;
    if (x.compare(q,4,"<!--")==0){q=x.find("-->",q+4); if (q==string::npos){return false;} q+=3;}
    else { //(a literal ends at the first > outside its brackets: <!DOCTYPE x [ <!ENTITY ...> ]>)
      int d=0; for (q+=2;q<x.size();q++){if (x[q]=='['){d++;} else if (x[q]==']'){d--;} else if ((x[q]=='>') && (d<=0)){break;}}
      if (q>=x.size()){return false;} q++;
    }
    if (x.compare(q,1,"\n")!=0){return false;}
    p=q+1; return true;
  }
  bool GXElement(const string &x,size_t &p,const size_t ind,gx_node &n,string &why)
  {
    const size_t at=p;
    const bool pi=(p+ind+2<x.size()) && (x.compare(p+ind,2,"<?")==0); //(a processing instruction: <?name attributes?>)
    if ((x.compare(p,ind,string(ind,' '))!=0) || (x.compare(p+ind,1,"<")!=0) || (p+ind+1>=x.size()) || !GXNameChar(x[p+ind+(pi?2:1)])){
      why="byte "+to_string((unsigned long long)at)+": an element at an indentation of "+to_string((unsigned long long)ind)+" blanks expected"; return false;
    }
    p+=ind+(pi?2:1); size_t q=p; while ((q<x.size()) && GXNameChar(x[q])){q++;}
    n.name=x.substr(p,q-p); p=q;
    while ((p+1<x.size()) && (x[p]==' ') && GXNameChar(x[p+1])){
      p++; q=p; while ((q<x.size()) && GXNameChar(x[q])){q++;}
      const string a=x.substr(p,q-p); p=q;
      if (x.compare(p,2,"=\"")!=0){why="element "+n.name+", attribute "+a+": =\"value\" expected"; return false;}
      p+=2; q=x.find('"',p); string v;
      if ((q==string::npos) || !GXUnescape(x.substr(p,q-p),true,v)){why="element "+n.name+", attribute "+a+": a value GDAL does not write"; return false;}
      for (size_t k=0;k<n.attr.size();k++){if (n.attr[k].first==a){why="element "+n.name+": attribute "+a+" given twice"; return false;}}
      n.attr.push_back(make_pair(a,v)); p=q+1;
    }
    if (pi){
      if (x.compare(p,3,"?>\n")!=0){why="processing instruction "+n.name+": ?> and a line end expected at byte "+to_string((unsigned long long)p); return false;}
      n.name="?"+n.name; p+=3; return true;
    }
    if (x.compare(p,4," />\n")==0){p+=4; return true;}
    if (x.compare(p,1,">")!=0){why="element "+n.name+": > expected at byte "+to_string((unsigned long long)p); return false;}
    p++;
    const string close="</"+n.name+">\n";
    const size_t lt=x.find('<',p);
    const bool textOnly=(lt!=string::npos) && (x.compare(lt,close.size(),close)==0); //(a text may start with a line end: no '<' in it)
    if (!textOnly){ //(child elements, comments and literals at the next indentation, a text before, between or after them,
                    // then the closing tag at this indentation)
      const string in0(ind,' '),in2(ind+2,' '); string tx;
      for (bool first=true;;first=false){
        const size_t lq=x.find('<',p);
        if (lq==string::npos){why="element "+n.name+" is not closed"; return false;}
        if (!first && (lq>=p+ind) && (x.compare(lq-ind,ind,in0)==0) && (x.compare(lq,close.size(),close)==0)){
          if (!GXUnescape(x.substr(p,lq-ind-p),false,tx)){why="element "+n.name+": a text holds a character or entity GDAL does not write"; return false;}
          p=lq+close.size(); return true;
        }
        if ((lq<p+ind+2) || (x.compare(lq-ind-2,ind+2,in2)!=0) || (first && ((lq==p+ind+2) || (x[lq-ind-3]!='\n')))){
          why="byte "+to_string((unsigned long long)p)+": an element at an indentation of "+to_string((unsigned long long)(ind+2))+" blanks expected"; return false;
        }
        if (!GXUnescape(x.substr(p,lq-ind-2-p-(first?1:0)),false,tx)){why="element "+n.name+": a text holds a character or entity GDAL does not write"; return false;}
        p=lq-ind-2;
        if (x.compare(lq,2,"<!")==0){
          if (!GXSkip(x,p,ind+2)){why="byte "+to_string((unsigned long long)lq)+": a comment or literal of element "+n.name+" not in the form GDAL writes"; return false;}
          continue;
        }
        gx_node c; if (!GXElement(x,p,ind+2,c,why)){return false;}
        n.kids.push_back(c);
      }
    }
    q=x.find('<',p);
    if ((q==string::npos) || (x.compare(q,close.size(),close)!=0)){why="element "+n.name+": its text is not followed by </"+n.name+"> and a line end"; return false;}
    if (!GXUnescape(x.substr(p,q-p),false,n.value)){why="element "+n.name+": its text holds a character or entity GDAL does not write"; return false;}
    n.text=true; p=q+close.size(); return true;
  }
  //the whole text is one element in GDAL's layout (a missing line end after the last tag is allowed: the TIFF tag has none)
  bool GXDocument(const string &x0,gx_node &root,string &why)
  {
    const string x=(!x0.empty() && (x0[x0.size()-1]!='\n'))?x0+"\n":x0;
    size_t p=0; if (!GXElement(x,p,0,root,why)){return false;}
    if (p!=x.size()){why="text after the closing </"+root.name+">"; return false;}
    return true;
  }
  //a name GDAL reads in any case (CPLGetXMLNode, EQUAL) written in another case than GDAL writes it
  bool GXOtherCase(const string &a,const string &gdal)
  {
    if ((a.size()!=gdal.size()) || (a==gdal)){return false;}
    for (size_t i=0;i<a.size();i++){if (tolower((unsigned char)a[i])!=tolower((unsigned char)gdal[i])){return false;}}
    return true;
  }
  bool GXInteger(const string &s){return !s.empty() && (s.find_first_not_of("0123456789")==string::npos) && ((s=="0") || (s[0]!='0'));} //(as GDAL writes %d)
  const gx_node *GXChild(const gx_node &n,const string &name){for (size_t i=0;i<n.kids.size();i++){if (n.kids[i].name==name){return &n.kids[i];}} return NULL;}
  bool GXHolds(const string &x,const char **words,const int nw)
  {
    string l=x; for (size_t i=0;i<l.size();i++){l[i]=(char)tolower((unsigned char)l[i]);}
    if (l.find('\0')!=string::npos){return true;}
    for (int k=0;k<nw;k++){if (l.find(words[k])!=string::npos){return true;}}
    return false;
  }
  //the scale, offset and unit of band 1 in GDAL's metadata (tag 42112), as GDAL writes them:
  //  <GDALMetadata>
  //    <Item name="SCALE" sample="0" role="scale">0.1</Item>
  //  </GDALMetadata>
  //(an item's attributes name, sample, role, domain in this order; its value escaped twice, as GDAL does). Only items of
  //sample 0 with the role scale, offset or unittype count (as GDAL reads them; any other item is ordinary metadata):
  //upper-case role -> text. Returns 0 when read, 1 when the text is not in GDAL's form but names no scale, offset or unit
  //(GDAL takes none from it), 2 when it is not in GDAL's form and names one (why in e)
  int GDALBandItems(const string &x,map<string,string> &items,string &e)
  {
    gx_node r; string why; bool ok=GXDocument(x,r,why);
    if (!ok && (x.find("\r\n")!=string::npos)){why="it has CRLF line ends; GDAL writes LF ("+why+")";}
    if (ok && (r.name!="GDALMetadata")){ok=false; why="its element is "+r.name+", not GDALMetadata";}
    for (size_t i=0;ok && (i<r.kids.size());i++){
      const gx_node &it=r.kids[i];
      if ((it.name!="Item") || !it.kids.empty() || it.attr.empty() || (it.attr[0].first!="name")){ok=false; why="an element that is not <Item name=...>"; break;}
      const char *order[4]={"name","sample","role","domain"}; int o=0; string sample,role;
      for (size_t k=0;ok && (k<it.attr.size());k++){
        while ((o<4) && (it.attr[k].first!=order[o])){o++;}
        if (o==4){ok=false; why="item "+it.attr[0].second+": attribute "+it.attr[k].first+" (GDAL writes name, sample, role and domain, in this order)"; break;}
        if (o==1){sample=it.attr[k].second; if (!GXInteger(sample)){ok=false; why="item "+it.attr[0].second+": sample=\""+sample+"\" is not an integer as GDAL writes it";}}
        if (o==2){role=it.attr[k].second;}
        o++;
      }
      if (ok && (GXOtherCase(role,"scale") || GXOtherCase(role,"offset") || GXOtherCase(role,"unittype"))){ok=false; why="item "+it.attr[0].second+": role=\""+role+"\" (GDAL writes it in lower case)"; break;}
      if (!ok || (sample!="0") || ((role!="scale") && (role!="offset") && (role!="unittype"))){continue;}
      string v; if (!GXUnescape(it.value,true,v)){ok=false; why="item "+it.attr[0].second+": its value is not escaped as GDAL escapes it"; break;}
      string R=role; for (size_t k=0;k<R.size();k++){R[k]=(char)toupper((unsigned char)R[k]);}
      if (items.count(R)){ok=false; why="two items of band 1 with the role "+role; break;}
      items[R]=v;
    }
    if (ok){return 0;}
    items.clear();
    const char *w[4]={"scale","offset","unittype","&#"};
    if (!GXHolds(x,w,4)){return 1;}
    e=why; return 2;
  }
  //GDAL's side file <raster>.aux.xml (PAM: what GDAL keeps outside the raster file), as GDAL writes it:
  //  <PAMDataset>
  //    <PAMRasterBand band="1">
  //      <UnitType>ft</UnitType>
  //      <Offset>5</Offset>
  //      <Scale>0.01</Scale>
  //      <NoDataValue le_hex_equiv="...">-9999</NoDataValue>
  //    </PAMRasterBand>
  //  </PAMDataset>
  //(PAMRasterBand elements only under PAMDataset, one per band, band="n" an integer; in band 1 each of UnitType, Offset,
  //Scale and NoDataValue at most once, with text). Band 1's items -> upper-case name (SCALE, OFFSET, UNITTYPE, NODATA) ->
  //text (NODATA from le_hex_equiv when it is given, as GDAL takes it). Returns 0 when read, 1 when the file is not in
  //GDAL's form but names no band and no band value (GDAL takes no band value from such a file), 2 when it is not in GDAL's
  //form and names one (why in e)
  int PAMBandItems(const string &x,map<string,string> &items,string &e)
  {
    gx_node r; string why; bool ok=GXDocument(x,r,why);
    if (!ok && (x.find("\r\n")!=string::npos)){why="it has CRLF line ends; GDAL writes LF, so convert them (e.g. dos2unix) or let GDAL write the file again";}
    if (ok && (r.name!="PAMDataset")){ok=false; why="its element is "+r.name+", not PAMDataset";}
    set<string> bands; const gx_node *b1=NULL;
    for (size_t i=0;ok && (i<r.kids.size());i++){
      const gx_node &b=r.kids[i];
      if (GXOtherCase(b.name,"PAMRasterBand")){ok=false; why="an element "+b.name+" (GDAL writes PAMRasterBand)"; break;}
      if (b.name!="PAMRasterBand"){ //(no band element below another element)
        vector<const gx_node*> st(1,&b);
        while (ok && !st.empty()){
          const gx_node *c=st.back(); st.pop_back();
          if ((c->name=="PAMRasterBand") || GXOtherCase(c->name,"PAMRasterBand")){ok=false; why="a PAMRasterBand inside "+b.name;}
          for (size_t k=0;k<c->kids.size();k++){st.push_back(&c->kids[k]);}
        }
        continue;
      }
      if ((b.attr.size()!=1) || (b.attr[0].first!="band") || !GXInteger(b.attr[0].second)){ok=false; why="a PAMRasterBand without band=\"n\" as GDAL writes it"; break;}
      if (!bands.insert(b.attr[0].second).second){ok=false; why="two PAMRasterBand elements of band "+b.attr[0].second; break;}
      if (b.attr[0].second=="1"){b1=&b;}
    }
    if (ok && (b1!=NULL)){
      const char *nm[4][2]={{"Scale","SCALE"},{"Offset","OFFSET"},{"UnitType","UNITTYPE"},{"NoDataValue","NODATA"}};
      for (int k=0;ok && (k<4);k++){
        int cnt=0;
        for (size_t i=0;i<b1->kids.size();i++){
          if (b1->kids[i].name==nm[k][0]){cnt++;}
          if (GXOtherCase(b1->kids[i].name,nm[k][0])){ok=false; why="band 1's element "+b1->kids[i].name+" (GDAL writes "+string(nm[k][0])+")";}
        }
        if (!ok){break;}
        if (cnt==0){continue;}
        const gx_node *c=GXChild(*b1,nm[k][0]);
        const bool hex=(k==3) && (c->attr.size()==1) && (c->attr[0].first=="le_hex_equiv") && (c->attr[0].second.size()==16) &&
                       (c->attr[0].second.find_first_not_of("0123456789ABCDEFabcdef")==string::npos);
        if ((cnt>1) || !c->text || c->value.empty() || (!c->attr.empty() && !hex)){ok=false; why="band 1's "+string(nm[k][0])+" is not written as GDAL writes it"; break;}
        items[nm[k][1]]=c->value;
        if (hex){ //(GDAL takes the value's IEEE bits from le_hex_equiv: little-endian bytes)
          const string &hx=c->attr[0].second; uint64_t u=0;
          for (int j=7;j>=0;j--){u=(u<<8)|(uint64_t)strtol(hx.substr(2*j,2).c_str(),NULL,16);}
          double dv; memcpy(&dv,&u,8); ostringstream os; os<<setprecision(17)<<dv; items[nm[k][1]]=os.str();
        }
      }
    }
    if (ok){return 0;}
    items.clear();
    const char *w[5]={"pamrasterband","scale","offset","unittype","nodatavalue"};
    if (!GXHolds(x,w,5)){return 1;}
    e=why; return 2;
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

//---------------------------------------------------------------
// band unit (UnitType: GDAL_METADATA or the .aux.xml side file): read as metres only when the whole text is a metre form
//---------------------------------------------------------------
namespace
{
  //The band unit is free text, so its rule is closed: the raster is read as metres only when the whole text, normalised,
  //is one listed metre form (below); a text without a unit keeps the raster as metres; any other text stops the run
  //(GWGeoTIFF::BandUnit). Normalising: full-width characters read as ASCII, no-break and ideographic spaces as blanks,
  //blanks trimmed and runs of them read as one, lower case (ASCII, Latin-1, Latin Extended-A, Greek with its accents
  //removed, Cyrillic), and surrounding brackets ( ), [ ], { } removed. The lists of other length units only choose the
  //message: a text with another length unit and no metre gives "its elevations are in ...", which no confirmation lifts.
  //
  //metre forms: the metre written out (English, French, Spanish, Portuguese, Italian, German, Dutch, Scandinavian, Polish,
  //Czech, Hungarian, Romanian, Finnish, Russian, Ukrainian, Greek), its symbol, the metre-above-sea-level abbreviations,
  //and EPSG 9001
  const char *BU_METRE_WORD[]={"metre","metres","meter","meters","metern","m\xC3\xA8tre","m\xC3\xA8tres","m\xC3\xA9ter","m\xC3\xA9terben",
    "metro","metros","metri","metr","metra","metry","metr\xC3\xB3w","metrach","metrami","metrze","metrem","metr\xC5\xAF",
    "metrech","metru","metrin","metri\xC3\xA4","metrit",
    "\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80","\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80\xD0\xB0","\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80\xD0\xBE\xD0\xB2",
    "\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80\xD1\x8B","\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80\xD0\xB0\xD1\x85","\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80\xD0\xB0\xD0\xBC\xD0\xB8",
    "\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80\xD0\xB5","\xD0\xBC\xD0\xB5\xD1\x82\xD1\x80\xD1\x96\xD0\xB2",       //(метр, метра, метров, метры, метрах, метрами, метре, метрів)
    "\xCE\xBC\xCE\xAD\xCF\x84\xCF\x81\xCE\xBF","\xCE\xBC\xCE\xAD\xCF\x84\xCF\x81\xCE\xB1"};                       //(μέτρο, μέτρα)
  const char *BU_METRE_ABBR[]={"m","masl","m asl","m a.s.l.","m.a.s.l.","mamsl","m amsl","m a.m.s.l.","moh","m.o.h.","m\xC3\xB6h","m \xC3\xB6.h.",
    "m\xC3\xBC","m\xC3\xBCnn","m \xC3\xBC. nn","m \xC3\xBC.nn","m \xC3\xBC. m.","m \xC3\xBC.m.","\xD0\xBC","\xD0\xBC\xD0\xB1\xD1\x81","\xD0\xBC \xD0\xB1\xD1\x81", //(..., möh, mü, müNN, m ü. M., м, мБС: Baltic system)
    "m.\xC3\xB6.h.","m o.h.","mt","mts","mtr","mtrs",                                                            //(m.ö.h., m o.h., mt, mts, mtr, mtrs)
    "m.","epsg:9001","epsg::9001","epsg 9001","urn:ogc:def:uom:epsg::9001","http://www.opengis.net/def/uom/epsg/0/9001",
    "https://www.opengis.net/def/uom/epsg/0/9001"};
  //vertical datums and references that may follow a metre form ("m MSL", "m (NAVD88)", "m, AHD", "Meter Ã¼ber NN", fused
  //"mAHD", "msnm", "m.s.n.m."; a datum may end in "geoid": "m (EGM2008 geoid)")
  const char *BU_DATUM[]={"msl","amsl","asl","a.s.l.","a.m.s.l.","above sea level","above mean sea level","ahd","aod","odn","nap","+nap","taw",
    "nn","nhn","\xC3\xBC. nn","\xC3\xBC.nn","\xC3\xBC nn","\xC3\xBC\x62\x65r nn","\xC3\xBC. nhn","\xC3\xBC.nhn","\xC3\xBC\x62\x65r nhn","\xC3\xBC. m.","\xC3\xBC.m.",
    "nn2000","rh2000","dvr90","navd88","navd 88","ngvd29","ngvd 29","cgvd2013","cgvd28","egm96","egm2008","s.l.m.","slm","n.p.m.","npm",
    "s.n.m.","snm","n. m.","n.m.","s. n. m.","ellipsoidal","above ellipsoid","above the ellipsoid","ellipsoid","wgs84 ellipsoid","the wgs84 ellipsoid",
    "egm08","od","mllw","mlw","mhw","mhhw","lat","hae","orthometric","ngf","ign69","n2000","nzvd2016","t.p.",
    "above sea-level","above mean sea-level","navd","ngvd","evrf2007","evrf2000","dhhn2016","dhhn92","nhn2016","\xC3\xBC. nhn2016","\xC3\xBC.nhn2016",
    "ln02","\xC3\xBC. a.","\xC3\xBC.a.","bpv","n. m. bpv","n.m. bpv","ngf-ign69","ngf ign69","t.o.v. nap","n60","dnn","ahd71","cd","chart datum",
    "osgm15","ordnance datum","geoid","wgs84","tp"};                                                          //(..., m ü. NHN2016, m ü. A.)
  //quantities that may precede a metre form in brackets or after "in": "elevation (m)", "Elevation [m]", "height in metres"
  const char *BU_QUANT[]={"elevation","elevations","height","heights","altitude","elevation above sea level","elevation above mean sea level",
    "height above sea level","height above mean sea level","height above ellipsoid","ellipsoidal height","orthometric height","height_above_ellipsoid",
    "elev","elev.","z","dem","dtm","dsm","h\xC3\xB6he","altitud","altura","cota","hauteur","\xC3\xA9l\xC3\xA9vation","quota","altezza","altitudine","hoogte",
    "wysoko\xC5\x9B\xC4\x87","v\xC3\xBD\xC5\xA1ka","nadmo\xC5\x99sk\xC3\xA1 v\xC3\xBD\xC5\xA1ka","magass\xC3\xA1g","korkeus","h\xC3\xB6jd","h\xC3\xB8yde","h\xC3\xB8jde",
    "\xC3\xAEn\xC4\x83l\xC8\x9Bime","\xD0\xB2\xD1\x8B\xD1\x81\xD0\xBE\xD1\x82\xD0\xB0","\xD0\xB2\xD0\xB8\xD1\x81\xD0\xBE\xD1\x82\xD0\xB0",
    "\xCF\x85\xCF\x88\xCF\x8C\xCE\xBC\xCE\xB5\xCF\x84\xCF\x81\xCE\xBF","\xCF\x8D\xCF\x88\xCE\xBF\xCF\x82"};
    //(..., Höhe, élévation, wysokość, výška, nadmořská výška, magasság, höjd, høyde, højde, înălțime, высота, висота, υψόμετρο, ύψος)
  //parts and multiples, which with a metre make a scaled metre ("m (in thousands)", "tenths of a metre")
  const char *BU_PART[]={"tenth","tenths","hundredth","hundredths","thousandth","thousandths","tens","hundred","hundreds","thousand","thousands"};
  const char *BU_CJK_METRE[]={"\xE7\xB1\xB3","\xE5\x85\xAC\xE5\xB0\xBA","\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB",
    "\xE3\x83\xA1\xE3\x83\xBC\xE3\x82\xBF\xE3\x83\xBC","\xEB\xAF\xB8\xED\x84\xB0"};                       //(米, 公尺, メートル, メーター, 미터)
  const char *BU_CJK_NONE[]={"\xE7\xB1\xB3\xE5\x9B\xBD","\xE7\xB1\xB3\xE8\xBB\x8D"};                      //(米国, 米軍: the United States, no unit)
  //other lengths (for the message only): feet, inches, yards, miles, fathoms, chains in English and other languages and
  //case forms; the metre's multiples and parts as symbols (also Cyrillic); feet fused with a qualifier (BUFusedFeet)
  const char *BU_OTHER[]={"ft","foot","feet","fasl","fmsl","inch","inches","yard","yards","yd","yds","mile","miles","kmi","nmi",
    "fathom","fathoms","chain","chains","furlong","furlongs",
    "pie","pies","pied","pieds","piede","piedi","p\xC3\xA9","p\xC3\xA9s","voet","voeten","fu\xC3\x9F","f\xC3\xBC\xC3\x9F\x65","fuss","fot","fod","stopa","stopy",
    "st\xC3\xB3p","stopach","stop\xC3\xA1\x63h","jalka","jalkaa","l\xC3\xA1\x62","l\xC3\xA1\x62\x62\x61n",
    "zoll","pouce","pouces","pulgada","pulgadas","pollice","pollici","duim","yarda","yardas","milla","millas","meile","meilen",
    "mijl","miglio","miglia",
    "polegada","polegadas","milha","milhas","jarda","jardas","cal","cale","cali","mila","mile","mil","jard","jardy","jard\xC3\xB3w","palec","palce","palc\xC5\xAF",
    "m\xC3\xADle","h\xC3\xBCvelyk","m\xC3\xA9rf\xC3\xB6ld","tuuma","tuumaa","maili","mailia","jaardi","jaardia","tum","tomme","tommer","picior","picioare",
    "\xC8\x9Boli","\xC5\xA3oli","toli","inci","yarzi","mil\xC4\x83","milles","verge","verges","iarda","iarde","vara","varas",
    //(..., jardów, palců, míle, hüvelyk, mérföld, țoli, ţoli, milă)
    "braza","brazas","bra\xC3\xA7\x61","bra\xC3\xA7\x61s","brasse","brasses","toise","toises","faden","klafter","vadem","vademen","famn","famnar","favn","favne",
    "s\xC4\x85\xC5\xBC\x65\xC5\x84","s\xC3\xA1h","syli","st\xC3\xA2njen",                                   //(fathoms: ..., braça(s), sążeń, sáh, stânjen)
    "\xCE\xAF\xCE\xBD\xCF\x84\xCF\x83\xCE\xB1","\xCE\xAF\xCE\xBD\xCF\x84\xCF\x83\xCE\xB5\xCF\x82","\xCE\xBC\xCE\xAF\xCE\xBB\xCE\xB9","\xCE\xBC\xCE\xAF\xCE\xBB\xCE\xB9\xCE\xB1",
    "\xCE\xB3\xCE\xB9\xCE\xAC\xCF\x81\xCE\xB4\xCE\xB1","\xCE\xB3\xCE\xB9\xCE\xAC\xCF\x81\xCE\xB4\xCE\xB5\xCF\x82","\xCE\xBF\xCF\x81\xCE\xB3\xCE\xB9\xCE\xAC",
    //(ίντσα, ίντσες, μίλι, μίλια, γιάρδα, γιάρδες, οργιά)
    "\xD1\x81\xD0\xB0\xD0\xB6\xD0\xB5\xD0\xBD\xD1\x8C","\xD1\x81\xD0\xB0\xD0\xB6\xD0\xB5\xD0\xBD\xD0\xB8","\xD1\x81\xD0\xB0\xD0\xB6\xD0\xB5\xD0\xBD\xD0\xB5\xD0\xB9",
    "\xD0\xB2\xD0\xB5\xD1\x80\xD1\x81\xD1\x82\xD0\xB0","\xD0\xB2\xD0\xB5\xD1\x80\xD1\x81\xD1\x82\xD1\x8B","\xD0\xB2\xD1\x91\xD1\x80\xD1\x81\xD1\x82","\xD0\xB2\xD0\xB5\xD1\x80\xD1\x81\xD1\x82",
    "\xD0\xB0\xD1\x80\xD1\x88\xD0\xB8\xD0\xBD","\xD0\xB0\xD1\x80\xD1\x88\xD0\xB8\xD0\xBD\xD0\xB0","\xD0\xB0\xD1\x80\xD1\x88\xD0\xB8\xD0\xBD\xD0\xBE\xD0\xB2",
    //(сажень, сажени, саженей, верста, версты, вёрст, верст, аршин, аршина, аршинов)
    "micron","microns","rod","rods","perch","perches","league","leagues","kft",
    "\xD1\x84\xD1\x83\xD1\x82","\xD1\x84\xD1\x83\xD1\x82\xD0\xB0","\xD1\x84\xD1\x83\xD1\x82\xD0\xBE\xD0\xB2","\xD1\x84\xD1\x83\xD1\x82\xD1\x8B", //(фут, фута, футов, футы,
    "\xD1\x84\xD1\x83\xD1\x82\xD0\xB0\xD1\x85","\xD1\x84\xD1\x83\xD1\x82\xD1\x96\xD0\xB2",                                                   // футах, футів)
    "\xD0\xB4\xD1\x8E\xD0\xB9\xD0\xBC","\xD0\xB4\xD1\x8E\xD0\xB9\xD0\xBC\xD0\xB0\xD1\x85","\xD1\x8F\xD1\x80\xD0\xB4","\xD1\x8F\xD1\x80\xD0\xB4\xD0\xB0\xD1\x85",
    "\xD0\xBC\xD0\xB8\xD0\xBB\xD1\x8F","\xD0\xBC\xD0\xB8\xD0\xBB\xD0\xB8","\xD0\xBC\xD0\xB8\xD0\xBB\xD1\x8C","\xD0\xBC\xD0\xB8\xD0\xBB\xD1\x8F\xD1\x85", //(дюйм, дюймах, ярд, ярдах, миля, мили, миль, милях)
    "\xCF\x80\xCF\x8C\xCE\xB4\xCE\xB9","\xCF\x80\xCF\x8C\xCE\xB4\xCE\xB9\xCE\xB1",                               //(πόδι, πόδια)
    "km","cm","mm","dm","nm","kms","cms","mms","\xC2\xB5m","\xCE\xBCm",                                     //(..., µm, μm)
    "\xD0\xBA\xD0\xBC","\xD1\x81\xD0\xBC","\xD0\xBC\xD0\xBC","\xD0\xB4\xD0\xBC","\xD0\xBC\xD0\xBA\xD0\xBC","\xD0\xBD\xD0\xBC"}; //(км, см, мм, дм, мкм, нм)
  //short ids that are also ordinary words (PROJ's ch, mi, fath, link, in; the decametre and hectometre dam, hm): another
  //length only as the whole text, or after us, ind, int or intl ("us-in")
  const char *BU_AMBIG[]={"in","ch","mi","fath","link","links","dam","hm"};
  const char *BU_QUAL[]={"us","ind","int","intl"};
  //SI prefixes written as words: another length joined to a metre word or right before one ("kilometres", "Kilo Meter")
  const char *BU_PREFIX[]={"kilo","kil\xC3\xB3","quilo","quil\xC3\xB4","quil\xC3\xB3","chilo","hecto","hekto","hect\xC3\xB3","hect\xC3\xB4","deca","deka","d\xC3\xA9\x63\x61","d\xC3\xA9ka",
    "dec\xC3\xA1","deci","d\xC3\xA9\x63i","dec\xC3\xAD","dezi","decy","desi",
    "centi","cent\xC3\xAD","centy","zenti","sentti","santi","senti","etto","\xCF\x87\xCE\xB9\xCE\xBB\xCE\xB9\xCE\xBF","\xCE\xB5\xCE\xBA\xCE\xB1\xCF\x84\xCE\xBF\xCF\x83\xCF\x84\xCE\xBF","\xCF\x87\xCE\xB9\xCE\xBB\xCE\xB9\xCE\xBF\xCF\x83\xCF\x84\xCE\xBF","\xCE\xB4\xCE\xB5\xCE\xBA\xCE\xB1\xCF\x84\xCE\xBF","\xCE\xB4\xCE\xB5\xCE\xBA\xCE\xB1","\xCE\xB5\xCE\xBA\xCE\xB1\xCF\x84\xCE\xBF","milli","mili","mil\xC3\xAD","micro","mikro","nano","mega","giga",
    "\xD0\xBA\xD0\xB8\xD0\xBB\xD0\xBE","\xD0\xB3\xD0\xB5\xD0\xBA\xD1\x82\xD0\xBE","\xD0\xB4\xD0\xB5\xD0\xBA\xD0\xB0","\xD0\xB4\xD0\xB5\xD1\x86\xD0\xB8",
    "\xD1\x81\xD0\xB0\xD0\xBD\xD1\x82\xD0\xB8","\xD0\xBC\xD0\xB8\xD0\xBB\xD0\xBB\xD0\xB8","\xD0\xBC\xD0\xB8\xD0\xBA\xD1\x80\xD0\xBE",
    "\xD0\xBD\xD0\xB0\xD0\xBD\xD0\xBE","\xD0\xBC\xD0\xB5\xD0\xB3\xD0\xB0","\xD0\xBA\xD1\x96\xD0\xBB\xD0\xBE","\xD0\xBC\xD1\x96\xD0\xBB\xD1\x96"}; //(кило, гекто, дека, деци, санти, милли, микро, нано, мега, кіло, мілі)
  //Chinese, Japanese and Korean other lengths (searched in a run of these scripts, the longest entry first at each place)
  const char *BU_CJK_OTHER[]={"\xE5\x8E\x98\xE7\xB1\xB3","\xE6\xAF\xAB\xE7\xB1\xB3","\xE5\x88\x86\xE7\xB1\xB3","\xE5\x8D\x83\xE7\xB1\xB3",
    "\xE5\x85\xAC\xE9\x87\x8C","\xE5\xBE\xAE\xE7\xB1\xB3","\xE7\xBA\xB3\xE7\xB1\xB3","\xE7\xB4\x8D\xE7\xB1\xB3","\xE5\x85\xAC\xE5\x88\x86",
    "\xE8\x8B\xB1\xE5\xB0\xBA","\xE5\x91\x8E","\xE8\x8B\xB1\xE5\xAF\xB8","\xE5\x90\x8B","\xE8\x8B\xB1\xE9\x87\x8C","\xE7\xA0\x81","\xE7\xA2\xBC",
    "\xE5\xB0\xBA","\xE5\xB8\x82\xE5\xB0\xBA","\xE5\xAF\xB8","\xE4\xB8\x88","\xE9\x87\x8C","\xE9\x96\x93",                            //(尺, 市尺, 寸, 丈, 里, 間)
    "\xE6\xB5\xB7\xE9\x87\x8C","\xE9\x87\x90\xE7\xB1\xB3","\xE5\x85\xAC\xE9\x87\x90","\xE7\x99\xBE\xE7\xB1\xB3","\xE5\x8D\x81\xE7\xB1\xB3","\xE5\x85\xAC\xE5\x8E\x98",
    //(厘米, 毫米, 分米, 千米, 公里, 微米, 纳米, 納米, 公分, 英尺, 呎, 英寸, 吋, 英里, 码, 碼, 海里, 釐米, 公釐, 百米, 十米)
    "\xE3\x82\xAD\xE3\x83\xAD\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB","\xE3\x82\xBB\xE3\x83\xB3\xE3\x83\x81\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB","\xE3\x83\x9F\xE3\x83\xAA\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB",
    "\xE3\x83\x87\xE3\x82\xB7\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB","\xE3\x83\x87\xE3\x82\xAB\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB","\xE3\x83\x98\xE3\x82\xAF\xE3\x83\x88\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB",
    "\xE3\x83\x9E\xE3\x82\xA4\xE3\x82\xAF\xE3\x83\xAD\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB","\xE3\x83\x8A\xE3\x83\x8E\xE3\x83\xA1\xE3\x83\xBC\xE3\x83\x88\xE3\x83\xAB","\xE3\x82\xBB\xE3\x83\xB3\xE3\x83\x81","\xE3\x83\x9F\xE3\x83\xAA",
    "\xE3\x83\x95\xE3\x82\xA3\xE3\x83\xBC\xE3\x83\x88","\xE3\x82\xA4\xE3\x83\xB3\xE3\x83\x81","\xE3\x83\x9E\xE3\x82\xA4\xE3\x83\xAB","\xE3\x83\xA4\xE3\x83\xBC\xE3\x83\x89",
    //(キロ-, センチ-, ミリ-, デシ-, デカ-, ヘクト-, マイクロ-, ナノメートル, センチ, ミリ, フィート, インチ, マイル, ヤード)
    "\xED\x82\xAC\xEB\xA1\x9C\xEB\xAF\xB8\xED\x84\xB0","\xEC\x84\xBC\xED\x8B\xB0\xEB\xAF\xB8\xED\x84\xB0","\xEB\xB0\x80\xEB\xA6\xAC\xEB\xAF\xB8\xED\x84\xB0","\xEB\x8D\xB0\xEC\x8B\x9C\xEB\xAF\xB8\xED\x84\xB0",
    "\xEB\x8D\xB0\xEC\xB9\xB4\xEB\xAF\xB8\xED\x84\xB0","\xED\x97\xA5\xED\x86\xA0\xEB\xAF\xB8\xED\x84\xB0","\xEB\xA7\x88\xEC\x9D\xB4\xED\x81\xAC\xEB\xA1\x9C\xEB\xAF\xB8\xED\x84\xB0","\xEB\x82\x98\xEB\x85\xB8\xEB\xAF\xB8\xED\x84\xB0",
    "\xED\x94\xBC\xED\x8A\xB8","\xEC\x9D\xB8\xEC\xB9\x98","\xEB\xA7\x88\xEC\x9D\xBC","\xEC\x95\xBC\xEB\x93\x9C"};
    //(킬로-, 센티-, 밀리-, 데시-, 데카-, 헥토-, 마이크로-, 나노미터, 피트, 인치, 마일, 야드)
  //EPSG unit codes of length other than the metre (PROJ's unit_of_measure table): another length only in a unit form
  //(EPSG:n or urn:ogc:def:uom:EPSG::n as the whole text), since several of them are also coordinate-system codes
  const int BU_EPSG_OTHER[]={1025,1033,9002,9003,9005,9014,9030,9031,9033,9034,9035,9036,9037,9038,9039,9040,9041,9042,9043,
    9050,9051,9052,9053,9060,9061,9062,9063,9070,9080,9081,9082,9083,9084,9085,9086,9087,9093,9094,9095,9096,9097,9098,9099,
    9300,9301,9302};
  const char *BU_FEET[]={"ft","feet","foot"};                            //(feet fused with a qualifier: ftus, footus, ...
  const char *BU_FEETQ[]={"us","s","m","i","int","intl","ind","uk","br"};// ... and usft, sft, usfeet, intfoot, ...)

  template<size_t N> bool BUIn(const char *(&tab)[N],const string &w);

  //UTF-8: the code point at s[i] and its length (a malformed byte is U+FFFD of length 1)
  uint32_t BUDecode(const string &s,const size_t i,size_t &len){
    const unsigned char c=(unsigned char)s[i]; len=1;
    if (c<0x80){return c;}
    int n=0; uint32_t cp=0;
    if      ((c&0xE0)==0xC0){n=1; cp=c&0x1F;}
    else if ((c&0xF0)==0xE0){n=2; cp=c&0x0F;}
    else if ((c&0xF8)==0xF0){n=3; cp=c&0x07;}
    else {return 0xFFFD;}
    for (int k=1;k<=n;k++){
      if (i+k>=s.size()){return 0xFFFD;}
      const unsigned char d=(unsigned char)s[i+k]; if ((d&0xC0)!=0x80){return 0xFFFD;}
      cp=(cp<<6)|(d&0x3F);
    }
    len=n+1; return cp;
  }
  void BUEncode(string &s,const uint32_t c){
    if      (c<0x80)   {s+=(char)c;}
    else if (c<0x800)  {s+=(char)(0xC0|(c>>6)); s+=(char)(0x80|(c&0x3F));}
    else if (c<0x10000){s+=(char)(0xE0|(c>>12)); s+=(char)(0x80|((c>>6)&0x3F)); s+=(char)(0x80|(c&0x3F));}
    else               {s+=(char)(0xF0|(c>>18)); s+=(char)(0x80|((c>>12)&0x3F)); s+=(char)(0x80|((c>>6)&0x3F)); s+=(char)(0x80|(c&0x3F));}
  }
  //lower case, and Greek without its accents (Greek capitals are often written without them: ΜΕΤΡΑ is μέτρα)
  uint32_t BULower(uint32_t c){
    if ((c>='A') && (c<='Z')){return c+32;}
    if ((c>=0xC0) && (c<=0xDE) && (c!=0xD7)){return c+32;}                          //(Latin-1)
    if ((c>=0x100) && (c<=0x137) && (c%2==0)){return c+1;}                          //(Latin Extended-A)
    if ((c>=0x139) && (c<=0x148) && (c%2==1)){return c+1;}
    if ((c>=0x14A) && (c<=0x177) && (c%2==0)){return c+1;}
    if ((c==0x179) || (c==0x17B) || (c==0x17D)){return c+1;}
    if ((c>=0x391) && (c<=0x3A9) && (c!=0x3A2)){c+=32;}                             //(Greek capitals, then the accents removed)
    else if (c==0x386){c=0x3B1;} else if (c==0x388){c=0x3B5;} else if (c==0x389){c=0x3B7;} else if (c==0x38A){c=0x3B9;}
    else if (c==0x38C){c=0x3BF;} else if (c==0x38E){c=0x3C5;} else if (c==0x38F){c=0x3C9;}
    if (c==0x3AC){return 0x3B1;} if (c==0x3AD){return 0x3B5;} if (c==0x3AE){return 0x3B7;} if ((c==0x3AF) || (c==0x3CA) || (c==0x390)){return 0x3B9;}
    if (c==0x3CC){return 0x3BF;} if ((c==0x3CD) || (c==0x3CB) || (c==0x3B0)){return 0x3C5;} if (c==0x3CE){return 0x3C9;}
    if ((c>=0x410) && (c<=0x42F)){return c+32;}                                     //(Cyrillic)
    if ((c>=0x400) && (c<=0x40F)){return c+80;}
    return c;
  }
  bool BULetter(const uint32_t c){
    if (c<0x80){return isalpha((int)c)!=0;}
    if (c==0xB5){return true;}                                                      //(the micro sign)
    if (c<0xC0){return false;}                                                      //(Latin-1 signs)
    if ((c==0xD7) || (c==0xF7)){return false;}                                      //(multiplication and division signs)
    if ((c>=0x2000) && (c<=0x2BFF)){return false;}                                  //(punctuation, symbols, arrows, ...)
    if ((c>=0x3000) && (c<=0x303F)){return false;}                                  //(CJK punctuation)
    if ((c>=0xFE30) && (c<=0xFE4F)){return false;}
    if ((c>=0xFF00) && (c<=0xFF20)){return false;}                                  //(full-width signs and digits)
    if ((c>=0xFF3B) && (c<=0xFF40)){return false;}
    if ((c>=0xFF5B) && (c<=0xFF65)){return false;}
    if ((c==0xFFFD) || (c>=0xFFF0 && c<=0xFFFF)){return false;}
    return true;
  }
  bool BUCJK(const uint32_t c){return (c>=0x2E80) && BULetter(c);}
  //the normalised text (above)
  string BUNormal(const string &t0,const bool brackets=true){
    string t;
    for (size_t i=0;i<t0.size();){
      size_t len=1; uint32_t c=BUDecode(t0,i,len); i+=len;
      if ((c>=0xFF01) && (c<=0xFF5E)){c-=0xFEE0;}                                   //(full-width ASCII)
      if ((c==0x3000) || (c==0xA0) || (c==0x2007) || (c==0x202F) || (c=='\t') || (c=='\r') || (c=='\n')){c=' ';}
      if ((c==' ') && (t.empty() || (t[t.size()-1]==' '))){continue;}
      BUEncode(t,BULower(c));
    }
    while (!t.empty() && (t[t.size()-1]==' ')){t.erase(t.size()-1);}
    while (brackets && (t.size()>=2)){                                              //(surrounding brackets)
      const char a=t[0],b=t[t.size()-1];
      if (!(((a=='(') && (b==')')) || ((a=='[') && (b==']')) || ((a=='{') && (b=='}')))){break;}
      t=t.substr(1,t.size()-2);
      while (!t.empty() && (t[0]==' ')){t.erase(0,1);}
      while (!t.empty() && (t[t.size()-1]==' ')){t.erase(t.size()-1);}
    }
    return t;
  }
  template<size_t N> bool BUIn(const char *(&tab)[N],const string &w){for (size_t i=0;i<N;i++){if (w==BUNormal(tab[i])){return true;}} return false;}

  bool BUMetre(const string &w){return BUIn(BU_METRE_WORD,w) || BUIn(BU_METRE_ABBR,w);}
  bool BUMetreForm(const string &w){
    if (BUMetre(w)){return true;}
    if ((w.size()>1) && (w.find('.')!=string::npos) && BUMetre((w[w.size()-1]=='.')?w.substr(0,w.size()-1):w+".")){return true;} //(a dotted abbreviation
                                                                                 // with or without its final dot: m.a.s.l, moh.)
    if ((w.size()>1) && (w[w.size()-1]=='.') && BUMetre(w.substr(0,w.size()-1))){return true;}
    for (size_t q=0;q<sizeof(BU_CJK_METRE)/sizeof(char*);q++){if (w==BU_CJK_METRE[q]){return true;}}
    return false;
  }
  bool BUDatum1(string r){
    for (size_t i=1;i+1<r.size();i++){                                           //(a blank before a datum's year: rh 2000, egm 96)
      if ((r[i]==' ') && isalpha((unsigned char)r[i-1]) && isdigit((unsigned char)r[i+1])){r.erase(i,1);}
    }
    if (BUIn(BU_DATUM,r)){return true;}
    if ((r.size()>1) && (r.find('.')!=string::npos)){                            //(a dotted abbreviation with or without its final dot)
      if (BUIn(BU_DATUM,(r[r.size()-1]=='.')?r.substr(0,r.size()-1):r+".")){return true;}
    }
    if ((r.size()>1) && (r[r.size()-1]=='.') && BUIn(BU_DATUM,r.substr(0,r.size()-1))){return true;}
    return false;
  }
  bool BUDatum(string r){
    { const size_t b=r.find_last_of(' ');                                       //(a geoid model after the datum: "navd88 geoid18",
      if ((b!=string::npos) && (b>0)){                                            // "navd88 (geoid12b)", "egm96 geoid")
        string g=r.substr(b+1);
        if ((g.size()>2) && (((g[0]=='(') && (g[g.size()-1]==')')) || ((g[0]=='[') && (g[g.size()-1]==']')))){g=g.substr(1,g.size()-2);}
        if (g.compare(0,5,"geoid")==0){r=r.substr(0,b);}
      }
    }
    return BUDatum1(r) || ((r.size()>6) && (r.compare(0,6,"above ")==0) && BUDatum1(r.substr(6))) || //("m above msl", "m over msl")
           ((r.size()>5) && (r.compare(0,5,"over ")==0) && BUDatum1(r.substr(5)));
  }
  //(a text in one pair of brackets: its inside)
  bool BUBracketed(const string &r,string &in){
    if ((r.size()<2) || !(((r[0]=='(') && (r[r.size()-1]==')')) || ((r[0]=='[') && (r[r.size()-1]==']')))){return false;}
    in=r.substr(1,r.size()-2);
    if (in.find_first_of("()[]")!=string::npos){return false;}
    while (!in.empty() && (in[0]==' ')){in.erase(0,1);}
    while (!in.empty() && (in[in.size()-1]==' ')){in.erase(in.size()-1);}
    return !in.empty();
  }
  //a metre form, alone or followed by a datum ("m msl", "m, msl", "m (msl)", fused "mahd", "m.s.n.m.") or by a metre form or
  //EPSG 9001 in brackets ("meter (m)", "metre (EPSG:9001)")
  bool BUMetreDatum(const string &t){
    if (BUMetreForm(t)){return true;}
    if ((t.size()>1) && (t[0]=='m')){
      const string r=t.substr(1);
      if ((r.find(' ')==string::npos) && BUDatum(r)){return true;}
      if ((r[0]=='.') && (r.size()>1) && BUDatum(r.substr(1))){return true;}
    }
    for (size_t p=1;p<t.size();p++){                                            //(after a blank, a comma, ; / - _ or a bracket:
      if (strchr(" ,;/-_([",t[p])==NULL){continue;}                               // "m, navd88", "m-navd88", "m(navd88)")
      if (!BUMetreForm(t.substr(0,p))){continue;}
      string r=t.substr(p); while (!r.empty() && (strchr(" ,;/-_",r[0])!=NULL)){r.erase(0,1);}
      if (r.empty()){continue;}
      string in;
      if (BUDatum(r) || (BUBracketed(r,in) && (BUDatum(in) || BUMetreDatum(in)))){return true;}
    }
    return false;
  }
  //the whole text is metres: a metre form with a datum (above), or a quantity with one in brackets or after "in"
  bool BUMetreText(const string &t){
    if (BUMetreDatum(t)){return true;}
    for (size_t q=0;q<sizeof(BU_QUANT)/sizeof(char*);q++){
      const string qn=BUNormal(BU_QUANT[q]); const size_t n=qn.size();
      if ((t.size()<=n) || (t.compare(0,n,qn)!=0) || ((t[n]!=' ') && (t[n]!='(') && (t[n]!='[') && (t[n]!=',') && (t[n]!='_'))){continue;}
      if ((t[n]=='_') && BUMetreDatum(t.substr(n+1))){return true;}             //("elevation_m")
      string r=t.substr(n); while (!r.empty() && (r[0]==' ')){r.erase(0,1);}
      if (!r.empty() && (r[0]==',')){r.erase(0,1); while (!r.empty() && (r[0]==' ')){r.erase(0,1);} if (BUMetreDatum(r)){return true;} continue;} //("Elevation, m")
      string in;
      if ((r.compare(0,3,"in ")==0) && BUMetreDatum(r.substr(3))){return true;}
      if (BUBracketed(r,in) && BUMetreDatum(in)){return true;}
    }
    return false;
  }
  //(a token that is a number: 0.1, 0,01, 1e-3, 1/1000, 1:1000, x0.1, 1000s)
  bool BUNumberToken(string k){
    if (!k.empty() && ((k[0]=='x') || (k[0]=='+') || (k[0]=='-'))){k.erase(0,1);}
    if (!k.empty() && (k[k.size()-1]=='s')){k.erase(k.size()-1);}
    if (k.empty() || !(isdigit((unsigned char)k[0]) || ((k[0]=='.') && (k.size()>1) && isdigit((unsigned char)k[1])))){return false;} //(.1 too)
    for (size_t i=0;i<k.size();i++){
      const char c=k[i];
      if (isdigit((unsigned char)c) || (c=='.') || (c==',') || (c==':')){continue;}
      if ((c=='e') && (i>0) && (i+1<k.size()) && (isdigit((unsigned char)k[i+1]) || (k[i+1]=='-') || (k[i+1]=='+'))){continue;}
      if (((c=='-') || (c=='+')) && (i>0) && (k[i-1]=='e')){continue;}
      return false;
    }
    return true;
  }
  bool BUFusedFeet(const string &w){
    for (int i=0;i<3;i++){const size_t n=strlen(BU_FEET[i]); if ((w.size()>n) && (w.compare(0,n,BU_FEET[i])==0) && BUDatum(w.substr(n))){return true;}} //(ftmsl, ftasl, feetmsl)
    for (int i=0;i<3;i++){
      const size_t n=strlen(BU_FEET[i]);
      if ((w.size()>n) && (w.compare(0,n,BU_FEET[i])==0) && BUIn(BU_FEETQ,w.substr(n))){return true;}  //(ftus, footus, ...)
      if ((w.size()>n) && (w.compare(w.size()-n,n,BU_FEET[i])==0)){
        if (i>0){return true;}                                                      //(...feet, ...foot)
        if (BUIn(BU_FEETQ,w.substr(0,w.size()-n))){return true;}                    //(usft, sft, ...)
      }
    }
    return false;
  }
  //feet fused with a datum, its digits kept ("ftnavd88", "ftnavd 88", "usftngvd29", "ftnavd88 geoid"): read on a whole
  //text or on a word between blanks, commas and brackets, since words of letters alone part the datum from its year
  bool BUFeetDatum(const string &w){
    const char *f[6]={"usft","usfeet","usfoot","ft","feet","foot"};
    for (int i=0;i<6;i++){const size_t n=strlen(f[i]); if ((w.size()>n) && (w.compare(0,n,f[i])==0) && isalpha((unsigned char)w[n]) && BUDatum(w.substr(n))){return true;}}
    return false;
  }
  bool BUPrefixed(const string &w){                                                 //(a prefix joined to a metre word)
    for (size_t i=0;i<sizeof(BU_PREFIX)/sizeof(BU_PREFIX[0]);i++){
      const string p=BUNormal(BU_PREFIX[i]); const size_t n=p.size();
      if ((w.size()>n) && (w.compare(0,n,p)==0) && BUIn(BU_METRE_WORD,w.substr(n))){return true;}
    }
    return false;
  }
}

//the class of a band unit text (see above): BU_METRE, BU_NONE (no unit), BU_LENGTH (another length unit and no metre,
//for the message) or BU_TEXT (any other text)
int GWGeoTIFF::BandUnit(const string &ut0)
{
  const string t=BUNormal(ut0);
  if (t.empty()){return BUNormal(ut0,false).empty()?BU_NONE:BU_TEXT;}             //(empty brackets, "()", are a text)
  if (BUMetreText(t)){return BU_METRE;}
  //-- another length unit, for the message: an EPSG length code in a unit form, or words of another length and none of the metre
  { size_t r=string::npos;
    const string urn="urn:ogc:def:uom:epsg::",uri="http://www.opengis.net/def/uom/epsg/0/",uris="https://www.opengis.net/def/uom/epsg/0/";
    if ((t.compare(0,4,"epsg")==0) && (t.size()>4) && ((t[4]==':') || (t[4]==' '))){r=4; while ((r<t.size()) && ((t[r]==':') || (t[r]==' '))){r++;}} //(epsg:n, epsg n, epsg: n)
    else if (t.compare(0,urn.size(),urn)==0){r=urn.size();}
    else if (t.compare(0,uri.size(),uri)==0){r=uri.size();} else if (t.compare(0,uris.size(),uris)==0){r=uris.size();}
    if ((r!=string::npos) && (r<t.size()) && (t.find_first_not_of(":0123456789",r)==string::npos)){
      const long code=atol(t.substr(t.find_first_not_of(':',r)==string::npos?r:t.find_first_not_of(':',r)).c_str());
      for (size_t q=0;q<sizeof(BU_EPSG_OTHER)/sizeof(int);q++){if (code==BU_EPSG_OTHER[q]){return BU_LENGTH;}}
    }
  }
  if ((t.compare(0,5,"link_")==0) || (t.compare(0,6,"links_")==0) || (t.compare(0,6,"chain_")==0) || (t.compare(0,4,"rod_")==0)){return BU_LENGTH;} //(ESRI
                                                                                    // unit names: Link_Clarke, Chain_Benoit_1895_B, Rod_US)
  bool ftd=BUFeetDatum(t);
  { string tok; for (size_t i=0;i<=t.size();i++){if ((i==t.size()) || (strchr(" ,;()[]{}",t[i])!=NULL)){if (BUFeetDatum(tok)){ftd=true;} tok="";} else {tok+=t[i];}} }
  vector<string> wd; vector<bool> cj; string cur; bool curCJK=false,gapBlank=false;  //(words: runs of letters; Chinese, Japanese and
  for (size_t i=0;i<=t.size();){                                                    // Korean runs apart, joined across blanks: 千 米, キロ メートル)
    size_t len=1; const uint32_t c=(i<t.size())?BUDecode(t,i,len):(uint32_t)' ';
    const bool let=BULetter(c),cjk=BUCJK(c);
    if (!cur.empty() && (!let || (cjk!=curCJK))){
      if (curCJK && gapBlank && !wd.empty() && cj.back()){wd.back()+=cur;} else {wd.push_back(cur); cj.push_back(curCJK);}
      cur=""; gapBlank=true;
    }
    if (let){if (cur.empty()){curCJK=cjk;} cur.append(t,i,len);}
    else if (c!=' '){gapBlank=false;}
    i+=len;
  }
  bool other=false,metre=false,part=false;
  for (size_t j=0;j<wd.size();j++){
    const string &w=wd[j];
    if (cj[j]){                                                                     //(the longest entry at each place)
      for (size_t p=0;p<w.size();){
        size_t best=0; int bc=0;
        for (int k=0;k<3;k++){
          const char **tab=(k==0)?BU_CJK_OTHER:((k==1)?BU_CJK_METRE:BU_CJK_NONE);
          const size_t n=(k==0)?sizeof(BU_CJK_OTHER)/sizeof(char*):((k==1)?sizeof(BU_CJK_METRE)/sizeof(char*):sizeof(BU_CJK_NONE)/sizeof(char*));
          for (size_t q=0;q<n;q++){const size_t L=strlen(tab[q]); if ((L>best) && (w.compare(p,L,tab[q])==0)){best=L; bc=k;}}
        }
        if (best==0){size_t l2=1; BUDecode(w,p,l2); p+=l2; continue;}
        if (bc==0){other=true;} else if (bc==1){metre=true;}
        p+=best;
      }
      continue;
    }
    if (BUMetre(w)){metre=true; continue;}
    if (BUIn(BU_OTHER,w) || BUFusedFeet(w) || BUPrefixed(w)){other=true; continue;}
    if (BUIn(BU_PREFIX,w) && (j+1<wd.size()) && BUIn(BU_METRE_WORD,wd[j+1])){other=true; j++; continue;} //("Kilo Meter": one unit)
    if (BUIn(BU_AMBIG,w) && ((wd.size()==1) || ((wd.size()==2) && (j==1) && BUIn(BU_QUAL,wd[0])))){other=true;}
    if (BUIn(BU_PART,w)){part=true;}
  }
  if (ftd){other=true;}
  if (other && !metre){return BU_LENGTH;}
  if (other && metre && !t.empty() && ((t[t.size()-1]==')') || (t[t.size()-1]==']'))){ //-- another length named with its metre conversion:
    const size_t ob=t.find_last_of("([");                                          //   "ft (0.3048 m)" or a metre form and the length, "m (feet)"
    if ((ob!=string::npos) && (ob>0)){
      string head=t.substr(0,ob),in=t.substr(ob+1,t.size()-ob-2);
      while (!head.empty() && (head[head.size()-1]==' ')){head.erase(head.size()-1);}
      while (!in.empty() && (in[0]==' ')){in.erase(0,1);} while (!in.empty() && (in[in.size()-1]==' ')){in.erase(in.size()-1);}
      size_t i=0; while ((i<in.size()) && (isdigit((unsigned char)in[i]) || (in[i]=='.') || (in[i]==','))){i++;}
      string rest=in.substr(i); while (!rest.empty() && (rest[0]==' ')){rest.erase(0,1);}
      const bool conv=(i>0) && BUNumberToken(in.substr(0,i)) && BUMetreForm(rest);
      if ((in.find_first_of("()[]")==string::npos) && ((conv && (GWGeoTIFF::BandUnit(head)==BU_LENGTH)) || (BUMetreForm(head) && (GWGeoTIFF::BandUnit(in)==BU_LENGTH)))){return BU_LENGTH;}
    }
  }
  if (metre && !other){ //-- a scaled metre: a part or multiple ("m (in thousands)", "tenths of a metre"), a number with a factor
                        //   sign or word ("m*0.1", "Meter (x0.1)", "m (scaled by 0.1)", "0.1 of a metre", "m^2"), or a number and a
                        //   metre form as the whole text ("0.1 m", "1000 m"); other numbers (datum years, resolutions) are not scales
    if (part){return BU_SCALED;}
    string s=t; size_t q=0;
    while ((q=s.find("+/-"))!=string::npos){s.replace(q,3," ");}
    const char *sign[]={"*","/","^","\xC3\x97","\xC3\xB7","\xC2\xB7","\xE2\x8B\x85"};   //(*, /, ^, the signs of multiplication and division, dots)
    bool hasSign=false;
    for (int k=0;k<7;k++){q=0; while ((q=s.find(sign[k],q))!=string::npos){hasSign=true; s.replace(q,strlen(sign[k])," ");}}
    //(punctuation parts the text into groups, "|": a unit word or an ordinal counts only next to a number of its own group)
    { string g; for (q=0;q<s.size();q++){if (strchr(";()[]{}=,:",s[q])!=NULL){g+=" | ";} else {g+=s[q];}} s=g; }
    istringstream is(s); string k; bool num=false,word=false,prevNum=false; vector<string> tk;
    while (is>>k){tk.push_back(k);}
    for (size_t i=0;i<tk.size();i++){
      k=tk[i]; if (k=="|"){continue;}
      const string pv=(i>0)?tk[i-1]:string("|"),nx=(i+1<tk.size())?tk[i+1]:string("|"),nx2=(i+2<tk.size())?tk[i+2]:string("|");
      const string pv2=(i>1)?tk[i-2]:string("|");
      const bool n=BUNumberToken(k); if (n || isdigit((unsigned char)k[0])){num=true;}       //(also with superscripts: 10⁻²)
      if (n && (k[k.size()-1]=='s')){word=true;}                                  //("in 1000s": a multiple)
      if (n && (k.find('e')!=string::npos)){word=true;}                           //(a power of ten: "m 1e-2")
      if ((k.size()>1) && (k[0]=='e') && BUNumberToken("1"+k)){num=true; word=true;} //("m e-2")
      if (isdigit((unsigned char)k[0]) && (k.size()>2)){                            //(an ordinal as a part: "100ths of a metre", "in 100ths";
        const bool th=(k.compare(k.size()-2,2,"th")==0),ths=(k.size()>3) && (k.compare(k.size()-3,3,"ths")==0); // not "19th century")
        if ((th || ths) && ((nx=="of") || (ths && (pv=="in")))){word=true;}
      }
      if ((k.size()>1) && (k[0]=='x') && isdigit((unsigned char)k[1])){hasSign=true;}
      if ((k=="x") || (k=="scale") || (k=="scaled") || (k=="times") || (k=="multiplier") || (k=="factor") || (k=="divide") || (k=="divided") ||
          (k=="multiply") || (k=="multiplied") || (k=="per") || ((k=="of") && prevNum) || (k=="faktor") || (k=="facteur") || (k=="fattore") ||
          (k=="fator")){word=true;}
      if ((k=="unit") || (k=="units")){ //(next to a number of its group: "units of 0.1 m", "0.1 m units", "0.1 units"; not "NAVD 1988 (units)")
        if (BUNumberToken(pv) || BUNumberToken(nx) || ((nx=="of") && BUNumberToken(nx2)) || (BUMetreForm(pv) && BUNumberToken(pv2))){word=true;}
      }
      prevNum=n;
    }
    if (num && (hasSign || word)){return BU_SCALED;}
    { const size_t b=t.find_first_of(" "); //(a number, then a metre form, as the whole text)
      string nt=t,rest; while (!nt.empty() && ((nt[0]=='(') || (nt[0]=='['))){nt.erase(0,1);} while (!nt.empty() && ((nt[nt.size()-1]==')') || (nt[nt.size()-1]==']'))){nt.erase(nt.size()-1);}
      (void)b; size_t i=0; while ((i<nt.size()) && (isdigit((unsigned char)nt[i]) || strchr(".,/:e",nt[i]) || (((nt[i]=='-') || (nt[i]=='+')) && (i>0) && (nt[i-1]=='e')))){i++;} //(1e-2 m)
      rest=nt.substr(i); while (!rest.empty() && (rest[0]==' ')){rest.erase(0,1);}
      if ((i>0) && BUNumberToken(nt.substr(0,i)) && BUMetreForm(rest)){return BU_SCALED;}
    }
  }
  return BU_TEXT;
}

//GDAL's side file <file>.aux.xml, as GDAL looks for it (the raster's name with .aux.xml added, in that case): has is set
//when it exists
bool GWGeoTIFF::ReadSideFile(const string &filename,string &pam,bool &has,string &err)
{
  pam=""; has=false;
  ifstream X((filename+".aux.xml").c_str(),ios::binary);
  if (!X){return true;}
  X.seekg(0,ios::end); const streamoff xs=X.tellg(); X.seekg(0,ios::beg);
  if ((xs<0) || ((unsigned long long)xs>(16ULL<<20))){err="raster "+filename+": its GDAL side file "+filename+".aux.xml is larger than the 16 MB Raven reads, so its band scale, offset and unit cannot be checked"; return false;}
  ostringstream s; s<<X.rdbuf(); pam=s.str(); has=true;
  return true;
}

//the band's scale, offset, unit and side-file no-data value as GDAL reads them, and the band-unit rule (see GWGeoTIFF.h)
bool GWGeoTIFF::BandMeta(const string &F,const string &meta,const string *pam,const string &side,const bool metres,const string &cmd,
                         double &bscale,double &boffset,bool &pamHasNodata,string &pamNodata,string &note,string &err,const bool georef)
{
  bscale=1.0; boffset=0.0; pamHasNodata=false; pamNodata=""; note="";
  const string rewrite="rewrite it with GDAL: copy the raster with gdal_translate, or set the values with gdal_edit.py (-scale, -offset, -units, -a_nodata), "
                       "which write them in GDAL's form";
  map<string,string> items; string ge;
  if ((meta!="") && (GDALBandItems(meta,items,ge)==2)){
    err=F+"its GDAL_METADATA tag is not in the form GDAL writes ("+ge+") and names a band scale, offset or unit, which GDAL may read differently: "+rewrite; return false;
  }
  for (int k=0;k<2;k++){
    const char *nm=(k==0)?"SCALE":"OFFSET"; map<string,string>::const_iterator it=items.find(nm); if (it==items.end()){continue;}
    char *endp; const double v=strtod(it->second.c_str(),&endp);
    if (it->second.empty() || (*endp!='\0') || !std::isfinite(v)){err=F+"its band "+string(nm)+" in GDAL_METADATA ("+it->second+") is not a number"; return false;}
    if (k==0){bscale=v;} else {boffset=v;}
  }
  string bsrc="in GDAL_METADATA";
  map<string,string>::const_iterator iu=items.find("UNITTYPE");
  string unit=(iu!=items.end())?iu->second:string(""),usrc="UNITTYPE in GDAL_METADATA";
  if (pam!=NULL){ //GDAL's side file, which GDAL reads before the file's own metadata: a Scale or Offset there replaces both of
                  //the file's (the other is then 1 or 0); a UnitType replaces the file's; a NoDataValue replaces GDAL_NODATA
    map<string,string> pi; string pe;
    if (georef){ //georeferencing in the side file (GeoTransform, SRS, GCPs): GDAL may take it in place of a GeoTIFF's own, in any form
      string l=*pam; for (size_t i=0;i<l.size();i++){l[i]=(char)tolower((unsigned char)l[i]);}
      const char *g[4]={"geotransform","<srs","gcplist","<gcp"}; string found;
      for (int k=0;k<4;k++){if (l.find(g[k])!=string::npos){found=g[k];}}
      for (size_t i=0;found.empty() && ((i=l.find('<',i))!=string::npos);i++){ //("< SRS", a blank after the bracket)
        size_t j=i+1; while ((j<l.size()) && isspace((unsigned char)l[j])){j++;}
        if ((l.compare(j,3,"srs")==0) || (l.compare(j,3,"gcp")==0)){found="srs";}
      }
      if (!found.empty()){
        err=F+"its GDAL side file "+side+" holds georeferencing (GeoTransform, SRS or GCPs), which GDAL may take in place of the raster's own: rewrite it with "
            "gdal_translate (gdal_translate in.tif out.tif writes the georeferencing into the raster), or remove the side file"; return false;
      }
    }
    if (PAMBandItems(*pam,pi,pe)==2){
      err=F+"its GDAL side file "+side+" is not in the form GDAL writes ("+pe+") and names a band or a band value (Scale, Offset, UnitType, NoDataValue), "
          "which GDAL may read differently: "+rewrite+", or remove the side file"; return false;
    }
    if (pi.count("SCALE") || pi.count("OFFSET")){
      bscale=1.0; boffset=0.0; bsrc="in "+side;
      for (int k=0;k<2;k++){
        const char *nm=(k==0)?"SCALE":"OFFSET"; map<string,string>::const_iterator it=pi.find(nm); if (it==pi.end()){continue;}
        char *endp; const double v=strtod(it->second.c_str(),&endp);
        if (it->second.empty() || (*endp!='\0') || !std::isfinite(v)){err=F+"its band "+string((k==0)?"Scale":"Offset")+" in "+side+" ("+it->second+") is not a number"; return false;}
        if (k==0){bscale=v;} else {boffset=v;}
      }
    }
    if (pi.count("UNITTYPE")){unit=pi["UNITTYPE"]; usrc="UnitType in "+side;}
    if (pi.count("NODATA")){pamNodata=pi["NODATA"]; pamHasNodata=true;}
  }
  if (bscale==0.0){err=F+"its band scale "+bsrc+" is 0 (every value would be the offset): save the raster with its elevations (e.g. gdal_translate -unscale)"; return false;}
  //the band unit (BandUnit): read as metres when the whole text is a metre form, or when there is no unit; another length and
  //a scaled metre stop, also with METRES; any other text stops unless METRES follows the file name: the user's statement that
  //the values are metres, after which the file's own unit is reported every run (also for a metre form or no unit)
  const int bu=BandUnit(unit);
  const string shown=Shown(unit),hint=cmd.empty()?string(":LandSurface or :BedrockSurface <file>"):cmd;
  if (bu==BU_LENGTH){
    err=F+"its elevations are in "+shown+" (the band unit, "+usrc+"); Raven reads elevations in metres: convert the values to metres and set the band unit to m "
        "(e.g. gdal_calc.py, then gdal_edit.py -units m)"+(metres?string("; METRES after the file name does not apply to a unit of another length"):string("")); return false;
  }
  if (bu==BU_SCALED){
    err=F+"its band unit is '"+shown+"' ("+usrc+"), a metre with a number or a part (a scale, factor or exponent), so the values are not plain metres: "
        "write the elevations unscaled in metres, or give the scale as the band's scale and offset (e.g. gdal_edit.py -scale 0.1 -units m)"+
        (metres?string("; METRES after the file name does not apply to a scaled unit"):string("")); return false;
  }
  if (metres){
    note=F+"read as metres because METRES follows the file name ("+hint+" METRES); "+((bu==BU_NONE)?string("the file gives no band unit"):
         "the file's own band unit is '"+shown+"' ("+usrc+")")+": check that its values are metres";
    return true;
  }
  if (bu==BU_TEXT){
    err=F+"its band unit is '"+shown+"' ("+usrc+"), which is not one of the metre forms Raven reads as metres (m, metre, metres, meter, meters, m MSL, EPSG:9001 and the "
        "others in the manual's table of band units): if its elevations are in metres, state it by writing METRES after the file name ("+hint+" METRES); otherwise "
        "convert the values to metres and set the band unit to m (e.g. gdal_calc.py, then gdal_edit.py -units m)"; return false;
  }
  return true;
}

bool GWGeoTIFF::IsTIFFName(const string &filename)
{
  size_t p=filename.find_last_of('.'); if (p==string::npos){return false;}
  string e=filename.substr(p+1); for (size_t i=0;i<e.size();i++){e[i]=(char)tolower((unsigned char)e[i]);}
  size_t sl=filename.find_last_of("/\\"); if ((sl!=string::npos) && (sl>p)){return false;}
  return (e=="tif") || (e=="tiff") || (e=="gtiff");
}

bool GWGeoTIFF::Read(const string &filename,gw_geotiff &R,string &err,const bool metres,const string &cmd)
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
  string pam; bool hasPam=false; //(GDAL's side file, read before the file's own metadata as GDAL reads it)
  if (!ReadSideFile(filename,pam,hasPam,err)){return false;}
  return ReadMemory(data.empty()?NULL:&data[0],data.size(),filename,R,err,metres,hasPam?&pam:NULL,cmd);
}

bool GWGeoTIFF::ReadMemory(const unsigned char *data,const size_t n,const string &name,gw_geotiff &R,string &err,const bool metres,const string *pam,const string &cmd)
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
  //-- GDAL's band scale, offset and unit (GDAL_METADATA, and the side file <raster>.aux.xml, which GDAL reads first): the
  //   values are read unscaled, as GDAL's -unscale and GIS programs show them; the unit must be a metre form (BandMeta)
  double bscale=1.0,boffset=0.0; string pamNodata; bool pamHasNodata=false;
  { const uint64_t MAXMETA=16ULL<<20;
    const TEntry *tm=T.Tag(42112);
    if ((tm!=NULL) && (tm->type==2) && (tm->count>MAXMETA)){
      err=F+"its GDAL_METADATA tag holds "+to_string((unsigned long long)(tm->count>>20))+" MB, more than the 16 MB Raven reads, so its band scale, offset and unit cannot be checked: "
          "save the raster without the extra metadata (e.g. gdal_translate -unscale -co PROFILE=GeoTIFF in.tif out.tif)"; return false;
    }
    string note;
    if (!BandMeta(F,T.Ascii(42112,MAXMETA),pam,name+".aux.xml",metres,cmd,bscale,boffset,pamHasNodata,pamNodata,note,err)){return false;}
    if (note!=""){R.notes+=string(R.notes.empty()?"":"\n")+note;}
  }
  if (R.rasterType==2){xul-=0.5*sx; yul+=0.5*sy;}    //PixelIsPoint: the tie point is the centre of the pixel
  //-- NODATA
  string nd=pamHasNodata?pamNodata:T.Ascii(42113); //(the side file's NoDataValue first, as GDAL reads it)
  { size_t a=nd.find_first_not_of(" \t\r\n"),b=nd.find_last_not_of(" \t\r\n");
    nd=(a==string::npos)?"":nd.substr(a,b-a+1); }
  if (nd!=""){
    string l=nd; for (size_t i=0;i<l.size();i++){l[i]=(char)tolower((unsigned char)l[i]);}
    if ((l=="nan") || (l=="-nan") || (l=="+nan")){R.hasNodata=true; R.nodata=NAN;}
    else {
      char *endp; double v=strtod(nd.c_str(),&endp);
      if ((endp==nd.c_str()) || (*endp!='\0')){err=F+(pamHasNodata?"its NoDataValue in "+name+".aux.xml":string("its GDAL_NODATA tag"))+" ("+nd+") is not a number"; return false;}
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
