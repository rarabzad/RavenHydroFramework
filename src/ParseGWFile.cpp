/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Parsing of the groundwater (.rvg) file for the Raven-MODFLOW 6 coupling
----------------------------------------------------------------*/
#include <cerrno>
#include <climits>
#include "RavenInclude.h"
#include "GWGrid.h"
#include "GWUtil.h"
#include "Model.h"
#include "ParseLib.h"
#include "GroundwaterModel.h"

//////////////////////////////////////////////////////////////////
/// \brief Parses the groundwater file, model.rvg
/// \details Geometry, grid and coupling settings. Hydrogeology (aquifer classes and profiles)
///  lives in the .rvp file; which HRUs have groundwater is set by AQUIFER_PROFILE in the .rvh file.
///
///  :HRUGeometry      [file.geojson] {ID field, default HRU_ID}   (REQUIRED)
///  :LandSurface      [ESRI ASCII raster, lon/lat]               (optional; else HRU ELEVATION)
///  :BedrockSurface   [ESRI ASCII raster, lon/lat]               (optional; else DEFAULT_BEDROCK_DEPTH)
///  :InputCRS         [projection of the geometry and raster files] (optional; else lon/lat)
///  :GridCellSize     [m]                                        (default 1000)
///  :MinCellCoverage  [fraction]                                 (default 0.25)
///  :GridType REGULAR|QUADTREE|HRU_MESH|FILE, :GridRotation, :GridRefinement, :GridFile, :QuadtreeMaxLevel,
///  :FlowCorrection, :LayerConnection  (grid options, see the manual)
///  :MF6Library       [path to libmf6]                           (optional; else RAVEN_MF6_LIB)
///  :SeepageReturnTo  [state variable]                           (default: lowest soil layer)
///  :HeadSaveFrequency [time steps]                              (default 30)
///  :WriteGWBudgetFile                                           (optional; MODFLOW budget file and porosity)
///
/// \param *&pModel [in] The input model object
/// \param &Options [in] Global model options information
/// \return true if parsing succeeded
//
//reads an :Attributes table block until the given end tag; rows are {ID, values...}
static void ReadAttributeTable(CParser *p,const char *endtag,vector<string> &attrs,vector<vector<string> > &rows)
{
  char *s[MAXINPUTITEMS]; int Len;
  while (!p->Tokenize(s,Len)){
    if (IsComment(s[0],Len)){continue;}
    if (!strcmp(s[0],endtag)){break;}
    if (!strcmp(s[0],":Units")){continue;}
    if (!strcmp(s[0],":Attributes")){attrs.clear(); for (int i=1;i<Len;i++){attrs.push_back(s[i]);} continue;}
    ExitGracefullyIf(Len!=(int)attrs.size()+1,(string("ParseGWFile: wrong number of values in a row of the :")+(endtag+4)+" table (one ID and one value per :Attributes entry)").c_str(),BAD_DATA);
    vector<string> r; for (int i=0;i<Len;i++){r.push_back(s[i]);} rows.push_back(r);
  }
}
static string GetAttr(const vector<string> &attrs,const vector<string> &row,const string &name,const string &def)
{
  for (size_t i=0;i<attrs.size();i++){if (attrs[i]==name){return row[i+1];}}
  ExitGracefullyIf(def=="REQUIRED",("ParseGWFile: missing attribute "+name).c_str(),BAD_DATA);
  return def;
}

//a malformed groundwater command changes the physics if ignored: stop with a clear message
static void GWImproperFormat(char **s,const optStruct &Options)
{
  ExitGracefully(("ParseGWFile: improper format of command "+string(s[0])+" in "+Options.rvg_filename).c_str(),BAD_DATA);
}
//numbers are read strictly: the whole token must be a finite number (no units, no thousands separators; a comma
//separates tokens, so 2,000.0 is read as the two tokens 2 and 000.0)
static double GWNumber(const string &cmd,const char *tok,const optStruct &Options)
{
  double v=0.0;
  if (!gwutil::StrictNumber(tok,v)){
    ExitGracefully(("ParseGWFile: command "+cmd+" in "+Options.rvg_filename+": '"+string(tok)+"' is not a number (give a plain number "
                    "in the units of the manual, without units or thousands separators)").c_str(),BAD_DATA);
  }
  return v;
}
static int GWInteger(const string &cmd,const char *tok,const optStruct &Options)
{
  double v=GWNumber(cmd,tok,Options);
  if ((v!=floor(v)) || (fabs(v)>2.0e9)){
    ExitGracefully(("ParseGWFile: command "+cmd+" in "+Options.rvg_filename+": '"+string(tok)+"' is not a whole number").c_str(),BAD_DATA);
  }
  return (int)v;
}
//a whole number that names something (role: "HRU ID", "cell number", "well ID", ...), within the long long range
static long long GWID(const string &cmd,const char *tok,const optStruct &Options,const string &role="ID")
{
  char *e=NULL; errno=0; long long v=strtoll(tok,&e,10);
  if (errno==ERANGE){
    ExitGracefully(("ParseGWFile: command "+cmd+" in "+Options.rvg_filename+": the "+role+" '"+string(tok)+"' is out of range (beyond "+(v<0?"-9223372036854775808":"9223372036854775807")+")").c_str(),BAD_DATA);
  }
  if ((e==tok) || (*e!='\0')){
    ExitGracefully(("ParseGWFile: command "+cmd+" in "+Options.rvg_filename+": the "+role+" '"+string(tok)+"' is not a whole number").c_str(),BAD_DATA);
  }
  return v;
}
static double GWNum(char **s,const int i,const optStruct &Options){return GWNumber(s[0],s[i],Options);}
static int    GWInt(char **s,const int i,const optStruct &Options){return GWInteger(s[0],s[i],Options);}
//a command with more tokens than it takes stops (a thousands separator splits a number into two tokens)
static void GWMaxTokens(char **s,const int Len,const int nmax,const optStruct &Options)
{
  if (Len>nmax){
    string glued; //(a # inside a value is part of it: a comment starts at a # after a blank)
    for (int i=1;i<nmax;i++){if (strchr(s[i],'#')!=NULL){glued=string(s[i]); break;}}
    ExitGracefully(("ParseGWFile: command "+string(s[0])+" in "+Options.rvg_filename+" has the surplus value '"+string(s[nmax])+"' (it takes "+
                    to_string(nmax-1)+" value"+((nmax==2)?"":"s")+" at most; "+(glued.empty()?string("a comma separates values, so write 2000.0, not 2,000.0)"):
                    "a # glued to a value is part of it, as in '"+glued+"': write a blank before # to start a comment)")).c_str(),BAD_DATA);
  }
}
//the most tokens a command takes (with the command itself; 0: a list of any length, or a form the command's own reader
//refuses): one check for every command, so that no value after the last one a command takes is dropped without a message
static int GWCommandArity(char **s,const int Len)
{
  const string c=s[0],a=(Len>=2)?string(s[1]):string("");
  if ((c==":StreamPackage") && (Len>=3)){return !strcmp(s[2],"AUX")?4:(!strcmp(s[2],"DOMINANT_HRU")?3:0);}
  if (c==":GridRefinement"){return ((a=="RIVERS") || (a=="WELLS") || (a=="LAKES"))?3:(((a=="LINES") || (a=="POLYGONS"))?4:0);}
  if (c==":SoilZone"){return (a=="CONCEPTUAL")?2:0;}
  const char *one[]={":OverlapWeights",":Wells",":ObservationWells",":ReservoirExchange",":WriteNetCDFHeads",":WriteGWBudgetFile"};
  const char *two[]={":GridCellSize",":MF6Simulation",":MF6Model",":StreamPackage",":SeepageToRaven",":GridType",":GridRotation",":QuadtreeMaxLevel",
    ":FlowCorrection",":LayerConnection",":MinCellCoverage",":MF6Library",":SoilZoneSubgrid",":SoilZoneVadoseClasses",":SoilZoneVadoseInitialFlux",
    ":SoilZoneCapillaryRise",":SoilZoneExchangeTiming",":SoilZoneClasses",":SoilZoneAquiferClass",":SoilZoneOutput",":SeepageReturnTo",":RiverBedDepth",
    ":LinkageCache",":SolverHeadTolerance",":SolverMaxOuterIterations",":HeadSaveFrequency"};
  const char *three[]={":HRUGeometry",":RiverGeometry",":LandSurface",":BedrockSurface",":MF6GridFile",":MF6StartDate",":GridFile",":GWInitialization"};
  const char *four[]={":SoilZoneIteration",":SoilZoneVadose"};
  for (size_t i=0;i<sizeof(one)/sizeof(one[0]);i++){if (c==one[i]){return 1;}}
  for (size_t i=0;i<sizeof(two)/sizeof(two[0]);i++){if (c==two[i]){return 2;}}
  for (size_t i=0;i<sizeof(three)/sizeof(three[0]);i++){if (c==three[i]){return 3;}}
  for (size_t i=0;i<sizeof(four)/sizeof(four[0]);i++){if (c==four[i]){return 4;}}
  if (c==":SurfaceWaterExchange"){return 5;}
  return 0;
}
//////////////////////////////////////////////////////////////////
/// \brief true if the .rvg file asks for :SoilZone PHYSICAL. Read ahead of the .rvp file: Raven discards soil
///  parameters that no process needs, and the soil zone needs HYDRAUL_COND, POROSITY and FIELD_CAPACITY
//
bool GWFileUsesPhysicalSoilZone(const optStruct &Options)
{
  if ((Options.modeltype!=MODELTYPE_COUPLED) || (Options.rvg_filename=="")){return false;}
  ifstream RVG(Options.rvg_filename.c_str());
  if (RVG.fail()){return false;}
  CParser *p=new CParser(RVG,Options.rvg_filename,0);
  char *s[MAXINPUTITEMS]; int Len; bool found=false;
  while (!p->Tokenize(s,Len)){                                    //as ParseGWFile: the last :SoilZone counts, :End stops
    if ((Len>=1) && !strcmp(s[0],":End")){break;}
    if ((Len>=2) && !strcmp(s[0],":SoilZone")){found=!strcmp(s[1],"PHYSICAL");}
  }
  delete p; RVG.close();
  return found;
}

bool ParseGWFile(CModel*& pModel, const optStruct& Options)
{
  ExitGracefullyIf(pModel==NULL,"ParseGWFile: model has not yet been created.",BAD_DATA);
  CGroundwaterModel *pGW=pModel->GetGroundwaterModel();
  char *s[MAXINPUTITEMS];
  int   Len,line(0);

  ifstream RVG;
  RVG.open(Options.rvg_filename.c_str());
  if (RVG.fail()){ExitGracefully(("Cannot find or read .rvg file "+Options.rvg_filename).c_str(),BAD_DATA); return false;}
  CParser *p=new CParser(RVG,Options.rvg_filename,line);

  if (Options.noisy){
    cout<<"============================================================="<<endl;
    cout<<"Parsing Groundwater File "<<Options.rvg_filename<<"..."<<endl;
    cout<<"============================================================="<<endl;
  }

  string szOpts="";                                               //soil-zone options given (they need :SoilZone PHYSICAL)
  bool end_of_file=p->Tokenize(s,Len);
  while (!end_of_file)
  {
    if (!IsComment(s[0],Len)){const int nmax=GWCommandArity(s,Len); if (nmax>0){GWMaxTokens(s,Len,nmax,Options);}}
    if      (IsComment(s[0],Len)){}
    else if (!strcmp(s[0],":End")){break;}
    else if ((!strcmp(s[0],":FileType")) || (!strcmp(s[0],":Application")) || (!strcmp(s[0],":Version")) ||
             (!strcmp(s[0],":WrittenBy")) || (!strcmp(s[0],":CreationDate")) || (!strcmp(s[0],":Name"))){}
    else if (!strcmp(s[0],":HRUGeometry"))
    {
      if (Len<2){GWImproperFormat(s,Options);}
      else {pGW->SetHRUGeometry(CorrectForRelativePath(s[1],Options.rvg_filename),(Len>=3)?string(s[2]):string("HRU_ID"));}
    }
    else if (!strcmp(s[0],":RiverGeometry"))
    {
      if (Len<2){GWImproperFormat(s,Options);}
      else {pGW->SetRiverGeometry(CorrectForRelativePath(s[1],Options.rvg_filename),(Len>=3)?string(s[2]):string(""));}
    }
    else if (!strcmp(s[0],":LandSurface") || !strcmp(s[0],":BedrockSurface"))
    {/*:LandSurface [file] [METRES], :BedrockSurface [file] [METRES]: METRES states that the raster's elevations are in
       metres for a band unit text Raven cannot place (another length or a scaled metre still stops) */
      if (Len<2){GWImproperFormat(s,Options);}
      else {
        string m=(Len>=3)?string(s[2]):string(""); for (size_t i=0;i<m.size();i++){m[i]=(char)toupper((unsigned char)m[i]);}
        ExitGracefullyIf((Len>=3) && (m!="METRES") && (m!="METERS"),("ParseGWFile: command "+string(s[0])+" in "+Options.rvg_filename+": '"+string(s[2])+
          "' after the file name is not METRES (or METERS, the only word it takes, which confirms that the raster's elevations are in metres)").c_str(),BAD_DATA);
        if (!strcmp(s[0],":LandSurface")){pGW->SetLandSurfaceFile(CorrectForRelativePath(s[1],Options.rvg_filename),Len>=3,s[1]);}
        else                             {pGW->SetBedrockFile    (CorrectForRelativePath(s[1],Options.rvg_filename),Len>=3,s[1]);}
      }
    }
    else if (!strcmp(s[0],":GridCellSize"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetCellSize(GWNum(s,1,Options));}
    }
    //-- existing MODFLOW 6 model ------------------------------------------------------------------
    else if (!strcmp(s[0],":MF6Simulation"))
    {/*:MF6Simulation [mfsim.nam of an existing model] - couples that model instead of generating one */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMF6Simulation(CorrectForRelativePath(s[1],Options.rvg_filename));}
    }
    else if (!strcmp(s[0],":MF6Model"))
    {/*:MF6Model [GWF model name] - required when the simulation holds several GWF models */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMF6Model(s[1]);}
    }
    else if (!strcmp(s[0],":MF6CRS"))
    {/*:MF6CRS EPSG:code | TMERC lat0 lon0 k0 FE FN {WGS84|GRS80} | ALBERS|LCC lat1 lat2 lat0 lon0 FE FN {..} - projection of the
       model's coordinates; Raven then takes the cell outlines from the model itself */
      if (Len<2){GWImproperFormat(s,Options);}
      else {
        string spec; for (int i=1;i<Len;i++){spec+=string(s[i])+" ";}
        string err; if (!pGW->SetMF6CRS(spec,err)){ExitGracefully(("ParseGWFile: :MF6CRS - "+err).c_str(),BAD_DATA);}
      }
    }
    else if (!strcmp(s[0],":InputCRS"))
    {/*:InputCRS EPSG:code | TMERC ... | ALBERS ... | LCC ...  (grammar of :MF6CRS) - projection of the coordinates of every
       GeoJSON and raster input file (and of wells given by X and Y); without it they are longitude/latitude */
      if (Len<2){GWImproperFormat(s,Options);}
      else {
        string spec; for (int i=1;i<Len;i++){spec+=string(s[i])+" ";}
        string err; if (!pGW->SetInputCRS(spec,err)){ExitGracefully(("ParseGWFile: :InputCRS - "+err).c_str(),BAD_DATA);}
      }
    }
    else if (!strcmp(s[0],":MF6GridFile"))
    {/*:MF6GridFile [GeoJSON of the model's top-layer cells, lon/lat] {ID field, default CELL_ID} */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMF6GridFile(CorrectForRelativePath(s[1],Options.rvg_filename),(Len>=3)?s[2]:"CELL_ID");}
    }
    else if (!strcmp(s[0],":MF6StartDate"))
    {/*:MF6StartDate [yyyy-mm-dd] {hh:mm:ss} - calendar date of the model's time zero, with its time of day (also written
       yyyy-mm-ddThh:mm:ss; default 00:00:00); Raven may start then or later: restarts */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMF6StartDate((Len>=3)?string(s[1])+" "+string(s[2]):string(s[1]));}
    }
    else if (!strcmp(s[0],":RavenTakesOver"))
    {/*:RavenTakesOver [package name] {package name ...} - recharge/ET/UZF packages removed and replaced by Raven */
      if (Len<2){GWImproperFormat(s,Options);} else {for (int i=1;i<Len;i++){pGW->AddTakeOver(s[i]);}}
    }
    else if (!strcmp(s[0],":KeepPackages"))
    {/*:KeepPackages [package name] {package name ...} - recharge/ET/UZF packages that stay beside Raven */
      if (Len<2){GWImproperFormat(s,Options);} else {for (int i=1;i<Len;i++){pGW->AddKeepPackage(s[i]);}}
    }
    else if (!strcmp(s[0],":StreamPackage"))
    {/*:StreamPackage [RIV/DRN/GHB package name] {AUX [aux name] | DOMINANT_HRU} - its flows go to Raven's reaches */
      if (Len==2){pGW->SetStreamPackage(s[1],"AUX","SUBBASIN");}
      else if ((Len>=4) && !strcmp(s[2],"AUX")){pGW->SetStreamPackage(s[1],"AUX",s[3]);}
      else if ((Len>=3) && !strcmp(s[2],"DOMINANT_HRU")){pGW->SetStreamPackage(s[1],"DOMINANT_HRU","");}
      else {GWImproperFormat(s,Options);}
    }
    else if (!strcmp(s[0],":SeepageToRaven"))
    {/*:SeepageToRaven ON | OFF - DRN_RAVEN at the top of the top active cells (default ON) */
      if ((Len<2) || (strcmp(s[1],"ON") && strcmp(s[1],"OFF"))){GWImproperFormat(s,Options);} else {pGW->SetExtSeepage(!strcmp(s[1],"ON"));}
    }
    else if (!strcmp(s[0],":OverlapWeights"))
    {/*:OverlapWeights  {[HRU ID] [cell] [w = A_overlap/A_cell]} x n  :EndOverlapWeights  (MODFLOW-USG coupling format) */
      while (!p->Tokenize(s,Len)){
        if (IsComment(s[0],Len)){continue;}
        if (!strcmp(s[0],":EndOverlapWeights")){break;}
        if (!strcmp(s[0],":Attributes") || !strcmp(s[0],":Units")){continue;}
        if (Len<3){GWImproperFormat(s,Options); break;}
        if (Len>3){ExitGracefully(("ParseGWFile: :OverlapWeights row "+string(s[0])+" "+string(s[1])+" in "+Options.rvg_filename+" has the surplus value '"+string(s[3])+"' (each row is HRU ID, cell, weight)").c_str(),BAD_DATA);}
        double w=GWNumber(":OverlapWeights",s[2],Options);
        ExitGracefullyIf((w<0.0) || (w>1.0+1e-9),"ParseGWFile: :OverlapWeights weights must be between 0 and 1",BAD_DATA);
        long long hid=GWID(":OverlapWeights",s[0],Options,"HRU ID");
        if (w>0.0){ //(a row of weight 0 is skipped, whatever its cell number)
          long long lcell=GWID(":OverlapWeights",s[1],Options,"cell number");
          if ((lcell<1) || (lcell>INT_MAX)){
            ExitGracefully(("ParseGWFile: :OverlapWeights row "+string(s[0])+" "+string(s[1])+" in "+Options.rvg_filename+": the cell number "+string(s[1])+
                            " is not between 1 and the number of cells per layer").c_str(),BAD_DATA);
          }
          pGW->AddOverlapWeight(hid,(int)lcell,w);
        }
      }
    }
    else if (!strcmp(s[0],":GridType"))
    {/*:GridType REGULAR | QUADTREE | HRU_MESH | FILE  (default REGULAR) */
      if (Len<2){GWImproperFormat(s,Options);}
      else if (!strcmp(s[1],"REGULAR")) {pGW->SetGridType(GWGRID_REGULAR);}
      else if (!strcmp(s[1],"QUADTREE")){pGW->SetGridType(GWGRID_QUADTREE);}
      else if (!strcmp(s[1],"HRU_MESH")){pGW->SetGridType(GWGRID_HRUMESH);}
      else if (!strcmp(s[1],"FILE"))    {pGW->SetGridType(GWGRID_FILE);}
      else {GWImproperFormat(s,Options);}
    }
    else if (!strcmp(s[0],":GridRotation"))
    {/*:GridRotation [degrees counterclockwise] (regular and quadtree grids) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetGridRotation(GWNum(s,1,Options));}
    }
    else if (!strcmp(s[0],":GridRefinement"))
    {/*:GridRefinement RIVERS [size] | WELLS [size] | LAKES [size] | LINES [file] [size] | POLYGONS [file] [size] */
      if ((Len>=3) && (!strcmp(s[1],"RIVERS") || !strcmp(s[1],"WELLS") || !strcmp(s[1],"LAKES"))){pGW->AddGridRefinement(s[1],"",GWNum(s,2,Options));}
      else if ((Len>=4) && (!strcmp(s[1],"LINES") || !strcmp(s[1],"POLYGONS"))){pGW->AddGridRefinement(s[1],CorrectForRelativePath(s[2],Options.rvg_filename),GWNum(s,3,Options));}
      else {GWImproperFormat(s,Options);}
    }
    else if (!strcmp(s[0],":GridFile"))
    {/*:GridFile [GeoJSON of cell polygons] {ID field (default CELL_ID)} */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetGridFile(CorrectForRelativePath(s[1],Options.rvg_filename),(Len>=3)?s[2]:"CELL_ID");}
    }
    else if (!strcmp(s[0],":QuadtreeMaxLevel"))
    {/*:QuadtreeMaxLevel [halvings of the base cell size] (default: from the refinement sizes) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetQuadtreeMaxLevel(GWInt(s,1,Options));}
    }
    else if (!strcmp(s[0],":FlowCorrection"))
    {/*:FlowCorrection AUTO | XT3D | XT3D_RHS | NONE  (default AUTO: no correction; XT3D options for strongly non-orthogonal cells) */
      if ((Len<2) || (strcmp(s[1],"AUTO") && strcmp(s[1],"XT3D") && strcmp(s[1],"XT3D_RHS") && strcmp(s[1],"NONE"))){GWImproperFormat(s,Options);} else {pGW->SetFlowCorrection(s[1]);}
    }
    else if (!strcmp(s[0],":LayerConnection"))
    {/*:LayerConnection INDEX | ELEVATION  (default INDEX; ELEVATION connects overlapping cells of any layers, DISU) */
      if ((Len<2) || (strcmp(s[1],"INDEX") && strcmp(s[1],"ELEVATION"))){GWImproperFormat(s,Options);} else {pGW->SetLayerConnection(!strcmp(s[1],"ELEVATION"));}
    }
    else if (!strcmp(s[0],":MinCellCoverage"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMinCoverage(GWNum(s,1,Options));}
    }
    else if (!strcmp(s[0],":MF6Library"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetLibraryPath(CorrectForRelativePath(s[1],Options.rvg_filename));}
    }
    else if (!strcmp(s[0],":SoilZone"))
    {/*:SoilZone CONCEPTUAL                      (default: soil and aquifer linked only by recharge and the seepage drain)
       :SoilZone PHYSICAL SOIL[0] ... SOIL[M-1]   (the water table can rise into these soil layers; recharge must come
                                                   from the bottom one, SOIL[M-1]; see the manual, "The soil zone") */
      if ((Len<2) || (strcmp(s[1],"CONCEPTUAL") && strcmp(s[1],"PHYSICAL"))){GWImproperFormat(s,Options);}
      else if (!strcmp(s[1],"CONCEPTUAL")){pGW->SetSoilZone(false,0);}
      else {
        int M=Len-2;
        for (int j=0;j<M;j++){
          string want="SOIL["+to_string(j)+"]";
          ExitGracefullyIf(want!=string(s[2+j]),(":SoilZone PHYSICAL: list the soil-zone layers from the top, SOIL[0] SOIL[1] ...; found "+string(s[2+j])+" where "+want+" was expected").c_str(),BAD_DATA);
        }
        ExitGracefullyIf(M<1,":SoilZone PHYSICAL needs at least one soil layer (SOIL[0])",BAD_DATA);
        pGW->SetSoilZone(true,M);
      }
    }
    else if (!strcmp(s[0],":SoilZoneIteration"))
    {/*:SoilZoneIteration [max updates] {head tolerance m} {flow tolerance -}  (default 100 1E-4 1E-6)
       max updates: interface updates a step may add once MODFLOW itself has converged (0: no iteration, the explicit
       predictor only); tolerances not given keep their defaults */
      if (Len<2){GWImproperFormat(s,Options);}
      else {
        int n=GWInt(s,1,Options); double th=(Len>=3)?GWNum(s,2,Options):1e-4,tq=(Len>=4)?GWNum(s,3,Options):1e-6;
        szOpts+=" :SoilZoneIteration";
        ExitGracefullyIf((n<0) || !(th>0) || !(tq>0),":SoilZoneIteration: max updates >= 0, tolerances > 0",BAD_DATA);
        pGW->SetSoilZoneIteration(n,th,tq);
      }
    }
    else if (!strcmp(s[0],":SoilZoneSubgrid"))
    {/*:SoilZoneSubgrid [DEM samples per cell side]  (default 5; 1 = cell-mean land surface only) */
      if (Len<2){GWImproperFormat(s,Options);}
      else {int n=GWInt(s,1,Options); ExitGracefullyIf((n<1) || (n>20),":SoilZoneSubgrid: samples per cell side must be 1 to 20",BAD_DATA); pGW->SetSoilZoneSubgrid(n); szOpts+=" :SoilZoneSubgrid";}
    }
    else if (!strcmp(s[0],":SoilZoneVadose"))
    {/*:SoilZoneVadose NONE | RICHARDS {accuracy} {nodes} | KINEMATIC {accuracy} {nodes} | LAG {parameter} {nodes}
       (default NONE: recharge reaches the water table in the step). An unsaturated column between the soil-zone base and
       the water table delays and stores recharge: RICHARDS (Richards equation), KINEMATIC (gravity-only kinematic wave),
       LAG (linear storage cascade); accuracy = largest content change in a substep (default 0.02), nodes per column
       (default 16) */
      const int mode=(Len<2)?-1:(!strcmp(s[1],"NONE")?0:(!strcmp(s[1],"RICHARDS")?1:(!strcmp(s[1],"KINEMATIC")?2:(!strcmp(s[1],"LAG")?3:-1))));
      if (mode<0){GWImproperFormat(s,Options);}
      else {
        double acc=0.0,par=0.0; int n=(Len>=4)?GWInt(s,3,Options):0;
        if (mode==3){par=(Len>=3)?GWNum(s,2,Options):0.0; ExitGracefullyIf((Len>=3) && !(par>=0.0),":SoilZoneVadose LAG: the parameter may not be negative",BAD_DATA);}
        else {
          acc=(Len>=3)?GWNum(s,2,Options):0.0;
          ExitGracefullyIf((Len>=3) && !((acc>0.0) && (acc<=0.2)),":SoilZoneVadose: accuracy (largest content change in a substep) must lie in (0, 0.2]",BAD_DATA);
        }
        ExitGracefullyIf((Len>=4) && ((n<3) || (n>200)),":SoilZoneVadose: nodes per column must be 3 to 200",BAD_DATA);
        pGW->SetSoilZoneVadose(mode,acc,n,par); szOpts+=" :SoilZoneVadose";
      }
    }
    else if (!strcmp(s[0],":SoilZoneVadoseClasses"))
    {/*:SoilZoneVadoseClasses K  (default 1: one unsaturated column per soil-zone column, at its mean base; K > 1: up to K
       columns, one per depth class of its DEM samples) */
      if (Len<2){GWImproperFormat(s,Options);}
      else {int K=GWInt(s,1,Options); ExitGracefullyIf((K<1) || (K>10),":SoilZoneVadoseClasses: classes must be 1 to 10",BAD_DATA); pGW->SetSoilZoneVadoseClasses(K); szOpts+=" :SoilZoneVadoseClasses";}
    }
    else if (!strcmp(s[0],":SoilZoneVadoseInitialFlux"))
    {/*:SoilZoneVadoseInitialFlux q  [mm/d] (default 0: the columns start at rest, at their reference profile): on a cold start
       every unsaturated column starts at its steady profile for the downward flux q (not on a groundwater hotstart) */
      if (Len<2){GWImproperFormat(s,Options);}
      else {double q=GWNum(s,1,Options); ExitGracefullyIf(!(q>=0.0),":SoilZoneVadoseInitialFlux: the flux may not be negative",BAD_DATA); pGW->SetSoilZoneVadoseInitialFlux(q); szOpts+=" :SoilZoneVadoseInitialFlux";}
    }
    else if (!strcmp(s[0],":SoilZoneCapillaryRise"))
    {/*:SoilZoneCapillaryRise ON | OFF  (default ON with an unsaturated column, :SoilZoneVadose other than NONE; else OFF) */
      if ((Len<2) || (strcmp(s[1],"ON") && strcmp(s[1],"OFF"))){GWImproperFormat(s,Options);}
      else {pGW->SetSoilZoneCapillaryRise(strcmp(s[1],"ON")?0:1); szOpts+=" :SoilZoneCapillaryRise";}
    }
    else if (!strcmp(s[0],":SoilZoneExchangeTiming"))
    {/*:SoilZoneExchangeTiming AFTER | LAGGED  (default LAGGED: the last step's exchange is a source before Raven's
       processes and only its change follows; AFTER: the whole exchange follows the processes) */
      if ((Len<2) || (strcmp(s[1],"AFTER") && strcmp(s[1],"LAGGED"))){GWImproperFormat(s,Options);}
      else {pGW->SetSoilZoneTiming(strcmp(s[1],"LAGGED")?0:1); szOpts+=" :SoilZoneExchangeTiming";}
    }
    else if (!strcmp(s[0],":SoilZoneClasses"))
    {/*:SoilZoneClasses K  (default 1: Raven's processes act on the HRU's mean soil water; K > 1: on up to K classes of
       its soil-zone columns, grouped by wetness) */
      if (Len<2){GWImproperFormat(s,Options);}
      else {int K=GWInt(s,1,Options); ExitGracefullyIf((K<1) || (K>50),":SoilZoneClasses: 1 to 50 classes",BAD_DATA); pGW->SetSoilZoneClasses(K); szOpts+=" :SoilZoneClasses";}
    }
    else if (!strcmp(s[0],":SoilZoneAquiferClass"))
    {/*:SoilZoneAquiferClass [aquifer class]  (an existing model: CLAPP_B, AIR_ENTRY_PRESSURE and POROSITY of its cells for
       the unsaturated columns; their conductivity and specific yield come from the model) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetSoilZoneAquiferClass(s[1]); szOpts+=" :SoilZoneAquiferClass";}
    }
    else if (!strcmp(s[0],":SoilZoneOutput"))
    {/*:SoilZoneOutput [HRU group]  (HRUs written to GWSoilExchange.csv; default all soil-zone HRUs) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetSoilZoneOutput(s[1]); szOpts+=" :SoilZoneOutput";}
    }
    else if (!strcmp(s[0],":SeepageReturnTo"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetSeepageTarget(s[1]);}
    }
    else if ((!strcmp(s[0],":Wells")) || (!strcmp(s[0],":ObservationWells")))
    {/*:Wells
         :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_TOP, SCREEN_BOT
         {ID, values} x nWells
       :EndWells
       :ObservationWells
         :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_ELEV
       :EndObservationWells
       with :InputCRS, X and Y (projected, m) may replace LATITUDE and LONGITUDE */
      bool obs=!strcmp(s[0],":ObservationWells");
      vector<string> attrs; vector<vector<string> > rows;
      ReadAttributeTable(p,obs?":EndObservationWells":":EndWells",attrs,rows);
      for (size_t r=0;r<rows.size();r++){
        if (obs){
          const string c=":ObservationWells";
          gw_obswell O; O.id=GWID(c,rows[r][0].c_str(),Options); O.name=GetAttr(attrs,rows[r],"NAME",rows[r][0]);
          O.xy=(GetAttr(attrs,rows[r],"X","")!="") && (GetAttr(attrs,rows[r],"LATITUDE","")=="");
          if (O.xy){O.lat=GWNumber(c,GetAttr(attrs,rows[r],"Y","REQUIRED").c_str(),Options); O.lon=GWNumber(c,GetAttr(attrs,rows[r],"X","REQUIRED").c_str(),Options);}
          else {O.lat=GWNumber(c,GetAttr(attrs,rows[r],"LATITUDE","REQUIRED").c_str(),Options); O.lon=GWNumber(c,GetAttr(attrs,rows[r],"LONGITUDE","REQUIRED").c_str(),Options);}
          O.z=GWNumber(c,GetAttr(attrs,rows[r],"SCREEN_ELEV","REQUIRED").c_str(),Options); O.node=-1;
          pGW->AddObservationWell(O);
        }
        else {
          const string c=":Wells";
          gw_well W; W.id=GWID(c,rows[r][0].c_str(),Options); W.name=GetAttr(attrs,rows[r],"NAME",rows[r][0]);
          W.xy=(GetAttr(attrs,rows[r],"X","")!="") && (GetAttr(attrs,rows[r],"LATITUDE","")=="");
          if (W.xy){W.lat=GWNumber(c,GetAttr(attrs,rows[r],"Y","REQUIRED").c_str(),Options); W.lon=GWNumber(c,GetAttr(attrs,rows[r],"X","REQUIRED").c_str(),Options);}
          else {W.lat=GWNumber(c,GetAttr(attrs,rows[r],"LATITUDE","REQUIRED").c_str(),Options); W.lon=GWNumber(c,GetAttr(attrs,rows[r],"LONGITUDE","REQUIRED").c_str(),Options);}
          W.screen_top=GWNumber(c,GetAttr(attrs,rows[r],"SCREEN_TOP","REQUIRED").c_str(),Options); W.screen_bot=GWNumber(c,GetAttr(attrs,rows[r],"SCREEN_BOT","REQUIRED").c_str(),Options);
          W.pRate=NULL;
          ExitGracefullyIf(W.screen_top<W.screen_bot,("ParseGWFile: well "+W.name+" has SCREEN_TOP below SCREEN_BOT").c_str(),BAD_DATA);
          pGW->AddWell(W);
        }
      }
    }
    else if ((!strcmp(s[0],":GeneralHeadBoundary")) || (!strcmp(s[0],":SpecifiedHeadBoundary")))
    {/*:GeneralHeadBoundary   [name] [file.geojson] HEAD [m] CONDUCTANCE [m2/d per cell] {LAYERS ALL|TOP}
       :SpecifiedHeadBoundary [name] [file.geojson] HEAD [m] {LAYERS ALL|TOP}
       polygons select cells whose centre lies inside; lines select cells they cross.
       A :BoundaryHead [name] time series in the .rvt file overrides HEAD. */
      gw_boundary B; B.isCHD=!strcmp(s[0],":SpecifiedHeadBoundary");
      if (Len<3){GWImproperFormat(s,Options);}
      else {
        B.name=s[1]; B.file=CorrectForRelativePath(s[2],Options.rvg_filename);
        B.head=0.0; B.cond=0.0; B.allLayers=true; B.pHead=NULL; bool hasHead=false;
        ExitGracefullyIf((Len-3)%2!=0,("ParseGWFile: boundary "+B.name+" in "+Options.rvg_filename+": the values after the file must come in "
                         "keyword-value pairs (HEAD, CONDUCTANCE, LAYERS); '"+string(s[Len-1])+"' has no partner (a comma separates values: write 1146.85, not 1,146.85)").c_str(),BAD_DATA);
        for (int i=3;i+1<Len;i+=2){
          if      (!strcmp(s[i],"HEAD"       )){B.head=GWNum(s,i+1,Options); hasHead=true;}
          else if (!strcmp(s[i],"CONDUCTANCE")){B.cond=GWNum(s,i+1,Options);}
          else if (!strcmp(s[i],"LAYERS"     )){
            ExitGracefullyIf(strcmp(s[i+1],"TOP") && strcmp(s[i+1],"ALL"),("ParseGWFile: boundary "+B.name+": LAYERS must be ALL or TOP, not "+string(s[i+1])).c_str(),BAD_DATA);
            B.allLayers=strcmp(s[i+1],"TOP")!=0;}
          else {ExitGracefully(("ParseGWFile: boundary "+B.name+" in "+Options.rvg_filename+": unknown keyword '"+string(s[i])+"' (HEAD, CONDUCTANCE or LAYERS)").c_str(),BAD_DATA);}
        }
        ExitGracefullyIf(!hasHead,("ParseGWFile: boundary "+B.name+" needs a HEAD value").c_str(),BAD_DATA);
        ExitGracefullyIf((!B.isCHD) && (B.cond<=0),("ParseGWFile: general-head boundary "+B.name+" needs a positive CONDUCTANCE").c_str(),BAD_DATA);
        pGW->AddBoundary(B);
      }
    }
    else if (!strcmp(s[0],":GWInitialization"))
    {/*:GWInitialization STEADY_STATE [uniform recharge, mm/d] */
      if ((Len<3) || strcmp(s[1],"STEADY_STATE")){GWImproperFormat(s,Options);}
      else {pGW->SetSteadyStateInit(GWNum(s,2,Options));}
    }
    else if (!strcmp(s[0],":RiverBedDepth"))
    {/*:RiverBedDepth [m | REFERENCE_FLOW]  channel bed depth below the DEM along rivers (default REFERENCE_FLOW) */
      if (Len<2){GWImproperFormat(s,Options);}
      else {pGW->SetRiverBedDepth(strcmp(s[1],"REFERENCE_FLOW")?GWNum(s,1,Options):-1.0);}
    }
    else if (!strcmp(s[0],":SurfaceWaterExchange"))
    {/*:SurfaceWaterExchange [state variable] [HRU group] LEAKANCE [1/d]   e.g. DEPRESSION WetlandHRUs LEAKANCE 0.01 */
      if (Len<5){GWImproperFormat(s,Options);}
      else if (strcmp(s[3],"LEAKANCE")){GWImproperFormat(s,Options);}
      else {pGW->AddSurfaceWaterExchange(s[1],s[2],GWNum(s,4,Options));}
    }
    else if (!strcmp(s[0],":ReservoirExchange"))
    {/*:ReservoirExchange   couple Raven reservoirs (with :SeepageParameters, on lake HRUs that have an AQUIFER_PROFILE)
         to the aquifer: MODFLOW heads set each reservoir's local groundwater head; its seepage enters the aquifer */
      pGW->SetReservoirExchange(true);
    }
    else if (!strcmp(s[0],":LinkageCache"))
    {/*:LinkageCache [file]  where the HRU/grid intersection is cached (default <output>/mf6/linkage.cache) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetLinkageCache(CorrectForRelativePath(s[1],Options.rvg_filename));}
    }
    else if (!strcmp(s[0],":WriteNetCDFHeads"))
    {/*:WriteNetCDFHeads   full head field to GWHeads.nc every :HeadSaveFrequency steps (needs -Dnetcdf) */
      pGW->SetWriteNetCDFHeads(true);
    }
    else if (!strcmp(s[0],":WriteGWBudgetFile"))
    {/*:WriteGWBudgetFile   MODFLOW budget file mf6/gwf.cbc every :HeadSaveFrequency steps, with mf6/gwf.porosity (generated models) */
      pGW->SetWriteBudgetFile();
    }
    else if (!strcmp(s[0],":SolverHeadTolerance"))
    {/*:SolverHeadTolerance [m]  Newton head-change tolerance (default 0.00001) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetHeadTolerance(GWNum(s,1,Options));}
    }
    else if (!strcmp(s[0],":SolverMaxOuterIterations"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMaxOuterIterations(GWInt(s,1,Options));}
    }
    else if (!strcmp(s[0],":HeadSaveFrequency"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetHeadSaveFrequency(GWInt(s,1,Options));}
    }
    else if (!strcmp(s[0],":NameFile")         || !strcmp(s[0],":GWRiverConnection") || !strcmp(s[0],":Drains") ||
             !strcmp(s[0],":ETProperties")     || !strcmp(s[0],":RavenCBCUnitNumber")|| !strcmp(s[0],":Recharge") ||
             !strcmp(s[0],":RechargeOptionCode")|| !strcmp(s[0],":NumNodes")         || !strcmp(s[0],":NumDrains"))
    { //commands of Raven's MODFLOW-USG coupling (Raven 4.15 and earlier): refused, naming the MODFLOW 6 commands
      ExitGracefully(("ParseGWFile: "+string(s[0])+" belongs to the MODFLOW-USG coupling of Raven 4.15 and earlier, which this version "
                      "replaces with MODFLOW 6. To couple an existing MODFLOW 6 model use :MF6Simulation (an :OverlapWeights "
                      "table in the same format is accepted); see the manual chapter \"Coupling an existing MODFLOW 6 model\"").c_str(),BAD_DATA);
    }
    else
    {
      WriteWarning("ParseGWFile: unrecognized command "+string(s[0])+" in "+Options.rvg_filename+" (ignored)",Options.noisy);
    }
    end_of_file=p->Tokenize(s,Len);
  }
  RVG.close();
  delete p;
  pGW->ResolveInputCRS();
  if ((szOpts!="") && !pGW->UsesPhysicalSoilZone()){
    WriteWarning("ParseGWFile:"+szOpts+" ignored: these apply only with :SoilZone PHYSICAL",Options.noisy);
  }
  else if (szOpts!=""){pGW->WarnIgnoredSoilZoneOptions(szOpts,Options.noisy);}
  return true;
}
