/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
  Licensed under the Artistic License 2.0 (see LICENSE)
  SPDX-License-Identifier: Artistic-2.0
  Parsing of the groundwater (.rvg) file for the Raven-MODFLOW 6 coupling
----------------------------------------------------------------*/
#include "RavenInclude.h"
#include "GWGrid.h"
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
///  :GridCellSize     [m]                                        (default 1000)
///  :MinCellCoverage  [fraction]                                 (default 0.25)
///  :GridType REGULAR|QUADTREE|HRU_MESH|FILE, :GridRotation, :GridRefinement, :GridFile, :QuadtreeMaxLevel,
///  :FlowCorrection, :LayerConnection  (grid options, see the manual)
///  :MF6Library       [path to libmf6]                           (optional; else RAVEN_MF6_LIB)
///  :SeepageReturnTo  [state variable]                           (default: lowest soil layer)
///  :HeadSaveFrequency [time steps]                              (default 30)
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
    ExitGracefullyIf(Len!=(int)attrs.size()+1,(string("ParseGWFile: wrong number of values in ")+endtag+" table").c_str(),BAD_DATA);
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
bool ParseGWFile(CModel*& pModel, const optStruct& Options)
{
  ExitGracefullyIf(pModel==NULL,"ParseGWFile: model has not yet been created.",BAD_DATA);
  CGroundwaterModel *pGW=pModel->GetGroundwaterModel();
  char *s[MAXINPUTITEMS];
  int   Len,line(0);

  ifstream RVG;
  RVG.open(Options.rvg_filename.c_str());
  if (RVG.fail()){cout<<"ERROR opening *.rvg file: "<<Options.rvg_filename<<endl; return false;}
  CParser *p=new CParser(RVG,Options.rvg_filename,line);

  if (Options.noisy){
    cout<<"============================================================="<<endl;
    cout<<"Parsing Groundwater File "<<Options.rvg_filename<<"..."<<endl;
    cout<<"============================================================="<<endl;
  }

  bool end_of_file=p->Tokenize(s,Len);
  while (!end_of_file)
  {
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
    else if (!strcmp(s[0],":LandSurface"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetLandSurfaceFile(CorrectForRelativePath(s[1],Options.rvg_filename));}
    }
    else if (!strcmp(s[0],":BedrockSurface"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetBedrockFile(CorrectForRelativePath(s[1],Options.rvg_filename));}
    }
    else if (!strcmp(s[0],":GridCellSize"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetCellSize(s_to_d(s[1]));}
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
    else if (!strcmp(s[0],":MF6GridFile"))
    {/*:MF6GridFile [GeoJSON of the model's top-layer cells, lon/lat] {ID field, default CELL_ID} */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMF6GridFile(CorrectForRelativePath(s[1],Options.rvg_filename),(Len>=3)?s[2]:"CELL_ID");}
    }
    else if (!strcmp(s[0],":MF6StartDate"))
    {/*:MF6StartDate [yyyy-mm-dd] - calendar date of the model's time zero (Raven may start then or later: restarts) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMF6StartDate(s[1]);}
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
        double w=s_to_d(s[2]);
        ExitGracefullyIf((w<0.0) || (w>1.0+1e-9),"ParseGWFile: :OverlapWeights weights must be between 0 and 1",BAD_DATA);
        if (w>0.0){pGW->AddOverlapWeight(s_to_ll(s[0]),s_to_i(s[1]),w);}
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
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetGridRotation(s_to_d(s[1]));}
    }
    else if (!strcmp(s[0],":GridRefinement"))
    {/*:GridRefinement RIVERS [size] | WELLS [size] | LAKES [size] | LINES [file] [size] | POLYGONS [file] [size] */
      if ((Len>=3) && (!strcmp(s[1],"RIVERS") || !strcmp(s[1],"WELLS") || !strcmp(s[1],"LAKES"))){pGW->AddGridRefinement(s[1],"",s_to_d(s[2]));}
      else if ((Len>=4) && (!strcmp(s[1],"LINES") || !strcmp(s[1],"POLYGONS"))){pGW->AddGridRefinement(s[1],CorrectForRelativePath(s[2],Options.rvg_filename),s_to_d(s[3]));}
      else {GWImproperFormat(s,Options);}
    }
    else if (!strcmp(s[0],":GridFile"))
    {/*:GridFile [GeoJSON of cell polygons] {ID field (default CELL_ID)} */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetGridFile(CorrectForRelativePath(s[1],Options.rvg_filename),(Len>=3)?s[2]:"CELL_ID");}
    }
    else if (!strcmp(s[0],":QuadtreeMaxLevel"))
    {/*:QuadtreeMaxLevel [halvings of the base cell size] (default: from the refinement sizes) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetQuadtreeMaxLevel(s_to_i(s[1]));}
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
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMinCoverage(s_to_d(s[1]));}
    }
    else if (!strcmp(s[0],":MF6Library"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetLibraryPath(CorrectForRelativePath(s[1],Options.rvg_filename));}
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
       :EndObservationWells */
      bool obs=!strcmp(s[0],":ObservationWells");
      vector<string> attrs; vector<vector<string> > rows;
      ReadAttributeTable(p,obs?":EndObservationWells":":EndWells",attrs,rows);
      for (size_t r=0;r<rows.size();r++){
        if (obs){
          gw_obswell O; O.id=s_to_ll(rows[r][0].c_str()); O.name=GetAttr(attrs,rows[r],"NAME",rows[r][0]);
          O.lat=atof(GetAttr(attrs,rows[r],"LATITUDE","REQUIRED").c_str()); O.lon=atof(GetAttr(attrs,rows[r],"LONGITUDE","REQUIRED").c_str());
          O.z=atof(GetAttr(attrs,rows[r],"SCREEN_ELEV","REQUIRED").c_str()); O.node=-1;
          pGW->AddObservationWell(O);
        }
        else {
          gw_well W; W.id=s_to_ll(rows[r][0].c_str()); W.name=GetAttr(attrs,rows[r],"NAME",rows[r][0]);
          W.lat=atof(GetAttr(attrs,rows[r],"LATITUDE","REQUIRED").c_str()); W.lon=atof(GetAttr(attrs,rows[r],"LONGITUDE","REQUIRED").c_str());
          W.screen_top=atof(GetAttr(attrs,rows[r],"SCREEN_TOP","REQUIRED").c_str()); W.screen_bot=atof(GetAttr(attrs,rows[r],"SCREEN_BOT","REQUIRED").c_str());
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
        for (int i=3;i+1<Len;i+=2){
          if      (!strcmp(s[i],"HEAD"       )){B.head=s_to_d(s[i+1]); hasHead=true;}
          else if (!strcmp(s[i],"CONDUCTANCE")){B.cond=s_to_d(s[i+1]);}
          else if (!strcmp(s[i],"LAYERS"     )){B.allLayers=strcmp(s[i+1],"TOP")!=0;}
          else {WriteWarning("ParseGWFile: unknown boundary keyword "+string(s[i]),Options.noisy);}
        }
        ExitGracefullyIf(!hasHead,("ParseGWFile: boundary "+B.name+" needs a HEAD value").c_str(),BAD_DATA);
        ExitGracefullyIf((!B.isCHD) && (B.cond<=0),("ParseGWFile: general-head boundary "+B.name+" needs a positive CONDUCTANCE").c_str(),BAD_DATA);
        pGW->AddBoundary(B);
      }
    }
    else if (!strcmp(s[0],":GWInitialization"))
    {/*:GWInitialization STEADY_STATE [uniform recharge, mm/d] */
      if ((Len<3) || strcmp(s[1],"STEADY_STATE")){GWImproperFormat(s,Options);}
      else {pGW->SetSteadyStateInit(s_to_d(s[2]));}
    }
    else if (!strcmp(s[0],":RiverBedDepth"))
    {/*:RiverBedDepth [m | REFERENCE_FLOW]  channel bed depth below the DEM along rivers (default REFERENCE_FLOW) */
      if (Len<2){GWImproperFormat(s,Options);}
      else {pGW->SetRiverBedDepth(strcmp(s[1],"REFERENCE_FLOW")?s_to_d(s[1]):-1.0);}
    }
    else if (!strcmp(s[0],":SurfaceWaterExchange"))
    {/*:SurfaceWaterExchange [state variable] [HRU group] LEAKANCE [1/d]   e.g. DEPRESSION WetlandHRUs LEAKANCE 0.01 */
      if (Len<5){GWImproperFormat(s,Options);}
      else if (strcmp(s[3],"LEAKANCE")){GWImproperFormat(s,Options);}
      else {pGW->AddSurfaceWaterExchange(s[1],s[2],s_to_d(s[4]));}
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
    else if (!strcmp(s[0],":SolverHeadTolerance"))
    {/*:SolverHeadTolerance [m]  Newton head-change tolerance (default 0.00001) */
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetHeadTolerance(s_to_d(s[1]));}
    }
    else if (!strcmp(s[0],":SolverMaxOuterIterations"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetMaxOuterIterations(s_to_i(s[1]));}
    }
    else if (!strcmp(s[0],":HeadSaveFrequency"))
    {
      if (Len<2){GWImproperFormat(s,Options);} else {pGW->SetHeadSaveFrequency(s_to_i(s[1]));}
    }
    else if (!strcmp(s[0],":NameFile")         || !strcmp(s[0],":GWRiverConnection") || !strcmp(s[0],":Drains") ||
             !strcmp(s[0],":ETProperties")     || !strcmp(s[0],":RavenCBCUnitNumber")|| !strcmp(s[0],":Recharge") ||
             !strcmp(s[0],":RechargeOptionCode")|| !strcmp(s[0],":NumNodes")         || !strcmp(s[0],":NumDrains"))
    { //commands of the earlier MODFLOW-USG coupling (preserved, outside the build, in legacy/mfusg/)
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
  return true;
}
