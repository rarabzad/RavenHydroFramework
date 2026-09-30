#-----------------------------------------------------------------
# DataType         Raven RVP file
# For Liard Watershed Model, Canada
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# Soil Classes
#-----------------------------------------------------------------
:SoilClasses
  :Attributes,
  :Units,  
   TOPSOIL,      1.0,    0.0,       0
   SLOW_RES,     1.0,    0.0,       0
   FAST_RES,     1.0,    0.0,       0
   SLOW_RES2,    1.0,    0.0,       0
   FAST_RES2,    1.0,    0.0,       0
   SLOW_RES3,    1.0,    0.0,       0
   FAST_RES3,    1.0,    0.0,       0
   TOP_WETLAND,  1.0,    0.0,       0
   RES_WETLAND,  1.0,    0.0,       0
   RES2_WETLAND, 1.0,    0.0,       0
   TOPSOIL2,     1.0,    0.0,       0
   VADOSE,       1.0,    0.0,       0
:EndSoilClasses

#-----------------------------------------------------------------
# Land Use Classes
#-----------------------------------------------------------------
:LandUseClasses, 
  :Attributes,        IMPERM,    FOREST_COV, 
       :Units,          frac,          frac, 
  FOREST,               0.05,           1.0,
  BARREN,               0.15,           0.0, 
  SHRUBLAND,            0.10,           0.0, 
  GRASSLAND,             0.0,           0.0,
  GLACIER,               0.0,           0.0,
  WETLAND,               0.0,           0.0,
  WATER,                 0.0,           0.0,
:EndLandUseClasses

#-----------------------------------------------------------------
# Vegetation Classes
#-----------------------------------------------------------------
:VegetationClasses, 
  :Attributes,        MAX_HT,       MAX_LAI, MAX_LEAF_COND, 
       :Units,             m,          none,      mm_per_s, 
       FOREST,            17,            6,            5.0, 
       BARREN,           0.0,          0.0,            0.0,
    SHRUBLAND,           4.0,            5,            5.0,
    GRASSLAND,           0.6,            2,            5.0,
      WETLAND,           4.0,            5,            5.0,
        WATER,           0.0,          0.0,            0.0,
      GLACIER,           0.0,          0.0,            0.0,
:EndVegetationClasses

#-----------------------------------------------------------------
# Soil Profiles
#-----------------------------------------------------------------
:SoilProfiles
         LAKE, 0
         ROCK, 0
      GLACIER, 0
      WETLAND, 4, TOP_WETLAND, 1.0, VADOSE, 1.0 , RES_WETLAND, 100.0, RES2_WETLAND,   1.0,
	DEFAULT_P, 4, TOPSOIL    , 1.0, VADOSE, 1.0 , FAST_RES   , 100.0, SLOW_RES    , 100.0, 
	  SOIL_P2, 4, TOPSOIL    , 1.0, VADOSE, 1.0 , FAST_RES2  , 100.0, SLOW_RES2   , 100.0, 
	  SOIL_P3, 4, TOPSOIL2   , 1.0, VADOSE, 10.0, FAST_RES3  , 100.0, SLOW_RES3   , 100.0,  
	P_WETLAND, 4, TOP_WETLAND, 1.0, VADOSE, 1.0 , RES_WETLAND, 100.0, RES2_WETLAND,   1.0, 
:EndSoilProfiles


#-----------------------------------------------------------------
# Global Parameters
#-----------------------------------------------------------------
:GlobalParameter        RAINSNOW_TEMP 0.0 
:GlobalParameter       RAINSNOW_DELTA 2.0 
:GlobalParameter      ADIABATIC_LAPSE 5.810948E+00
:GlobalParameter             SNOW_SWI 0.055
:GlobalParameter     AVG_ANNUAL_RUNOFF 280
:GlobalParameter          PRECIP_LAPSE 1.002140E+00

#-----------------------------------------------------------------
# Soil Parameters
#-----------------------------------------------------------------
:SoilParameterList
:Parameters ,        POROSITY    ,        HBV_BETA    ,  FIELD_CAPACITY    ,        SAT_WILT, MAX_CAP_RISE_RATE,   MAX_PERC_RATE    ,  BASEFLOW_COEFF    ,      BASEFLOW_N    , 
    :Units ,            none    ,            none    ,            none    ,            none,              mm/d,            mm/d    ,             1/d    ,            none    ,
    [DEFAULT],        1.013835E-01,        6.750959E-01,        7.741520E-01,             0.0,               0.0,             0.0    ,             0.0    ,             0.0    , 
    FAST_RES ,        _DEFAULT    ,        _DEFAULT    ,        _DEFAULT    ,             0.0,          _DEFAULT,        2.676131E+00,        6.928413E-02,        1.143030E+00, 
    SLOW_RES ,        _DEFAULT    ,        _DEFAULT    ,        _DEFAULT    ,             0.0,          _DEFAULT,        _DEFAULT    ,        2.258343E-02,             1.0    , 
    FAST_RES2,       _DEFAULT     ,        _DEFAULT    ,        _DEFAULT    ,             0.0,          _DEFAULT,             1.5    ,            0.08    ,            1.25    , 
    SLOW_RES2,       _DEFAULT     ,        _DEFAULT    ,        _DEFAULT    ,             0.0,          _DEFAULT,        _DEFAULT    ,            0.02    ,             1.0    ,  
TOP_WETLAND  ,        2.310000E-01,        5.590000E-01,        5.000000E-01,             0.0,               0.0,        _DEFAULT    ,             0.0    ,             0.0    , 
RES_WETLAND  ,        _DEFAULT    ,        _DEFAULT    ,        _DEFAULT    ,             0.0,          _DEFAULT,        _DEFAULT    ,        1.995250E-01,        1.489806E+00, 
RES2_WETLAND ,       _DEFAULT     ,        _DEFAULT    ,        _DEFAULT    ,             0.0,          _DEFAULT,        _DEFAULT    ,             0.0    ,             1.0    , 
    TOPSOIL2 ,        3.000000E-01,        5.000000E-01,        2.000000E-01,             0.0,               0.0,        1.000000E+01,             0.0    ,             1.0    , 
    FAST_RES3,        _DEFAULT    ,       _DEFAULT     ,        _DEFAULT    ,             0.0,          _DEFAULT,        _DEFAULT    ,        7.655112E-02,        1.502327E+00, 
    SLOW_RES3,        _DEFAULT    ,       _DEFAULT     ,        _DEFAULT    ,             0.0,          _DEFAULT,        _DEFAULT    ,        2.064129E-02,             1.0    ,    
:EndSoilParameterList

:SoilParameterList
  :Parameters,    MAX_BASEFLOW_RATE,      BASEFLOW_THRESH,
      :Units,                 mm/d,                 none,
   [DEFAULT],                    0,                    0,
   TOPSOIL2,              1.999336E+01,               3.149430E-01,
:EndSoilParameterList

#-----------------------------------------------------------------
# Land Use Parameters
#-----------------------------------------------------------------
:LandUseParameterList
  :Parameters,     MELT_FACTOR, MIN_MELT_FACTOR, HBV_MELT_ASP_CORR,  HBV_MELT_FOR_CORR, REFREEZE_FACTOR,
         :Units,        mm/d/K,          mm/d/K,              none,               none,          mm/d/K,
[DEFAULT],       2.010079E+00,       1.007513E+00,          3.001893E-01,           4.520247E-01,             2.0,
   FOREST,       2.834272E+00,       1.697670E+00,          5.993902E-01,           3.279659E-01,        _DEFAULT, 
  GLACIER,      2.007555E+00 ,       1.002627E+00,          _DEFAULT    ,           5.199134E-01,        _DEFAULT,
   BARREN,       1.292391E+00,       1.070814E+00,          _DEFAULT    ,           4.625021E-01,             3.0,
SHRUBLAND,     2.124607E+00  ,       2.012541E+00,          _DEFAULT    ,           4.236297E-01,        _DEFAULT,
  WETLAND,       5.496716E+00,       3.999620E+00,          _DEFAULT    ,           8.603706E-01,            1.22,  
:EndLandUseParameterList

:LandUseParameterList
  :Parameters,     HBV_MELT_GLACIER_CORR, HBV_GLACIER_KMIN, GLAC_STORAGE_COEFF,  HBV_GLACIER_AG,
  :Units,                           none,              1/d,                1/d,            1/mm,
   [DEFAULT],                          1,             0.35,               0.35,            0.05,
:EndLandUseParameterList

:LandUseParameterList
  :Parameters,   DEP_THRESHOLD,  DEP_MAX_FLOW,   DEP_N,  DEP_MAX,    DEP_SEEP_K,   
  :Units,                    mm,          mm/d,    none,       mm,           1/d,              
   [DEFAULT],          1.997862E+02,      1.039893E+02,    1.995070E+00,   4.868635E+02,    7.138236E-03,           
:EndLandUseParameterList

:LandUseParameterList
  :Parameters, OW_PET_CORR,
  :Units,             none,
   [DEFAULT],          0.8,
:EndLandUseClasses

#-----------------------------------------------------------------
# Vegetation Parameters
#-----------------------------------------------------------------
:VegetationParameterList
  :Parameters,    SAI_HT_RATIO,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,    MAX_CAPACITY, MAX_SNOW_CAPACITY,
       :Units,            none,            none,            none,              mm,                mm,
    [DEFAULT],             2.0,            0.05,            0.05,              10,             10000,
       FOREST,        _DEFAULT,        _DEFAULT,        _DEFAULT,           10000,             10000,
        WATER,               0,               0,               0,               0,                 0,
    GRASSLAND,             4.0,        _DEFAULT,        _DEFAULT,        _DEFAULT,          _DEFAULT,
       BARREN,               0,               0,               0,               0,                 0,
      GLACIER,               0,               0,               0,               0,                 0,
:EndVegetationParameterList

#
# ----Channel profiles-------------------------
:RedirectToFile channel_sections.rvp


#-----------------------------------------------------------------
# Groundwater (Raven-MODFLOW 6): hypothetical aquifer for testing
#-----------------------------------------------------------------
:AquiferClasses
  :Attributes, K_HORIZ, K_VERT, SPEC_STORAGE, SPEC_YIELD, POROSITY
  :Units,      m/d,     m/d,    1/m,          none,       none
  GRAVEL,      par_K_gravel, 5.0, 1.0E-5, par_Sy_gravel, 0.30
  CLAY_TILL,   0.001,   0.0001, 1.0E-4,       0.03,       0.40
  SAND,        10.0,    1.0,    1.0E-5,       0.20,       0.35
  UPLAND_TILL, par_K_till, 0.02, 1.0E-4, 0.08, 0.30
:EndAquiferClasses

:AquiferProfiles
  VALLEY_ALLUVIUM, 3, GRAVEL, 15.0, AQUIFER, CLAY_TILL, 5.0, AQUITARD, SAND, TO_BEDROCK, CONFINED_AQUIFER
  UPLAND_TILL,     1, UPLAND_TILL, TO_BEDROCK, AQUIFER
:EndAquiferProfiles

:AquiferProfileParameters
  :Attributes,     INITIAL_HEAD_DEPTH, DEFAULT_BEDROCK_DEPTH, SEEPAGE_LEAKANCE, SOIL_ZONE_DEPTH
  :Units,          m,                  m,                     1/d,              m
  VALLEY_ALLUVIUM, 3.0,                60.0,                  1.0,              0.0
  UPLAND_TILL,     8.0,                15.0,                  1.0,              0.0
:EndAquiferProfileParameters
