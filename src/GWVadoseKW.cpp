/*----------------------------------------------------------------
  Raven Library Source Code
  Copyright (c) 2008-2026 the Raven Development Team
  ----------------------------------------------------------------
  GWVadoseKW.cpp
  Not implemented: see GWVadose.h; the column passes its inflow straight through.
  ----------------------------------------------------------------*/
#include "GWVadose.h"

namespace gwvz
{
  double StepKinematic(vz_column &c,const vz_strat &S,const double dt,const double qtop,double &rejected,vz_step_info *info)
  {
    rejected=0.0; if (info){info->substeps=1; info->newton=0; info->fails=0;}
    if (qtop<0.0){rejected=qtop*dt; return 0.0;}
    return qtop;
  }
}
