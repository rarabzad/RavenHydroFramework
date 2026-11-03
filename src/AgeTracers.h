/*----------------------------------------------------------------
Raven Library Source Code
Copyright (c) 2008-2026 the Raven Development Team
----------------------------------------------------------------*/

#ifndef AGE_TRACERS_H
#define AGE_TRACERS_H

#include "RavenInclude.h"
#include "Transport.h"

///////////////////////////////////////////////////////////////////
/// \brief Class for coordinating age tracer simulation (child of ConstituentModel)
/// \details Implemented in AgeTracers.cpp
//
class CAgeTracer :public CConstituentModel
{
private:
   double _model_time; //local storage of tt.model_time, measured from the datum of the entry times (_t_diff + tt.model_time)
   double _t_diff;     //[d] start of this run after the datum of the entry times (the start of the run that wrote the
                       //hotstart file; 0 on a cold start): the masses keep that datum, the clock is shifted

public:/*-------------------------------------------------------*/
  CAgeTracer(CModel *pMod,CTransportModel *pTMod,string name,const int c);
  ~CAgeTracer();

  double CalculateReportingConcentration(const double &mass, const double &vol) const; //returns age [d]
  double ConvertConcentration(const double &delta) const;
  double GetOutflowConcentration(const int p) const;

  double GetAdvectionCorrection(const CHydroUnit* pHRU,const int iFromWater,const int iToWater,const double& mass, const double &vol, const double &Q) const;

  void   SetModelTime(const double &t);
  void   SetInitialTimeOffset(const double &tdiff);
  double GetTimeOffset() const {return _t_diff;} ///< [d] datum of the entry times before this run's start

  //Manipulators (inherited from CConstitModel)
  void   Initialize                 (const optStruct &Options);
  void   UpdateInitialConditions    (const optStruct &Options);
  void   WriteEnsimOutputFileHeaders(const optStruct &Options);
  void   WriteEnsimMinorOutput      (const optStruct &Options,const time_struct &tt);
};
#endif
