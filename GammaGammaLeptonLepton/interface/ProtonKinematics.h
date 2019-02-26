#ifndef ProtonKinematics_h
#define ProtonKinematics_h

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

/**
 * \date Jun 2016
 * \author Jan Kaspar <jan.kaspar@cern.ch>
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
class ProtonKinematics
{
 public:
  /**
   * \param[in] arm_id 0 (F), 1 (N)
   * \param[in] side 0 (L), 1 (R)
   */
  inline ProtonKinematics(unsigned int run_id, unsigned short arm_id, unsigned short side, const TotemRPLocalTrack& lt) :
    fX0(lt.getX0()), fY0(lt.getY0()), fIsValid(lt.isValid()), fArm(arm_id), fSide(side)
  {
    if (run_id<274244) fRunPeriod = 0;
    else fRunPeriod = 1;
  }
  inline bool isValid() const { return fIsValid; }
  inline double getXi() const {
    const double x_corr = fX0+getDeX();
    return x_corr*1e-3/getD();
  }

 private:
  // alignment corrections (in m)
  inline double getDeX() const {
    switch (fSide) {
      case 0: default: { return getDeX_L(); }
      case 1:          { return getDeX_R(); }     
    }
  }
  inline double getDeX_L() const {
    switch (fRunPeriod) {
      case 0:          { return (fArm==0) ? -3.40 : -0.90; }
      case 1: default: { return (fArm==0) ? -1.45 : -3.90; }
    }
  }
  inline double getDeX_R() const {
    switch (fRunPeriod) {
      case 0:          { return (fArm==0) ? -2.75 : -2.40; }
      case 1: default: { return (fArm==0) ? -2.85 : -3.25; }
    }
  }
  // optics properties (in m)
  inline double getD() const {
    switch (fSide) {
      case 0: default: { return getD_L(); }
      case 1:          { return getD_R(); }
    }
  }
  inline double getD_L() const { return (fArm==0) ? 9.22e-2 : 9.26e-2; }
  inline double getD_R() const { return (fArm==0) ? 5.81e-2 : 5.16e-2; }

  double fX0, fY0;
  bool fIsValid;
  unsigned short fRunPeriod, fArm, fSide;
 
};

#endif
