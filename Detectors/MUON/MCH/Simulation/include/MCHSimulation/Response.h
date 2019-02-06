// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MCH_SIMULATION_RESPONSE_H_
#define O2_MCH_SIMULATION_RESPONSE_H_

#include "MCHSimulation/Digit.h"
#include "MCHSimulation/Detector.h"
#include "MCHSimulation/Hit.h"

namespace o2
{
namespace mch
{

class Response
{
 public:
  Response(int station = 0): mStation(station) {

    if (mStation==0)
      {
	//parameter for Mathieson station 1
	mK2x = 1.021026;
	mSqrtK3x = 0.7000;
	mK4x = 0.40934890;
	mK2y = 0.9778207;
	mSqrtK3y = 0.7550;
	mK4y = 0.38658194;
	//inverse anode-cathode Pitch in 1/cm, station 1
	mInversePitch = 1. / 0.21;
      } else if((mStation>0)&&(mStation<5))
    {
      //parameter for Mathieson station 2-5
      mK2x = 1.010729;
      mSqrtK3x = 0.7131;
      mK4x = 0.40357476;
      mK2y = 0.970595;
      mSqrtK3y = 0.7642;
      mK4y = 0.38312571;
      //inverse anode-cathode Pitch in 1/cm, station 2-5
      mInversePitch = 1. / 0.25;
    } else
    {
      LOG(ERROR) << "station number for response out of bound for MCH ";
      return;
    }
  };
  
  ~Response() = default;

  float getQspreadX() { return mQspreadX; };
  float getQspreadY() { return mQspreadY; };
  float getChargeSat() { return mChargeSat; };
  float getChargeThreshold() { return mChargeThreshold; };
  float etocharge(float edepos);
  double chargePadfraction(float xmin, float xmax, float ymin, float ymax);
  double chargefrac1d(float min, float max, double k2, double sqrtk3, double k4);
  double response(float charge);
  float getAnod(float x);
  float chargeCorr();

 private:
  //parameter for station number
  int mStation = 0;
  //proper parameter in aliroot in AliMUONResponseFactory.cxx
  const float mQspreadX = 0.144; //charge spread in cm
  const float mQspreadY = 0.144;
  
  //ChargeSlope for Station 2-5
  const float mChargeSlope = 25;  //why float in Aliroot?
  const float mChargeCorr = 0.11; // number from line 122
  //of AliMUONResponseFactory.cxx

  const float mChargeThreshold = 1e-4;
  //AliMUONResponseV0.cxx constr.
  const float mChargeSat = 0.61 * 1.25 * 0.2;
  //from AliMUONResponseV0.cxx
  //equals AliMUONConstants::DefaultADC2MV()*AliMUONConstants::DefaultA0()*AliMUONConstants::DefaultCapa()
  //Mathieson parameter: NIM A270 (1988) 602-603
  //should be a common place for MCH
  // Mathieson parameters from L.Kharmandarian's thesis, page 190
  //  fKy2 = TMath::Pi() / 2. * (1. - 0.5 * fSqrtKy3);//AliMUONMathieson::SetSqrtKx3AndDeriveKx2Kx4(Float_t SqrtKx3)
  //  Float_t cy1 = fKy2 * fSqrtKy3 / 4. / TMath::ATan(Double_t(fSqrtKy3));
  //  fKy4 = cy1 / fKy2 / fSqrtKy3; //this line from AliMUONMathieson::SetSqrtKy3AndDeriveKy2Ky4
  //why this multiplicitation before again division? any number small compared to Float precision?
  double mK2x = 0.0;
  double mSqrtK3x = 0.0;
  double mK4x = 0.0;
  double mK2y = 0.0;
  double mSqrtK3y = 0.0;
  double mK4y = 0.0;

  //anode-cathode Pitch in 1/cm
  float mInversePitch = 0.0;
};
} // namespace mch
} // namespace o2
#endif
