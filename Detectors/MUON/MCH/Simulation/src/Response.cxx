// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/** @file Response.cxx
 * C++ MCH charge induction and signal generation incl. Mathieson.
 * constants and functions taken from Aliroot.
 * @author Michael Winn, Laurent Aphecetche
 */

#include "MCHSimulation/Response.h"

#include "TMath.h"
#include "TRandom.h"

using namespace o2::mch;

//_____________________________________________________________________
void Response::init()
{
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
  
}


//_____________________________________________________________________
float Response::etocharge(float edepos)
{
  //Todo convert in charge in number of electrons
  //equivalent if IntPH in AliMUONResponseV0 in Aliroot
  //to be clarified:
  //1) why effective parameterisation with Log?
  //2) any will to provide random numbers
  //3) Float in aliroot, Double needed?
  //with central seed to be reproducible?
  //TODO: dependence on station
  //TODO: check slope meaning in thesis
  int nel = int(edepos * 1.e9 / 27.4);
  float charge = 0;
  if (nel == 0)
    nel = 1;
  for (int i = 1; i <= nel; i++) {
    float arg = 0.;
    while (!arg)
      arg = gRandom->Rndm();
    charge -= mChargeSlope * TMath::Log(arg);
  }
  //translate to fC roughly, equivalent to AliMUONConstants::DefaultADC2MV()*AliMUONConstants::DefaultA0()*AliMUONConstants::DefaultCapa() multiplication in aliroot
  charge *= 0.61 * 1.25 * 0.2; // put this in header as constants?
  return charge;
}
//_____________________________________________________________________
double Response::chargePadfraction(float xmin, float xmax, float ymin, float ymax)
{
  //see AliMUONResponseV0.cxx (inside DisIntegrate)
  // and AliMUONMathieson.cxx (IntXY)
  //see: https://edms.cern.ch/ui/file/1054937/1/ALICE-INT-2009-044.pdf
  // normalise w.r.t. Pitch

  xmin *= mInversePitch;
  xmax *= mInversePitch;
  ymin *= mInversePitch;
  ymax *= mInversePitch;

  // The Mathieson function integral
  double ux1 = mSqrtK3x * TMath::TanH(mK2x * xmin);
  double ux2 = mSqrtK3x * TMath::TanH(mK2x * xmax);
  double uy1 = mSqrtK3y * TMath::TanH(mK2y * ymin);
  double uy2 = mSqrtK3y * TMath::TanH(mK2y * ymax);

  return 4. * mK4x * (TMath::ATan(ux2) - TMath::ATan(ux1)) *
         mK4y * (TMath::ATan(uy2) - TMath::ATan(uy1));
}
//______________________________________________________________________
double Response::response(float charge)
{
  //to be done: calculate from induced charge signal
  return charge;
}
//______________________________________________________________________
float Response::getAnod(float x)
{
  int n = Int_t(x / mInversePitch);
  float wire = (x > 0) ? n + 0.5 : n - 0.5;
  return mInversePitch * wire;
}
//______________________________________________________________________
float Response::chargeCorr()
{
  //taken from AliMUONResponseV0
  //conceptually not at all understood why this should make sense
  //mChargeCorr not taken
  return TMath::Exp(gRandom->Gaus(0.0, mChargeCorr / 2.0));
}
