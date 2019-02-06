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

  return chargefrac1d(xmin, xmax, mK2x, mSqrtK3x, mK4x) * chargepad1d(ymin, ymax, mK2y, mSqrtK3y, mK4y); 
}
//______________________________________________________________________
double Response::chargefrac1d(float min, float max, double k2, double sqrtk3, double k4){
  // The Mathieson function integral (1D)
  double u1 = sqrtk3 * TMath::TanH( k2 * xmin );
  double u2 = sqrtk3 * TMath::TanH( k2 * xmax );

  return 2. * k4 * (TMath::ATan(u2) - TMath::ATan(u2));
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
