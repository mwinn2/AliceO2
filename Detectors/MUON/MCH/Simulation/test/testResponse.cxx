// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \file testResponse.cxx                                                                                                                    
/// \brief This task tests of the Response of the MCH digitization
/// \author Michael Winn, DPhN/IRFU/CEA, michael.winn@cern.ch

#define BOOST_TEST_MODULE Test MCHSimulation Response
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <memory>
#include <vector>
#include "MCHSimulation/Response.h"

#include "TH1D.h"
#include "TF1.h"

namespace o2
{
namespace MCH
{

BOOST_AUTO_TEST_CASE(Response_test)
{
  //Problem: function is not deterministic!
  //check transition between energy and charge?
  //check integration via Mathieson
  //check FEE response, station 1 and 2-5 only different for Mathieson integral for the time being
  Response r_stat1 = Response(0);//station 1
  Response r_stat2 = Response(1);//station 2-5
  //check threshold
  float threshold = r_stat1.getChargeThreshold();
  float threshold_target = 1e-4;
  float threshold_precision = threshold_target/10.f;
  BOOST_CHECK_CLOSE(threshold, threshold, threshold_target, threshold_precision);
  
  //check conversion energy to charge
  float eloss = 1e-6; //typical example from stepper for one hit, average
  TH1D hTest("hTest", "", 10000, 0, 1000);
  TF1 gaus("gaus", "gaus");
  for(Int_t = 0; i<100000; i++) {
    hTest.Fill(r_stat1.etatocharge(eloss));
  }

  hTest.Fit("gaus","Q");
  
  float charge_target = 134.f;//taken from one run with Aliroot parameters, tbc
  float charge_precision = charge_target/10.f; //tbc
  BOOST_CHECK_CLOSE(gaus.GetParameter(1),charge_target,charge_precision);
  //TODO, if needed
  float charge_resolution_target = 1.f*charge_target;
  float charge_resolution_precision = charge_resolution*2.f;
  BOOST_CHECK_CLOSE(gauss.GetParameter(2), charge_resolution_target,charge_resolution_precision);
  
  //test Mathieson integration on Pad
  //total charge to be distributed
  float charge_on_plane = 100.f;
  //borders, one wire pitch distance to hit
  //smallest pad dimension for  stat.1 (bending plane)
  float xmin = 0.21;
  float xmax = 0.0;
  float ymin =0.0;
  float ymax =0.0;

  //
  float expected_chargeonpad = 0.0;//hard coded
  //or from some alternative implementation of Mathieson?
  float chargeonpad_precision = expected_chargeonpad/100.f;
  float result_chargeonpad =  r_stat1.chargePad(xmin,xmax,ymin,ymax,charge_on_plane);
  BOOST_CHECK_CLOSE(result_chargeonpad,expected_chargeonpad,chargeonpad_precision);
  //todo test r_stat2
  
  //todo some test of chargeCorr? function
  
  
}

