// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \file testDigitisation.cxx
/// \brief This task tests the Digitizer and the Response of the MCH digitization
/// \author Michael Winn, DPhN/IRFU/CEA, michael.winn@cern.ch

#define BOOST_TEST_MODULE Test MCH Digitization
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <memory>
#include <vector>
#include "MCHSimulation/Digit.h"
#include "MCHSimulation/Digitizer.h"
#include "MCHSimulation/Hit.h"


namespace o2
{
namespace mch
{

/// \brief Test of the Digitization
/// A couple of values are filled into a Hits and we check whether we get reproducible output in terms of digits
/// and MClabels

BOOST_AUTO_TEST_CASE(Response_test)
{
  //Problem: function is not deterministic!
  //check transition between energy and charge?
  //check integration via Mathieson
  //check FEE response

  
}
  
BOOST_AUTO_TEST_CASE(Digitizer_test1)
{
  
  Digitizer digitizer;
  //test amplification
  
 
  //need to produce hits with proper MC labels
  int trackId1 = 0;
  int trackId2 = 1;
  short detElemId1 = 101;//check what to put
  short detElemId2 = 102;
  Point3D<float> entrancePoint1(0.f, 1.f, 2.f); //x,y,z coordinates
  Point3D<float> exitPoint1(0.f, 1.f, 2.f);
  Point3D<float> entrancePoint2(0.f, 1.f, 2.f);
  Point3D<float> exitPoint2(0.f, 1.f, 2.f);
  float eloss1 = 1e-6;
  float eloss2 =1e-6;
  float length = 0.f;//no ida what it is good for
  float tof = 0.0;//not used
  //could also check to give the same input and see whether I get the same output as well
    
  std::vector<Hit> hits(2);
  vector.at(0) = Hit(trackId1, detElemId1, entrancePoint1, exitPoint1, eloss1, length, tof );//put some values
  vector.at(1) = Hit(trackId2, detElemId2, entrancePoint2, exitPoint2, eloss2, length, tof);//put some values
  // one hit per station, if feasible and energy deposition such that from 1 to 4 pad digits all included
  //
  //test first only single processHit
  
  std::vector<Digit> digits;
  process(hits, &digits);
  
  //what to test:
  //1) compare label of hit and of MCtruthContainer, certainly makes sense
  //1) check condition that equal or more than 1 digit per hit: ok, but not very stringent
  //2) compare position of digit with hit position within a certain precision?
  // (would rely on segmentation to map back from pad number to position, not really this unit)
  //3) something with charge size? Not really useful.
  //could also check charge summation for different hits acting on same pad...

  //should one introduce member constants and getters to display intermediate steps?
} 

}
}
