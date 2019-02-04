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


#include "MCHMappingInterface/Segmentation.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"


namespace o2
{
namespace mch
{

/// \brief Test of the Digitization
/// A couple of values are filled into a Hits and we check whether we get reproducible output in terms of digits
/// and MClabels

  
BOOST_AUTO_TEST_CASE(Digitizer_test1)
{
  
  Digitizer digitizer;
  //test amplification
  
 
  //need to produce hits with proper MC labels
  int trackId1 = 0;
  int trackId2 = 1;
  short detElemId1 = 101;//check what to put
  short detElemId2 = 1012;
  Point3D<float> entrancePoint1(-17.7993, 8.929883, -522.201); //x,y,z coordinates in cm
  /*
     hit.GetEnergyLoss() 1.17596e-06
        [22200]: charge 157.925
[22200]: detID 101
[22200]: indexID 1
[22200]: mSeg[indexID].detElemId() 101
[22200]: lpos.X() 17.9843 lpos.Y() 9.44883 lpos.Z() -0.209473
[22200]: fraction charge 0.977048
[22200]:  mSeg[indexID].isBendingPad(padid) 1
[22200]: mSeg[indexID].bending().nofPads() 14392
[22200]:  padid 2598
[22200]: bending: xmin -0.343334 xmax 0.286666 ymin -0.211174 ymax 0.208826
[22200]: padid 2598
[22200]: signal 114.19
[22200]: labelIndex 0
[22200]:  mSeg[indexID].isBendingPad(padid) 0
[22200]: mSeg[indexID].bending().nofPads() 14392
[22200]:  padid 16981
[22200]: nonbending: xmin -0.0283337 xmax 0.601666 ymin -0.001174 ymax 0.418826
[22200]: padid 16981
[22200]: signal 45.4551
[22200]: labelIndex 0
[22200]:  mSeg[indexID].isBendingPad(padid) 0
[22200]: mSeg[indexID].bending().nofPads() 14392
[22200]:  padid 16982
[22200]: nonbending: xmin -0.0283337 xmax 0.601666 ymin -0.421174 ymax -0.001174
[22200]: padid 16982
[22200]: signal 44.8909
[22200]: labelIndex 0
[22200]:  mSeg[indexID].isBendingPad(padid) 0
[22200]: mSeg[indexID].bending().nofPads() 14392
[22200]:  padid 16997
[22200]: nonbending: xmin -0.658334 xmax -0.0283337 ymin -0.421174 ymax -0.001174
[22200]: padid 16997
[22200]: signal 32.7726
[22200]: labelIndex 0
[22200]:  mSeg[indexID].isBendingPad(padid) 0
[22200]: mSeg[indexID].bending().nofPads() 14392
[22200]:  padid 16996
[22200]: nonbending: xmin -0.658334 xmax -0.0283337 ymin -0.001174 ymax 0.418826
[22200]: padid 16996
[22200]: signal 33.1845
[22200]: labelIndex 0
   */
  Point3D<float> exitPoint1(-17.8136, 8.93606, -522.62);
  Point3D<float> entrancePoint2(-49.2793, 28.8673, -1441.25);
  Point3D<float> exitPoint2(-49.2965, 28.8806, -1441.75);
  float eloss1 = 1e-6;
  float eloss2 =1e-6;
  float length = 0.f;//no ida what it is good for
  float tof = 0.0;//not used
  /*
hit.GetX() -49.2793 hit.GetY() 28.8673hit.GetZ() -1441.25 hit.GetEnergyLoss() 1.12173e-06
[22200]: hit.GetTrackID() 0
[22200]: hit.entrancePoint().x() -49.2793 hit.entrancePoint().y() 28.8673 hit.entrancePoint().z() -1441.25
[22200]: hit.exitPoint().x() -49.2965 hit.exitPoint().y() 28.8806 hit.exitPoint().z() -1441.75
[22200]: charge 145.878
[22200]: detID 1012
[22200]: indexID 142
[22200]: mSeg[indexID].detElemId() 1012
[22200]: lpos.X() -71.9707 lpos.Y() -9.28491 lpos.Z() -0.249413
[22200]: fraction charge 1.01618
[22200]:  mSeg[indexID].isBendingPad(padid) 1
[22200]: mSeg[indexID].bending().nofPads() 4096
[22200]:  padid 1301
[22200]: bending: xmin -10 xmax 0 ymin -0.284906 ymax 0.215094
[22200]: padid 1301
[22200]: signal 58.5768
[22200]: labelIndex 9
[22200]:  mSeg[indexID].isBendingPad(padid) 1
[22200]: mSeg[indexID].bending().nofPads() 4096
[22200]:  padid 1237
[22200]: bending: xmin 0 xmax 10 ymin -0.284906 ymax 0.215094
[22200]: padid 1237
[22200]: signal 58.5768
[22200]: labelIndex 9
[22200]:  mSeg[indexID].isBendingPad(padid) 0
[22200]: mSeg[indexID].bending().nofPads() 4096
[22200]:  padid 4931
[22200]: nonbending: xmin -2.38419e-07 xmax 0.714285 ymin -9.28491 ymax 0.715094
[22200]: padid 4931
[22200]: signal 71.2336
[22200]: labelIndex 9
[22200]:  mSeg[indexID].isBendingPad(padid) 0
[22200]: mSeg[indexID].bending().nofPads() 4096
[22200]:  padid 4933
[22200]: nonbending: xmin -0.714286 xmax -2.38419e-07 ymin -9.28491 ymax 0.715094
[22200]: padid 4933
[22200]: signal 71.2334
[22200]: labelIndex 9
   */
  //could also check to give the same input and see whether I get the same output as well
    
  std::vector<Hit> hits(2);
  vector.at(0) = Hit(trackId1, detElemId1, entrancePoint1, exitPoint1, eloss1, length, tof);//put some values
  vector.at(1) = Hit(trackId2, detElemId2, entrancePoint2, exitPoint2, eloss2, length, tof);//put some values
  // one hit per station, if feasible and energy deposition such that from 1 to 4 pad digits all included
  //
  //test first only single processHit
  
  std::vector<Digit> digits;
  mapping::Segmentation seg1{ detElemId1 };
  mapping:: Segmentation seg2{ detElemId2 };
  process(hits, &digits);
  //digit members: 
  //retrieve information from digits: getPad(), getADC(), getLabel()
  //compare Hit
  int digitcounter1 = 0;
  int digitcounter2 = 1;
  
  for (auto& digit : digits) {
    
    int padid = digit.getPad();
    int adc = digit.getADC();
    int label = digit.getLabel();
    
    if(label == trackId1)
      {
	bool check = seg1.isValid(digit.getPad());// is pad ID unique across full detector?
	//check true
	//adc
	BOOST_CHECK_CLOSE();
	digitcounter1++;
      } else if (label == trackId2)
      {
	digitcounter2++;
      } else
      {
	//some boost functionality, if not  one of two values
      };
    
  }
  //check both digitcounters being between 1 and 3

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
