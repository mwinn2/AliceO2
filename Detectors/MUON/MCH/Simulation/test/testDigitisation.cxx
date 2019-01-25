// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/// \file testTPCDigitContainer.cxx
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
#include "MCHSimulation/Response.h"

namespace o2
{
namespace MCH
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
  //need to produce hits with proper MC labels
  std::vector<Hit> hits = ;
  // one hit per station, if feasible and energy deposition such that from 1 to 4 pad digits all included
  //                            
  
  //need to produce 

  //what to test:
  //1) compare label of hit and of MCtruthContainer, certainly makes sense
  //1) check condition that equal or more than 1 digit per hit
  //2) compare position of digit with hit position within a certain precision?
  // (would rely on segmentation to map back from pad number to position)
  //3) something with charge size? Not really useful
  //could also check charge summation for different hits acting on same pad...

  //should one introduce member constants and getters to display intermediate steps?
  //1) response translation to charge? Can be done already now with Response 
  //2) size of charge deposition?
  //3) something else
  } 

}
}
