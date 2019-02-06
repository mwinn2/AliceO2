// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHSimulation/Digitizer.h"

#include "TMath.h"
#include "TProfile2D.h"
#include "TRandom.h"

#include <algorithm>
#include <cassert>

using namespace o2::mch;

namespace
{
std::map<int, int> createDEMap()
{
  std::map<int, int> m;
  int i{ 0 };
  o2::mch::mapping::forEachDetectionElement([&m, &i](int deid) {
    m[deid] = i++;
  });
  return m;
}

std::vector<o2::mch::mapping::Segmentation> createSegmentations()
{
  std::vector<o2::mch::mapping::Segmentation> segs;

  o2::mch::mapping::forEachDetectionElement([&segs](int deid) {
    segs.emplace_back(deid);
  });
  return segs;
}
} // namespace

Digitizer::Digitizer(int) : mdetID{ createDEMap() }, mSeg{ createSegmentations() }
{
}

void Digitizer::init()
{
}

//______________________________________________________________________
void Digitizer::process(const std::vector<Hit> hits, std::vector<Digit>& digits)
{
  digits.clear();
  mDigits.clear();
  mMCTruthContainer.clear();
  
  //array of MCH hits for a given simulated event
  for (auto& hit : hits) {
    int labelIndex = mMCTruthContainer.getIndexedSize();
    //index for this hit
    int detID = hit.GetDetectorID();
    if(isStation1(detID)){
      processHit(hit, detID, mMuonresponse_stat1, mEventTime, labelIndex);
    } else {
      processHit(hit, detID, mMuonresponse_stat2, mEventTime, labelIndex);
    }
 
    MCCompLabel label(hit.GetTrackID(), mEventID, mSrcID);
    mMCTruthContainer.addElementRandomAccess(labelIndex,label);
    auto labels = mMCTruthContainer.getLabels(labelIndex);
    std::sort(labels.begin(), labels.end());
  } //loop over hits

  fillOutputContainer(digits);
}
//______________________________________________________________________
int Digitizer::processHit(const Hit& hit, int detID, Response response, double event_time, int labelIndex)
{
  //hit position(cm),hit has global coordinates
  Point3D<float> pos(hit.GetX(), hit.GetY(), hit.GetZ());
  
  //convert energy to charge
  float charge = response.etocharge(hit.GetEnergyLoss());
  //time information
  float time = hit.GetTime();

  //get index for this detID
  int indexID = mdetID.at(detID);

  //# digits for hit
  int ndigits = 0;

  //transformation from global to local
  auto t = o2::mch::getTransformation(detID, *gGeoManager);
  Point3D<float> lpos;
  t.MasterToLocal(pos, lpos);
  float anodpos = response.getAnod(lpos.X());
  //  float anoddis = TMath::Abs(pos.X() - anodpos);
  float fracplane = response.chargeCorr();
  //should become a function of anoddis
  float chargebend = fracplane * charge;
  float chargenon = charge / fracplane;
  float signal = 0.0;

  //borders of charge gen.
  double xMin = anodpos - response.getQspreadX() * 0.5;
  double xMax = anodpos + response.getQspreadX() * 0.5;
  double yMin = lpos.Y() - response.getQspreadY() * 0.5;
  double yMax = lpos.Y() + response.getQspreadY() * 0.5;

  //pad-borders
  float xmin = 0.0;
  float xmax = 0.0;
  float ymin = 0.0;
  float ymax = 0.0;

  //use DetectorID to get area for signal induction
  //single pad as check
  int padidbendcent = 0;
  int padidnoncent = 0;
  bool padexists = mSeg[indexID].findPadPairByPosition(anodpos, lpos.Y(), padidbendcent, padidnoncent);
  if (!padexists)
    return 0;

  //pad-id vecotr
  std::vector<int> padIDs;

  //retrieve pads with signal
  mSeg[indexID].forEachPadInArea(xMin, yMin, xMax, yMax, [&padIDs](int padid) { padIDs.emplace_back(padid); });

  //induce signal pad-by-pad: bending
  for (auto& padid : padIDs) {
    //retrieve coordinates for each pad
    xmin = (anodpos - mSeg[indexID].padPositionX(padid)) - mSeg[indexID].padSizeX(padid) * 0.5;
    xmax = xmin + mSeg[indexID].padSizeX(padid);
    ymin = (lpos.Y() - mSeg[indexID].padPositionY(padid)) - mSeg[indexID].padSizeY(padid) * 0.5;
    ymax = ymin + mSeg[indexID].padSizeY(padid);
    // 1st step integrate induced charge for each pad
    if (mSeg[indexID].isBendingPad(padid))
      {
	signal = response.chargePadfraction(xmin, xmax, ymin, ymax) * chargebend;
      } else
      {
	signal = response.chargePadfraction(xmin, xmax, ymin, ymax) * chargenon;
      }
    // if(signal>mMuonresponse.getChargeThreshold()
    //&&     signal<mMuonresponse.getChargeSat()
    //  ) {
    //translate charge in signal
    signal = response.response(signal);
    //write digit
    mDigits.emplace_back(padid, signal, labelIndex);
     ++ndigits;
    //  }
  }
  return ndigits;
}
//______________________________________________________________________
void Digitizer::fillOutputContainer(std::vector<Digit>& digits)
{
  // filling the digit container
  if (mDigits.empty())
    return;

  auto itBeg = mDigits.begin();
  auto iter = itBeg;
  for (; iter != mDigits.end(); ++iter) {
    digits.emplace_back(*iter);
  }
  mDigits.erase(itBeg, iter);
  mMCTruthOutputContainer.clear();
  for (int index =0; index < mMCTruthContainer.getIndexedSize(); ++index) {
    mMCTruthOutputContainer.addElements(index, mMCTruthContainer.getLabels(index));
  }
}
//______________________________________________________________________
void Digitizer::setSrcID(int v)
{
  //set current MC source ID
  if (v > MCCompLabel::maxSourceID()) {
    LOG(FATAL) << "MC source id " << v << " exceeds max storable in the label " << MCCompLabel::maxSourceID() << FairLogger::endl;
  }
  mSrcID = v;
}
//______________________________________________________________________
void Digitizer::setEventID( int v )
{
  // set current MC event ID
  if (v > MCCompLabel::maxEventID()) {
    LOG(FATAL) << "MC event id " << v << " exceeds max storabel in the label " << MCCompLabel::maxEventID() << FairLogger::endl;
  }
  mEventID = v; 
}
//______________________________________________________________________
void Digitizer::provideMC(o2::dataformats::MCTruthContainer<o2::MCCompLabel>& mcContainer)
{
  //fill MCtruth info
  mcContainer.clear();
  if (mMCTruthOutputContainer.getNElements()==0)
    return;
  
  for (int index =0; index < mMCTruthOutputContainer.getIndexedSize(); ++index) {
    mcContainer.addElements(index, mMCTruthOutputContainer.getLabels(index));
  }
  
  mMCTruthOutputContainer.clear();
}
