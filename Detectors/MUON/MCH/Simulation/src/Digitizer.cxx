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
  // To be done:
  //1) add a test
  //2) check read-out chain efficiency handling in current code: simple efficiency values? masking?
  //different strategy?
  //3) handling of time dimension: what changes w.r.t. aliroot check HMPID
  //4) handling of MCtruth information
}

//______________________________________________________________________

void Digitizer::process(const std::vector<Hit> hits, std::vector<Digit>& digits)
{
  digits.clear();
  mDigits.clear();
  mMCTruthContainer.clear();
  
  //array of MCH hits for a given simulated event
  //to be added: exception handling
  for (auto& hit : hits) {
    int labelIndex = mMCTruthContainer.getIndexedSize();//index for this hit
    processHit(hit, mEventTime, labelIndex);

    MCCompLabel label(hit.GetTrackID(), mEventID, mSrcID);
    //Todo:: check/update DPL workflow how EventID and srcID set
    mMCTruthContainer.addElementRandomAccess(labelIndex,label);//ToDo: check if needed like this

    auto labels = mMCTruthContainer.getLabels(labelIndex);
    std::sort(labels.begin(), labels.end());
  } // end loop over hits
  std::cout << "mMCTruthContainer.getIndexedSize() " << mMCTruthContainer.getIndexedSize() << std::endl;
  std::cout << "MMCTruthContainer.getNElements() " << mMCTruthContainer.getNElements() << std::endl;
  //TODO: merge (new member function, add charge) of digits that are on same pad:
  //things to think about in terms of time costly

  fillOutputContainer(digits);
}

//______________________________________________________________________
int Digitizer::processHit(const Hit& hit, double event_time, int labelIndex)
{
  //Question: what happens if several tracks producing signal on same pad up to now? -> should be additive? is this the case?
  //labelIndex: to be done:
  //assign the labelIndex to every digit that corresponds to
  //this hit
  //if multiple hits on producing same digit or multiple
  //hit position(cm)
  Point3D<float> pos(hit.GetX(), hit.GetY(), hit.GetZ());
  //hit has global coordinates
  //convert energy to charge, float enough?
  float charge = mMuonresponse.etocharge(hit.GetEnergyLoss());
  //std::cout << "charge " << charge << std::endl;

  //time information
  float time = hit.GetTime(); //to be used for pile-up
  int detID = hit.GetDetectorID();
  //get index for this detID
  int indexID = mdetID.at(detID);
  //# digits for hit
  int ndigits = 0;

  //transformation from global to local
  auto t = o2::mch::getTransformation(detID, *gGeoManager);
  Point3D<float> lpos;
  t.MasterToLocal(pos, lpos);

  float anodpos = mMuonresponse.getAnod(lpos.X(), detID);

  //TODO/Questions:
  //- possibility to do random seeding in controlled way?
  //distance of impact point and anode, needed for charge sharing
  float anoddis = TMath::Abs(pos.X() - anodpos);
  //should be related to electrons fluctuating out/in one/both halves (independent x)
  float fracplane = mMuonresponse.chargeCorr();
  //should become a function of anoddis
  //std::cout << "fracplane " << fracplane << std::endl;
  float chargebend = fracplane * charge;
  float chargenon = charge / fracplane;
  //last line  from Aliroot, not understood why
  //since charge = charchbend+chargenon and not multiplication
  float signal = 0.0;

  //borders of charge gen.
  double xMin = anodpos - mMuonresponse.getQspreadX() * 0.5;
  double xMax = anodpos + mMuonresponse.getQspreadX() * 0.5;
  double yMin = lpos.Y() - mMuonresponse.getQspreadY() * 0.5;
  double yMax = lpos.Y() + mMuonresponse.getQspreadY() * 0.5;

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
    return 0; //to be decided if needed

  //need to keep both electros separated since afterwards signal generation on each plane
  //how misalignment enters?
  std::vector<int> padIDs;

  //retrieve pads with signal
  mSeg[indexID].forEachPadInArea(xMin, yMin, xMax, yMax, [&padIDs](int padid) { padIDs.emplace_back(padid); });

  //induce signal pad-by-pad: bending
  for (auto& padid : padIDs) {
    //squeeze both loop, in one if statement for signal
    //retrieve coordinates for each pad
    xmin = (anodpos - mSeg[indexID].padPositionX(padid)) - mSeg[indexID].padSizeX(padid) * 0.5;
    xmax = xmin + mSeg[indexID].padSizeX(padid);
    ymin = (lpos.Y() - mSeg[indexID].padPositionY(padid)) - mSeg[indexID].padSizeY(padid) * 0.5;
    ymax = ymin + mSeg[indexID].padSizeY(padid);

    // 1st step integrate induced charge for each pad
    if (mSeg[indexID].isBendingPad(padid))
      {
	signal = mMuonresponse.chargePad(xmin, xmax, ymin, ymax, detID, chargebend);
      } else
      {
	signal = mMuonresponse.chargePad(xmin, xmax, ymin, ymax, detID, chargenon);
      }
    // if(signal>mMuonresponse.getChargeThreshold()
    //&&     signal<mMuonresponse.getChargeSat()//not yet tuned
    //  ) {
    //translate charge in signal
    signal = mMuonresponse.response(signal, detID);
    //write digit
    mDigits.emplace_back(padid, signal, labelIndex);
    //to do check if label is really making sense
    //i.e. what happens if multiple hits for one pad
    ++ndigits;
    //  }
  } 
  return ndigits;
}
//______________________________________________________________________
//not clear if needed for DPL or modifications required
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
  //why needed and not directly in MCTruthContainer for EMCal?
  //check how this is written into the output in DPL workflow
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
  //works only if called after process (ugly!!!!!)
  mcContainer.clear();
  if (mMCTruthOutputContainer.getNElements()==0)//tbc if working...
      return;
  
  for (int index =0; index < mMCTruthOutputContainer.getIndexedSize(); ++index) {
    mcContainer.addElements(index, mMCTruthOutputContainer.getLabels(index));
  }
  
  mMCTruthOutputContainer.clear();
}
