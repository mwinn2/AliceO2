// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Dimuon task
//
// Author: Michael Winn

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "Analysis/Dimuon.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DerivedMuonTask{
  Produces<o2::aod::DerivedMuons> derivedmuons; 

  void process(aod::Muons const& muonTracks){

     for (auto it0 = muonTracks.begin(); it0 != muonTracks.end(); ++it0) {
       derivedmuons(it0.inverseBendingMomentum(), it0.thetaX(), it0.thetaY());
       //Todo: add MFT information
     }
  }

};


struct DiMuonTask {
  Produces<o2::aod::Dimuons> dimuons;

  void process(aod::Collision const& collision,
               aod::DerivedMuons const& derivedmuonTracks)
  {
    float mMuonmassSquared = 0.10566*0.10566;
    //to run twice the loop makes very little sense just for deriving...
    //would be better to calculate things on the fly
    for (auto it0 = derivedmuonTracks.begin(); it0 != derivedmuonTracks.end(); ++it0) {
      auto& muon_0 = *it0;

      for(auto it1 = it0 +1; it1 != derivedmuonTracks.end(); ++it1){
	
	auto& muon_1 = *it1;
	dimuons(collision,
		sqrtf(it0.px()*it0.px() + it0.py()*it0.py() + it0.pz()*it0.pz() - mMuonmassSquared)
		+sqrtf(it1.px()*it1.px() + it1.py()*it1.py() + it1.pz()*it1.pz() - mMuonmassSquared),
		it0.px()+it1.px(), it0.py()+it1.py(), it0.pz()+it1.pz(),
		//it0.globalIndex(),
		it0.pt(), it0.eta(), it0.phi(), it0.charge(),
		//it1.globalIndex(),
		it1.pt(), it1.eta(), it1.phi(), it1.charge()
		);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<DerivedMuonTask>("derivedMuonsbuilder"),
    adaptAnalysisTask<DiMuonTask>("diMuonBuilder")};
}
