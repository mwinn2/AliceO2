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

struct DiMuonTask {
  Produces<o2::aod::Dimuons> dimuons;

  // TODO: use abs, eventually use pt
  // TODO: use values from configurables
  // TODO: add eta cuts
  Filter trackCuts = (aod::track::signed1Pt < 10.f) &&
                     (aod::track::signed1Pt > -10.f);

  std::unique_ptr<fastjet::BackgroundEstimatorBase> bge;
  std::unique_ptr<fastjet::Subtractor> sub;

  std::vector<fastjet::PseudoJet> pJets;

  void process(aod::Collision const& collision,
               aod::Muons const& muonTracks)
  {

    for (auto it0 = muonTracks.begin(); it0 != tracks.end(); ++it0) {
      auto& muon_0 = *it0;
      for(auto it1 = it0 +1; it1 != tracks.end(); ++it1){
	auto& muon_1 = *it1;
	
	float mass = invmassdimuon();//invmassdimuon to be implemented
	//https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODDimuon.h
	//https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODDimuon.cxx
      }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<diMuonTask>("diMuonTask")};
}
