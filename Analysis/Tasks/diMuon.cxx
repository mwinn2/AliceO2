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
  Produces<o2::aod::DerivedMuons> derivedmuons;

  void process(aod::Collision const& collision,
               aod::Muons const& muonTracks)
  {
    float muonmass = 0.10566;
    for (auto it0 = muonTracks.begin(); it0 != tracks.end(); ++it0) {
      auto& muon_0 = *it0;
      aod::DerivedMuons muondev_0 = DerivedMuons(it0.InverseBendingMomentum(), it0.ThetaX(), it0.ThetaY());
      derivedmuons(it0.InverseBendingMomentum(), it0.ThetaX(), it0.ThetaY());
      for(auto it1 = it0 +1; it1 != tracks.end(); ++it1){
	auto& muon_1 = *it1;
	aod::DerivedMuons muondev_1 = DerivedMuons(it1.InverseBendingMomentum(), it1.ThetaX(), it1.ThetaY());
	//https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODDimuon.h
	//https://github.com/alisw/AliRoot/blob/master/STEER/AOD/AliAODDimuon.cxx
	dimuons(collision,
		TMath::Sqrt(TMath::Sqrt(muondev_0.Px()*muondev_0.Px() + muondev_0.Py()*muondev_0.Py() + muondev_0.Pz()*muondev_0.Pz() - muonmass*muonmass )
			    + TMath::Sqrt(muondev_1.Px()*muondev_1.Px() + muondev_1.Py()*muondev_1.Py() + muondev_1.Pz()*muondev_1.Pz() - muonmass*muonmass)
			    ),
		muondev_0.Px()+muondev_1.Px(), muondev_0.Py()+muondev_1.Py(), muondev_0.Pz()+muondev_1.Pz(),
		muondev_0.globalIndex(), muondev_0.Pt(), muondev_0.Eta(), muondev_0.Phi(),
		muondev_1.globalIndex(), muondev_1.Pt(), muondev_1.Eta(), muondev_1.Phi());
      }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<diMuonTask>("diMuonTask")};
}
