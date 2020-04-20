// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "Analysis/Dimuon.h"



using namespace o2;
using namespace o2::framework;

struct DimuonAnalysis {
  OutputObj<TH1F> hDimuonMass{"mass"};

  void init(InitContext const&)
  {
    hDimuonMass.setObject(new TH1F("mass", "dimuon mass;mass (GeV/#it{c^{2}})",
				   1000, 0.,12.));
  }

  void process(aod::dimuon const& dimuon)
  {
    //    hDimuonMass->Fill(dimuon.mass());
  }
  
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<DimuonAnalysis>("dimuon-analysis")};

}
  
