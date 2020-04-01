// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// table definitions of dimuon-candidate for quarkonium analyses
//
// Author: Michael Winn
#ifndef O2_ANALYSIS_DIMUON_H_
#define O2_ANALYSIS_DIMUON_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace dimuon
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Eta, eta, float, "fEta");
DECLARE_SOA_COLUMN(Y, y, float, "fY");
DECLARE_SOA_COLUMN(Phi, phi, float, "fPhi");
DECLARE_SOA_COLUMN(Pt, pt, float, "fPt");
DECLARE_SOA_COLUMN(Mass, mass, float, "fMass");
DECLARE_SOA_INDEX_COLUMN_FULL(Index0, index0, int, MuonTracks, "fIndex0");
DECLARE_SOA_COLUMN(Pt0, pt0, float, "fPt0");
DECLARE_SOA_COLUMN(Eta0, eta0, float, "fEta0");
DECLARE_SOA_COLUMN(Phi0, phi0, float, "fPhi0");
DECLARE_SOA_COLUMN(Charge0, charge0, bool, "fCharge0");
DECLARE_SOA_INDEX_COLUMN_FULL(Index1, index1, int, MuonTracks, "fIndex1");
DECLARE_SOA_COLUMN(Pt1, pt1, float, "fPt1");
DECLARE_SOA_COLUMN(Eta1, eta1, float, "fEta1");
DECLARE_SOA_COLUMN(Phi1, phi1, float, "fPhi1");
DECLARE_SOA_COLUMN(Charge1, charge1, bool, "fCharge1");
//TODO add:
//MTR/MID information
//option for cluster information (see Z-analysis)
//2ndary vertexing compatibiliy
//dielectron analysis compatibility: same derived table from different inputs
} // namespace dimuon

DECLARE_SOA_TABLE(Dimons, "AOD", "DIMUON",
                  o2::soa::Index<>, dimuon::CollisionId,
                  dimuon::Eta, dimuon::Y, dimuon::Phi, dimuon::Pt, dimuon::Mass,
		  dimuon::Index0, dimuon::Pt0, dimuon::Eta0, dimuon::Phi0,
		  dimuon::Index1, dimuon::Pt1, dimuon::Eta1, dimuon::Phi0);

} // namespace o2::aod

#endif // O2_ANALYSIS_SECONDARYVERTEX_H_
