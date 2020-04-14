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

#include <TMath.h>

namespace o2::aod
{
namespace derivedmuon
{
DECLARE_SOA_COLUMN(InverseBendingMomentum, inverseBendingMomentum, float, "fInverseBendingMomentum");
DECLARE_SOA_COLUMN(ThetaX, thetaX, float, "fThetaX");
DECLARE_SOA_COLUMN(ThetaY, thetaY, float, "fThetaY");
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float inverseBendingMomentum, float thetaX, float thetaY) {
     float nonBendingSlope = TMath::Tan(thetaX);
     float bendingSlope = TMath::Tan(thetaY);
     float pYZ = (inverseBendingMomentum != 0.) ? TMath::Abs(1. / inverseBendingMomentum) : - FLT_MAX;
     float pZ = -pYZ / TMath::Sqrt(1.0 + bendingSlope * bendingSlope);
     return pZ * nonBendingSlope; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float inverseBendingMomentum, float thetaY) {
     float bendingSlope = TMath::Tan(thetaY);
     float pYZ = (inverseBendingMomentum != 0.) ? TMath::Abs(1. / inverseBendingMomentum) : - FLT_MAX;
     float pZ = -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);
     return pZ * bendingSlope;
   });
 DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float inverseBendingMomentum, float thetaY) {
     float bendingSlope = TMath::Tan(thetaY);
     float pYZ = (inverseBendingMomentum != 0.) ? TMath::Abs(1. / inverseBendingMomentum) : - FLT_MAX;
     return -pYZ / TMath::Sqrt(1.0 + bendingSlope*bendingSlope);
   });
 DECLARE_SOA_COLUMN(Phi, phi, [](float px, float py) {
     return TMath::AcosH(py/TMath::Sqrt(px*px+py*py));
   });
 DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) {
     return TMath::AtanH(pz/TMath::Sqrt(px*px+py*py+z*pz));
   });
  DECLARE_SOA_COLUMN(Pt, pt, [](float px, float py) {
     return TMath::Sqrt(px*px+py*py);
 });
 
}
 
 
DECLARE_SOA_TABLE(DerivedMuons, "AOD", "DERIVEDMUON",
		  derivedmuon::InverseBendingMomentum,
		  derivedmuon::ThetaX, derivedmuon::ThetaY,
                  derivedmuon::Px<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaX, derivedmuon::ThetaY>,
		  derivedmuon::Py<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaY>,
		  derivedmuon::Pz<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaY>,
		  derivedmuon::Eta<derivedmuon::Px, derivedmuon::Py, dimuon::Pz>,
                  derivedmuon::Phi<derivedmuon::Px, derivedmuon::Py>,
                  derivedmuon::Pt<derivedmuon::Px, derivedmuon::Py>
		  );

  
namespace dimuon
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(E, e, float, "fE");
DECLARE_SOA_INDEX_COLUMN(Px, px, float, "fPx");
DECLARE_SOA_INDEX_COLUMN(Py, py, float, "fPy");
DECLARE_SOA_INDEX_COLUMN(Pz, pz, float, "fPz");
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
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) {
    return TMath::AtanH(pz/TMath::Sqrt(px*px+py*py+z*pz));
  });
 DECLARE_SOA_DYNAMIC_COLUMN(Y, y, [](float e, float px, float py, float pz) {
     return TMath::AtanH(pz/TMath::Sqrt(e*e+px*px+py*py+pz*pz));
   });
 DECLARE_SOA_COLUMN(Phi, phi, [](float px, float py) {
     return TMath::AcosH(py/TMath::Sqrt(px*px+py*py));
   });
 DECLARE_SOA_COLUMN(Pt, pt, [](float px, float py) {
     return TMath::Sqrt(px*px+py*py);
 });
DECLARE_SOA_DYNAMIC_COLUMN(Mass, mass, [](float e, float px, float py, float pz) {
    return Sqrt(e*e-px*px-py*py-pz*pz);//numerically, very bad way to calculate mass...
 });

//TODO add:
//MTR/MID information
//option for cluster information (see Z-analysis)
//2ndary vertexing compatibiliy
//dielectron analysis compatibility: same derived table from different inputs
} // namespace dimuon

DECLARE_SOA_TABLE(Dimuons, "AOD", "DIMUON",
                  //o2::soa::Index<>,
		  dimuon::CollisionId,
		  dimuon::E,
		  dimuon::Px, dimuon::Py, dimuon::Pz,
		  dimuon::Index0, dimuon::Pt0, dimuon::Eta0, dimuon::Phi0,
		  dimuon::Index1, dimuon::Pt1, dimuon::Eta1, dimuon::Phi0,
		  dimuon::Eta<dimuon::Px, dimuon::Py, dimuon::Pz>,
                  dimuon::Y<dimuon::Px, dimuon::Py, dimuon::Pz>,
                  dimuon::Phi<dimuon::Px, dimuon::Py>,
                  dimuon::Pt<dimuon::Px, dimuon::Py>,
                  dimuon::Mass<dimuon::E, dimuon::Px, dimuon::Py, dimuon::Pz>
		  );

} // namespace o2::aod

using namespace o2;
using namespace o2::framework;



#endif // O2_ANALYSIS_SECONDARYVERTEX_H_
