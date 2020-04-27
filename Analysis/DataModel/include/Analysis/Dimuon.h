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
namespace derivedmuon
{
DECLARE_SOA_COLUMN(InverseBendingMomentum, inverseBendingMomentum, float, "fInverseBendingMomentum");
DECLARE_SOA_COLUMN(ThetaX, thetaX, float, "fThetaX");
DECLARE_SOA_COLUMN(ThetaY, thetaY, float, "fThetaY");
DECLARE_SOA_DYNAMIC_COLUMN(Charge, charge, [] (float inverseBendingMomentum) {
     if (inverseBendingMomentum > 0) return 1;
     if (inverseBendingMomentum < 0) return -1;
     return 0;
   });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float inverseBendingMomentum, float thetaX, float thetaY) {
     float nonBendingSlope = tanf(thetaX);
     float bendingSlope = tanf(thetaY);
     float pYZ = (inverseBendingMomentum != 0.) ? fabsf(1. / inverseBendingMomentum) : - 100000;//TODO: put proper max resolvable momentum
     float pZ = -pYZ / sqrtf(1.0 + bendingSlope * bendingSlope);
     return pZ * nonBendingSlope; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float inverseBendingMomentum, float thetaY) {
     float bendingSlope = tanf(thetaY);
     float pYZ = (inverseBendingMomentum != 0.) ? fabsf(1. / inverseBendingMomentum) : - 100000;//TODO: see above
     float pZ = -pYZ / sqrtf(1.0 + bendingSlope*bendingSlope);
     return pZ * bendingSlope; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float inverseBendingMomentum, float thetaY) {
     float bendingSlope = tanf(thetaY);
     float pYZ = (inverseBendingMomentum != 0.) ? fabsf(1. / inverseBendingMomentum) : - 100000;//TODO: see above
     return -pYZ / sqrtf(1.0 + bendingSlope*bendingSlope);  });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float thetaX, float thetaY) {
    return acoshf(tanf(thetaX)/tanf(thetaY));//sign gets lost
   });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float inverseBendingMomentum, float thetaX, float thetaY) {
    return atanhf(1.0/sqrtf(1+1./(tanf(thetaX)*tanf(thetaX))+1./(tanf(thetaY)*tanf(thetaY))));
   });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float inverseBendingMomentum, float thetaX, float thetaY) {
    return sqrtf((tanf(thetaX)*tanf(thetaX)+tanf(thetaY)*tanf(thetaY))*(sqrtf(1.0 + tanf(thetaY)*tanf(thetaY)))/(inverseBendingMomentum*inverseBendingMomentum));
    });
}
 
 
DECLARE_SOA_TABLE(DerivedMuons, "AOD", "DERIVEDMUON",
		  derivedmuon::InverseBendingMomentum,
		  derivedmuon::ThetaX, derivedmuon::ThetaY,
		  derivedmuon::Charge<derivedmuon::InverseBendingMomentum>,
                  derivedmuon::Px<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaX, derivedmuon::ThetaY>,
		  derivedmuon::Py<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaY>,
		  derivedmuon::Pz<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaY>,
		  derivedmuon::Phi<derivedmuon::ThetaX, derivedmuon::ThetaY>,
		  derivedmuon::Eta<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaX, derivedmuon::ThetaY>,
                  derivedmuon::Pt<derivedmuon::InverseBendingMomentum, derivedmuon::ThetaX, derivedmuon::ThetaY>
		  );

  
namespace dimuon
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(E, e, float, "fE");
DECLARE_SOA_COLUMN(Px, px, float, "fPx");
DECLARE_SOA_COLUMN(Py, py, float, "fPy");
DECLARE_SOA_COLUMN(Pz, pz, float, "fPz");
//DECLARE_SOA_INDEX_COLUMN_FULL(Index0, index0, int, DerivedMuons, "fIndex0");
DECLARE_SOA_COLUMN(Pt0, pt0, float, "fPt0");
DECLARE_SOA_COLUMN(Eta0, eta0, float, "fEta0");
DECLARE_SOA_COLUMN(Phi0, phi0, float, "fPhi0");
DECLARE_SOA_COLUMN(Charge0, charge0, int, "fCharge0");
//DECLARE_SOA_INDEX_COLUMN_FULL(Index1, index1, int, DerivedMuons, "fIndex1");
DECLARE_SOA_COLUMN(Pt1, pt1, float, "fPt1");
DECLARE_SOA_COLUMN(Eta1, eta1, float, "fEta1");
DECLARE_SOA_COLUMN(Phi1, phi1, float, "fPhi1");
DECLARE_SOA_COLUMN(Charge1, charge1, int, "fCharge1");
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) {
    return atanhf(pz/sqrtf(px*px+py*py+pz*pz));
  });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, [](float e, float px, float py, float pz) {
     return atanhf(pz/sqrtf(e*e+px*px+py*py+pz*pz));
   });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) {
     return acoshf(py/sqrtf(px*px+py*py));
   });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) {
     return sqrtf(px*px+py*py);
 });
DECLARE_SOA_DYNAMIC_COLUMN(Mass, mass, [](float e, float px, float py, float pz) {
    return sqrtf(e*e-px*px-py*py-pz*pz);//numerically, very bad way to calculate mass...
 });

//TODO add:
//MTR/MID information
//option for cluster information (see Z-analysis)
//2ndary vertexing compatibiliy
//dielectron analysis compatibility: same derived table from different inputs
} // namespace dimuon

DECLARE_SOA_TABLE(Dimuons, "AOD", "DIMUON",
		  o2::soa::Index<>,
		  dimuon::CollisionId,
		  dimuon::E,
		  dimuon::Px, dimuon::Py, dimuon::Pz,
		  //dimuon::Index0Id,
		  dimuon::Pt0, dimuon::Eta0, dimuon::Phi0, dimuon::Charge0,
		  //dimuon::Index1Id,
		  dimuon::Pt1, dimuon::Eta1, dimuon::Phi1, dimuon::Charge1,
		  dimuon::Eta<dimuon::Px, dimuon::Py, dimuon::Pz>,
                  dimuon::Y<dimuon::Px, dimuon::Py, dimuon::Pz>,
                  dimuon::Phi<dimuon::Px, dimuon::Py>,
                  dimuon::Pt<dimuon::Px, dimuon::Py>,
                  dimuon::Mass<dimuon::E, dimuon::Px, dimuon::Py, dimuon::Pz>
		  );

 using Dimuon = Dimuons::iterator;
 
} // namespace o2::aod

#endif 
