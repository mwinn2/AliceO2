// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisHelpers.h"

#include <ROOT/RDataFrame.hxx>

using namespace ROOT::RDF;
using namespace o2::framework;

// A dummy workflow which creates a few of the tables proposed by Ruben,
// using ARROW
WorkflowSpec defineDataProcessing(ConfigContext const& specs)
{
  // Workflow definition. A workflow can be one or more DataProcessors
  // each implementing (part of) an analysis. Each DataProcessor has
  // (at least) a name, some Inputs, some Outputs and they get arranged
  // together accordingly.
  WorkflowSpec workflow{
    //  Multiple DataProcessor specs
    // can be provided per workflow
    DataProcessorSpec{
      // The name of my analysis
      "dimuon-analysis",
      Inputs{
        // Dangling inputs of type AOD will be automatically picked up
        // by DPL and an extra reader device will be instanciated to
        // read them. In this particular case the signature
        // AOD/DZEROFLAGGED is associated to Gianmichele's
        // D0 candidates schema. The first string is just a label
        // so that the algorithm can be in principle be reused for different
        // kind of candidates.
	//        InputSpec{ "Muon", "RN2", "MUON" },//tobe modified for Dimuons, still wrong
	//	InputSpec{ "tracks", "RN2", "TRACKPAR" }, //works
	InputSpec{ "tracks", "RN2", "MUON" },
	  },//muons instead of tracks not working (together with input)
      // No outputs for the time being.TODO tuples!
      Outputs{},
      AlgorithmSpec{
        // This is the actual per "message" loop, where a message could
        // be the contents of a file or part of it.
        // FIXME: Too much boilerplate.
        adaptStateless([](InputRecord& inputs) {
	    auto input = inputs.get<TableConsumer>("tracks");//tbc //has to match first parameter of InputSpec
	    //muons not working instead of tracks (in sync)
	    //muons instead of tracks not working
          // This does a single loop on all the candidates in the input message
          // using a simple mask on the cand_type_ML column and does
          // a simple 1D histogram of the filtered entries.
	    auto tracks = o2::analysis::doSingleLoopOn(input);//need to introduce "CombineParticle"

	    auto h1 = tracks.Filter("fInverseBendingMomentum < 1.00000").Histo1D("fInverseBendingMomentum");
          // A lambda function subtracting two quantities. This defines
          // a function "delta" which can be invoked with
          //
          // delta(1,2)
          //
          // and will return 1 - 2.
	    // auto delta = [](float x, float y) { return x - y; };
	    auto massmuon = [](float invBendMom1, float invBendMom2, float thetax1, float thetax2, float thetay1, float thetay2 //, float mass1 = 0.14, float mass2 = 0.14
			       )
	    {
	      float mass1 = 0.14;
	      float mass2 = 0.14;
	      float e1 = sqrt( mass1 * mass1 + 1/ (invBendMom1 * invBendMom1) * ( 1.0 + sin(thetay1)*sin(thetay1) * (1.0 + 1.0 / (sin(thetax1)*sin(thetax1)) ) ));
	      float e2 = sqrt( mass2 * mass2 + 1/ (invBendMom2 * invBendMom2) * ( 1.0 + sin(thetay2)*sin(thetay2) * (1.0 + 1.0 / (sin(thetax2)*sin(thetax2)) ) ));
	      float scalprod = 1.0 / (invBendMom1 * invBendMom2) * (1.0 + sin(thetay1) * sin(thetay1) * (1.0 + 1.0 / sin(thetax1) * sin(thetax2) ) );
	      return sqrt( mass1 * mass1 + mass2 * mass2 + 2.0 * e1 * e2 - 2.0 * scalprod);
	  };
	  //for correlation analysis for Olltraut's paper
          // This does all the combinations for all the candidates which have
          // the same value for cand_evtID_ML (the Event ID).
          // d0_ is the prefix assigned to the outer variable of the double loop.
          // d0bar_ is the prefix assigned to the inner variable of the double loop.
          //
          // The lines below will:
          // * Filter the combinations according to some mask
          // * Define a column delta_phi with the difference in phi between d0 and d0bar phi
          // * Define a column delta_eta with the difference in phi between d0 and d0bar eta
          // * Do two histograms with delta_phi, delta_eta
	  //ok, this can be used
	  auto combinations = o2::analysis::doSelfCombinationsWith(input, "mu", "fID4mu");//ivent index fID4mu"
	  auto massmuons = combinations
	    .Define("massmuons", massmuon, { "mu_fInverseBendingMomentum", "mubar_fInverseBendingMomentum", "mu_fThetaX", "mubar_fThetaX", "mu_fThetaY", "mubar_fThetaY"});//need to think how it knows which variable belongs to which particle

	  //write a function to extract px,py,pz and to add

	  auto h2 = massmuons.Histo1D("massmuons");
	  
          // FIXME: For the moment we hardcode saving the histograms.
          // In reality it should send the results as outputs to a downstream merger
          // process which merges them as wished.
          TFile f("result.root", "RECREATE");
          h1->SetName("InverseBendingMomentum");
          h1->Write();
	  h2->SetName("massmuons");
	  h2->Write();
        }) } }
  };
  return workflow;
}
