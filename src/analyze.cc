#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>

// ROOT includes
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TLorentzVector.h>

// fastjet includes
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/SoftKiller.hh" // In external code, this should be fastjet/contrib/SoftKiller.hh
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh" 

// my classes
#include "RecHit.hh"
#include "RecHitCalibration.hh"
#include "RecHitCollection.hh"
#include "Cluster.hh"
#include "ClusterCollection.hh"
#include "GenParticle.hh"
#include "GenParticleCollection.hh"
#include "Jet.hh"
#include "JetCollection.hh"
#include "JetTree.hh"
#include "JetPlots.hh"
#include "JetEnergyResolutionPlots.hh"
#include "JetAngularResolutionPlots.hh"

bool debug = false;
//bool debug = true;

// -- Test                                                                                                                                                   
TH1D *h1 = new TH1D("h1","",100,0,10000);
TH1D *h2 = new TH1D("h2","",100,0,10000);

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;
//---------------------------------
std::vector<int> mapIDs = {111,211,-211,2212,2112};

//----------------------------------------------------------------------
struct selection {
   double ptmin;
   double ptmax;
   double absetamin;
   double absetamax;
};

//----------------------------------------------------------------------

template <class Collection>
void produceJets(Collection& input_particles, JetCollection & jets, const float& r, const selection cuts, const bool doPuSubtraction = false, const bool doSubstructure = false, const double beta_sd=0.);

template <class Sequence>
void convertJets(Sequence seq, vector<PseudoJet> pseudojets, const float& r, JetCollection& jets, const bool doSubstructure = false, const double beta_sd=0.);

void matchJets(JetCollection& genjets, JetCollection& recojets, float dr);
void matchTrackJets(JetCollection genjets, JetCollection recojets, JetCollection& outputCollection, ClusterCollection clusters, GenParticleCollection tracks, ClusterCollection& outputClusterCollection, float dr, float mDr, float ptThr);

void printJets(JetCollection & jets);
void computePuOffset(RecHitCollection& rechits);
void sumOverCone(JetCollection& newjets, JetCollection& recojets, RecHitCollection& rechits,float dr);

void addHitsToJet(JetCollection& recojets, RecHitCollection& rec,float minR,float dr);

void sumEnergiesAroundJet(JetCollection& recojets, RecHitCollection& input_collection, GenParticleCollection& genPart,float alpha,float dr);
void sumClusterEnergiesAroundJet(JetCollection& recojets, ClusterCollection& input_collection, GenParticleCollection& genPart,float alpha,float dr);
void sumTrackEnergiesAroundJet(JetCollection& recojets, GenParticleCollection& input_collection,GenParticleCollection& genPart,float alpha,float dr);

void setMultiSD(JetCollection& recojets, JetCollection& recojets02, JetCollection& recojets04, JetCollection& recojets08);

//----------------------------------------------------------------------
int main(int argc, char* argv[]){


 // Check the number of parameters
  if (argc < 9) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " [input.root] " << " [output.root] " <<" [Nevts] [DeltaR] [maxEta] " << " [checkJetCone] " << " [beta] "<< " [alpha] "<< std::endl;
    return 1;
  }

  double DeltaR = atof(argv[4]);
  double matchDeltaR = 0.3;
  double EF_alpha = atof(argv[9]);

  // if cone size smaller 0.3, cone size for matching is set to the total cone size
  if ( DeltaR < matchDeltaR)
    matchDeltaR = DeltaR;

  double maxEta = atof(argv[5]);
  int check = atof(argv[7]);
  double DeltaR_pfa_match = 0.1;
  double beta_softDrop = atof(argv[8]);
  std::cout << "Soft-drop parameter beta set to:  " << beta_softDrop << std::endl;

  bool doConeCheck = false;
  bool useSmearedParts = false;
  bool useClusters = false;
  bool useRechits = false;
  
  if (check==1)
    doConeCheck = true;
  else if (check==2) {
    useSmearedParts = true;
    useRechits = false;
    std::cout << "Use smeared generated particles, as tracks. " << std::endl;
  }
  else if (check==3) {
    useSmearedParts = true;
    useClusters = true;
    std::cout << "Use smeared generated particles, as tracks, and clusters for dummy PFA. " << std::endl;
  }
  else {
    useRechits = true;
    useClusters = true;
  }

  // ---   Tree stuff declarations
  TFile *f = new TFile(argv[1]);
  TTree *t = (TTree*)f->Get("events");
  
  vector<Float_t> *rechit_pt        = 0;
  vector<Float_t> *rechit_eta       = 0;
  vector<Float_t> *rechit_phi       = 0;
  vector<Float_t> *rechit_energy    = 0;
  vector<Float_t> *rechit_x         = 0;
  vector<Float_t> *rechit_y         = 0;
  vector<Float_t> *rechit_z         = 0;
  vector<Float_t> *rechit_layer     = 0;
  vector<int> *rechit_detid     = 0;

  vector<Float_t> *cluster_pt        = 0;
  vector<Float_t> *cluster_eta       = 0;
  vector<Float_t> *cluster_phi       = 0;
  vector<Float_t> *cluster_energy    = 0;
  vector<Float_t> *cluster_x         = 0;
  vector<Float_t> *cluster_y         = 0;
  vector<Float_t> *cluster_z         = 0;

  vector<Float_t> *genpart_pt       = 0;
  vector<Float_t> *genpart_eta      = 0;
  vector<Float_t> *genpart_phi      = 0;
  vector<Float_t> *genpart_energy   = 0;
  vector<Float_t> *genpart_status   = 0;
  vector<Float_t> *genpart_pdgid    = 0;

  vector<Float_t> *track_pt       = 0;
  vector<Float_t> *track_eta      = 0;
  vector<Float_t> *track_phi      = 0;
  vector<Float_t> *track_energy   = 0;
  vector<Float_t> *track_status   = 0;
  vector<Float_t> *track_pdgid    = 0;

  t->SetBranchAddress("rechit_eta", &rechit_eta);
  t->SetBranchAddress("rechit_phi", &rechit_phi);
  t->SetBranchAddress("rechit_pt", &rechit_pt);
  t->SetBranchAddress("rechit_energy", &rechit_energy);
  t->SetBranchAddress("rechit_x", &rechit_x);
  t->SetBranchAddress("rechit_y", &rechit_y);
  t->SetBranchAddress("rechit_z", &rechit_z);
  t->SetBranchAddress("rechit_layer", &rechit_layer);
  t->SetBranchAddress("rechit_detid", &rechit_detid);

  t->SetBranchAddress("cluster_eta",    &cluster_eta);
  t->SetBranchAddress("cluster_phi",    &cluster_phi);
  t->SetBranchAddress("cluster_pt",     &cluster_pt);
  t->SetBranchAddress("cluster_energy", &cluster_energy);
  t->SetBranchAddress("cluster_x",      &cluster_x);
  t->SetBranchAddress("cluster_y",      &cluster_y);
  t->SetBranchAddress("cluster_z",      &cluster_z);
  
  t->SetBranchAddress("gen_eta", &genpart_eta);
  t->SetBranchAddress("gen_phi", &genpart_phi);
  t->SetBranchAddress("gen_pt", &genpart_pt);
  t->SetBranchAddress("gen_energy", &genpart_energy);
  t->SetBranchAddress("gen_status", &genpart_status);
  t->SetBranchAddress("gen_pdgid", &genpart_pdgid);

  t->SetBranchAddress("track_eta",    &track_eta);
  t->SetBranchAddress("track_phi",    &track_phi);
  t->SetBranchAddress("track_pt",     &track_pt);
  t->SetBranchAddress("track_energy", &track_energy);
  t->SetBranchAddress("track_status", &track_status);
  t->SetBranchAddress("track_pdgid",  &track_pdgid);

  // declare histograms
  vector<float> ptvals;
  ptvals = {10., 20., 30.,50., 75., 100., 150., 200., 300., 500., 750., 1000., 1500., 2000., 3500., 5000., 7500., 15000.};
  
  JetPlots gen_plots  = JetPlots("gen", ptvals);
  JetPlots reco_plots = JetPlots("reco", ptvals);
  JetPlots track_plots = JetPlots("track", ptvals);

  JetAngularResolutionPlots jet_angular_plots = JetAngularResolutionPlots("angular", ptvals);
  JetEnergyResolutionPlots jet_energy_plots = JetEnergyResolutionPlots("energy", ptvals);

  JetTree reco_tree = JetTree();
  
  // calibration 
  RecHitCalibration recHitCalibration;

  // generic declarations
  TLorentzVector cluster_p4, cluster_pos;
  TLorentzVector genpart_p4, track_p4;

  //read all entries and fill the histograms
  Long64_t nentries = t->GetEntries();
  Long64_t nmax = atoi(argv[3]);
  Int_t nrun = TMath::Min(nentries, nmax);


  //store plots in output file                                                                                                                              
  TFile outfile(argv[2],"RECREATE");

  for (Long64_t i=0;i<nrun;i++) {
    t->GetEntry(i);

    cout<<" ---- processing event : "<<i<<endl;

    // ---  prepare genparts
    GenParticleCollection genparts;
    GenParticleCollection clean_genparts;

    GenParticleCollection smearedparts;
    GenParticleCollection clean_smearedparts;

    TString findGen("gen_pt");
    TBranch* useGen = t->FindBranch(findGen);
    unsigned genpart_size = 0;
    
    if ( useGen ){
      genpart_size = genpart_pt->size();
      std::cout << "Generated particles size: " << genpart_size << std::endl;
      
      for (unsigned i = 0; i < genpart_size; i++) {
	// initialize genpart
	genpart_p4.SetPtEtaPhiE(genpart_pt->at(i), genpart_eta->at(i), genpart_phi->at(i), genpart_energy->at(i));
	genparts.AddGenParticle(genpart_p4, genpart_pdgid->at(i), genpart_status->at(i));
      }  
      for (unsigned i = 0; i < genparts.size(); i++) {
        GenParticle *g = genparts.at(i);
	if (!g->isClusterable()) continue;
        clean_genparts.Add(new GenParticle(*g));
      }
    }
    if ( useSmearedParts ){
      uint smearedpart_size = track_pt->size();
      std::cout << "Smeared particles size: " << smearedpart_size << std::endl;
      for (unsigned i = 0; i < smearedpart_size; i++) {
        // initialize genpart
	track_p4.SetPtEtaPhiE(track_pt->at(i), track_eta->at(i), track_phi->at(i), track_energy->at(i));
	smearedparts.AddGenParticle(track_p4, track_pdgid->at(i), track_status->at(i));
      }
      std::cout << "smeared particles size: " << smearedparts.size() << std::endl;
      for (unsigned i = 0; i < smearedparts.size(); i++) {
	GenParticle *g = smearedparts.at(i);
        if (!g->isClusterable()) continue;
        clean_smearedparts.Add(new GenParticle(*g));
      }
    }
    
    // ---  prepare rechits ----------------------------------------------
    
    TString findRechit("rechit_pt");
    bool haveRechits = t->FindBranch(findRechit);
    
    RecHitCollection rechits; 
    
    if ( haveRechits && !useSmearedParts ){

      unsigned rechit_size = rechit_pt->size(); 
      if (rechit_size!=0){
	useClusters = false;
	std::cout << "rechit size before: " << rechit_size << std::endl;
	
	for (unsigned i = 0; i < rechit_size; i++) {
	  TLorentzVector rechit_p4, rechit_pos;
	  // initialize rechit
	  rechit_p4.SetPtEtaPhiE(rechit_pt->at(i), rechit_eta->at(i), rechit_phi->at(i), rechit_energy->at(i));
	  rechit_pos.SetX(rechit_x->at(i));
	  rechit_pos.SetY(rechit_y->at(i));
	  rechit_pos.SetZ(rechit_z->at(i));
	  
	  if ( std::isnan(rechit_p4.Pz())){
	    std::cout << "rechit : det:           " <<rechit_detid->at(i)<< std::endl;
	    std::cout << "rechit : eta:           " << rechit_eta->at(i) << std::endl;
	    std::cout << "rechit : x,y,z:           " << rechit_x->at(i)<< ", " <<rechit_y->at(i)<<", "<< rechit_z->at(i) << std::endl;
	    std::cout << "rechit : px,py,pz,energy: " << rechit_p4.Px()<< ", " <<rechit_p4.Py()<<", "<< rechit_p4.Pz()<<", "<< rechit_energy->at(i) << std::endl; 	
	  }
	  // select only hits in Calorimeter
	  if ( rechit_detid->at(i) ==5 || rechit_detid->at(i) == 8 ){ //or rechit_detid->at(i) == 9)){
	    rechits.AddRecHit(rechit_p4, rechit_pos, rechit_energy->at(i));
	  }
	}
	std::cout << "rechit size: " << rechits.size() << std::endl;
      }
    }
    // ---  prepare clusters ----------------------------------------------

    ClusterCollection clusters;
    // make sure the clusters are written in file
    TString findCluster("cluster_pt");
    TBranch* br = t->FindBranch(findCluster);
    
    if ( br && useClusters ){ 
      useRechits=false;
      unsigned cluster_size = cluster_pt->size();
      for (unsigned i = 0; i < cluster_size; i++) {
	// clean clusters from either close to 0 energetic and wrong positions
	if  (std::fabs(cluster_eta->at(i))>2)
	  continue;
	// initialize cluster 
	cluster_p4.SetPtEtaPhiE(cluster_pt->at(i), cluster_eta->at(i), cluster_phi->at(i), cluster_energy->at(i));
	cluster_pos.SetXYZT(cluster_x->at(i), cluster_y->at(i), cluster_z->at(i), 0.0);
	clusters.AddCluster(cluster_p4, cluster_pos);
	//	if (debug && i<10){
	//  std::cout << "pt cluster: " << cluster_pt->at(i) << "\n";
	//}
	h1->Fill(cluster_pt->at(i));
      }
      //      if(debug) 
      cout<<"cluster size: "<<clusters.size()<<endl;
    }

    // ---------- Produce jets ------------------------------------------------
    
    // declare jet collections
    JetCollection genjets;
    JetCollection recojets;
    JetCollection trackjets;
    JetCollection newPFAjets, PFAjets;
    ClusterCollection newPFAs;
    JetCollection recojets_02,recojets_04,recojets_08;
    
    // produce jet collections (anti-kT R = 0.4)
    bool doSubstructure  = true;
    bool doPuSubtraction = false;
    bool useMulti = true;

    selection cuts;
    cuts.ptmin  = 2.5;
    cuts.ptmax  = 20000.;
    cuts.absetamin = 0.0;
    cuts.absetamax = maxEta;
      
    produceJets(clean_genparts, genjets, DeltaR, cuts, false,false, beta_softDrop);
    cout<<"Number of gen jets: " << genjets.size() <<endl;

    if ( useRechits ){
      cout<<"Rechits used for jet production. "<<endl;
      produceJets(rechits, recojets, DeltaR, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
      if (useMulti){
	produceJets(rechits, recojets_02, 0.2, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
	produceJets(rechits, recojets_04, 0.4, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
	produceJets(rechits, recojets_08, 0.8, cuts, doPuSubtraction, doSubstructure, beta_softDrop);

	setMultiSD(recojets, recojets_02, recojets_04, recojets_08);
      }
    } 
    if ( useSmearedParts ){
      cout<<"Tracks used for jet production. "<<endl;
      produceJets(clean_smearedparts, trackjets, DeltaR, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
      if (useMulti){
        produceJets(clean_smearedparts, recojets_02, 0.2, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
        produceJets(clean_smearedparts, recojets_04, 0.4, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
        produceJets(clean_smearedparts, recojets_08, 0.8, cuts, doPuSubtraction, doSubstructure, beta_softDrop);

        setMultiSD(trackjets, recojets_02, recojets_04, recojets_08);
      }
    } 
    if ( useClusters ){
      cout<<"Clusters used for jet production. "<<endl;
      produceJets(clusters, recojets, DeltaR, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
      if (useMulti){
        produceJets(clusters, recojets_02, 0.2, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
        produceJets(clusters, recojets_04, 0.4, cuts, doPuSubtraction, doSubstructure, beta_softDrop);
        produceJets(clusters, recojets_08, 0.8, cuts, doPuSubtraction, doSubstructure, beta_softDrop);

        setMultiSD(recojets, recojets_02, recojets_04, recojets_08);
      }
   }
    if (debug) { 
      cout<<" ------  gen jets ------"<<endl;
      printJets(genjets);
      cout<<" ------  reco jets ------ "<<endl;
      printJets(recojets);
      cout<<" ------  track jets ------ "<<endl;
      printJets(trackjets);
    }
    
    // if you check the anti-kt altogithm for basicaly summing up the cells within a cone
    // !!!! NEEDED FOR ANGULAR RESOLUTION 
    if (doConeCheck){
      cout<<" ------ rechits summed around anti-kt jet axis' ------"<<endl;
      
      JetCollection newjets;
      sumOverCone(newjets, recojets, rechits, DeltaR);
      
      // match reco to gen (need this in order to make resolution plots) 
      // matching aroud 0.3 (ATLAS-CONF-2015-037)
      matchJets(genjets, newjets, matchDeltaR);
      //addHitsToJet(newjets,clean_rechits, 0., DeltaR);
      if (debug) { 
	cout<<" ------  jets in cone ------ "<<endl;
	printJets(newjets);
      }
      reco_tree.fill(newjets);
      gen_plots.fill(genjets);  
      reco_plots.fill(newjets);
      jet_energy_plots.fill(newjets);
      jet_angular_plots.fill(newjets);
    }
    // run normally
    else {
      // match reco to gen (need this in order to make resolution plots)                                                                                            
      if ( useRechits ){
	matchJets(genjets, recojets, matchDeltaR);
	cout<<"Rechits summed around jet axis. "<<endl;
	sumEnergiesAroundJet(recojets, rechits, genparts, EF_alpha, DeltaR);
	addHitsToJet(recojets,rechits, 0., DeltaR);
      }
      else{
	if ( useSmearedParts && !useClusters ){
	  matchJets(genjets, trackjets, matchDeltaR);
	  cout<<"Smeared particles summed around jet axis. "<<endl;
	  sumTrackEnergiesAroundJet(trackjets, clean_smearedparts, genparts, EF_alpha, DeltaR);
	}
	if ( useClusters ){
	  matchJets(genjets, recojets, matchDeltaR);
          cout<<"Cluster summed around jet axis. "<<endl;
	  sumClusterEnergiesAroundJet(recojets, clusters, genparts, EF_alpha, DeltaR);
	}
	if (useSmearedParts && useClusters) {
	  // match reco and trackjets
	  matchTrackJets(genjets, recojets, newPFAjets, clusters, smearedparts, newPFAs, DeltaR, DeltaR_pfa_match, 50);
	}
      }
      // fill plots
      gen_plots.fill(genjets);
      reco_plots.fill(recojets);
      track_plots.fill(trackjets);
      // merged track-cluster pfa objects
      if (useSmearedParts && useClusters){
	reco_tree.fill(newPFAjets);
	jet_energy_plots.fill(newPFAjets);
      }
      else if (useSmearedParts && !useClusters){
	reco_tree.fill(trackjets);
	jet_energy_plots.fill(trackjets);
      }
      else{
	reco_tree.fill(recojets);
	jet_energy_plots.fill(recojets);
      }
      jet_angular_plots.fill(recojets);
    }
  } // end event loop
  
    // store plots in output root tree
  reco_tree.write();
  gen_plots.write();
  reco_plots.write();
  jet_energy_plots.write();
  jet_angular_plots.write();
  
  outfile.Close();

  if (debug){
    std::cout << "cluster inputs: " << h1->GetEntries() << "\n";
    std::cout << "jet inputs:     " << h2->GetEntries() << "\n";
  }
  TCanvas can("can","can");
  can.SetLogy();
  h1->SetLineWidth(2);
  h1->SetLineColor(1);
  h1->Draw();
  h2->SetLineColor(2);
  h2->Draw("same");
  can.Print("test_pt_input.pdf");

  return 0;
}

//------------------------------------------------------------------------------------------------------
void computePuOffset(RecHitCollection& rechits) {
   
   const int nLayers = 53;
   const int etaSlices = 5;
   const float etamin = 1.479;
   const float etamax = 3.;
   float step = float((etamax - etamin)/etaSlices);
   
   vector<double> layer_energy_vector[nLayers][etaSlices]; 
   
   for (unsigned i = 0; i < rechits.size(); i++) {
       for (unsigned j = 0; j < etaSlices; j++) {
           float eta1 = etamin + float(j)*step;
           float eta2 = etamin + float(j+1)*step;
           RecHit r = *rechits.at(i);
           if (fabs(r.eta()) < eta2 && fabs(r.eta()) > eta1 )
               layer_energy_vector[r.layer()][j].push_back(r.energy());
       }
   }

   // now compute median for each layer
   double medians[nLayers][etaSlices]; 
   for (unsigned i = 0; i < nLayers ; i++) { 
       for (unsigned j = 0; j < etaSlices; j++) {
           size_t size = layer_energy_vector[i][j].size();
           sort(layer_energy_vector[i][j].begin(), layer_energy_vector[i][j].end());
           double median = 0;
           if (size > 0){
               if (size  % 2 == 0) median = (layer_energy_vector[i][j][size / 2 - 1] + layer_energy_vector[i][j][size / 2]) / 2;
               else median = layer_energy_vector[i][j][size / 2];
           }
           medians[i][j] = median;
       }
   }

   for (unsigned i = 0; i < rechits.size(); i++) {
       for (unsigned j = 0; j < etaSlices; j++) {
           float eta1 = etamin + float(j)*step;
           float eta2 = etamin + float(j+1)*step;
           RecHit r = *rechits.at(i);
           if (fabs(r.eta()) < eta2 && fabs(r.eta()) > eta1 )
              r.setPuOffset(medians[r.layer()][j]);
       }
   }
}

//------------------------------------------------------------------------------------------------------
template <class Collection>
void produceJets(Collection& constituents, JetCollection& jets, const float& r, const selection cuts, const bool doPuSubtraction = false, const bool doSubstructure = false, const double beta_sd) {
   // first convert constituents into fastjet pseudo-jets
   vector <PseudoJet> input_particles;
   for (unsigned i = 0; i < constituents.size(); i++) {
     
     input_particles.push_back( PseudoJet( constituents.at(i)->px(), constituents.at(i)->py(), constituents.at(i)->pz(), constituents.at(i)->energy()));
     if (doSubstructure)
       h2->Fill(PseudoJet( constituents.at(i)->px(), constituents.at(i)->py(), constituents.at(i)->pz(), constituents.at(i)->energy() ).perp());

     if (debug && i<10){
       std::cout << "mass of pseudojet : " << PseudoJet( constituents.at(i)->px(), constituents.at(i)->py(), constituents.at(i)->pz(), constituents.at(i)->energy() ).m() << "\n"; 
       std::cout << "pt of pseudojet : " << PseudoJet( constituents.at(i)->px(), constituents.at(i)->py(), constituents.at(i)->pz(), constituents.at(i)->energy() ).perp() << " \n";
     }
   }
   
   // Initial clustering with anti-kt algorithm
   JetAlgorithm algorithm = antikt_algorithm;
   double jet_rad = r; // jet radius for anti-kt algorithm
   JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
   
   vector<PseudoJet> akjets;
   
   Selector select_eta   = SelectorAbsEtaRange(cuts.absetamin, cuts.absetamax);
   Selector select_pt    = SelectorPtRange(cuts.ptmin, cuts.ptmax);
   Selector select_jets  = select_eta && select_pt;
   
   if(doPuSubtraction) {
       
       // jet area correction parameters, are used only if doPU = true
       float etamin = 0;
       float etamax = cuts.absetamax;
       float spacing = 0.50;
       Selector selector = SelectorAbsRapRange(etamin, etamax);

       AreaDefinition areaDef(active_area, GhostedAreaSpec(selector));
       ClusterSequenceArea clust_seq(input_particles,jetDef, areaDef);
       
       vector<PseudoJet> antikt_jets  = clust_seq.inclusive_jets();

       // now apply PU subtraction
       RectangularGrid grid(-etamax, etamax, spacing, spacing, selector);
       GridMedianBackgroundEstimator gmbge(grid);
       gmbge.set_particles(input_particles);
       Subtractor subtractor(&gmbge);
       akjets = subtractor(antikt_jets);
       
       // apply cuts
       akjets = select_jets(sorted_by_pt(akjets));
       
       // eventually apply substructure and store in custom dataformat
       convertJets(clust_seq, akjets, r, jets, doSubstructure, beta_sd);
   }
   else {
     
     ClusterSequence clust_seq(input_particles,jetDef);
     akjets  = sorted_by_pt(clust_seq.inclusive_jets());

     // apply cuts
     akjets = select_jets(akjets);
     
     // eventually apply substructure and store in custom dataformat
     convertJets(clust_seq, akjets, r, jets, doSubstructure, beta_sd);
   }
   input_particles.clear();
   akjets.clear();
}

//------------------------------------------------------------------------------------------------------
template <class Sequence>
void convertJets(Sequence seq, vector<PseudoJet> pseudojets, const float& r, JetCollection& jets, const bool doSubstructure = false, const double beta_sd){

   TLorentzVector p4;
   for (unsigned j = 0; j < pseudojets.size() ; j++) { 

       // get the jet for analysis
       PseudoJet this_jet = pseudojets[j];

       p4.SetPtEtaPhiM(this_jet.pt(), this_jet.eta(), this_jet.phi(), std::max(this_jet.m(),0.0));
       Jet jet(p4);

       if (doSubstructure) {
       
           // N-subjettiness
           double beta = 1.0;

           Nsubjettiness         nSub1_beta1(1,   OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           Nsubjettiness         nSub2_beta1(2,   OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           Nsubjettiness         nSub3_beta1(3,   OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           NsubjettinessRatio    nSub21_beta1(2,1, OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           NsubjettinessRatio    nSub31_beta1(3,1, OnePass_KT_Axes(), NormalizedMeasure(beta, r));
           NsubjettinessRatio    nSub32_beta1(3,2, OnePass_KT_Axes(), NormalizedMeasure(beta, r));

           jet.setTau1(nSub1_beta1(this_jet));
           jet.setTau2(nSub2_beta1(this_jet));
           jet.setTau3(nSub3_beta1(this_jet));
           jet.setTau21(nSub21_beta1(this_jet));
           jet.setTau31(nSub31_beta1(this_jet));
           jet.setTau32(nSub32_beta1(this_jet));

           // soft drop
           beta = beta_sd;
           double zcut = 0.1;
           SoftDrop softDrop(beta,zcut);
           PseudoJet softdrop_jet = softDrop(this_jet);
           p4.SetPtEtaPhiM(softdrop_jet.pt(), softdrop_jet.eta(), softdrop_jet.phi(), std::max(softdrop_jet.m(),0.0));
           jet.setSDmass(p4.M());

           // store jet in collection
       }
       jets.Add(new Jet(jet));
   }   
}

//------------------------------------------------------------------------------------------------------                                                                     
 void matchJets(JetCollection& genjets, JetCollection& recojets, float dr){

  for (unsigned i = 0; i < genjets.size() ; i++) {
    float dr0 = 999.;
    Jet *gj = genjets.at(i);
    // will be best matching recojet                                                                                                                                                
    Jet *rj0;
    for (unsigned j = 0; j < recojets.size() ; j++) {
      Jet *rj = recojets.at(j);
      float dr_gr = rj->p4().DeltaR(gj->p4());
      if( dr_gr < dr0 ){
	rj0 = rj;
	dr0 = dr_gr;
      }
    }
    // assign genjet ref. to best matching recojet (and vice versa)                                                                                                                 
    if (dr0 < dr){
      rj0->setRef(gj);
      gj->setRef(rj0);
    }
  }
}

//------------------------------------------------------------------------------------------------------
void setMultiSD(JetCollection& recojets, JetCollection& recojets02, JetCollection& recojets04, JetCollection& recojets08){
  // match jets with 3 three reco jets.
  
  // find jets belonging to same gen jet
  matchJets(recojets02, recojets, 0.2);
  for (unsigned i = 0; i < recojets.size() ; i++) {
    Jet* rj = recojets.at(i);
    if (rj->ref()){
      Jet* rj1 = rj->ref();
      //std::cout << "Fill massSD for 0.2 : "<< rj1->massSD() << std::endl;
      recojets.at(i)->setSDmass_02(rj1->massSD());
    }
  }
  matchJets(recojets04, recojets, 0.4);
  for (unsigned i = 0; i < recojets.size() ; i++) {
    Jet* rj = recojets.at(i);
    if (rj->ref()){
      Jet* rj2 = rj->ref();
      //      std::cout << "Fill massSD for 0.4 : "<< rj2->massSD() << std::endl;
      recojets.at(i)->setSDmass_04(rj2->massSD());
    }
  }
  matchJets(recojets08, recojets, 0.8);
  for (unsigned i = 0; i < recojets.size() ; i++) {
    Jet* rj = recojets.at(i);
    if (rj->ref()){
      Jet* rj3 = rj->ref();
      // std::cout << "Fill massSD for 0.8 : "<< rj3->massSD() << std::endl;
      recojets.at(i)->setSDmass_08(rj3->massSD());
    }
  }
}

//------------------------------------------------------------------------------------------------------
void matchTrackJets(JetCollection genjets, JetCollection recojets, JetCollection& outputCollection, ClusterCollection clusters, GenParticleCollection tracks, ClusterCollection& outputClusterCollection, float dr, float mDr, float ptThr){
  
  for (unsigned i = 0; i < genjets.size() ; i++) { 
    Jet *gj = genjets.at(i);
    // get ref jet build from calo clusters
    Jet *rj = gj->ref();
 
    float eta=0.;
    float phi=0.;
    float energy=0.;
    float x=0.;
    float y=0.;
    float z=0.;

    if (rj){
      // new jet consisting of clusters from neutrals and low pt tracks
      
      // fill vectors of tracks and clusters per jet
      std::vector<std::pair<GenParticle*,float> > includedTracks;
      std::vector<Cluster*> includedClusters;
      std::vector<Cluster*> matched;

      std::map<GenParticle*, std::vector<Cluster*> > matchMap;
      TLorentzVector new_p4(0.,0.,0.,0.);

      // collect the tracks/cluster of the jet:
      for (uint t = 0; t < tracks.size() ; t++) {
	GenParticle *track = tracks.at(t);
	float dr_gr = track->p4().DeltaR(rj->p4());
	if( dr_gr < dr ){
	  includedTracks.push_back(std::make_pair(track, track->p4().Pt()));
	}
      } 
      for (uint c = 0; c < clusters.size() ; c++) {
	Cluster *cluster = clusters.at(c);
	float dr_gr = cluster->p4().DeltaR(rj->p4());
	if( dr_gr < dr && std::find(includedClusters.begin(), includedClusters.end(), cluster) == includedClusters.end() ){
	  includedClusters.push_back(cluster);
	}
      }
      uint allCl = includedClusters.size();
      uint allTr = includedTracks.size();

      std::cout << "Found "<< allCl <<" cluster." << std::endl;
      // sort tracks by pt
      std::sort(includedTracks.begin(), includedTracks.end(), [](auto &left, auto &right) {
	  return left.second > right.second;   
	});
      // match tracks and clusters:
      for (auto t : includedTracks){
	GenParticle* charged = t.first;
	bool matchedSuccess=false;
	for (auto c : includedClusters){
	  Cluster *cl = c;
	  float drTC = charged->p4().DeltaR(cl->p4());
	  
	  // if cluster was matched with track already, continue
	  if ( std::find(matched.begin(), matched.end(), cl) != matched.end() )
	    continue;
	  
	  // find match of track and cluster, if within DeltaR
	  if ( drTC < mDr) {
	    if (debug){
	      std::cout << "Found matching track-cluster." << std::endl;}
	    matchMap[charged].push_back(cl);
	    matched.push_back(cl);
	    includedClusters.erase(std::remove(includedClusters.begin(), includedClusters.end(), cl), includedClusters.end());
	    matchedSuccess=true;
	  }
	}
	if (matchedSuccess)
	  includedTracks.erase(std::remove(includedTracks.begin(), includedTracks.end(), t), includedTracks.end());
      }
      
      int counter=0;
      for (auto t : matchMap){
	GenParticle* charged = t.first;
	if (debug){
	  std::cout << "Track with pt:       "<< charged->p4().Pt() << std::endl;
	}
	if ( charged->p4().Pt() < ptThr ){
	  //	  new_p4+=charged->p4();
	  x += charged->px()*charged->energy();
	  y += charged->py()*charged->energy();
	  z += charged->pz()*charged->energy();
	  energy += charged->energy();
	  eta += charged->eta() * charged->energy();
	  phi += charged->phi() * charged->energy();
	  TLorentzVector p4;
	  p4.SetPtEtaPhiE(charged->pt(), charged->eta(), charged->phi(), charged->energy());
	  TLorentzVector pos(charged->px(), charged->py(), charged->pz(), charged->energy());
	  outputClusterCollection.Add(new Cluster(p4, pos));
	  if (debug){
	    std::cout << "gets added to jet.       " << std::endl;
	  }
	}
	else{
	  if (debug){
	    std::cout << "gets the connected " << t.second.size() << " clusters added to jet.  " << std::endl;
	  }
	  for (auto cl : t.second){
	    //  new_p4+=cl->p4();
	    x += cl->px()*cl->energy();
	    y += cl->py()*cl->energy();
	    z += cl->pz()*cl->energy();
	    energy += cl->energy();
	    eta += cl->eta() * cl->energy();
	    phi += cl->phi() * cl->energy();
	    outputClusterCollection.Add(new Cluster(cl->p4(), cl->pos()));
	    
	    counter++;
	  }
	}
      }
      
      if ((matched.size()+includedClusters.size())!=allCl){
	std::cout << "Matched       "<< matched.size() <<" cluster." << std::endl;
	std::cout << "Not matched   "<< includedClusters.size() <<" cluster." << std::endl;
	std::cout << "Originally    "<< allCl <<" cluster." << std::endl;
      }

      if ((matchMap.size()+includedTracks.size())!=allTr){
	std::cout << "Matched       "<< matchMap.size() <<" tracks." << std::endl;
	std::cout << "Not matched   "<< includedTracks.size() <<" tracks." << std::endl;
	std::cout << "Originally    "<< allTr <<" cluster." << std::endl;
      }

      // remaining clusters are counted as neutral and added to cluster
      std::cout << "Left over clusters:    "<< includedClusters.size() <<" are added as neutrals." << std::endl;
      std::cout << "Left over tracks:    "<< includedTracks.size() <<" are added as low pt. " << std::endl;

      for (auto left : includedClusters){
	//	new_p4+=left->p4();
	x += left->px()*left->energy();
        y += left->py()*left->energy();
        z += left->pz()*left->energy();
        energy += left->energy();
        eta += left->eta() * left->energy();
        phi += left->phi() * left->energy();
	outputClusterCollection.Add(new Cluster(left->p4(), left->pos()));

      }
      for (auto left : includedTracks){
        //new_p4+=left.first->p4();
	x += left.first->px()*left.first->energy();
	y += left.first->py()*left.first->energy();
	z += left.first->pz()*left.first->energy();
	energy += left.first->energy();
	eta += left.first->eta() * left.first->energy();
	phi += left.first->phi() * left.first->energy();
	TLorentzVector p4;
	p4.SetPtEtaPhiE(left.first->pt(), left.first->eta(), left.first->phi(), left.first->energy());
	TLorentzVector pos(left.first->px(), left.first->py(), left.first->pz(), left.first->energy());
	outputClusterCollection.Add(new Cluster(p4, pos));
	if(debug){
	  std::cout << "track pt: " << left.first->p4().Pt() << std::endl;
	}
      }
      
      TVector3 vec;
      vec.SetXYZ(x/energy,y/energy,z/energy);
      new_p4.SetPtEtaPhiE(energy*vec.Unit().Perp(), eta/energy, phi/energy,energy);

      Jet theJet(new_p4);
      theJet.setTau1(rj->tau1());
      theJet.setTau2(rj->tau2());
      theJet.setTau3(rj->tau3());
      theJet.setTau21(rj->tau21());
      theJet.setTau32(rj->tau32());
      theJet.setEf1(rj->ef1());
      theJet.setEf2(rj->ef2());
      theJet.setEf3(rj->ef3());
      theJet.setEf4(rj->ef4());
      theJet.setEf5(rj->ef5());
      theJet.setEnergy_clusters(rj->energy_clusters());
      theJet.setEnergy_gens(rj->energy_gens());
      theJet.setEnergy_hadrons(rj->energy_hadrons());
      theJet.setNum_tracks(allTr);
      theJet.setNum_pfas(includedTracks.size());
      theJet.setNum_neutrals(includedClusters.size());
      theJet.setNum_clusters(rj->num_clusters());
      theJet.setNum_gens(rj->num_gens());
      theJet.setNum_hadrons(rj->num_hadrons());
      theJet.setRef(gj);

      std::cout << "True jet pt: " << gj->p4().Pt() <<  ", old jet pt: " <<  rj->p4().Pt() << ", and new jet pt: " << new_p4.Pt() << std::endl;
      std::cout << "True jet energy: " << gj->p4().E() <<  ", old jet energy: " <<  rj->p4().E() << ", and new jet energy: " << new_p4.E() << std::endl;

      // clear vectors/maps
      includedTracks.clear();
      includedClusters.clear();
      matched.clear();
      matchMap.clear();
      
      outputCollection.Add(new Jet(theJet));
    }
  }
}

//------------------------------------------------------------------------------------------------------
void printJets(JetCollection& jets){
    cout<<" -- Print Jet collection -- "<<endl; 
    for (unsigned j = 0; j < jets.size() ; j++){
       jets.at(j)->print();
       if(jets.at(j)->ref()) {
           (jets.at(j)->ref())->print();
       }
    }
}

//------------------------------------------------------------------------------------------------------
void sumOverCone(JetCollection& newjets, JetCollection& recojets, RecHitCollection& recHits,float dr){
 
 for (unsigned i = 0; i < recojets.size() ; i++) { 
   Jet *rj = recojets.at(i);
   float eta;
   float phi;
   float energy;
   float x;
   float y;
   float z;
   TLorentzVector p4; // use SetPtEtaPhiM

   // sums up all reHits aroud DeltaR=0.4
   for (unsigned j = 0; j < recHits.size() ; j++) { 
     RecHit *hit = recHits.at(j);
     float dr_gh = hit->p4().DeltaR(rj->p4());
     if( dr_gh < dr ){
       // Add hit to new jet
       // Calculate transverse momentum
       x += hit->pos().X() * hit->energy();
       y += hit->pos().Y() * hit->energy();
       z += hit->pos().Z() * hit->energy();
       energy += hit->energy();
       eta += hit->pos().Eta() * hit->energy();
       phi += hit->pos().Phi() * hit->energy();
     }
   }
   TVector3 vec;
   vec.SetXYZ(x/energy,y/energy,z/energy);
   p4.SetPtEtaPhiE(energy*vec.Unit().Perp(), eta/energy, phi/energy,energy);
   Jet newjet(p4);
   newjets.Add(new Jet(newjet));
 }
}

//------------------------------------------------------------------------------------------------------
void addHitsToJet(JetCollection& recojets, RecHitCollection& rec, float minR, float dr){
  
  float maxR = minR + dr;

  for (unsigned i = 0; i < recojets.size() ; i++) {
    Jet *rj = recojets.at(i);
    
    vector<RecHit*> hits;
    // sums up all reHits aroud DeltaR                          
    for (unsigned j = 0; j < rec.size() ; j++) {
      auto hit = rec.at(j);
      float dr_gh = hit->p4().DeltaR(rj->p4());
      if( dr_gh > minR && dr_gh < maxR ){
	// Add hit to new jet
	hits.push_back(hit);
      }
    }
    // if (typeid(hits).name() == "RecHitCollection"){
    // rj->setHits(hits);
    //}
  }
}
//------------------------------------------------------------------------------------------------------                                          
void sumEnergiesAroundJet(JetCollection& recojets, RecHitCollection& inputCollection, GenParticleCollection& genParts, float alpha, float dr){

  // Five rings around jet axis with max = 5*alpha
  // Alpha can effectively be maximal 0.4/5 = 0.08
  // Use 0.05 as default!
  // And add the reference to the true generated particles within the jet.
  uint n = 5;
  std::vector<RecHit*> vecHits;

  for (uint i = 0; i < recojets.size() ; i++) {
    Jet *rj = recojets.at(i);
    std::vector<GenParticle*> particles;
    float energies[n]={0.};
    if (debug){
      std::cout << "Particle collection of size   : " << genParts.size() << endl;
    }

    // sums up all energies within rings
    for (uint j = 0; j < inputCollection.size() ; j++) {
      auto hit = inputCollection.at(j);
      float dr_gh = hit->p4().DeltaR(rj->p4());
      
      for (uint k=1; k<=n; k++){
        float min = alpha*(k-1)/n;
        float max = alpha*k/n;
        if( dr_gh >= min && dr_gh < max ){
          energies[k-1] += hit->pt();
	  vecHits.push_back(hit);
	}
      }
    }  // find all generated particles within the cone
    for (uint k = 0; k < genParts.size(); k++){
      auto part = genParts.at(k);
      float dr_gh = part->p4().DeltaR(rj->p4());
      if (dr_gh<dr)
	// add particle to jet          
	particles.push_back(part);
    }
    //rj->setParticles(particles);
    //rj->setHits(vecHits);

    double total = energies[0] + energies[1] + energies[2] + energies[3] + energies[4];
    //if (debug){
    //  std::cout << "Test ring energy n=1: " << energies[0] << endl;
    //  std::cout << "Test ring energy n=2: " << energies[1] << endl;
    //  std::cout << "Test ring energy n=3: " << energies[2] << endl;
    //  std::cout << "Test ring energy n=4: " << energies[3] << endl;
    //  std::cout << "Test ring energy n=5: " << energies[4] << endl;
    //}
    rj->setPt_005(total);
    rj->setEf1(energies[0]/total);
    rj->setEf2(energies[1]/total);
    rj->setEf3(energies[2]/total);
    rj->setEf4(energies[3]/total);
    rj->setEf5(energies[4]/total);

    vecHits.clear();
  }
}
//------------------------------------------------------------------------------------------------------
void sumTrackEnergiesAroundJet(JetCollection& recojets, GenParticleCollection& constituents2, GenParticleCollection& genParts, float alpha, float dr){

  // Five rings around jet axis with max = 5*alpha
  // Alpha can effectively be maximal 0.4/5 = 0.08
  // Use 0.05 as default!
  // And add the reference to the true generated particles within the jet.
  uint n = 5;

  for (unsigned i = 0; i < recojets.size() ; i++) {
    Jet *rj = recojets.at(i);
    float energies[n]={0.};
    
    if (debug){
      std::cout << "Particle collection of size   : " << genParts.size() << endl;
    }

    // sums up all energies within rings
    for (unsigned j = 0; j < constituents2.size() ; j++) {
      auto hit = constituents2.at(j);
      float dr_gh = hit->p4().DeltaR(rj->p4());

      for (uint k=1; k<=n; k++){
        float min = alpha*(k-1)/n;
        float max = alpha*k/n;
        if( dr_gh >= min && dr_gh < max ){
          energies[k-1] += hit->pt();
        }
      }
    }
   
    double total = energies[0] + energies[1] + energies[2] + energies[3] + energies[4];
    if (debug){
      std::cout << "Test ring energy n=1: " << energies[0] << endl;
      std::cout << "Test ring energy n=2: " << energies[1] << endl;
      std::cout << "Test ring energy n=3: " << energies[2] << endl;
      std::cout << "Test ring energy n=4: " << energies[3] << endl;
      std::cout << "Test ring energy n=5: " << energies[4] << endl;
    }
    rj->setPt_005(total);
    rj->setEf1(energies[0]/total);
    rj->setEf2(energies[1]/total);
    rj->setEf3(energies[2]/total);
    rj->setEf4(energies[3]/total);
    rj->setEf5(energies[4]/total);
    
  }
}

//------------------------------------------------------------------------------------------------------
void sumClusterEnergiesAroundJet(JetCollection& recojets, ClusterCollection& inputCollection, GenParticleCollection& genParts, float alpha, float dr){

  // Five rings around jet axis with max = 5*alpha
  // Alpha can effectively be maximal 0.4/5 = 0.08
  // Use 0.05 as default!
  // And add the reference to the true generated particles within the jet.
  uint n = 5;
  
  if (debug){
    std::cout << "Cluster collection of size    : " << inputCollection.size() << endl;
    std::cout << "Particle collection of size   : " << genParts.size() << endl;
  }

  for (unsigned i = 0; i < recojets.size(); i++) {
    Jet *rj = recojets.at(i);
    float energy_cluster=0., energy_gens=0., energy_hadrons=0.; //, energy_hadrons_1GeV=0., energy_gens_1GeV=0.;
    int num_cluster=0, num_gens=0, num_hadrons=0; //, num_hadrons_1GeV=0, num_gens_1GeV=0;

    float energies[n]={0.};
    
    // sums up all energies within rings
    for (unsigned j = 0; j < inputCollection.size() ; j++) {
      auto hit = inputCollection.at(j);
      float dr_gh = hit->p4().DeltaR(rj->p4());
      if (dr_gh<dr){
	num_cluster++;
	energy_cluster += hit->pt();
      }
      for (uint k=1; k<=n; k++){
	float min = alpha*(k-1)/n;
	float max = alpha*k/n;
	
	if( dr_gh >= min && dr_gh < max ){
	  energies[k-1] += hit->pt();
	}
      }
    }
    // find all generated particles within the cone
    for (uint k = 0; k < genParts.size(); k++){
      auto part = genParts.at(k);
      float dr_gh = part->p4().DeltaR(rj->p4());
      if (dr_gh<dr){
	// add particle to jet 
	energy_gens += part->pt();
	num_gens++;
	auto it = find(mapIDs.begin(), mapIDs.end(), part->pdgid());
	if (it!=mapIDs.end()){
	  //its a hadron
	  energy_hadrons += part->pt();
	  num_hadrons++;
//	  if (part->energy()>1){
//	    energy_hadrons_1GeV += part->energy();
//	    num_hadrons_1GeV++;
//	  }
	}
//	if (part->energy()>1){
//	  energy_gens_1GeV += part->energy();
//	  num_gens_1GeV++;
//	}
      }
    }
    rj->setEnergy_clusters(energy_cluster);
    rj->setEnergy_gens(energy_gens);
    //    rj->setEnergy_gens_leq_1GeV(energy_gens_1GeV);
    rj->setEnergy_hadrons(energy_hadrons);
    //rj->setEnergy_hadrons_leq_1GeV(energy_hadrons_1GeV);

    rj->setNum_clusters(num_cluster);
    rj->setNum_gens(num_gens);
    //rj->setNum_gens_leq_1GeV(num_gens_1GeV);
    rj->setNum_hadrons(num_hadrons);
    //rj->setNum_hadrons_leq_1GeV(num_hadrons_1GeV);

    double total = energies[0] + energies[1] + energies[2] + energies[3] + energies[4];
    if (debug){
      std::cout << "Test ring energy n=1: " << energies[0] << endl;
      std::cout << "Test ring energy n=2: " << energies[1] << endl;
      std::cout << "Test ring energy n=3: " << energies[2] << endl;
      std::cout << "Test ring energy n=4: " << energies[3] << endl;
      std::cout << "Test ring energy n=5: " << energies[4] << endl;
    }
    rj->setPt_005(total);
    rj->setEf1(energies[0]/total);
    rj->setEf2(energies[1]/total);
    rj->setEf3(energies[2]/total);
    rj->setEf4(energies[3]/total);
    rj->setEf5(energies[4]/total);
  }
}
