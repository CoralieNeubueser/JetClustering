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

// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
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

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;


//----------------------------------------------------------------------
struct selection {
   double ptmin;
   double ptmax;
   double absetamin;
   double absetamax;
};

//----------------------------------------------------------------------

template <class Collection>
void produceJets(Collection& input_particles, JetCollection & jets, const float& r, const selection cuts, const bool doPuSubtraction = false, const bool doSubstructure = false);

template <class Sequence>
void convertJets(Sequence seq, vector<PseudoJet> pseudojets, const float& r, JetCollection& jets, const bool doSubstructure = false);
void matchJets(JetCollection& genjets, JetCollection& recojets, float dr);
void printJets(JetCollection & jets);
void computePuOffset(RecHitCollection& rechits);
void sumOverCone(JetCollection& newjets, JetCollection& recojets, RecHitCollection& rechits,float dr);
void addHitsToJet(JetCollection& recojets, RecHitCollection& rechits,float minR,float dr);
void sumEnergiesAroundJet(JetCollection& recojets, RecHitCollection& rechits,float alpha);

//----------------------------------------------------------------------
int main(int argc, char* argv[]){


 // Check the number of parameters
  if (argc < 6) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " [input.root] " << " [output.root] " <<" [Nevts] [DeltaR] [maxEta] " << " [checkJetCone] " << std::endl;
    return 1;
  }

  double DeltaR = atof(argv[4]);
  double maxEta = atof(argv[5]);
   
  TString doConeCheck = argv[6];
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
  vector<Float_t> *rechit_thickness = 0;
  vector<Float_t> *rechit_layer     = 0;

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

  t->SetBranchAddress("rechit_eta", &rechit_eta);
  t->SetBranchAddress("rechit_phi", &rechit_phi);
  t->SetBranchAddress("rechit_pt", &rechit_pt);
  t->SetBranchAddress("rechit_energy", &rechit_energy);
  t->SetBranchAddress("rechit_x", &rechit_x);
  t->SetBranchAddress("rechit_y", &rechit_y);
  t->SetBranchAddress("rechit_z", &rechit_z);
  t->SetBranchAddress("rechit_layer", &rechit_layer);

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

  // declare histograms
  vector<float> ptvals;
  ptvals = {10., 20., 30.,50., 75., 100., 150., 200., 300., 500., 750., 1000., 1500., 2000., 3500., 5000., 7500., 15000.};
  
  JetPlots gen_plots  = JetPlots("gen", ptvals);
  JetPlots reco_plots = JetPlots("reco", ptvals);

  JetAngularResolutionPlots jet_angular_plots = JetAngularResolutionPlots("angular", ptvals);
  JetEnergyResolutionPlots jet_energy_plots = JetEnergyResolutionPlots("energy", ptvals);

  // write jets into tree
  JetTree reco_tree = JetTree();

  // calibration 
  RecHitCalibration recHitCalibration;

  // generic declarations
  TLorentzVector rechit_p4, rechit_pos;
  TLorentzVector cluster_p4, cluster_pos;
  TLorentzVector genpart_p4;

  //read all entries and fill the histograms
  Long64_t nentries = t->GetEntries();
  Long64_t nmax = atoi(argv[3]);
  Int_t nrun = TMath::Min(nentries, nmax);

  for (Long64_t i=0;i<nrun;i++) {
    t->GetEntry(i);

    cout<<" ---- processing event : "<<i<<endl;

    // ---  prepare genparts
    GenParticleCollection genparts;
    unsigned genpart_size = genpart_pt->size();

    for (unsigned i = 0; i < genpart_size; i++) {
       // initialize genpart
       genpart_p4.SetPtEtaPhiE(genpart_pt->at(i), genpart_eta->at(i), genpart_phi->at(i), genpart_energy->at(i));
       genparts.AddGenParticle(genpart_p4, genpart_pdgid->at(i), genpart_status->at(i));
    }

    GenParticleCollection clean_genparts;
    for (unsigned i = 0; i < genparts.size(); i++) {
        GenParticle *g = genparts.at(i);
        if (!g->isClusterable()) continue;
        clean_genparts.Add(new GenParticle(*g));
    }

    // ---  prepare rechits ----------------------------------------------
    
    RecHitCollection rechits;
    unsigned rechit_size = rechit_pt->size();
    for (unsigned i = 0; i < rechit_size; i++) {
      // initialize rechit
      rechit_p4.SetPtEtaPhiE(rechit_pt->at(i), rechit_eta->at(i), rechit_phi->at(i), rechit_energy->at(i));
      rechit_pos.SetXYZT(rechit_x->at(i), rechit_y->at(i), rechit_z->at(i), 0.0);
      rechits.AddRecHit(rechit_p4, rechit_pos, rechit_layer->at(i));
    }
    
    if(debug) cout<<"rechit size: "<<rechits.size()<<endl;
    RecHitCollection clean_rechits;
    for (unsigned i = 0; i < rechits.size(); i++) {
      RecHit *r = rechits.at(i);
      clean_rechits.Add(new RecHit(*r));
    }
    if(debug) cout<<"clean rechit size: "<<clean_rechits.size()<<endl;
    
    // ---  prepare clusters ----------------------------------------------                                                                                                                                                                    
    ClusterCollection clusters;
    // make sure the clusters are written in file
    TString findCluster("cluster_pt");
    TBranch* br = t->FindBranch(findCluster);
    if (br){ 
      unsigned cluster_size = cluster_pt->size();
      for (unsigned i = 0; i < cluster_size; i++) {
	// initialize cluster 
	cluster_p4.SetPtEtaPhiE(cluster_pt->at(i), cluster_eta->at(i), cluster_phi->at(i), cluster_energy->at(i));
	cluster_pos.SetXYZT(cluster_x->at(i), cluster_y->at(i), cluster_z->at(i), 0.0);
	clusters.AddCluster(cluster_p4, cluster_pos);
      }
      if(debug) cout<<"cluster size: "<<clusters.size()<<endl;
    }
    // ---------- Produce jets ------------------------------------------------
    
    // declare jet collections
    JetCollection genjets;
    JetCollection recojets;
    
    // produce jet collections (anti-kT R = 0.4)
    bool doSubstructure  = true;
    bool doPuSubtraction = false;
    
    selection cuts;
    cuts.ptmin  = 2.5;
    cuts.ptmax  = 20000.;
    cuts.absetamin = 0.0;
    cuts.absetamax = maxEta;
    
    produceJets(clean_genparts, genjets, DeltaR, cuts, false, doSubstructure);
    if (clean_rechits.size() > 0)
      produceJets(clean_rechits, recojets, DeltaR, cuts, doPuSubtraction, doSubstructure);
    else if (clusters.size() > 0)
      produceJets(clusters, recojets, DeltaR, cuts, doPuSubtraction, doSubstructure);

    if (debug) { 
        cout<<" ------  gen jets ------"<<endl;
        printJets(genjets);
        cout<<" ------  reco jets ------ "<<endl;
        printJets(recojets);
    }

    // if you check the anti-kt altogithm for basicaly summing up the cells within a cone
    // !!!! NEEDED FOR ANGULAR RESOLUTION 
    if (doConeCheck == "1"){
      cout<<" ------ rechits summed around anti-kt jet axis' ------"<<endl;
	
      JetCollection newjets;
      sumOverCone(newjets, recojets, clean_rechits, DeltaR);
    
      // match reco to gen (need this in order to make resolution plots) 
      // matching aroud 0.3 (ATLAS-CONF-2015-037)
      matchJets(genjets, newjets, 0.3);
      addHitsToJet(newjets,clean_rechits, 0., DeltaR);
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
      matchJets(genjets, recojets, 0.3);
      addHitsToJet(recojets,clean_rechits, 0., DeltaR);
      sumEnergiesAroundJet(recojets,clean_rechits, 0.05);
      // fill plots
      reco_tree.fill(recojets);
      gen_plots.fill(genjets);
      reco_plots.fill(recojets);
      jet_energy_plots.fill(recojets);
      jet_angular_plots.fill(recojets);
    }
  } // end event loop
  
  //store plots in output file
  TFile outfile(argv[2],"RECREATE");
  
  // store plots in output root tree
  reco_tree.write();
  gen_plots.write();
  reco_plots.write();
  jet_energy_plots.write();
  jet_angular_plots.write();

  outfile.Close();

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
           RecHit *r = rechits.at(i);
           if (fabs(r->eta()) < eta2 && fabs(r->eta()) > eta1 )
               layer_energy_vector[r->layer()][j].push_back(r->energy());
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
           RecHit *r = rechits.at(i);
           if (fabs(r->eta()) < eta2 && fabs(r->eta()) > eta1 )
              r->setPuOffset(medians[r->layer()][j]);
       }
   }
}


//------------------------------------------------------------------------------------------------------
template <class Collection>
void produceJets(Collection& constituents, JetCollection& jets, const float& r, const selection cuts, const bool doPuSubtraction = false, const bool doSubstructure = false) {
      
   // first convert constituents into fastjet pseudo-jets
   vector <PseudoJet> input_particles;
   for (unsigned i = 0; i < constituents.size(); i++) {
      //Constituent pj = constituents.at(i);
      input_particles.push_back( PseudoJet( constituents.at(i)->px(), constituents.at(i)->py(), constituents.at(i)->pz(), constituents.at(i)->energy()));
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
       convertJets(clust_seq, akjets, r, jets, doSubstructure);

       // soft killer PU (could be used later)

       /*double grid_size = 0.4;
       SoftKiller soft_killer(-3.0, 3.0, grid_size, grid_size) ;
       double pt_threshold;
       vector<PseudoJet> soft_killed_event;
       soft_killer.apply(input_particles, soft_killed_event, pt_threshold);

       ClusterSequence clust_seq_kill(soft_killed_event, jetDef);  
       vector<PseudoJet> kill_jets = clust_seq_kill.inclusive_jets();

       kill_jets = sel_jets(kill_jets);

       if(debug) cout << setprecision(4);
       if(debug) cout << "Soft Killer applied a pt threshold of " << pt_threshold << endl;

       // run things and print the result
       //----------------------------------------------------------
       if(debug) cout << "# original hard jets" << endl;
       for (unsigned int i=0; i<antikt_jets.size(); i++){
         const PseudoJet &jet = antikt_jets[i];
         if(debug) cout << "pt = " << jet.pt()
              << ", rap = " << jet.rap()
              << ", mass = " << jet.m() << endl;
       }
       if(debug) cout << endl;

       if(debug) cout << "# jets after applying the soft killer" << endl;
       for (unsigned int i=0; i<kill_jets.size(); i++){
         const PseudoJet &jet = kill_jets[i];
         if(debug) cout << "pt = " << jet.pt()
              << ", rap = " << jet.rap()
              << ", mass = " << jet.m() << endl;
       }
       if(debug) cout << endl;
       }
       */
   }
   else {
       
       ClusterSequence clust_seq(input_particles,jetDef);
       akjets  = sorted_by_pt(clust_seq.inclusive_jets());

       // apply cuts
       akjets = select_jets(akjets);

       // eventually apply substructure and store in custom dataformat
       convertJets(clust_seq, akjets, r, jets, doSubstructure);
   }
   
}



//------------------------------------------------------------------------------------------------------
template <class Sequence>
void convertJets(Sequence seq, vector<PseudoJet> pseudojets, const float& r, JetCollection& jets, const bool doSubstructure = false){

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
           NsubjettinessRatio    nSub32_beta1(3,2, OnePass_KT_Axes(), NormalizedMeasure(beta, r));

           jet.setTau1(nSub1_beta1(this_jet));
           jet.setTau2(nSub2_beta1(this_jet));
           jet.setTau3(nSub3_beta1(this_jet));
           jet.setTau21(nSub21_beta1(this_jet));
           jet.setTau32(nSub32_beta1(this_jet));

           // soft drop
           beta = 0.0;
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
   float pt;
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
void addHitsToJet(JetCollection& recojets, RecHitCollection& recHits, float minR, float dr){
  
  float maxR = minR + dr;

  for (unsigned i = 0; i < recojets.size() ; i++) {
    Jet *rj = recojets.at(i);
    
    std::vector<RecHit*> hits;
    // sums up all reHits aroud DeltaR                          
    for (unsigned j = 0; j < recHits.size() ; j++) {
      RecHit *hit = recHits.at(j);
      float dr_gh = hit->p4().DeltaR(rj->p4());
      if( dr_gh > minR && dr_gh < maxR ){
	// Add hit to new jet
	hits.push_back(hit);
      }
    }
    rj->setHits(hits);
  }
}
//------------------------------------------------------------------------------------------------------
 
void sumEnergiesAroundJet(JetCollection& recojets, RecHitCollection& recHits, float alpha){

  // Five rings around jet axis with max = 5*alpha
  // Alpha can effectively be maximal 0.4/5 = 0.08
  // Use 0.05 as default!
  int n = 5;

  for (unsigned i = 0; i < recojets.size() ; i++) {
    Jet *rj = recojets.at(i);

    float energies[n]={0.};
    // sums up all energies within rings
    for (unsigned j = 0; j < recHits.size() ; j++) {
      RecHit *hit = recHits.at(j);
      float dr_gh = hit->p4().DeltaR(rj->p4());
      
      for (unsigned int k=1; k<=n; k++){
	float min = alpha*(k-1)/n;
	float max = alpha*k/n;

	if( dr_gh > min && dr_gh < max ){
	  energies[k-1] += hit->energy();
	}
      }
    }

//    std::cout << "Test ring energy n=1: " << energies[0] << endl;
//    std::cout << "Test ring energy n=2: " << energies[1] << endl;
//    std::cout << "Test ring energy n=3: " << energies[2] << endl;
//    std::cout << "Test ring energy n=4: " << energies[3] << endl;
//    std::cout << "Test ring energy n=5: " << energies[4] << endl;

    rj->setEf1(energies[0]);
    rj->setEf2(energies[1]);
    rj->setEf3(energies[2]);
    rj->setEf4(energies[3]);
    rj->setEf5(energies[4]);
  
  }
}
