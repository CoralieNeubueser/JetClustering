#include "JetTree.hh"
#include<iostream>
#include<string>

JetTree::JetTree(){

   events->Branch("pt",     &_pt     ); 
   events->Branch("pt_005",     &_pt_005     );
   events->Branch("eta",    &_eta    );
   events->Branch("phi",    &_phi    );
   events->Branch("mass",   &_mass   ); 
   events->Branch("energy", &_energy ); 
   events->Branch("tau1",   &_tau1   );  
   events->Branch("tau2",   &_tau2   );
   events->Branch("tau3",   &_tau3   );
   events->Branch("tau21",  &_tau21  ); 
   events->Branch("tau31",  &_tau31  ); 
   events->Branch("tau32",  &_tau32  );
   events->Branch("massSD", &_massSD ); 
   events->Branch("massSD_02", &_massSD_02 );
   events->Branch("massSD_04", &_massSD_04 );
   events->Branch("massSD_08", &_massSD_08 );
   events->Branch("ef1", &_ef1 );
   events->Branch("ef2", &_ef2 );
   events->Branch("ef3", &_ef3 );
   events->Branch("ef4", &_ef4 );
   events->Branch("ef5", &_ef5 );
   events->Branch("energy_clusters",    &_energy_clusters );
   events->Branch("energy_gens",        &_energy_gens );
   events->Branch("energy_hadrons",     &_energy_hadrons );
   events->Branch("num_clusters",    &_num_clusters );
   events->Branch("num_tracks",    &_num_tracks );
   events->Branch("num_pfas",    &_num_pfas );
   events->Branch("num_neutrals",    &_num_neutrals );
   events->Branch("num_gens",        &_num_gens );
   events->Branch("num_hadrons",     &_num_hadrons );

   events->Branch("gen_pt",     &_gen_pt     );
   events->Branch("gen_eta",    &_gen_eta    );
   events->Branch("gen_phi",    &_gen_phi    );
   events->Branch("gen_mass",   &_gen_mass   );
   events->Branch("gen_energy", &_gen_energy );   
}


JetTree::JetTree(double rad){
  std::string r = "0"+std::to_string(int(rad*10));
  std::cout << "name: pt_" << r << std::endl;
  events->Branch(("pt_"+r).c_str(),     &_pt     );
  events->Branch(("eta_"+r).c_str(),    &_eta    );
  events->Branch(("phi_"+r).c_str(),    &_phi    );
  events->Branch(("mass_"+r).c_str(),   &_mass   );
  events->Branch(("energy_"+r).c_str(), &_energy );
  events->Branch(("tau1_"+r).c_str(),   &_tau1   );
  events->Branch(("tau2_"+r).c_str(),   &_tau2   );
  events->Branch(("tau3_"+r).c_str(),   &_tau3   );
  events->Branch(("tau21_"+r).c_str(),  &_tau21  );
  events->Branch(("tau31_"+r).c_str(),  &_tau31  );
  events->Branch(("tau32_"+r).c_str(),  &_tau32  );
  events->Branch(("massSD_"+r).c_str(), &_massSD );
  events->Branch("ef1", &_ef1 );
  events->Branch("ef2", &_ef2 );
  events->Branch("ef3", &_ef3 );
  events->Branch("ef4", &_ef4 );
  events->Branch("ef5", &_ef5 );
  events->Branch(("energy_clusters_"+r).c_str(),    &_energy_clusters );
  events->Branch(("energy_gens_"+r).c_str(),        &_energy_gens );
  events->Branch(("energy_hadrons_"+r).c_str(),     &_energy_hadrons );
  events->Branch(("num_clusters_"+r).c_str(),    &_num_clusters );
  events->Branch(("num_tracks_"+r).c_str(),    &_num_tracks );
  events->Branch(("num_pfas_"+r).c_str(),    &_num_pfas );
  events->Branch(("num_neutrals_"+r).c_str(),    &_num_neutrals );
  events->Branch(("num_gens_"+r).c_str(),        &_num_gens );
  events->Branch(("num_hadrons_"+r).c_str(),     &_num_hadrons );
  events->Branch(("gen_pt_"+r).c_str(),     &_gen_pt     );
  events->Branch(("gen_eta_"+r).c_str(),    &_gen_eta    );
  events->Branch(("gen_phi_"+r).c_str(),    &_gen_phi    );
  events->Branch(("gen_mass_"+r).c_str(),   &_gen_mass   );
  events->Branch(("gen_energy_"+r).c_str(), &_gen_energy );
}


void JetTree::fill(JetCollection& coll){
  
  for (unsigned i = 0; i < coll.size(); i++) {

    _pt    =0.;
    _eta   =0.;
    _phi   =0.;
    _mass  =0.;
    _energy=0.;
    _tau1  =0.;
    _tau2  =0.;
    _tau3  =0.;
    _tau21 =0.;
    _tau31 =0.;
    _tau32 =0.;
    _massSD=0.;
    _massSD_02=0.;
    _massSD_04=0.;
    _massSD_08=0.;
    _pt_005=0.;
    _ef1=0.;
    _ef2=0.;
    _ef3=0.;
    _ef4=0.;
    _ef5=0.;
    _energy_clusters=0.;
    _energy_gens=0.;
    _energy_hadrons=0.;
    _num_tracks=0;
    _num_pfas=0;
    _num_neutrals=0;
    _num_clusters=0;
    _num_gens=0;
    _num_hadrons=0;
 
    _gen_pt    =0.;
    _gen_eta   =0.;
    _gen_phi   =0.;
    _gen_mass  =0.;
    _gen_energy=0.;

    // fill first two only for now
    if(i > 1) break;

    Jet jet = *(coll.at(i));
    _pt      = jet.pt();
    _pt_005      = jet.pt_005();
    _eta     = jet.eta();
    _phi     = jet.phi();
    _mass    = jet.mass();
    _energy  = jet.energy();
    _tau1    = jet.tau1();
    _tau2    = jet.tau2();
    _tau3    = jet.tau3();
    _tau21   = jet.tau21();
    _tau31   = jet.tau31();
    _tau32   = jet.tau32();
    _massSD  = jet.massSD();
    _massSD_02  = jet.massSD_02();
    _massSD_04  = jet.massSD_04();
    _massSD_08  = jet.massSD_08();
    _ef1     = jet.ef1();
    _ef2     = jet.ef2();
    _ef3     = jet.ef3();
    _ef4     = jet.ef4();
    _ef5     = jet.ef5();
    _energy_clusters = jet.energy_clusters();
    _energy_gens = jet.energy_gens();
    _energy_hadrons = jet.energy_hadrons();
    _num_clusters         = jet.num_clusters();
    _num_tracks         = jet.num_tracks();
    _num_pfas         = jet.num_pfas();
    _num_neutrals         = jet.num_neutrals();
    _num_gens         = jet.num_gens();
    _num_hadrons      = jet.num_hadrons();

    Jet *genjet = jet.ref();
    if(genjet){
      _gen_pt      = genjet->pt();
      _gen_eta     = genjet->eta();
      _gen_phi     = genjet->phi();
      _gen_mass    = genjet->mass();
      _gen_energy  = genjet->energy();
    }

    events->Fill();
  }
}

void JetTree::write(){
   
  events->Write();

}
