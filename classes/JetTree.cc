#include "JetTree.hh"

JetTree::JetTree(){

   events->Branch("pt",     &_pt     ); 
   events->Branch("eta",    &_eta    );
   events->Branch("phi",    &_phi    );
   events->Branch("mass",   &_mass   ); 
   events->Branch("energy", &_energy ); 
   events->Branch("tau1",   &_tau1   );  
   events->Branch("tau2",   &_tau2   );
   events->Branch("tau3",   &_tau3   );
   events->Branch("tau21",  &_tau21  ); 
   events->Branch("tau32",  &_tau32  ); 
   events->Branch("massSD", &_massSD ); 

}

void JetTree::fill(JetCollection& coll){
   
  _pt    .clear();
  _eta   .clear();
  _phi   .clear();
  _mass  .clear();
  _energy.clear();
  _tau1  .clear();
  _tau2  .clear();
  _tau3  .clear();
  _tau21 .clear();
  _tau32 .clear();
  _massSD.clear();
  
  for (unsigned i = 0; i < coll.size(); i++) {

      // fill first two only for now
       if(i > 1) break;

       Jet jet = *(coll.at(i));
       _pt      .push_back(jet.pt());
       _eta     .push_back(jet.eta());
       _phi     .push_back(jet.phi());
       _mass    .push_back(jet.mass());
       _energy  .push_back(jet.energy());
       _tau1    .push_back(jet.tau1());
       _tau2    .push_back(jet.tau2());
       _tau3    .push_back(jet.tau3());
       _tau21   .push_back(jet.tau21());
       _tau32   .push_back(jet.tau32());
       _massSD  .push_back(jet.massSD());
  }
  
  events->Fill();
}

void JetTree::write(){
   
  events->Write();

}
