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
   events->Branch("ef1", &_ef1 );
   events->Branch("ef2", &_ef2 );
   events->Branch("ef3", &_ef3 );
   events->Branch("ef4", &_ef4 );
   events->Branch("ef5", &_ef5 );

   //events->Branch("jets_massSD", &_jets_massSD ); 
   //events->Branch("jets_mass",   &_jets_mass );   
   //events->Branch("jets_pt",     &_jets_pt     ); 
   
}

void JetTree::fill(JetCollection& coll){

  //_jets_pt    .clear();   
  //_jets_mass  .clear();   
  //_jets_massSD.clear();   
  //
  //for (unsigned i = 0; i < coll.size(); i++) {
  //  // fill first two only for now 
  //  if(i > 1) break;
  //  Jet jet = *(coll.at(i));
  //  _jets_pt      .push_back(jet.pt());
  //  _jets_mass    .push_back(jet.mass());
  //  _jets_massSD  .push_back(jet.massSD());
  //}
  
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
    _tau32 =0.;
    _massSD=0.;
    _ef1=0.;
    _ef2=0.;
    _ef3=0.;
    _ef4=0.;
    _ef5=0.;

    // fill first two only for now
    if(i > 1) break;

    Jet jet = *(coll.at(i));
    _pt      = jet.pt();
    _eta     = jet.eta();
    _phi     = jet.phi();
    _mass    = jet.mass();
    _energy  = jet.energy();
    _tau1    = jet.tau1();
    _tau2    = jet.tau2();
    _tau3    = jet.tau3();
    _tau21   = jet.tau21();
    _tau32   = jet.tau32();
    _massSD  = jet.massSD();
    _ef1     = jet.ef1();
    _ef2     = jet.ef2();
    _ef3     = jet.ef3();
    _ef4     = jet.ef4();
    _ef5     = jet.ef5();

    events->Fill();
  }
}

void JetTree::write(){
   
  events->Write();

}
