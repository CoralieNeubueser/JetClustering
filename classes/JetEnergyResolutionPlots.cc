#include "JetEnergyResolutionPlots.hh"
#include "fastjet/PseudoJet.hh"

JetEnergyResolutionPlots::JetEnergyResolutionPlots(const TString name, const vector<float> ptvals){
   
   _name    = name;   

   //create dictionary for resolution plots                                                                                                                                                                                
   for(vector<float>::const_iterator it = ptvals.begin(); it != ptvals.end()-1; it++) {
     float ptmin = *it;
     float ptmax = *(it+1);
     TString energy_str;
     energy_str.Form("_energy_%.0f_%.0f",ptmin,ptmax);
     energy_str = name + energy_str;
     _energybins.push_back(make_pair(make_pair(ptmin, ptmax), TH1F(energy_str,energy_str,2000,0.,4.0)));  
   }
}

void JetEnergyResolutionPlots::fill(JetCollection& col){
   for (unsigned i = 0; i < col.size(); i++) {
     Jet jet = *(col.at(i));
     Jet *genjet = jet.ref();
     
     // fill resolution plots                                                  
     for(vector<pair<pair<float,float>, TH1F>>::iterator it = _energybins.begin(); it != _energybins.end(); it++) {
       float ptmin=(it->first).first;
       float ptmax=(it->first).second;
       TH1F histo = it->second;
       if(genjet){
	 if(genjet->pt() > ptmin && genjet->pt() < ptmax){
	   (it->second).Fill(jet.energy() / genjet->energy());
	 }
       }
     }    
   }
}

void JetEnergyResolutionPlots::write(){
   
   gDirectory->mkdir(_name);
   gDirectory->cd(_name);
   
   // now fill resolution plots                                                                                                                                                                                                             
   gDirectory->mkdir("reso");
   gDirectory->cd("reso");

   for(vector<pair<pair<float,float>, TH1F>>::const_iterator it = _energybins.begin(); it != _energybins.end(); it++) {
     (it->second).Write();
   }
   gDirectory->cd("../");
}
