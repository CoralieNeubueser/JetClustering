#include "JetAngularResolutionPlots.hh"
#include "fastjet/PseudoJet.hh"

JetAngularResolutionPlots::JetAngularResolutionPlots(const TString name, const vector<float> ptvals){
   
   _name    = name;   

   //create dictionary for resolution plots                                                                                                                                                                                
   for(vector<float>::const_iterator it = ptvals.begin(); it != ptvals.end()-1; it++) {
     float ptmin = *it;
     float ptmax = *(it+1);
     TString eta_str;
     TString theta_str;
     TString phi_str;
     eta_str.Form("_eta_%.0f_%.0f",ptmin,ptmax);
     theta_str.Form("_theta_%.0f_%.0f",ptmin,ptmax);
     phi_str.Form("_phi_%.0f_%.0f",ptmin,ptmax);
     eta_str = name + eta_str;
     theta_str = name + theta_str;
     phi_str = name + phi_str;
     _etabins.push_back(make_pair(make_pair(ptmin, ptmax), TH1F(eta_str,eta_str,6000,-1.,1.)));
     _thetabins.push_back(make_pair(make_pair(ptmin, ptmax), TH1F(theta_str,theta_str,6000,-1.,1.)));
     _phibins.push_back(make_pair(make_pair(ptmin, ptmax), TH1F(phi_str,phi_str,6000,-1.,1.)));
   }
}

void JetAngularResolutionPlots::fill(JetCollection& col){
   for (unsigned i = 0; i < col.size(); i++) {
     Jet jet = *(col.at(i));
     Jet *genjet = jet.ref();
     
     double etaW = 0.;
     double thetaW = 0.;
     double phiW = 0.;
     double energy =0.;

     for(const auto& hit : jet.hits()){
       
       etaW += hit->eta()*hit->energy();
       thetaW += (2. * atan( exp(-hit->eta()))) *hit->energy();
       phiW += hit->phi()*hit->energy(); 
       energy += hit->energy();
     
     }
    
     // fill resolution plots                                                                                                                                                                                                              
     for(vector<pair<pair<float,float>, TH1F>>::iterator it = _etabins.begin(); it != _etabins.end(); it++) {
       float ptmin=(it->first).first;
       float ptmax=(it->first).second;
       TH1F histo = it->second;
       if(genjet){
	 if(genjet->pt() > ptmin && genjet->pt() < ptmax){
	   (it->second).Fill((etaW / energy) - genjet->eta());
	 }
       }
     }
     for(vector<pair<pair<float,float>, TH1F>>::iterator it = _thetabins.begin(); it != _thetabins.end(); it++) {
       float ptmin=(it->first).first;
       float ptmax=(it->first).second;
       TH1F histo = it->second;
       if(genjet){
         if(genjet->pt() > ptmin && genjet->pt() < ptmax){
           (it->second).Fill((thetaW / energy) - (2. * atan( exp(-genjet->eta()))));
         }
       }
     }

     for(vector<pair<pair<float,float>, TH1F>>::iterator it = _phibins.begin(); it != _phibins.end(); it++) {
       float ptmin=(it->first).first;
       float ptmax=(it->first).second;
       TH1F histo = it->second;
       if(genjet){
         if(genjet->pt() > ptmin && genjet->pt() < ptmax){
           (it->second).Fill((phiW / energy) - genjet->phi());
         }
       }
     }
     
   }
}

void JetAngularResolutionPlots::write(){
   
   gDirectory->mkdir(_name);
   gDirectory->cd(_name);
   
   // now fill resolution plots                                                                                                                                                                                                             
   gDirectory->mkdir("reso");
   gDirectory->cd("reso");

   for(vector<pair<pair<float,float>, TH1F>>::const_iterator it = _etabins.begin(); it != _etabins.end(); it++) {
     (it->second).Write();
   }
   for(vector<pair<pair<float,float>, TH1F>>::const_iterator it = _thetabins.begin(); it != _thetabins.end(); it++) {
     (it->second).Write();
   }
   for(vector<pair<pair<float,float>, TH1F>>::const_iterator it = _phibins.begin(); it != _phibins.end(); it++) {
     (it->second).Write();
   }
   gDirectory->cd("../");
}
