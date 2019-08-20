#ifndef JETANGULARRESOLUTIONPLOTS_H
#define JETANGULARRESOLUTIONPLOTS_H

#include "Jet.hh"
#include "JetCollection.hh"
#include <vector>
#include <TH1F.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TDirectory.h>
#include <iostream>

using namespace std;

class JetAngularResolutionPlots{

    private:

        TString _name;

  vector<pair<pair<float,float>, TH1F>> _etabins;
  vector<pair<pair<float,float>, TH1F>> _thetabins;
  vector<pair<pair<float,float>, TH1F>> _phibins;

    public:
        
        // constructors
        JetAngularResolutionPlots(const TString, const vector<float> ptvals);
        void fill(JetCollection&);
        void write();

};

#endif 
