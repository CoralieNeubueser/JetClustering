#ifndef JETENERGYRESOLUTIONPLOTS_H
#define JETENERGYRESOLUTIONPLOTS_H

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

class JetEnergyResolutionPlots{

    private:

        TString _name;

  vector<pair<pair<float,float>, TH1F>> _energybins;


    public:
        
        // constructors
        JetEnergyResolutionPlots(const TString, const vector<float> ptvals);
        void fill(JetCollection&);
        void write();

};

#endif 
