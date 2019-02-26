#ifndef JETTREE_H
#define JETTREE_H

#include "Jet.hh"
#include "JetCollection.hh"
#include <vector>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TDirectory.h>
#include <iostream>

using namespace std;

class JetTree{

    private:

        TTree* events = new TTree("events", "a Tree with data from a fake Geant3");

        float _pt;
        float _eta;
        float _phi;
        float _mass;
        float _energy;
        float _tau1;
        float _tau2;
        float _tau3;
        float _tau21;
        float _tau32;
        float _massSD;
        float _ef1;
        float _ef2;
        float _ef3;
        float _ef4;
        float _ef5;

    public:
        
        // constructors
        JetTree();
        void fill(JetCollection&);
        void write();

};

#endif 
