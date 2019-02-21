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

        std::vector<double> _pt;
        std::vector<double> _eta;
        std::vector<double> _phi;
        std::vector<double> _mass;
        std::vector<double> _energy;
        std::vector<double> _tau1;
        std::vector<double> _tau2;
        std::vector<double> _tau3;
        std::vector<double> _tau21;
        std::vector<double> _tau32;
        std::vector<double> _massSD;

    public:
        
        // constructors
        JetTree();
        void fill(JetCollection&);
        void write();

};

#endif 
