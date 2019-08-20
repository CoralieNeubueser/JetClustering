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
  float _pt_005;
        float _eta;
        float _phi;
        float _mass;
        float _energy;
        float _tau1;
        float _tau2;
        float _tau3;
        float _tau21;
        float _tau31;
        float _tau32;
        float _massSD;
  float _massSD_02;
  float _massSD_04;
  float _massSD_08;
        float _ef1;
        float _ef2;
        float _ef3;
        float _ef4;
        float _ef5;
        int _num_clusters=0;
        int _num_tracks=0;
        int _num_pfas=0;
        int _num_neutrals=0;
        int _num_gens=0;
  //    int _num_gens_1GeV=0;
        int _num_hadrons=0;
  //     int _num_hadrons_1GeV=0;
        float _energy_clusters=0.;
        float _energy_gens=0.;
  //      float _energy_gens_1GeV=0.;
        float _energy_hadrons=0.;
  //        float _energy_hadrons_1GeV=0.;

        float _gen_pt;
        float _gen_eta;
        float _gen_phi;
        float _gen_mass;
        float _gen_energy;

    public:
        
        // constructors
        JetTree();
        JetTree(double rad);
        void fill(JetCollection&);
        void write();

};

#endif 
