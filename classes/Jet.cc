#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include "Jet.hh"
#include "GenParticle.hh"
#include "GenParticleCollection.hh"

Jet::Jet(){};

Jet::Jet(TLorentzVector p4){
   _mom = p4;
}

Jet::Jet(Jet& j){
    _mom = j._mom;
    _tau1 = j._tau1;
    _tau2 = j._tau2;
    _tau3 = j._tau3;
    _tau21 = j._tau21;
    _tau31 = j._tau31;
    _tau32 = j._tau32;
    _sdmass = j._sdmass;
    _sdmass_02 = j._sdmass_02;
    _sdmass_04 = j._sdmass_04;
    _sdmass_08 = j._sdmass_08;
    _pt_005 = j._pt_005;
    _ref = j._ref;
    //    _hits = j._hits;
    _gens = j._gens;
    _ef1 = j._ef1;
    _ef2 = j._ef2;
    _ef3 = j._ef3;
    _ef4 = j._ef4;
    _ef5 = j._ef5;
    _num_tracks = j._num_tracks;
    _num_pfas = j._num_pfas;
    _num_neutrals = j._num_neutrals;
    _num_clusters = j._num_clusters;
    _num_gens = j._num_gens;
    //    _num_gens_leq_1GeV = j._num_gens_leq_1GeV;
    _num_hadrons = j._num_hadrons;
    //_num_hadrons_leq_1GeV = j._num_hadrons_leq_1GeV;
    _energy_clusters = j._energy_clusters;
    _energy_gens = j._energy_gens;
    //_energy_gens_leq_1GeV = j._energy_gens_leq_1GeV;
    _energy_hadrons = j._energy_hadrons;
    //_energy_hadrons_leq_1GeV = j._energy_hadrons_leq_1GeV;

}


void Jet::setTau1(float tau1){_tau1 = tau1;}
void Jet::setTau2(float tau2){_tau2 = tau2;}
void Jet::setTau3(float tau3){_tau3 = tau3;}
void Jet::setTau21(float tau21){_tau21 = tau21;}
void Jet::setTau31(float tau31){_tau31 = tau31;}
void Jet::setTau32(float tau32){_tau32 = tau32;}

void Jet::setRef(Jet* jet){_ref = jet;}
//void Jet::setHits(std::vector<RecHit*> hits){_hits = hits;}
void Jet::setParticles(GenParticleCollection* gens){_gens = gens;}

void Jet::setEf1(float ef1){_ef1 = ef1;}
void Jet::setEf2(float ef2){_ef2 = ef2;}
void Jet::setEf3(float ef3){_ef3 = ef3;}
void Jet::setEf4(float ef4){_ef4 = ef4;}
void Jet::setEf5(float ef5){_ef5 = ef5;}

void Jet::setNum_tracks(int num_tracks){_num_tracks = num_tracks;}
void Jet::setNum_pfas(int num_pfas){_num_pfas = num_pfas;}
void Jet::setNum_neutrals(int num_neutrals){_num_neutrals = num_neutrals;}
void Jet::setNum_clusters(int num_clusters){_num_clusters = num_clusters;}
void Jet::setNum_gens(int num_gens){_num_gens = num_gens;}
//void Jet::setNum_gens_leq_1GeV(int num_gens_leq_1GeV){_num_gens_leq_1GeV = num_gens_leq_1GeV;}
void Jet::setNum_hadrons(int num_hadrons) {_num_hadrons=num_hadrons;}
//void Jet::setNum_hadrons_leq_1GeV(int num_hadrons_leq_1GeV){_num_hadrons_leq_1GeV=num_hadrons_leq_1GeV;}

void Jet::setEnergy_clusters(float energy_clusters){_energy_clusters=energy_clusters;}
void Jet::setEnergy_gens(float energy_gens){_energy_gens=energy_gens;}
//void Jet::setEnergy_gens_leq_1GeV(float energy_gens_leq_1GeV){_energy_gens_leq_1GeV=energy_gens_leq_1GeV;}
void Jet::setEnergy_hadrons(float energy_hadrons){_energy_hadrons=energy_hadrons;}
//void Jet::setEnergy_hadrons_leq_1GeV(float energy_hadrons_leq_1GeV){_energy_hadrons_leq_1GeV=energy_hadrons_leq_1GeV;}

void Jet::setSDmass(float sdmass){_sdmass = sdmass;}
void Jet::setSDmass_02(float sdmass02){_sdmass_02 = sdmass02;}
void Jet::setSDmass_04(float sdmass04){_sdmass_04 = sdmass04;}
void Jet::setSDmass_08(float sdmass08){_sdmass_08 = sdmass08;}

void Jet::setPt_005(float pt_005){_pt_005 = pt_005;}

float Jet::eta() {return _mom.Eta();}
float Jet::phi() {return _mom.Phi();}
float Jet::pt() {return _mom.Pt();}
float Jet::energy() {return _mom.Energy();}
float Jet::px() {return _mom.Px();}
float Jet::py() {return _mom.Py();}
float Jet::pz() {return _mom.Pz();}
float Jet::mass() {return _mom.M();}
TLorentzVector Jet::p4() {return _mom;}
float Jet::tau1() {return _tau1;}
float Jet::tau2() {return _tau2;}
float Jet::tau3() {return _tau3;}
float Jet::tau21() {return _tau21;}
float Jet::tau31() {return _tau31;}
float Jet::tau32() {return _tau32;}
float Jet::massSD() {return _sdmass;}
float Jet::massSD_02() {return _sdmass_02;}
float Jet::massSD_04() {return _sdmass_04;}
float Jet::massSD_08() {return _sdmass_08;}
float Jet::pt_005() {return _pt_005;}
float Jet::ef1() {return _ef1;}
float Jet::ef2() {return _ef2;}
float Jet::ef3() {return _ef3;}
float Jet::ef4() {return _ef4;}
float Jet::ef5() {return _ef5;}
int Jet::num_tracks() {return _num_tracks;}
int Jet::num_pfas() {return _num_pfas;}
int Jet::num_neutrals() {return _num_neutrals;}
int Jet::num_clusters() {return _num_clusters;}
int Jet::num_gens() {return _num_gens;}
//int Jet::num_gens_leq_1GeV() {return _num_gens_leq_1GeV;}
int Jet::num_hadrons() {return _num_hadrons;}
//int Jet::num_hadrons_leq_1GeV() {return _num_hadrons_leq_1GeV;}
float Jet::energy_clusters() {return _energy_clusters;}
float Jet::energy_gens() {return _energy_gens;}
//float Jet::energy_gens_leq_1GeV() {return _energy_gens_leq_1GeV;}
float Jet::energy_hadrons() {return _energy_hadrons;}
//float Jet::energy_hadrons_leq_1GeV() {return _energy_hadrons_leq_1GeV;}

Jet* Jet::ref() {return _ref;}
//std::vector<RecHit*> Jet::hits() {return _hits;}
GenParticleCollection* Jet::gens() {return _gens;}

void Jet::print(){
   std::string commentStr = "";
   std::cout << commentStr << std::setprecision(4) <<  std::setw(5) << _mom.Pt()
   << std::setprecision(4) <<  std::setw(10) << _mom.Eta()
   << std::setprecision(4) <<  std::setw(11) << _mom.Phi()
   << std::setprecision(4) <<  std::setw(11) << _mom.M()
   << std::setprecision(4) <<  std::setw(11) << _mom.E();
   std::cout << std::setprecision(4) <<  std::setw(11) << _tau21;
   std::cout << std::setprecision(4) <<  std::setw(11) << _tau32;
   std::cout << std::setprecision(4) <<  std::setw(11) << _sdmass;
   std::cout << std::endl;
   
}
