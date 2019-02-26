#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h>
#include "Jet.hh"
#include "RecHit.hh"

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
    _tau32 = j._tau32;
    _sdmass = j._sdmass;
    _ref = j._ref;
    _hits = j._hits;
    _ef1 = j._ef1;
    _ef2 = j._ef2;
    _ef3 = j._ef3;
    _ef4 = j._ef4;
    _ef5 = j._ef5;

}


void Jet::setTau1(float tau1){_tau1 = tau1;}
void Jet::setTau2(float tau2){_tau2 = tau2;}
void Jet::setTau3(float tau3){_tau3 = tau3;}
void Jet::setTau21(float tau21){_tau21 = tau21;}
void Jet::setTau32(float tau32){_tau32 = tau32;}

void Jet::setRef(Jet* jet){_ref = jet;}
void Jet::setHits(std::vector<RecHit*> hits){_hits = hits;}

void Jet::setEf1(float ef1){_ef1 = ef1;}
void Jet::setEf2(float ef2){_ef2 = ef2;}
void Jet::setEf3(float ef3){_ef3 = ef3;}
void Jet::setEf4(float ef4){_ef4 = ef4;}
void Jet::setEf5(float ef5){_ef5 = ef5;}

void Jet::setSDmass(float sdmass){_sdmass = sdmass;}

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
float Jet::tau32() {return _tau32;}
float Jet::massSD() {return _sdmass;}
float Jet::ef1() {return _ef1;}
float Jet::ef2() {return _ef2;}
float Jet::ef3() {return _ef3;}
float Jet::ef4() {return _ef4;}
float Jet::ef5() {return _ef5;}


Jet* Jet::ref() {return _ref;}
std::vector<RecHit*> Jet::hits() {return _hits;}

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
