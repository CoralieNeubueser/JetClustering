#ifndef JET_H
#define JET_H

#include "TLorentzVector.h"
#include "RecHit.hh"
#include "GenParticle.hh"
#include "GenParticleCollection.hh"

class Jet{

    private:
        TLorentzVector _mom;
        float _tau1;
        float _tau2;
        float _tau3;
        float _tau21;
        float _tau31;
        float _tau32;
  float _sdmass, _sdmass_02, _sdmass_04, _sdmass_08;
  float _pt_005;
        float _ef1;
        float _ef2;
        float _ef3;
        float _ef4;
        float _ef5;
        int _num_tracks=0;
        int _num_pfas=0;
        int _num_neutrals=0;
        int _num_clusters=0;
        int _num_gens=0;
  //        int _num_gens_leq_1GeV=0;
        int _num_hadrons=0;
  //      int _num_hadrons_leq_1GeV=0;
        float _energy_clusters=0.;
        float _energy_gens=0.;
  //     float _energy_gens_leq_1GeV=0.;
        float _energy_hadrons=0.;
  //    float _energy_hadrons_leq_1GeV=0.;

        Jet* _ref = NULL;
  //        std::vector<RecHit*> _hits;
        GenParticleCollection* _gens = NULL;

    public:
        
        // constructors
        Jet();
        Jet(TLorentzVector p4);
        Jet(Jet& j);

        void setTau1(float tau1);
        void setTau2(float tau2);
        void setTau3(float tau3);
        void setTau21(float tau21);
        void setTau31(float tau31);
        void setTau32(float tau32);
        void setSDmass(float sdmass);
  void setSDmass_02(float sdmass02);
  void setSDmass_04(float sdmass04);
  void setSDmass_08(float sdmass08);
  void setPt_005(float pt_005);      
  void setEf1(float ef1);
        void setEf2(float ef2);
        void setEf3(float ef3);
        void setEf4(float ef4);
        void setEf5(float ef5);
        void setNum_tracks(int num_tracks);
        void setNum_pfas(int num_pfas);
        void setNum_neutrals(int num_neutrals);
        void setNum_clusters(int num_clusters);
        void setNum_gens(int num_gens);
  //    void setNum_gens_leq_1GeV(int num_gens_leq_1GeV);
        void setNum_hadrons(int num_hadrons);
  //    void setNum_hadrons_leq_1GeV(int num_hadrons_leq_1GeV);      
        void setEnergy_clusters(float energy_clusters);
        void setEnergy_gens(float energy_gens);
  //    void setEnergy_gens_leq_1GeV(float energy_gens_leq_1GeV);
        void setEnergy_hadrons(float energy_hadrons);
  //    void setEnergy_hadrons_leq_1GeV(float energy_hadrons_leq_1GeV);

        void setRef(Jet*);
  //    void setHits(std::vector<RecHit*>);
        void setParticles(GenParticleCollection*);
  
        float eta();
        float phi();
        float pt();
  float pt_005();
        float energy();
        float px();
        float py();
        float pz();
        float mass();

        float tau1();
        float tau2();
        float tau3();
        float tau21();
        float tau31();
        float tau32();

        float ef1();
        float ef2();
        float ef3();
        float ef4();
        float ef5();

        int num_tracks();
        int num_pfas();
        int num_neutrals();
        int num_clusters();
        int num_gens();
  //    int num_gens_leq_1GeV();
        int num_hadrons();
  //    int num_hadrons_leq_1GeV();
      
        float energy_clusters();
        float energy_gens();
  //    float energy_gens_leq_1GeV();
        float energy_hadrons();
  //    float energy_hadrons_leq_1GeV();

        float massSD();
  float massSD_02();
  float massSD_04();
  float massSD_08();

        TLorentzVector p4();

        Jet* ref();
  //    std::vector<RecHit*> hits();
        GenParticleCollection* gens();

        void print();

};

#endif 
