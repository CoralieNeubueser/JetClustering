#ifndef JET_H
#define JET_H

#include "TLorentzVector.h"
#include "RecHit.hh"

class Jet{

    private:
        TLorentzVector _mom;
        float _tau1;
        float _tau2;
        float _tau3;
        float _tau21;
        float _tau32;
        float _sdmass;
        float _ef1;
        float _ef2;
        float _ef3;
        float _ef4;
        float _ef5;

        Jet* _ref = NULL;
        std::vector<RecHit*> _hits;

    public:
        
        // constructors
        Jet();
        Jet(TLorentzVector p4);
        Jet(Jet& j);

        void setTau1(float tau1);
        void setTau2(float tau2);
        void setTau3(float tau3);
        void setTau21(float tau21);
        void setTau32(float tau32);
        void setSDmass(float sdmass);
        void setEf1(float ef1);
        void setEf2(float ef2);
        void setEf3(float ef3);
        void setEf4(float ef4);
        void setEf5(float ef5);

        void setRef(Jet*);
        void setHits(std::vector<RecHit*>);

        float eta();
        float phi();
        float pt();
        float energy();
        float px();
        float py();
        float pz();
        float mass();

        float tau1();
        float tau2();
        float tau3();
        float tau21();
        float tau32();

        float ef1();
        float ef2();
        float ef3();
        float ef4();
        float ef5();

        float massSD();

        TLorentzVector p4();

        Jet* ref();
        std::vector<RecHit*> hits();
      
        void print();

};

#endif 
