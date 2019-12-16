#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/20GeV/reco/positions/ -n 50 -o 20GeVljets --njobs 171 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/50GeV/reco/positions/ -n 50 -o 50GeVljets --njobs 553 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/100GeV/reco/positions/ -n 50 -o 100GeVljets --njobs 794 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/200GeV/reco/positions/ -n 50 -o 200GeVljets --njobs 203 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/500GeV/reco/positions/ -n 50 -o 500GeVljets --njobs 670 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/1000GeV/reco/positions/ -n 50 -o 1000GeVljets --njobs 711 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/2000GeV/reco/positions/ -n 50 -o 2000GeVljets --njobs 353 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/5000GeV/reco/positions/ -n 50 -o 5000GeVljets --njobs 665 -q 1nh
#python submitJetClustering.py -i /eos/experiment/fcc/hh/simulation/samples/v02_pre/physics/ljets/bFieldOn/etaTo1.5/10000GeV/reco/positions/ -n 50 -o 10000GeVljets --njobs 381 -q 1nh

find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/20GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/50GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/100GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/200GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/500GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/1000GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/2000GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/5000GeVljets/out/* -size 0 -delete
find /eos/experiment/fcc/users/c/cneubuse/JetClustering/antiKt/10000GeVljets/out/* -size 0 -delete

python submitJetClustering.py --collect antiKt --collectPt 20
python submitJetClustering.py --collect antiKt --collectPt 50
python submitJetClustering.py --collect antiKt --collectPt 100
python submitJetClustering.py --collect antiKt --collectPt 200
python submitJetClustering.py --collect antiKt --collectPt 500
python submitJetClustering.py --collect antiKt --collectPt 1000
python submitJetClustering.py --collect antiKt --collectPt 2000
python submitJetClustering.py --collect antiKt --collectPt 5000
python submitJetClustering.py --collect antiKt --collectPt 10000

python submitJetClustering.py --collect antiKt --allPts
