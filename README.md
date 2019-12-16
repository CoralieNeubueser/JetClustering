[]() Package installation
--------------------------
First load the required environment on lxplus:
```
source init.sh
```
Install fastjet (to be done only once):
```
./install-fastjet.sh
```

[]() Run instructions
----------------------


First load the required environment on lxplus:
```
source init.sh
```
Compile analysis code (can be found and modified in ```src/analyze.cc```):

```
make -j 4
```
Run analysis:
```
./analyze [input_file] [output_file] [number of events] [jet cone size DeltaR] [maximum eta range] [run a check of jet pt versus sum over collected objects] [beta (sub-jettiness)] [alpha (energy flow cone size)]
```
e.g:
```
./analyze /eos/experiment/fcc/hh/simulation/samples/tracker_calobarrel/ditop/500GeV/NTUP/output_helsens_20171011151211690.root test.root 10 0.4 0.5 0 0 0.05
```


[]() Condor submission
--------------------

The script ```batch/submitJetClusteringEos.py``` allows to run this script on condor queues.
A new output directory is automatically build from the input string, following the logic of FCCSimJobs version v04 by default.
Specify algorithm: "anti-kt" (if calo cells have been used in full granularity), "anti-kt_resegmentedHCal" (if input default calo cell granularity) or "anti-kt_split" (if cluster were split) and if input are cluster add: -t
``
example:
``
```
python submitJetClusteringEos.py --version v04 --etaMax 0.5 --etaCut 0.5  --algorithm anti-kt_resegmentedHCal -i /eos/experiment/fcc/hh/simulation/samples/v04/physics/ljets/bFieldOn/etaTo0.5/10000GeV/ntup/topoClusters/electronicsNoise/resegmentedHCal/calibrated/benchmark/420/ -p 10000 -n 50 --njobs 4000 --deltaR 0.4 --alpha 0.05 --topoCluster --process ljets --beta 0.
```
Jobs will be collected for certain jet alorithm and pt, in directory /eos/experiment/fcc/hh/simulation/samples/v04/physics/ljets/bFieldOn/etaTo0.5/10000GeV/ana/JetClustering/anti-kt_resegmentedHCal/
An example to send all possible reconstructions of a certain process and pt can be specified in batch/send.sh
To merge all files of a certain jet reconstruction, use batch/merge.sh script by specifying the reco that you want to collect in one file e.g. to be used as input for the TMVA analysis..

