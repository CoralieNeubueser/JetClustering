deltaR=0.2
alpha=0.25
energy=10000
pro="ljets" #, "Wqq" or "top"
beta=0. # default

# rec: cluster
python submitJetClusteringEos.py --version v04 --etaMax 0.5 --etaCut 0.5  --algorithm anti-kt_resegmentedHCal -i /eos/experiment/fcc/hh/simulation/samples/v04/physics/${pro}/bFieldOn/etaTo0.5/${energy}GeV/ntup/topoClusters/electronicsNoise/resegmentedHCal/calibrated/benchmark/420/ -p ${energy} -n 50 --njobs 4000 --deltaR ${deltaR} --alpha ${alpha} --topoCluster --process ${pro} --beta ${beta} 
# rec: split cluster
python submitJetClusteringEos.py --version v04 --etaMax 0.5 --etaCut 0.5  --algorithm anti-kt_resegmentedHCal_split -i /eos/experiment/fcc/hh/simulation/samples/v04/physics/${pro}/bFieldOn/etaTo0.5/${energy}GeV/ntup/topoClusters/electronicsNoise/resegmentedHCal/splitted/calibrated/benchmark/420/ -p ${energy} -n 50 --njobs 8000 --deltaR ${deltaR} --alpha ${alpha} --topoCluster --process ${pro} --beta ${beta}
# rec: cells
python submitJetClusteringEos.py --version v04 --etaMax 0.5 --etaCut 0.5  --algorithm anti-kt_resegmentedHCal -i /eos/experiment/fcc/hh/simulation/samples/v04/physics/${pro}/bFieldOn/etaTo0.5/${energy}GeV/ntup/positions/resegmentedHCal/ -p ${energy} -n 50 --njobs 8000 --deltaR ${deltaR} --alpha ${alpha} --process ${pro} --beta ${beta}
# rec: full granularity cells
python submitJetClusteringEos.py --version v04 --etaMax 0.5 --etaCut 0.5  --algorithm anti-kt -i /eos/experiment/fcc/hh/simulation/samples/v04/physics/${pro}/bFieldOn/etaTo0.5/${energy}GeV/ntup/positions/ -p ${energy} -n 50 --njobs 8000 --deltaR ${deltaR} --alpha ${alpha} --process ${pro} --beta ${beta}
# rec: tracks (tracks input saved in same file as cluster)
python submitJetClusteringEos.py --version v04 --etaMax 0.5 --etaCut 0.5  --algorithm anti-kt -i /eos/experiment/fcc/hh/simulation/samples/v04/physics/${pro}/bFieldOn/etaTo0.5/${energy}GeV/ntup/topoClusters/electronicsNoise/resegmentedHCal/calibrated/benchmark/420/ -p ${energy} -n 50 --njobs 8000  --deltaR ${deltaR} --alpha ${alpha} --track --process ${pro} --beta ${beta}
