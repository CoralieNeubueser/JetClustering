delta=0.2
algo="anti-kt_track" #resegmentedHCal_split_cluster" #resegmentedHCal_split_cluster" #resegmentedHCal_split_cluster"
bfield="bFieldOn"
version="v04"
process="Wqq" #"Wqq" #"ljets" #top"
alpha=0.25

echo ${delta}
echo ${algo}

for en in 10000;
do
    echo ${en}
    for i in {1..9}; 
    do 
	echo $i;
	if [ $alpha == 0.05 ]
	then
	    hadd -f -k /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/all_${i}.root /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/${en}GeV_*_${i}*.root;
	else
	    hadd -f -k /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/alpha_${alpha}/all_${i}.root /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/alpha_${alpha}/${en}GeV_*_${i}*.root;
	fi
    done
    if [ $alpha == 0.05 ]
    then
	hadd -f -k /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/all.root /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/all_*.root;
    else
	hadd -f -k /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/alpha_${alpha}/all.root /eos/experiment/fcc/hh/simulation/samples/${version}/physics/${process}/${bfield}/etaTo0.5/${en}GeV/ana/JetClustering/${algo}/DeltaR_${delta}/maxEta0.5/alpha_${alpha}/all_*.root;
    fi
done
#rm /eos/experiment/fcc/hh/simulation/samples/v03/physics/Wqq/bFieldOff/etaTo0.5/500GeV/ana/JetClustering/anti-kt_resegmentedHCal/DeltaR_0.4/maxEta0.5/500GeV_Wqq*
#rm /eos/experiment/fcc/hh/simulation/samples/v03/physics/Wqq/bFieldOff/etaTo0.5/1000GeV/ana/JetClustering/anti-kt_resegmentedHCal/DeltaR_0.4/maxEta0.5/1000GeV_Wqq*
#rm /eos/experiment/fcc/hh/simulation/samples/v03/physics/Wqq/bFieldOff/etaTo0.5/2000GeV/ana/JetClustering/anti-kt_resegmentedHCal/DeltaR_0.4/maxEta0.5/2000GeV_Wqq*
#rm /eos/experiment/fcc/hh/simulation/samples/v03/physics/Wqq/bFieldOff/etaTo0.5/5000GeV/ana/JetClustering/anti-kt_resegmentedHCal/DeltaR_0.4/maxEta0.5/5000GeV_Wqq*
#rm /eos/experiment/fcc/hh/simulation/samples/v03/physics/Wqq/bFieldOff/etaTo0.5/10000GeV/ana/JetClustering/anti-kt_resegmentedHCal/DeltaR_0.4/maxEta0.5/10000GeV_Wqq*
