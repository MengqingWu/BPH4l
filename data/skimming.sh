#!/bin/sh

# selection (FIXME: if you want to use)
zzPair="(nllnunu>0)"
zWindow="(abs(llnunu_l1_mass-91.1876)<20.0)"
zptCut="(llnunu_l1_pt>100.0)"
metCut="(llnunu_l2_pt>0.0)"

selection="(1)"
#selection=$subLep"&&"$signalLep"&&"$zzPair"&&"$zWindow
#selection=$subLep"&&"$zzPair"&&"$zWindow

# compile
g++ skimming.cc -o skimming.exe `root-config --cflags` `root-config --libs`

#inputs
inputdir=./bph4l_8025_2017Feb15
#outputdir=AnalysisRegion
outputdir=$inputdir"/treeOut"
mkdir -p ${outputdir}

n=0

for infile in $inputdir/*/twoLeptonTreeProducerOnia/tree.root ;
do
  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"
  outfile="${outfile/\/twoLeptonTreeProducerOnia\/tree/}"

  inSkimFile=${infile/twoLeptonTreeProducerOnia\/tree.root/skimAnalyzerCount\/SkimReport.txt}

  #echo $inSkimFile
  AllEvents=`grep "All Events" ${inSkimFile} | awk {'print $3'}`
  SumWeights=`grep "Sum Weights" ${inSkimFile} | awk {'print $3'}`

  if [ -z $AllEvents ]; then
    echo Fail to get All Events from file ${inSkimFile}
    continue
  fi
  if [ -z $SumWeights ]; then
    SumWeights=$AllEvents
  fi

  echo -- Input file: $infile
  echo -- Output file: $outfile
  echo -- AllEvents: $AllEvents , SumWeights: $SumWeights
  echo -- Selection: $selection
  
  ./skimming.exe $infile $outfile $AllEvents $SumWeights $selection &
  
  n=$(( n + 1 ))

  # if [ "$n" -gt "0" ]; then
  #     break
  # fi # to test with only one sample run

  # if [ "$n" -gt "4" ]; then
  #     wait
  #     n=0
  # fi # to have 4 processes run in parallel

done

