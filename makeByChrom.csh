#!/bin/csh -efv

set faFilter = $GOBIN/faFilter
set faBin = $GOBIN/faBin
foreach line (`cat /work/cf189/runPairwiseAlignments/ref.list`)
set genome=/work/cf189/runPairwiseAlignments/genomes/$line/$line.fa

set dir = pairwise/ 
#output directory
mkdir -p $dir/$line.byChrom

zcat $genome.gz | $faFilter -minSize=20000 stdin stdout | $faBin -minSize=100000000 -assembly=$line stdin $dir/$line.byChrom/
echo $genome
end
echo DONE
