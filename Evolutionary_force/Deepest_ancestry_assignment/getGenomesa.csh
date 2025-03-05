#!/bin/csh -evf

mkdir -p genomes

foreach name (`cat species.list`)

wget -nc -P /work/cf189/runPairwiseAlignments/genomes/$name/ https://hgdownload.soe.ucsc.edu/goldenPath/$name/bigZips/$name.fa.gz
wget -nc -P /work/cf189/runPairwiseAlignments/genomes/$name/ https://hgdownload.soe.ucsc.edu/goldenPath/$name/bigZips/$name.chrom.sizes
wget -nc -P /work/cf189/runPairwiseAlignments/genomes/$name/ https://hgdownload.soe.ucsc.edu/goldenPath/$name/bigZips/$name.2bit

echo $name

end

echo DONE
