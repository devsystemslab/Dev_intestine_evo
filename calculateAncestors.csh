#!/bin/csh -exf
set db = "hg38"
#db is our reference
set threshold = 0.33
mkdir -p coverageCalls
rm -rf coverageCalls
mkdir coverageCalls

#you will need an installation of kentutils for this version, and instead of this first part where we give our absolute path to overlapSelect you will need to use yours, same with the rest of the paths to files. Keep the coverageCalls directory creation above and its usage throughout. 
 /hpc/group/vertgenlab/softwareShared/kent/kent.Jul.30.2021/overlapSelect -mergeOutput /work/cf189/birthData/allSpecies.bed hg38OrganoidAtac.campLab.4col.bed stdout \
| /work/cf189/birthData/punchHolesInBed.pl /dev/stdin coverageCalls/tmp.a.bed

 foreach x (`cat ../speciesNode.list`)

	set species = $x:r
        set node = $x:e
        echo $species $node

        /hpc/group/vertgenlab/softwareShared/kent/kent.Jul.30.2021/overlapSelect -aggregate -overlapThreshold=$threshold /work/cf189/birthData/alignmentData/$db.wholeGenomes/$db.$species.bed coverageCalls/tmp.a.bed coverageCalls/tmp.out.bed

        cat coverageCalls/tmp.out.bed >> coverageCalls/node"$node".bed
        cat coverageCalls/tmp.out.bed | cut -f 4 | awk -v "NODE=$node" '{ print $1"\t"NODE }' >> coverageCalls/cneNode.tsv 

        /hpc/group/vertgenlab/softwareShared/kent/kent.Jul.30.2021/overlapSelect -nonOverlapping coverageCalls/node"$node".bed coverageCalls/tmp.a.bed coverageCalls/tmp.b.bed

        mv coverageCalls/tmp.b.bed coverageCalls/tmp.a.bed

        rm coverageCalls/tmp.out.bed
end
  
mv coverageCalls/tmp.a.bed coverageCalls/refSpecific.bed
wc -l coverageCalls/refSpecific.bed

echo DONE
