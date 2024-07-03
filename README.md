# Dev_intestine_evo
Scripts and resources for the manuscript titled as Recent evolution of the developing human small intestinal epithelium

makeByChrom.csh
this will take each genome and split it into bins of similar size so computational resources can be used efficiently. This will filter out smaller contigs and then group the fastas into bins.

getGenomes.csh
this function will take the species list and grab the genomes for the analysis from UCSC as well as their chrom.sizes and 2bit files which will be necessary in later steps.

datedHumanGenome/
contains a file for each node of the tree which has human regions that are associated with that node. 

allSpecies.bed 
all of our alignable regions in human referenced alignments (this was the input file for our calculation of ancestor for human) ** this now includes zebrafish**.

Node info:
node0= human/Chimp/bonobo ancestor
node1= human/gorilla ancestor
node2= human/crab eating macaque ancestor
node3= human/marmoset ancestor
node4= human/mouse ancestor
node5= human/cow/dog ancestor
node6= placental mammal ancestor
node7= human/opossum ancestor
node8= human/platypus ancestor
node9= human/bird/alligator/turtle ancestor
node10= tetrapod ancestor (human/frog)
node11= human/ceolacanth ancestor
node12= human/fish ancestor
node13= vertebrate ancestor

alignments/
this directory contains all genome.genome alignments with the target/reference listed first. 

punchHolesInBed.pl
this is a perl script that requires the bed inputs (allSpecies.bed and whatever data its being compared to) to have 4 columns that are all unique. It will error if that isn't true (typically it'll say something like this if that's the problem: "chroms don't match:" followed by the line info that it failed on.) To make things unique we typically will just add the line number that it appears on. The resources listed here already have this formatting but for new data users will need to ensure this requirement is met. You can do that a few ways, I like to use awk:

"cat <file.bed> | cut -f1-4 | awk '{print $0"."NR}' > <file.4colUnique.bed>"

you can also use perl:

"cat <file.bed> | cut -f1-4 | perl -ne 'chomp($_);$x+=1;print("$_.$x\n");' > <file.4colUnique.bed>"

as a side note, when you copy files on the terminal they tend to translate tabs to spaces, so if that happens you can use this to correct it:

"cat <file> | tr " " "\t" > <newFile>"

campTree.png
Tree file with the nodes labelled with the corresponding number of their ancestor which is labelled as node above.

campLabTree.nh
Newick formatted file for running phast/all_dist if user wants to regenerate the phylogenetic distances. This was used as input for gonomics/drawNewickTree.go function which oututs a png file. documentation for that function can be found here: https://pkg.go.dev/github.com/vertgenlab/gonomics@v1.0.1/cmd/drawNewickTree
this tree is a adaptation of a tree avialable here: https://github.com/joelarmstrong/vertebrate-phylogeny/blob/master/vertebrates.commonNames.nh


speciesNode.list
contains information for calculating ancestors program on which species are on which nodes in relation to the human branch ad which species share nodes. This file must bein descending distance to human tarting wit the mostdistant node.

calculateAncestors.csh
This program takes in the final alignments, the regions of interest and the allSpecies.bed file. It also requires the speciesNode.list file. This outputs dated regions into node specific files in a directory called coverageCalls/ and will list the original coordinates of each region in its file.

species.list
contains the assembly name of every species included (order isn't important)

dists.comma.txt
contains the comma separated version of dists.txt

dists.txt
describes the phylogenetic distance between all pairwise sets of species in the newick. Calculated using PHAST/all_dists: http://compgen.cshl.edu/phast/

axtToNetBed.csh
this function processes axt files output by lastz into final BED files for input into calculateAncestors.csh. It requires an installation of kentutils and that the pairwise alignment outputs and other required files are available in the correct locations (see paths in script for locations required). This also requires alignment matrices.

hoxD55.mat
default.mat
human_chimp_v2.mat
these are alignment matrices needed for lastz and for chaining and neeting in netToAxtBed.csh


