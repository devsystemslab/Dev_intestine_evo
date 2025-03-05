#!/bin/csh -efv

set kentutils=<location of kentutils functions>
set axtToPsl=$kentutils/axtToPsl
set chainToPsl=$kentutils/chainToPslBasic
set pslToBed=$kentutils/pslToBed
set netToBed=$kentutils/netToBed
set axtChain=$kentutils/axtChain
set chainSort=$kentutils/chainSort
set chainPreNet=$kentutils/chainPreNet
set chainFilter=$kentutils/chainFilter
set chainNet=$kentutils/chainNet
set netSyntenic=$kentutils/netSyntenic
set netFilter=$kentutils/netFilter
set netClass=$kentutils/netClass
set netToBed=$kentutils/netToBed
set netToAxt=$kentutils/netToAxt
set temp=temp/
set lowerLim=0.2
set upperLim=0.7
set genomes=<location of chrom.sizes and 2bit files>

foreach ref (hg38)
set tSize=$genomes/$ref/$ref.chrom.sizes
set tBit=$genomes/$ref/$ref.2bit

echo $ref

foreach species (`cat species.list`)
set qSize=$genomes/$species/$species.chrom.sizes
set qBit=$genomes/$species/$species.2bit

echo $species

if (`echo $species` == `echo $ref`) then
continue
endif

foreach dist (`cat dists.comma.txt`)
set word=($dist:as/,/ /)

if ($word[1] != "(total)") then

if (`echo $ref` == `echo $species`) then
set matrix="none"

else if (`echo $species` == `echo $word[1]` && `echo $ref` == `echo $word[2]`) then
set outDir=chainingNetting/$ref.$species

set actual=$word[3]

echo $actual

set resultClose=`awk -v lowerLim="$lowerLim" -v upperLim="$upperLim" -v actual="$actual" 'BEGIN { if (actual <= lowerLim) print "true"; else print "false" }'`
set resultDefault=`awk -v lowerLim="$lowerLim" -v upperLim="$upperLim" -v actual="$actual" 'BEGIN { if (actual >= lowerLim && actual <= upperLim) print "true"; else print "false" }'`
echo $resultClose
echo $resultDefault

if ($resultClose == "true") then
set matrix=human_chimp_v2.mat
else if ($resultsDefault == "true") then
set matrix=default.mat
else
set matrix=hoxD55.mat
endif

else if (`echo $species` == `echo $word[2]` && `echo $ref` == `echo $word[1]`) then
set outDir=chainingNetting/$ref.$species

set actual=$word[3]

echo $actual

set resultClose=`awk -v lowerLim="$lowerLim" -v upperLim="$upperLim" -v actual="$actual" 'BEGIN { if (actual <= lowerLim) print "true"; else print "false" }'`
set resultDefault=`awk -v lowerLim="$lowerLim" -v upperLim="$upperLim" -v actual="$actual" 'BEGIN { if (actual >= lowerLim && actual <= upperLim) print "true"; else print "false" }'`
echo $resultClose
echo $resultDefault


if ($resultClose == "true") then
set matrix=human_chimp_v2.mat
else if ($resultDefault == "true") then
set matrix=default.mat
else
set matrix=hoxD55.mat

endif

else
continue

endif

echo $matrix

if ($matrix:t == hoxD55.mat) then
set linearGap="loose"
set minSynScore=50000
set minChainScore=20000
#I'm guessing at some decent values, I'm going to leave minSynSize at default value for all

foreach file (`find chainingNetting/$ref.$species -name "*.axt"`)
set name=$file:t:r
set tName=$file:t:r:r
set outChain=$temp/$name.filteredScore.chain
set syntenic=$temp/$name.unfilteredSyntenic.net
set outSynTarget=$outDir/$name.filteredSyn.net
set outPsl=$temp/$name.filteredNet.psl
set bed=$outDir/$name.syntenic.bed
set outAxt=$temp/$name.filteredNet.axt

$axtChain -linearGap=$linearGap -scoreScheme=$matrix $file $tBit $qBit stdout \
| $chainSort stdin stdout \
| $chainPreNet stdin $tSize $qSize stdout \
| $chainFilter -minScore=$minChainScore stdin > $outChain
$chainNet $outChain $tSize $qSize stdout /dev/null \
| $netSyntenic stdin stdout > $syntenic
$netFilter -syn -minSynScore=$minSynScore $syntenic > $outSynTarget
$netToBed -maxGap=0 $outSynTarget $bed
$netToAxt $outSynTarget $outChain $tBit $qBit $outAxt
$axtToPsl $outAxt $tSize $qSize $outPsl

end

else if ($matrix:t == human_chimp_v2.mat) then
set linearGap="medium"
set minChainScore=1000000

foreach file (`find chainingNetting/$ref.$species -name "*.axt"`)
set name=$file:t:r
set tName=$file:t:r:r
set outChain=$temp/$name.filteredScore.chain
set syntenic=$temp/$name.unfilteredSyntenic.net
set outSynTarget=$outDir/$name.filteredSyn.net
set outPsl=$temp/$name.filteredNet.psl
set bed=$outDir/$name.syntenic.bed
set outAxt=$temp/$name.filteredNet.axt

$axtChain -linearGap=$linearGap -scoreScheme=$matrix $file $tBit $qBit stdout \
| $chainSort stdin stdout \
| $chainPreNet stdin $tSize $qSize stdout \
| $chainFilter -minScore=$minChainScore stdin > $outChain
$chainNet $outChain $tSize $qSize stdout /dev/null \
| $netSyntenic stdin stdout > $syntenic
$netFilter -chimpSyn $syntenic > $outSynTarget
$netToBed -maxGap=0 $outSynTarget $bed
$netToAxt $outSynTarget $outChain $tBit $qBit $outAxt
$axtToPsl $outAxt $tSize $qSize $outPsl

end

else if ($matrix:t == default.mat) then
set linearGap="medium"
set minChainScore=100000
foreach file (`find chainingNetting/$ref.$species -name "*.axt"`)
set name=$file:t:r
set tName=$file:t:r:r
set outChain=$temp/$name.filteredScore.chain
set syntenic=$temp/$name.unfilteredSyntenic.net
set outSynTarget=$outDir/$name.filteredSyn.net
set outPsl=$temp/$name.filteredNet.psl
set bed=$outDir/$name.syntenic.bed
set outAxt=$temp/$name.filteredNet.axt

$axtChain -linearGap=$linearGap -scoreScheme=$matrix $file $tBit $qBit stdout \
| $chainSort stdin stdout \
| $chainPreNet stdin $tSize $qSize stdout \
| $chainFilter -minScore=$minChainScore stdin > $outChain
$chainNet $outChain $tSize $qSize stdout /dev/null \
| $netSyntenic stdin stdout > $syntenic
$netFilter -syn $syntenic > $outSynTarget
$netToBed -maxGap=0 $outSynTarget $bed
$netToAxt $outSynTarget $outChain $tBit $qBit $outAxt
$axtToPsl $outAxt $tSize $qSize $outPsl

end

else
continue
endif
endif
end
end
end
echo "DONE"

