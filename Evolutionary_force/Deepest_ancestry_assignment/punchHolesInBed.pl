#!/usr/bin/env perl

use Getopt::Long;
use warnings FATAL => 'all';
use strict;

#hash
my %chromHash = ();
my %startHash = ();
my %endHash = ();
my %blockCountHash = ();
my %blockSizesHash = ();
my %blockStartsHash = ();

#options
my $optVerbose;


sub usage
{
    print STDERR "
usage:   punchHolesInBed.pl overlapSelectMergedOutput.tsv out.bed 
options:
     -verbose
";
    exit(1);
}

sub checkOptions
{
    # Make sure command line options are valid/supported.
    my $ok = GetOptions("verbose" => \$optVerbose);
    if(! $ok){&usage();}
    if(!(defined $optVerbose)){$optVerbose = 0;}

    if($optVerbose){print(STDERR "Options look ok\n");}
}

sub abort
{
    my ($message) = @_;
    print(STDERR "$message\n");
    exit(1);
}

sub min
{
	my ($a,$b) = @_;
	if($a <= $b){return($a);}
	else{return($b);}
}

sub max
{
        my ($a,$b) = @_;
        if($a >= $b){return($a);}
        else{return($b);}
}

sub readData
{
	my ($inFile) = @_;

	my ($size,$start,$end,$inChrom,$inStart,$inEnd,$inName,$selectChrom,$selectStart,$selectEnd,$selectName,$line);

	open(Fin, "$inFile") || &abort("$inFile could not be opened");
	while(defined($line = <Fin>))
	{
		chomp($line);
		($inChrom,$inStart,$inEnd,$inName,$selectChrom,$selectStart,$selectEnd,$selectName) = split('\t',$line);
		#print("inChrom = $inChrom\n");
		#print("inStart = $inStart\n");
		#print("inEnd = $inEnd\n");
		#print("inName = $inName\n");
		#print("selectChrom = $selectChrom\n");
		#print("selectStart = $selectStart\n");
		#print("selectEnd = $selectEnd\n");
		#print("selectName = $selectName\n");
		$start = &max($inStart,$selectStart);
		$end = &min($inEnd,$selectEnd);
		$size = $end - $start;
		if(!(defined($chromHash{$inName}))){$chromHash{$inName} = $inChrom;}
		elsif($inChrom ne $chromHash{$inName}){&abort("chroms don't match: $line");}
		if(!(defined($startHash{$inName}))){$startHash{$inName} = $start;}
		else{$startHash{$inName} = &min($start,$startHash{$inName});}
		if(!(defined($endHash{$inName}))){$endHash{$inName} = $end;}
		else{$endHash{$inName} = &max($end,$endHash{$inName});}
		if(!(defined($blockCountHash{$inName}))){$blockCountHash{$inName} = 1;}
                else{$blockCountHash{$inName}+=1;}
		if(!(defined($blockSizesHash{$inName}))){$blockSizesHash{$inName} = "$size";}
                else{$blockSizesHash{$inName} = "$blockSizesHash{$inName},$size";}
		if(!(defined($blockStartsHash{$inName}))){$blockStartsHash{$inName} = "$start";}
                else{$blockStartsHash{$inName} = "$blockStartsHash{$inName},$start";}
	}
	close(Fin);
}

sub adjustBlockStarts
{
	my ($start,$stringList) = @_;

	my $block;
	my @blocks = split(",",$stringList);
	foreach $block (@blocks)
	{
		$block -= $start;
	}
	return(join(",",@blocks));
}

sub writeData
{
        my ($outFile) = @_;

        my ($name);

	open(Fout, ">$outFile") || &abort("could not open $outFile");
	foreach $name (keys(%chromHash))
	{
		$blockStartsHash{$name} = &adjustBlockStarts($startHash{$name},$blockStartsHash{$name});
		print(Fout "$chromHash{$name}\t$startHash{$name}\t$endHash{$name}\t$name\t0\t+\t$startHash{$name}\t$endHash{$name}\t255,255,255\t$blockCountHash{$name}\t$blockSizesHash{$name}\t$blockStartsHash{$name}\n");
	}
	close(Fout);
}

######################################################
# main
######################################################


&checkOptions();

if (($#ARGV + 1) != 2){&usage();}

my $inFile = $ARGV[0];
my $outFile = $ARGV[1];

&readData($inFile);
&writeData($outFile);


