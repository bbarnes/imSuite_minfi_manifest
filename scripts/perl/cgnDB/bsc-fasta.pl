#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);
use File::Basename;
use PerlIO::gzip;
use Compress::Zlib;
use Storable;
use IPC::Open2;
use FileHandle;

#use Data::Dumper::Simple;
#use Compress::BGZF::Writer;
#use Compress::BGZF::Reader;

my @date = ();
my @runTimes = ();
my ($endTime, $runTime, $curTime, $curDate);
my $startTime=time; $endTime=$startTime;
use constant FALSE => 0;
use constant TRUE  => 1;
use constant WIDTH => 2;

## improbe Format:
use constant SEQ_ID                              => 0;
use constant FORWARD_SEQUENCE                    => 1;
use constant GENOME_BUILD                        => 2;
use constant CHROMOSOME                          => 3;
use constant COORDINATE                          => 4;
use constant DESIGN_STATE                        => 5;
use constant SEQ_LENGTH                          => 6;
use constant FORWARD_CPG_COORD                   => 7;
use constant TB_STRAND                           => 8;
use constant TOP_SEQUENCE                        => 9;
use constant TOP_CPG_COORD                       => 10;
use constant PROBE_TYPE                          => 11;
use constant PROBESET_ID                         => 12;
use constant PROBESET_SCORE                      => 13;
use constant METHYL_PROBE_ID                     => 14;
use constant METHYL_PROBE_SEQUENCE               => 15;
use constant METHYL_PROBE_LENGTH                 => 16;
use constant METHYL_START_COORD                  => 17;
use constant METHYL_END_COORD                    => 18;
use constant METHYL_PROBE_COVERED_TOP_SEQUENCE   => 19;
use constant METHYL_ALLELE_FR_STRAND             => 20;
use constant METHYL_ALLELE_TB_STRAND             => 21;
use constant METHYL_ALLELE_CO_STRAND             => 22;
use constant METHYL_ALLELE_TYPE                  => 23;
use constant METHYL_FINAL_SCORE                  => 24;
use constant METHYL_TM                           => 25;
use constant METHYL_TM_SCORE                     => 26;
use constant METHYL_GC_PERCENT                   => 27;
use constant METHYL_GC_SCORE                     => 28;
use constant METHYL_13MER_COUNT                  => 29;
use constant METHYL_13MER_SCORE                  => 30;
use constant METHYL_ADDRESS_COUNT                => 31;
use constant METHYL_ADDRESS_SCORE                => 32;
use constant METHYL_SELF_COMPLEMENTARITY         => 33;
use constant METHYL_SELF_COMPLEMENTARITY_SCORE   => 34;
use constant METHYL_MONO_RUN                     => 35;
use constant METHYL_MONO_RUN_SCORE               => 36;
use constant METHYL_ECTOPIC_COUNT                => 37;
use constant METHYL_ECTOPIC_SCORE                => 38;
use constant METHYL_UNDERLYING_CPG_COUNT         => 39;
use constant METHYL_UNDERLYING_CPG_MIN_DIST      => 40;
use constant METHYL_UNDERLYING_CPG_SCORE         => 41;
use constant METHYL_IN_CPG_ISLAND_RELAXED        => 42;
use constant METHYL_CPG_ISLAND_SCORE             => 43;
use constant METHYL_NEXT_BASE                    => 44;
use constant METHYL_NEXT_BASE_SCORE              => 45;
use constant UNMETHYL_PROBE_ID                   => 46;
use constant UNMETHYL_PROBE_SEQUENCE             => 47;
use constant UNMETHYL_PROBE_LENGTH               => 48;
use constant UNMETHYL_START_COORD                => 49;
use constant UNMETHYL_END_COORD                  => 50;
use constant UNMETHYL_PROBE_COVERED_TOP_SEQUENCE => 51;
use constant UNMETHYL_ALLELE_FR_STRAND           => 52;
use constant UNMETHYL_ALLELE_TB_STRAND           => 53;
use constant UNMETHYL_ALLELE_CO_STRAND           => 54;
use constant UNMETHYL_ALLELE_TYPE                => 55;
use constant UNMETHYL_FINAL_SCORE                => 56;
use constant UNMETHYL_TM                         => 57;
use constant UNMETHYL_TM_SCORE                   => 58;
use constant UNMETHYL_GC_PERCENT                 => 59;
use constant UNMETHYL_GC_SCORE                   => 60;
use constant UNMETHYL_13MER_COUNT                => 61;
use constant UNMETHYL_13MER_SCORE                => 62;
use constant UNMETHYL_ADDRESS_COUNT              => 63;
use constant UNMETHYL_ADDRESS_SCORE              => 64;
use constant UNMETHYL_SELF_COMPLEMENTARITY       => 65;
use constant UNMETHYL_SELF_COMPLEMENTARITY_SCORE => 66;
use constant UNMETHYL_MONO_RUN                   => 67;
use constant UNMETHYL_MONO_RUN_SCORE             => 68;
use constant UNMETHYL_ECTOPIC_COUNT              => 69;
use constant UNMETHYL_ECTOPIC_SCORE              => 70;
use constant UNMETHYL_UNDERLYING_CPG_COUNT       => 71;
use constant UNMETHYL_UNDERLYING_CPG_MIN_DIST    => 72;
use constant UNMETHYL_UNDERLYING_CPG_SCORE       => 73;
use constant UNMETHYL_IN_CPG_ISLAND_RELAXED      => 74;
use constant UNMETHYL_CPG_ISLAND_SCORE           => 75;
use constant UNMETHYL_NEXT_BASE                  => 76;
use constant UNMETHYL_NEXT_BASE_SCORE            => 77;


##
## Common Methylation && Format Variables:
##
my ($II, $TOP, $BOT, $UNK, $GRN, $RED, $NEG);
my (@TB, @FR, @CO, @III, @MU, @RG, @AB, @ACTG, %ACTG, @MUD);
my ($SL, $TAB, $COM, $COL, $RET, $SEM, $PLUS, $MINUS, $DASH, $DOT);
my ($CG, $TG, $UM, $MA, $FO);
my ($A, $B, $C, $D, $E, $F, $G, $H, $I, $J, $K, $L, $M, $N,
    $O, $P, $Q, $R, $S, $T, $U, $V, $W, $X, $Y, $Z);
my ($d);
my ($iU, $iM, $iD, $iN, $iF, $iR, $iT, $iB, $iC, $iO);
my ($iFC, $iFO, $iRC, $iRO, $iBC, $iBO, $iTC, $iTO);

my %base2color=();
my %expandBase=();
my @iSrd2bsp=();
my (@idx2mu, %mu2idx, %iSrd, %bases2iupac, %iupac2bases);
my (%mapM, %mapD, %bspFR, %bspCO, %bspTB);
my (%cmplNucs, %iMU, %baseIdx, @idxBase);
my (%iupacTuple, %iupacPair, %iupacNucs, %cmplChr);
my ($lev5, $forcastLim, $forcastFlag, $forcastCnt, $cpgPadLen);
my @iStrands=(0 .. 3);
my $iupacCG = "[CYSMBHV][GRSKBDV]";
my ($True,$False,$Ref,$Alt,$Red,$Grn,$Mix,$Nan,$chunkSize);

my $improbeExe  = "/illumina/scratch/darkmatter/Projects/NZT/multiUnique/bin/improbe";
my $bwbbleExe   = "/illumina/scratch/darkmatter/bscCGN/bin/bwbble";
my $bsmapExe    = "/illumina/scratch/darkmatter/Software/bsmap-2.90/bsmap";
my $ntthalExe   ="/Users/bbarnes/Documents/Projects/programs/primer3/src/ntthal";
my $ntthalParams=""; ## Use defaults...

my $impTango    = "/illumina/scratch/methylation/Magellan_2010/dat/Tango_A_or_B_11mer_s1.dat";
my $imp13mer    = "/illumina/scratch/methylation/Magellan_2010/dat/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat";
my $impExe      = "/illumina/scratch/darkmatter/Projects/NZT/multiUnique/bin/improbe";
my $impLog      = "/home/bbarnes/tmp/improbe.log";
my $impCmd      = "$impExe -oASPE -tBOTH -cBoth -n$imp13mer -a$impTango -V - 2>$impLog";

my $annHeader   = "Seq_ID,Forward_Sequence,Genome_Build,Chromosome,Coordinate,Design_State,Seq_Length,Forward_CpG_Coord,TB_Strand,Top_Sequence,Top_CpG_Coord,Probe_Type,Probeset_ID,Probeset_Score,Methyl_Probe_ID,Methyl_Probe_Sequence,Methyl_Probe_Length,Methyl_Start_Coord,Methyl_End_Coord,Methyl_Probe_Covered_Top_Sequence,Methyl_Allele_FR_Strand,Methyl_Allele_TB_Strand,Methyl_Allele_CO_Strand,Methyl_Allele_Type,Methyl_Final_Score,Methyl_Tm,Methyl_Tm_Score,Methyl_GC_Percent,Methyl_GC_Score,Methyl_13mer_Count,Methyl_13mer_Score,Methyl_Address_Count,Methyl_Address_Score,Methyl_Self_Complementarity,Methyl_Self_Complementarity_Score,Methyl_Mono_Run,Methyl_Mono_Run_Score,Methyl_Ectopic_Count,Methyl_Ectopic_Score,Methyl_Underlying_CpG_Count,Methyl_Underlying_CpG_Min_Dist,Methyl_Underlying_CpG_Score,Methyl_In_CpG_Island_Relaxed,Methyl_CpG_Island_Score,Methyl_Next_Base,Methyl_Next_Base_Score,UnMethyl_Probe_ID,UnMethyl_Probe_Sequence,UnMethyl_Probe_Length,UnMethyl_Start_Coord,UnMethyl_End_Coord,UnMethyl_Probe_Covered_Top_Sequence,UnMethyl_Allele_FR_Strand,UnMethyl_Allele_TB_Strand,UnMethyl_Allele_CO_Strand,UnMethyl_Allele_Type,UnMethyl_Final_Score,UnMethyl_Tm,UnMethyl_Tm_Score,UnMethyl_GC_Percent,UnMethyl_GC_Score,UnMethyl_13mer_Count,UnMethyl_13mer_Score,UnMethyl_Address_Count,UnMethyl_Address_Score,UnMethyl_Self_Complementarity,UnMethyl_Self_Complementarity_Score,UnMethyl_Mono_Run,UnMethyl_Mono_Run_Score,UnMethyl_Ectopic_Count,UnMethyl_Ectopic_Score,UnMethyl_Underlying_CpG_Count,UnMethyl_Underlying_CpG_Min_Dist,UnMethyl_Underlying_CpG_Score,UnMethyl_In_CpG_Island_Relaxed,UnMethyl_CpG_Island_Score,UnMethyl_Next_Base,UnMethyl_Next_Base_Score,UnMethyl_assay_id";

##
## Basic Program Variables
##
my ($topDir, $statusCnt, $qsub, $qsubCnt, $local,
    $clean, $verbose, $silent, $mode, $usage, $submit);

my $validate=TRUE;
#$validate=FALSE;

# Set Global Parameter::
#
setParams();

$mode="parse";

##
## Input Variables and Data Structures
##
my %opts=();
my %stats = ();

$opts{'srcMax'}=0;
$opts{'chrMax'}=0;
$opts{'lenMax'}=0;

# FAI Fasta Index Data Structures::
#
my %faiRef=();
my %fai2idx=();
my %idx2fai=();
my %isGenomic=();

# CGN Database Data Structures::
#
# my %top2cgn=();
# my @cgn2top=();

#
#
# Missing Global Data Structures::
#
#

# Input Variables and Data Structures:: dbSNP()
my %usrDef = ();
my @srcNames = ();
my @bedFiles = (); my $bedMax; my @bedSrc = ();
my @manFiles = (); my $manMax; my @manSrc = ();
my @desFiles = (); my $datMax; my @datSrc = ();
my @fasFiles = (); my $fasMax; my @fasSrc = ();
my @vcfFiles = (); my $vcfMax;
my @binFiles = ();

my %impFhs = (); my %binFhs = (); my %snpFhs = (); my %bedFhs = ();
my %impFns = (); my %binFns = (); my %snpFns = (); my %bedFns = ();

my %fasFhsFCU = (); my %fasFhsFCM = (); my %fasFhsFCD = ();
my %fasFnsFCU = (); my %fasFnsFCM = (); my %fasFnsFCD = ();

my %fasFhsRCU = (); my %fasFhsRCM = (); my %fasFhsRCD = ();
my %fasFnsRCU = (); my %fasFnsRCM = (); my %fasFnsRCD = ();

my ($outName);

my %snpCgns = ();
my %srcCgns = ();
my %cgnDB   = ();
my %cgnQC   = ();
my $cpgIdx  = 1;
my $cgnMax  = $cpgIdx;
my $srcLenMax = 6;

#
# NOTE: These steps are no logner needed, but should be tracked...
#
# my @badCGNs = qw(3458191 11037148 11718090 15092802 15408454 16390856 23509027);
# foreach my $cgn (@badCGNs) { $cgnDB{$cgn} = $cgn; }
#
# my $nnnSeq = "TTATTTATTAGCCACAAGGGAACTCTTTTTTTTCAAGTTCTTAATCAGAGCACTGGTCATCGTTCCCTGGAGGTGAATCCTGATTATTCATAAGACAAACCTGAATTCNNNNNNNNNNNNNN";
# $nnnSeq =~ s/N/A/gi;
# $cgnDB{"13279585"} = $nnnSeq;
# $cgnDB{$nnnSeq} = "13279585";
#
# my ($nnnCgn, $nnnSeq);
# # 718858; umanMethylation450_15017482_v.1.2
# $nnnCgn = 718858;
# $nnnSeq = "NATTTTAAAATACCCAGCTCCACCCCTTCCTGTTAGGCTTTCGCGTGTCGCAGCTGTGCACGCTGATTGGTCCTCTGCTGGCCAATCACCACTGCACTTCATGACGGCTGTAGTTTTCAAAA";
# $nnnSeq =~ s/N/A/gi;
# $cgnDB{$nnnCgn} = $nnnSeq;
# $cgnDB{$nnnSeq} = $nnnCgn;
#
# END NOTE
#

my $DB_NEW="new"; my $DB_MAT="mat"; my $DB_NNN="nnn";

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Program Functions:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Program Main Subroutines:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub isPoorTemplate {
    my $func=(caller(0))[3];
    my ($seq, $dieOnLen) = @_;

    my $len = length($seq);
    die("[ERROR]: SeqLen($len) != 122! Seq=Seq\n")
        if ($dieOnLen && $len != 122);

    return(TRUE) if ($seq =~ m/[^ACTG]/ || $len != 122);
    return(FALSE);
}

sub addSeq($) {
    my $func=(caller(0))[3];
    my ($fwdSeq) = @_;

    my ($topSeq, $botSeq, $topTB, $botTB) = returnTopBot($fwdSeq);

    my ($cpg, $isNew);
    if ($topTB eq $S && $botTB eq $S) {
        if (defined($cgnDB{uc($topSeq)})) {
            die("[$func]: ERROR: $S Failed to match mate:\n".
		"[$func]: ERROR: fwd=".addBrac($fwdSeq)."\n".
		"[$func]: ERROR: top=".addBrac($topSeq)."\n".
		"[$func]: ERROR: bot=".addBrac($botSeq)."\n") if (! defined($cgnDB{uc($botSeq)}));
        } else {
            $cgnDB{uc($topSeq)} = $botSeq;
            $cgnDB{uc($botSeq)} = $topSeq;
        }
        $cpg = 0;
    } elsif (defined($cgnDB{uc($topSeq)})) {
        $cpg = $cgnDB{uc($topSeq)};
        die("[$func]: ERROR: top defined and bot defined! cpg=$cpg, cpg(bot)=$cgnDB{uc($botSeq)}\n".
            "[$func]: ERROR: fwd=".addBrac($fwdSeq)."\n".
            "[$func]: ERROR: top=".addBrac($topSeq)."\n".
            "[$func]: ERROR: bot=".addBrac($botSeq)."\n") if ($topTB ne $U && defined($cgnDB{uc($botSeq)}));
        die("[$func]: ERROR: Reported topSeq doesn't match lookup cpg=$cpg\n".
            "[$func]: ERROR: cpg=".addBrac($cgnDB{$cpg})."\n".
            "[$func]: ERROR: fwd=".addBrac($fwdSeq)."\n".
            "[$func]: ERROR: top=".addBrac($topSeq)."\n".
            "[$func]: ERROR: bot=".addBrac($botSeq)."\n") if (uc($cgnDB{$cpg}) ne uc($topSeq));
        $isNew = FALSE;
        $cgnQC{$DB_MAT}++;
    } else {
        die("[$func]: ERROR: top defined and bot defined! cpg=$cpg, cpg(bot)=$cgnDB{uc($botSeq)}\n".
            "[$func]: ERROR: fwd=".addBrac($fwdSeq)."\n".
            "[$func]: ERROR: top=".addBrac($topSeq)."\n".
            "[$func]: ERROR: bot=".addBrac($botSeq)."\n") if (defined($cgnDB{uc($botSeq)}));

        $cpgIdx++ while(defined($cgnDB{$cpgIdx}));

	$cgnDB{uc($topSeq)} = $cpgIdx;
        $cgnDB{$cpgIdx} = $topSeq;
        $isNew = TRUE;
        $cgnQC{$DB_NEW}++;
        $cpg = $cpgIdx;
    }
    $cgnMax = $cpg if ($cpg > $cgnMax);
    return($cpg, $topSeq, $topTB, $isNew);
}

sub addCpgSeq($$) {
    my $func=(caller(0))[3];
    my ($cpg, $fwdSeq) = @_;

    my ($topSeq, $botSeq, $topTB, $botTB) = returnTopBot($fwdSeq);

    my ($isNew);
    if (defined($cgnDB{uc($topSeq)})) {
        die("[$func]: ERROR: top defined and cpg is NOT defined! cpg=$cpg, cpg(top)=$cgnDB{uc($topSeq)}\n".
            "[$func]: ERROR: fwd=$fwdSeq\n".
            "[$func]: ERROR: top=$topSeq\n".
            "[$func]: ERROR: cpg=$cgnDB{$cpg}\n".
            "[$func]: ERROR: bot=$botSeq\n")
	    if (! defined($cgnDB{$cpg}) || $cpg ne $cgnDB{uc($topSeq)} ||
		$cgnDB{$cpg} ne $topSeq);
	
        die("[$func]: ERROR: top defined and bot defined! cpg=$cpg, cpg(bot)=$cgnDB{uc($botSeq)}\n".
            "[$func]: ERROR: fwd=$fwdSeq\n".
            "[$func]: ERROR: top=$topSeq\n".
            "[$func]: ERROR: bot=$botSeq\n")
	    if (defined($cgnDB{uc($botSeq)}));
	
        $isNew = FALSE;
        $cgnQC{$DB_MAT}++;
    } else {
        die("[$func]: ERROR: top NOT defined, but cpg is! cpg=$cpg\n".
            "[$func]: ERROR: fwd=$fwdSeq\n".
            "[$func]: ERROR: top=$topSeq\n".
            "[$func]: ERROR: cpg=$cgnDB{$cpg}\n".
            "[$func]: ERROR: bot=$botSeq\n")
	    if (defined($cgnDB{$cpg}));

        ## Make assignment
        $cgnDB{uc($topSeq)} = $cpg;
        $cgnDB{$cpg} = $topSeq;
        $isNew = TRUE;
        $cgnQC{$DB_NEW}++;
    }
    $cgnMax = $cpg if ($cpg > $cgnMax);
    return($cpg, $topSeq, $topTB, $isNew);
}

sub parseCpg($$$) {
    my $func=(caller(0))[3];
    my ($chr, $chrSeq, $source) = @_;

    my $chrLen = length($chrSeq);
    my $chrTen = int($chrLen/1000);
    $chrTen++ if ($chrTen==0);
    
    mssg("[$func]: chr=$chr, seqLen=$chrLen; chrTen=$chrTen\n",3);
    # $chrSeq .= 'N' x 60;
    $chrSeq .= 'N' x 124;

    my $isSNP = FALSE;
    $isSNP = TRUE if ($source =~ m/snp/);

    my $iupacCG = "[CYSMBHV][GRSKBDV]";
    my $chrStr = "chr$chr";
    my $posIdx;
    my $iter = 0;
    mssg("[$func]: Starting: chr=$chr, iter=$iter\n");
    while($chrSeq =~ m/[CYSMBHV][GRSKBDV]/g) {
	$posIdx = pos($chrSeq);

	if ($posIdx < 60) {
	    warn("[$func]: WARNING: posIdx=$posIdx < 60; Skipping...\n",3);
	    next;
	}

	my $fwdSeq = substr($chrSeq, $posIdx-62, 122);
	my $midFwd = substr($fwdSeq, 60, 2);
        my $seqLen = length($fwdSeq);
	my $seqSrd = "+";

	if ($seqLen != 122) {
	    mssg("[$func]: Warning: chr=$chr; posIdx=$posIdx; seqLen=$seqLen != 122; idx=$posIdx; ".
		 "fwdSeq=$fwdSeq, chrLen=$chrLen; Skipping...\n\n");
	    # print "\n\n";
	    # print "CHR_SEQ='$chrSeq'\n\n";
	    next;
	}
        # die("[$func]: ERROR: posIdx=$posIdx; seqLen=$seqLen != 122; idx=$posIdx; ".
	#     "fwdSeq=$fwdSeq, chrLen=$chrLen\n") if ($seqLen != 122);
		
	if ($fwdSeq =~ m/N/) {
	    $cgnDB{$DB_NNN}++;
	    #mssg("[NNN]: $fwdSeq\n\n");
	    #pos($chrSeq) = $posIdx+2;
	    next;
	}

	my $pos = $posIdx-1;
	my @expSeqs = expandSeq($fwdSeq);
	if ($isSNP && @expSeqs == 1) {
	    my $revSeq = revCmp($fwdSeq);
	    die("[$func]: ERROR: Non-SNP without definitions:\n".
		"[$func]: ERROR: fwd=".addBrac($fwdSeq)."\n".
		"[$func]: ERROR: rev=".addBrac($revSeq)."\n")
		if (! defined($cgnDB{$fwdSeq}) && ! defined($cgnDB{$revSeq}));
	    $cgnQC{"snpFree"}++;
	    next;
	}

	foreach my $expSeq (@expSeqs) {
	    my $midExp = substr($expSeq, 60, 2);
	    if (uc($midExp) ne "CG") {
		#mssg("[SkipFW]: $midFwd ".addBrac($fwdSeq)."\t$pos\n");
		#mssg("[SkipNN]: $midExp ".addBrac($expSeq)."\t$pos\n\n");
		$cgnQC{"nonCpG"}++;
		next;
	    }

	    my ($cpg, $topSeq, $topTB, $isNew) = addSeq($expSeq);
	    my $cpgStr = padCGN($cpg);
	    $seqSrd = "-" if ($fwdSeq eq $topSeq);
	    
	    #mssg("[Adding]: $midFwd ".addBrac($fwdSeq)."\t$pos\t$cpg\n");
	    #mssg("[Adding]: $midExp ".addBrac($topSeq)."\t$pos\t$cpg\n");
	    if ($topTB ne $S) {
		# OLD CODE that writes TOP not FWD to improbe input::
		# print {$impFhs{$source}} join($TAB, padCGN($cpg), addBrac($topSeq), $source, $chr, $pos, "FALSE").$RET;

		# Should Write FWD Sequence to improbe input and TOP Sequence to BED OR not at all...
		#  - Adding Top Strand (+/-) to Genome Build, this can be used for verification after running improbe...
		#    If seqSrd == + then fwd is top, if seqSrd == - then rev is top
		print {$impFhs{$source}} join($TAB, $cpgStr, addBrac($fwdSeq), $source.$seqSrd, $chr, $pos, "FALSE").$RET;

		# Removing TOP Sequence from BED file::
		# print {$bedFhs{$source}} join($TAB, $chrStr, $pos-1, $pos, $cpgStr.$L.$topTB, $source, $seqSrd, $topSeq).$RET;
		print {$bedFhs{$source}} join($TAB, $chrStr, $pos-1, $pos, $cpg, $source, $seqSrd).$RET;
		
		$srcCgns{$cpg}{$source}++;
		$snpCgns{$cpg}{$fwdSeq}++ if ($isSNP);
	    } else {
		#print CGN_S padCGN($cpg).$TAB.addBrac($fwdSeq).$TAB.$topTB.$RET;
	    }
	}
	mssg("[$func]:\t Current posIdx=$pos; iter=$iter\n") if ($iter % $chrTen==0);

	$iter++;
    }
    mssg("[$func]: cgnDB=".scalar(keys %cgnDB)."\n",3);
    mssg("[$func]: cgnQC=\n".hashStr(\%cgnQC)."\n",3);
    mssg("\n",3);
}

sub loadFasta($$) {
    my $func=(caller(0))[3];
    my ($file, $source) = @_;

    die("[$func]: ERROR: loadFasta: file or source are not defined!\n") if (! defined($file) || ! defined($source));
    die("[$func]: ERROR: loadFasta: Unable to find ($source) file=$file!\n") if (! -e $file);

    mssg("[$func]: src=$source\n",2);
    mssg("[$func]: fas=$file\n\n",2);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    trimSuffix(\$fileName);

    if ($fileSuffix eq ".gz") { open(FAS, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(FAS, "<$file") or die("Failed to open $file: $!") }
    mssg("[$func]: Loading genome ($source) fasta=$file...\n", 3);

    my $fasCnt = 0;
    my $nucTotCnt = 0;

    my $tag = "";
    $tag=$opts{'runName'} if (defined($opts{'runName'}));

    #
    # Define all Forward (F) Converted (C) BSC (U,M,D) Files::
    #
    my $fasFileFCU = "$opts{'outDir'}/$source.genome$tag.FCU.fa.gz";
    if (! defined($fasFhsFCU{$source})) {
        $fasFnsFCU{$source} = $fasFileFCU;
        open($fasFhsFCU{$source}, ">:gzip", $fasFileFCU) or die("Failed to open $fasFileFCU file for writing: $!");
    }
    my $fasFileFCM = "$opts{'outDir'}/$source.genome$tag.FCM.fa.gz";
    if (! defined($fasFhsFCM{$source})) {
        $fasFnsFCM{$source} = $fasFileFCM;
        open($fasFhsFCM{$source}, ">:gzip", $fasFileFCM) or die("Failed to open $fasFileFCM file for writing: $!");
    }
    my $fasFileFCD = "$opts{'outDir'}/$source.genome$tag.FCD.fa.gz";
    if (! defined($fasFhsFCD{$source})) {
        $fasFnsFCD{$source} = $fasFileFCD;
        open($fasFhsFCD{$source}, ">:gzip", $fasFileFCD) or die("Failed to open $fasFileFCD file for writing: $!");
    }

    #
    # Define all Reverse (R) Converted (C) BSC (U,M,D) Files::
    #
    my $fasFileRCU = "$opts{'outDir'}/$source.genome$tag.RCU.fa.gz";
    if (! defined($fasFhsRCU{$source})) {
        $fasFnsRCU{$source} = $fasFileRCU;
        open($fasFhsRCU{$source}, ">:gzip", $fasFileRCU) or die("Failed to open $fasFileRCU file for writing: $!");
    }
    my $fasFileRCM = "$opts{'outDir'}/$source.genome$tag.RCM.fa.gz";
    if (! defined($fasFhsRCM{$source})) {
        $fasFnsRCM{$source} = $fasFileRCM;
        open($fasFhsRCM{$source}, ">:gzip", $fasFileRCM) or die("Failed to open $fasFileRCM file for writing: $!");
    }
    my $fasFileRCD = "$opts{'outDir'}/$source.genome$tag.RCD.fa.gz";
    if (! defined($fasFhsRCD{$source})) {
        $fasFnsRCD{$source} = $fasFileRCD;
        open($fasFhsRCD{$source}, ">:gzip", $fasFileRCD) or die("Failed to open $fasFileRCD file for writing: $!");
    }

    #
    # These files below are not used::
    #
    my $impFile = "$opts{'outDir'}/$source.improbeDesignInput.tsv.gz";
    if (! defined($impFhs{$source})) {
        # $impFns{$source} = $impFile;
        # open($impFhs{$source}, ">:gzip", $impFile) or die("Failed to open $impFile file for writing: $!");
	# print {$impFhs{$source}} "Seq_ID\tSequence\tGenome_Build\tChromosome\tCoordinate\tCpG_Island\n";
    }

    my $bedFile = "$opts{'outDir'}/$source.cgn.bed.gz";
    if (! defined($bedFhs{$source})) {
	# $bedFns{$source} = $bedFile;
	# open($bedFhs{$source}, ">:gzip", $bedFile) or die("Failed to open $bedFile file for writing: $!");
    }
    
    my $curChr;
    my $curSeq = "";
    my $curLen = 0;
    my $seqLen = 0;
    my $chrCnt = 0;
    while(<FAS>) {
        s/^\s+//;
        s/\s+$//;

        my $line = $_;
        if ($line =~ m/^>([^\s]+)/) {
            if (defined($curChr)) {
		$curSeq=substr($curSeq,0,$opts{'lenMax'})
		    if (defined($opts{'lenMax'}) && $opts{'lenMax'} != 0);
		$curLen=length($curSeq);

		mssg("[$func]: Parsing CpGs chr=$curChr; curLen=$curLen...\n",3);

		# Convert/Write all Forward (F) Converted (C) BSC (U,M,D) Sequences::
		#
		my $curSeqU = $curSeq;
		my $curSeqM = $curSeq;
		my $curSeqD = $curSeq;
		
		bscU(\$curSeqU);
		print {$fasFhsFCU{$source}} ">$curChr\n$curSeqU\n";
		mssg("[$func]: Wrote bsc=FCU chr=$curChr, chrCnt=$chrCnt...\n",3);
		
		bscM(\$curSeqM);
		print {$fasFhsFCM{$source}} ">$curChr\n$curSeqM\n";
		mssg("[$func]: Wrote bsc=FCM chr=$curChr, chrCnt=$chrCnt...\n",3);
		
		bscD(\$curSeqD);
		print {$fasFhsFCD{$source}} ">$curChr\n$curSeqD\n";
		mssg("[$func]: Wrote bsc=FCD chr=$curChr, chrCnt=$chrCnt...\n",3);

		# Need to generate the reverse strand...
		# Convert/Write all Reverse (R) Converted (C) BSC (U,M,D) Sequences::
		#
		$curSeq = revCmpl($curSeq);
		$curSeqU = $curSeq;
		$curSeqM = $curSeq;
		$curSeqD = $curSeq;
		
		bscU(\$curSeqU);
		$curSeqU = revCmpl($curSeqU);
		print {$fasFhsRCU{$source}} ">$curChr\n$curSeqU\n";
		mssg("[$func]: Wrote bsc=RCU chr=$curChr, chrCnt=$chrCnt...\n",3);
		
		bscM(\$curSeqM);
		$curSeqM = revCmpl($curSeqM);
		print {$fasFhsRCM{$source}} ">$curChr\n$curSeqM\n";
		mssg("[$func]: Wrote bsc=RCM chr=$curChr, chrCnt=$chrCnt...\n",3);
		
		bscD(\$curSeqD);
		$curSeqD = revCmpl($curSeqD);
		print {$fasFhsRCD{$source}} ">$curChr\n$curSeqD\n";
		mssg("[$func]: Wrote bsc=RCD chr=$curChr, chrCnt=$chrCnt...\n",3);
		
            }
	    $chrCnt++;
            $curChr = $1;
            $curSeq = "";
	    $curLen = 0;
            $seqLen = 0;
            $curChr =~ s/^chr//;

	    last if (defined($opts{'chrMax'}) &&
		     $opts{'chrMax'} != 0 &&
		     $opts{'chrMax'} <= $chrCnt-1);
	    
            mssg("[$func]: Loading chr=$curChr, chrCnt=$chrCnt...\n",3);
        } else {
            $seqLen += length($line);
            $curSeq .= uc($line);
        }

        mssg("[$func]: statusCnt=$nucTotCnt, chr=$curChr, chrCnt=$chrCnt, seqlen=$seqLen...\n", 3)
            if (defined($statusCnt) && $nucTotCnt % $statusCnt == 0);
        $nucTotCnt++;
	
        # last if (defined($opts{'lenMax'}) &&
	# $opts{'lenMax'} != 0 &&
	# $opts{'lenMax'} <= $seqLen);
	# $opts{'lenMax'} <= $nucTotCnt);
    }
    close(FAS);
    
    if (defined($curChr)) {
	$curSeq=substr($curSeq,0,$opts{'lenMax'})
	    if (defined($opts{'lenMax'}) && $opts{'lenMax'} != 0);
	$curLen=length($curSeq);

        mssg("[$func]: Parsing CpGs chr=$curChr; curLen=$curLen...\n",3);

	# Convert/Write all Forward (F) Converted (C) BSC (U,M,D) Sequences::
	#
	my $curSeqU = $curSeq;
	my $curSeqM = $curSeq;
	my $curSeqD = $curSeq;
		
	bscU(\$curSeqU);
	print {$fasFhsFCU{$source}} ">$curChr\n$curSeqU\n";
	mssg("[$func]: Wrote bsc=FCU chr=$curChr, chrCnt=$chrCnt...\n",3);
		
	bscM(\$curSeqM);
	print {$fasFhsFCM{$source}} ">$curChr\n$curSeqM\n";
	mssg("[$func]: Wrote bsc=FCM chr=$curChr, chrCnt=$chrCnt...\n",3);
		
	bscD(\$curSeqD);
	print {$fasFhsFCD{$source}} ">$curChr\n$curSeqD\n";
	mssg("[$func]: Wrote bsc=FCD chr=$curChr, chrCnt=$chrCnt...\n",3);

	# Need to generate the reverse strand...
	# Convert/Write all Reverse (R) Converted (C) BSC (U,M,D) Sequences::
	#
	$curSeq = revCmpl($curSeq);
	$curSeqU = $curSeq;
	$curSeqM = $curSeq;
	$curSeqD = $curSeq;
		
	bscU(\$curSeqU);
	$curSeqU = revCmpl($curSeqU);
	print {$fasFhsRCU{$source}} ">$curChr\n$curSeqU\n";
	mssg("[$func]: Wrote bsc=RCU chr=$curChr, chrCnt=$chrCnt...\n",3);
		
	bscM(\$curSeqM);
	$curSeqM = revCmpl($curSeqM);
	print {$fasFhsRCM{$source}} ">$curChr\n$curSeqM\n";
	mssg("[$func]: Wrote bsc=RCM chr=$curChr, chrCnt=$chrCnt...\n",3);
		
	bscD(\$curSeqD);
	$curSeqD = revCmpl($curSeqD);
	print {$fasFhsRCD{$source}} ">$curChr\n$curSeqD\n";
	mssg("[$func]: Wrote bsc=RCD chr=$curChr, chrCnt=$chrCnt...\n",3);
    }
    close($fasFhsFCU{$source});
    close($fasFhsFCM{$source});
    close($fasFhsFCD{$source});

    close($fasFhsRCU{$source});
    close($fasFhsRCM{$source});
    close($fasFhsRCD{$source});

    
    # close($impFhs{$source});
    # close($bedFhs{$source});
    
    runTimes();
    mssg("[$func]: runTime=$runTimes[0]\n");
    mssg("[$func]: Loaded Genome Data=$nucTotCnt, candidatCnt=$fasCnt\n",2);
    mssg("[$func]: cgnDB=".scalar(keys %cgnDB)."\n",3);
    mssg("[$func]: cgnQC=\n".hashStr(\%cgnQC)."\n",3);
    mssg("\n",3);
}

sub loadBinary($) {
    my $func=(caller(0))[3];
    my ($file) = @_;

    die("[$func]: ERROR: loadBinary: Unable to find file=$file!\n") if (! -e $file);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    trimSuffix(\$fileName);

    #if ($fileSuffix eq ".gz") { open(BIN, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    #else { open(BIN, "<$file") or die("Failed to open $file: $!") }
    open(BIN, "<:raw","$file") or die("Failed to open $file for reading: $!");
    #binmode(BIN);
    mssg("[Binary]: Loading BinaryInput=$file...\n", 3);

    my $byteSize = 126;
    my ($buffer, $bytesRead);
    while( $bytesRead = read( BIN, $buffer, $byteSize ) ) {
	die("[$func]: ERROR: Unexpected bytesRead=$bytesRead, expecting $byteSize!\n") if ($bytesRead != $byteSize);

	my ($srcSeq, $srcCgn) = unpack("A122" ."N", $buffer);

	mssg("[$func]: $srcCgn, $srcSeq\n");
        my ($cgn, $topSeq, $topTB, $isNew) = addCpgSeq($srcCgn, $srcSeq);
	mssg("[$func]: $cgn, $topSeq, $topTB, $isNew\n\n");
    }
    close(BIN);

}

sub loadCanonical($) {
    my $func=(caller(0))[3];
    my ($inpFile) = @_;

    die("[$func]: ERROR: loadBinary: Unable to find file=$inpFile!\n") if (! -e $inpFile);

    # my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    # trimSuffix(\$fileName);

    my ($inpFH, $outFH);
    my ($inpFileName, $inpFilePath, $inpFileSuffix, $isBin)=openFile(\$inpFH,$inpFile);
    mssg("[$func]: Loading file=$inpFile...\n", 2);

    my %colToIdx=();
    
    setHeader(\%colToIdx, $inpFH, "CGN", $COM, 3);
    mssg("[$func]: Loading inpFile=$inpFile...\n", 3);
    
    my @lines=<$inpFH>;
    chomp(@lines);
    
    my ($totCnt, $tenCnt, $curCnt) = tenCnt(scalar(@lines));
    mssg("[$func]: TenCnt=$tenCnt, TotalCnt=$totCnt\n", 3);

    while(@lines) {
        my $line=shift(@lines);

	my ($srcCgn,$srcSeq,$srcGrp)=split(/,/, $line);

	my ($cgn, $topSeq, $topTB, $isNew) = addCpgSeq($srcCgn, $srcSeq);
	# mssg("[$func]: $cgn, $topSeq, $topTB, $isNew\n\n");
	
	if ($curCnt % $tenCnt == 0) {
            runTimes();
            mssg("[$func]: ".sprintf("%2d", 100*$curCnt/$totCnt).
		 "\% ($curCnt/$totCnt); cgn=$cgn; tb=$topTB; new=$isNew; ".
		 "RunTime=$runTimes[0]\n", 3);
        }
        $curCnt++;
        if (defined($opts{'canMax'}) && $curCnt >= $opts{'canMax'}) {
            mssg("[$func]: Readched Max $curCnt >= $opts{'canMax'}. Last...\n");
            last;
        }
    }
    close($inpFH);

    runTimes();
    mssg("[$func]: runTime=$runTimes[0]\n",2);
    # mssg("[$func]: Hack Percent = $hacPer\% ($hacCnt/$sumCnt)!"."\n", 2);    
    mssg("\n",2);
}

sub loadDesign($) {
    my $func=(caller(0))[3];
    my ($file) = @_;

    die("[$func]: ERROR: loadDesignInputs: Unable to find file=$file!\n") if (! -e $file);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    trimSuffix(\$fileName);

    if ($fileSuffix eq ".gz") { open(DES, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(DES, "<$file") or die("Failed to open $file: $!") }
    mssg("[$func]: Loading designInput=$file...\n", 3);

    my $desBinFile;

    my $desTotCnt = 0;
    my %printed = ();
    while(<DES>) {
        s/^\s+//;
        s/\s+$//;
        next if /^Seq_ID/;

        my $line = $_;
        my @data = split(/\t/, $line);

        my $orgCpg = $data[SEQ_ID];
        my $srcSeq = shearBrac($data[1]);
        my $build  = $data[2];
        my $chr    = $data[3];
        my $pos    = $data[4];
        $orgCpg =~ s/^cg//; $orgCpg =~ s/^0+//;

        if ($desTotCnt == 0 && $mode =~ m/w/) {
            $desBinFile = "$opts{'outDir'}/manifest.GRCh$build.cgn.bin";
	    push @binFiles, $desBinFile;
            open(BIN, ">$desBinFile") or die("Failed to open $desBinFile for writing: $!");
	    binmode(BIN);	    
	}
        if ($orgCpg eq "13279585") { $srcSeq =~ s/N/A/gi; print "$orgCpg\t$srcSeq\n"; }

        my ($cpg, $topSeq, $topTB, $isNew) = addCpgSeq($orgCpg, $srcSeq);

        if ($mode =~ m/w/) {
	    print BIN pack("A122", $topSeq);
	    print BIN pack("N", $cpg);
	    last if ($desTotCnt >= 3);
        }

	## Status:
        mssg("[$func]: statusCnt=$desTotCnt, curSeqCnt=".scalar(keys %cgnDB).", cpg=$cpg...\n", 3)
            if (defined($statusCnt) && $desTotCnt % $statusCnt == 0);
        $desTotCnt++;
        last if (defined($datMax) && $desTotCnt >= $datMax);
    }
    close(DES);
    close(BIN) if ($mode =~ m/w/);
    runTimes();
    mssg("[$func]: runTime=$runTimes[0]\n",2);
    mssg("[$func]: Total Designs Processed=$desTotCnt, lastCpgIdx=$cpgIdx\n",2);
    mssg("[$func]: cgnDB=".scalar(keys %cgnDB)."\n",3);
    mssg("[$func]: cgnQC=\n".hashStr(\%cgnQC)."\n",3);
    mssg("\n",2);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Function Calls:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub loadBins($) {
    my $func=(caller(0))[3];
    my ($fileRefs) = @_;

    foreach my $file (@{$fileRefs}) {
        loadBinary($file);
    }
}

sub loadCans($) {
    my $func=(caller(0))[3];
    my ($fileRefs) = @_;

    foreach my $file (@{$fileRefs}) {
        # loadBinary($file);
	loadCanonical($file)
    }
}

sub loadDesigns($) {
    my $func=(caller(0))[3];
    my ($fileRefs) = @_;

    foreach my $file (@{$fileRefs}) {
        loadDesign($file);
    }
}

sub loadFastas($$) {
    my $func=(caller(0))[3];
    my ($fileRefs, $nameRefs) = @_;

    ## Validate srcNames:
    for (my $ii=0; $ii<@{$fileRefs}; $ii++) {
	my $src = $$nameRefs[$ii];
	if (length($src) > $srcLenMax) {
	    die("[$func]: ERROR: Source Name=$src(".length($src).") is over limit($srcLenMax)!\n");
	} elsif (length($src) < $srcLenMax) {
	    $$nameRefs[$ii] .= "_" x ($srcLenMax - length($src));
	}

	last if (defined($opts{'srcMax'}) &&
		 $opts{'srcMax'} != 0 &&
		 $opts{'srcMax'} <= $ii+1);
    }

    for (my $ii=0; $ii<@{$fileRefs}; $ii++) {
	mssg("[$func]:\tCurrent src=$$nameRefs[$ii], fas=$$fileRefs[$ii]\n",2);
	
	loadFasta($$fileRefs[$ii], $$nameRefs[$ii]);

	last if (defined($opts{'srcMax'}) &&
		 $opts{'srcMax'} != 0 &&
		 $opts{'srcMax'} <= $ii+1);
    }

}


## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Main:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

GetOptions("o|out=s"      => \$opts{'outDir'},
	   
	   "s|src=s"      => \@{$opts{'srcName'}},
	   "f|fas=s"      => \@{$opts{'fasFile'}},
	   "c|can=s"      => \@{$opts{'canFile'}},
	   "d|des=s"      => \@{$opts{'desFile'}},
	   "b|bin=s"      => \@{$opts{'binFile'}},

	   "n|name=s"     => \$opts{'runName'},
	   
	   # "u|usr=s"      => \%usrDef,
	   
	   "m|man=s"      => \$opts{'manFile'},
	   "h|hac=s"      => \$opts{'hacFile'},

	   "srcMax=i"     => \$opts{'srcMax'},
	   "chrMax=i"     => \$opts{'chrMax'},
	   "lenMax=i"     => \$opts{'lenMax'},

	   "impExe=s"     => \$opts{'impExe'},
	   "impMer=s"     => \$opts{'impMer'},
	   "impTan=s"     => \$opts{'impTan'},
	   
	   "mod=s"        => \$mode,
           "vcfMax=i"     => \$vcfMax,
           "clean=i"      => \$clean,

           "v|verbose=i"  => \$verbose,
    );

@date = getTime(); runTimes();
mssg("\n[Starting]: $date[0], $date[2], runTime=$runTimes[0]\n\n");

$usage="[Usage]: $0 -o|out outputDir\n"
    . " -f|fas fasFile(s)\n"
    . " -s|src sourceNames(s)]\n"
    . " -d|des designInputFile(s) ]\n"
    . "\n"
    . " [ -mod mode ]\n"
    . " [ -srcMax max chromosome soruces to process ]\n"
    . " [ -chrMax max chromosomes to process ]\n"
    . " [ -lenMax max chromosome length to process ]\n"
    . "\n"
    . " [ -u|user userDefined key=value pair(s) ]\n"
    . " [ w|writeInput ]\n"
    . "\n"
    . " [ v|verbose verbosityLevel ]\n"
    . " [ options ]\n\n"
    # . "[Usage]: Will read data files from STDIN if not specified.\n"
    . "[Usage]: Missing arguments! Exiting...\n\n";

if (! defined($opts{'outDir'})  ||
    ! defined($opts{'fasFile'}) ||
    ! defined($opts{'srcName'})
    ) {
    print STDERR $usage;
    print STDERR opt2str(\%opts);
    
    exit(1);
}

if ($verbose) { print STDERR "Options::\n".opt2str(\%opts)."\n"; }

$opts{'outDir'}=~ s/\/+$//;
make_path($opts{'outDir'}) if (! -e $opts{'outDir'});
$opts{'outDir'}=abs_path($opts{'outDir'});

# buildDirs();
# loadDesigns(\@{$opts{'desFile'}});
# loadBins(\@{$opts{'binFile'}});
# loadCans(\@{$opts{'canFile'}});

loadFastas(\@{$opts{'fasFile'}}, \@{$opts{'srcName'}});

##
## Done!
##
@date = getTime(); runTimes();
mssg("[Done]: $date[0], $date[2], runTime=$runTimes[0]\n\n");
exit(0);

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Manifest Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Standard File Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

#
# Load Single Chromosome Sequence from Fasta File
#  chrs{chr}=seq
sub fas2chrSeq {
    my $func=(caller(0))[3];
    my ($datRef,$inpChr,$fasFile)=@_;

    die("[$func]: ERROR Missing paramaters (datRef, chr*, fasFile)\n") if (! defined($inpChr));
    die("[$func]: ERROR Missing paramaters (datRef, chr, fasFile*)\n") if (! defined($fasFile));

    die("[$func]: ERROR File='$fasFile' does not exist!\n") if (! -e $fasFile);
    $fasFile=abs_path($fasFile);

    my ($fasFH, %fasCols);
    my ($fasFileName, $fasFilePath, $fasFileSuffix) = fileparse( abs_path($fasFile), qr/\.[^.]*/);
    trimSuffix(\$fasFileName);
    mssg("[$func]: Loading file=$fasFile...\n",2);

    # Load fai index::
    #
    fai2idx($fasFile,$opts{'source'});

    my $cmd="samtools faidx $fasFile $inpChr";
    open($fasFH, "$cmd | ") or die("[$func]: Failed to open fasta stream=$cmd: $!");

    local $/ = ">";
    my $first = <$fasFH>;
    my $matChr;
    while(my $record = <$fasFH>){
        chomp $record;
        my $newline_loc = index($record,"\n");
        $matChr = substr($record,0,$newline_loc);
	warn("[$func]: Fasta Chromsome header=$matChr does not match input chr=$inpChr")
	    if ($inpChr ne $matChr);
        $$datRef{$inpChr}= uc(substr($record,$newline_loc+1));
        $$datRef{$inpChr}=~ tr/\n//d;
        $$datRef{$inpChr}=~ s/\s+.*$//;
        $$datRef{$inpChr}=substr($$datRef{$inpChr},0,$opts{'seqMax'}) if (defined($opts{'seqMax'}));
    } local $/ = "\n";
    close($fasFH);
    mssg("[$func]: Chromosome=$inpChr sequence=".length($$datRef{$inpChr})."\n",2,1);
    mssg("\n",2);
}

sub seq2iupacIdx {
    my $func=(caller(0))[3];
    my ($seq)=@_;

    die("[$func]: ERROR Function call with missing paramaters (seq)\n") if (! defined($seq)||length($seq)==0);

    my @cigs=();
    my $iupacMat="[RYSWKMBDHV]";
    my $posIdx;
    my $preIdx=0;
    while($seq =~ m/($iupacMat)/gi) {
        $posIdx = pos($seq);
	push @cigs,[$posIdx,$1];
    }
    return(@cigs);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Basic Subroutines:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

sub opt2str {
    my $func=(caller(0))[3];
    my ($opt_ref)=@_;

    my $ret_str = "";
    foreach my $key (sort keys %{$opt_ref}) {
	if (defined($$opt_ref{$key})) {
	    if (ref($$opt_ref{$key}) eq "ARRAY") {
		$ret_str .= "[Usage]: $key(s) => \n\t".join("\n\t",@{$$opt_ref{$key}})."\n";
	    } else {
		$ret_str .= "[Usage]: $key => '$$opt_ref{$key}'\n";
	    }
	} else {
	    $ret_str .= "[Usage]: $key => ''\n";
	}
    }

    return($ret_str);
}

sub setParams {
    my $func=(caller(0))[3];
    $lev5 = 5;

    $verbose   = TRUE;
    $topDir    = abs_path(getcwd);
    $statusCnt = 100000;
    $qsub = "qsub -cwd -pe threaded 16 -l excl=true -N";
    $qsubCnt=0;

    $forcastLim  = 16384;
    $forcastFlag = TRUE;
    $forcastCnt  = 0;
    $cpgPadLen   = 8;
    $chunkSize   = 100;

    $True='True'; $False='False'; $Ref='Ref'; $Alt='Alt'; $Red='Red'; $Grn='Grn'; $Mix='Mix'; $Nan='Nan';

    $A='A'; $B='B'; $C='C'; $D='D'; $E='E'; $F='F'; $G='G'; $H='H'; $I='I'; $J='J'; $K='K'; $L='L'; $M='M'; $N='N';
    $O='O'; $P='P'; $Q='Q'; $R='R'; $S='S'; $T='T'; $U='U'; $V='V'; $W='W'; $X='X'; $Y='Y'; $Z='Z';

    $d='d';

    $I='I'; $II='II'; $GRN='Grn'; $RED='Red'; $NEG='Neg';
    $L='_'; $SL='/';

    $CG='CG'; $TG='TG'; $UM='UM'; $MA='MA'; $FO='FO';
    @ACTG = ($A, $C, $T, $G); @FR  = ($F, $R); @CO  = ($C, $O); @MU  = ($M, $U); @MUD = ($M, $U, $D);
    @RG = ($RED, $GRN); @AB = ($A, $B); @III = ($I, $II); @TB=($T,$B);
    foreach my $nuc (@ACTG) { $ACTG{$nuc} = $nuc; }

    $TOP='TOP'; $BOT='BOT'; $UNK='UNK';
    $TAB="\t"; $COM=','; $COL=':'; $RET="\n"; $SEM=';'; 
    $PLUS='+'; $MINUS='-'; $DASH='-'; $DOT='.';

    $base2color{$A}=$Red; $base2color{lc($A)}=$Red;
    $base2color{$T}=$Red; $base2color{lc($T)}=$Red;

    $base2color{$C}=$Grn; $base2color{lc($C)}=$Grn;
    $base2color{$G}=$Grn; $base2color{lc($G)}=$Grn;

    $cmplChr{$F} = $R; $cmplChr{$R} = $F; ## F/R
    $cmplChr{$C} = $O; $cmplChr{$O} = $C; ## C/O
    $cmplChr{$M} = $U; $cmplChr{$U} = $M; ## M/U
    $cmplChr{$T} = $B; $cmplChr{$B} = $T; ## T/B
    $cmplChr{$D} = $D; ## Degnerate???

    $baseIdx{$T}=0; $idxBase[0]=$T;
    $baseIdx{$G}=1; $idxBase[1]=$G;
    $baseIdx{$A}=2; $idxBase[2]=$A;
    $baseIdx{$C}=3; $idxBase[3]=$C;
    
    $iF=0; $iR=1;
    $iT=0; $iB=1;
    $iC=0; $iO=1;
    $iU=0; $iM=1; $iD=2; $iN=3;

    $iTC=0; $iFC=0;
    $iTO=1; $iFO=1;
    $iBC=2; $iRC=2;
    $iBO=3; $iRO=3;

    $mu2idx{$U}=$iU; $idx2mu[$iU]=$U;
    $mu2idx{$M}=$iM; $idx2mu[$iM]=$M;
    $mu2idx{$D}=$iD; $idx2mu[$iD]=$D;
    $mu2idx{$N}=$iN; $idx2mu[$iN]=$N;

    $iSrd{$F}=$iF; $iSrd{$R}=$iR;
    $iSrd{$T}=$iT; $iSrd{$B}=$iB;
    $iSrd{$C}=$iC; $iSrd{$O}=$iO;
    $iSrd{$U}=$iU; $iSrd{$M}=$iM;

    $iSrd{$T.$C}=$iC+$iT;
    $iSrd{$T.$O}=$iO+$iT;
    $iSrd{$B.$C}=$iC+($iB*2);
    $iSrd{$B.$O}=$iO+($iB*2);

    $iSrd{'++'}=join($L,$F,$C);
    $iSrd{'+-'}=join($L,$F,$O);
    $iSrd{'-+'}=join($L,$R,$C);
    $iSrd{'--'}=join($L,$R,$O);

    $iSrd{$iFC}=join($L, $T, $C);
    $iSrd{$iFO}=join($L, $T, $O);
    $iSrd{$iBC}=join($L, $B, $C);
    $iSrd{$iBO}=join($L, $B, $O);
#    foreach my $i (sort keys %iSrd) { mssg("[$func]: iSrd[$i]=>$iSrd{$i}\n"); }

    $iSrd2bsp[$iFC]=$PLUS.$PLUS;
    $iSrd2bsp[$iFO]=$PLUS.$MINUS;
    $iSrd2bsp[$iRC]=$MINUS.$PLUS;
    $iSrd2bsp[$iRO]=$MINUS.$MINUS;

    $bspFR{$PLUS}  = $F; $bspFR{$F} = $PLUS;
    $bspFR{$MINUS} = $R; $bspFR{$R} = $MINUS;

    $bspTB{$PLUS}  = $T; $bspTB{$T} = $PLUS;
    $bspTB{$MINUS} = $B; $bspTB{$B} = $MINUS;

    $bspCO{$PLUS}  = $C; $bspCO{$C} = $PLUS;
    $bspCO{$MINUS} = $O; $bspCO{$O} = $MINUS;

    ## [CYSMBHV] [GRSKBDV]/[ATCYWMH]
    ## D -> tr/
    ##      U M D
    ## C -> T C Y/T
    ## Y -> T Y Y/T
    ## S -> K S B/K
    ## M -> W M H/W
    ## B -> K B B/K
    ## H -> W H H/W
    ## V -> D V N/D
    ## 
    %mapM = ();
    %mapD = ();
    my @dNucsM = qw(Y Y B H B H N);
    my @dNucsU = qw(T T K W K W D);
    my @cNucs  = qw(C Y S M B H V);
    my @gNucs  = qw(G R S K B D V);
    my @hNucs  = qw(A T C Y W M H N);
    my ($inp, $out);
    for (my $ii=0; $ii <@cNucs; $ii++) {
	my $cc = $cNucs[$ii];
	my $dd = $dNucsM[$ii];
	foreach my $gg (@gNucs) {
	    $inp = $cc.$gg;
	    $mapM{$inp} = lc($cc).$gg;
	    #mssg("[mapM]: $inp => ".lc($cc).$gg." [mapD]: $inp => ".lc($dd).$gg."\n");
	    $mapD{$inp} = lc($dd).$gg;
	}
    }

    my @letters = ($A, $B, $C, $D, $E, $F, $G, $H, $I, $J, $K, $L, $M, $N,
		   $O, $P, $Q, $R, $S, $T, $U, $V, $W, $X, $Y, $Z);

    foreach my $let (@letters) { $cmplNucs{$let} = cmpl(\$let); }

    push @{$iupacNucs{$R}}, $A; push @{$iupacNucs{$Y}}, $T;
    push @{$iupacNucs{$R}}, $G; push @{$iupacNucs{$Y}}, $C;

    # $bases2iupac{$A}=$A;
    # $bases2iupac{$C}=$C;
    # $bases2iupac{$T}=$T;
    # $bases2iupac{$G}=$G;

    $bases2iupac{$A.$G}=$R;     $bases2iupac{$C.$T}=$Y;
    $bases2iupac{$G.$A}=$R;     $bases2iupac{$T.$C}=$Y;

    push @{$iupacNucs{$S}}, $C; push @{$iupacNucs{$W}}, $A;
    push @{$iupacNucs{$S}}, $G; push @{$iupacNucs{$W}}, $T;

    $bases2iupac{$C.$G}=$S;     $bases2iupac{$A.$T}=$W;
    $bases2iupac{$G.$C}=$S;     $bases2iupac{$T.$A}=$W;

    push @{$iupacNucs{$K}}, $G; push @{$iupacNucs{$M}}, $A;
    push @{$iupacNucs{$K}}, $T; push @{$iupacNucs{$M}}, $C;

    $bases2iupac{$T.$G}=$K;     $bases2iupac{$C.$A}=$M;
    $bases2iupac{$G.$T}=$K;     $bases2iupac{$A.$C}=$M;

    push @{$iupacNucs{$B}}, $C; push @{$iupacNucs{$D}}, $A; push @{$iupacNucs{$H}}, $A; push @{$iupacNucs{$V}}, $A;
    push @{$iupacNucs{$B}}, $G; push @{$iupacNucs{$D}}, $G; push @{$iupacNucs{$H}}, $C; push @{$iupacNucs{$V}}, $C;
    push @{$iupacNucs{$B}}, $T; push @{$iupacNucs{$D}}, $T; push @{$iupacNucs{$H}}, $T; push @{$iupacNucs{$V}}, $G;

    push @{$iupacNucs{$N}}, $A;
    push @{$iupacNucs{$N}}, $C;
    push @{$iupacNucs{$N}}, $T;
    push @{$iupacNucs{$N}}, $G;

    %expandBase=();
    foreach my $refNuc (sort keys %iupacNucs) {
	foreach my $canNuc (sort @{$iupacNucs{$refNuc}}) {
	    push @{$expandBase{uc($refNuc)}}, uc($canNuc);
	    push @{$expandBase{lc($refNuc)}}, uc($canNuc);
	}
    }
    foreach my $refNuc (sort keys %ACTG) {
	push @{$expandBase{uc($refNuc)}}, uc($refNuc);
	push @{$expandBase{lc($refNuc)}}, uc($refNuc);
    }

    $bases2iupac{$C.$G.$T}=$B;    $bases2iupac{$A.$C.$G.$T}=$N;
    $bases2iupac{$C.$T.$G}=$B;    $bases2iupac{$A.$C.$T.$G}=$N;
    $bases2iupac{$G.$C.$T}=$B;    $bases2iupac{$A.$G.$C.$T}=$N;
    $bases2iupac{$G.$T.$C}=$B;    $bases2iupac{$A.$G.$T.$C}=$N;
    $bases2iupac{$T.$C.$G}=$B;    $bases2iupac{$A.$T.$C.$G}=$N;
    $bases2iupac{$T.$G.$C}=$B;    $bases2iupac{$A.$T.$G.$C}=$N;

    $bases2iupac{$A.$G.$T}=$D;    $bases2iupac{$C.$A.$G.$T}=$N;
    $bases2iupac{$A.$T.$G}=$D;    $bases2iupac{$C.$A.$T.$G}=$N;
    $bases2iupac{$G.$A.$T}=$D;    $bases2iupac{$C.$G.$A.$T}=$N;
    $bases2iupac{$G.$T.$A}=$D;    $bases2iupac{$C.$G.$T.$A}=$N;
    $bases2iupac{$T.$A.$G}=$D;    $bases2iupac{$C.$T.$A.$G}=$N;
    $bases2iupac{$T.$G.$A}=$D;    $bases2iupac{$C.$T.$G.$A}=$N;

    $bases2iupac{$A.$C.$T}=$H;    $bases2iupac{$G.$A.$C.$T}=$N;
    $bases2iupac{$A.$T.$C}=$H;    $bases2iupac{$G.$A.$T.$C}=$N;
    $bases2iupac{$C.$A.$T}=$H;    $bases2iupac{$G.$C.$A.$T}=$N;
    $bases2iupac{$C.$T.$A}=$H;    $bases2iupac{$G.$C.$T.$A}=$N;
    $bases2iupac{$T.$A.$C}=$H;    $bases2iupac{$G.$T.$A.$C}=$N;
    $bases2iupac{$T.$C.$A}=$H;    $bases2iupac{$G.$T.$C.$A}=$N;

    $bases2iupac{$A.$G.$C}=$V;    $bases2iupac{$T.$A.$G.$C}=$N;
    $bases2iupac{$A.$C.$G}=$V;    $bases2iupac{$T.$A.$C.$G}=$N;
    $bases2iupac{$G.$A.$C}=$V;    $bases2iupac{$T.$G.$A.$C}=$N;
    $bases2iupac{$G.$C.$A}=$V;    $bases2iupac{$T.$G.$C.$A}=$N;
    $bases2iupac{$C.$A.$G}=$V;    $bases2iupac{$T.$C.$A.$G}=$N;
    $bases2iupac{$C.$G.$A}=$V;    $bases2iupac{$T.$C.$G.$A}=$N;

    foreach my $deg (sort keys %iupacNucs) {
	my $nucs="";
	foreach my $nuc (sort @{$iupacNucs{$deg}}) {
	    $nucs.=$nuc;
	    $iupacTuple{lc($nuc)}{lc($deg)} = 1;
	    $iupacTuple{lc($deg)}{lc($nuc)} = 1;

	    $iupacTuple{lc($nuc)}{$deg} = 1;
	    $iupacTuple{$deg}{lc($nuc)} = 1;
	    $iupacTuple{$nuc}{lc($deg)} = 1;
	    $iupacTuple{lc($deg)}{$nuc} = 1;

	    $iupacTuple{$nuc}{$deg} = 1;
	    $iupacTuple{$deg}{$nuc} = 1;
	    # mssg("[$func]: iupacTuple{$deg}{$nuc}=$iupacTuple{$deg}{$nuc}\n");
	}
	$bases2iupac{$nucs}=$deg;
	$bases2iupac{$nucs}=$deg;
    }

    foreach my $str (sort keys %bases2iupac) {
	my $deg=$bases2iupac{$str};
	$iupac2bases{$deg}{$str}=@{$iupacNucs{$deg}};
    }
    if (FALSE) {
	foreach my $deg (sort keys %iupac2bases) {
	    foreach my $str (sort keys %{$iupac2bases{$deg}}) {
		mssg("$deg=$str\n");
	    }
	    mssg("\n");
	}
	mssg("\n");
    }

}

sub buildDir {
    my $func=(caller(0))[3];
    my ($dirName,$outDir)=@_;
    $outDir //=$opts{'outDir'};

    die("[$func]: ERROR Missing paramaters (dirName*)\n") if (! defined($dirName));
    die("[$func]: ERROR No output directory defined:$usage\n") if (! defined($outDir));

    $outDir=~ s/\/+$//;
    if (! -e $outDir) {
	make_path($outDir);
	$opts{'outDir'}=abs_path($opts{'outDir'}) if (! defined($opts{'outDir'}));
    }
    $dirName=~ s/\/+$//;
    $dirName=~ s/Dir$//;

    my $dirKey=$dirName."Dir";
    my $subDir="$outDir/$dirName";
    if (defined($opts{$dirKey}) && $opts{$dirKey} ne $subDir) {
	warn("[$func]: Multiple definitions in options{$dirKey}=$subDir=$opts{$dirKey}. Skipping...\n");
	return();
    }
    $opts{$dirKey}=$subDir;
    if (! -e $opts{$dirKey}) {
	make_path($opts{$dirKey});
	mssg("[$func]:\tBuilding dir($dirKey)=$opts{$dirKey}\n",3);
    }
    return($opts{$dirKey});
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Common File Format Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub fai2idx {
    my $func=(caller(0))[3];
    my ($inpFile,$srcName)=@_;

    die("[$func]: ERROR Function call with missing paramaters (head, [delimiter])\n")
	if (! defined($inpFile) || length($srcName)==0);

    die("[$func]: ERROR Fasta File does not exist: '$inpFile'!\n") if (! -e $inpFile);

    warn("[Warning]: Overwritting faiRef{$srcName}=$faiRef{$srcName} with $inpFile!\n")
	if (defined($faiRef{$srcName}));
    $faiRef{$srcName}=$inpFile;

#    $inpFile=~ s/\.gz/.fai/;
    $inpFile.=".fai";
    die("[$func]: ERROR Fasta Index File does not exist: '$inpFile'!\n") if (! -e $inpFile);

    my ($inpFH, $outFH, $binFH);
    my ($inpFileName, $inpFilePath, $inpFileSuffix, $isBin)=openFile(\$inpFH,$inpFile);
    mssg("[$func]: Loading file=$inpFile...\n", 3);

    my @lines=<$inpFH>;
    chomp(@lines);

    my $faiIdx=0;
    while(@lines) {
	my $line=shift(@lines);
        my @data=split(/\t/,$line);
	my $faiChr=shift(@data);

        if (! defined($fai2idx{$srcName}{$faiChr})) {
            $fai2idx{$srcName}{$faiChr}=$faiIdx;
            $idx2fai{$srcName}[$faiIdx]=$faiChr;
	    # Mark the standard 23 chromosomes
	    #
	    if ($faiChr =~ m/^[1-9]$/ || $faiChr =~ m/^[1-9][0-9]$/ || $faiChr =~ m/^[XYM]$/) {
		$isGenomic{$faiChr}=1;
		$isGenomic{$faiIdx}=1;		      
	    }
            $faiIdx++;
        }
    }
    close($inpFH);
    mssg("[$func]: Loaded faiIdx=$faiIdx, chromIdx=".scalar(keys %{$fai2idx{$srcName}}).
         ", idxChrom=".scalar(keys @{$idx2fai{$srcName}})."\n\n",2);
}

sub tabixStream {
    my $func=(caller(0))[3];
    my ($fh, $file, $range)=@_;
    $range //= "";

    die("[$func]: Function call with missing arguments (fh, file*, [range])\n") if (! defined($file));
    die("[$func]: file does not exist='$file'\n") if (! -e $file);
    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    trimSuffix(\$fileName);

    my $cmd="tabix -h $file $range";
    open($$fh, "$cmd | ") or die("[$func]: Failed to open stream='$cmd': $!");

    return($fileName, $filePath, $fileSuffix);
}

sub setHeader {
    my $func=(caller(0))[3];
    my ($datREF,$fh,$prefix,$del,$verbosity)=@_;
    $del //='\t';
    $verbosity //=0;

    die("[$func]: ERROR Function call with missing paramaters (colsREF,FH,prefix,[delimiter=tab],[verbosity])\n")
	if (! defined($prefix)||length($prefix)==0);

    my @headers=getHeader($fh,$prefix,$verbosity);
    my $header=$headers[$#headers];
    parseHeader(\%{$datREF},$header,$del,$verbosity);
    # foreach my $key (sort keys %{$datREF}) { mssg("[$func]: $key=>$$datREF{$key}\n") if ($verbosity > 0); }

    return(@headers);
}

sub getHeader {
    my $func=(caller(0))[3];
    my ($fh,$prefix,$verbosity)=@_;
    $verbosity //=0;

    my $header;
    my @lines=();
    while(<$fh>) {
	s/^\s+//;
	s/\s+$//;
	$header=$_;
	push @lines,$header;
	last if ($header =~m/^$prefix/);
    }
    die("[$func]: ERROR Failed to find header='^$prefix'\n") if (! defined($header));
    # mssg("[$func]: Prefix=$prefix\n".
    #     "[$func]: Header=$header\n") if ($verbosity > 0);
    return(@lines);
#    return($header);
}

sub parseHeader {
    my $func=(caller(0))[3];
    my ($datREF,$head,$del,$verbosity)=@_;
    $del //='\t';
    $verbosity //=0;

    die("[$func]: ERROR Function call with missing paramaters (colsREF,head,[delimiter=tab],[verbosity])\n")
	if (! defined($head) || length($head)==0);

    $head=~ s/\s+$//;
    my @dat=split(/$del/, $head);
    my $idx=0;
    # mssg("[$func]: HeaderLine=$head\n".
    # 	 "[$func]: Delimeter='$del'\n") if ($verbosity > 0);
    foreach my $key (@dat) { 
	# mssg("[$func]:\t$key=>$idx\n") if ($verbosity > 0);
	$$datREF{$key}=$idx;
	$idx++;
    }
}

sub parseHeaderIndex {
    my $func=(caller(0))[3];
    my ($idxREF,$colREF,$line,$del)=@_;
    $del //=$COM;

    die("[$func]: ERROR Missing arguments(idxREF,$colREF,line*, [del])\n") if (! defined($line));
    $line=~ s/\s+$//;
    my @dat=split(/$del/,$line);
    foreach my $i (0 .. $#dat) {
	my ($chip,$type)=split(/\./,$dat[$i]);
	$$idxREF[$i]=$dat[$i];
	$$colREF{$dat[$i]}=$i;
	#mssg("[$func]: $$idxREF[$i]=$dat[$i]\t$$colREF{$dat[$i]}=$i\n",3);
    }
}

sub parseHeaderGS {
    my $func=(caller(0))[3];
    my ($idxREF,$colREF,$line,$del)=@_;
    $del //=',';

    die("[$func]: ERROR Missing arguments(idxREF,$colREF,line*, [del])\n") if (! defined($line));
    $line=~ s/\s+$//;
    my @dat=split(/$del/, $line);
    foreach my $i (4 .. $#dat) {
	my ($chip,$type)=split(/\./,$dat[$i]);
	$$idxREF[$i]=[$chip,$type];
	$$colREF{$chip}{$type}=$i;
	#mssg("[$func]: $$idxREF[$i]=$dat[$i]\t$$colREF{$dat[$i]}=$i\n",3);
    }
}

sub parseHeaderVcf {
    my $func=(caller(0))[3];
    my ($vcfFH,$colsRef,$del,$com,$dbl)=@_;

    $del //='\t';
    $com //='#';
    $dbl //=$com.$com;

    my ($vcfSrc,$vcfBld,$refBld);
    die("[$func]: ERROR Function call with missing paramaters (vcfFH*,colRef)\n") if (! defined($vcfFH));

    my $line=<$vcfFH>;
    $line=~ s/\s+$//;
    while($line=~ m/^$dbl/) {
	$line=~ s/^$dbl//;
	$vcfSrc=$1 if ($line=~ m/^source=(.*)$/);
	$vcfBld=$1 if ($line=~ m/^dbSNP_BUILD_ID=([0-9]+)$/);
	$refBld=$1 if ($line=~ m/^reference=([^.]+)/);

	$line=<$vcfFH>;
	chomp($line);
    }
    die("[$func]: Failed to match VCF header=$line\n") if ($line !~ m/^$com/);

    my @dat=split(/$del/, $line);
    my $idx=0;
    foreach my $key (@dat) { $$colsRef{$key}=$idx; $idx++; }

    return($vcfSrc,$vcfBld,$refBld);
}

sub tenCnt {
    my $func=(caller(0))[3];
    my ($totCnt,$num)=@_;
    $num //=10;

    my $tenCnt=int($totCnt/10);
    $tenCnt++ if (! $tenCnt);
    return(($totCnt,$tenCnt,0));
}

sub loadList {
    my $func=(caller(0))[3];
    my ($inpFile, $datRef)=@_;

    die("[$func]: ERROR Function call with missing paramaters (listFile, listRef)\n") if (! defined($inpFile));
    $inpFile=abs_path($inpFile);

    ## Load Inprobe Design Data::
    ##
    my ($inpFH);
    my ($inpFileName, $inpFilePath, $inpFileSuffix, $isInpBin)=openFile(\$inpFH,$inpFile);
    mssg("[$func]: Loading file=$inpFile...\n", 2);

    my @inpFiles=<$inpFH>;
    close($inpFH);
    chomp(@inpFiles);
    my $fileCnt=scalar(@inpFiles);
    runTimes();

    while(@inpFiles) {
	my $file=shift(@inpFiles);
	die("[$func]: Listed file does not exists='$file'\n") if (! -e $file);
	next if ($file !~ m/\.22\./);

	push @{$datRef}, $file;
    }
    mssg("[$func]: Added $fileCnt to list size=".scalar(@{$datRef})."\n\n",3);
}

sub loadSingleFasta {
    my $func=(caller(0))[3];
    my $file = shift;
    my ($fasFH);
    die("[$func]: ERROR fasta file is not defined!\n") if (! defined($file) || ! -e $file);
    if ($file =~ m/\.gz/) { open($fasFH, "<:gzip", $file) or die("Failed to open $file stream: $!") }
    else { open($fasFH, "<$file") or die("Failed to open $file: $!") }
    mssg("[$func]: Loading single fasta=$file...\n", 3);

    local $/ = ">";
    my $first = <$fasFH>;
    while(my $record = <$fasFH>){
	chomp $record;
	my $newline_loc = index($record,"\n");
	my $chr = substr($record,0,$newline_loc);
	my $seq = substr($record,$newline_loc+1);
	$seq =~ tr/\n//d;
	$chr =~ s/\s+.*$//;
	close($fasFH);
	local $/ = "\n";
	return($chr, uc($seq));
    }
    die("[$func]: ERROR loadSingleFasta: Failed to find chromsome!\n");
}

sub file2sortedStream {
    die("[ERROR]: Function under construction!!!\n");
    my $func=(caller(0))[3];
    my ($fh,$idx,$file)=@_;
    
    die("[$func]: ERROR File does not exist '$file'\n")
	if (! defined($file) || ! -e $file);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    open($fh,"gzip -dc $file | sort -k $idx,$idx | ") or die("[$func]: Failed to open sorted stream from $file: $!")
	if ($fileSuffix eq '.gz');

    exit(1);
}

sub openFile {
    my $func=(caller(0))[3];
    my $fh=shift;
    my $file=shift;
    
    die("[$func]: ERROR File does not exist '$file'\n")
	if (! defined($file) || ! -e $file);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $isBin = ($fileName =~ /\.bin$/) ? 1 : 0;
    trimSuffix(\$fileName);

    $$fh = gzopen($file, "rb") or die("Failed to open binary gzip stream=$file for reading: $!")
	if ($isBin && $fileSuffix eq ".gz");

    open($$fh, "<:gzip", $file) or die("Failed to open gzip stream=$file for reading: $!")
	if (! $isBin && $fileSuffix eq ".gz");

    open($$fh, "<:raw", $file) or die("Failed to open binary file=$file for reading: $!")
	if ($isBin && $fileSuffix ne ".gz");

    open($$fh, "<$file") or die("Failed to open regular file=$file for reading: $!")
	if (! $isBin && $fileSuffix ne ".gz");

    return(($fileName, $filePath, $fileSuffix, $isBin));
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Execute External Program Functions:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

# Run bsmap on all improbe designs
#
sub bsmap {
    my $func=(caller(0))[3];
    my ($rSource,$rFasta,$qSource,$qFasta,$checkPrev)=@_;
    $checkPrev //=0;

    die("[$func]: ERROR Missing paramaters (refSrc*,refFas,querySrc,queryFas)\n") if (! defined($rSource));
    die("[$func]: ERROR Missing paramaters (refSrc,refFas*,querySrc,queryFas)\n") if (! defined($rFasta));
    die("[$func]: ERROR Missing paramaters (refSrc,refFas,querySrc*,queryFas)\n") if (! defined($qSource));
    die("[$func]: ERROR Missing paramaters (refSrc,refFas,querySrc,queryFas*)\n") if (! defined($qFasta));

    die("[$func]: ERROR File='$rFasta' does not exist!\n") if (! -e $rFasta);
    $rFasta=abs_path($rFasta);
    die("[$func]: ERROR File='$qFasta' does not exist!\n") if (! -e $qFasta);
    $qFasta=abs_path($qFasta);

    buildDir('bsp');
    my $bspFile=$opts{'bspDir'}."/".join($L,$rSource,$qSource).".bsp";
    my $sortedFile=$opts{'bspDir'}."/".join($L,$rSource,$qSource).".cgn-sorted.bsp.gz";
    if ($checkPrev) { mssg("[$func]:\tChecking($checkPrev) for previous alignment='$sortedFile'\n"); }
    else { mssg("[$func]:\tWill not check($checkPrev) for previous alignment='$sortedFile'\n"); }

    if ($checkPrev && -e $sortedFile) { mssg("[$func]:\tReturning previous results='$sortedFile'\n\n"); return($sortedFile); }

    my $cmd="$bsmapExe -s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R -u -a $qFasta -d $rFasta -o $bspFile";
    mssg("[$func]:\tcmd='$cmd'\n",3);
    my $sysCall=system($cmd);
    #$cmd="sort -k 4,4 $bspFile | gzip -c -> $sortedFile";
    $cmd="sort -k 1,1 $bspFile | gzip -c -> $sortedFile";
    mssg("[$func]: sorting='$cmd'\n");
    system($cmd);
    system("rm -f $bspFile");

    mssg("[$func]: bsmap alignment complete! Returning alignment sorted-results='$sortedFile'\n\n");
    return($sortedFile);
}

sub launchBsmap {
    my $func=(caller(0))[3];
    my ($rName, $rFile, $qName, $qFile, $fileIdx) = @_;

    die("[$func]: ERROR fasta queryFile=$qFile, name=$qName is null or does not exist!\n")
	if (! defined($qFile) || ! defined($qName) || ! -e $qFile);
    die("[$func]: ERROR fasta referenceFile=$rFile, name=$rName is null does not exist!\n")
	if (! defined($rFile) || ! defined($rName) || ! -e $rFile);
    $fileIdx  //= 1;

    my $binDir="$opts{'outDir'}/bin";
    my $bspDir="$opts{'outDir'}/bsp";
    make_path($bspDir) if (! -e $bspDir);
    
    my $binFH;
    my $outName=join($L, $rName, $qName);
    my $bspFile="$bspDir/$outName.idx$fileIdx.bsp";
    my $sortChrFile="$bspDir/$outName.idx$fileIdx.chr-sorted.bed";
    my $sortSeqFile="$bspDir/$outName.idx$fileIdx.seq-sorted.bed.gz";

    my $logFile="$bspDir/$outName.idx$fileIdx.log";
    my $binFile="$binDir/$outName.launch-bsmap.idx$fileIdx.sh";

    my $exe=$bsmapExe;

    system("rm -f bspFile.gz") if (-e "$bspFile.gz");
    my $cmd="$exe -a $qFile -d $rFile -o $bspFile ".
	"-s 12 -v 5 -g 0 -p 16 -n 1 -r 2 -R -u &> $logFile\n\n";

    $cmd .= "sort -k 4,4 $bspFile | gzip -c -> $sortSeqFile\n\n";

    open($binFH, ">$binFile") or die("Failed to open $binFile for writing: $!");
    print {$binFH} "$cmd";
    close($binFH);

    system("chmod 777 $binFile");
    $qsubCnt++;
    sleep(1);
    system("$qsub bsmap$qsubCnt $binFile") if (! $silent);
    sleep(1);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## File Conversion Executable Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub bwbbleIdx {
    my $func=(caller(0))[3];
    my ($file) = @_;

    if ($file =~ m/\.gz$/) { warn("[Warning]: Can't create bwwble index form gzip file=$file!\n"); return(); }
    my $cmd = "$bwbbleExe index $file";
    die("[$func]: ERROR Failed bwbble index: '$cmd'") if (system($cmd));
}

sub gzip {
    my $func=(caller(0))[3];
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);

    return($file) if ($fileSuffix eq ".gz");
    system("rm -f $file.gz") if (-e "$file.gz");
    die("[$func]: Failed bgzip $file") if (system("bgzip $file"));
    system("rm -f $file") if ($clean);
    return("$file.gz");
}

sub faiIdx {
    my $func=(caller(0))[3];
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $sortedBed = "$filePath$fileName.sorted.bed";

    if ($fileSuffix eq ".gz") { mssg("[Warning]: bgzipIdx file already compressed: $file!\n"); return(); }
    system("rm -f $file.gz") if (-e "$file.gz");
    system("rm -f $file.fai") if (-e "$file.fai");
    die("[$func]: Failed samtools faidx $file") if (system("samtools faidx $file"));
    die("[$func]: Failed gzip $file") if (system("gzip $file"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub bgzipIdx {
    my $func=(caller(0))[3];
    my ($file, $clean,$header)=@_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $sortedBed = "$filePath$fileName.sorted.bed";

    if ($fileSuffix eq ".gz") { mssg("[Warning]: bgzipIdx file already compressed: $file!\n"); return(); }
    system("rm -f $file.gz") if (-e "$file.gz");
    system("rm -f $file.gz.tbi") if (-e "$file.gz.tbi");
    die("[$func]: Failed bgzip $file") if (system("bgzip $file"));
    die("[$func]: Failed tabix $file.gz") if (system("tabix $file.gz"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub sortBgzipIdx {
    my $func=(caller(0))[3];
    my ($file,$clean,$header) = @_;
    $clean //=0; $header //="";

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    trimSuffix(\$fileName);
    my $sortedBed = "$filePath$fileName.chr-sorted.bed";

    $header="| grep -v $header" if (defined($header));
    my $cmd = "sort -k 1,1V -k 2,2n -k 3,3n $file $header > $sortedBed";
    if ($fileSuffix eq ".gz") { $cmd = "gzip -dc $file $header | sort -k 1,1V -k 2,2n -k 3,3n > $sortedBed"; }
    die("[$func]: Failed sort: '$cmd'") if (system("$cmd"));
    die("[$func]: Failed bgzip $sortedBed") if (system("bgzip -f -i $sortedBed"));
    $sortedBed .= ".gz";
    die("[$func]: Failed tabix $sortedBed") if (system("tabix -f $sortedBed"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub sortColZip {
    my $func=(caller(0))[3];
    my ($file,$colIdx,$name,$clean,$header) = @_;
    $clean //=0; $header //="";

    die("[$func]: ERROR Function call with missing paramaters (file, columnIndex)\n")
	if (! defined($file) || ! defined($colIdx));

    $name="col$colIdx" if (! defined($name));

    die("[$func]: ERROR Improbe file='$file' does not exist!\n") if (! -e $file);

    $header="| grep -v $header" if (defined($header));
    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    trimSuffix(\$fileName);
    my $sortedBed = "$filePath$fileName.$name-sorted.bed";

    my $cmd = "sort -k $colIdx,$colIdx $file $header > $sortedBed";
    if ($fileSuffix eq ".gz") { $cmd = "gzip -dc $file $header | sort -k $colIdx,$colIdx > $sortedBed"; }
    die("[$func]: Failed sort: '$cmd'") if (system("$cmd"));
    die("[$func]: Failed bgzip $sortedBed") if (system("bgzip -f $sortedBed"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub samToBam {
    my $func=(caller(0))[3];
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    trimSuffix(\$fileName);
    my $bam = "$filePath$fileName.bam";
    my $sortedBam = "$filePath$fileName.sorted.bam";

    my $cmd = "samtools view -S -b $file > $bam";
    die("[$func]: Failed sam <-> bam: $cmd") if (system($cmd));
    $cmd = "samtools sort $bam > $sortedBam";
    die("[$func]: Failed sort: $cmd") if (system($cmd));
    $cmd = "samtools index $sortedBam";
    die("[$func]: Failed index: $cmd") if (system($cmd));
    system("rm -f $file") if ($clean);
    system("rm -f $bam") if ($clean);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## BSC Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub cpgCnt {
    my $func=(caller(0))[3];
    my ($seq,$srdCO,$isDeg)=@_;
    $srdCO //=$C;
    $isDeg //=0;

    $seq=uc($seq);
    my $iupacCG="CG";
    if ($srdCO eq $C) { $iupacCG="R"; }
    elsif ($srdCO eq $O) { $iupacCG="Y"; }
    else { die("[$func]: ERROR Invalid strand=$srdCO (valid=$C|$O)\n\n"); }
    my $cpgCnt=0;
    while($seq =~ m/($iupacCG)/g) { $cpgCnt++; }
    return($cpgCnt);
}

sub bsc {
    my $func=(caller(0))[3];
    my ($seq, $type) = @_;

    die("\n[$func]: ERROR: bsc(seq, type): sequence not defined!\n") if (! defined($seq) || length($seq) == 0);
    die("\n[$func]: ERROR: bsc(seq, type): type not defined or lenggth 1!\n") if (! defined($type) || length($type) != 1);

    if ($type eq $D) {
        return(bscNewD($seq));
    } elsif ($type eq $M) {
	return(bscNewM($seq));
    } elsif ($type eq $U) {
        return(bscNewU($seq));
    } else {
        die("[ERROR]: Unsupported bsc type=$type!\n");
    }
}

sub bscBasic {
    my $func=(caller(0))[3];
    my ($seq,$type,$retUC)=@_;
    $type //=$U; $retUC //=0;

    die("[$func]: ERROR bsc missing parameters: (seq*,[type=U*MD],[returnUC])!\n") if (! defined($seq)||length($seq)==0);

    $type=$idx2mu[$type] if ($type =~ m/^[0-2]$/);
    
    $seq=uc($seq);
    if ($type eq $D) {
	#$seq=~ s/([CYSMBHV][GRSKBDV])/
	$seq=~ s/CG/yG/gi; }
    elsif ($type eq $M) {
	$seq=~ s/CG/cG/gi;
    } elsif ($type eq $U) { }
    else { die("[$func]: Unknown type=$type!\n"); }
    $seq=~ tr/CYSMBHV/ttkwkwd/;
    $seq=uc($seq) if ($retUC);
    return($seq);
}

sub bscU {
    my $func=(caller(0))[3];
    my $seq = shift;
    my $retUC = shift;
    die("[$func]: ERROR Seq=Null\n") if (! defined($$seq) || length($$seq) == 0);
    $$seq = uc($$seq);
    if (!$retUC) { $$seq =~ tr/CYSMBHV/ttkwkwd/; }
    if ( $retUC) { $$seq =~ tr/CYSMBHV/TTKWKWD/; }
}

sub bscNewU {
    my $func=(caller(0))[3];
    my $seq = shift;
    my $retUC = shift;
    die("[$func]: ERROR bscSeqU(): Seq=Null\n") if (! defined($seq) || length($seq) == 0);
    $seq = uc($seq);
    if (!$retUC) { $seq =~ tr/CYSMBHV/ttkwkwd/; }
    if ( $retUC) { $seq =~ tr/CYSMBHV/TTKWKWD/; }
    return($seq);
}

sub bscM {
    my $func=(caller(0))[3];
    my $seq = shift;
    my $retUC = shift;
    die("[$func]: ERROR bscSeqM(): Seq=Null\n") if (! defined($$seq) || length($$seq) == 0);
    $$seq = uc($$seq);
    $$seq =~ s/([CYSMBHV][GRSKBDV])/$mapM{$1}/g;
    $$seq =~ tr/CYSMBHV/ttkwkwd/;
    if ($retUC) { $$seq = uc($$seq); }
}

sub bscNewM {
    my $func=(caller(0))[3];
    my $seq = shift;
    my $retUC = shift;
    die("[$func]: ERROR bscSeqM(): Seq=Null\n") if (! defined($seq) || length($seq) == 0);
    $seq = uc($seq);
    $seq =~ s/([CYSMBHV][GRSKBDV])/$mapM{$1}/g;
    $seq =~ tr/CYSMBHV/ttkwkwd/;
    if ($retUC) { return(uc($seq)); }
    return($seq);
}

sub bscD {
    my $func=(caller(0))[3];
    my $seq = shift;
    my $retUC = shift;
    die("[$func]: ERROR bscSeqD(): Seq=Null\n") if (! defined($$seq) || length($$seq) == 0);
    $$seq = uc($$seq);
    $$seq =~ s/([CYSMBHV][GRSKBDV])/$mapD{$1}/g;
    $$seq =~ tr/CYSMBHV/ttkwkwd/;
    if ($retUC) { $$seq = uc($$seq); }
}

sub bscNewD {
    my $func=(caller(0))[3];
    my $seq = shift;
    my $retUC = shift;
    die("[$func]: ERROR bscSeqD(): Seq=Null\n") if (! defined($seq) || length($seq) == 0);
    $seq = uc($seq);
    $seq =~ s/([CYSMBHV][GRSKBDV])/$mapD{$1}/g;
    $seq =~ tr/CYSMBHV/ttkwkwd/;
    if ($retUC) { return(uc($seq)); }
    return($seq);
}

sub bscRcmD {
    my $func=(caller(0))[3];
    my $seq = shift;
    my $retUC = shift;
    die("[$func]: ERROR bscSeqD(): Seq=Null\n") if (! defined($seq) || length($seq) == 0);
    $seq = uc(reverse($seq));
    $seq =~ tr/ACTGRYSWKMBDHV/TGACYRSWMKVHDB/;
    $seq =~ s/([CYSMBHV][GRSKBDV])/$mapD{$1}/g;
    $seq =~ tr/CYSMBHV/ttkwkwd/;
    if ($retUC) { return(uc($seq)); }
    return($seq);
}

sub bscNew {
    my $func=(caller(0))[3];
    my ($type,$seq,$retUC)=@_;

    die("[$func]: ERROR bsc missing parameters: (type*,seq,[returnUC])!\n") if (! defined($type)||$type!~m/[UMD012]/);
    die("[$func]: ERROR bsc missing parameters: (type,seq*,[returnUC])!\n") if (! defined($seq)||length($seq)==0);

    $type=$idx2mu[$type] if ($type =~ m/^[0-2]$/);

    if ($type eq $D) { return(bscNewD($seq,$retUC)); }
    elsif ($type eq $M) { return(bscNewM($seq,$retUC)); }
    elsif ($type eq $U) { return(bscNewU($seq,$retUC)) }
    die("[$func]: ERROR Unsupported bsc type=$type!\n");
}

sub cmpl {
    my $func=(caller(0))[3];
    my $seq = shift;
    return("") if (! defined($seq)||length($seq)==0);
    $$seq =~ tr/ACTGRYSWKMBDHV/TGACYRSWMKVHDB/;
    $$seq =~ tr/actgryswkmbdhv/tgacyrswmkvhdb/;
}

sub revCmpl {
    my $seq = shift;
    return("") if (! defined($seq)||length($seq)==0);
    $seq = reverse($seq);
    cmpl(\$seq);
    return($seq);
}

sub cmplSeq {
    my $func=(caller(0))[3];
    my ($seq) = @_;

    die("[$func]: ERROR Cmpl: Seq=Null\n") if (! defined($seq));
    #return("") if (! defined($seq) || length($seq) == 0);

    $seq =~ tr/ACTGRYSWKMBDHV/TGACYRSWMKVHDB/;
    $seq =~ tr/actgryswkmbdhv/tgacyrswmkvhdb/;
    return($seq);
}

sub des2prbs {
    my $func=(caller(0))[3];
    my ($srd, $desSeq, $prbRef,$retUC) = @_;
    $retUC //=0;

#    my @desSetU=();
#    my @desSetM=();

    my $bscNewU = bscNewU($desSeq,$retUC);
    my $bscNewM = bscNewM($desSeq,$retUC);
    my $bscNewD = bscNewD($desSeq,$retUC);

    my @desSetU = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewU);
    my @desSetM = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewM);
    my @desSetD = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewD);

#    my @desSetU=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewU($desSeq,$retUC));
#    my @desSetM=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewM($desSeq,$retUC));

    my ($desNxbU, $desCpgU, $desSeqU, $desEndU);
    my ($desNxbM, $desCpgM, $desSeqM, $desEndM);
    my ($desNxbD, $desCpgD, $desSeqD, $desEndD);

#         (  $desNxbU,    $desCpgU,               $desSeqU,                $desEndU) = 
#   return( $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7],
#           $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7]) if ($desCO eq $C);
#          ( $desNxbM,    $desCpgM,               $desSeqM,                $desEndM) = 

#                  ( $desNxbU,            $desCpgU,                    $desSeqU,                  $desEndU) = 
#   return( revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]),
#           revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0])) if ($desCO eq $O);
#                  ( $desNxbM,            $desCpgM,                    $desSeqM,                  $desEndM) = 

    $$prbRef[$srd][$iU]=[ $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7] ];
    $$prbRef[$srd][$iM]=[ $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7] ];
    $$prbRef[$srd][$iD]=[ $desSetD[2], $desSetD[3], $desSetD[4].$desSetD[5].$desSetD[6], $desSetD[7] ];

    $srd++;
    $$prbRef[$srd][$iU]=[ revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]) ];
    $$prbRef[$srd][$iM]=[ revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0]) ];
    $$prbRef[$srd][$iD]=[ revCmpl($desSetD[4]), revCmpl($desSetD[3]), revCmpl($desSetD[1].$desSetD[2]), revCmpl($desSetD[0]) ];
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Top/Bot BSC Sequence Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub ref2bsc {
    my $func=(caller(0))[3];
    my ($bspSrd,$refSeq,$retUC,$subIdx,$subLen)=@_;

    die("[$func]: ERROR Function call with missing paramaters (srd*,refSeq,[idx,len])\n") if (! defined($bspSrd));
    die("[$func]: ERROR Function call with missing paramaters (srd,refSeq*,[idx,len])\n") if (! defined($refSeq));

    my @refSeqs=();
    $refSeqs[$iN]=$refSeq;
    if ($bspSrd eq "--") {
#	$refSeqs[$iN]=cmplSeq($refSeqs[$iN]);
	foreach my $i ($iU .. $iD) {
	    $refSeqs[$i]=bscNew($i,$refSeqs[$iN],$retUC);
#	    $refSeqs[$i]=cmplSeq($refSeqs[$i]);
	}
#	$refSeqs[$iN]=cmplSeq($refSeqs[$iN]);
    } elsif ($bspSrd eq "-+") {
	## FIXED!!!
	$refSeqs[$iN]=revCmpl($refSeqs[$iN]);
	foreach my $i ($iU .. $iD) { $refSeqs[$i]=bscNew($i,$refSeqs[$iN],$retUC); }
    } elsif ($bspSrd eq "+-") {
	## Looks FIXED???
	$refSeqs[$iN]=revCmpl($refSeqs[$iN]);
	foreach my $i ($iU .. $iD) { $refSeqs[$i]=bscNew($i,$refSeqs[$iN],$retUC); }
    } elsif ($bspSrd eq "++") {
	## FIXED!!!!
	foreach my $i ($iU .. $iD) { $refSeqs[$i]=bscNew($i,$refSeqs[$iN],$retUC); }
    } else {
	die("[$func]: ERROR Unrecognized strand='$bspSrd'\n");
    }
    if (defined($subIdx)) {
	if (defined($subLen)) {
	    foreach my $i ($iU .. $iN) { $refSeqs[$i]=substr($refSeqs[$i],$subIdx,$subLen); } }	    
	else {
	    foreach my $i ($iU .. $iN) { $refSeqs[$i]=substr($refSeqs[$i],$subIdx); } 
	}
    }
    return(@refSeqs);
}

sub srd2bed {
    my $func=(caller(0))[3];
    my ($srd,$pos,$len)=@_;
    $len //=50;

    die("[$func]: ERROR Function call with missing params (srd*,pos)\n") if (! defined($srd));
    die("[$func]: ERROR Function call with missing params (srd,pos*)\n") if (! defined($pos));

    my $start=$pos; my $end=$pos;
    if ($srd eq "++") { $start--; $end+=$len-1; }
    elsif ($srd eq "+-") { $start-=$len-1; }
    elsif ($srd eq "-+") { $start-=$len; $end++; }
    elsif ($srd eq "--") {$end+=$len+1-1; }
    else { die("[$func]: Unrecognized strand: $srd!\n"); }
    return(($start,$end));
}

sub returnTopBot_OLD {
    my $func=(caller(0))[3];
    my ($fwdSeq) = @_;

    my ($topSeq, $botSeq);
    my $revSeq = revCmp($fwdSeq);
    my $fwdTB  = substr(seqTopBotCG($fwdSeq), 0, 1);
    my $revTB  = substr(seqTopBotCG($revSeq), 0, 1);
    my $passTB = TRUE;
    $cgnQC{$fwdTB}++;
    $cgnQC{$revTB}++;

    my @ret = ($fwdSeq, $revSeq, $fwdTB, $revTB);
    if ($fwdTB eq $U && $revTB eq $U) {
        warn("[$func]: $U fwdTB=$fwdTB, revTB=$revTB\n".
             "[$func]: fwd=".addBrac($fwdSeq).", fwdTB=$fwdTB\n".
             "[$func]: rev=".addBrac($revSeq).", revTB=$revTB\n");
        return(@ret);
    } elsif ($fwdTB eq $S && $revTB eq $S) {
        warn("[$func]: $S fwdTB=$fwdTB, revTB=$revTB\n".
             "[$func]: fwd=".addBrac($fwdSeq).", fwdTB=$fwdTB\n".
             "[$func]: rev=".addBrac($revSeq).", revTB=$revTB\n");
        return(@ret);
    } elsif ($fwdTB eq uc($T) && $revTB eq uc($B)) {
        return(@ret);
    } elsif ($fwdTB eq uc($B) && $revTB eq uc($T)) {
        @ret = ($revSeq, $fwdSeq, $revTB, $fwdTB);
	return(@ret);
    } else {
        die("[$func]: ERROR: fwdTB($fwdTB) != revTB($revTB)\n".
            "[$func]: ERROR: fwd=".addBrac($fwdSeq).", fwdTB=$fwdTB\n".
            "[$func]: ERROR: rev=".addBrac($revSeq).", revTB=$revTB\n");
    }
}

sub returnTopBot {
    my $func=(caller(0))[3];
    my ($fwdSeq)=@_;

    die("[$func]: ERROR: Function call with missing paramaters (seq)\n")
	if (! defined($fwdSeq) || length($fwdSeq)==0);

    my $revSeq = revCmpl($fwdSeq);
    my $fwdTB  = substr(seq2topBotCG($fwdSeq), 0, 1);
    my $revTB  = substr(seq2topBotCG($revSeq), 0, 1);
    # mssg("[seqTB]:\tfwdTB=$fwdTB, revTB=$revTB\n\n");

    if ($fwdTB eq uc($T) && $revTB eq uc($B)) {
        return($fwdSeq, $revSeq, $fwdTB, $revTB, $F);
    } elsif ($fwdTB eq uc($B) && $revTB eq uc($T)) {
        return($revSeq, $fwdSeq, $revTB, $fwdTB, $R);
    } else {
        warn("[$func]: ERROR: NonTB: fwdTB=$fwdTB, revTB=$revTB\n".
             "[$func]: ERROR: NonTB: fwd=".addBrac($fwdSeq).", fwdTB=$fwdTB\n".
             "[$func]: ERROR: NonTB: rev=".addBrac($revSeq).", revTB=$revTB\n");
        return($fwdSeq, $revSeq, $fwdTB, $revTB, $N);
    }
}

sub seq2topBotCG {
    my $func=(caller(0))[3];
    my ($seq,$retName)=@_;

    die("[$func]: ERROR Function call with missing paramaters (seq)\n")
	if (! defined($seq) || length($seq)==0);

    $seq=uc($seq);
    shearBrac(\$seq);
    my $len = length($seq);
    my $mid = ($len/2) - 1;
    die("[$func]: ERROR Sequence is not even in length: mid=$mid; seq=$seq\n") if ($mid !~ m/^[0-9]+$/);

    my $top=$T; $top=$TOP if ($retName);
    my $bot=$B; $bot=$BOT if ($retName);
    my $unk=$U; $bot=$UNK if ($retName);

    my @nucs = split(//, uc($seq));
#    mssg("[Seq]:\t".addBrac($seq)."\n");
    foreach my $ii (1 .. $mid-1) {
	#for (my $ii = 1; $ii < $mid; $ii++) {
        my $nucA = $nucs[$mid-$ii];
        my $nucB = $nucs[$mid+$ii+1];
#	mssg("[$ii]:\t"." " x ($mid-$ii).$nucA." " x ($ii*2+2).$nucB."\n");

	next if ($nucA eq $nucB);
	next if ((defined($iupacTuple{$W}{$nucA}) && defined($iupacTuple{$W}{$nucB})) ||
		(defined($iupacTuple{$S}{$nucA}) && defined($iupacTuple{$S}{$nucB})));

        return($top) if (defined($iupacTuple{$W}{$nucA}) && defined($iupacTuple{$S}{$nucB}));
	return($bot) if (defined($iupacTuple{$S}{$nucA}) && defined($iupacTuple{$W}{$nucB}));
	return($U);
    }
    return($U);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Degenerate and CIGAR Sequence Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub doubleArray {
    my $func=(caller(0))[3];
    my ($datRef, $x) = @_;

    die("[$func]: ERROR Function call with missing paramaters (datRef*, x)\n") if (! defined($datRef));
    die("[$func]: ERROR Function call with missing paramaters (datRef, x*)\n") if (! defined($x));
    die("[$func]: ERROR Function call with non-numeric value x=$x\n") if ($x !~ m/^[0-9]+$/);

    ## FIXME::Not sure why we can just add onto the current array???
    my @ret = ();
    foreach my $seq (@{$datRef}) {
        for (my $ii = 0; $ii < $x; $ii++) {
            push @ret, $seq;
	}
    }
    return(@ret);
}

sub expandSeqs {
    my $func=(caller(0))[3];
    my ($seqREF)=@_;

    my @expStack=();
    foreach my $seq (@{$seqREF}) {
	my @expSeqs=expandSeq($seq);
	push @expStack, @expSeqs;
    }
    return(@expStack)
}

sub expandSeq {
    my $func=(caller(0))[3];
    my ($orgSeq,$expSEQS,$defVal)=@_;

    die("[$func]: ERROR Function call with missing paramaters (seq)\n") if (! defined($orgSeq)||length($orgSeq)==0);

    my @seqs   = ($orgSeq);
    my @seqDat = split(//, $orgSeq);

    my @degIdx = grep { $seqDat[$_] =~ /[^actgACTG]+/ } (0..$#seqDat);
    if (! @degIdx) {
	if (defined($defVal)) { foreach my $seq (@seqs) { $$expSEQS{uc($seq)}=$defVal; } }
	return(@seqs); }

    my $fc = 1;
    foreach my $dIdx (@degIdx) { $fc *= scalar(@{$iupacNucs{uc($seqDat[$dIdx])}}); }
    if ($fc > $forcastLim) { $forcastFlag = $fc; $forcastCnt++; die("[Forcast]: $fc>$forcastLim; $orgSeq!\n"); return(@seqs); }
#    if ($fc > $forcastLim) { $forcastFlag = $fc; $forcastCnt++; warn("[Forcast]: $fc>$forcastLim; $orgSeq!\n"); return(@seqs); }

    foreach my $ii (@degIdx) {
        my $seqNuc = uc($seqDat[$ii]);
        #mssg("[expandSeq]: seqNuc=$seqDat[$ii] <=> $seqNuc\n",$lev5) if (! defined($ACTG{$seqNuc}));                                                                            
        die("[ERROR]: Only checking IUPAC codes now!\n") if (! defined($iupacNucs{$seqNuc}));

        my $x = scalar(@{$iupacNucs{$seqNuc}});
        @seqs = doubleArray(\@seqs, $x);
        for (my $seqIdx = 0; $seqIdx < @seqs; $seqIdx++) {
            my $degNuc = lc($iupacNucs{$seqNuc}[$seqIdx % $x]);
            substr($seqs[$seqIdx], $ii, 1, $degNuc);
        }
    }
    if (defined($defVal)) { foreach my $seq (@seqs) { $$expSEQS{uc($seq)}=$defVal; } }
    return(@seqs);
}

sub cig2seq($$) {
    my $func=(caller(0))[3];
    my $cig=shift;
    my $seq=shift;

    die("[$func]: ERROR Function call with missing paramaters (cigar*, seq)\n") if (! defined($cig));
    die("[$func]: ERROR Function call with missing paramaters (cigar, seq*)\n") if (! defined($seq));

    $cig=uc($cig);
    my $idx=0;
    my $nuc;
    my $retSeq=$seq;
    while(length($cig)) {
	if ($cig =~ m/^([0-9]+)$/) {
	    return($retSeq);
	} elsif ($cig =~ m/^([0-9]+)([A-Z])/) {
	    $idx+=$1;
	    $nuc= $2;
	    substr($retSeq,$idx,1,$nuc);
	    $idx++;
	    $cig =~ s/^([0-9]+)([A-Z])//;
	} elsif ($cig=~ m/^([A-Z])/) {
	    $nuc=$1;
	    substr($retSeq,$idx,1,$nuc);
	    $idx++;
	    $cig=~ s/^([A-Z])//;
	} else {
	    die("[$func]: ERROR Unrecognized cigar=$cig, seq=$seq, retSeq=$retSeq\n");
	}
    }
    return($retSeq);
}

sub seqs2iupacCigar {
    my $func=(caller(0))[3];
    my ($seqA,$seqB)=@_;

    die("[$func]: ERROR Function call with missing paramaters (seqA*,seqB)\n")
	if (! defined($seqA) || length($seqA)==0);
    die("[$func]: ERROR Function call with missing paramaters (seqA,seqB*)\n")
	if (! defined($seqB) || length($seqB)==0);

    $seqA=uc($seqA); $seqB=uc($seqB);

    my $lenA=length($seqA);
    my $lenB=length($seqB);
    my $maxLen=max($lenA,$lenB);
    $seqA.=" " x ($maxLen-$lenA) if ($lenA!=$maxLen);
    $seqB.=" " x ($maxLen-$lenB) if ($lenB!=$maxLen);
    die("[$func]: Should never be here: $seqA/$seqB\n") if (length($seqA)!=length($seqB));

    my $bscMis=0; my $bscMat=0;
    my $iupMis=0; my $iupMat=0;
    my $misCnt=0; my $matCnt=0;
    my $snpCnt=0; my $witCnt=0;
    my $cigCnt=0;

    my $cigStr="";
    my $misStr=' ' x $maxLen;
#    return(($iupMis,$bscMis,$snpCnt,$witCnt,$misStr,$maxLen)) if ($seqA eq $seqB);
#    return(($iupMis,$bscMis,$snpCnt,($lenB-$lenA),$misStr,$maxLen)) if ($lenA<$lenB && $seqA =~m/^$seqB/);
#    return(($iupMis,$bscMis,$snpCnt,($lenA-$lenB),$misStr,$maxLen)) if ($lenB<$lenA && $seqB =~m/^$seqA/);

    my @charsA=split(//,$seqA);
    my @charsB=split(//,$seqB);


    foreach my $idx(0 .. $#charsA) {
	if ($charsA[$idx] eq ' ' || $charsB[$idx] eq ' ') { $witCnt++; substr($misStr, $idx, 1, 's'); next }
	if ($charsA[$idx] ne $charsB[$idx]) {
	    $snpCnt++;
	    if (defined($iupacTuple{$charsA[$idx]}{$charsB[$idx]})) {
		if ($cigCnt != 0) { $cigStr.=$cigCnt; $cigCnt=0; }
		die("[$func]: ERROR Undef chars(idx=$idx): charA='$charsA[$idx]', charB='$charsB[$idx]'\n".
		    "[$func]: ERROR  seqA=$seqA\n".
		    "[$func]: ERROR  seqB=$seqB\n") 
		    if (! defined($idx) || ! defined($cigStr) || ! defined($charsA[$idx]) || ! defined($charsB[$idx]));
		$cigStr.=lc(intersectIupac($charsA[$idx],$charsB[$idx]));
		$cigCnt++;

		$iupMat++; substr($misStr, $idx, 1, '*'); }
	    else {
		if ($cigCnt != 0) { $cigStr.=$cigCnt; $cigCnt=0; }
		$cigStr.=$charsB[$idx];
		$cigCnt++;

		$iupMis++;
		if ($charsA[$idx] eq $C && $charsB[$idx] eq $T) { $bscMat++; substr($misStr, $idx, 1, '^'); }
		elsif ($charsA[$idx] eq $A && $charsB[$idx] eq $G) { $bscMat++; substr($misStr, $idx, 1, '^'); }
		else { $bscMis++; substr($misStr, $idx, 1, '|'); }
	    }
	}
	else { $matCnt++; $cigCnt++; }
    }
    if ($cigCnt != 0) { $cigCnt.=$cigCnt; }
    return(($iupMis,$bscMis,$snpCnt,$cigStr,$misStr,$witCnt,$maxLen));
}

sub intersectIupac {
    my $func=(caller(0))[3];
    my ($degA,$degB)=@_;

    $degA=uc($degA);
    $degB=uc($degB);
    my @nucsA=($degA);
    my @nucsB=($degB);

    @nucsA=@{$iupacNucs{$degA}} if (defined($iupacNucs{$degA}));
    @nucsB=@{$iupacNucs{$degB}} if (defined($iupacNucs{$degB}));

    my %int=();
    foreach my $nucA (sort @nucsA) {
	foreach my $nucB (sort @nucsB) {
	    $int{$nucA}=1 if ($nucA eq $nucB);
	}
    }
    my $nucStr="";
    foreach my $nuc (sort keys %int) { $nucStr.=$nuc; }
#    mssg("[$func]: degA=$degA, degB=$degB, nucsA=@nucsA, nucsB=@nucsB, nucStr=$nucStr\n\n");
    return($nucStr) if (defined(lc($ACTG{$nucStr})));

    die("[$func]: ERROR Failed to find bases2iupac{$nucStr}!\n\n") if (! defined($bases2iupac{$nucStr}));
    my $intDeg=$bases2iupac{$nucStr};
    return($intDeg);
}

sub seqs2basicCigar {
    my $func=(caller(0))[3];
    my ($seqA,$seqB)=@_;

    die("[$func]: ERROR Missing agruments (seqA,seqB)\n")
	if (! defined($seqA) || length($seqA) == 0 ||
	    ! defined($seqB) || length($seqB) == 0);

    $seqA=uc($seqA);
    $seqB=uc($seqB);
    my $lenA=length($seqA);
    my $lenB=length($seqB);
    $seqB.=" " x ($lenA-$lenB) if ($lenA-$lenB>0);
    my $minLen=min($lenA,$lenB);

    return(($minLen,0)) if ($seqA eq $seqB || $seqA=~ m/^$seqB/);

    my @nucA=split(//,$seqA,$minLen);
    my @nucB=split(//,$seqB,$minLen);

    my $str=""; my $matCnt=0; my $misCnt=0;
    foreach my $idx (0 .. $#nucA) {
	if ($nucA[$idx] eq $nucB[$idx]) { $matCnt++; }
	else {
	    if ($matCnt != 0) { $str.=$matCnt; $matCnt=0; }
	    $str.=$nucB[$idx];
	    $misCnt++;
	}
    }
    if ($matCnt != 0) { $str.=$matCnt; }

#    mssg("[$func]: seqA=$seqA\tmisCnt\n".
#	 "[$func]:     =$str\n".
#	 "[$func]: seqB=$seqB\n\n") if ($str eq "0");
#    exit(1) if ($str eq "0");
    return(($str,$misCnt));
}

sub seqs2deg2 {
    my $func=(caller(0))[3];
    my ($prbA, $prbB) = @_;

    die("[$func]: ERROR prbA/B not defined!") if (! defined($prbA) || ! defined($prbB));
    die("[$func]: ERROR PrbA/B not equal length or zero: $prbA/$prbB!\n") if (length($prbA) != length($prbB) || length($prbA) == 0);

    my @mat = ();
    my $prbLen = length($prbA);
    my $sLen = 2;
    my @matA = split(//, uc($prbA));
    my @matB = split(//, uc($prbB));

    my $degPrb = "";
    for (my $ii = 0; $ii < $prbLen; $ii++) {
	die("[$func]: ERROR Unknown nuc=$matA[$ii], not found in basic or iupac!\n" .
            "[$func]: ERROR MatA=@matA\n") if (! defined($ACTG{$matA[$ii]}) && ! defined($iupacNucs{$matA[$ii]}));
        die("[$func]: ERROR Unknown nuc=$matB[$ii], not found in basic or iupac!\n" .
            "[$func]: ERROR MatB=@matB\n") if (! defined($ACTG{$matB[$ii]}) && ! defined($iupacNucs{$matB[$ii]}));

	# Sum individual nucleotides::A & B
	#
        my %nucs = ();
        if (defined($ACTG{$matA[$ii]})) { $nucs{$matA[$ii]} = 1; }
	else { map($nucs{$_}++, @{$iupacNucs{$matA[$ii]}}); }

        if (defined($ACTG{$matB[$ii]})) { $nucs{$matB[$ii]} = 1; }
	else { map($nucs{$_}++, @{$iupacNucs{$matB[$ii]}}); }

	# Build sorted ACCT String::
	#
        my $nucStr = "";
	foreach my $nuc (sort keys %nucs) { $nucStr .= $nuc; }

	# Add IUPAC Base::
	#
        if (defined($ACTG{$nucStr})) { $degPrb .= $ACTG{$nucStr}; next; }
        die("[$func]: ERROR Unable to find iupac for $nucStr!") if (! defined($bases2iupac{$nucStr}));
	$degPrb .= $bases2iupac{$nucStr};
    }
    die("[$func]: ERROR degPrb not equal to original length: $prbLen != ".length($degPrb)."\n") if (length($degPrb) ne $prbLen);
    return($degPrb);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Basic String Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub addBrac {
    my $func=(caller(0))[3];
    my ($seq) = @_;

    return($seq) if ($seq =~ m/\[.*\]/);

    my $idx = int(length($seq)/2);
    my $pre = substr($seq, 0, $idx-1);
    my $mid = substr($seq, $idx-1,2);
    my $suf = substr($seq, $idx+1);
    my $bra = $pre."[".$mid."]".$suf;
    die("[$func]: ERROR addBrac flanks are not 60!\n".
	"[$func]: ERROR seq=$seq\n".
	"[$func]: ERROR pre=$pre (len=".length($pre).")\n".
	"[$func]: ERROR mid=$mid (len=".length($mid).")\n".
	"[$func]: ERROR pre=$suf (len=".length($suf).")\n")
	if (length($pre) != 60 || length($suf) != 60);

    return($bra);
}

sub addBracs {
    die("FIXME\n");
    my ($seq, $swap) = @_;

    $swap = "CG" if (! defined($swap));
    $seq =~ s/($swap)/[$1]/gi;

    return($seq);
}

sub cleanTango {
    my $func=(caller(0))[3];
    my ($tango,$padLen,$trim)=@_;
    $padLen //= 8;

    return('') if (! defined($tango)||length($tango)==0);
    $tango =s/^$trim// if (defined($trim));
    my $tanLen=length($tango);
    return($tango) if ($tanLen==$padLen);

    $tango =~ s/^0+//;
    $tanLen=length($tango);
    die("[$func]: ERROR tango=$tango longer than $padLen; length=$tanLen\n") if ($tanLen>$padLen);

#    mssg("[$func]: padLen=$padLen, tango=$tango, tanLen=$tanLen\n");
    return(sprintf("%0".$padLen."d",$tango));
}

sub shearPadCG {
    my $seq = shift;
    $$seq =~ s/^[a-zA-Z]+//;
    $$seq =~ s/^0+//;
    return(0) if (length($$seq)==0);

    return($$seq);
}

sub shearBrac {
    my $seq = shift;
    $$seq =~ s/\[//; $$seq =~ s/\]//;
    return($$seq);
}

sub shearBracSeq {
    my ($seq) = @_;
    $seq =~ s/\[//; $seq =~ s/\]//;
    return($seq);
}

sub padCGN {
    my ($cpg,$pre) = @_;
    $pre //= "cg";

    my $len = length($cpg);
    $cpg =~ s/^cg//;
    $cpg = $pre. '0' x ($cpgPadLen - $len).$cpg;
    return($cpg);
}

sub iupacMisCnt {
    my $func=(caller(0))[3];
    my ($seqA,$seqB)=@_;

    die("[$func]: ERROR Function call with missing paramaters (seqA*,seqB)\n")
	if (! defined($seqA) || length($seqA)==0);
    die("[$func]: ERROR Function call with missing paramaters (seqA,seqB*)\n")
	if (! defined($seqB) || length($seqB)==0);

    $seqA=uc($seqA); $seqB=uc($seqB);

    my $lenA=length($seqA);
    my $lenB=length($seqB);
    my $maxLen=max($lenA,$lenB);
    $seqA.=" " x ($maxLen-$lenA) if ($lenA!=$maxLen);
    $seqB.=" " x ($maxLen-$lenB) if ($lenB!=$maxLen);
    die("[$func]: Should never be here: $seqA/$seqB\n") if (length($seqA)!=length($seqB));

    my $bscMis=0; my $bscMat=0;
    my $iupMis=0; my $iupMat=0;
    my $misCnt=0; my $matCnt=0;
    my $snpCnt=0; my $witCnt=0;

    my $misStr=' ' x $maxLen;
    return(($iupMis,$bscMis,$snpCnt,$witCnt,$misStr,$maxLen)) if ($seqA eq $seqB);
    return(($iupMis,$bscMis,$snpCnt,($lenB-$lenA),$misStr,$maxLen)) if ($lenA<$lenB && $seqA =~m/^$seqB/);
    return(($iupMis,$bscMis,$snpCnt,($lenA-$lenB),$misStr,$maxLen)) if ($lenB<$lenA && $seqB =~m/^$seqA/);

    my @charsA=split(//,$seqA);
    my @charsB=split(//,$seqB);

    foreach my $idx(0 .. $#charsA) {
	if ($charsA[$idx] eq ' ' || $charsB[$idx] eq ' ') { $witCnt++; substr($misStr, $idx, 1, 's'); next }
	if ($charsA[$idx] ne $charsB[$idx]) {
	    $snpCnt++;
	    if (defined($iupacTuple{$charsA[$idx]}{$charsB[$idx]})) { $iupMat++; substr($misStr, $idx, 1, '*'); }
	    else {
		$iupMis++;
		if ($charsA[$idx] eq $C && $charsB[$idx] eq $T) { $bscMat++; substr($misStr, $idx, 1, '^'); }
		elsif ($charsA[$idx] eq $A && $charsB[$idx] eq $G) { $bscMat++; substr($misStr, $idx, 1, '^'); }
		else { $bscMis++; substr($misStr, $idx, 1, '|'); }
	    }
	}
	else { $matCnt++; }
    }
    return(($iupMis,$bscMis,$snpCnt,$witCnt,$misStr,$maxLen));
}

sub diffStringsBSC {
    my $seqA=shift;
    my $seqB=shift;

    $seqA=uc($seqA);
    $seqB=uc($seqB);
    $seqA=~ s/\s+$//;
    $seqB=~ s/\s+$//;
    my $lenA=length($seqA);
    my $lenB=length($seqB);
    my $lenMax=max($lenA,$lenB);
    $seqA .= " " x ($lenMax-$lenA);
    $seqB .= " " x ($lenMax-$lenB);
    $lenA=length($seqA);
    $lenB=length($seqB);
#    mssg("seqA($lenA)='$seqA' max=$lenMax\n");
#    mssg("seqB($lenB)='$seqB' max=$lenMax\n");

    return(-1,-1,-1,' 'x $lenA) if ($lenA != $lenB);
    my $misStr=" " x $lenA;
    return(0,0,0,' 'x $lenA) if (uc($seqA) eq uc($seqB));

    my @charsA=split(//,$seqA);
    my @charsB=split(//,$seqB);
    my $misCnt=0; my $ctCnt=0; my $agCnt=0;
    foreach my $idx(0 .. $lenA-1) {
	#mssg("[$idx]: $misCnt\t$charsA[$idx]/$charsB[$idx]\n");
	if ($charsA[$idx] eq " " || $charsB[$idx] eq " ") {
	} elsif ($charsA[$idx] ne $charsB[$idx]) {
	    if ($charsA[$idx] eq $C && $charsB[$idx] eq $T) {
		$ctCnt++;
		substr($misStr, $idx, 1, '*');
	    } elsif ($charsA[$idx] eq $A && $charsB[$idx] eq $G) {
		$agCnt++;
		substr($misStr, $idx, 1, '*');
	    } else {
		$misCnt++;
		substr($misStr, $idx, 1, '|');
	    }
	}
    }
    return(($misCnt,$ctCnt,$agCnt,$misStr));
}

sub diffStrings {
    my $seqA=shift;
    my $seqB=shift;

    $seqA=uc($seqA);
    $seqB=uc($seqB);
    $seqA=~ s/\s+$//;
    $seqB=~ s/\s+$//;
    my $lenA=length($seqA);
    my $lenB=length($seqB);
    my $lenMax=max($lenA,$lenB);
    $seqA .= " " x ($lenMax-$lenA);
    $seqB .= " " x ($lenMax-$lenB);
    $lenA=length($seqA);
    $lenB=length($seqB);

    return(-1,-1,-1,' 'x $lenA) if ($lenA != $lenB);
    my $misStr=" " x $lenA;
    return(0,0,0,' 'x $lenA) if (uc($seqA) eq uc($seqB));

    my @charsA=split(//,$seqA);
    my @charsB=split(//,$seqB);
    my $misCnt=0;
    foreach my $idx(0 .. $lenA-1) {
	#mssg("[$idx]: $misCnt\t$charsA[$idx]/$charsB[$idx]\n");
	if ($charsA[$idx] eq " " || $charsB[$idx] eq " ") {
	} elsif ($charsA[$idx] ne $charsB[$idx]) {
	    $misCnt++;
	    substr($misStr, $idx, 1, '|');
	}
    }
    return(($misCnt,$misStr));
}

sub hashStr {
    my ($datRef, $cnt, $sep) = @_;

    $cnt=2 if (! defined($cnt));
    $sep=$RET if (! defined($sep));

    my $str = "";
    foreach my $key (sort keys %{$datRef}) {
	$str .=$TAB x $cnt ."$key => $$datRef{$key}$sep";
    }
    #$str =~ s/$sep$//;
    return($str);
}

sub hashCntStr {
    my $func=(caller(0))[3];
    my ($hash)=@_;

    my $str="";
    my %cnts=();
    foreach my $key (sort keys %{$hash}) {
	my $cnt = keys %{$$hash{$key}};
	$cnts{$key}=$cnt;
    }
    foreach my $key (sort {$$hash{$b} <=> $$hash{$a}} keys %cnts) {
	$str.="$key=$cnts{$key},";
    }
    $str=~ s/,$//;

    return($str);
}

sub sortedStr {
    my $func=(caller(0))[3];
    my ($hash)=@_;

    my $str="";
    foreach my $key (sort {$$hash{$b} <=> $$hash{$a}} keys %{$hash}) {
	$str.="$key=$$hash{$key},";
    }
    $str=~ s/,$//;

    return($str);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Cluster and File I/O Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub readStream {
    my ($fh) = @_;

    my @ret = ();
    while(<$fh>) {
	s/^\s+//;
	s/\s+$//;
	next if /^(\s)*$/;
	next if /^\#/;

	my $line = $_;
	push @ret, $line;
    }
    return(@ret);
}

sub listFiles {
    my ($dir, $key) = @_;

    my $cmd = "find $dir -name $key";
    mssg("[find]: $cmd\n");
    open(LIST, "$cmd | ") or die("Failed to open stream: $cmd for reading: $!");
    my @list = readStream(\*LIST);
    close(LIST);
    return(@list);
}

sub trimSuffix {
    my $str=shift;

    for (my $ii = 0; $ii < 5; $ii++) {
	$$str=~ s/\.gz$//;  $$str=~ s/\.fa$//;
	$$str=~ s/\.bed$//; $$str=~ s/\.txt$//; $$str=~ s/\.csv$//; $$str=~ s/\.tsv$//;
	$$str=~ s/\.bam$//; $$str=~ s/\.sam$//; $$str=~ s/\.vcf$//; $$str=~ s/\.fas$//;
	$$str=~ s/\.bsp$//; $$str=~ s/\.manifest$//; $$str=~ s/\.manifests$//;
	$$str=~ s/\.sorted$//; $$str =~ s/\.merged$//; $$str =~ s/\.intersect$//;
	$$str=~ s/\.[^\.]-sorted$//;
    }
}

sub qsubJob {
    my $func="qsubJob";
    my ($jobCmd,$jobDir,$jobName,$jobSub,$launch)=@_;

    die("[$func]: ERROR Function call missing arguments(cmd*, dir, name, subName, execute)\n") if (! defined($jobCmd));
    die("[$func]: ERROR Function call missing arguments(cmd, dir*, name, subName, execute)\n") if (! defined($jobDir));
    die("[$func]: ERROR Function call missing arguments(cmd, dir, name*, subName, execute)\n") if (! defined($jobName));
    die("[$func]: ERROR Function call missing arguments(cmd, dir, name, subName*, execute)\n") if (! defined($jobSub));
    die("[$func]: ERROR Function call missing arguments(cmd, dir, name, subName, execute*)\n") if (! defined($launch));

    make_path($jobDir) if (! -e $jobDir);

    my $shellFile = "$jobDir/$jobName.sh";
    open(SHELL, ">$shellFile") or die("Failed to open $shellFile for writing: $!");
    print SHELL "$jobCmd\n";
    close(SHELL);
    system("chmod 777 $shellFile");
    mssg("[qsub]: $shellFile\n",3);
    return(system("$qsub $jobSub $shellFile")) if ($launch);
    return(system($shellFile)) if ($local);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Time Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub runTimes
{
    my ($curTime, $perTime, $segTime, $totTime);
    $curTime = time;
    $segTime = sprintf("%.2f", ($curTime - $endTime)/60);
    $endTime = $curTime;
    $totTime = $endTime - $startTime;
    $totTime = sprintf("%.2f", ($totTime/60));
    unshift @runTimes, "$segTime(m)\t$totTime(m)";
}

sub getTime {
    my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
    my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    $yday += 1875;

    my $dateStr   = "$mon/$mday/$yday";
    my $dateStamp = "$sec.$min.$hour.$mday.$mon.$yday";

    my $tStr = "$mday $months[$mon] $days[$wday] $yday";
    return(($tStr, $dateStr, $dateStamp));
}


sub mssg {
    my $str=shift;
    my $level=shift;
    my $pTime=shift;

    if ($pTime) {
	runTimes();
	if ($str=~ m/\n$/) { $str=~ s/\n+$//; $str.=" time=$runTimes[0]\n"; }
	else { $str.=" time=$runTimes[0]"; }
    }
    if (defined($level)) {
	print STDERR $str if ($verbose >= $level);
    } else {
	print STDERR $str if ($verbose);
    }
}

sub status {
    my $func=(caller(0))[3];
    my ($name,$cur,$tot,$ver);
    $name=shift;
    $cur=shift;
    $tot=shift;
    $ver=shift;
    my @args=@_;

    return() if (! defined($verbose) || $ver > $verbose);

    my $tab="";
    $tab=$TAB x ($ver-2) if ($ver-2>0);
    my $str="[$name]:$tab".sprintf("%2d",100*$cur/$tot)."\% ($cur/$tot). ";
    $str.=join(", ",@args) if (@args);
    $str.="." if (@args);
    runTimes();
    $str.=" time=$runTimes[0]\n";
    print STDERR $str;
}

sub statsRecursion {
    my $func=(caller(0))[3];
    my ($cnt,$statsREF,$preStr)=@_;
    $cnt //=0;

    die("\n[$func]: ERROR Parameters(cnt,statsREF,key*)\n\n") if (! defined($preStr)||length($preStr)==0);

    my $tabA=$TAB x $cnt;
    my $tabB=$TAB x ($cnt+1);
    print STDERR "$tabA $preStr\n";
    foreach my $key (sort keys %{$statsREF}) {
	my $refType=ref($$statsREF{$key});
	unless ($refType) {
	    print STDERR $tabB." $key => $$statsREF{$key}\n"; }
	elsif ($refType eq 'HASH') {
	    statsRecursion($cnt+1,\%{$$statsREF{$key}},$key); }
	else {
	    die("\n[$func]: ERROR Unsupported stats reference type=$refType!\n\n"); }
    }
}

sub statusStats {
    my $func=(caller(0))[3];
    my ($ver)=@_;
    $ver //=$verbose;

    die("\n[$func]: ERROR Parameters(version=$ver) not numeric!\n\n") unless ($ver =~ m/^[0-9]+$/);
    return() if (! defined($verbose) || $ver > $verbose);

    my $cnt=$ver-1;
    $cnt=0 if ($cnt<=0);
    my @keys=sort keys %stats;
    foreach my $key (sort keys %stats) {
	my $refType=ref($stats{$key});
	unless ($refType) {
	    print STDERR $TAB x $cnt." $key => $stats{$key}\n"; }
	elsif ($refType eq 'HASH') {
	    statsRecursion($cnt+1,\%{$stats{$key}},$key); }
	else {
	    die("\n[$func]: ERROR Unsupported stats reference type=$refType!\n\n"); }
	print STDERR "\n";
    }
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Basic Math Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub isValidRange {
    my ($start, $end) = @_;

    return(1) if ($start =~ m/^[0-9]+$/ &&
		  $end =~ m/^[0-9]+$/ &&
		  $start > 0 &&
		  $start < $end);
    return(0);
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
#    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}
sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}

sub minIdx {
    my ($datRef)=@_;
    my $min=shift(@{$datRef});

    my $idx=0;
    return(($idx,$min)) if ($min !~ m/^[-0-9.]+$/ || $min eq '.');

    while(@{$datRef}) {
	my $val=shift(@{$datRef});
	$idx++;
	if ($val !~ m/^[-0-9.]+$/ || $val eq '.') { ($idx,$min)=($idx,$val); last; }
#	mssg("\t$idx,$min,$val\n");
	($idx,$min)=($idx,$val) if ($val<$min);
    }
    return(($idx,$min));
}

sub min {
    my ($a, $b) = @_;
    return($a) if (! defined($b) || $a < $b);
    return($b);
}

sub max {
    my ($a, $b) = @_;
    return($a) if (! defined($b) || $a > $b);
    return($b);
}

sub getMax {
    my ($ref) = @_;

    my @dat = sort {$b <=> $a} @{$ref};
    return($dat[0])
}

sub getMin {
    my ($ref) = @_;

    my @dat = sort {$a <=> $b} @{$ref};
    return($dat[0])
}

sub getSum {
    my ($ref) = @_;

    my $sum = 0;
    foreach my $val (@{$ref}) {
        $sum += $val;
    }
    return($sum);
}

sub getAvg {
    my ($ref) = @_;

    my @dat = sort @{$ref};
    my $mid=int(scalar(@dat)/2);
    my $med=$dat[$mid];
    my $sum = 0;
    foreach my $val (@dat) {
        next if (! defined($val) || length($val) == 0 || $val =~ m/^(\s)*$/);
        $sum += $val;
    }
    my $avg = $sum / scalar(@dat);
    return($avg,$med);
}

sub getStd {
    my ($ref, $avg) = @_;

    my @dat = @{$ref};
    my $sum = 0;
    foreach my $val (@dat) {
        next if (! defined($val) || length($val) == 0 || $val =~ m/^(\s)*$/);
        $sum += ($val - $avg) * ($val - $avg);
    }
    my $std = $sum / (scalar(@dat) - 1);
    $std = sqrt($std);
    return($std);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Old/Slow BSC Functions::
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub parseProbes {
    my $seq = shift;
    my $kmer = shift;
    $kmer = 48 if (! defined($kmer));
    my $seqLen = (length($seq)-2)/2;
    my $endLen = $seqLen-$kmer;

    shearBrac(\$seq);
    my ($endOpp, $prbOpp, $midCon, $prbCon, $endCon) = unpack("A$endLen"."A$kmer"."A2"."A$kmer"."A$endLen",$seq);
    my $nxbCon = substr($prbOpp, length($prbOpp)-1);
    my $midOpp = substr($midCon, 0, 1);
    my $nxbOpp = substr($midCon, 1, 1);

    my @opp = (revCmpl($nxbOpp), revCmpl($midOpp), revCmpl($prbOpp), revCmpl($endOpp));
    my @con = ($nxbCon, $midCon, $prbCon, $endCon);

    return(\@con, \@opp);
}

sub prbSet {
    my ($set, $seq, $idx, $kmer) = @_;
    $idx //= 0;
    $kmer //= 48;
    $idx*=2;

    my $seqLen = (length($seq)-2)/2;
    my $endLen = $seqLen-$kmer;

    shearBrac(\$seq);
    my ($endOpp, $prbOpp, $midCon, $prbCon, $endCon) = unpack("A$endLen"."A$kmer"."A2"."A$kmer"."A$endLen",$seq);
    my $nxbCon = substr($prbOpp, length($prbOpp)-1);
    my $midOpp = substr($midCon, 0, 1);
    my $nxbOpp = substr($midCon, 1, 1);

    ## Converted goes first
    $$set[0+$idx] = [ $nxbCon, $midCon, $prbCon, $endCon ];
    $$set[1+$idx] = [ revCmpl($nxbOpp), revCmpl($midOpp), revCmpl($prbOpp), revCmpl($endOpp) ];
}

sub bscSeqD {
    my $func=(caller(0))[3];
    my ($snvSeq, $retUP, $addBrac) = @_;

    die("[$func]: ERROR bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) == 0);
    my $bscSeq = ""; $snvSeq = uc($snvSeq);

    for (my $ii = 0; $ii < length($snvSeq)-1; $ii++) {
	my $snvCurr = substr($snvSeq, $ii, 1);
	my $snvNext = substr($snvSeq, $ii+1, 1);

	if (uc($snvCurr) =~ m/[CY]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($Y);
	} elsif (uc($snvCurr) =~ m/[S]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($B);
	} elsif (uc($snvCurr) =~ m/[M]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($H);
	} elsif (uc($snvCurr) =~ m/[B]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($B);
	} elsif (uc($snvCurr) =~ m/[H]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($H);
	} elsif (uc($snvCurr) =~ m/[V]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($N);


	} elsif (uc($snvCurr) =~ m/[CY]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $T;
	} elsif (uc($snvCurr) =~ m/[S]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/[M]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/[B]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/[H]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/[V]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $D;
	} elsif (uc($snvCurr) =~ m/[CYSMBHV]/ && uc($snvNext) =~ m/[ATYWMH]/) {
	    $bscSeq .= $T;
	} else {
	    $bscSeq .= $snvCurr;
	}
    }
    my $snvLast = substr($snvSeq, length($snvSeq)-1, 1);
    if ($snvLast =~ m/[CY]/) {
	$bscSeq .= $T;
    } elsif ($snvLast =~ m/[S]/) {
	$bscSeq .= $K;
    } elsif ($snvLast =~ m/[M]/) {
	$bscSeq .= $W;
    } elsif ($snvLast =~ m/[B]/) {
	$bscSeq .= $K;
    } elsif ($snvLast =~ m/[H]/) {
	$bscSeq .= $W;
    } elsif ($snvLast =~ m/[V]/) {
	$bscSeq .= $D;
    } else {
	$bscSeq .= $snvLast;
    }
    $bscSeq = uc($bscSeq) if ($retUP);
    $bscSeq = addBrac($bscSeq) if ($addBrac);
    return($bscSeq);
}

sub bscSeqM {
    my $func=(caller(0))[3];
    my ($snvSeq, $retUP, $addBrac) = @_;

    die("[$func]: ERROR bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) == 0);
    my $bscSeq = ""; $snvSeq = uc($snvSeq);

    for (my $ii = 0; $ii < length($snvSeq)-1; $ii++) {
	my $snvCurr = substr($snvSeq, $ii, 1);
	my $snvNext = substr($snvSeq, $ii+1, 1);

	if (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/[G]/) {
	    $bscSeq .= lc($C);
	} elsif (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/[RSKBDV]/) {
	    $bscSeq .= lc($Y);
	} elsif (uc($snvCurr) =~ m/[Y]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($Y);
	} elsif (uc($snvCurr) =~ m/[S]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($S);
	} elsif (uc($snvCurr) =~ m/[M]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($M);
	} elsif (uc($snvCurr) =~ m/[B]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($B);
	} elsif (uc($snvCurr) =~ m/[H]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($H);
	} elsif (uc($snvCurr) =~ m/[V]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($V);


	} elsif (uc($snvCurr) =~ m/[CY]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $T;
	} elsif (uc($snvCurr) =~ m/[S]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/[M]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/[B]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/[H]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/[V]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $D;

	} elsif (uc($snvCurr) =~ m/[CYSMBHV]/ && uc($snvNext) =~ m/[ATYWMH]/) {
	    $bscSeq .= $T;
	} else {
	    $bscSeq .= $snvCurr;
	}
    }
    my $snvLast = substr($snvSeq, length($snvSeq)-1, 1);
    if ($snvLast =~ m/[CY]/) {
	$bscSeq .= $T;
    } elsif ($snvLast =~ m/[S]/) {
	$bscSeq .= $K;
    } elsif ($snvLast =~ m/[M]/) {
	$bscSeq .= $W;
    } elsif ($snvLast =~ m/[B]/) {
	$bscSeq .= $K;
    } elsif ($snvLast =~ m/[H]/) {
	$bscSeq .= $W;
    } elsif ($snvLast =~ m/[V]/) {
	$bscSeq .= $D;
    } else {
	$bscSeq .= $snvLast;
    }
    $bscSeq = uc($bscSeq) if ($retUP);
    $bscSeq = addBrac($bscSeq) if ($addBrac);
    return($bscSeq);
}

sub bscSeqU {
    my $func=(caller(0))[3];
    my ($snvSeq, $retUP, $addBrac) = @_;

    die("[$func]: ERROR bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) == 0);
    my $bscSeq = ""; $snvSeq = uc($snvSeq);

    for (my $ii = 0; $ii < length($snvSeq)-1; $ii++) {
	my $snvCurr = substr($snvSeq, $ii, 1);
	my $snvNext = substr($snvSeq, $ii+1, 1);

	if (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/[G]/) {
	    $bscSeq .= lc($T);
	} elsif (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/[RSKBDV]/) {
	    $bscSeq .= lc($T);
	} elsif (uc($snvCurr) =~ m/[Y]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($T);
	} elsif (uc($snvCurr) =~ m/[S]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($K);
	} elsif (uc($snvCurr) =~ m/[M]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($W);
	} elsif (uc($snvCurr) =~ m/[B]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($K);
	} elsif (uc($snvCurr) =~ m/[H]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($W);
	} elsif (uc($snvCurr) =~ m/[V]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $bscSeq .= lc($D);


	} elsif (uc($snvCurr) =~ m/[CY]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $T;
	} elsif (uc($snvCurr) =~ m/[S]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/[M]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/[B]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/[H]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/[V]/ && uc($snvNext) =~ m/[ATCYWMH]/) {
	    $bscSeq .= $D;

	} elsif (uc($snvCurr) =~ m/[CYSMBHV]/ && uc($snvNext) =~ m/[ATYWMH]/) {
	    $bscSeq .= $T;
	} else {
	    $bscSeq .= $snvCurr;
	}
    }
    my $snvLast = substr($snvSeq, length($snvSeq)-1, 1);
    if ($snvLast =~ m/[CY]/) {
	$bscSeq .= $T;
    } elsif ($snvLast =~ m/[S]/) {
	$bscSeq .= $K;
    } elsif ($snvLast =~ m/[M]/) {
	$bscSeq .= $W;
    } elsif ($snvLast =~ m/[B]/) {
	$bscSeq .= $K;
    } elsif ($snvLast =~ m/[H]/) {
	$bscSeq .= $W;
    } elsif ($snvLast =~ m/[V]/) {
	$bscSeq .= $D;
    } else {
	$bscSeq .= $snvLast;
    }
    $bscSeq = uc($bscSeq) if ($retUP);
    $bscSeq = addBrac($bscSeq) if ($addBrac);
    return($bscSeq);
}

sub bscSeq {
    my $func=(caller(0))[3];
    my ($seq, $type, $retUP, $addBrac) = @_;

    die("[$func]: ERROR bsc(seq, type): sequence not defined!\n") if (! defined($seq) || length($seq) == 0);
    die("[$func]: ERROR bsc(seq, type): type not defined or lenggth 1!\n") if (! defined($type) || length($type) != 1);

    $retUP = FALSE if (! $retUP);
    $addBrac = FALSE if (! $addBrac);

    if ($type eq $D) {
	return(bscD($seq, $retUP, $addBrac));
    } elsif ($type eq $M) {
	return(bscM($seq, $retUP, $addBrac));
    } elsif ($type eq $U) {
	return(bscU($seq, $retUP, $addBrac));
    } else {
	die("[$func]: ERROR Unsupported bsc type=$type!\n");
    }

    return($seq);
}

## End of file
