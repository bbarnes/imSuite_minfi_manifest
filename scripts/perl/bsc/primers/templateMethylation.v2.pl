#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);
use File::Basename;

use constant FALSE => 0;
use constant TRUE  => 1;

## BED Format:
use constant BED_CHR    => 0;
use constant BED_START  => 1;
use constant BED_END    => 2;
use constant BED_NAME   => 3;
use constant BED_SCORE  => 4;
use constant BED_STRAND => 5;

## VCF Format:
use constant VCF_CHROM   => 0;
use constant VCF_POS     => 1;
use constant VCF_ID      => 2;
use constant VCF_REF     => 3;
use constant VCF_ALT     => 4;
use constant VCF_QUAL    => 5;
use constant VCF_FILTER  => 6;
use constant VCF_INFO    => 7;
use constant VCF_FORMAT  => 8;
use constant VCF_SAMPLE  => 9;
use constant VCF_GT      => 10;
use constant VCF_QC      => 11;
use constant VCF_REF_CNT => 12;
use constant VCF_ALT_CNT => 13;

## Order Format
use constant ORD_KEY     => 0;
use constant ORD_KEY_A   => 1;
use constant ORD_SEQ_A   => 2;
use constant ORD_KEY_B   => 3;
use constant ORD_SEQ_B   => 4;
use constant ORD_BIN     => 5;
use constant ORD_INFTYPE => 6;

## BAM Format:
use constant BAM_QNAME    => 0;
use constant BAM_FLAG     => 1;
use constant BAM_RNAME    => 2;
use constant BAM_POS      => 3;
use constant BAM_MAPQ     => 4;
use constant BAM_CIGAR    => 5;
use constant BAM_RNEXT    => 6;
use constant BAM_PNEXT    => 7;
use constant BAM_TLEN     => 8;
use constant BAM_SEQ      => 9;
use constant BAM_QUAL     => 10;
use constant BAM_SAMPLE   => 11;

## Methylation Manifest File
my @manHeader = ();
use constant MAN_ILMNID               => 0;     push @manHeader, "IlmnID";
use constant MAN_NAME                 => 1;     push @manHeader, "Name";
use constant MAN_ADDRESSA_ID          => 2;     push @manHeader, "AddressA_ID";
use constant MAN_ALLELEA_PROBESEQ     => 3;     push @manHeader, "AlleleA_ProbeSeq";
use constant MAN_ADDRESSB_ID          => 4;     push @manHeader, "AddressB_ID";
use constant MAN_ALLELEB_PROBESEQ     => 5;     push @manHeader, "AlleleB_ProbeSeq";
use constant MAN_INFINIUM_DESIGN_TYPE => 6;     push @manHeader, "Infinium_Design_Type";
use constant MAN_NEXT_BASE            => 7;     push @manHeader, "Next_Base";
use constant MAN_COLOR_CHANNEL        => 8;     push @manHeader, "Color_Channel";
use constant MAN_FORWARD_SEQUENCE     => 9;     push @manHeader, "Forward_Sequence";
use constant MAN_GENOME_BUILD         => 10;    push @manHeader, "Genome_Build";
use constant MAN_CHR                  => 11;    push @manHeader, "CHR";
use constant MAN_MAPINFO              => 12;    push @manHeader, "MAPINFO";
use constant MAN_SOURCESEQ            => 13;    push @manHeader, "SourceSeq";
use constant MAN_STRAND               => 14;    push @manHeader, "Strand";

## Probe Parameters:
use constant MAN_PROBE_CPG            => 15;    push @manHeader, "probeCG";
use constant MAN_STRAND_TB            => 16;    push @manHeader, "probeTB";
use constant MAN_STRAND_FR            => 17;    push @manHeader, "probeFR";
use constant MAN_STRAND_CO            => 18;    push @manHeader, "probeCO";

use constant MAN_SCR_M                => 19;    push @manHeader, "probeScoreM";
use constant MAN_SCR_U                => 20;    push @manHeader, "probeScoreU";
use constant MAN_SCR_MIN              => 21;    push @manHeader, "probeScoreMin";

use constant MAN_NUM_CPGS             => 22;    push @manHeader, "probeUnderlyingCpgsCount";
use constant MAN_NUM_SNVS             => 23;    push @manHeader, "probeUnderlyingSnvsCount";

use constant MAN_SNV_POS              => 24;    push @manHeader, "snvCoordinates";
use constant MAN_SNV_REF              => 25;    push @manHeader, "snvRefs";
use constant MAN_SNV_ALT              => 26;    push @manHeader, "snvAlts";
use constant MAN_VAR_POS              => 27;    push @manHeader, "iupacSnvCoordinates";
use constant MAN_VAR_REF              => 28;    push @manHeader, "iupacSnvRefs";
use constant MAN_VAR_ALT              => 29;    push @manHeader, "iupacSnvAlts";

## Order vs. Design Comparison
use constant MAN_OFF_TARGET           => 30;    push @manHeader, "offTargetMatch";
use constant MAN_MATCH_CALL           => 31;    push @manHeader, "designMatch";

use constant MAN_VAR_WASH_POS         => 32;    push @manHeader, "probeAbsentPosition";
use constant MAN_VAR_WASH_CNT         => 33;    push @manHeader, "probeAbsentCount";
use constant MAN_VAR_MISS_POS         => 34;    push @manHeader, "sampleAbsentPosition";
use constant MAN_VAR_MISS_CNT         => 35;    push @manHeader, "sampleAbsentCount";

## Ordered Probe
use constant MAN_ORD_KEY_S            => 36;    push @manHeader, "orderSnvKey";
use constant MAN_ORD_KEY_G            => 37;    push @manHeader, "orderGenKey";
use constant MAN_ORD_CIRGAR_S         => 38;    push @manHeader, "orderCigarSnv";
use constant MAN_ORD_CIRGAR_G         => 39;    push @manHeader, "orderCigarGeneral";
use constant MAN_ORD_PRB_U            => 40;    push @manHeader, "orderProbeU";
use constant MAN_ORD_PRB_M            => 41;    push @manHeader, "orderProberM";
use constant MAN_ORD_PRB_D            => 42;    push @manHeader, "orderProbeD";

use constant MAN_ORD_CPG              => 43;    push @manHeader, "orderCpgCount";
use constant MAN_ORD_SNV              => 44;    push @manHeader, "orderSnvCount";
use constant MAN_ORD_COST_U           => 45;    push @manHeader, "orderCostU";
use constant MAN_ORD_COST_M           => 46;    push @manHeader, "orderCostM";
use constant MAN_ORD_COST_D           => 47;    push @manHeader, "orderCostD";

## Design Probe
use constant MAN_DES_KEY_S            => 48;	push @manHeader, "designSnvKey";
use constant MAN_DES_KEY_G            => 49;	push @manHeader, "designGenKey";
use constant MAN_DES_CIRGAR_S         => 50;	push @manHeader, "designCigarSnv";
use constant MAN_DES_CIRGAR_G         => 51;	push @manHeader, "designCigarGeneral";
use constant MAN_DES_PRB_U            => 52;	push @manHeader, "designProbeU";
use constant MAN_DES_PRB_M            => 53;	push @manHeader, "designProbeM";
use constant MAN_DES_PRB_D            => 54;	push @manHeader, "designProbeD";

use constant MAN_DES_CPG              => 55;	push @manHeader, "designCpgCount";
use constant MAN_DES_SNV              => 56;	push @manHeader, "designSnvCount";
use constant MAN_DES_COST_U           => 57;	push @manHeader, "designCostU";
use constant MAN_DES_COST_M           => 58;	push @manHeader, "designCostM";
use constant MAN_DES_COST_D           => 59;	push @manHeader, "designCostD";
use constant MAN_SRC_TYPE             => 60;	push @manHeader, "orderSourceType";

use constant MAN_LAST_KEY             => 60;

## SNP Chip (Omni5/GSA) Format:
use constant MAN_SNP_ILMNID           => 0;
use constant MAN_SNP_NAME             => 1;
use constant MAN_SNP_ILMNSTRAND       => 2;
use constant MAN_SNP_SNP              => 3;
use constant MAN_SNP_ADDRESSA_ID      => 4;
use constant MAN_SNP_ALLELEA_PROBESEQ => 5;
use constant MAN_SNP_ADDRESSB_ID      => 6;
use constant MAN_SNP_ALLELEB_PROBESEQ => 7;
use constant MAN_SNP_GENOMEBUILD      => 8;
use constant MAN_SNP_CHR              => 9;
use constant MAN_SNP_MAPINFO          => 10;
use constant MAN_SNP_PLOIDY           => 11;
use constant MAN_SNP_SPECIES          => 12;
use constant MAN_SNP_SOURCE           => 13;
use constant MAN_SNP_SOURCEVERSION    => 14;
use constant MAN_SNP_SOURCESTRAND     => 15;
use constant MAN_SNP_SOURCESEQ        => 16;
use constant MAN_SNP_TOPGENOMICSEQ    => 17;
use constant MAN_SNP_BEADSETID        => 18;
use constant MAN_SNP_EXP_CLUSTERS     => 19;
use constant MAN_SNP_REFSTRAND        => 20;

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

## UCLA Design Format:
use constant UCLA_PROBESET_ID              => 0;
use constant UCLA_SEQ_ID                   => 1;
use constant UCLA_FORWARD_SEQUENCE         => 2;
use constant UCLA_GENOME_BUILD             => 3;
use constant UCLA_CHROMOSOME               => 4;
use constant UCLA_COORDINATE               => 5;
use constant UCLA_DESIGN_STATE             => 6;
use constant UCLA_SEQ_LENGTH               => 7;
use constant UCLA_FORWARD_CPG_COORD        => 8;
use constant UCLA_TB_STRAND                => 9;
use constant UCLA_TOP_SEQUENCE             => 10;
use constant UCLA_TOP_CPG_COORD            => 11;
use constant UCLA_PROBE_TYPE               => 12;
use constant UCLA_PROBESET_SCORE           => 13;
use constant UCLA_METHYL_PROBE_SEQUENCE    => 14;
use constant UCLA_ALLELE_FR_STRAND         => 15;
use constant UCLA_ALLELE_TB_STRAND         => 16;
use constant UCLA_ALLELE_CO_STRAND         => 17;
use constant UCLA_METHYL_PROBE_SCORE       => 18;
use constant UCLA_UNDERLYING_CPG_COUNT     => 19;
use constant UCLA_UNMETHYL_PROBE_SEQUENCE  => 20;
use constant UCLA_UNMETHYL_PROBE_SCORE     => 21;
use constant UCLA_INFINIUM_SCORESEQ_ID     => 22;
use constant UCLA_NUM_SPECIES              => 23;
use constant UCLA_SPECIES                  => 24;
use constant UCLA_PROBE_START_COORD        => 25;
use constant UCLA_PROBE_END_COORD          => 26;
use constant UCLA_REFERENCE_PROBE_SEQUENCE => 27;
use constant UCLA_SNV_LOCATION             => 28;
use constant UCLA_SNV_ORIGINAL             => 29;
use constant UCLA_SNV_CHANGE               => 30;
use constant UCLA_INFINIUM_TYPE            => 31;
use constant UCLA_IS_EPIC_SITE             => 32;
use constant UCLA_IS_EPIC_DESIGN           => 33;
use constant UCLA_NVARIANTS                => 34;
use constant UCLA_BABOON                   => 35;
use constant UCLA_CAT                      => 36;
use constant UCLA_CHIMP                    => 37;
use constant UCLA_COW                      => 38;
use constant UCLA_DOG                      => 39;
use constant UCLA_GIBBON                   => 40;
use constant UCLA_GREENMONKEY              => 41;
use constant UCLA_HORSE                    => 42;
use constant UCLA_HUMAN                    => 43;
use constant UCLA_MACACQUE                 => 44;
use constant UCLA_MARMOSET                 => 45;
use constant UCLA_MOUSE                    => 46;
use constant UCLA_RABBIT                   => 47;
use constant UCLA_RAT                      => 48;
use constant UCLA_RHESUSMONKEY             => 49;
use constant UCLA_SHEEP                    => 50;
use constant UCLA_ALL_MAPPABLE             => 51;
use constant UCLA_NONPRIMATE_MAPPABLE      => 52;
use constant UCLA_RRBS                     => 53;
use constant UCLA_ISBABOON                 => 54;
use constant UCLA_ISCAT                    => 55;
use constant UCLA_ISCHIMP                  => 56;
use constant UCLA_ISCOW                    => 57;
use constant UCLA_ISDOG                    => 58;
use constant UCLA_ISGIBBON                 => 59;
use constant UCLA_ISGREENMONKEY            => 60;
use constant UCLA_ISHORSE                  => 61;
use constant UCLA_ISHUMAN                  => 62;
use constant UCLA_ISMACACQUE               => 63;
use constant UCLA_ISMARMOSET               => 64;
use constant UCLA_ISMOUSE                  => 65;
use constant UCLA_ISRABBIT                 => 66;
use constant UCLA_ISRAT                    => 67;
use constant UCLA_ISRHESUSMONKEY           => 68;
use constant UCLA_ISSHEEP                  => 69;
use constant UCLA_NUMMAPS                  => 70;
use constant UCLA_NUMGENOMESTESTING        => 71;
use constant UCLA_MAPFILTER                => 72;

## Input format:
use constant IN_SEQ_ID       => 0;
use constant IN_SEQUENCE     => 1;
use constant IN_GENOME_BUILD => 2;
use constant IN_CHROMOSOME   => 3;
use constant IN_COORDINATE   => 4;
use constant IN_CPG_ISLAND   => 5;

##
## Common Methylation && Format Variables:
##
my ($II);
my @TB = (); my @FR = (); my @CO = (); my @III = (); my @MU = (); my @ACTG = ();
my ($SL, $TAB, $COM, $COL, $RET, $SEM, $PLUS, $MINUS, $DASH);
my ($A, $B, $C, $D, $E, $F, $G, $H, $I, $J, $K, $L, $M, $N,
    $O, $P, $Q, $R, $S, $T, $U, $V, $W, $X, $Y, $Z);
#my ($a, $b, $c, $d, $e, $f, $g, $h, $i, $j, $k, $l, $m, $n,
#    $o, $p, $q, $r, $s, $t, $u, $v, $w, $x, $y, $z);
my %iupacNucs = (); my %cmplNucs = ();
my %iupacStrs = (); my %cmplChar = ();
my %basicNucs = (); my %ACTG = ();

##
## Basic Program Variables
##
my ($topDir, $statusCnt, $inpMax, $desMax, $snvMax, $manMax, $qsub);
my $verbose   = TRUE;

my $lev5 = 5;
my $pLen = 60;
my $forcastLim  = 16384;
my $forcastFlag = TRUE;
my $forcastCnt  = 0;

##
## Input Variables and Structures
##
my ($outDir, $logDir, $bedDir, $datDir, $fasDir,
    $bamDir, $ordDir, $desDir, $manDir, $vizDir);

my @desFiles = ();  ## Probe Designs
my @snvFiles = ();

my @ordFiles = ();  ## Manufacturing
my @matFiles = ();
my @aqpFiles = ();
my @pqcFiles = ();

my @manFiles = ();  ## Manifests

##
## Output Variables
##
my %vizFhsI  = ();
my %vizFhsII = ();
my %vizFnsI  = ();
my %vizFnsII = ();

##
## Variables and Data Structures:
##
my %manDat = ();
my %snvDat = ();
my %desDat = ();
my %desMap = ();

my %ordDat = ();
my %ordMap = ();
my %dupDat = ();
my %dupMap = ();
my %ordKeys = ();

my %pqcMap = ();

my %matDat = ();
my %aqpDat = ();
my %pqcDat = ();

my $projName = "HorvathMammalMethyl";

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Program Main Subroutines:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub loadDesign($) {
    my ($file) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    $fileName = trimSuffix($fileName);

    if ($fileSuffix eq ".gz") { open(DES, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(DES, "<$file") or die("Failed to open $file: $!") }
    mssg("[Design]: Loading designFile=$file...\n", 3);

    my $inpCnt = 0;
    while(<DES>) {
        s/^\s+//;
        s/\s+$//;
	next if /^Seq_ID/;

        my $line = $_;
        my @data = split(/\t/, $line);

        my $prbCG  = $data[SEQ_ID];
        my $prbFR  = $data[METHYL_ALLELE_FR_STRAND];
        my $prbTB  = substr($data[METHYL_ALLELE_TB_STRAND], 0, 1);
        my $prbCO  = $data[METHYL_ALLELE_CO_STRAND];

	if (defined($desDat{$prbCG}{$prbTB}{$prbCO})) {
	    die("[ERROR]: Design already defined: {$prbCG}{$prbTB}{$prbCO} and not equal!\n" .
		"[ERROR]: Prev=@{$desDat{$prbCG}{$prbTB}{$prbCO}}\n" .
		"[ERROR]: Data=@data\n")
		if (! compareArrays(\@{$desDat{$prbCG}{$prbTB}{$prbCO}}, \@data));
	    next;
	}
	die("[ERROR]: Design map already defined: {$prbCG}{$prbFR}{$prbCO}!\n")
	    if (defined($desMap{$prbCG}{$prbFR}{$prbCO}) && $desMap{$prbCG}{$prbFR}{$prbCO} ne $prbTB);

	$desDat{$prbCG}{$prbTB}{$prbCO} = [ @data ];
	$desMap{$prbCG}{$prbFR}{$prbCG} = $prbTB;

        ## Status
        mssg("[Bed]: statusCnt=$inpCnt...\n", 4)
            if (defined($statusCnt) && $inpCnt % $statusCnt == 0);
        $inpCnt++;
        last if (defined($desMax) && $inpCnt >= $desMax);
    }
    close(DES);
    mssg("[Designs]: Loaded $inpCnt designs, total designs=" . scalar(keys %desDat) . "\n\n", 2);
}

sub loadSnv($) {
    my ($file) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    $fileName = trimSuffix($fileName);

    if ($fileSuffix eq ".gz") { open(DES, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(DES, "<$file") or die("Failed to open $file: $!") }
    mssg("[Design]: Loading SNV designFile=$file...\n", 3);

    my $inpCnt = 0;
    while(<DES>) {
        s/^\s+//;
        s/\s+$//;
	next if /^Probeset_ID/;

        my $line = $_;
        my @data = split(/\t/, $line);

        my $prbCG  = $data[UCLA_SEQ_ID];
	my $prbTB  = substr($data[UCLA_ALLELE_TB_STRAND], 0, 1);
        my $prbFR  = $data[UCLA_ALLELE_FR_STRAND];
        my $prbCO  = $data[UCLA_ALLELE_CO_STRAND];
	my $prbTP  = infType($data[UCLA_INFINIUM_TYPE]);

	die("[ERROR]: Multiple definitions for UCLA SNV data={$prbCG}{$prbFR}{$prbCO}{$prbTP}!\n")
	    if (defined($snvDat{$prbCG}{$prbTB}{$prbCO}{$prbTP}));

	$snvDat{$prbCG}{$prbTB}{$prbCO}{$prbTP} = [ @data ];

	if ($prbCG eq "cg26121053") { die("[ERROR]: @{$snvDat{$prbCG}{$prbTB}{$prbCO}{$prbTP}}\n"); }

        ## Status
        mssg("[Bed]: statusCnt=$inpCnt...\n", 4)
            if (defined($statusCnt) && $inpCnt % $statusCnt == 0);
        $inpCnt++;
        last if (defined($desMax) && $inpCnt >= $desMax);
    }
    close(DES);
    mssg("[Designs]: Loaded $inpCnt SNV designs, total SNV designs=" . scalar(keys %desDat) . "\n\n", 2);
}

sub addOrder($$$) {
    my ($ordKeyD, $ordSeqD, $datRef) = @_;

    my ($cg, $tb, $co, $ii, $mu, $cig ) = split($L, $ordKeyD);
    if (defined($ordMap{$ordSeqD}) && $ordMap{$ordSeqD} ne $ordKeyD) {
	$dupMap{$ordSeqD} = $ordKeyD;
	$dupMap{$ordKeyD} = $ordSeqD;
	$dupDat{$ordKeyD} = [ @{$datRef} ];
	mssg("[Warning]: Probe Map Key: ordKeyD already defined:\n" .
	     "[Warning]: $ordMap{$ordSeqD}\n" .
	     "[Warning]: $ordKeyD\n" .
	     "[Warning]: Data=@{$datRef}\n\n", 3);
	return(FALSE);
    }
    if (defined($ordMap{$ordKeyD}) && $ordMap{$ordKeyD} ne $ordSeqD) {
	
	mssg("[Warning]: Probe Map Seq: ordSeqD already defined:\n" .
	     "[Warning]: $ordMap{$ordKeyD}\n" .
	     "[Warning]: $ordSeqD\n" .
	     "[Warning]: Data=@{$datRef}\n\n");
	return(FALSE);
    }
    
    die("[ERROR]: Multiple different order defintions: $ordKeyD\n" .
	"[ERROR]: Prev=@{$ordDat{$ordKeyD}}\n" .
	"[ERROR]: Data=@{$datRef}\n\n")
	if (defined($ordDat{$ordKeyD}) && ! compareArrays(\@{$ordDat{$ordKeyD}}, \@{$datRef}));

    die("[ERROR]: Order data key map already defined: {$cg}{$tb}{$co}{$ii}{$mu}{$cig}\n" .
	"[ERROR]: Prev=$ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig}\n" .
	"[ERROR]: Data=$ordKeyD\n\n")
	if (defined($ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig}) && $ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig} ne $ordKeyD);

    if (! defined($ordDat{$ordKeyD}) && ! defined($ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig})) {
	$ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig} = $ordKeyD;

	$ordMap{$ordSeqD} = $ordKeyD;
	$ordMap{$ordKeyD} = $ordSeqD;
	$ordDat{$ordKeyD} = [ @{$datRef} ];
	return(TRUE);
    } elsif (! defined($ordDat{$ordKeyD}) && defined($ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig}) ||
	     defined($ordDat{$ordKeyD}) && ! defined($ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig})) {

	mssg("[ERROR]: ordDat Not defined!\n") if (! defined($ordDat{$ordKeyD}));
	mssg("[ERROR]: ordDat Defined!\n") if (defined($ordDat{$ordKeyD}));
	mssg("[ERROR]: ordKeys Not defined!\n") if (! defined($ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig}));
	mssg("[ERROR]: ordKeys Defined!\n") if (defined($ordKeys{$cg}{$tb}{$co}{$ii}{$mu}{$cig}));

	die("[ERROR]: Somehow here: ($cg, $tb, $co, $ii, $mu, $cig)\n" .
	    "[ERROR]: ordKeys: {$cg}{$tb}{$co}{$ii}{$mu}{$cig} = $ordKeyD\n" .
	    "[ERROR]: ordMap: $ordMap{$ordSeqD} = $ordKeyD\n" .
	    "[ERROR]: odrMap: $ordMap{$ordKeyD} = $ordSeqD\n" .
	    "[ERROR]: ordDat: \@{$ordDat{$ordKeyD}} = [ @{$datRef} ];\n\n");
    }
    return(-1);
}

sub loadOrder($) {
    my ($file) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    $fileName = trimSuffix($fileName);

    if ($fileSuffix eq ".gz") { open(ORD, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(ORD, "<$file") or die("Failed to open $file: $!") }
    mssg("[Order]: Loading orderFile=$file...\n", 3);

    my $prbCnt = 0;
    my $dupCnt = 0;
    my $repCnt = 0;
    my $totCnt = 0;
    my $inpCnt = 0;
    my $readingData = FALSE;
    while(<ORD>) {
        s/^\s+//;
        s/\s+$//;

        my $line = $_;
        my @data = split(/,/, $line);

        if ($line =~ m/^Assay_Design_Id/) { $readingData = TRUE; next; }

        if ($readingData) {
            $data[0] =~ s/_deg1$//;
            $data[1] =~ s/_deg1_A$/_A/;
            $data[3] =~ s/_deg1_B$/_B/;

            my $orgKey  = $data[ORD_KEY];
            my $orgKeyU = $data[ORD_KEY_A];
            my $ordSeqU = $data[ORD_SEQ_A];
            my $orgKeyM = $data[ORD_KEY_B];
            my $ordSeqM = $data[ORD_SEQ_B];
            my $normBin = $data[ORD_BIN];

            $orgKey =~ s/-37//;

            my ($ordCG, $ordTB, $ordFR, $ordCO);
            if ($orgKey =~ m/^([^_]+)_([^_]+)_([^_]+)_([^_]+)/) {
                $ordCG = $1;
                $ordTB = $2;
                $ordFR = $3;
                $ordCO = $4;
            } else {
                mssg("[Warning]: Unable to parse order key; will skip for now: $orgKey! Skipping...\n", 4);
                next;
            }
	    my $ordTP  = $II;
            my $ordCIG = getCigarDeg($ordSeqU);
            $ordTP = $I if (defined($ordSeqM) && length($ordSeqM) != 0);
	    my $ordKey = join($L, $ordCG, $ordTB, $ordCO, $ordTP);

	    my ($addRet);
	    if ($ordTP eq $II) {
		my $ordSeqD = $ordSeqU;
		my $cigarD  = getCigarDeg($ordSeqD);
		my $ordKeyD = join($L, $ordKey, $D, $cigarD);

		$addRet = addOrder($ordKeyD, $ordSeqD, \@data);
		$prbCnt++ if ($addRet == TRUE);
		$dupCnt++ if ($addRet == FALSE);
		$repCnt++ if ($addRet == -1);
		$totCnt++;
	    } else {
		my $cigarU  = getCigarDeg($ordSeqU);
		my $ordKeyU = join($L, $ordKey, $U, $cigarU);
		my $cigarM  = getCigarDeg($ordSeqM);
		my $ordKeyM = join($L, $ordKey, $M, $cigarM);

		$addRet = addOrder($ordKeyU, $ordSeqU, \@data);
		$prbCnt++ if ($addRet == TRUE);
		$dupCnt++ if ($addRet == FALSE);
		$repCnt++ if ($addRet == -1);
		$totCnt++;

		$addRet = addOrder($ordKeyM, $ordSeqM, \@data);
		$prbCnt++ if ($addRet == TRUE);
		$dupCnt++ if ($addRet == FALSE);
		$repCnt++ if ($addRet == -1);
		$totCnt++;
	    }

	    ## Status
	    mssg("[Bed]: statusCnt=$inpCnt...\n", 4)
		if (defined($statusCnt) && $inpCnt % $statusCnt == 0);
	    $inpCnt++;
	    last if (defined($desMax) && $inpCnt >= $desMax);
	}
    }
    close(ORD);
    my $prbPer = 0; my $dupPer = 0;
    $prbPer = sprintf("%.2f", (100*$prbCnt/$totCnt)) if ($totCnt != 0);
    $dupPer = sprintf("%.2f", (100*$dupCnt/$totCnt)) if ($totCnt != 0);
    mssg("[Order]: Unique Probes Add Rate   = $prbPer ($prbCnt/$totCnt)\t" .
	 "ordDat=".scalar(keys %ordDat).$TAB."ordMap=".scalar(keys %ordMap).$RET, 2);
    mssg("[Order]: Duplicate Probes Add Rate= $dupPer ($dupCnt/$totCnt)\t" .
	 "dupDat=".scalar(keys %dupDat).$TAB."dupMap=".scalar(keys %dupMap).$RET, 2);
    mssg("[Order]: Order Paired-Design Count=$inpCnt\n\n",2);
}

sub loadMatch($) {
    my ($file) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    $fileName = trimSuffix($fileName);

    if ($fileSuffix eq ".gz") { open(ORD, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(ORD, "<$file") or die("Failed to open $file: $!") }
    mssg("[Match]: Loading matchFile=$file...\n", 3);

    my $prbCnt = 0;
    my $inpCnt = 0;
    my $readingData = FALSE;
    while(<ORD>) {
        s/^\s+//;
        s/\s+$//;
	next if /^probe_id/;

        my $line = $_;
        my @data = split(/\t/, $line);

	$data[0] =~ s/_deg1$//;
	$data[1] =~ s/_deg1_A$/_A/;
	$data[3] =~ s/_deg1_B$/_B/;

	my $orgKey  = $data[0];
	my $prbSeq  = $data[1];
	my $address = $data[3];

	$orgKey =~ s/-37//;

	my %ordPrbs = ();
	my ($ordCG, $ordTB, $ordFR, $ordCO);
	if ($orgKey =~ m/^([^_]+)_([^_]+)_([^_]+)_([^_]+)/) {
	    $ordCG = $1;
	    $ordTB = $2;
	    $ordFR = $3;
	    $ordCO = $4;
	} else {
	    mssg("[Warning]: Unable to parse order key; will skip for now: $orgKey! Skipping...\n", 4);
	    next;
	}

	if (defined($ordMap{$prbSeq})) {
	    my $prbKey = $ordMap{$prbSeq};

	    $pqcMap{$address} = $prbKey;
	    $pqcMap{$prbKey} = $address;
	} elsif (defined($dupMap{$prbSeq})) {
	    my $prbKey = $dupMap{$prbSeq};

	    $pqcMap{$address} = $prbKey;
	    $pqcMap{$prbKey} = $address;
	} else {
	    die("[ERROR]: Failed to find match probe in order data: {$ordCG}{$ordTB}{$ordCO}; probe=$prbSeq\n")
	}

	## Status
	mssg("[Bed]: statusCnt=$inpCnt...\n", 4)
	    if (defined($statusCnt) && $inpCnt % $statusCnt == 0);
	$inpCnt++;
	last if (defined($desMax) && $inpCnt >= $desMax);
    }
    close(ORD);
    mssg("[Match]: Loaded $inpCnt match data!\n\n", 2);
}

sub loadPQC($) {
    my ($file) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    $fileName = trimSuffix($fileName);

    if ($fileSuffix eq ".gz") { open(PQC, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(PQC, "<$file") or die("Failed to open $file: $!") }
    mssg("[PQC]: Loading PQC File=$file...\n", 3);

    my $passCnt = 0;
    my $failCnt = 0;
    my $codeCnt = 0;

    my $matCnt = 0;
    my $misCnt = 0;

    my $prbCnt = 0;
    my $inpCnt = 0;
    my $readingData = FALSE;
    while(<PQC>) {
        s/^\s+//;
        s/\s+$//;

        my $line = $_;
        my @data = split(/\t/, $line);

	if ($line =~ m/^Address/) { $readingData = TRUE; next; }

	if ($readingData) {
	    my $code = $data[0];
	    my $pass = $data[1];

	    $pass = 1 if ($data[1] == 0);
	    $pass = 0 if ($data[1] != 0);
	    $passCnt++ if ($pass);
	    $failCnt++ if (! $pass);
	    $codeCnt++;

	    if (! defined($pqcMap{$code})) {
		$matCnt++;
		$pqcDat{$code} = $pass;
	    } else {
		$misCnt++;
		#mssg("[Warning]: PQC result, but no match record! Address=$code, Value=$pass\n");
	    }
	}
    }
    close(PQC);
    my $passPer = 0; my $failPer = 0; my $matPer = 0; my $misPer = 0;
    $passPer = sprintf("%.2f", (100*$passCnt/$codeCnt)) if ($codeCnt != 0);
    $failPer = sprintf("%.2f", (100*$failCnt/$codeCnt)) if ($codeCnt != 0);
    $matPer  = sprintf("%.2f", (100*$matCnt/$codeCnt))  if ($codeCnt != 0);
    $misPer  = sprintf("%.2f", (100*$misCnt/$codeCnt))  if ($codeCnt != 0);
    mssg("[PQC]: Pass  Rate = ".pad($passPer,1,6)." ($passCnt/$codeCnt)\n", 3);
    mssg("[PQC]: Fail  Rate = ".pad($failPer,1,6)." ($failCnt/$codeCnt)\n", 3);
    mssg("[PQC]: Match Rate = ".pad($matPer,1,6)." ($matCnt/$codeCnt)\n", 3);
    mssg("[PQC]: Miss  Rate = ".pad($misPer,1,6)." ($misCnt/$codeCnt)\n", 3);
    mssg("\n",3);
}

sub loadDesignFiles($) {
    my ($filesRef) = @_;
    foreach my $file (@{$filesRef}) {
	loadDesign($file);
    }
}

sub loadSnvFiles($) {
    my ($filesRef) = @_;
    foreach my $file (@{$filesRef}) {
	loadSnv($file);
    }
}

sub loadOrderFiles($) {
    my ($filesRef) = @_;
    foreach my $file (@{$filesRef}) {
	loadOrder($file);
    }
}

sub loadMatchFiles($) {
    my ($filesRef) = @_;
    foreach my $file (@{$filesRef}) {
	loadMatch($file);
    }
}

sub loadQCFiles($$) {
    my ($aqpRef, $pqcRef) = @_;

    if (scalar(@{$pqcRef}) != 0) {
	foreach my $file (@{$pqcRef}) {
	    loadPQC($file);
	}
    } elsif (scalar(@{$aqpRef}) != 0) {
	foreach my $file(@{$aqpRef}) {
	    loadAQP($file);
	}
    } else {
	die("[ERROR]: No AQP or PQC data provided!");
    }
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Main:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

GetOptions("o|out=s"     => \$outDir,
	   "d|des=s@"    => \@desFiles,
	   "s|snv=s@"    => \@snvFiles,
	   "m|mat=s@"    => \@matFiles,
	   "p|pqc=s@"    => \@pqcFiles,
	   "man=s@"      => \@manFiles,


	   "ord=s@"      => \@ordFiles,

	   "desMax=i"    => \$desMax,
	   "snvMax=i"    => \$snvMax,
	   "manMax=i"    => \$manMax,

	   "v|verbose=i" => \$verbose,
    );

my @curDate = getTime();
print STDERR "[Date]: $curDate[0]\n\n";

if (! defined($outDir) || @desFiles == 0 || @snvFiles == 0) {
    print STDERR "[Usage]: $0 -o outDir\n"
	. " -d|dat designFiles(s)\n"
	. " -s|snv taretFile(s)\n"
	. " -p|pqc pqcFile(s)\n"
	. " [ -m|man manifest(s) ]\n"
	. " [ -desMax designMax ]\n"
	. " [ -snvMax snvMax ]\n"
	. " [ -manMax manifestMax ]\n"
	. " [ v|verbose verbosityLevel ]\n"
	. " [ options ]\n\n";
    print STDERR "[Usage]: Will read inptus from STDIN if not specified.\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}

setParams();
buildDirs();

loadDesignFiles(\@desFiles);
loadSnvFiles(\@snvFiles);

loadOrderFiles(\@ordFiles);
loadMatchFiles(\@matFiles);
loadQCFiles(\@aqpFiles, \@pqcFiles);



##
## Done!
##
@curDate = getTime();
mssg("[Done]: $curDate[0], $curDate[2]\n\n");


## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Basic Subroutines:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub setParams {
    $topDir    = abs_path(getcwd);
    $statusCnt = 10000;
    $qsub = "qsub -cwd -pe threaded 16 -l excl=true -N";

    $A="A"; $B="B"; $C="C"; $D="D"; $E="E"; $F="F"; $G="G"; $H="H"; $I="I"; $J="J"; $K="K"; $L="L"; $M="M"; $N="N";
    $O="O"; $P="P"; $Q="Q"; $R="R"; $S="S"; $T="T"; $U="U"; $V="V"; $W="W"; $X="X"; $Y="Y"; $Z="Z";

#    $a="a"; $b="b"; $c="c"; $d="d"; $e="e"; $f="f"; $g="g"; $h="h"; $i="i"; $j="j"; $k="k"; $l="l"; $m="m"; $n="n";
#    $o="o"; $p="p"; $q="q"; $r="r"; $s="s"; $t="t"; $u="u"; $v="v"; $w="w"; $x="x"; $y="y"; $z="z";

    $I="I"; $II="II";
    $L="_"; $SL="/";

    @ACTG = ($A, $C, $T, $G); @FR  = ($F, $R); @CO  = ($C, $O); @MU  = ($M, $U); @III = ($I, $II);
    foreach my $nuc (@ACTG) { $ACTG{$nuc} = $nuc; }

    $cmplChar{$F} = $R; $cmplChar{$R} = $F; ## F/R
    $cmplChar{$C} = $O; $cmplChar{$O} = $C; ## C/O
    $cmplChar{$M} = $U; $cmplChar{$U} = $M; ## M/U
    $cmplChar{$T} = $B; $cmplChar{$B} = $T; ## T/B
    $cmplChar{$D} = $D; ## Degnerate???

    my @letters = ($A, $B, $C, $D, $E, $F, $G, $H, $I, $J, $K, $L, $M, $N,
		   $O, $P, $Q, $R, $S, $T, $U, $V, $W, $X, $Y, $Z);

#		   $a, $b, $c, $d, $e, $f, $g, $h, $i, $j, $k, $l, $m, $n,
#		   $o, $p, $q, $r, $s, $t, $u, $v, $w, $x, $y, $z);

    foreach my $let (@letters) { $cmplNucs{$let} = cmpl($let); }

    push @{$iupacNucs{$R}}, $A; push @{$iupacNucs{$Y}}, $T;
    push @{$iupacNucs{$R}}, $G; push @{$iupacNucs{$Y}}, $C;

    push @{$iupacNucs{$S}}, $C; push @{$iupacNucs{$W}}, $A;
    push @{$iupacNucs{$S}}, $G; push @{$iupacNucs{$W}}, $T;

    push @{$iupacNucs{$K}}, $G; push @{$iupacNucs{$M}}, $A;
    push @{$iupacNucs{$K}}, $T; push @{$iupacNucs{$M}}, $C;

    push @{$iupacNucs{$B}}, $C; push @{$iupacNucs{$D}}, $A; push @{$iupacNucs{$H}}, $A; push @{$iupacNucs{$V}}, $A;
    push @{$iupacNucs{$B}}, $G; push @{$iupacNucs{$D}}, $G; push @{$iupacNucs{$H}}, $C; push @{$iupacNucs{$V}}, $C;
    push @{$iupacNucs{$B}}, $T; push @{$iupacNucs{$D}}, $T; push @{$iupacNucs{$H}}, $T; push @{$iupacNucs{$V}}, $G;

    push @{$iupacNucs{$N}}, $A;
    push @{$iupacNucs{$N}}, $C;
    push @{$iupacNucs{$N}}, $T;
    push @{$iupacNucs{$N}}, $G;

    $TAB="\t"; $COM="\t"; $COL=":"; $RET="\n"; $SEM=";"; 
    $PLUS="+"; $MINUS="-"; $DASH="-";

    foreach my $nuc (@ACTG) { $basicNucs{$nuc} = $nuc; }
    foreach my $deg (sort keys %iupacNucs) { my $nucs = ""; foreach my $nuc (sort @{$iupacNucs{$deg}}) { $nucs .= $nuc; } $iupacStrs{$nucs} = $deg; }
}

sub buildDirs {
    my $cleanAll = shift;

    $outDir =~ s/\/+$//;
    system("rm -rf $outDir") if (defined($cleanAll) && $cleanAll && -e $outDir);

    make_path($outDir);
    $outDir = abs_path($outDir);
    mssg("[Input]: outDir=$outDir\n\n");

    $logDir = "$outDir/log";
    make_path($logDir) if (! -e $logDir);
    system("rm -f $logDir/* &> $outDir/cleanLog.log");

#    $bedDir = "$outDir/bed";
#    make_path($bedDir) if (! -e $bedDir);
#    system("rm -f $bedDir/* &> $logDir/cleanBed.log");

    $vizDir = "$outDir/vizifest";
    make_path($vizDir) if (! -e $vizDir);
    system("rm -f $vizDir/* &> $logDir/cleanVizifest.log");

    $manDir = "$outDir/manifest";
    make_path($manDir) if (! -e $manDir);
    system("rm -f $manDir/* &> $logDir/cleanManifest.log");

#    $ordDir = "$outDir/order";
#    make_path($ordDir) if (! -e $ordDir);
#    system("rm -f $ordDir/* &> $logDir/cleanOrd.log");

    $desDir = "$outDir/design";
    make_path($desDir) if (! -e $desDir);
    system("rm -f $desDir/* &> $logDir/cleanDes.log");

    $fasDir = "$outDir/fas";
    make_path($fasDir) if (! -e $fasDir);
    system("rm -f $fasDir/* &> $logDir/cleanFasta.log");

    #$bamDir = "$outDir/bam";
    #make_path($bamDir) if (! -e $bamDir);
    #system("rm -f $bamDir/* &> $logDir/cleanBam.log");
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
    my $str = shift;
    my $level = shift;

    if (defined($level)) {
	print STDERR $str if ($verbose >= $level);
    } else {
	print STDERR $str if ($verbose);
    }
}

sub arrayIterStr {
    my ($datRef, $prefix) = @_;

    $prefix = "" if (! defined($prefix));
    my $str .= "";
    for (my $ii = 0; $ii < @{$datRef}; $ii++) {
	$str .= $prefix ."$ii='$$datRef[$ii]'\n";
    }
    return($str);
}

##
## Quick Math
##

sub min {
    my ($a, $b) = @_;
    return($a) if ($a < $b);
    return($b);
}

sub max {
    my ($a, $b) = @_;
    return($a) if ($a > $b);
    return($b);
}

sub getAvg {
    my ($ref) = @_;

    my @dat = @{$ref};
    my $sum = 0;
    foreach my $val (@dat) {
        next if (! defined($val) || length($val) == 0 || $val =~ m/^(\s)*$/);
        $sum += $val;
    }
    my $avg = $sum / scalar(@dat);
    return($avg);
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

##
## Quick I/O:
##

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

    my $cmd = "ls $dir/*$key";
    open(LIST, "$cmd | ") or die("Failed to open stream: $cmd for reading: $!");
    my @list = readStream(\*LIST);
    close(LIST);
    return(@list);
}

sub trimSuffix {
    my ($str) = @_;

    for (my $ii = 0; $ii < 5; $ii++) {
    $str =~ s/\.gz//; $str =~ s/\.fa//;
    $str =~ s/\.bed//; $str =~ s/\.txt//; $str =~ s/\.csv//; $str =~ s/\.tsv//;
    $str =~ s/\.bam//; $str =~ s/\.sam//; $str =~ s/\.vcf//; $str =~ s/\.fas//;
    $str =~ s/\.sorted//; $str =~ s/\.merged//; $str =~ s/\.intersect//;
    }
    return($str);
}

##
## Common file formats:
##

sub bgzipIdx {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $sortedBed = "$filePath$fileName.sorted.bed";

    my $cmd = "sort -k 1,1V -k 2,2n -k 3,3n $file > $sortedBed";
    die("[bgzipIdx]: Failed sort: '$cmd'") if (system("$cmd"));
    system("rm -f $sortedBed.gz") if (-e "$sortedBed.gz");
    system("rm -f $sortedBed.gz.tbi") if (-e "$sortedBed.gz.tbi");
    die("[BgzipIdx]: Failed bgzip $sortedBed") if (system("bgzip $sortedBed"));
    $sortedBed .= ".gz";
    die("[BgzipIdx]: Failed tabix $sortedBed") if (system("tabix $sortedBed"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub samToBam {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $bam = "$filePath$fileName.bam";
    my $sortedBam = "$filePath$fileName.sorted.bam";

    my $cmd = "samtools view -S -b $file > $bam";
    die("[samToBam]: Failed sam <-> bam: $cmd") if (system($cmd));
    $cmd = "samtools sort $bam > $sortedBam";
    die("[samToBam]: Failed sort: $cmd") if (system($cmd));
    $cmd = "samtools index $sortedBam";
    die("[samToBam]: Failed index: $cmd") if (system($cmd));
    system("rm -f $file") if ($clean);
    system("rm -f $bam") if ($clean);
}

sub compareArrays {
    my ($a, $b) = @_;

    return(FALSE) if (@{$a} != @{$b});
    for (my $ii = 0; $ii < @{$a}; $ii++) {
	return(FALSE) if ($$a[$ii] ne $$b[$ii]);
    }
    return(TRUE);
}

##
## Program subroutines
##
sub cmpl($) {
    my ($seq) = @_;
    $seq =~ tr/ACTGRYSWKMBDHV/TGACYRSWMKVHDB/;
    $seq =~ tr/actgryswkmbdhv/tgacyrswmkvhdb/;
    return($seq);
}

sub revCmp($) {
    my ($seq) = @_;
    return(reverse(cmpl($seq)));
}

sub bscGeneral {
    my ($snvSeq) = @_;

    die("[ERROR]: bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) != 0);
    my $bscSeq = ""; $snvSeq = uc($snvSeq);

    for (my $ii = 0; $ii < length($snvSeq)-1; $ii++) {
	my $snvCurr = substr($snvSeq, $ii+0, 1);
	my $snvNext = substr($snvSeq, $ii+1, 1);

	if (uc($snvCurr) =~ m/[CY]/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($Y);
	} elsif (uc($snvCurr) =~ m/([S])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($B);
	} elsif (uc($snvCurr) =~ m/([M])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($H);
	} elsif (uc($snvCurr) =~ m/([B])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($B);
	} elsif (uc($snvCurr) =~ m/([H])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($H);
	} elsif (uc($snvCurr) =~ m/([V])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($N);


	} elsif (uc($snvCurr) =~ m/([CY])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $T;
	} elsif (uc($snvCurr) =~ m/([S])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/([M])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/([B])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/([H])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/([V])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $D;
	} elsif (uc($snvCurr) =~ m/([CYSMBHV])/ && uc($snvNext) =~ m/([ATYWMH])/) {
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
    return($bscSeq);
}

sub bscM {
    my ($snvSeq) = @_;

    die("[ERROR]: bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) != 0);
    my $bscSeq = ""; $snvSeq = uc($snvSeq);

    for (my $ii = 0; $ii < length($snvSeq)-1; $ii++) {
	my $snvCurr = substr($snvSeq, $ii+0, 1);
	my $snvNext = substr($snvSeq, $ii+1, 1);

	if (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/([G])/) {
	    $bscSeq .= lc($C);
	} elsif (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/([RSKBDV])/) {
	    $bscSeq .= lc($Y);
	} elsif (uc($snvCurr) =~ m/[Y]/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($Y);
	} elsif (uc($snvCurr) =~ m/([S])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($S);
	} elsif (uc($snvCurr) =~ m/([M])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($M);
	} elsif (uc($snvCurr) =~ m/([B])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($B);
	} elsif (uc($snvCurr) =~ m/([H])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($H);
	} elsif (uc($snvCurr) =~ m/([V])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($V);


	} elsif (uc($snvCurr) =~ m/([CY])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $T;
	} elsif (uc($snvCurr) =~ m/([S])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/([M])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/([B])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/([H])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/([V])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $D;

	} elsif (uc($snvCurr) =~ m/([CYSMBHV])/ && uc($snvNext) =~ m/([ATYWMH])/) {
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
    return($bscSeq);
}

sub bscU {
    my ($snvSeq, $refSeq, $prbFR, $prbCO) = @_;

    die("[ERROR]: bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) != 0);
    my $bscSeq = ""; $snvSeq = uc($snvSeq);

    for (my $ii = 0; $ii < length($snvSeq)-1; $ii++) {
	my $snvCurr = substr($snvSeq, $ii+0, 1);
	my $snvNext = substr($snvSeq, $ii+1, 1);

	if (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/([G])/) {
	    $bscSeq .= lc($T);
	} elsif (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/([RSKBDV])/) {
	    $bscSeq .= lc($T);
	} elsif (uc($snvCurr) =~ m/[Y]/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($T);
	} elsif (uc($snvCurr) =~ m/([S])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($K);
	} elsif (uc($snvCurr) =~ m/([M])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($W);
	} elsif (uc($snvCurr) =~ m/([B])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($K);
	} elsif (uc($snvCurr) =~ m/([H])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($W);
	} elsif (uc($snvCurr) =~ m/([V])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $bscSeq .= lc($D);


	} elsif (uc($snvCurr) =~ m/([CY])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $T;
	} elsif (uc($snvCurr) =~ m/([S])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/([M])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/([B])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $K;
	} elsif (uc($snvCurr) =~ m/([H])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $W;
	} elsif (uc($snvCurr) =~ m/([V])/ && uc($snvNext) =~ m/([ATCYWMH])/) {
	    $bscSeq .= $D;

	} elsif (uc($snvCurr) =~ m/([CYSMBHV])/ && uc($snvNext) =~ m/([ATYWMH])/) {
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
    return($bscSeq);
}

sub bsc($) {
    my ($seq) = @_;

    $seq = uc($seq);
    $seq =~ s/CG/yg/gi;
    $seq =~ s/C/T/gi;
    return($seq);
}

sub shearBrac($) {
    my ($seq) = @_;
    $seq =~ s/\[//; $seq =~ s/\]//;
    return($seq);
}

sub doubleArray($$) {
    my ($datRef, $x) = @_;

    my @ret = ();
    foreach my $seq (@{$datRef}) {
        for (my $ii = 0; $ii < $x; $ii++) {
            push @ret, $seq;
        }
    }
    return(@ret);
}

sub expandSeq($) {
    my ($orgSeq) = @_;

    my @seqDat = split(//, $orgSeq);

    my @seqs = ($orgSeq);
    my @subDat = ();
    my $forcast = 1;
    for (my $ii = 0; $ii < @seqDat; $ii++) {
	my $seqNuc = uc($seqDat[$ii]);
	mssg("[expandSeq]: seqNuc=$seqDat[$ii] <=> $seqNuc\n",$lev5) if (! defined($ACTG{$seqNuc}));
        if (defined($iupacNucs{$seqNuc})) {
            my $x = scalar(@{$iupacNucs{$seqNuc}});
	    $forcast *= $x;
	}
    }
    if ($forcast > $forcastLim) { $forcastFlag = $forcast; $forcastCnt++; return(@seqs); }

    for (my $ii = 0; $ii < @seqDat; $ii++) {
	my $seqNuc = uc($seqDat[$ii]);
	mssg("[expandSeq]: seqNuc=$seqDat[$ii] <=> $seqNuc\n",$lev5) if (! defined($ACTG{$seqNuc}));
        if (defined($iupacNucs{$seqNuc})) {
            my $x = scalar(@{$iupacNucs{$seqNuc}});
            @seqs = doubleArray(\@seqs, $x);
            for (my $seqIdx = 0; $seqIdx < @seqs; $seqIdx++) {
                my $degIdx = $seqIdx % $x;
                my $degNuc = lc($iupacNucs{$seqNuc}[$degIdx]);
                substr($seqs[$seqIdx], $ii, 1, $degNuc);
            }
        }
    }
    return(@seqs);
}

sub expandSet {
    my ($prbKey, $degPrb, $orgPrb, $datRef) = @_;

    die("[ERROR]: expandSet: prbKey not defined!") if (! defined($prbKey) || length($prbKey) == 0);
    die("[ERROR]: expandSet: degPrb not defined!") if (! defined($degPrb) || length($degPrb) == 0);
    die("[ERROR]: expandSet: orgPrb not defined!") if (! defined($orgPrb) || length($orgPrb) == 0);
    die("[ERROR]: expandSet: datRef not defined!") if (! defined($datRef));

    my @degPrbs = expandSeq($degPrb);
    foreach my $prb (@degPrbs) {
	my $degKey = join($L, $prbKey, getCigarDeg($prb));
	my $snvKey = join($L, $prbKey, getCigarRef($orgPrb, $prb));
	$$datRef{uc($prb)} = [ $prb, $degKey, $snvKey ];
    }
}

sub getCigarDeg {
    my ($ref) = @_;

    die("[ERROR]: getCigarDeg: ref not defined!") if (! defined($ref) || length($ref) == 0);

    my @refDat = split(//, uc($ref));
    my $matCnt = 0;
    my $cigar = "";
    for (my $ii = 0; $ii < @refDat; $ii++) {
        if (defined($ACTG{$refDat[$ii]})) {
            $matCnt++;
        } else {
            $cigar .= $matCnt ."M" if ($matCnt != 0);
            $cigar .= "$refDat[$ii]";
            $matCnt = 0;
        }
    }
    $cigar .= $matCnt . "M" if ($matCnt != 0);
    return($cigar);
}

sub getCigarRef {
    my ($ref, $can) = @_;

    die("[ERROR]: getCigarRef: ref not defined!") if (! defined($ref) || length($ref) == 0);
    die("[ERROR]: getCigarRef: can not defined!") if (! defined($can) || length($can) == 0);

    my @refDat = split(//, uc($ref));
    my @canDat = split(//, uc($can));

    die("[ERROR]: getCigarRef: Candidate and reference length are not of equal length!\n" .
	"[ERROR]: CAN=$can\n" .
	"[ERROR]: REF=$ref\n") if (@canDat != @refDat);

    my $matCnt = 0;
    my $cigar = "";
    for (my $ii = 0; $ii < @canDat; $ii++) {
        if ($refDat[$ii] eq $canDat[$ii]) { $matCnt++; next; }
        if ($refDat[$ii] ne $canDat[$ii]) {
            $cigar .= $matCnt. "M" if ($matCnt != 0);
            $cigar .= "$canDat[$ii]";
            $matCnt = 0;
        }
    }
    $cigar .= $matCnt ."M" if ($matCnt != 0);
    return($cigar);
}

sub getCigarsRef {
    my ($ref, $canRef, $prefix) = @_;

    die("[ERROR]: getCigarRef: ref not defined!") if (! defined($ref) || length($ref) == 0);
    die("[ERROR]: getCigarRef: canRef not defined!") if (! defined($canRef));

    my @ret = ();
    foreach my $can (@{$canRef}) {
	push @ret, getCigarRef($ref, $can) if (! defined($prefix));
	push @ret, join($L, $prefix, getCigarRef($ref, $can)) if (defined($prefix));
    }
    return(@ret);
}

sub infType {
    my ($str) = @_;

    die("[ERROR]: infType: str not defined!") if (! defined($str) || length($str) == 0);

    return($II) if ($str eq "Inf2");
    return($I) if ($str eq "Inf1");
    die("[ERROR]: Unknown Infinium Type during conversions: $str!\n");
}


sub pad {
    my $str  = shift;
    my $side = shift;
    my $tLen = shift;

    $tLen = $pLen if (! defined($tLen));

    if (length($str) > $tLen) { mssg("[Warning]: too long (" .length($str) ."): $str!\n", 5); return($str); }
    $str .= " " x ($tLen - length($str)) if (! $side);
    $str = " " x ($tLen - length($str)) . $str if ($side);
    return($str);
}

sub getDegPrb {
    my ($prbs) = @_;

    die("[ERROR]: getDegPrb: prbs not defined!") if (! defined($prbs));

    my @mat = ();
    my $pLen;
    my $sLen = 0;
    foreach my $prb (@{$prbs}) {
	my @seqDat = split(//, uc($prb));
	push @mat, [ @seqDat ];

	die("[ERROR]: Failed matching lengths! $prb: " . length($prb) . " != $pLen!") if (defined($pLen) && $pLen != length($prb));
	$pLen = length($prb) if (! defined($pLen));
	$sLen++;
    }
    mssg("[SLen]: $sLen\n",$lev5);
    my $degPrb = "";
    for (my $ii = 0; $ii < $pLen; $ii++) {
	my %nucs = ();
	for (my $jj = 0; $jj < $sLen; $jj++) {
	    my $nuc = $mat[$jj][$ii];
	    die("[ERROR]: Unknown nuc=$nuc, not found in basic or iupac!\n" .
		"[ERROR]: Mat=$mat[$ii][$jj+1]\n") if (! defined($basicNucs{$nuc}) && ! defined($iupacNucs{$nuc}));
	    if (defined($basicNucs{$nuc})) { $nucs{$nuc} = 1; next; }
	    foreach my $iNuc (@{$iupacNucs{$nuc}}) { $nucs{$iNuc}++; };
	}
	my $nucStr = ""; foreach my $nuc (sort keys %nucs) { $nucStr .= $nuc; }

	if (defined($basicNucs{$nucStr})) { $degPrb .= $basicNucs{$nucStr}; next; }
	die("[ERROR]: Unable to find iupac for $nucStr!") if (! defined($iupacStrs{$nucStr}));
	$degPrb .= $iupacStrs{$nucStr};
    }
    return($degPrb);
}

sub missMatchStr {
    my ($can, $ref) = @_;

    die("[ERROR]: misMatchStr: can not defined!") if (! defined($can));
    die("[ERROR]: misMatchStr: ref not defined!") if (! defined($ref));

    my @canDat = split(//, uc($can));
    my @refDat = split(//, uc($ref));

    my $matStr = "";
    for (my $ii = 0; $ii < @canDat; $ii++) {
	if ($canDat[$ii] ne $refDat[$ii]) {
	    $matStr .= "X";
	} else {
	    if (! defined($ACTG{$refDat[$ii]})) {
		$matStr .= "*";
	    } else {
		$matStr .= " ";
	    }
	}
    }
    return($matStr);
}

sub getDegCost {
    my ($seq) = @_;

    die("[ERROR]: getDegCost: seq not defined!") if (! defined($seq) || length($seq) == 0);

    my @seqDat = split(//, uc($seq));

    my $cost = 0;
    for (my $ii = 0; $ii < @seqDat; $ii++) {
	my $nuc = $seqDat[$ii];
	$cost += scalar(@{$iupacNucs{$nuc}}) if (defined($iupacNucs{$nuc}));
    }
    return($cost);
}

sub getPrbStats {
    my ($seq) = @_;

    die("[ERROR]: getPrbStats: seq not defined!") if (! defined($seq) || length($seq) == 0);

    my $cpgCnt = getCpgCnt($seq);
    my ($degCnt, $cost) = getDegCnt($seq);
    my @ret = ($cpgCnt, $degCnt, $cost);
    return(@ret);
}

sub getDegCnt {
    my ($snvSeq) = @_;

    die("[ERROR]: getDegCnt: snvSeq not defined!") if (! defined($snvSeq) || length($snvSeq) == 0);

    $snvSeq = uc($snvSeq);
    $snvSeq =~ s/A//gi;
    $snvSeq =~ s/C//gi;
    $snvSeq =~ s/T//gi;
    $snvSeq =~ s/G//gi;
    my $degCnt = length($snvSeq);

    my @nucs = split(//, $snvSeq);
    my $cost = 1;
    foreach my $nuc (@nucs) {
	die("[ERROR]: getDegCnt: Unknown nuc=$nuc in iupac look up! snvSeq=$snvSeq!") if (! defined($iupacNucs{$nuc}));
	$cost *= scalar(@{$iupacNucs{$nuc}});
    }
    my @ret = ($degCnt, $cost);

    return(@ret);
}

sub getCpgCnt {
    my ($snvSeq) = @_;

    die("[ERROR]: getCpgCnt: snvSeq not defined!") if (! defined($snvSeq) || length($snvSeq) == 0);

    my $cpgCnt = 0;
    $snvSeq = uc($snvSeq);
    for (my $ii = 0; $ii < length($snvSeq)-1; $ii++) {
	my $snvCurr = substr($snvSeq, $ii+0, 1);
	my $snvNext = substr($snvSeq, $ii+1, 1);

	if (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/([G])/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/([RSKBDV])/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[Y]/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/([S])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/([M])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/([B])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/([H])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/([V])/ && uc($snvNext) =~ m/([GRSKBDV])/) {
	    $cpgCnt++;
	}
    }
    return($cpgCnt);
}

sub collapseVariants {
    my ($expCrds, $expRefs, $expAlts,
	$snvCrds, $snvRefs, $snvAlts) = @_;

    die("[ERROR]: collapseVariants: expCrds not defined!") if (! defined($expCrds));
    die("[ERROR]: collapseVariants: expRefs not defined!") if (! defined($expRefs));
    die("[ERROR]: collapseVariants: expAlts not defined!") if (! defined($expAlts));

    die("[ERROR]: collapseVariants: snvCrds not defined!") if (! defined($expCrds));
    die("[ERROR]: collapseVariants: snvRefs not defined!") if (! defined($snvRefs));
    die("[ERROR]: collapseVariants: snvAlts not defined!") if (! defined($snvAlts));

    my %crds = (); my %refs = (); my %alfs = ();
    for (my $ii = 0; $ii < @{$expCrds}; $ii++) {
	my $crd = $$expCrds[$ii];
	my $ref = $$expRefs[$ii];
	my $alt = $$expAlts[$ii];
	next if ($crd eq "-");
	next if ($ref eq $N || $alt eq $N);

	## Make sure we're not adding an N: ACTG
	next if (scalar(keys %{$crds{$crd}}) >= 3);

	$refs{$crd} = $ref if (! defined($refs{$crd}));
	$crds{$crd}{$ref} = 1;
	$crds{$crd}{$alt} = 1;
    }

    foreach my $crd (sort {$a <=> $b} keys %crds) {
	push @{$snvCrds}, $crd;
	push @{$snvRefs}, $refs{$crd};
	my $nucs = ""; foreach my $nuc (sort keys %{$crds{$crd}}) { $nucs .= $nuc; }

	die("[ERROR]: Unable to find degenerate base for $nucs!\n"
	    . "\tCrds=@{$expCrds}\n"
	    . "\tRefs=@{$expRefs}\n"
	    . "\tAlts=@{$expAlts}\n") if (! defined($iupacStrs{$nucs}));
	my $deg = $iupacStrs{$nucs};
	push @{$snvAlts}, $deg;
    }
}

sub addVariants {
    my ($fwdSeq, $snvPos, $snvLocs, $snvRefs, $snvAlts, $varDat) = @_;

    die("[ERROR]: addVariants: fwdSeq not defined!")  if (! defined($fwdSeq));
    die("[ERROR]: addVariants: snvPos not defined!")  if (! defined($snvPos));
    die("[ERROR]: addVariants: snvLocs not defined!") if (! defined($snvLocs));
    die("[ERROR]: addVariants: snvRefs not defined!") if (! defined($snvRefs));
    die("[ERROR]: addVariants: snvAlts not defined!") if (! defined($snvAlts));
    die("[ERROR]: addVariants: varDat not defined!")  if (! defined($varDat));

    $fwdSeq  = shearBrac($fwdSeq);
    my @expCrds = split(/,/, $snvLocs);
    my @expRefs = split(/,/, $snvRefs);
    my @expAlts = split(/,/, $snvAlts);

    my @snvCrds = (); my @snvRefs = (); my @snvAlts = ();
    collapseVariants(\@expCrds, \@expRefs, \@expAlts,
		     \@snvCrds, \@snvRefs, \@snvAlts);
    my $fwdStart = $snvPos - 60;

    my @snvIdxs = ();
    my $altStr = " " x 122; my $refStr = " " x 122;
    for (my $ii = 0; $ii < @snvCrds; $ii++) {
        next if ($snvCrds[$ii] eq "-");
        my $snvIdx = $snvCrds[$ii] - $fwdStart;

        die("[ERROR]: snvIdx=$snvIdx < 0, $snvCrds[$ii] - $fwdStart, snvPositions=@snvCrds\n" .
            "[ERROR]: line=@{$snvRefs}!") if ($snvIdx < 0);
        push @snvIdxs, $snvIdx;

	my $srcNuc = substr($fwdSeq, $snvIdx, 1);
        if ($srcNuc ne $snvRefs[$ii]) {
            die("[ERROR]: srcNuc=$srcNuc != snvRef=$snvRefs[$ii], snvIdx=$snvIdx\n"
                . "[ERROR]: line=@{$snvRefs}!\n"
		. "\tCrds: @snvCrds\n\tRefs: @snvRefs\n\tAlts: @snvAlts\n");
        }
        substr($fwdSeq, $snvIdx, 1, $snvAlts[$ii]);
        substr($altStr, $snvIdx, 1, $snvAlts[$ii]);
        substr($refStr, $snvIdx, 1, $snvRefs[$ii]);
    }
    push @{$varDat}, join($SEM, @snvCrds);
    push @{$varDat}, join($SEM, @snvRefs);
    push @{$varDat}, join($SEM, @snvAlts);
    push @{$varDat}, $refStr;
    push @{$varDat}, $altStr;

    return($fwdSeq);
}

sub getProbes {
    my ($desType, $state, $desSeq, $snvSeq, $orgSeq, $bsuSeq, $bsmSeq, 
	$prbKey, $desCG, $desTB, $desFR, $desCO, $datRef) = @_;

    my $prbKeyI  = join($L, $desCG, $desTB, $desFR, $desCO, $I);
    my $prbKeyII = join($L, $desCG, $desTB, $desFR, $desCO, $II);
    my ($degPrbI, $degPrbII, $snvPrbI, $snvPrbII, $orgPrbI, $orgPrbII);
    my ($bsuPrbI, $bsuPrbII, $bsmPrbI, $bsmPrbII);

    if ($desType eq "snp") {
	##
	## TBD: only one of these works:
	##
	$degPrbI  = revCmp(substr($desSeq, 60, 50)) if ($desFR eq $F && $desCO eq $C);
	$degPrbI  = revCmp(substr($desSeq, 61, 50)) if ($desFR eq $F && $desCO eq $O);
	$degPrbI  = revCmp(substr($desSeq, 61, 50)) if ($desFR eq $R && $desCO eq $C); 
	$degPrbI  = revCmp(substr($desSeq, 62, 50)) if ($desFR eq $R && $desCO eq $O);

	$degPrbII = revCmp(substr($desSeq, 60, 50)) if ($desFR eq $F && $desCO eq $C);
	$degPrbII = revCmp(substr($desSeq, 61, 50)) if ($desFR eq $F && $desCO eq $O);
	$degPrbII = revCmp(substr($desSeq, 61, 50)) if ($desFR eq $R && $desCO eq $C);
	$degPrbII = revCmp(substr($desSeq, 62, 50)) if ($desFR eq $R && $desCO eq $O);
    } elsif ($desType eq "cpg") {
	$degPrbI  = revCmp(substr($desSeq, 60, 50)) if ($desCO eq $C);
	$degPrbI  = revCmp(substr($desSeq, 61, 50)) if ($desCO eq $O);
	$degPrbII = revCmp(substr($desSeq, 61, 50)) if ($desCO eq $C);
	$degPrbII = revCmp(substr($desSeq, 62, 50)) if ($desCO eq $O);

	$snvPrbI   = revCmp(substr($snvSeq, 60, 50)) if ($desCO eq $C);
	$snvPrbI   = revCmp(substr($snvSeq, 61, 50)) if ($desCO eq $O);
	$snvPrbII  = revCmp(substr($snvSeq, 61, 50)) if ($desCO eq $C);
	$snvPrbII  = revCmp(substr($snvSeq, 62, 50)) if ($desCO eq $O);

	$orgPrbI   = revCmp(substr($orgSeq, 60, 50)) if ($desCO eq $C);
	$orgPrbI   = revCmp(substr($orgSeq, 61, 50)) if ($desCO eq $O);
	$orgPrbII  = revCmp(substr($orgSeq, 61, 50)) if ($desCO eq $C);
	$orgPrbII  = revCmp(substr($orgSeq, 62, 50)) if ($desCO eq $O);

	$bsuPrbI   = revCmp(substr($bsuSeq, 60, 50)) if ($desCO eq $C);
	$bsuPrbI   = revCmp(substr($bsuSeq, 61, 50)) if ($desCO eq $O);
	$bsuPrbII  = revCmp(substr($bsuSeq, 61, 50)) if ($desCO eq $C);
	$bsuPrbII  = revCmp(substr($bsuSeq, 62, 50)) if ($desCO eq $O);

	$bsmPrbI   = revCmp(substr($bsmSeq, 60, 50)) if ($desCO eq $C);
	$bsmPrbI   = revCmp(substr($bsmSeq, 61, 50)) if ($desCO eq $O);
	$bsmPrbII  = revCmp(substr($bsmSeq, 61, 50)) if ($desCO eq $C);
	$bsmPrbII  = revCmp(substr($bsmSeq, 62, 50)) if ($desCO eq $O);
    } else {
	die("[ERROR]: Unknown Design Type=$desType!");
    }
}

## End of file
