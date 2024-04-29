#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Path qw(make_path remove_tree);
use File::Basename;
#use Data::Dumper::Simple;
use Storable;

my %bin = (
    "T" => '00',
    "G" => '01',
    "A" => '10',
    "C" => '11',
    );
foreach my $key (keys %bin) { $bin{$bin{$key}} = $key; }

my @runTimes = ();
my ($endTime, $runTime, $curTime, $curDate);
my $startTime=time; $endTime=$startTime;
use constant FALSE => 0;
use constant TRUE  => 1;

## Set Format:
use constant SET_REF_SEQS => 0;
use constant SET_REF_PRBS => 1;
use constant SET_REF_EXPS => 2;
use constant SET_SNV_SEQS => 3;
use constant SET_SNV_PRBS => 4;
use constant SET_SNV_EXPS => 5;
use constant SET_OUT_STR  => 6;
use constant SET_REF_KEY  => 7;
use constant SET_DEG_KEY  => 8;

## Variant Format:
use constant SNV_REF_FWDSEQ => 0;
use constant SNV_ALT_FWDSEQ => 1;
use constant SNV_REF_MODSTR => 2;
use constant SNV_ALT_MODSTR => 3;
use constant SNV_CRDS_STR   => 4;
use constant SNV_ALTS_STR   => 5;
use constant SNV_REFS_STR   => 6;

## BSP Format:
use constant BSP_ILMNID       => 0;
use constant BSP_SEQ          => 1;
use constant BSP_QUAL         => 2;
use constant BSP_FLAG         => 3;
use constant BSP_CHR          => 4;
use constant BSP_POS          => 5;
use constant BSP_STRAND       => 6;
use constant BSP_INS_LEN      => 7;
use constant BSP_REF_SEQ      => 8;
use constant BSP_NUM_MISMATCH => 9;
use constant BSP_MIS_STR      => 10;

## Manifest Map File:
use constant MM_KEY                   => 0;
use constant MM_BSC_KEY               => 1;
use constant MM_PRB_M                 => 2;
use constant MM_PRB_U                 => 3;
use constant MM_PRB_D                 => 4;
use constant MM_ILMNID                => 5;
use constant MM_NAME                  => 6;
use constant MM_ADDRESSA_ID           => 7;
use constant MM_ALLELEA_PROBESEQ      => 8;
use constant MM_ADDRESSB_ID           => 9;
use constant MM_ALLELEB_PROBESEQ      => 10;
use constant MM_INFINIUM_DESIGN_TYPE  => 11;
use constant MM_NEXT_BASE             => 12;
use constant MM_COLOR_CHANNEL         => 13;
use constant MM_FORWARD_SEQUENCE      => 14;
use constant MM_GENOME_BUILD          => 15;
use constant MM_CHR                   => 16;
use constant MM_MAPINFO               => 17;
use constant MM_SOURCESEQ             => 18;
use constant MM_STRAND                => 19;

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

## Methylation Manifest General Format:
use constant MAN_ILMNID               => 0;
use constant MAN_NAME                 => 1;
use constant MAN_ADDRESSA_ID          => 2;
use constant MAN_ALLELEA_PROBESEQ     => 3;
use constant MAN_ADDRESSB_ID          => 4;
use constant MAN_ALLELEB_PROBESEQ     => 5;
use constant MAN_INFINIUM_DESIGN_TYPE => 6;
use constant MAN_NEXT_BASE            => 7;
use constant MAN_COLOR_CHANNEL        => 8;
use constant MAN_FORWARD_SEQUENCE     => 9;
use constant MAN_GENOME_BUILD         => 10;
use constant MAN_CHR                  => 11;
use constant MAN_MAPINFO              => 12;
use constant MAN_SOURCESEQ            => 13;
use constant MAN_STRAND               => 14;

## Methylation Manifest 27k:
use constant MAN27_ILMNID               => 0;
use constant MAN27_NAME                 => 1;
use constant MAN27_ILMNSTRAND           => 2;
use constant MAN27_ADDRESSA_ID          => 3;
use constant MAN27_ALLELEA_PROBESEQ     => 4;
use constant MAN27_ADDRESSB_ID          => 5;
use constant MAN27_ALLELEB_PROBESEQ     => 6;
use constant MAN27_GENOMEBUILD          => 7;
use constant MAN27_CHR                  => 8;
use constant MAN27_MAPINFO              => 9;
use constant MAN27_PLOIDY               => 10;
use constant MAN27_SPECIES              => 11;
use constant MAN27_SOURCE               => 12;
use constant MAN27_SOURCEVERSION        => 13;
use constant MAN27_SOURCESTRAND         => 14;
use constant MAN27_SOURCESEQ            => 15;
use constant MAN27_TOPGENOMICSEQ        => 16;
use constant MAN27_NEXT_BASE            => 17;
use constant MAN27_COLOR_CHANNEL        => 18;
use constant MAN27_TSS_COORDINATE       => 19;
use constant MAN27_GENE_STRAND          => 20;
use constant MAN27_GENE_ID              => 21;
use constant MAN27_SYMBOL               => 22;
use constant MAN27_SYNONYM              => 23;
use constant MAN27_ACCESSION            => 24;
use constant MAN27_GID                  => 25;
use constant MAN27_ANNOTATION           => 26;
use constant MAN27_PRODUCT              => 27;
use constant MAN27_DISTANCE_TO_TSS      => 28;
use constant MAN27_CPG_ISLAND           => 29;
use constant MAN27_CPG_ISLAND_LOCATIONS => 30;
use constant MAN27_MIR_CPG_ISLAND       => 31;
use constant MAN27_MIR_NAMES            => 32;

## Methylation Manifest 450k:
use constant MAN450_ILMNID                      => 0;
use constant MAN450_NAME                        => 1;
use constant MAN450_ADDRESSA_ID                 => 2;
use constant MAN450_ALLELEA_PROBESEQ            => 3;
use constant MAN450_ADDRESSB_ID                 => 4;
use constant MAN450_ALLELEB_PROBESEQ            => 5;
use constant MAN450_INFINIUM_DESIGN_TYPE        => 6;
use constant MAN450_NEXT_BASE                   => 7;
use constant MAN450_COLOR_CHANNEL               => 8;
use constant MAN450_FORWARD_SEQUENCE            => 9;
use constant MAN450_GENOME_BUILD                => 10;
use constant MAN450_CHR                         => 11;
use constant MAN450_MAPINFO                     => 12;
use constant MAN450_SOURCESEQ                   => 13;
use constant MAN450_CHROMOSOME_36               => 14;
use constant MAN450_COORDINATE_36               => 15;
use constant MAN450_STRAND                      => 16;
use constant MAN450_PROBE_SNPS                  => 17;
use constant MAN450_PROBE_SNPS_10               => 18;
use constant MAN450_RANDOM_LOCI                 => 19;
use constant MAN450_METHYL27_LOCI               => 20;
use constant MAN450_UCSC_REFGENE_NAME           => 21;
use constant MAN450_UCSC_REFGENE_ACCESSION      => 22;
use constant MAN450_UCSC_REFGENE_GROUP          => 23;
use constant MAN450_UCSC_CPG_ISLANDS_NAME       => 24;
use constant MAN450_RELATION_TO_UCSC_CPG_ISLAND => 25;
use constant MAN450_PHANTOM                     => 26;
use constant MAN450_DMR                         => 27;
use constant MAN450_ENHANCER                    => 28;
use constant MAN450_HMM_ISLAND                  => 29;
use constant MAN450_REGULATORY_FEATURE_NAME     => 30;
use constant MAN450_REGULATORY_FEATURE_GROUP    => 31;
use constant MAN450_DHS                         => 32;

## New Methyatlion Manifest CORE:
use constant MANCORE_ILMNID               => 0;
use constant MANCORE_NAME                 => 1;
use constant MANCORE_ADDRESSA_ID          => 2;
use constant MANCORE_ALLELEA_PROBESEQ     => 3;
use constant MANCORE_ADDRESSB_ID          => 4;
use constant MANCORE_ALLELEB_PROBESEQ     => 5;
use constant MANCORE_INFINIUM_DESIGN_TYPE => 6;
use constant MANCORE_NEXT_BASE            => 7;
use constant MANCORE_COLOR_CHANNEL        => 8;
use constant MANCORE_FORWARD_SEQUENCE     => 9;
use constant MANCORE_GENOME_BUILD         => 10;
use constant MANCORE_CHR                  => 11;
use constant MANCORE_MAPINFO              => 12;
use constant MANCORE_SOURCESEQ            => 13;
use constant MANCORE_STRAND               => 14;
use constant MANCORE_STRANDTB             => 15;
use constant MANCORE_SCOREMIN             => 16;
use constant MANCORE_SCOREU               => 17;
use constant MANCORE_SCOREM               => 18;
use constant MANCORE_SCORES               => 19;
use constant MANCORE_CPGCNT               => 20;
use constant MANCORE_PRBSEQU              => 21;
use constant MANCORE_PRBSEQM              => 22;
use constant MANCORE_PRBSEQD              => 23;
#use constant 

## Methylation Manifest EPIC:
use constant MANEPIC_ILMNID                                => 0;
use constant MANEPIC_NAME                                  => 1;
use constant MANEPIC_ADDRESSA_ID                           => 2;
use constant MANEPIC_ALLELEA_PROBESEQ                      => 3;
use constant MANEPIC_ADDRESSB_ID                           => 4;
use constant MANEPIC_ALLELEB_PROBESEQ                      => 5;
use constant MANEPIC_INFINIUM_DESIGN_TYPE                  => 6;
use constant MANEPIC_NEXT_BASE                             => 7;
use constant MANEPIC_COLOR_CHANNEL                         => 8;
use constant MANEPIC_FORWARD_SEQUENCE                      => 9;
use constant MANEPIC_GENOME_BUILD                          => 10;
use constant MANEPIC_CHR                                   => 11;
use constant MANEPIC_MAPINFO                               => 12;
use constant MANEPIC_SOURCESEQ                             => 13;
use constant MANEPIC_STRAND                                => 14;
use constant MANEPIC_UCSC_REFGENE_NAME                     => 15;
use constant MANEPIC_UCSC_REFGENE_ACCESSION                => 16;
use constant MANEPIC_UCSC_REFGENE_GROUP                    => 17;
use constant MANEPIC_UCSC_CPG_ISLANDS_NAME                 => 18;
use constant MANEPIC_RELATION_TO_UCSC_CPG_ISLAND           => 19;
use constant MANEPIC_PHANTOM4_ENHANCERS                    => 20;
use constant MANEPIC_PHANTOM5_ENHANCERS                    => 21;
use constant MANEPIC_DMR                                   => 22;
use constant MANEPIC_450K_ENHANCER                         => 23;
use constant MANEPIC_HMM_ISLAND                            => 24;
use constant MANEPIC_REGULATORY_FEATURE_NAME               => 25;
use constant MANEPIC_REGULATORY_FEATURE_GROUP              => 26;
use constant MANEPIC_GENCODEBASICV12_NAME                  => 27;
use constant MANEPIC_GENCODEBASICV12_ACCESSION             => 28;
use constant MANEPIC_GENCODEBASICV12_GROUP                 => 29;
use constant MANEPIC_GENCODECOMPV12_NAME                   => 30;
use constant MANEPIC_GENCODECOMPV12_ACCESSION              => 31;
use constant MANEPIC_GENCODECOMPV12_GROUP                  => 32;
use constant MANEPIC_DNASE_HYPERSENSITIVITY_NAME           => 33;
use constant MANEPIC_DNASE_HYPERSENSITIVITY_EVIDENCE_COUNT => 34;
use constant MANEPIC_OPENCHROMATIN_NAME                    => 35;
use constant MANEPIC_OPENCHROMATIN_EVIDENCE_COUNT          => 36;
use constant MANEPIC_TFBS_NAME                             => 37;
use constant MANEPIC_TFBS_EVIDENCE_COUNT                   => 38;
use constant MANEPIC_METHYL27_LOCI                         => 39;
use constant MANEPIC_METHYL450_LOCI                        => 40;
use constant MANEPIC_CHROMOSOME_36                         => 41;
use constant MANEPIC_COORDINATE_36                         => 42;
use constant MANEPIC_SNP_ID                                => 43;
use constant MANEPIC_SNP_DISTANCE                          => 44;
use constant MANEPIC_SNP_MINORALLELEFREQUENCY              => 45;
use constant MANEPIC_RANDOM_LOCI                           => 46;

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
my ($II, $TOP, $BOT, $UNK, $GRN, $RED, $NEG);
my @TB = (); my @FR = (); my @CO = (); my @III = (); my @MU = ();
my @RG = (); my @AB = (); my @ACTG = (); my @MUD = ();
my ($SL, $TAB, $COM, $COL, $RET, $SEM, $PLUS, $MINUS, $DASH);
my ($CG, $TG, $UM, $MA, $FO);
my ($A, $B, $C, $D, $E, $F, $G, $H, $I, $J, $K, $L, $M, $N,
    $O, $P, $Q, $R, $S, $T, $U, $V, $W, $X, $Y, $Z);
my ($d);
#my ($a, $b, $c, $d, $e, $f, $g, $h, $i, $j, $k, $l, $m, $n,
#    $o, $p, $q, $r, $s, $t, $u, $v, $w, $x, $y, $z);
my %iupacNucs = (); my %cmplNucs = ();
my %iupacStrs = (); my %cmplChar = ();
my %basicNucs = (); my %ACTG = ();
my %bspFR = (); my %bspCO = (); my %iupac = ();
my ($lev5, $forcastLim, $forcastFlag, $forcastCnt, $cpgPadLen);

##
## Basic Program Variables
##
my ($topDir, $statusCnt, $qsub, $clean, $verbose, $mode);
my ($outDir, $logDir, $bedDir, $datDir, $fasDir,
    $bamDir, $ordDir, $desDir, $manDir, $vizDir,
    $inpDir, $pngDir, $binDir);
setParams();
$mode = "parse";
##
## Input Variables and Data Structures
##
my %usrDef = ();
my @srcNames = ();
my @bedFiles = (); my $bedMax; my @bedSrc = ();
my @manFiles = (); my $manMax; my @manSrc = ();
my @datFiles = (); my $datMax; my @datSrc = ();
my @fasFiles = (); my $fasMax; my @fasSrc = ();
my @vcfFiles = (); my $vcfMax;

my $vtLine = "##INFO=<ID=VT,Number=.,Type=String,Description=\"indicates what type of variant the line represents\">";

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Program Main Subroutines:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub loadVcfFile($) {
    my ($file) = @_;

    die("[ERROR]: loadVcfFile: Unable to find vcf file=$file!\n") if (! -e $file);

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    $fileName = trimSuffix($fileName);

    if ($fileSuffix eq ".gz") { open(IN, "gzip -dc $file | ") or die("Failed to open $file stream: $!") }
    else { open(IN, "<$file") or die("Failed to open $file: $!") }
    mssg("[Vcf]: Loading reference variants from vcf file=$file...\n", 3);

    my ($vcfFile, $header, $curChr, $date, $dbSrc, $dbVer, $refVer);
    my $GT="1|1";

    my $rsnUncCnt = 0;
    my $rsnComCnt = 0;
    my $rsnTotCnt = 0;
    while(<IN>) {
        s/^\s+//;
        s/\s+$//;

        my $line = $_;
        my @data = split(/\t/, $line);

        if ($line =~ m/^\#\#/) {
            if ($line =~ m/^##fileDate=(.*)$/) {
                $date = $1;
		} elsif ($line =~ m/^##source=(.*)$/) {
			 $dbSrc = $1;
		    } elsif ($line =~ m/^##dbSNP_BUILD_ID=(.*)$/) {
                $dbVer = $1;
            } elsif ($line =~ m/^##reference=([^.]+)/) {
                $refVer = $1;
            }
            $header .= $line."\n";
            die("[ERROR]: Line already contains 'VT='\n") if ($line =~ m/VT=/);
            next;
        }
        if ($line =~ m/^\#CHROM/) {
            die("[ERROR]: Failed to find dbSrc, dbVer, date or refVer!\n")
                if (! defined($dbSrc) || ! defined($dbVer) || ! defined($date) || ! defined($refVer));
            $header .= "$vtLine\n";
            $header .= $line."\t$refVer\n";
            $vcfFile = "$outDir/".join($L, "common", "$dbSrc.$dbVer", $refVer. $date).".vcf";
            mssg("[vcf]: Writing output to: $vcfFile\n");
            open(VCF, ">$vcfFile") or die("Faield to open $vcfFile for writing: $!");
            print VCF $header;
            next;
	}

        my ($chr, $pos, $rsn, $ref, $alt, $qual, $filter, $info, $format) = @data;
	$chr =~ s/^chr//;

        my $isCommon;
        $isCommon = TRUE if ($info =~ m/COMMON=1/);
        if (! $isCommon) { $rsnUncCnt++; next; }

        my $isSNP = TRUE;
	my @refDat = split(/,/, $ref);
        my @altDat = split(/,/, $alt);
        foreach my $refVal (@refDat) {
            $isSNP = FALSE if (length($refVal) != 1);
        }
        foreach my $altVal (@altDat) {
            $isSNP = FALSE if (length($altVal) != 1);
        }
        $info = "COMMON=1;"; $format = "GT";
        $info .= "VT=INDEL" if (! $isSNP);
        $info .= "VT=SNP" if ($isSNP);
        $rsnComCnt++;

        print VCF join($TAB, $chr, $pos, $rsn, $ref, $alt, $qual, $filter, $info, $format, $GT).$RET;

        mssg("[vcf]: statusCnt=$rsnTotCnt, common=$rsnComCnt, uncommon=$rsnUncCnt, chr=$chr...\n", 3)
            if (defined($statusCnt) && $rsnTotCnt % $statusCnt == 0);
        $rsnTotCnt++;
        last if (defined($vcfMax) && $rsnTotCnt >= $vcfMax);
    }
    close(IN);
    close(VCF);
    mssg("[vcf]: status=DONE; total=$rsnTotCnt, common=$rsnComCnt, uncommon=$rsnUncCnt\n",2);
    mssg("\n",2);
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Function Calls:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub loadVcfFiles($) {
    my ($fileRefs) = @_;

    foreach my $file (@{$fileRefs}) {
	loadVcfFile($file);
    }
}

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Main:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

GetOptions("o|out=s"      => \$outDir,
           "vcf=s@"       => \@vcfFiles,

           "s|src=s@"     => \@srcNames,

           "u|usr=s"      => \%usrDef,

           "mod=s"        => \$mode,
           "max=i"        => \$datMax,
           "manMax=i"     => \$manMax,
           "vcfMax=i"     => \$vcfMax,
           "clean=i"      => \$clean,
           "v|verbose=i"  => \$verbose,
    );

my @date = getTime();
mssg("[Date]: $date[0] $date[2]\n\n");

my @curDate = getTime();
print STDERR "[Date]: $curDate[0]\n\n";

if (! defined($outDir) || @vcfFiles == 0) {
    print STDERR "[Usage]: $0 -o|out outputDir\n"
        . " -v|vcf vcfFile(s)\n"
        . "\n"
        . " [ -d|dat designInputFile(s) ]\n"
        . " [ -n|name outputName ]\n"
        . " [ -s|src sourceNames(s)]\n"
        . " [ -m|man manMapFiles(s) ]\n"
        . " [ -tm tmDataFile ]\n"
        . " [ -c|chr targetChromosome ]\n"
        . " [ -b|bed bedFiles(s) ]\n"
        . " [ -u|user userDefined key=value pair(s) ]\n"
        . " [ -mod mode ]\n"
        . " [ -max dataMax ]\n"
        . " [ -manMax manifestMax ]\n"
        . " [ w|writeInput ]\n"
        . " [ v|verbose verbosityLevel ]\n"
        . " [ options ]\n\n";
    print STDERR "[Usage]: Will read data files from STDIN if not specified.\n";
    print STDERR "[Usage]: Missing arguments! Exiting...\n\n";
    exit(1);
}
buildDirs();

loadVcfFiles(\@vcfFiles);


##
## Done!
##
@date = getTime();
runTimes();
mssg("[Done]: $date[0], $date[2], runTime=$runTimes[0]\n\n");
foreach my $runTime (reverse(@runTimes)) {
    mssg("[Time]:\t$runTime\n");
}
exit(0);

## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##
## Basic Subroutines:
##
## ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

sub seqTopBotCG {
    my ($seq, $retName) = @_;

    die("[ERROR]: seqToptBotCG sequence is not defined or length zero!\n") if (! defined($seq) || length($seq) == 0);
    $seq = shearBrac($seq);
    my $len = length($seq);
    my $mid = ($len/2) - 1;
    die("[ERROR]: seqTopBotCG sequence is not even in length: mid=$mid; seq=$seq\n") if ($mid !~ m/^[0-9]+$/);

    my $top=$T; $top=$TOP if ($retName);
    my $bot=$B; $bot=$BOT if ($retName);
    my $unk=$U; $bot=$UNK if ($retName);

    my $rev = revCmp($seq);
    return($unk) if ($seq eq $rev);

    my $bracSeq = addBrac($seq);
    my @nucs = split(//, uc($seq));

    for (my $ii = 1; $ii < $mid; $ii++) {
        my $nucA = $nucs[$mid-$ii];
        my $nucB = $nucs[$mid+$ii+1];
	next if ($nucA eq $nucB);

	my $skip = FALSE;
	$skip = TRUE if ((defined($iupac{$W}{$nucA}) && defined($iupac{$W}{$nucB})) ||
			 (defined($iupac{$S}{$nucA}) && defined($iupac{$S}{$nucB})));

        next if ($skip);
        return($top) if (defined($iupac{$W}{$nucA}) && defined($iupac{$S}{$nucB}));
	return($bot) if (defined($iupac{$S}{$nucA}) && defined($iupac{$W}{$nucB}));
	die("[ERROR]: Should not be here!!! Failed above:\n".
	    "[ERROR]: ! defined(iupac{W/S}{nucA=$nucA})=!defined($iupac{$W}{$nucA}/$iupac{$S}{$nucA})'\n".
	    "[ERROR]: ! defined(iupac{W/S}{nucB=$nucB})=!defined($iupac{$W}{$nucB}/$iupac{$S}{$nucB})\n".
	    "[ERROR]: Seq=$seq\n");
    }
    return($S);
    die("[ERROR]: seqBotTopCG: Should not be here!!! Failed above:\n".
	"[ERROR]: fwd=".addBrac($seq)."\n".
	"[ERROR]: rev=".addBrac($rev)."\n");
}

sub setParams {
    $lev5 = 5;

    $verbose   = TRUE;
    $topDir    = abs_path(getcwd);
    $statusCnt = 100000;
    $qsub = "qsub -cwd -pe threaded 16 -l excl=true -N";

    $forcastLim  = 16384;
    $forcastFlag = TRUE;
    $forcastCnt  = 0;
    $cpgPadLen   = 8;

    $A="A"; $B="B"; $C="C"; $D="D"; $E="E"; $F="F"; $G="G"; $H="H"; $I="I"; $J="J"; $K="K"; $L="L"; $M="M"; $N="N";
    $O="O"; $P="P"; $Q="Q"; $R="R"; $S="S"; $T="T"; $U="U"; $V="V"; $W="W"; $X="X"; $Y="Y"; $Z="Z";

    $d="d";
#    $a="a"; $b="b"; $c="c"; $d="d"; $e="e"; $f="f"; $g="g"; $h="h"; $i="i"; $j="j"; $k="k"; $l="l"; $m="m"; $n="n";
#    $o="o"; $p="p"; $q="q"; $r="r"; $s="s"; $t="t"; $u="u"; $v="v"; $w="w"; $x="x"; $y="y"; $z="z";

    $I="I"; $II="II"; $GRN="Grn"; $RED="Red"; $NEG="Neg";
    $L="_"; $SL="/";

    $CG="CG"; $TG="TG"; $UM="UM"; $MA="MA"; $FO="FO";
    @ACTG = ($A, $C, $T, $G); @FR  = ($F, $R); @CO  = ($C, $O); @MU  = ($M, $U); @MUD = ($M, $U, $D);
    @RG = ($RED, $GRN); @AB = ($A, $B); @III = ($I, $II); @TB=($T,$B);
    foreach my $nuc (@ACTG) { $ACTG{$nuc} = $nuc; }

    $cmplChar{$F} = $R; $cmplChar{$R} = $F; ## F/R
    $cmplChar{$C} = $O; $cmplChar{$O} = $C; ## C/O
    $cmplChar{$M} = $U; $cmplChar{$U} = $M; ## M/U
    $cmplChar{$T} = $B; $cmplChar{$B} = $T; ## T/B
    $cmplChar{$D} = $D; ## Degnerate???

    $TOP="TOP"; $BOT="BOT"; $UNK="UNK";
    $TAB="\t"; $COM=","; $COL=":"; $RET="\n"; $SEM=";"; 
    $PLUS="+"; $MINUS="-"; $DASH="-";

    $bspFR{$PLUS}  = $F; $bspFR{$F} = $PLUS;
    $bspFR{$MINUS} = $R; $bspFR{$R} = $MINUS;

    $bspCO{$PLUS}  = $O; $bspCO{$O} = $PLUS;
    $bspCO{$MINUS} = $C; $bspCO{$C} = $MINUS;


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

    foreach my $nuc (@ACTG) { $basicNucs{$nuc} = $nuc; }
    foreach my $deg (sort keys %iupacNucs) {
	my $nucs = "";
	foreach my $nuc (sort @{$iupacNucs{$deg}}) {
	    $nucs .= $nuc;
	    $iupac{$nuc}{$deg} = 1;
	    $iupac{$deg}{$nuc} = 1;
	}
	$iupacStrs{$nucs} = $deg;
    }
}

sub runTimes
{
    my ($curTime, $perTime, $segTime, $totTime);
    $curTime = time;
    $segTime = sprintf("%.2f", ($curTime - $endTime)/60);
    $endTime = $curTime;
    $totTime = $endTime - $startTime;
    #$perTime = sprintf("%.2f", (100*$segTime/$totTime));
    $totTime = sprintf("%.2f", ($totTime/60));
    unshift @runTimes, "$segTime(m)\t$totTime(m)";
    #return(reverse(@runTimes));
}

sub buildDirs {
    my $cleanAll = shift;

    return(0) if (! defined($outDir));

    $outDir =~ s/\/+$//;
    system("rm -rf $outDir") if (defined($outDir) && defined($cleanAll) && 
				 $cleanAll && -e $outDir);

    make_path($outDir) if (! -e $outDir);
    $outDir = abs_path($outDir);
    mssg("[Input]: outDir=$outDir\n\n",2);

#    $logDir = "$outDir/log";
#    make_path($logDir) if (! -e $logDir);
#    system("rm -f $logDir/* &> $outDir/cleanLog.log");

#    $binDir = "$outDir/bin";
#    make_path($binDir) if (! -e $binDir);
#    system("rm -f $bedDir/* &> $logDir/cleanBed.log");

#    $bedDir = "$outDir/bed";
#    make_path($bedDir) if (! -e $bedDir);
#    system("rm -f $bedDir/* &> $logDir/cleanBed.log");

#    $fasDir = "$outDir/fas";
#    make_path($fasDir) if (! -e $fasDir);
#    system("rm -f $fasDir/* &> $logDir/cleanFas.log");

#    $manDir = "$outDir/man";
#    make_path($manDir) if (! -e $manDir);
#    system("rm -f $manDir/* &> $logDir/cleanMan.log");

#    $datDir = "$outDir/dat";
#    make_path($datDir) if (! -e $datDir);
#    system("rm -f $datDir/* &> $logDir/cleanDat.log");

#    $pngDir = "$outDir/png";
#    make_path($pngDir) if (! -e $pngDir);
#    system("rm -f $pngDir/* &> $logDir/cleanPng.log");
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

sub padCGN {
    my ($cpg) = @_;

    my $len = length($cpg);
    $cpg =~ s/^cg//;
    $cpg = "cg". '0' x ($cpgPadLen - $len).$cpg;
    return($cpg);
}

sub hashStr {
    my ($datRef, $sep) = @_;

    $sep = $RET if (! defined($sep));

    my $str = "";
    foreach my $key (sort keys %{$datRef}) {
	$str .= "\t$key => $$datRef{$key}$sep";
    }
    $str =~ s/$sep$//;
    return($str);
}

sub keysStr {
    my ($datRef, $sep) = @_;

    $sep = $COM if (! defined($sep));

    my $str = "";
    foreach my $key (sort keys %{$datRef}) {
	$str .= $key.$sep;
    }
    $str =~ s/$sep$//;
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

sub getBasicStats {
    my ($ref) = @_;

    my @dat = sort {$a <=> $b} @{$ref};
    my $len = scalar(@dat);
    my $sum = 0; my $avg = 0;
    my $mid = int(scalar(@dat)/2);
    my $median = $dat[$mid];
    foreach my $val (@dat) {
        next if (! defined($val) || length($val) == 0 || $val =~ m/^(\s)*$/);
        $sum += $val;
    }
    $avg = $sum / $len if ($len != 0);

    #$sum = 0;
    my $muMAD = 0;
    my $meMAD = 0;
    foreach my $val (@dat) {
        next if (! defined($val) || length($val) == 0 || $val =~ m/^(\s)*$/);
        #$sum += ($val - $avg) * ($val - $avg);
	$muMAD += abs($val-$avg);
	$meMAD += abs($val-$median);
    }
    #my $std = $sum / (scalar(@dat) - 1);
    #$std = sqrt($std);
    $muMAD = $muMAD/$len if ($len != 0);
    $meMAD = $meMAD/$len if ($len != 0);

    return($avg, $median, $muMAD, $meMAD);
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

    my $cmd = "find $dir -name $key";
    #mssg("[find]: $cmd\n");
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
    $str =~ s/\.bsp//;
    $str =~ s/\.sorted//; $str =~ s/\.merged//; $str =~ s/\.intersect//;
    }
    return($str);
}

sub subJobs {
    my ($fileRefs) = @_;

    my $jobTotCnt = 0;
    foreach my $file (@{$fileRefs}) {
        die("[ERROR]: loadGenFile: Unable to find file=$file!\n") if (! -e $file);

        my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
        $fileName = trimSuffix($fileName);
        #my $chr = $fileName;
        #$chr =~ s/^([^_]+).*$/$1/;

        my $binFile = "$binDir/bscConvert$fileName.sh";
        open(OUT, ">$binFile") or die("Failed to open $binFile for writing: $!");
        #print OUT "$bscExe -g $file -o $outDir -v 3\n\n";

        close(OUT);
        system("chmod 777 $binFile");

        #my $cmd = "$qsub bsc$chr $binFile";
        my $cmd = "$qsub $jobTotCnt $binFile";
        mssg("[launching]: $cmd\n");
        system($cmd);

	##
	## Status:
	##
        mssg("[jobs]: statusCnt=$jobTotCnt...\n", 3)
            if (defined($statusCnt) && $jobTotCnt % $statusCnt == 0);
        $jobTotCnt++;
        last if (defined($manMax) && $jobTotCnt >= $manMax);
    }
    mssg("[jobs]: Done launching jobs\n",2);
    mssg("\n",2);
}

##
## Common file formats:
##
sub gzip {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);

    return($file) if ($fileSuffix eq ".gz");
    system("rm -f $file.gz") if (-e "$file.gz");
    die("[bgzip]: Failed bgzip $file") if (system("bgzip $file"));
    system("rm -f $file") if ($clean);
    return("$file.gz");
}

sub faiIdx {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $sortedBed = "$filePath$fileName.sorted.bed";

    if ($fileSuffix eq ".gz") { mssg("[Warning]: bgzipIdx file already compressed: $file!\n"); return(); }
    system("rm -f $file.gz") if (-e "$file.gz");
    system("rm -f $file.gz.tbi") if (-e "$file.gz.tbi");
    system("rm -f $file.gz.fai") if (-e "$file.gz.fai");
    die("[faiIdx]: Failed bgzip $file") if (system("bgzip $file"));
    die("[faiIdx]: Failed samtools faidx $file.gz") if (system("samtools faidx $file.gz"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub bgzipIdx {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $sortedBed = "$filePath$fileName.sorted.bed";

    if ($fileSuffix eq ".gz") { mssg("[Warning]: bgzipIdx file already compressed: $file!\n"); return(); }
    system("rm -f $file.gz") if (-e "$file.gz");
    system("rm -f $file.gz.tbi") if (-e "$file.gz.tbi");
    die("[BgzipIdx]: Failed bgzip $file") if (system("bgzip $file"));
    die("[BgzipIdx]: Failed tabix $file.gz") if (system("tabix $file.gz"));
    system("rm -f $file") if ($clean);
    return($sortedBed);
}

sub sortBgzipIdx {
    my ($file, $clean) = @_;

    my ($fileName, $filePath, $fileSuffix) = fileparse( abs_path($file), qr/\.[^.]*/);
    my $sortedBed = "$filePath$fileName.sorted.bed";

    if ($fileSuffix eq ".gz") { mssg("[Warning]: bgzipIdx file already compressed: $file!\n"); return(); }
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

#    if (@{$a} != @{$b}) {
#	mssg("Different Sizes: " . scalar(@{$a}) . ", " . scalar(@{$b}) . "\n");
#    }

    return(FALSE) if (@{$a} != @{$b});
    my $misCnt = 0;
    for (my $ii = 0; $ii < @{$a}; $ii++) {
#	if ($$a[$ii] ne $$b[$ii]) { mssg("[Diff]: $ii: $$a[$ii] ne $$b[$ii]\n"); }
	$misCnt++ if ($$a[$ii] ne $$b[$ii]);
	#return(FALSE) if ($$a[$ii] ne $$b[$ii]);
    }
    my $per = $misCnt / scalar(@{$a});
#    mssg("[Mis]: misPer=$per\n\n");
    return(FALSE) if ($per > 0.95);

    return(TRUE);
}

##
## Program subroutines
##
sub cmpl {
    my ($seq) = @_;

    die("[ERROR]: Cmpl: Seq=Null\n") if (! defined($seq));
    return("") if (! defined($seq) || length($seq) == 0);

    $seq =~ tr/ACTGRYSWKMBDHV/TGACYRSWMKVHDB/;
    $seq =~ tr/actgryswkmbdhv/tgacyrswmkvhdb/;
    return($seq);
}

sub revCmp {
    my ($seq) = @_;
    return(reverse(cmpl($seq)));
}

sub bscD {
    my ($snvSeq) = @_;

    die("[ERROR]: bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) == 0);
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
    return($bscSeq);
}

sub bscM {
    my ($snvSeq) = @_;

    die("[ERROR]: bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) == 0);
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
    return($bscSeq);
}

sub bscU {
    my ($snvSeq, $refSeq, $prbFR, $prbCO) = @_;

    die("[ERROR]: bscGeneral: snvSeq == 0!") if (! defined($snvSeq) || length($snvSeq) == 0);
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
    return($bscSeq);
}

sub bsc {
    my ($seq, $type) = @_;

    die("[ERROR]: bsc(seq, type): sequence not defined!\n") if (! defined($seq) || length($seq) == 0);
    die("[ERROR]: bsc(seq, type): type not defined or lenggth 1!\n") if (! defined($type) || length($type) != 1);

    if ($type eq $D) {
	return(bscD($seq));
    } elsif ($type eq $M) {
	return(bscM($seq));
    } elsif ($type eq $U) {
	return(bscU($seq));
    } else {
	die("[ERROR]: Unsupported bsc type=$type!\n");
    }

    return($seq);
}

sub addBrac {
    my ($seq) = @_;

    my $idx = int(length($seq)/2);
    my $pre = substr($seq, 0, $idx-1);
    my $mid = substr($seq, $idx-1,2);
    my $suf = substr($seq, $idx+1);
    my $bra = $pre."[".$mid."]".$suf;
    die("[ERROR]: addBrac flanks are not 60!\n".
	"[ERROR]: seq=$seq\n".
	"[ERROR]: pre=$pre (len=".length($pre).")\n".
	"[ERROR]: mid=$mid (len=".length($mid).")\n".
	"[ERROR]: pre=$suf (len=".length($suf).")\n")
	if (length($pre) != 60 || length($suf) != 60);

    return($bra);
}

sub addBracs {
    my ($seq, $swap) = @_;

    $swap = "CG" if (! defined($swap));
    $seq =~ s/($swap)/[$1]/gi;

    return($seq);
}

sub shearBrac {
    my ($seq) = @_;
    $seq =~ s/\[//; $seq =~ s/\]//;
    return($seq);
}

sub doubleArray {
    my ($datRef, $x) = @_;

    die("[ERROR]: doubleArray: datRef not defined!\n") if (! defined($datRef));
    die("[ERROR]: doubleArray: x not defined!\n") if (! defined($x) || $x !~ m/^[0-9]+$/);

    my @ret = ();
    foreach my $seq (@{$datRef}) {
        for (my $ii = 0; $ii < $x; $ii++) {
            push @ret, $seq;
        }
    }
    return(@ret);
}

sub expandSeqHash {
    my ($orgSeq) = @_;

    my @expSeqs = expandSeq($orgSeq);
    my %seqHash = ();
    foreach my $seq (@expSeqs) { $seqHash{$seq} = 1; }
    return(%seqHash);
}

sub expandSeq {
    my ($orgSeq, $retHash) = @_;

    if ($retHash) { mssg("[expandSeq]: Will return hash!\n"); }

    die("[ERROR]: epxandSeq: orgSeq not defined!") if (! defined($orgSeq) || length($orgSeq) == 0);
    my @seqs   = ($orgSeq);
    my @seqDat = split(//, $orgSeq);

    my @degIdx = grep { $seqDat[$_] =~ /[^actgACTG]+/ } (0..$#seqDat);
    return(@seqs) if (! @degIdx);

    my $fc = 1;
    foreach my $dIdx (@degIdx) { $fc *= scalar(@{$iupacNucs{uc($seqDat[$dIdx])}}); }
    if ($fc > $forcastLim) { $forcastFlag = $fc; $forcastCnt++; return(@seqs); mssg("[Forcast]: $fc>$forcastLim; $orgSeq!\n"); }

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
    return(@seqs);
}

sub expandSeqOld {
    my ($orgSeq) = @_;

    die("[ERROR]: epxandSeq: orgSeq not defined!") if (! defined($orgSeq) || length($orgSeq) == 0);


    my @seqs = ($orgSeq);
    my @isDeg = grep { /[^actgACTG]+/ } split(//, $orgSeq); 
    return(@seqs) if (! @isDeg);

    my $fc = 1;
    foreach my $deg (@isDeg) { $fc *= scalar(@{$iupacNucs{$deg}}); }
    mssg("FC=$fc\n");

    my @seqDat = split(//, $orgSeq);
    my @subDat = ();
    my $forcast = 1;
    $forcastFlag = FALSE;
    for (my $ii = 0; $ii < @seqDat; $ii++) {
	my $seqNuc = uc($seqDat[$ii]);
	mssg("[expandSeq]: seqNuc=$seqDat[$ii] <=> $seqNuc\n",$lev5) if (! defined($ACTG{$seqNuc}));
        if (defined($iupacNucs{$seqNuc})) {
            my $x = scalar(@{$iupacNucs{$seqNuc}});
	    $forcast *= $x;
	}
    }
    die("[ERROR]: Forcasts do not match: $fc != $forcast: @isDeg, orgSeq=$orgSeq, S=".scalar(@{$iupacNucs{$isDeg[0]}})."\n") if ($fc != $forcast);

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
    my ($prbKey, $degPrb, $datRef) = @_;

    die("[ERROR]: expandSet: prbKey not defined!") if (! defined($prbKey) || length($prbKey) == 0);
    die("[ERROR]: expandSet: degPrb not defined!") if (! defined($degPrb) || length($degPrb) == 0);

    my @degPrbs = expandSeq($degPrb);
    die("[ERROR]: Failed foracst limit for expansion=$forcastFlag, prbKey=$prbKey, degPrb=$degPrb!\n") if ($forcastFlag != FALSE);

    foreach my $prb (@degPrbs) {
	my $degKey = join($L, $prbKey, getCigarDeg($prb));
	$$datRef{uc($prb)} = [ $prb, $degKey ];
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

sub padRC {
    my ($str, $side, $tLen) = @_;

    $str  = revCmp($str);
    return(pad($str, $side, $tLen));
}

sub pad {
    my ($str, $side, $tLen) = @_;

    $tLen = 110 if (! defined($tLen));

    if (length($str) > $tLen) { mssg("[Warning]: too long (" .length($str) ."): $str!\n", 5); return($str); }
    $str .= " " x ($tLen - length($str)) if (! $side);
    $str = " " x ($tLen - length($str)) . $str if ($side);
    return($str);
}

sub getDegPrb {
    my ($prbs) = @_;

    die("[ERROR]: getDegPrb: prbs not defined!") if (! defined($prbs));

    my @mat = ();
    my $prbLen;
    my $sLen = 0;
    foreach my $prb (@{$prbs}) {
#    foreach my $prb (sort keys %{$prbs}) {
	my @seqDat = split(//, uc($prb));
	push @mat, [ @seqDat ];

	die("[ERROR]: Failed matching lengths! $prb: " . length($prb) . " != $prbLen!") if (defined($prbLen) && $prbLen != length($prb));
	$prbLen = length($prb) if (! defined($prbLen));
	$sLen++;
    }
    mssg("[SLen]: $sLen\n",$lev5);
    my $degPrb = "";
    for (my $ii = 0; $ii < $prbLen; $ii++) {
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

sub getDegPrbHash {
    my ($prbs) = @_;

    die("[ERROR]: getDegPrb: prbs not defined!") if (! defined($prbs));

    my @mat = ();
    my $prbLen;
    my $sLen = 0;
    foreach my $prb (sort keys %{$prbs}) {
        my @seqDat = split(//, uc($prb));
        push @mat, [ @seqDat ];

        die("[ERROR]: Failed matching lengths! $prb: " . length($prb) . " != $prbLen!") if (defined($prbLen) && $prbLen != length($prb));
        $prbLen = length($prb) if (! defined($prbLen));
        $sLen++;
    }
    mssg("[SLen]: $sLen\n",$lev5);
    my $degPrb = "";
    for (my $ii = 0; $ii < $prbLen; $ii++) {
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

sub getDegPrb2 {
    my ($prbA, $prbB) = @_;

    die("[ERROR]: getDegPrb2: prbA/B not defined!") if (! defined($prbA) || ! defined($prbB));
    die("[ERROR]: PrbA/B not equal length or zero: $prbA/$prbB!\n") if (length($prbA) != length($prbB) || length($prbA) == 0);

    my @mat = ();
    my $prbLen = length($prbA);
    my $sLen = 2;
    my @matA = split(//, uc($prbA));
    my @matB = split(//, uc($prbB));

    my $degPrb = "";
    for (my $ii = 0; $ii < $prbLen; $ii++) {

	die("[ERROR]: Unknown nuc=$matA[$ii], not found in basic or iupac!\n" .
	    "[ERROR]: MatA=@matA\n") if (! defined($basicNucs{$matA[$ii]}) && ! defined($iupacNucs{$matA[$ii]}));
	die("[ERROR]: Unknown nuc=$matB[$ii], not found in basic or iupac!\n" .
	    "[ERROR]: MatB=@matB\n") if (! defined($basicNucs{$matB[$ii]}) && ! defined($iupacNucs{$matB[$ii]}));

	my %nucs = ();
	if (defined($basicNucs{$matA[$ii]})) {
	    $nucs{$matA[$ii]} = 1;
	} else {
	    map($nucs{$_}++, @{$iupacNucs{$matA[$ii]}});
	}
	if (defined($basicNucs{$matB[$ii]})) {
	    $nucs{$matB[$ii]} = 1;
	} else {
	    map($nucs{$_}++, @{$iupacNucs{$matB[$ii]}});
	}
        my $nucStr = ""; foreach my $nuc (sort keys %nucs) { $nucStr .= $nuc; }

        if (defined($basicNucs{$nucStr})) { $degPrb .= $basicNucs{$nucStr}; next; }
        die("[ERROR]: Unable to find iupac for $nucStr!") if (! defined($iupacStrs{$nucStr}));
        $degPrb .= $iupacStrs{$nucStr};
    }
    die("[ERROR]: degPrb2: degPrb not equal to original length: $prbLen != ".length($degPrb)."\n") if (length($degPrb) ne $prbLen);
    return($degPrb);
}

sub missMatchStr {
    my ($can, $ref) = @_;

    die("[ERROR]: misMatchStr: can not defined!") if (! defined($can));
    die("[ERROR]: misMatchStr: ref not defined!") if (! defined($ref));

    my @canDat = split(//, uc($can));
    my @refDat = split(//, uc($ref));

    my $misCnt = 0;
    my $matStr = "";
    for (my $ii = 0; $ii < @canDat; $ii++) {
#       if ($canDat[$ii] ne $refDat[$ii]) {
        if (! defined($iupac{$canDat[$ii]}{$refDat[$ii]})) {
            $matStr .= "X";
            $misCnt++;
        } else {
            if (! defined($ACTG{$refDat[$ii]}) || ! defined($ACTG{$canDat[$ii]})) {
                $matStr .= "*";
            } else {
                die("[ERROR]: What is this: $refDat[$ii] ~ $canDat[$ii], pos=$ii\n".
                    "[ERROR]: can=$can\n".
                    "[ERROR]: ref=$ref\n") if ($refDat[$ii] ne $canDat[$ii]);
                $matStr .= " ";
            }
        }
    }
    return($matStr, $misCnt);
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

	if (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/[G]/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[C]/ && uc($snvNext) =~ m/[RSKBDV]/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[Y]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[S]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[M]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[B]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[H]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
	    $cpgCnt++;
	} elsif (uc($snvCurr) =~ m/[V]/ && uc($snvNext) =~ m/[GRSKBDV]/) {
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

    mssg("[collapse]: expCrds=@{$expCrds}, expRefs=@{$expRefs}, expAlts=@{$expAlts}\n", 5);

    my %crds = (); my %refs = (); my %alfs = ();
    for (my $ii = 0; $ii < @{$expCrds}; $ii++) {
	my $crd = $$expCrds[$ii];
	my $ref = $$expRefs[$ii];
	my $alt = $$expAlts[$ii];
	mssg("[collapse]: $ii: crd=$crd, ref=$ref, alt=$alt\n", 5);
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

sub addVariantsOld {
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

sub getProbesOld {
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
