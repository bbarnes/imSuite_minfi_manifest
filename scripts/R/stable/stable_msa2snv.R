
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Script for Screen Probe Analytical Performance
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Options, Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )

# Missing installation:: install.packages("RcppParallel")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Set Run Environment:: RStudio/Command-Line
#                            Source All Functions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par  <- list()
args <- commandArgs(trailingOnly = FALSE)

par$src_path <- NULL
par$run_mode <- args[1]
par$date_str <- Sys.Date() %>% as.character()
par$prgm_dir <- 'stable'
par$prgm_tag <- 'stable_msa2snv'
par$verbose  <- 3
local_paths  <- c( 
  # "/Users/bbarnes/Documents/tools/imSuite/scripts/R",
  "/Users/bbarnes/mylocal/Documents/tools/imSuite/scripts/R",
  NULL
)

if ( par$run_mode == "RStudio" ) {
  for (path in local_paths) if ( dir.exists(path) ) par$src_path <- path
} else {
  par$exe_path <- 
    base::normalizePath( base::substring( args[grep("--file=", args)], 8 ) )
  par$scr_path <- base::dirname( par$exe_path )
  par$src_path <- base::dirname( par$scr_path )
  par$run_mode <- 'Command_Line'
}
stopifnot( length(par$src_path) > 0 )
stopifnot( dir.exists(par$src_path) )

par$fun_path <- file.path( par$src_path, "functions" )
stopifnot( dir.exists(par$fun_path) )

src_ret <- base::lapply( base::list.files( 
  par$fun_path, pattern = "\\.R$", full.names = TRUE ), source )
par <- source_functions( pars = par, rcpp = FALSE, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = FALSE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$version <- 0

# par$run_name <- "EPICv1"
# par$run_name <- "Embarkv1"
# par$run_name <- "FAILv1"
# par$run_name <- "COREv1"
#
# par$run_name <- "MSAv03"
par$run_name <- "MSAv10"

opt <- NULL
opt <- imProbeQC_options( pars = par, args = args, vb = par$verbose )
vb  <- opt$verbose
vt  <- 0
tc  <- 0

p0  <- vb > vt + 0
p1  <- vb > vt + 1
p2  <- vb > vt + 2
p4  <- vb > vt + 4
p8  <- vb > vt + 8

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c( 'run_mode', 
               'src_path', 'scr_path', 'exe_path', 'prgm_dir', 'prgm_tag' )
opt_reqs <- c( 'out_path', 
               'Rscript', 'verbose' )

opt$rcpp <- 0
opt$rcpp <- 2
opt$rcpp <- 3
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = opt$rcpp,
                          vb = opt$verbose, vt=3, tc=0 )

opt <- prgm_dat$opt
par <- prgm_dat$par
opt_tib <- prgm_dat$opt_tib
par_tib <- prgm_dat$par_tib

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Initialize Run Objects
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tt <- timeTracker$new()

pmssg <- glue::glue("[{par$prgm_tag}]:")
pwarn <- glue::glue("{RET}[{par$prgm_tag}]: Warning:")
perrs <- glue::glue("{RET}[{par$prgm_tag}]: ERROR:")

success = TRUE;

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Pre-processing:: Wanding Cuts (MSA.v.1.0) Alpha
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

chop_cuts_tib <- NULL
chop_cuts_rds <- file.path( opt$top_path, "Projects.new/MSA/data/CHOP_CUTS_Wanding/20231018_suggested_masking.rds" )
chop_cuts_tib <- readr::read_rds( file = chop_cuts_rds ) %>% tibble::as_tibble()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         Pre-processing:: Load Trifecta Probes (MSA.v.1.0) Alpha
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

Trif_prb_tib <- NULL
Trif_prb_csv <- file.path( opt$top_path, "Projects.new/MSA/trifecta/GSIBIOINFO-668.cluster/20230117-trifecta_strand_specific_designs-rsID_map.prbID.csv.gz" )
Trif_prb_tib <- readr::read_csv( file = Trif_prb_csv, show_col_types = FALSE ) %>%
  dplyr::rename( Probe_ID = Assay_Design_Id ) %>%
  dplyr::mutate( Probe_ID = paste0(Probe_ID,"1" ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         Pre-processing:: Load Previous ILMN Cuts (MSA.v.1.0) Alpha
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

Ilmn_man_tib <- NULL
Ilmn_man_csv <- file.path( opt$top_path, "Projects.new/MSA/data/MSA-Alpha-CHOP/MSA-48v1-0-Post-PQC_3E-DoNotDistribute-Zhou.body.csv.gz" )
Ilmn_man_tib <- readr::read_csv( file = Ilmn_man_csv, show_col_types = FALSE )

Ilmn_cuts_tib <- NULL
Ilmn_cuts_tib <- Ilmn_man_tib %>% 
  dplyr::select( IlmnID, 'ProbeSelection_04-OCT-2023_280K' ) %>% 
  magrittr::set_names( value = c("Probe_ID","Ilmn_Cut_Status") ) %>%
  dplyr::mutate(
    Chop_Cut =  dplyr::case_when(
      Probe_ID %in% chop_cuts_tib$Probe_ID ~ "Y",
      TRUE ~ "N" ),
    Ilmn_Cut = dplyr::case_when(
      Ilmn_Cut_Status == "Bret-Fail-280"  ~ "Y",
      Ilmn_Cut_Status == "Bret-Pass-280"  ~ "N",
      Ilmn_Cut_Status == "Pass-Must-Have" ~ "M",
      TRUE ~ NA_character_ ),
    Trif_Cut = dplyr::case_when(
      Probe_ID %in% Trif_prb_tib$Probe_ID ~ "3",
      TRUE ~ "1" )
  ) %>% dplyr::select( -Ilmn_Cut_Status )
Ilmn_cuts_tib %>% print_sum( vec=c("Chop_Cut","Ilmn_Cut","Trif_Cut") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#     Pre-processing:: Set Training Data Path and Sample Sheet: Alpha
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

local_dat_path <- file.path( opt$top_path, "scratch/stable_imSesameCpp_MSAv1_Alpha/MSAv10-NCBI-v0" )
train_sdf_path <- file.path( local_dat_path, "sdf" )

ssh_tib <- NULL
ssh_csv <- file.path( local_dat_path, "sample_sheets/Samplesheet_48x1_preDVT_072023-SMG-Verbose.formatted.csv.gz" )
ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE ) %>% 
  dplyr::select( -Sentrix_Path )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#           Pre-processing:: Search for Idats (MSA v.1.0) Alpha
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# [TBD]: Add Alpha+DVT
#
#  - Projects.new/MSA/data/MSA-Alpha-CHOP/AlphaData-07122023_PreDVT_SGrun
#    - Samplesheet_48x1_preDVT_072023-SMG-Verbose.csv.gz  
#  - data/idats/idats_MSA_v1-48x1-04072023_MSAEX01V
#    - Methylation_Sample_Sheet_04072023_corrected.csv.gz 
#  - data/idats/idats_MSA_v1-48x1-08302023_DVT1_Run3
#    - Methylation_Sample_Sheet_DVT_Run3_08302023.csv.gz
#

# A tibble: 192 × 13
idat_dir <- file.path( opt$top_path, "Projects.new/MSA/data/MSA-Alpha-CHOP/AlphaData-07122023_PreDVT_SGrun" )
ssh_tib <- ssh_tib %>% dplyr::inner_join( 
  sesame::searchIDATprefixes( dir.name = idat_dir, recursive = TRUE ) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column( var = "Sentrix_Name") %>% 
    magrittr::set_names( c("Sentrix_Name", "Sentrix_Path") ) %>%
    tibble::as_tibble() %>%
    dplyr::distinct( Sentrix_Name, .keep_all = TRUE ),
  by=c("Sentrix_Name")
)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Pre-processing:: Manifest (MSA10):: Full
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# NOT NEEDED FOR ONLY MSA::
# neg_ctl_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_negative_ctls.rds" )

msa_man_dir <- safe_mkdir( dir = file.path( opt$out_path, "manifests" ) )
msa_man_csv <- file.path( msa_man_dir, "msa_man.csv.gz")
v03_man_csv <- file.path( opt$top_path, "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz" )
v10_man_csv <- file.path( opt$top_path, "data/pre-idats/MSA/MSA-48v1-0-Post-PQC_2B.body.csv.gz" )

# A tibble: 301,805 × 10
msa_man_tib <- NULL
msa_man_tib <- build_core_msa_man( neg_rds = NULL, 
                                   m03_csv = v03_man_csv,
                                   m10_csv = v10_man_csv,
                                   out_csv = msa_man_csv,
                                   vb=vb,vt=vt+3,tc=tc )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Pre-processing:: Manifest (MSA10):: Ctls
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# A tibble: 4,077 × 4
man_sub_ctl_tib <- NULL
man_sub_ctl_tib <- msa_man_tib %>% 
  dplyr::filter( stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
  dplyr::filter( Annotation == "NEGATIVE" ) %>%
  dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative") ) %>%
  dplyr::mutate( col = as.factor(col),
                 U = dplyr::case_when( U==0 ~ NA_real_, TRUE ~ U ),
                 M = dplyr::case_when( M==0 ~ NA_real_, TRUE ~ M ),
                 Address = U,
                 Sequence = "N",
                 Type = "NEGATIVE",
                 Color_Channel = "Red",
                 Name = paste0( "Negative_",dplyr::row_number() )
  ) %>% 
  dplyr::select( Address,Type,Color_Channel,Name )

# A tibble: 8,156 × 15
msa_all_ctls_tib <- NULL
msa_all_ctls_tib <- format_msa_ctls( 
  neg_man_csv = file.path( opt$top_path, "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz" ),
  
  out_dir    = file.path( opt$out_path, "manifests/inputs" ),
  run_tag    = opt$run_name,
  reload     = opt$reload,
  reload_min = 10,
  ret_data   = FALSE,
  parallel   = opt$parallel,
  write_out  = FALSE,
  
  vb=vb,vt=vt+1,tc=tc+1, tt=tt )

msa_neg_ctls_tib <- NULL
msa_neg_ctls_tib <- msa_ctls_to_negs( msa_all_ctls_tib ) %>%
  # dplyr::bind_rows( neg_ctls_tib ) %>%
  dplyr::arrange( Address ) %>%
  dplyr::distinct( Address, .keep_all = TRUE )

# A tibble: 178 × 3
# msa_neg_ctls_tib %>% dplyr::anti_join( man_sub_ctl_tib, by=c("Address") )
# A tibble: 0 × 4
# man_sub_ctl_tib %>% dplyr::anti_join( msa_neg_ctls_tib, by=c("Address") )

# Pre-compute neg-vector::
msa_neg_ctls_vec <- NULL
msa_neg_ctls_vec <- msa_neg_ctls_tib %>% 
  dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") ) %>% 
  dplyr::pull(Address) %>% as.integer() %>% unique() %>% as.vector()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Processing:: Single Sample Test
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

vb=0
vt=0
tc=0

outliers_val <- TRUE

opt$pre_pval <- 0.05

opt$min_beta <- 0.3
opt$max_beta <- 0.7
opt$min_perO <- 0.75
opt$min_perI <- 0.05

workflow_vec <- NULL
workflow_vec <- c( "id" )
workflow_vec <- c( "i" )

outliers_vec <- c( TRUE,FALSE )

# OLD SINGLE STATS::
# read_idat_pair_r  8.99      0.0270  0.220      0.00900   8.71    2023-10-18 01:39:32    NA
# A tibble: 301,803 × 22
# tmp_pair_tib <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Processing:: All MSA Data:: Alpha
#                               Parameters
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$one_at_a_time <- TRUE
opt$one_at_a_time <- FALSE

opt$parallel <- FALSE
opt$parallel <- TRUE

#
# Build Train/Tests Sets::
#
train_sentrix_tib <- NULL
train_sentrix_tib <- ssh_tib %>% 
  dplyr::filter( Sample_Input == 250 ) %>% 
  dplyr::filter( Sample_Group != "FFPE" )

tests_sentrix_tib <- NULL
tests_sentrix_tib <- ssh_tib %>%
  dplyr::anti_join( train_sentrix_tib, by=c("Sentrix_Name") )

# Validation:: 192
# dplyr::bind_rows( train_sentrix_tib, tests_sentrix_tib ) %>% dplyr::distinct( Sentrix_Name )

train_sentrix_list <- NULL
train_sentrix_list <- train_sentrix_tib %>% split( .$Sentrix_Name )
tests_sentrix_list <- NULL
tests_sentrix_list <- tests_sentrix_tib %>% split( .$Sentrix_Name )


# Change the naming below...
# msa_neg_ctls_tib => msa_neg_tib

cur_prefix <- NULL
cur_prefix <- ssh_tib %>% dplyr::filter( Sample_GroupKey == "B" ) %>% head(n=1) %>% dplyr::pull( Sentrix_Path )
cur_sename <- ssh_tib %>% dplyr::filter( Sample_GroupKey == "B" ) %>% head(n=1) %>% dplyr::pull( Sentrix_Name )

cur_sdf <- NULL
cur_sdf <- prefix_to_sdf( prefix = cur_prefix,
                          platform = "EPIC", 
                          manifest = msa_man_tib, 
                          # manifest = Ilmn_man_tib %>% 
                          #   dplyr::rename( Probe_ID=IlmnID, U=AddressA_ID, M=AddressB_ID ) %>%
                          #   dplyr::mutate(
                          #     U=as.integer(U), M=as.integer(M),col=as.factor(col)
                          #   ) %>%
                          #   dplyr::select( Probe_ID,U,M,col, dplyr::everything() ),
                          # controls = msa_neg_ctls_tib, 
                          out_dir  = opt$out_path,
                          run_tag  = paste( opt$run_name,cur_sename, sep="." ),
                          reload   = opt$reload, 
                          reload_min = 10,
                          parallel = opt$parallel, 
                          vb=vb, vt=vt, tc=0, tt = tt )

cur_tib <- NULL
cur_tib <- cur_sdf %>% as.data.frame() %>% tibble::as_tibble()

#
# Full CGN DB:
#    [bbarnes@ussd-prd-rdln03 tools]$ wc -l /illumina-sdcolo-01/scratch/darkmatter/data/cgnDB/canonical/20230412-GRCh38_GRCh37_GRCh36_GRCm10_GRCr00_GRCk01_GRGa02_CFAM31_GRCh13_GRCk03_GRLv04_GRSs05_GRSc06_GRBta1_GRCk07_GREL01.cgn-set-canonical.design-input.tsv
#    303987148 /illumina-sdcolo-01/scratch/darkmatter/data/cgnDB/canonical/20230412-GRCh38_GRCh37_GRCh36_GRCm10_GRCr00_GRCk01_GRGa02_CFAM31_GRCh13_GRCk03_GRLv04_GRSs05_GRSc06_GRBta1_GRCk07_GREL01.cgn-set-canonical.design-input.tsv
#
#   [bbarnes@ussd-prd-rdln03 tools]$ tail /illumina-sdcolo-01/scratch/darkmatter/data/cgnDB/canonical/20230412-GRCh38_GRCh37_GRCh36_GRCm10_GRCr00_GRCk01_GRGa02_CFAM31_GRCh13_GRCk03_GRLv04_GRSs05_GRSc06_GRBta1_GRCk07_GREL01.cgn-set-canonical.design-input.tsv
#   303987147	TCTCAGGACTCCCAAAGTTGGTCAATGCCATTTCATCAAAATCAAGCAGTCTTGAGAACT[CG]GCATCATTGTTCCCAATGTTCCCAGGTTTAGGACCCCATATTTCCATAATAAATGAAGGT	GRCk07	0	0	FALSE
#
#
#


# [TBD]: Standalone executable
# [TBD]: Wrapper Rcpp function (load ref-idats + loop/can-idats)
#
# [TBD]: Move TBD below into idat.h
#  - Basically allow the option to add a map of idats to compare against...
#  - Add this to the test function in idat_functions.cpp
# [TBD]: Direct idat comparison via raw intentisty: (%overlap, delta Intensity, r2 Intensity)
#  - Pre-load all individual idats based on Sample Sheet (Sentrix_Name, Sentrix_Path, Chip, Sample)
#  - Pick max manifest, then
#  - Pick max red/grn and ensure they agree
#

# [TBD]: Binary Files
#  - Manifest
#  - Reference Sample Sheet
#  - Manifest + Reference Sample Sheet???
#  - 
#
# [TBD]: Binary File Generation Functions
#

#
#
# Goal: Exapand the "col" column for SNVs...
#
#
# 1. Add Probe_Type (col+Probe_Type)
#
#  cg => col
# !cg2 => N
# !cg0 => col
#
cur_tib %>% 
  dplyr::mutate( 
    Probe_Type = Probe_ID %>% stringr:str_sub(1,2),
    Probe_Call = dplyr::case_when( 
      Probe_Type == "cg" & col == "2" ~ "2",  # 0
      Probe_Type == "cg" & col != "2" ~ col,  # 1
      
      Probe_Type == "ch" & col == "2" ~ "2",  # 0
      Probe_Type == "ch" & col != "2" ~ col,  # 1
      
      Probe_Type == "rs" & col == "2" ~ "D",  # 0
      Probe_Type == "rs" & col != "2" ~ col,  # 
      
      
      col != "2" & !stirngr::str_strats( Probe_ID,"rs") ~ "S",
      col != "2" & !stirngr::str_strats( Probe_ID,"rs") ~ "S",
      col != "2" & stirngr::str_strats( Probe_ID,"rs") ~ "S",
      col != "2" & stirngr::str_strats( Probe_ID,"rs") ~ "S",
      col != "2" & stirngr::str_strats( Probe_ID,"rs") ~ "S",
      
                           
# First target the potential Inf1/SNP probes
#
snv_sdf_tib <- NULL
snv_sdf_tib <- cur_tib %>% 
  dplyr::filter( col != "2" | stringr::str_starts(Probe_ID, "rs") | stringr::str_starts(Probe_ID, "nv") ) 

#
# Now Call SNV base on col/type...
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Data-processing:: Rcpp
#
# Need Alpha  (Singapore) + DVT
#  - Add JIRA Tickets...
#
# Wanding Order: QCDPB: https://zhou-lab.github.io/sesame/dev/sesame.html#Preprocessing_Function_Code
# [TBD]: Merge Codes: sdf_base_functions
# [TBD]: Update Rcpp code with as many functions as possible
#
# Rcpp Functions::
#
#   Update foreach in Rcpp read_idat_pair_rcpp( parallel = TRUE )
#    - Load all manifests
#    - Load all core samples
#
#   inferInfiniumIChannel_rcpp()
#
#   formatVCF_rcpp()
#    - Needs dbinom() => look up in boost
#    - 
#
#    predict_sample_rcpp()
#    - Sample_Name,Sample_Prefix => read_idat_pair_rcpp()
#
#
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Steps::
#
#  *** Update Developers Version of Sesame *** 
#
#  1. Review idat_pair code order of operations 
#  2. Update idat_pair_functions.cpp test code
#     - MSA Example
#  3. Add formatVCF()
#  4. Add inferInfiniumIChannel_rcpp












# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Data-processing:: Standard Sesame::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# opt$parallel <- FALSE
if ( FALSE ) {
  cur_sdfs <- NULL
  if ( opt$parallel ) {
    cur_sdfs <- mclapply( exp_idat_lst,
                          prefix_to_sdf, 
                          platform = "EPIC", 
                          manifest = exp_man_tib, 
                          controls = exp_ctl_tib, 
                          out_dir  = exp_out_dir,
                          run_tag  = opt$run_name,
                          reload   = opt$reload, 
                          reload_min = 10,
                          parallel = opt$parallel, 
                          vb=vb, vt=vt, tc=0, tt = tt )
  } else {
    cur_sdfs <- lapply( exp_idat_lst,
                        prefix_to_sdf, 
                        platform = "EPIC", 
                        manifest = exp_man_tib, 
                        controls = exp_ctl_tib, 
                        out_dir  = exp_out_dir,
                        run_tag  = opt$run_name,
                        reload   = opt$reload, 
                        reload_min = 10,
                        parallel = opt$parallel, 
                        vb=vb, vt=vt, tc=0, tt = tt )
  }
}











# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Reload Pre-Processed Data
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_data_rds <- NULL
all_data_rds <- file.path( opt$out_path, paste0( opt$run_name,"_all_data.rds" ) )
exp_data_rds <- NULL
exp_data_rds <- file.path( opt$out_path, paste0( opt$run_name,"_exp_data.rds" ) )

negs_sdf_rds <- NULL
negs_sdf_rds <- file.path( opt$out_path, paste0( opt$run_name,"_negs_sdf.rds" ) )
auto_ssh_rds <- NULL
auto_ssh_rds <- file.path( opt$out_path, paste0( opt$run_name, "_auto_ssh.rds") )

all_data_tab <- NULL
exp_data_tab <- NULL
negs_sdf_tab <- NULL
auto_ssh_tab <- NULL

pre_load_data <- FALSE
pre_load_data <- TRUE

if ( pre_load_data && (
  # file.exists( all_data_rds ) &&
  file.exists( exp_data_rds ) &&
  # file.exists( negs_sdf_rds ) &&
  file.exists( auto_ssh_rds ) ) ) {
  
  # all_data_tab <- readr::read_rds( file = all_data_rds )
  
  # A tibble: 17,203,683 × 34
  exp_data_tab <- readr::read_rds( file = exp_data_rds )
  
  # negs_sdf_tab <- readr::read_rds( file = negs_sdf_rds )
  
  # A tibble: 193 × 7
  auto_ssh_tab <- readr::read_rds( file = auto_ssh_rds )
  
} else {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Data-processing:: Titraiton Sample List
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  rm_outliers_vec <- NULL
  rm_outliers_vec <- c( TRUE, FALSE )
  # rm_outliers_vec <- c( TRUE )
  
  for ( rm_outliers in rm_outliers_vec ) {
    outliers_str <- "rm-outliers"
    if ( !rm_outliers ) outliers_str <- "ok-outliers"
    
    # cur_pre_file_tib <- NULL
    # cur_pre_file_tib <- tibble::tibble( 
    #   Pre_SDF_Path = list.files( path = file.path( opt$out_path, "pre_sdf",outliers_str ), pattern = paste0("_",outliers_str,".read_idat_pair_r.rds$"), recursive = TRUE, full.names = TRUE ),
    #   Pre_SSH_Path = list.files( path = file.path( opt$out_path, "pre_sdf",outliers_str ), pattern = paste0(".sample_sheet.csv$"), recursive = TRUE, full.names = TRUE ),
    #   Sentrix_Name = Pre_SDF_Path %>% stringr::str_remove("^.*\\/") %>%
    #     stringr::str_remove("_rm-outliers.read_idat_pair_r.rds$")
    # ) %>% dplyr::select( Sentrix_Name, dplyr::everything() )
    
    # Load File List::
    cur_pre_sdf_tib <- NULL
    cur_pre_sdf_tib <- file_list( path    = file.path( opt$out_path, "pre_sdf",outliers_str ),
                                  prefix  = file.path( opt$out_path, "pre_sdf",outliers_str ),
                                  pattern = paste0("_",outliers_str,".read_idat_pair_r.rds$"),
                                  suffix  = paste0("_",outliers_str,".read_idat_pair_r.rds"),
                                  recursive = TRUE,
                                  vb=vb,vt=vt+3,tc=tc ) %>%
      dplyr::bind_rows( .id = "Sentrix_Naem" ) %>% 
      t() %>% as.data.frame() %>% 
      tibble::rownames_to_column( var = "Sentrix_Name") %>% 
      magrittr::set_names( c("Sentrix_Name", "Pre_SDF_Path") ) %>%
      tibble::as_tibble()

    cur_pre_ssh_tib <- NULL
    cur_pre_ssh_tib <- file_list( path    = file.path( opt$out_path, "pre_sdf",outliers_str ),
                                  prefix  = file.path( opt$out_path, "pre_sdf",outliers_str ),
                                  pattern = paste0(".sample_sheet.csv$"),
                                  suffix  = paste0(".sample_sheet.csv"),
                                  recursive = TRUE,
                                  vb=vb,vt=vt+3,tc=tc ) %>%
      dplyr::bind_rows( .id = "Sentrix_Naem" ) %>% 
      t() %>% as.data.frame() %>% 
      tibble::rownames_to_column( var = "Sentrix_Name") %>% 
      magrittr::set_names( c("Sentrix_Name", "Pre_SSH_Path") ) %>%
      tibble::as_tibble()

    for ( exp_name in names(sel_ssh_list) ) {
      all_beta_tibs <- NULL
      all_poob_tibs <- NULL
      all_negs_tibs <- NULL
      
      # if ( exp_name == "T0_CC_250" ) next
      # exp_name <- "titration_0_100"
      # exp_name <- "uhm_0_50"
      # exp_name <- "T1_UH"
      
      cur_out_path <- NULL
      cur_out_path <- safe_mkdir( dir = file.path(opt$out_path, "analysis",outliers_str, paste0(exp_name,"_",outliers_str) ) )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                    Data-processing:: Sentrix Sample List
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # [TBD]: Group by Sample_Base
      ssh_sel_tibs <- NULL
      ssh_sel_tibs <- sel_ssh_list[[exp_name]] %>% 
        # dplyr::inner_join( cur_pre_file_tib, by=c("Sentrix_Name") ) %>% 
        dplyr::inner_join( cur_pre_sdf_tib, by=c("Sentrix_Name") ) %>%
        dplyr::inner_join( cur_pre_ssh_tib, by=c("Sentrix_Name") ) %>%
        dplyr::group_by( Sample_Base ) %>%
        dplyr::mutate( Rep_Idx = dplyr::row_number(),
                       Rep_Key = paste( Sample_Base,Rep_Idx, sep="_") ) %>%
        dplyr::ungroup()
      
      ssh_sel_list <- NULL
      ssh_sel_list <- ssh_sel_tibs %>% 
        split( .$Sentrix_Name ) # %>% head( n=opt$max_rep )
      
      # Load File List
      
      # sentrix_name <- names(ssh_sel_list)[1]
      # idat_path_tib
      
      # [TBD]: Screen this earlier...
      if ( ssh_sel_list %>% length() < 1 ) next
      
      for ( sentrix_name in names(ssh_sel_list) ) {
        
        pre_sdf_path <- ssh_sel_list[[sentrix_name]]$Pre_SDF_Path
        pre_ssh_path <- ssh_sel_list[[sentrix_name]]$Pre_SSH_Path
        pre_sdf_tib <- readr::read_rds( file = pre_sdf_path )
        pre_ssh_tib <- readr::read_csv( file = pre_ssh_path, show_col_types = FALSE )
        
        prefix      <- ssh_sel_list[[sentrix_name]]$Sentrix_Path
        src_key     <- ssh_sel_list[[sentrix_name]]$Source_Key
        ng_input    <- ssh_sel_list[[sentrix_name]]$Sample_Input
        # sample_base <- ssh_sel_list[[sentrix_name]]$Sample_Base
        sample_base <- ssh_sel_list[[sentrix_name]]$Sample_Name_rep
        if ( ssh_sel_list[[sentrix_name]]$Sample_Base_Key == "E" || ssh_sel_list[[sentrix_name]]$Sample_Base_Key == "Z" )
          sample_base <- ssh_sel_list[[sentrix_name]]$Sample_Name_uhm
        
        if ( p1 ) cat(glue::glue("{pmssg} Current prefix({exp_name}) = '{prefix}'{RET}"))
        if ( FALSE ) {
          idat_tib <- pre_sdf_tib
          
          # Negative Controls that fail:
          idat_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") ) %>% dplyr::filter( mask )
          
          sn <- 0.0001
          neg_beta_den_ggg <- NULL
          neg_beta_den_ggg <- idat_tib %>% 
            dplyr::filter( Probe_Type == "ct" ) %>% 
            dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") ) %>%
            ggplot2::ggplot( aes( x=Beta_1, fill=mask) ) +
            ggplot2::geom_density( alpha=0.2 )
          
          plot_tab <- NULL
          plot_tab <- idat_tib %>% 
            # dplyr::filter( Probe_Type == "ct" ) %>%
            dplyr::filter( Probe_Type == "nn" |
                             Probe_Type == "cg" | 
                             # Probe_Type == "ch" | 
                             # Probe_Type == "rs" |
                             Probe_ID %>% stringr::str_starts("ctl_Negative_") |
                             # Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_I_") |
                             # Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_II_") |
                             Probe_Type == "nn" ) %>%
            dplyr::mutate( 
              # Pval_Plot = as.integer( Pval_0 * 10 ),
              Pval_Plot = dplyr::case_when(
                Pval_0 <= 0.01 ~ "pval <= 01",
                Pval_0 <= 0.05 ~ "pval <= 05",
                Pval_0 <= 0.10 ~ "pval <= 10",
                Pval_0 <= 0.20 ~ "pval <= 20",
                Pval_0 <= 0.30 ~ "pval <= 30",
                Pval_0 >  0.30 ~ "pval >  30",
                
                TRUE ~ NA_character_ ),
              Plot_Name = dplyr::case_when(
                Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_I_") ~ "B1",
                Probe_ID %>% stringr::str_starts("ctl_BS_Conversion_II_") ~ "B2",
                Probe_Type == "ct" ~ Probe_ID %>% stringr::str_remove("ctl_") %>% stringr::str_sub(1,2),
                TRUE ~ Probe_Type ),
              int_sum = UG+UR+MG+MR
            ) %>%
            # dplyr::filter( int_sum < 10000 ) %>%
            # dplyr::filter( int_sum < 2000 ) %>%
            # dplyr::filter( int_sum < 1000 ) %>%
            tidyr::pivot_longer( cols = c(UG,UR,MG,MR), 
                                 names_to = "Color", 
                                 values_to = "Intensity", 
                                 values_drop_na = TRUE ) %>%
            dplyr::filter( Intensity != 0 ) 
          
          prb_ints_den_ggg <- plot_tab %>%
            # ggplot2::ggplot( aes( x=log(Intensity+sn), fill=Color) ) +
            ggplot2::ggplot( aes( x=Intensity, fill=Color ) ) +
            ggplot2::geom_density( alpha=0.2 ) +
            ggplot2::facet_grid( rows = vars(Pval_Plot), # cols = vars(mask),
                                 cols = vars(Plot_Name) )
          prb_ints_den_ggg
          
          plot_tab %>% dplyr::group_by( Plot_Name,Pval_Plot ) %>%
            dplyr::summarise( Count=n(), .groups = "drop" )
          
          # Count Stats::
          idat_sum1 <- idat_tib %>% print_sum( vec = c("Probe_Type"),
                                               vb=vb+1,vt=vt,tc=tc+1, tt=tt )
          idat_sum2 <- idat_tib %>% print_sum( vec = c("Probe_Type","mask"),
                                               vb=vb+1,vt=vt,tc=tc+1, tt=tt )
          idat_sum3 <- idat_tib %>% print_sum( vec = c("Probe_Type","mask_0"),
                                               vb=vb+1,vt=vt,tc=tc+1, tt=tt )
        }
        
        # Pretty Sure we can just call this from the pre_sdf files...
        if ( FALSE ) {
          # How to measure outlier removal::
          # 
          negs_sdf <- NULL
          negs_sdf <- idat_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% 
            dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") )
          readr::write_csv( x = negs_sdf, file = cur_negs_csv )
          
          negs_sdf_tab <- negs_sdf_tab %>% 
            dplyr::bind_rows(
              dplyr::mutate( negs_sdf,
                             Rm_Outliers = rm_outliers,
                             Sample = exp_name ) %>%
                dplyr::select( Rm_Outliers,Sample, dplyr::everything() )
            )
        }

        #
        # TBD:: Return select columns with appropriate data types in Rcpp...
        #   - This will allow the code above to be removed!
        #
        
        idat_tib <- pre_sdf_tib
        if ( is.null(idat_tib) ) {
          stop( glue::glue("{perrs} Failed to load idats: prefix='{prefix}'{RET}") )
          break
        }
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                       Data-processing:: Summary
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        idat_sum <- NULL
        idat_sum <- print_sum( tib = idat_tib, vec = c("mask","Probe_Type","col"),
                               vb=vb,vt=vt,tc=tc, tt=tt )
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Data-processing:: Parse Control Data
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        if ( FALSE ) {
          sdf_ctls_tib <- NULL
          sdf_ctls_tib <- pre_sdf_tib %>% 
            dplyr::filter(  stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
            dplyr::filter( !mask ) %>%
            dplyr::select( Probe_ID:mask ) # %>% as.data.frame()
          
          man_ctl_tib <- NULL
          man_ctl_tib <- msa_man_tib %>% 
            dplyr::filter( stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
            dplyr::filter( Annotation == "NEGATIVE" ) %>%
            dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative") ) %>%
            dplyr::mutate( col = as.factor(col),
                           U = dplyr::case_when( U==0 ~ NA_real_, TRUE ~ U ),
                           M = dplyr::case_when( M==0 ~ NA_real_, TRUE ~ M ),
                           Address = U,
                           Sequence = "N",
                           Type = "NEGATIVE",
                           Color_Channel = "Red",
                           Name = paste0( "Negative_",dplyr::row_number() )
            ) %>% 
            dplyr::select( Address,Type,Color_Channel,Name )
          
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          #                   Data-processing:: Workflow Pipeline
          # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
          
          # Need the Genome Studio Version of Controls...
          ses_sdf <- NULL
          ses_sdf <- sesame::readIDATpair( prefix.path = prefix, 
                                           manifest = msa_man_tib %>% dplyr::filter( !Probe_ID %>% stringr::str_starts("ct") ), 
                                           controls = man_ctl_tib ) # %>% dplyr::mutate( M = NA_real_) )
          
          ses_sdf %>% tibble::as_tibble() %>% dplyr::filter( col=="2")
          
          sesame::controls(ses_sdf) %>% tibble::as_tibble()
        }
        
        work_sdf <- NULL
        work_sdf <- pre_sdf_tib %>%
          dplyr::mutate( col = as.factor(col)
                         # U = dplyr::case_when( U==0 ~ NA_real_, TRUE ~ U ),
                         # M = dplyr::case_when( M==0 ~ NA_real_, TRUE ~ M )
          ) %>%
          # dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") ) %>% 
          dplyr::select( Probe_ID:mask ) %>%
          as.data.frame() %>%
          sesame::SigDF( # platform = "EPIC", 
                         ctl = pre_sdf_tib %>% 
                           dplyr::filter( stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
                           dplyr::mutate( col = as.factor(col) ) %>%
                           dplyr::filter( !mask ) %>% 
                             dplyr::select( Probe_ID:mask ) %>%
                           as.data.frame()
                         
                         # ctl = msa_all_ctls_tib
                         # ctl = ctls_tib
                         # ctl = NULL
          ) %>% dplyr::select( Probe_ID:mask )
        sesame::controls(work_sdf) %>% tibble::as_tibble()
        
        work_sdf <- work_sdf %>% 
          mutate_sdf_simple( 
            steps = "D", 
            negs_min = 1.0, 
            poob_min = 1.0, 
            vb=vb+3,vt=vt,tc=tc )
        work_sdf %>% tibble::as_tibble() %>% dplyr::filter( col=="2")
        
        work_sdf <- work_sdf %>% 
          mutate_sdf_simple( 
            steps = "B", 
            # steps = "b", 
            negs_min = 1.0, 
            poob_min = 1.0 )
        work_sdf %>% tibble::as_tibble() %>% dplyr::filter( col=="2")
        
        rep_key <- NULL
        rep_key <- ssh_sel_list[[sentrix_name]]$Rep_Key
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Data-processing:: Extract Beta Values
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        beta_tib <- NULL
        beta_tib <- mutate_sdf_simple( 
          sdf = work_sdf,
          # steps = "v", # No Masking...
          steps = "V", # With masking...
          negs_min = 1.0, 
          poob_min = 1.0, 
          vb=vb, vt=vt, tc=tc ) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Probe_ID") %>% 
          tibble::as_tibble() %>% 
          magrittr::set_names( value = c("Probe_ID",rep_key) )
        # magrittr::set_names( value = c("Probe_ID","Beta") )
        
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #                   Data-processing:: Extract pooBAH P-Values
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        poob_tib <- NULL
        poob_tib <- mutate_sdf_simple( 
          sdf = work_sdf,
          steps = "O", 
          negs_min = 1.0, 
          poob_min = 1.0, 
          vb=vb, vt=vt, tc=tc ) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Probe_ID") %>% 
          tibble::as_tibble() %>% 
          magrittr::set_names( value = c("Probe_ID",rep_key) )
        # magrittr::set_names( value = c("Probe_ID","Poob") )
        
        poob_mat <- NULL
        poob_mat <- as.matrix( poob_tib %>% tibble::column_to_rownames( var = "Probe_ID" ) )
        
        # [TBD]: Read Sample Sheet, but this should be done in an analysis function
        #   so everything is separated: build -> analyze...
        cur_auto_ssh_tib <- NULL
        # cur_auto_ssh_csv <- file.path( cur_out_path, paste0(sentrix_name,".sample_sheet.csv") )
        # cur_auto_ssh_tib <- readr::read_csv( file = cur_auto_ssh_csv, show_col_types = FALSE ) %>%
        cur_auto_ssh_tib <- cur_pre_ssh_tib %>% dplyr::filter( Sentrix_Name == sentrix_name ) %>%
          clean_tib() %>%
          dplyr::mutate( Remove_Outliers = rm_outliers, 
                         Experiment = exp_name,
                         Pass_Poob_Cnt = which(  poob_mat <= 0.05 ) %>% length(),
                         Fail_Poob_Cnt = which( !poob_mat <= 0.05 ) %>% length(),
                         Pass_Poob_Per = round( 100 * Pass_Poob_Cnt/ (Fail_Poob_Cnt+Pass_Poob_Cnt), 3 )
          ) %>% dplyr::select( Remove_Outliers,Experiment, 
                               dplyr::everything() )
        
        auto_ssh_tab <- auto_ssh_tab %>% 
          dplyr::bind_rows( cur_auto_ssh_tib )
        
        if ( cur_auto_ssh_tib$Pass_Poob_Per < opt$Pass_Pval_Min ) next
        
        #
        # Aggregate Data Sets if they pass...
        # [TBD]: Write all of this out and then select base on Call Rates...
        #
        
        if ( is.null( all_beta_tibs[[sample_base]] ) ) {
          all_beta_tibs[[sample_base]] <- beta_tib
        } else {
          all_beta_tibs[[sample_base]] <- all_beta_tibs[[sample_base]] %>% 
            inner_join( beta_tib, by = c("Probe_ID") )
        }
        
        if ( is.null( all_poob_tibs[[sample_base]] ) ) {
          all_poob_tibs[[sample_base]] <- poob_tib
        } else {
          all_poob_tibs[[sample_base]] <- all_poob_tibs[[sample_base]] %>% 
            inner_join( poob_tib, by = c("Probe_ID") )
        }
      }
      if ( is.null(all_beta_tibs[[sample_base]]) ) next
      # if ( all_beta_tibs[[sample_base]] %>% base::ncol() < 4 ) next
      
      if ( p0 ) cat(glue::glue("{pmssg} Finished Processing Idat Pairs sample='{exp_name}'.{RET2}"))
      
      #
      # LEFT OFF HERE::
      #  - [TBD]: Incorporate: {ng_input,rep_key}, could just pull from the sample sheet...
      #
      #  - [TBD]: Need to implement dual input matrices
      #     - This will allow Titration efforts
      #
      # [TBD]: Incorporate weight score
      # [TBD]: Add sample names
      # [TBD]: Reduce output
      # [TBD]: Add EPICv2
      # [TBD]: Create Comparison Groups Sheet
      #
      # [Done]: Add More Stats...
      #
      # [TBD]: Split all_beta_tibs/all_poob_tibs into titration if titration experiment...
      #
      
      # [TBD]: Add Poob Call Rate to sample sheets above...
      # pval_pas_cnt <- which( as.matrix( all_poob_tibs[[1]] %>% tibble::column_to_rownames( var = "Probe_ID" ) ) <= 0.05 ) %>% length()
      # pval_mis_cnt <- which( as.matrix( all_poob_tibs[[1]] %>% tibble::column_to_rownames( var = "Probe_ID" ) ) >  0.05 ) %>% length()
      # pval_pas_per <- base::round( 100 * pval_pas_cnt / ( pval_pas_cnt+pval_mis_cnt ), 3 )
      
      cur_min_dB  <- opt$min_dB
      cur_cmp_str <- "lte"
      cur_is_abs  <- TRUE
      
      exp_name1 <- exp_name
      
      all_beta_tibs2 <- NULL
      all_poob_tibs2 <- NULL
      all_negs_tibs2 <- NULL
      exp_name2 <- NULL
      
      # [TBD]: Find a better way of spliting this up...
      # Pretty Sure this isn't needed anymore...
      # [TBD]: Remove code below...
      if ( exp_name %>% stringr::str_starts("uhm_") ||
           exp_name %>% stringr::str_starts("T1_UM") ||
           exp_name %>% stringr::str_starts("T1_UH") ||
           exp_name %>% stringr::str_starts("T1_HM") ) {
        ssh_uhm_list <- NULL
        ssh_uhm_list <- ssh_sel_tibs %>% split( .$Sample_Titration )
        
        exp_name1 <- ssh_uhm_list[[1]] %>% dplyr::distinct( Sample_Name_uhm ) %>% dplyr::pull( Sample_Name_uhm )
        exp_name2 <- ssh_uhm_list[[2]] %>% dplyr::distinct( Sample_Name_uhm ) %>% dplyr::pull( Sample_Name_uhm )
        
        all_beta_tibs2 <- all_beta_tibs[[2]]
        all_poob_tibs2 <- all_poob_tibs[[2]]
        # all_negs_tibs2 <- all_negs_tibs[[2]]
        
        if ( exp_name == "uhm_0_100" ||
             exp_name == "T1_UM" ||
             exp_name %>% stringr::str_ends("_UM") ) {
          cur_min_dB  <- 0.5
          cur_cmp_str <- "gt"
          cur_is_abs  <- FALSE
        }
        if ( exp_name == "uhm_0_50" ||
             exp_name == "T1_UH" ||
             exp_name %>% stringr::str_ends("_UH") ) {
          cur_min_dB  <- 0.2
          cur_cmp_str <- "gt"
          cur_is_abs  <- FALSE
        }
        if ( exp_name == "uhm_50_100" ||
             exp_name == "T1_HM" ||
             exp_name %>% stringr::str_ends("_HM") ) {
          cur_min_dB  <- 0.1
          cur_cmp_str <- "gt"
          cur_is_abs  <- FALSE
        }
      }
      
      # [TBD]: Update the usage of the name_strA/B. Currently doesn't really do antyhing...
      cur_dB_tib <- NULL
      cur_dB_tib <- calc_dBs( beta_tibA = all_beta_tibs[[1]],
                              pval_tibA = all_poob_tibs[[1]],
                              name_strA = exp_name1, # names(all_beta_tibs[[1]])[1]
                              
                              beta_tibB = all_beta_tibs2,
                              pval_tibB = all_poob_tibs2,
                              name_strB = exp_name2, # names(all_beta_tibs[[2]])[1]
                              
                              min_dB  = cur_min_dB,
                              cmp_str = cur_cmp_str,
                              is_abs  = cur_is_abs,
                              
                              # [TBD]: Change the output directory to not have replicate,
                              #        use experiment output directory instead...
                              # out_dir    = file.path( opt$out_path, "replicate" ),
                              out_dir    = cur_out_path,
                              run_tag    = exp_name, 
                              reload     = opt$reload,
                              reload_min = 10, 
                              # ret_data   = TRUE,
                              ret_data   = FALSE,
                              parallel   = opt$parallel,
                              
                              vb=vb+2,vt=vt+1,tc=tc+1, tt=tt )
      
      if ( FALSE ) {
        
        cur_dB_tib2 <- cur_dB_tib
        
        cur_dB_tib %>% head() %>% as.data.frame()
        cur_dB_tib %>% dplyr::filter( dB_pas_call ) %>% base::nrow()
        
        cur_dB_tib %>% dplyr::filter( dB_med > 0.1 ) %>% base::nrow() /cur_dB_tib %>% base::nrow()
        # cur_dB_tib2 %>% dplyr::filter( dB_med > 0.1 ) %>% base::nrow() /cur_dB_tib %>% base::nrow()
        
        cur_dB_tib  %>% dplyr::filter( db_pas_per > 70 ) %>% base::nrow()
        cur_dB_tib  %>% dplyr::filter( dB_pas_call ) %>% base::nrow()
        # cur_dB_tib2 %>% dplyr::filter( dB_pas_call ) %>% base::nrow()
        
        den_wS_ggg <- NULL
        den_wS_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes( x=dB_med ) ) +
          ggplot2::geom_density()
        den_wS_ggg
        
      }
      
      #
      # General Test Case Plotting::
      #
      if ( FALSE ) {
        box_pas_wS_ggg <- NULL
        box_pas_wS_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes(x=db_pas_per,y=wS_per, group = db_pas_per) ) + 
          ggplot2::geom_boxplot( varwidth = TRUE )
        box_pas_wS_ggg
        
        pnt_pas_wS_ggg <- NULL
        pnt_pas_wS_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes(x=db_pas_per,y=wS_per) ) + 
          ggplot2::geom_point()
        pnt_pas_wS_ggg
        
        den_wS_ggg <- NULL
        den_wS_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes( x=wS_per ) ) +
          ggplot2::geom_density()
        den_wS_ggg
        
        pnt_wS_avg_ggg <- NULL
        pnt_wS_avg_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes( x=wS_per, y=dB_avg ) ) +
          ggplot2::geom_point()
        pnt_wS_avg_ggg
        
        d2d_wS_avg_ggg <- NULL
        d2d_wS_avg_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes( x=wS_per, y=dB_avg ) ) +
          ggplot2::geom_density2d()
        d2d_wS_avg_ggg
        
        pnt_wS_med_ggg <- NULL
        pnt_wS_med_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
          ggplot2::geom_point()
        pnt_wS_avg_ggg
        
        d2d_wS_med_ggg <- NULL
        d2d_wS_med_ggg <- cur_dB_tib %>% 
          ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
          ggplot2::geom_density2d()
        d2d_wS_med_ggg
        
        d2d_wS_med_ggg2 <- NULL
        d2d_wS_med_ggg2 <- cur_dB_tib %>% 
          ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
          ggplot2::geom_density2d() +
          ggplot2::facet_grid(
            rows = vars(db_pas_per)
          )
        d2d_wS_med_ggg2
      }
      
      # [TBD]: Need to plot EPIC correlation against db_pas_per & wS_per, which is better
      # cur_dat$ret_tib %>% dplyr::filter( is.na(wS_per) )
      
      # Bind Data
      # [TBD]: Add Annotation data...
      all_data_tab <- all_data_tab %>% dplyr::bind_rows( 
        cur_dB_tib %>% dplyr::mutate( 
          Remove_Outliers = rm_outliers, 
          Experiment = exp_name ) %>% 
          dplyr::select( Remove_Outliers,Experiment, dplyr::everything() ) )
      
      if ( opt$single_sam ) break
    }
    if ( opt$single_neg ) break
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               Write Data:
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$write_data_rds <- TRUE
opt$write_data_rds <- FALSE

if ( opt$write_data_rds ) {
  
  if ( FALSE ) {
    readr::write_rds( x = all_data_tab, file = all_data_rds, compress = "gz" )
    if ( p1 ) cat(glue::glue("{pmssg} Wrote All Data RDS = '{all_data_rds}'{RET}"))
    
    # readr::write_rds( x = negs_sdf_tab, file = negs_sdf_rds, compress = "gz" )
    # if ( p1 ) cat(glue::glue("{pmssg} Wrote All Data NEG = '{negs_sdf_rds}'{RET}"))
  }
  
  readr::write_rds( x = auto_ssh_tab, file = auto_ssh_rds, compress = "gz" )
  if ( p1 ) cat(glue::glue("{pmssg} Wrote All Data SSH = '{auto_ssh_rds}'{RET}"))
  
  exp_data_tab <- NULL
  exp_data_tab <- all_data_tab %>%
    tidyr::separate( Experiment, into=c("Set","Sam","Con"), sep="_", 
                     remove = FALSE, convert = TRUE ) %>%
    dplyr::mutate( 
      dB_pas_class = as.integer(db_pas_per),
      Grp = Sam %>% stringr::str_sub(2,2),
      Sam = Sam %>% stringr::str_sub(1,1) ) %>%
    dplyr::select( Remove_Outliers,Experiment,Set,Sam,Grp,Con,
                   dplyr::everything() )
  # auto_ssh_tab %>% dplyr::group_by( Experiment ) %>% dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  readr::write_rds( x = exp_data_tab, file = exp_data_rds, compress = "gz" )
  if ( p1 ) cat(glue::glue("{pmssg} Wrote All Data RDS = '{exp_data_rds}'{RET}"))
}

#
# Need these for plotting...
#
grp_uniq_vec <- NULL
grp_uniq_vec <- exp_data_tab$Grp %>% unique()

sam_uniq_vec <- NULL
sam_uniq_vec <- exp_data_tab$Sam %>% unique()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Summarize and Plot Results::
#                             Current Plots...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_basic_plots <- TRUE
run_basic_plots <- FALSE

if ( run_basic_plots ) {
  
  sdf_paths_vec <- NULL
  sdf_paths_vec <- list.files( path = file.path( opt$out_path, "pre_sdf","rm-outliers" ), 
                           pattern = paste0("_","rm-outliers",".read_idat_pair_r.rds$"), 
                           recursive = TRUE, full.names = TRUE )
  
  probe_map_tib <- NULL
  probe_map_tib <- readr::read_rds( file = sdf_paths_vec[length(sdf_paths_vec)] ) %>% 
    dplyr::mutate( Probe_Idx = dplyr::row_number() ) %>%
    dplyr::select( Probe_ID, Probe_Idx )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                        Negative Control Screening::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  neg1_screen_tib <- NULL
  neg1_screen_tib <- exp_data_tab %>% 
    dplyr::filter( Remove_Outliers ) %>%
    dplyr::filter( (Sam != "UM" & Sam != "UH" & Sam != "HM") ) %>%
    dplyr::group_by( Probe_ID ) %>%
    dplyr::summarise( 
      dB_max_med = median( dB_max, na.rm = TRUE ),
      mP_med_med = median( mP_med, na.rm = TRUE ),
      .groups = "drop" ) %>% 
    dplyr::rename( Probe_Idx = Probe_ID ) %>% 
    dplyr::mutate( Probe_Idx = Probe_Idx %>% as.integer() ) %>%
    dplyr::inner_join( probe_map_tib, by=c("Probe_Idx") ) %>%
    dplyr::select( Probe_ID, dplyr::everything() ) %>%
    dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative") )
  
  # neg1_screen_tib %>% dplyr::filter( dB_max_med > 0 & mP_med_med > 0 )
  # neg1_screen_tib %>% dplyr::filter( mP_med_med == -Inf )
  # neg1_screen_tib %>% dplyr::filter( dB_max_med > 0 | dB_max_med == -Inf )
  
  # Passing Negative Controls: 2325 / 4259 = 0.5459028
  #
  neg0_screen_ggg <- NULL
  neg0_screen_ggg <- neg1_screen_tib %>% 
    dplyr::filter( dB_max_med > 0 ) %>% 
    ggplot2::ggplot( aes(x=dB_max_med) ) + 
    ggplot2::geom_density( alpha=0.2 )
  
  neg1_screen_ggg <- NULL
  neg1_screen_ggg <- neg1_screen_tib %>% 
    dplyr::filter( dB_max_med > 0 & mP_med_med > 0 ) %>%
    ggplot2::ggplot( aes(x=mP_med_med) ) + 
    ggplot2::geom_density( alpha=0.2 )
  
  neg2_screen_ggg <- NULL
  neg2_screen_ggg <- neg1_screen_tib %>% 
    ggplot2::ggplot( aes(x=mP_med_med) ) + 
    ggplot2::geom_density( alpha=0.2 )
  
  #
  # Build Select Directory::
  #
  opt$sel_path <- safe_mkdir( file.path(opt$out_path, "selection") )
  
  neg1_sel_screen_tib <- NULL
  neg1_sel_screen_csv <- file.path( opt$sel_path, "negative_controls.addresses.csv.gz")
  neg1_sel_screen_tib <- neg1_screen_tib %>% 
    dplyr::filter( dB_max_med > 0 & mP_med_med > 0 ) %>%
    dplyr::inner_join( msa_man_tib, by=c("Probe_ID") ) %>%
    dplyr::distinct( U ) %>%
    dplyr::rename( Address = U ) # %>% dplyr::pull( Address )
  readr::write_csv( x = neg1_sel_screen_tib, file = neg1_sel_screen_csv )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                           Replicate Screening::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  rep1_screen_tib <- NULL
  rep1_screen_tib <- exp_data_tab %>% 
    dplyr::filter( Remove_Outliers ) %>%
    dplyr::filter( (Sam != "UM" & Sam != "UH" & Sam != "HM") ) %>%
    dplyr::group_by( Probe_ID ) %>%
    dplyr::summarise( 
      Full_Cnt = n(),
      Pass_Cnt = count( db_pas_per >  90 ),
      Fail_Cnt = count( db_pas_per <= 90 ),
      Pass_Per = base::round( 100 * Pass_Cnt/Full_Cnt, 3 ),
      .groups = "drop" ) %>% 
    dplyr::rename( Probe_Idx = Probe_ID ) %>% 
    dplyr::mutate( Probe_Idx = Probe_Idx %>% as.integer() ) %>%
    dplyr::inner_join( probe_map_tib, by=c("Probe_Idx") ) %>%
    dplyr::select( Probe_ID, dplyr::everything() )
  
  # Better be zero below::
  rep1_qc_cnt1 <- rep1_screen_tib %>% dplyr::filter( Pass_Cnt + Fail_Cnt != Full_Cnt ) %>% base::nrow()
  
  rep1_screen_sum <- NULL
  rep1_screen_sum <- rep1_screen_tib %>% 
    print_sum( vec = c("Pass_Per"), 
               vb=vb+3,vt=vt,tc=tc )

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                           Titration Screening::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  uhm1_screen_tib <- NULL
  uhm1_screen_tib <- exp_data_tab %>% 
    dplyr::filter( Remove_Outliers ) %>%
    dplyr::filter( !(Sam != "UM" & Sam != "UH" & Sam != "HM") ) %>%
    dplyr::group_by( Probe_ID ) %>%
    dplyr::summarise( 
      Full_Cnt = n(),
      Pass_Cnt = count( db_pas_per >  90 ),
      Fail_Cnt = count( db_pas_per <= 90 ),
      Pass_Per = base::round( 100 * Pass_Cnt/Full_Cnt, 3 ),
      .groups = "drop" ) %>% 
    dplyr::rename( Probe_Idx = Probe_ID ) %>% 
    dplyr::mutate( Probe_Idx = Probe_Idx %>% as.integer() ) %>%
    dplyr::inner_join( probe_map_tib, by=c("Probe_Idx") ) %>%
    dplyr::select( Probe_ID, dplyr::everything() )
  
  # Better be zero below::
  uhm1_qc_cnt1 <- rep1_screen_tib %>% dplyr::filter( Pass_Cnt + Fail_Cnt != Full_Cnt ) %>% base::nrow()
  
  uhm1_screen_sum <- NULL
  uhm1_screen_sum <- uhm1_screen_tib %>% 
    print_sum( vec = c("Pass_Per"), 
               vb=vb+3,vt=vt,tc=tc )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                            Combined Screening::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # (177787+91406) / 301819
  all1_screen_tib <- NULL
  all1_screen_tib <- dplyr::full_join( 
    rep1_screen_tib %>% dplyr::distinct( Probe_ID, .keep_all = TRUE ),
    uhm1_screen_tib %>% dplyr::distinct( Probe_ID, .keep_all = TRUE ), 
    by=c("Probe_ID"),
    suffix=c("_rep","_uhm") )
  
  # All: 239397 / 301815 = 0.7931912 [ Pass_Per_rep > 95 & Pass_Per_uhm > 60 ]
  # All: 249185 / 301815 = 0.8256217 [ Pass_Per_rep > 95 & Pass_Per_uhm > 60 ]
  #
  # CpG: 167537 / 283133 = 0.5917254 [ Pass_Per_rep > 95 & Pass_Per_uhm > 99 ]
  # CpG: 239297 / 283133 = 0.8451752 [ Pass_Per_rep > 95 & Pass_Per_uhm > 60 ]
  # CpG: 249075 / 283133 = 0.8797102 [ Pass_Per_rep > 95 & Pass_Per_uhm > 60 ]
  #
  all1_screen_tib %>% 
    dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) ) %>%
    dplyr::filter( Probe_Type == "cg" ) %>%
    dplyr::filter( Pass_Per_rep > 95 & Pass_Per_uhm > 60 ) %>%
    dplyr::group_by( Pass_Per_rep,Pass_Per_uhm ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  all1_sel95_screen_vec <- NULL
  all1_sel95_screen_vec <- all1_screen_tib %>% 
    dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) ) %>%
    # dplyr::filter( Probe_Type == "cg" ) %>%
    dplyr::filter( (Probe_Type == "cg" & Pass_Per_rep > 95 & Pass_Per_uhm > 60 ) |
                     Probe_Type != "cg" ) %>%
    dplyr::pull( Probe_ID ) %>% as.vector()
  
  man1_screen_csv <- file.path( opt$sel_path, "msa-v1.0_screened_manifest.v1.csv.gz" )
  man1_screen_tib <- NULL
  man1_screen_tib <- msa_man_tib %>% dplyr::filter( Probe_ID %in% all1_sel95_screen_vec )
  readr::write_csv( x = man1_screen_tib, file = man1_screen_csv )
  
  # Jira = 782
  

  
  
    
  
  #
  # Basic Summaries::
  #
  exp_data_tab %>%
    dplyr::mutate( dB_pas_class = as.integer(db_pas_per/10) ) %>%
    dplyr::group_by( dB_pas_class ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)

  exp_data_tab %>%
    dplyr::mutate( dB_pas_class = as.integer(db_pas_per/10) ) %>%
    dplyr::group_by( dB_pas_class,Experiment ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  # 13757074 / 15090950
  
  # Try Filtering and viewing pass percent...
  exp_data_tab %>% 
    dplyr::filter( Remove_Outliers ) %>%
    dplyr::filter( dB_pas_per > 75 )

  # Two variables wS and call rate. Color by wS and and split by experiment sam[1]
  # - just group by last character of Sam...
  
  # [TBD]: Look at percent passing scores...
  
  # Grp = {L,M,A,H}
  repL_den_wSper_dB_ggg <- NULL
  repL_den_wSper_dB_ggg <- exp_data_tab %>%
    dplyr::filter( (Sam != "UM" & Sam != "UH" & Sam != "HM") ) %>%
    dplyr::mutate( dB_pas_class = as.integer(db_pas_per),
                   Grp = Sam %>% stringr::str_sub(2,2),
                   Sam = Sam %>% stringr::str_sub(1,1) ) %>%
    dplyr::filter( Grp == "L" ) %>%
    ggplot2::ggplot( aes(x=wS_per,
                         # y=wS_per, 
                         color = dB_pas_call,
                         group = dB_pas_call) ) + 
    ggplot2::geom_density( alpha=0.2 ) +
    # ggplot2::facet_grid( rows = vars(Remove_Outliers,dB_pas_class, Con),
    ggplot2::facet_grid( rows = vars(Remove_Outliers, Con),
                         cols = vars(Grp,Sam) )
  repL_den_wSper_dB_ggg
  
  uhm1_den_wSper_dB_ggg <- NULL
  uhm1_den_wSper_dB_ggg <- exp_data_tab %>%
    dplyr::filter( !(Sam != "UM" & Sam != "UH" & Sam != "HM") ) %>%
    dplyr::mutate( dB_pas_class = as.integer(db_pas_per),
                   Grp = Sam %>% stringr::str_sub(2,2),
                   Sam = Sam %>% stringr::str_sub(1,1) ) %>%
    dplyr::filter( Grp == "L" ) %>%
    ggplot2::ggplot( aes(x=wS_per,
                         # y=wS_per, 
                         color = dB_pas_call,
                         group = dB_pas_call) ) + 
    ggplot2::geom_density( alpha=0.2 ) +
    # ggplot2::facet_grid( rows = vars(Remove_Outliers,dB_pas_class, Con),
    ggplot2::facet_grid( rows = vars(Remove_Outliers, Con),
                         cols = vars(Grp,Sam) )
  uhm1_den_wSper_dB_ggg
  
  
  exp_data_sum <- NULL
  exp_data_sum <- exp_data_tab %>% 
    print_sum( vec = c("Experiment","Remove_Outliers","dB_pas_call"),
               vb=vb+2,vt=vt+1,tc=tc+1, tt=tt )
  
  
  # Investigation::
  # exp_split_ssh_tib %>% dplyr::filter( Sentrix_Name == "207675480016_R12C02" ) %>% dplyr::distinct() %>% as.data.frame()
  
  # [TBD]: Rename this to exp_data_tab to be consistent...
  # [TBD]: Add this to above requriements...
  # [TBD]: Reload/Reboot everyhing and start analysis over...

  exp_split_ssh_tib <- NULL
  exp_split_ssh_tib <- auto_ssh_tab %>%
    tidyr::separate( Experiment, into=c("Set","Sam","Con"), sep="_", 
                     remove = FALSE, convert = TRUE ) 

  # Pass_Poob_Per
  exp_auto_ssh_sum1 <- NULL
  exp_auto_ssh_sum1 <- exp_split_ssh_tib %>%
    dplyr::group_by( Set,
                     Sam,
                     Con,
                     Experiment,
                     Remove_Outliers ) %>%
    dplyr::summarise(
      
      cr_cnt = n(),
      cr_min = min(Pass_Poob_Per),
      cr_avg = mean(Pass_Poob_Per),
      cr_med = median(Pass_Poob_Per),
      cr_sds = sd(Pass_Poob_Per),
      cr_mad = mad(Pass_Poob_Per),
      cr_max = max(Pass_Poob_Per),
      
      .groups = "drop" )
  exp_auto_ssh_sum1 %>% print( n=base::nrow(exp_auto_ssh_sum1) )
  
  exp_auto_ssh_sum1 %>% 
    dplyr::group_by( Experiment,Remove_Outliers ) %>%
    summarise(
      cr_min = min(cr_med, na.rm = TRUE ),
      .groups = "drop" ) %>% 
    dplyr::arrange( cr_min ) %>%
    print( n=base::nrow(exp_auto_ssh_sum1) )
  
  exp_auto_ssh_sum2 <- NULL
  exp_auto_ssh_sum2 <- auto_ssh_tab %>%
    # dplyr::filter( Pass_Poob_Per > 85 ) %>%
    dplyr::group_by( Experiment,Remove_Outliers ) %>%
    dplyr::summarise(
      
      cr_cnt = n(),
      cr_min = min(Pass_Poob_Per),
      cr_avg = mean(Pass_Poob_Per),
      cr_med = median(Pass_Poob_Per),
      cr_sds = sd(Pass_Poob_Per),
      cr_mad = mad(Pass_Poob_Per),
      cr_max = max(Pass_Poob_Per),
      
      .groups = "drop" )
  exp_auto_ssh_sum2 %>% print( n=base::nrow(exp_auto_ssh_sum2) )
  
  all_data_sum <- NULL
  all_data_sum <- all_data_tab %>% 
    print_sum( vec = c("Experiment","Remove_Outliers","dB_pas_call"),
               vb=vb+2,vt=vt+1,tc=tc+1, tt=tt )
  
  #
  # [TBD]: Add extra plotting layer by sample types
  # [TBD]: Write final calls
  #        - Get Actual Names BACK ASAP!!!
  # [TBD]: Reverse look up the filtered tango addresses
  #
  work_sdf %>% tibble::as_tibble()
  all_data_tab %>% 
    dplyr::filter( Experiment=="T1_UM" | Experiment=="T1_UM" | Experiment=="T1_HM" ) %>%
    dplyr::group_by( Remove_Outliers,Experiment,Probe_ID ) %>% 
    dplyr::summarise( Count=n(), 
                      Pass_Count = count( dB_pas_call, na.rm = TRUE ),
                      .groups = "drop" )
  
  # 276797 / 301819
  #
  # A tibble: 276,797 × 1
  all_data_tab %>% 
    dplyr::filter( Experiment=="T1_UM" | Experiment=="T1_UH" | Experiment=="T1_HM" ) %>%
    dplyr::filter( dB_pas_call == TRUE ) %>% dplyr::distinct( Probe_ID )
  
  #
  # Iterative Screening::
  #
  
  # A tibble: 301,819 × 1
  screen_tib <- NULL
  screen_tib <- all_data_tab %>% dplyr::distinct( Probe_ID )
  
  # A tibble: 276,575 × 1 [  Remove_Outliers ]
  # A tibble: 269,697 × 1 [ !Remove_Outliers ]
  screen_tib <- screen_tib %>% 
    dplyr::inner_join( 
      all_data_tab %>% 
        dplyr::filter( Remove_Outliers ) %>%
        dplyr::filter( Experiment=="T1_UM" ) %>%
        dplyr::filter( dB_pas_call == TRUE ) %>% 
        dplyr::distinct( Probe_ID ),
      by=c("Probe_ID") )
  
  if ( FALSE ) {
    # A tibble: 234,335 × 1 [  Remove_Outliers ]
    # A tibble: 229,919 × 1 [ !Remove_Outliers ]
    screen_tib <- screen_tib %>% 
      dplyr::inner_join( 
        all_data_tab %>% 
          dplyr::filter( Remove_Outliers ) %>%
          dplyr::filter( Experiment=="T1_UH" ) %>%
          dplyr::filter( dB_pas_call == TRUE ) %>% 
          dplyr::distinct( Probe_ID ),
        by=c("Probe_ID") )
  }
  
  if ( FALSE ) {
    # A tibble: 201,962 × 1 [  Remove_Outliers ]
    # A tibble: 197,971 × 1 [ !Remove_Outliers ]
    screen_tib <- screen_tib %>% 
      dplyr::inner_join( 
        all_data_tab %>% 
          dplyr::filter( Remove_Outliers ) %>%
          dplyr::filter( Experiment=="T1_HM" ) %>%
          dplyr::filter( dB_pas_call == TRUE ) %>% 
          dplyr::distinct( Probe_ID ),
        by=c("Probe_ID") )
  }
    
  # A tibble: 277,632 × 1 [  Remove_Outliers ]
  # A tibble: 277,116 × 1 [ !Remove_Outliers ]
  screen_tib <- screen_tib %>% 
    dplyr::inner_join( 
      all_data_tab %>% 
        dplyr::filter( Remove_Outliers ) %>%
        dplyr::filter( Experiment=="T1_HL_250" ) %>%
        dplyr::filter( dB_pas_call == TRUE ) %>% 
        dplyr::distinct( Probe_ID ),
      by=c("Probe_ID") )
  
  # A tibble: 283,868 × 1 [  Remove_Outliers ]
  # A tibble: 281,884 × 1 [ !Remove_Outliers ]
  screen_tib <- screen_tib %>% 
    dplyr::inner_join( 
      all_data_tab %>% 
        dplyr::filter( Remove_Outliers ) %>%
        dplyr::filter( Experiment=="T1_JL_250" ) %>%
        dplyr::filter( dB_pas_call == TRUE ) %>% 
        dplyr::distinct( Probe_ID ),
      by=c("Probe_ID") )
  
  # A tibble: 281,779 × 1 [  Remove_Outliers ]
  # A tibble: 279,599 × 1 [ !Remove_Outliers ]
  screen_tib <- screen_tib %>% 
    dplyr::inner_join( 
      all_data_tab %>% 
        dplyr::filter( Remove_Outliers ) %>%
        dplyr::filter( Experiment=="T1_RL_250" ) %>%
        dplyr::filter( dB_pas_call == TRUE ) %>% 
        dplyr::distinct( Probe_ID ),
      by=c("Probe_ID") )
  
  # > all_data_tab$Experiment %>% unique()
  # [1] "T1_HL_250"      "T1_HL_50"       "T1_HM"          "T1_JL_250"      "T1_JL_50"       "T1_KL_250"      "T1_KL_50"       "T1_ML_250"      "T1_ML_50"       "T1_NA11922_250" "T1_NA11922_50" 
  # [12] "T1_NA12752_250" "T1_NA12752_50"  "T1_NA12877_250" "T1_NA12877_50"  "T1_NA12878_250" "T1_NA12878_50"  "T1_NA12882_250" "T1_NA12882_50"  "T1_NA1879_250"  "T1_NA1879_50"   "T1_RL_250"     
  # [23] "T1_RL_50"       "T1_UH"          "T1_UM"      
  
  den_pas_dB_ggg <- NULL
  den_pas_dB_ggg <- exp_data_tab %>% 
    # dplyr::mutate( Group_Key = Experiment %>% stringr::str_sub(1,2) ) %>%
    ggplot2::ggplot( aes(x=db_pas_per,
                         # y=wS_per, 
                         group = db_pas_per) ) + 
    ggplot2::geom_density( alpha=0.2 ) +
    ggplot2::facet_grid( rows = vars(Remove_Outliers), 
                         cols = vars(Experiment) )
  den_pas_dB_ggg
  
  box_pas_wS_ggg <- NULL
  box_pas_wS_ggg <- exp_data_tab %>% 
    ggplot2::ggplot( aes(x=db_pas_per,
                         y=wS_per, 
                         group=db_pas_per) ) + 
    ggplot2::geom_boxplot( varwidth = TRUE ) +
    ggplot2::facet_grid( rows = vars(Remove_Outliers), 
                         cols = vars(Experiment) )
  box_pas_wS_ggg
  
  # auto_ssh_tab 
  #
  # [TBD]: Figure out how to add the size of each group to the plot
  # [TBD]: Try using the weight aes...
  if ( FALSE ) {
    
    
    box_pas_wS_ggg2 <- NULL
    box_pas_wS_ggg2 <- exp_data_tab %>% 
      ggplot2::ggplot( aes(x=db_pas_per,
                           y=wS_per,
                           group=db_pas_per) ) + 
      ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
      ggplot2::facet_grid( rows = vars(Sam,Remove_Outliers), 
                           cols = vars(Grp,Experiment) )
    
    # notch doesn't really do anything...
    # notch = TRUE ) +

    # This geom_text stuff doesn't work...
      # ggplot2::geom_text(
      #   stat = "db_pas_per" # aes(label = after_stat("db_pas_per") )
      # )
    box_pas_wS_ggg2
  }
  
  uhm_box_pas_wS_ggg2 <- NULL
  uhm_box_pas_wS_ggg2 <- exp_split_dat_tib %>% 
    ggplot2::ggplot( aes(x=db_pas_per,y=wS_per, group = db_pas_per) ) + 
    ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
    ggplot2::facet_grid( rows = vars(Remove_Outliers), 
                         cols = vars(Set, Sam) )
  
  uhm_box_pas_wS_ggg2
  
}
























if ( FALSE ) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #.                    Summarize and Plot Results::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  rep_data_sum <- NULL
  rep_data_sum <- rep_data_tib %>% 
    print_sum( vec = c("Sample","Rm_Outliers"), vb=vb+3,vt=vt,tc=tc )
  
  if ( par$run_name == "EPICv1" ||
       par$run_name == "COREv1" ||
       par$run_name == "FAILv1" ) {
    
    # rep_data_tib %>% dplyr::filter( Probe_ID %in% epic_ses_dat$mask )
    
    epic_mask_tib <- NULL
    epic_mask_tib <- tibble::tibble( Probe_ID = epic_ses_dat$mask, Masked = TRUE )
    
    rep_mask_tib <- NULL
    rep_mask_tib <- rep_data_tib %>% 
      dplyr::left_join( epic_mask_tib, by=c("Probe_ID") ) %>% 
      dplyr::mutate(
        Masked = dplyr::case_when(
          is.na(Masked) ~ FALSE,
          TRUE ~ TRUE )
      )
    
    rep_mask_sum <- NULL
    rep_mask_sum <- rep_mask_tib %>% 
      dplyr::group_by( Masked ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    
    # Quick Validation::
    #
    # rep_mask_tib %>% 
    #   dplyr::group_by( Masked ) %>% 
    #   dplyr::summarise( Count=n(), .groups = "drop" )
    
    sn <- 0.000001
    
    # Apply Cutoffs/Log
    rep_den_mask_gg <- NULL
    rep_den_mask_gg <- rep_mask_tib %>% 
      dplyr::filter( fin_scr >= 0.9 ) %>%
      # ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
      ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
      ggplot2::geom_density( alpha = 0.2 )  +
      ggplot2::facet_grid( rows = vars(Masked) )
    
    # Conclusion:: Score does show a difference for Replicate score vs. third
    #  party analysis...
    #
    rep_den_mask_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_den_mask.pdf", sep=".") )
    rep_den_mask_gg <- NULL
    rep_den_mask_gg <- rep_mask_tib %>% 
      dplyr::filter( fin_scr >= 0.9 ) %>%
      ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
      # ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
      ggplot2::geom_density( alpha = 0.2 )  +
      ggplot2::facet_grid( rows = vars(Masked),
                           cols = vars(Sample) )
    ggplot2::ggsave( filename = rep_den_mask_pdf, 
                     device = "pdf", width = 7, height = 7, dpi = 320 )
    
  } else {
    
    # Looks good::
    rep_den_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_density.pdf", sep=".") )
    rep_den_gg <- NULL
    rep_den_gg <- rep_data_tib %>% 
      ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
      ggplot2::geom_density( alpha = 0.2 ) +
      ggplot2::facet_grid( cols = vars(Sample),
                           rows = vars(Rm_Outliers) )
    
    ggplot2::ggsave( filename = rep_den_pdf, 
                     device = "pdf", 
                     width = 7, 
                     height = 7, 
                     dpi = 320 )
    
    #
    # Plot Negative w/w0 Outlier Removal
    #
    negs_sdf_ggg <- negs_sdf_tab %>% 
      ggplot2::ggplot( aes(x=UG, fill=Rm_Outliers) ) + 
      ggplot2::geom_density( alpha=0.2 ) +
      ggplot2::facet_grid( cols = vars(Sample),
                           rows = vars(Rm_Outliers) )
    
    # More or less worthless::
    rep_2den_gg <- NULL
    rep_2den_gg <- rep_data_tib %>% 
      ggplot2::ggplot( aes( x=fin_scr, y=med_dbs ) ) +
      ggplot2::geom_density2d() +
      ggplot2::facet_grid( rows = vars(Sample) )
    
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
