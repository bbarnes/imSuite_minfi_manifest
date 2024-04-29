
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
par$prgm_tag <- 'stable_imSesameCpp_MSAv1_Alpha_Analysis'
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
par$version <- 1

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
#
#         Qucik Extraction of Interesting Probes for Dream Designs...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_drm_dat <- TRUE
run_drm_dat <- FALSE

if ( run_drm_dat ) {
  dmr_cnt_tib <- NULL
  dmr_cnt_tsv <- "/Users/bbarnes/mylocal/Documents/data/manifests/methylation/GenomeStudio/EX100K_GRCh37_1F.body-section.cpg-cnts.tsv"
  dmr_cnt_tib <- readr::read_tsv( file = dmr_cnt_tsv, show_col_types = FALSE ) %>%
    magrittr::set_names( c("Probe_ID","Cpg_Cnt") ) %>%
    dplyr::mutate( 
      Probe_ID = Probe_ID %>% stringr::str_remove("_bsc_[A-Z]$"),
      Cpg_Cnt = Cpg_Cnt %>% as.integer()
    ) %>% 
    tidyr::separate( Probe_ID, into=c("Loci_ID","Srd_Str","Strand_FR"), sep="_" ) %>% 
    dplyr::mutate( Probe_ID = paste(Loci_ID,Srd_Str, sep="_") ) %>%
    tidyr::separate( 
      Srd_Str, into=c("Strand_TB","Strand_CO","Infinium_Design","Rep_Num"),
      sep=c(1,2,3,4), remove = TRUE, convert = TRUE ) %>%
    dplyr::select( Probe_ID,Loci_ID,Strand_TB,Strand_CO,Infinium_Design,Rep_Num,Strand_FR,Cpg_Cnt)
  
  dmr_all_tib <- NULL
  dmr_all_csv <- "/Users/bbarnes/mylocal/Documents/data/manifests/methylation/GenomeStudio/EX100K_GRCh37_1F.body-section.genomic-sorted.csv.gz"
  dmr_all_tib <- readr::read_csv( file = dmr_all_csv, show_col_types = FALSE )
  
  # Join Count Data::
  drm_all_sum1 <- NULL
  drm_all_sum1 <- dmr_all_tib %>% 
    dplyr::distinct( AlleleA_ProbeSeq, .keep_all = TRUE ) %>% 
    dplyr::filter( !IlmnID %in% dmr_cnt_tib$Probe_ID ) %>% 
    dplyr::group_by( CHR,Rep_Num ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  #
  # Build Select Annotation Groups::
  #
  drm_grp_tib <- NULL
  drm_grp_tib <- Ilmn_man_tib %>% 
    dplyr::select( Name,Probe_Type, `Trifecta?`, `Must-Have-Bucket`,TOPMED_MAF ) %>% 
    magrittr::set_names( c("Name","Probe_Type","isTrifecta","Must_Class","TOPMED_MAF") ) %>% 
    dplyr::filter( !is.na(Must_Class) ) %>% 
    dplyr::distinct( Name, .keep_all = TRUE ) %>%
    dplyr::arrange( Name )
  drm_grp_tib %>% print_sum( vec=c("Probe_Type") )
  
  dmr_grp_tib2 <- NULL
  dmr_grp_tib2 <- Ilmn_man_tib %>% 
    dplyr::select( Name,Probe_Type, `Trifecta?`, `Must-Have-Bucket`,TOPMED_MAF ) %>% 
    magrittr::set_names( c("Name","Probe_Type","isTrifecta","Must_Class","TOPMED_MAF") ) %>% 
    # dplyr::filter( !is.na(Must_Class) ) %>% 
    dplyr::distinct( Name, .keep_all = TRUE ) %>%
    dplyr::arrange( Name )
  
  # [Full]: A tibble: 59,309 × 26
  dmr_sel_tib0 <- NULL
  dmr_sel_tib0 <- dmr_all_tib %>% 
    dplyr::full_join( drm_grp_tib, by=c("Name","Probe_Type") )
  
  # [Inner:Full]: A tibble: 3,332 × 23
  dmr_sel_tib1 <- NULL
  dmr_sel_tib1 <- dmr_all_tib %>% 
    dplyr::inner_join( drm_grp_tib, by=c("Name","Probe_Type") )
  dmr_sel_tib1 %>% print_sum( vec=c("Probe_Type","isTrifecta","Must_Class") )
  dmr_sel_tib1 %>% print_sum( vec=c("CHR") ) %>% print(n=1000)
  dmr_sel_tib1 %>% print_sum( vec=c("Infinium_Design") ) %>% print(n=1000)
  dmr_sel_tib1 %>% dplyr::filter( !is.na(TOPMED_MAF) ) %>% as.data.frame()
  
  # [Inner:Name]: A tibble: 3,069 × 26
  dmr_sel_tib2 <- NULL
  dmr_sel_tib2 <- dmr_all_tib %>% 
    dplyr::distinct( Name, .keep_all = TRUE ) %>% 
    dplyr::inner_join( drm_grp_tib, by=c("Name","Probe_Type") )
  dmr_sel_tib2 %>% print_sum( vec=c("Probe_Type","isTrifecta","Must_Class") )
  dmr_sel_tib2 %>% print_sum( vec=c("CHR") ) %>% print(n=1000)
  dmr_sel_tib2 %>% print_sum( vec=c("Infinium_Design") ) %>% print(n=1000)
  dmr_sel_tib2 %>% dplyr::filter( !is.na(TOPMED_MAF) ) %>% as.data.frame()
  
  # [Inner:PrbSeqA]: A tibble: 3,135 × 26
  dmr_sel_tib3 <- NULL
  dmr_sel_tib3 <- dmr_all_tib %>% 
    dplyr::distinct( AlleleA_ProbeSeq, .keep_all = TRUE ) %>% 
    dplyr::inner_join( drm_grp_tib, by=c("Name","Probe_Type") )
  dmr_sel_tib3 %>% print_sum( vec=c("Probe_Type","isTrifecta","Must_Class") )
  dmr_sel_tib3 %>% print_sum( vec=c("CHR") ) %>% print(n=1000)
  dmr_sel_tib3 %>% print_sum( vec=c("Infinium_Design") ) %>% print(n=1000)
  dmr_sel_tib3 %>% dplyr::filter( !is.na(TOPMED_MAF) ) %>% as.data.frame()

  # [Inner:PrbSeqA]: A tibble: 3,135 × 26
  dmr_sel_tib4 <- NULL
  dmr_sel_tib4 <- dmr_all_tib %>% 
    dplyr::distinct( AlleleA_ProbeSeq, .keep_all = TRUE ) %>% 
    dplyr::inner_join( dmr_grp_tib2, by=c("Name","Probe_Type") )
  dmr_sel_tib4 %>% print_sum( vec=c("Probe_Type","isTrifecta","Must_Class") )
  dmr_sel_tib4 %>% print_sum( vec=c("CHR") ) %>% print(n=1000)
  dmr_sel_tib4 %>% print_sum( vec=c("Infinium_Design") ) %>% print(n=1000)
  dmr_sel_tib4 %>% dplyr::filter( !is.na(TOPMED_MAF) ) %>% as.data.frame()
  
  dmr_sel_tib5 <- NULL
  dmr_sel_tib5 <- dplyr::bind_rows(
    dmr_sel_tib3,
    head( dmr_sel_tib4, n=430 )
  ) %>%
    dplyr::distinct( AlleleA_ProbeSeq, .keep_all = TRUE )
  
  #
  # Intersecting with Counts::
  #
  
  # dmr_sel_tib1 %>% dplyr::filter( IlmnID %in% dmr_cnt_tib$Probe_ID )
  # dmr_sel_tib1 %>% dplyr::filter( Name %in% dmr_cnt_tib$Loci_ID )
  # dmr_sel_tib1 %>% dplyr::left_join( Name %in% dmr_cnt_tib$Loci_ID )
  
  #
  # Load Dream Order file and subset it...
  #
  drm_ord_cols <- NULL
  drm_ord_cols <- readr::cols(
    Assay_Design_Id        = readr::col_character(),
    AlleleA_Probe_Id       = readr::col_character(),
    AlleleA_Probe_Sequence = readr::col_character(),
    AlleleB_Probe_Id       = readr::col_character(),
    AlleleB_Probe_Sequence = readr::col_character(),
    Normalization_Bin      = readr::col_character()
  )
  
  drm_ord_all_tib <- NULL
  drm_ord_all_csv <- file.path( opt$top_path, "scratch/dream-designs/Homo_Sapiens.order.drm.csv" )
  drm_ord_all_tib <- readr::read_csv( file = drm_ord_all_csv, 
                                      col_names = names(drm_ord_cols$cols), 
                                      col_types = drm_ord_cols )
  
  drm_ord_fin_csv <- file.path( opt$top_path, "scratch/dream-designs/Homo_Sapiens.order.drm.selected.csv.gz" )
  drm_ord_fin_tib <- NULL
  drm_ord_fin_tib <- drm_ord_all_tib %>% 
    dplyr::mutate( Loci_ID = Assay_Design_Id %>% stringr::str_remove("_.*$") ) %>% 
    # dplyr::filter( Loci_ID %in% dmr_sel_tib1$Name ) %>%
    dplyr::filter( Loci_ID %in% dmr_sel_tib5$Name ) %>%
    dplyr::distinct( AlleleA_Probe_Sequence, .keep_all = TRUE )
  readr::write_csv( x = drm_ord_fin_tib %>% dplyr::select(-Loci_ID), 
                    file = drm_ord_fin_csv )
  
  drm_loc_fin_cnt = drm_ord_fin_tib %>% dplyr::distinct( Loci_ID ) %>% base::nrow()
  drm_cgn_fin_cnt = drm_ord_fin_tib %>% dplyr::distinct( Assay_Design_Id ) %>% base::nrow()
  drm_seq_fin_cnt = drm_ord_fin_tib %>% dplyr::distinct( AlleleA_Probe_Sequence ) %>% base::nrow()
  
  ( drm_cgn_fin_cnt / 2 ) + drm_cgn_fin_cnt
  
  drm_ord_fin_tib %>% dplyr::select( -Loci_ID )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#     Pre-processing:: Set Training Data Path and Sample Sheet: Alpha
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

local_dat_path <- file.path( opt$top_path, "scratch/stable_imSesameCpp_MSAv1_Alpha/MSAv10-NCBI-v0" )
train_sdf_path <- file.path( local_dat_path, "sdf" )

ssh_tib <- NULL
ssh_csv <- file.path( local_dat_path, "sample_sheets/Samplesheet_48x1_preDVT_072023-SMG-Verbose.formatted.csv.gz" )
ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Processing:: Load Top Training Data: Alpha
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# work_step_str <- paste("NCD","O","PB", sep="-" )
# train_sel_vec <- NULL
# train_sel_vec <- c( 
#   paste0("rm-outliers/pval-1/mask_1/",work_step_str), 
#   paste0("rm-outliers/pval-0.5/mask_1/",work_step_str),
#   NULL
# )
# dbs_file_tib <- NULL
# dbs_file_tib <- expand_grid( 
#   train_sdf_path, train_sel_vec ) %>% 
#   tidyr::unite( Path, dplyr::everything(), sep = "/", remove = FALSE )

outliers_vec <- NULL
outliers_vec <- c( 
  "rm-outliers",
  NULL
)
min_pval_vec <- NULL
min_pval_vec <- c( 
  "pval-1", 
  "pval-0.5",
  NULL
)
mask_vec <- NULL
mask_vec <- c(
  "mask_1",
  NULL
)
work_step_mat <- NULL
work_step_mat <- c( 
  paste("NCD","O","PB", sep="-" ),
  paste("ND", "O","PB", sep="-" ),
  NULL
)

dbs_file_tib <- NULL
dbs_file_tib <- expand_grid( 
  train_sdf_path, outliers_vec,min_pval_vec,mask_vec,work_step_mat ) %>% 
  tidyr::unite( Path, dplyr::everything(), sep = "/", remove = FALSE )

dbs_file_len <- dbs_file_tib %>% base::nrow()
if ( p1 ) cat(glue::glue("{pmssg} Experiment Counts: = '{dbs_file_len}'{RET}"))

all_dbs_tib <- NULL
all_dbs_sum <- NULL
all_sum_tab <- NULL
for ( ii in c(1: base::nrow(dbs_file_tib)) ) {
  if ( p1 ) cat(glue::glue("{pmssg} Current Path[{ii}] = '{dbs_file_tib$Path[ii]}'{RET}"))
  
  cur_dbs_csvs <- NULL
  cur_dbs_csvs <- file_list( path = dbs_file_tib$Path[ii],
                             # subs = c("rm-outliers","ok-outliers"),
                             # prefix = NULL,
                             suffix = ".calc_dBs.csv.gz",
                             pattern = ".calc_dBs.csv.gz$",
                             recursive = TRUE )
  
  params_tib <- NULL
  params_tib <- tibble::tibble(
    Outliers = dplyr:: case_when(
      dbs_file_tib$outliers_vec[ii] == "rm-outliers" ~ as.integer(0),
      dbs_file_tib$outliers_vec[ii] == "ok-outliers" ~ as.integer(1),
      TRUE ~ NA_integer_ ),
    Min_Pval = as.integer( dbs_file_tib$min_pval_vec[ii] %>% 
                             stringr::str_remove("^pval-") %>% 
                             as.double() * 100 ),
    Mask_Idx = as.integer( dbs_file_tib$mask_vec[ii] %>% 
                             stringr::str_remove("^mask_") %>% 
                             as.double() ),
    Work_Str = dbs_file_tib$work_step_mat[ii] )
  
  
  # [TBD]: Reduce Boolean Fields [DONE Below...]
  # [Done]: Us mcapply to multi-thread... [Testing]...
  #
  # [TBD]: In Previous (Alpha.R) script; clean up the output...
  cur_dbs_tib <- NULL
  cur_dbs_tib <- cur_dbs_csvs %>% # head() %>%
    mclapply( function(x) {
      readr::read_csv( file = x, show_col_types=FALSE ) %>%
        dplyr::mutate( Set = Set %>% as.character() %>% stringr::str_sub(1,1) ) %>%
        dplyr::filter( Probe_Type != "ct" ) %>%
        dplyr::filter( !(Set == "T" & Probe_Type == "rs") ) %>%
        dplyr::filter( !(Set == "T" & Probe_Type == "ch") ) %>%
        dplyr::mutate( 
          Work_Str = dbs_file_tib$work_step_mat[ii],
          Grp = dplyr::case_when( 
            Set != "T" ~ Set, 
            TRUE ~ as.character(Sam) )
        ) %>% dplyr::left_join( Ilmn_cuts_tib, by=c("Probe_ID") )
    }) %>% dplyr::bind_rows() %>% clean_tib() %>% 
    dplyr::select( Work_Str,Grp, dplyr::everything() )
  
  # [TBD]: Remove/Reduce Extra Columns::
  all_dbs_tib <- all_dbs_tib %>% dplyr::bind_rows( cur_dbs_tib )
  
}

if ( FALSE ) {
  sdf_path <- safe_mkdir( file.path( opt$out_path, "sdf") )
  sdf_rds  <- file.path( sdf_path, "all_dBs_sdf.rds" )
  readr::write_rds( x = all_dbs_tib, file = sdf_rds, compress = "gz" )
}
  
#
# Parse Selected Replicate/Titration Samples::
# Add MSA_280k Flag
#
if ( FALSE ) {
  sel_dbs_tib <- NULL
  sel_dbs_tib <- all_dbs_tib %>% 
    dplyr::filter( (Min_Pval == 0.5 & Set != "T") | (Min_Pval == 1 & Set == "T") )
  
  # [Sanity Check]: Validate all Probes Found in Previous Screening::
  chop_mis_cnt <- chop_cuts_tib %>% 
    dplyr::filter( !Probe_ID %in% sel_dbs_tib$Probe_ID ) %>% base::nrow()
  Ilmn_mis_cnt <- Ilmn_cuts_tib %>% 
    dplyr::filter( !Probe_ID %in% sel_dbs_tib$Probe_ID ) %>% base::nrow()
  
  if ( chop_mis_cnt != 0 || Ilmn_mis_cnt != 0 ) {
    stop(glue::glue("{pmssg} Failed chop={chop_mis_cnt} or Ilmn={Ilmn_mis_cnt} Counts!{RET}"))
  } else {
    cat(glue::glue("{pmssg} Passed chop={chop_mis_cnt} or Ilmn={Ilmn_mis_cnt} Counts!{RET}"))
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Processing:: Screen Top Training Data: Alpha
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# [JUNK]: 'ppp' Chop_Pass Beta Value (N/Y+N)
#
# [NOTE]: 'dBp' Delta Beta (dB) Passing Percent (Call Rate)
# [NOTE]: 'dBm' Delta Beta (dB) Median
# [NOTE]: 'mPm' Min Pval (mP) Median
# [NOTE]: 'wSp' Weight Score (wS) Percent

# [NOTE]: C = NA[0-9]+ CEPH
# [NOTE]: L = Cancer Cell Line (Hela,..)
# [NOTE]: T = Titration

# Core Metrics::
#
# dBp_Avg = mean(db_pas_per, na.rm=TRUE),
# dBm_Avg = mean(dB_med, na.rm=TRUE),
# mPm_Avg = mean(mP_med, na.rm=TRUE),
# wSp_Avg = mean(wS_per, na.rm=TRUE),

# across(starts_with('a') & where(is.numeric) ) ) %>% 
#  replace_na(list(ab = 999, ac = 999)

sel_dbs_tib <- NULL
sel_dbs_tib <- all_dbs_tib %>% 
  dplyr::select( -Remove_Outliers, -Mask_Index, -Con ) %>%
  dplyr::mutate(
    Trif_Cut = Trif_Cut %>% as.character(),
    Chop_Cut = dplyr::case_when(
      is.na(Chop_Cut) ~ "U", TRUE ~ Chop_Cut ) %>%
      as.factor(),
    Ilmn_Cut = dplyr::case_when(
      is.na(Ilmn_Cut) ~ "U", TRUE ~ Ilmn_Cut ) %>%
      as.factor(),
    Trif_Cut = dplyr::case_when(
      is.na(Trif_Cut) ~ "U", TRUE ~ Trif_Cut ) %>%
      as.factor(),
    
    Work_Str = Work_Str %>% as.factor(),
    Grp = Grp %>% as.factor(),
    Set = Set %>% as.factor(),
    
    Min_Pval = Min_Pval %>% as.factor(),
    Probe_Type = Probe_Type %>% as.factor(),
    Chop_Cut = Chop_Cut %>% as.factor(),
    Ilmn_Cut = Ilmn_Cut %>% as.factor(),
    Trif_Cut = Trif_Cut %>% as.factor(),
  )

# sel_dbs_tib %>% dplyr::distinct( Chop_Cut,Ilmn_Cut,Trif_Cut )

uhm_dbs_sum <- NULL
uhm_dbs_sum <- sel_dbs_tib %>%
  dplyr::filter( Probe_Type == "cg" ) %>%
  
  # For Set == T (Titration)
  dplyr::filter( Set == "T" ) %>%
  dplyr::filter( Work_Str == "ND-O-PB" ) %>%
  dplyr::filter( Min_Pval == 1 ) %>%
  # dplyr::filter( Grp == "UH" ) %>%

  dplyr::group_by( Work_Str,
                   
                   Experiment, 
                   Grp,Set,
                   
                   Min_Pval,
                   Chop_Cut,Trif_Cut, # Ilmn_Cut,
                   
                   Probe_Type ) %>%
  dplyr::summarise( 
    tot_cnt = n(),
    
    # dBp_cnt = sum( db_pas_per > 75, na.rm=TRUE ),
    # dBp_avg = mean( db_pas_per, na.rm=TRUE ),
    dBp_med = median( db_pas_per, na.rm=TRUE ),
    dBp_mad = mad( db_pas_per, na.rm=TRUE ),
    
    # dBm_cnt = sum( dB_med > 0.75, na.rm=TRUE ),
    # dBm_avg = mean( dB_med, na.rm=TRUE ),
    dBm_med = median( dB_med, na.rm=TRUE ),
    dBm_mad = mad( dB_med, na.rm=TRUE ),
    
    # mPm_cnt = sum( mP_med > 0.75, na.rm=TRUE ),
    # mPm_avg = mean( mP_med, na.rm=TRUE ),
    mPm_med = median( mP_med, na.rm=TRUE ),
    mPm_mad = mad( mP_med, na.rm=TRUE ),
    
    # wSp_cnt = sum( wS_per > 0.75, na.rm=TRUE ),
    # wSp_avg = mean( wS_per, na.rm=TRUE ),
    wSp_med = median( wS_per, na.rm=TRUE ),
    wSp_mad = mad( wS_per, na.rm=TRUE ),
    
    .groups = "drop"
  )

uhm_dbs_sum %>% dplyr::filter( Chop_Cut != "U" ) %>% 
  dplyr::mutate( 
    dBp = dplyr::case_when( 
      Chop_Cut == "N" ~ dBp_med - dBp_mad,
      Chop_Cut == "Y" ~ dBp_med + dBp_mad,
      TRUE ~ NA_real_ ),
    dBm = dplyr::case_when( 
      Chop_Cut == "N" ~ dBm_med - dBm_mad,
      Chop_Cut == "Y" ~ dBm_med + dBm_mad,
      TRUE ~ NA_real_ ),
    mPm = dplyr::case_when( 
      Chop_Cut == "N" ~ mPm_med - mPm_mad,
      Chop_Cut == "Y" ~ mPm_med + mPm_mad,
      TRUE ~ NA_real_ ),
    wSp = dplyr::case_when( 
      Chop_Cut == "N" ~ wSp_med - wSp_mad,
      Chop_Cut == "Y" ~ wSp_med + wSp_mad,
      TRUE ~ NA_real_ )
  ) %>% dplyr::select( Work_Str,Grp,
                       Experiment,
                       Chop_Cut,Trif_Cut,tot_cnt, dBp,dBm,mPm,wSp ) %>%
  print(n=1000)


rep_dbs_sum <- NULL
rep_dbs_sum <- sel_dbs_tib %>%
  dplyr::filter( Probe_Type == "cg" ) %>%
  
  # For Set != T (Replicate)
  dplyr::filter( Set != "T" ) %>%
  dplyr::filter( Work_Str == "NCD-O-PB" ) %>%
  # dplyr::filter( Work_Str == "ND-O-PB" ) %>%
  dplyr::filter( Min_Pval == 0.5 ) %>%

  dplyr::group_by( Work_Str,
                   
                   # Experiment, 
                   Grp,Set,
                   
                   Min_Pval,
                   Chop_Cut,Trif_Cut, #Ilmn_Cut,
                   
                   Probe_Type ) %>%
  dplyr::summarise( 
    tot_cnt = n(),
    
    # dBp_cnt = sum( db_pas_per > 75, na.rm=TRUE ),
    # dBp_avg = mean( db_pas_per, na.rm=TRUE ),
    dBp_med = median( db_pas_per, na.rm=TRUE ),
    dBp_mad = mad( db_pas_per, na.rm=TRUE ),
    
    # dBm_cnt = sum( dB_med > 0.75, na.rm=TRUE ),
    # dBm_avg = mean( dB_med, na.rm=TRUE ),
    dBm_med = median( dB_med, na.rm=TRUE ),
    dBm_mad = mad( dB_med, na.rm=TRUE ),
    
    # mPm_cnt = sum( mP_med > 0.75, na.rm=TRUE ),
    # mPm_avg = mean( mP_med, na.rm=TRUE ),
    mPm_med = median( mP_med, na.rm=TRUE ),
    mPm_mad = mad( mP_med, na.rm=TRUE ),
    
    # wSp_cnt = sum( wS_per > 0.75, na.rm=TRUE ),
    # wSp_avg = mean( wS_per, na.rm=TRUE ),
    wSp_med = median( wS_per, na.rm=TRUE ),
    wSp_mad = mad( wS_per, na.rm=TRUE ),
    
    .groups = "drop"
  )


# For Set != T (Replicate)
top_rep_tib <- NULL
top_rep_tib <- sel_dbs_tib %>%
  # dplyr::filter( Probe_Type == "cg" ) %>%
  dplyr::filter( Set != "T" ) %>%
  dplyr::filter( Work_Str == "NCD-O-PB" ) %>%
  # dplyr::filter( Work_Str == "ND-O-PB" ) %>%
  dplyr::filter( Min_Pval == 0.5 ) %>% 
  dplyr::group_by( Probe_ID,Chop_Cut,Trif_Cut ) %>%
  dplyr::summarise( wSp_max = max(wS_per, na.rm=TRUE), 
                    dBm_max = max(dB_med, na.rm=TRUE),
                    dBp_max = max(db_pas_per, na.rm=TRUE),
                    mPm_max = max(mP_med, na.rm=TRUE),
                    .groups = "drop" )

sWp_cuts <- c( 90,95,96,97,98,99 )
top_rep_sels <- NULL
for ( sWp_cut in sWp_cuts ) {
  if ( p1 ) cat(glue::glue("{pmssg} sWp_cut: = '{sWp_cut}'{RET}"))

  top_rep_sel <- NULL
  top_rep_sel <- top_rep_tib %>% 
    dplyr::group_by( Chop_Cut,Trif_Cut ) %>%
    dplyr::mutate( Flag = dplyr::case_when(
      wSp_max < sWp_cut ~ "Y",
      TRUE ~ "N" )
    ) %>% dplyr::ungroup()
  
  top_rep_sels <- top_rep_sels %>% dplyr::bind_rows(
    top_rep_sel %>%
    dplyr::group_by( Chop_Cut,Flag,Trif_Cut ) %>%
    dplyr::summarise( sWp_cut = sWp_cut,
                      tot_cnt = n(), 
                      wSp_avg = mean( wSp_max,na.rm=TRUE ),
                      wSp_med = median( wSp_max,na.rm=TRUE ),
                      
                      dBm_avg = mean( dBm_max,na.rm=TRUE ),
                      dBm_med = median( dBm_max,na.rm=TRUE ),
                      
                      dBp_avg = mean( dBp_max,na.rm=TRUE ),
                      dBp_med = median( dBp_max,na.rm=TRUE ),
                      
                      mPm_avg = mean( mPm_max,na.rm=TRUE ),
                      mPm_med = median( mPm_max,na.rm=TRUE ),
                      
                      # pas_cnt = sum( wSp_max < 92 ), 
                      .groups = "drop" )  
  )
}

top_rep_sels %>% split(.$sWp_cut)

fin_dbs_prb_ids_csv <- file.path( opt$out_path, "MSA_Removal_List.ProbeID.20102023.sorted.csv" )
fin_dbs_prb_ids_tib <- NULL
fin_dbs_prb_ids_tib <- dplyr::bind_rows( 
  chop_cuts_tib %>% dplyr::select( Probe_ID ) %>%
    dplyr::distinct(),
  top_rep_sel %>% dplyr::filter(Flag == "Y" & Chop_Cut == "N" ) %>% 
    dplyr::select( Probe_ID ) %>% dplyr::distinct()  ) %>%
  dplyr::distinct( Probe_ID ) %>%
  dplyr::arrange( Probe_ID )
readr::write_csv( x = fin_dbs_prb_ids_tib, file = fin_dbs_prb_ids_csv )


# A tibble: 5 × 5
# Chop_Cut Trif_Cut tot_cnt wSp_avg wSp_med
# <fct>    <fct>      <int>   <dbl>   <dbl>
# 1 N        1         271944    99.1    99.2
# 2 N        3          11520    99.1    99.3
# 3 U        U             12    99.1    99.3
# 4 Y        1           8850    96.4    96.5
# 5 Y        3           1323    96.7    96.9


# 8850+1323 = 10173
# 7340+509 =  7849
# 1679+128 = 1807

# 293649 - 18022 = 275627
# 293649 - 10173 = 283476

# 293649 - 10173 - 1807 = 281669

rep_dbs_sum %>% 
  dplyr::filter( Chop_Cut != "U" ) %>% 
  dplyr::select( Work_Str,Grp,
                 # Experiment,
                 Chop_Cut,Trif_Cut,tot_cnt, 
                 dBm_med, dBm_mad, mPm_med, mPm_mad, wSp_med, wSp_mad) %>% 
  split(.$Grp)
# split(.$Experiment)

# [TBD]: Only Filter on Chop_Cut
# [TBD]: Add Mad,Min,Max to Titration set properly for each vairable...
rep_dbs_sum %>% dplyr::filter( Chop_Cut != "U" ) %>%
  dplyr::mutate( 
    dBp = dplyr::case_when( 
      Chop_Cut == "N" ~ dBp_med - dBp_mad,
      Chop_Cut == "Y" ~ dBp_med + dBp_mad,
      TRUE ~ NA_real_ ),
    dBm = dplyr::case_when( 
      Chop_Cut == "N" ~ dBm_med - dBm_mad,
      Chop_Cut == "Y" ~ dBm_med + dBm_mad,
      TRUE ~ NA_real_ ),
    mPm = dplyr::case_when( 
      Chop_Cut == "N" ~ mPm_med - mPm_mad,
      Chop_Cut == "Y" ~ mPm_med + mPm_mad,
      TRUE ~ NA_real_ ),
    wSp = dplyr::case_when( 
      Chop_Cut == "N" ~ wSp_med - wSp_mad,
      Chop_Cut == "Y" ~ wSp_med + wSp_mad,
      TRUE ~ NA_real_ )
  ) %>% dplyr::select( Work_Str,Experiment,Chop_Cut,Trif_Cut,tot_cnt, dBp,dBm,mPm,wSp ) %>%
  print(n=1000)


# Filter for target experiments... Add Work_Str
# dplyr::filter( (Min_Pval == 0.5 & Set != "T") | (Min_Pval == 1 & Set == "T") )

# sel_dbs_sum %>% dplyr::filter( Probe_Type == "cg" )

# sel_dbs_tib %>% print_sum( vec=c("Grp") )








# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Pre-processing:: Load Training Data: Alpha
#
# [Description Note]: The code below validated the plotting used to select
#  the optimal training set for screening!
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #           Pre-processing:: Set Training Data Parameters: Alpha
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  outliers_vec <- NULL
  outliers_vec <- c( 
    "rm-outliers",
    "ok-outliers",
    NULL
  )
  
  min_pval_vec <- NULL
  min_pval_vec <- c( 
    "pval-1", 
    # "pval-0.1", 
    "pval-0.5",
    NULL
  )
  
  mask_vec <- NULL
  mask_vec <- c(
    "mask_0",
    "mask_1",
    NULL
  )
  
  work_step_mat <- NULL
  work_step_mat <- c( 
    paste("DC", "O","PB", sep="-" ),
    paste("D", "O","PB",  sep="-" ),
    
    paste("CD", "O","PB", sep="-" ),
    
    paste("NDC","O","PB", sep="-" ),
    paste("NCD","O","PB", sep="-" ),
    
    paste("ND","O","PB",  sep="-" ),
    paste("DC", "O","pb", sep="-" ),
    paste("CD", "O","pb", sep="-" ),
    NULL
  )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                 Processing:: Load Training Data: Alpha
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  dbs_file_tib <- NULL
  dbs_file_tib <- expand_grid( 
    train_sdf_path, outliers_vec,min_pval_vec,mask_vec,work_step_mat ) %>% 
    tidyr::unite( Path, dplyr::everything(), sep = "/", remove = FALSE )
  
  dbs_file_len <- dbs_file_tib %>% length()
  if ( p1 ) cat(glue::glue("{pmssg} Experiment Counts: = '{dbs_file_len}'{RET}"))
  
  all_dbs_tib <- NULL
  all_dbs_sum <- NULL
  all_sum_tab <- NULL
  for ( ii in c(1: base::nrow(dbs_file_tib)) ) {
    if ( p1 ) cat(glue::glue("{pmssg} Current Path[{ii}] = '{dbs_file_tib$Path[ii]}'{RET}"))
    
    cur_dbs_csvs <- NULL
    cur_dbs_csvs <- file_list( path = dbs_file_tib$Path[ii],
                               # subs = c("rm-outliers","ok-outliers"),
                               # prefix = NULL,
                               suffix = ".calc_dBs.csv.gz",
                               pattern = ".calc_dBs.csv.gz$",
                               recursive = TRUE )
    
    params_tib <- NULL
    params_tib <- tibble::tibble(
      Outliers = dplyr:: case_when(
        dbs_file_tib$outliers_vec[ii] == "rm-outliers" ~ as.integer(0),
        dbs_file_tib$outliers_vec[ii] == "ok-outliers" ~ as.integer(1),
        TRUE ~ NA_integer_ ),
      Min_Pval = as.integer( dbs_file_tib$min_pval_vec[ii] %>% 
                               stringr::str_remove("^pval-") %>% 
                               as.double() * 100 ),
      Mask_Idx = as.integer( dbs_file_tib$mask_vec[ii] %>% 
                               stringr::str_remove("^mask_") %>% 
                               as.double() ),
      Work_Str = dbs_file_tib$work_step_mat[ii] )
    
    
    # [TBD]: Reduce Boolean Fields [DONE Below...]
    # [Done]: Us mcapply to multi-thread... [Testing]...
    #
    # [TBD]: In Previous (Alpha.R) script; clean up the output...
    cur_dbs_tib <- NULL
    cur_dbs_tib <- cur_dbs_csvs %>% # head() %>%
      mclapply( function(x) {
        readr::read_csv( file = x, show_col_types=FALSE ) %>%
          dplyr::mutate( Set = Set %>% as.character() %>% stringr::str_sub(1,1) ) %>%
          dplyr::filter( Probe_Type != "ct" ) %>%
          dplyr::filter( !(Set == "T" & Probe_Type == "rs") ) %>%
          dplyr::filter( !(Set == "T" & Probe_Type == "ch") ) %>%
          dplyr::mutate( 
            Chop_Cut = dplyr::case_when(
              Probe_ID %in% chop_cuts_tib$Probe_ID ~ "Y",
              TRUE ~ "N" ),
            # Work_Str = dbs_file_tib$work_step_mat[ii],
            Grp = dplyr::case_when( 
              Set != "T" ~ Set, 
              TRUE ~ as.character(Sam) )
          )
      }) %>% dplyr::bind_rows() %>% clean_tib() # %>% dplyr::select( Work_Str,Grp, dplyr::everything() )
    
    # [TBD]: Remove/Reduce Extra Columns::
    # all_dbs_tib <- all_dbs_tib %>% dplyr::bind_rows( cur_dbs_tib )
    
    #
    # [TBD]: Get General Stats::
    #
    cur_sum_tib <- NULL
    cur_sum_tib <- cbind( 
      params_tib, cur_dbs_tib %>%
        dplyr::group_by( Probe_Type,
                         Set,
                         Grp,
                         Chop_Cut ) %>%
        summarise( Tot_Cnt = n(),
                   
                   # [NOTE]: 'dBp' Delta Beta (dB) Passing Percent (Call Rate)
                   #
                   # dBp_Min = min(db_pas_per, na.rm=TRUE),
                   # dBp_Sds = sd(db_pas_per, na.rm=TRUE),
                   dBp_Avg = mean(db_pas_per, na.rm=TRUE),
                   # dBp_Med = median(db_pas_per, na.rm=TRUE),
                   # dBp_Mad = mad(db_pas_per, na.rm=TRUE),
                   # dBp_Max = max(db_pas_per, na.rm=TRUE),
                   
                   # [NOTE]: 'dBm' Delta Beta (dB) Median
                   #
                   # dBm_Min = min(dB_med, na.rm=TRUE),
                   #dBm_Sds = sd(dB_med, na.rm=TRUE),
                   dBm_Avg = mean(dB_med, na.rm=TRUE),
                   # dBm_Med = median(dB_med, na.rm=TRUE),
                   # dBm_Mad = mad(dB_med, na.rm=TRUE),
                   # dBm_Max = max(dB_med, na.rm=TRUE),
                   
                   # [NOTE]: 'mPm' Min Pval (mP) Median
                   #
                   # mPm_Min = min(mP_med, na.rm=TRUE),
                   # mPm_Sds = sd(mP_med, na.rm=TRUE),
                   mPm_Avg = mean(mP_med, na.rm=TRUE),
                   # mPm_Med = median(mP_med, na.rm=TRUE),
                   # mPm_Mad = mad(mP_med, na.rm=TRUE),
                   # mPm_Max = max(mP_med, na.rm=TRUE),
                   
                   # [NOTE]: 'wSp' Weight Score (wS) Percent
                   #
                   # wSp_Min = min(wS_per, na.rm=TRUE),
                   # wSp_Sds = sd(wS_per, na.rm=TRUE),
                   wSp_Avg = mean(wS_per, na.rm=TRUE),
                   # wSp_Med = median(wS_per, na.rm=TRUE),
                   # wSp_Mad = mad(wS_per, na.rm=TRUE),
                   # wSp_Max = max(wS_per, na.rm=TRUE),
                   
                   .groups = "drop" ) ) %>% tibble::as_tibble()
    # cur_dbs_sum %>% print( n=base::nrow(cur_dbs_sum) )
    
    # Gather Summary Tibble::
    all_dbs_sum <- all_dbs_sum %>% dplyr::bind_rows(cur_sum_tib)
    
    # [TBD]: Get Min and Max {Y/N}? Maybe add more summary variables??>
    # [TBD]: Plot on individual datasets... by Grp...
    cur_sum_tab <- NULL
    cur_sum_tab <- cur_dbs_sum %>%
      tidyr::pivot_longer( cols = c(Tot_Cnt, dBp_Avg, dBm_Avg, mPm_Avg, wSp_Avg) ) %>% 
      dplyr::arrange( name ) %>%
      tidyr::pivot_wider( names_from = c(name,Chop_Cut), values_from = c(value) )
    # cur_sum_tab %>% print(n=1000)
    
    # Gather Summary Table::
    all_sum_tab <- all_sum_tab %>% dplyr::bind_rows(cur_sum_tab)
    
    # break
  }
  
  # [TBD]: Round variables or write rds...
  sum_path <- safe_mkdir( file.path(opt$out_path, "summary") )
  all_dBs_sum_tib_csv <- file.path( sum_path, "all_dBs_sum_tib.csv.gz" )
  all_dBs_sum_tab_csv <- file.path( sum_path, "all_dBs_sum_tab.csv.gz" )
  
  readr::write_csv( x = all_dbs_sum, file = all_dBs_sum_tib_csv )
  readr::write_csv( x = all_sum_tab, file = all_dBs_sum_tab_csv )
  
  # [TBD]: Rename all the data structures (all_dbs_tib,all_dbs_sum,all_sum_tab)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #.                  Post-processing:: Plotting:: Scratch
  #          Post-processing:: Compare/Plot Wanding vs. Training Data
  #
  # [Description Note]: The code below validated the plotting used to select
  #  the optimal training set for screening!
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # [TBD]: Pick Best (1-3) Methods based on Wanding Comparisons
  # [TBD]: Run Sesame Straight on all Alpha/DVT Data & EPIC Alpha Data with Best Methods
  #        - Generate Negative Control Selection
  #        - Generate Genotype Calls
  # [TBD]: Rate each Method based on dB and r2 inside (MSA) and accross (EPIC)
  
  # [NOTE]: 'ppp' Chop_Pass Beta Value (N/Y+N)
  # [NOTE]: 'dBp' Delta Beta (dB) Passing Percent (Call Rate)
  # [NOTE]: 'dBm' Delta Beta (dB) Median
  # [NOTE]: 'mPm' Min Pval (mP) Median
  # [NOTE]: 'wSp' Weight Score (wS) Percent
  
  # [NOTE]: C = NA[0-9]+ CEPH
  # [NOTE]: L = Cancer Cell Line (Hela,..)
  # [NOTE]: T = Titration
  
  fin_sum_tib <- NULL
  fin_sum_tib <- all_sum_tab %>% 
    dplyr::mutate( 
      ppp = round( Tot_Cnt_N / (Tot_Cnt_N + Tot_Cnt_Y), 4 ),
      dBp = dBp_Avg_N - dBp_Avg_Y,
      dBm = dBm_Avg_N - dBm_Avg_Y,
      mPm = mPm_Avg_N - mPm_Avg_Y,
      wSp = wSp_Avg_N - wSp_Avg_Y
    )
  
  fin_sum_tab <- NULL
  fin_sum_tab <- fin_sum_tib %>% 
    # dplyr::filter( Probe_Type == "cg" ) %>%
    dplyr::select( Outliers,Min_Pval,Mask_Idx,Work_Str,Probe_Type,Set,Grp, 
                   ppp, dBp, dBm, mPm, wSp
                   # Cnt_Rat, dBp_Dif, dBm_Dif, mPm_Dif, wSp_Dif
    ) %>%
    tidyr::pivot_longer( cols = c(ppp, dBp, dBm, mPm, wSp) )
  # c(Cnt_Rat, dBp_Dif, dBm_Dif, mPm_Dif, wSp_Dif) )
  
  #
  # Plots by Experiment Type::
  #
  grp_vec <- unique(fin_sum_tab$Grp)
  # [Done]:: Reduce Min_Pval, subtract min value for each group...
  for ( grp in grp_vec ) {
    if ( p1 ) cat(glue::glue("{pmssg} plotting '{grp}'{RET}"))
    
    plot_tib <- NULL
    plot_tib <- fin_sum_tab %>%
      # dplyr::mutate( PlotGrp = paste(Grp,Work_Str, sep="\n") ) %>%
      dplyr::mutate( 
        Work_Str = Work_Str %>% stringr::str_replace_all("-","\n"),
        Work_Str = paste0( Min_Pval,"\n",Work_Str )
      ) %>%
      dplyr::filter( Grp == grp ) %>%
      dplyr::filter( Probe_Type == "cg" )
    
    lim_tib <- NULL
    lim_tib <- plot_tib %>% dplyr::group_by( name ) %>% 
      dplyr::summarise( Min = min(value), Max=max(value) )
    
    ppp_max <- lim_tib %>% dplyr::filter(name == "ppp" ) %>% pull("Max")
    
    ppp_min <- lim_tib %>% dplyr::filter(name == "ppp" ) %>% pull("Min")
    dBp_min <- lim_tib %>% dplyr::filter(name == "dBp" ) %>% pull("Min")
    dBm_min <- lim_tib %>% dplyr::filter(name == "dBm" ) %>% pull("Min")
    mPm_min <- lim_tib %>% dplyr::filter(name == "mPm" ) %>% pull("Min")
    wSp_min <- lim_tib %>% dplyr::filter(name == "wSp" ) %>% pull("Min")
    
    plot_tib2 <- NULL
    plot_tib2 <- plot_tib %>% 
      dplyr::mutate( value = dplyr::case_when( 
        name == "ppp" ~ value - ppp_min, 
        name == "dBp" ~ value - dBp_min, 
        name == "dBm" ~ value - dBm_min, 
        name == "mPm" ~ value - mPm_min, 
        name == "wSp" ~ value - wSp_min, 
        TRUE ~ NA_real_ ) )
    # plot_tib2 %>% dplyr::filter( is.na(value) )
    
    plot_ggg <- plot_tib2 %>%
      ggplot2::ggplot( aes(Work_Str) ) + 
      ggplot2::geom_bar( aes(weight=value, group=Probe_Type, fill=Probe_Type) ) +
      # ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
      ggplot2::facet_grid(
        # rows = vars( Outliers,Mask_Idx,Work_Str),
        # rows = vars( Mask_Idx,Work_Str),
        rows = vars( name ),
        cols = vars( Grp,Mask_Idx,Outliers ),
        scales = "free_y"
        # scales = "fixed"
      ) +
      ggplot2::labs( title = paste0("Sample Group[3] = ",grp) )
    
    plot(plot_ggg)
  }
  
  #
  # NEW SINGLE PLOT
  #
  
  # [TBD]:: Reduce Min_Pval, subtract min value for each group...
  if ( TRUE ) {
    # if ( p1 ) cat(glue::glue("{pmssg} plotting '{grp}'{RET}"))
    
    fin_sum_tab2 <- NULL
    fin_sum_tab2 <- fin_sum_tab %>%
      dplyr::filter( Probe_Type == "cg" ) %>%
      dplyr::filter( name != "ppp" ) %>%
      dplyr::mutate( 
        Work_Str = Work_Str %>% stringr::str_replace_all("-","\n")
        # Work_Str = paste0( Min_Pval,"\n",Work_Str )
      )
    # dplyr::mutate( PlotGrp = paste(Grp,Work_Str, sep="\n") ) %>%
    
    #
    # Replicates::
    #
    plot_tibLC <- NULL
    plot_tibLC <- fin_sum_tab2 %>%
      dplyr::filter( Set != "T" ) %>%
      dplyr::filter( Min_Pval == 50 ) %>%
      dplyr::mutate( PlotGrp = "R" )
    
    lim_tibLC <- NULL
    lim_tibLC <- plot_tibLC %>% dplyr::group_by( name ) %>% 
      dplyr::summarise( Min = min(value), Max=max(value) )
    
    dBp_minLC <- lim_tibLC %>% dplyr::filter(name == "dBp" ) %>% pull("Min")
    dBm_minLC <- lim_tibLC %>% dplyr::filter(name == "dBm" ) %>% pull("Min")
    mPm_minLC <- lim_tibLC %>% dplyr::filter(name == "mPm" ) %>% pull("Min")
    wSp_minLC <- lim_tibLC %>% dplyr::filter(name == "wSp" ) %>% pull("Min")
    
    plot_tibLC <- plot_tibLC %>% 
      dplyr::mutate( value = dplyr::case_when( 
        name == "dBp" ~ value - dBp_minLC, 
        name == "dBm" ~ value - dBm_minLC, 
        name == "mPm" ~ value - mPm_minLC, 
        name == "wSp" ~ value - wSp_minLC, 
        TRUE ~ NA_real_ ) )
    # plot_tibLC %>% dplyr::filter( is.na(value) )
    
    #
    # Titration::
    #
    # [TBD]: Break up into { HM,UH,UM }
    #
    plot_tibTT <- NULL
    plot_tibTT <- fin_sum_tab2 %>%
      dplyr::filter( Set == "T" ) %>%
      dplyr::filter( Min_Pval == 100 ) %>%
      dplyr::mutate( PlotGrp = "T" )
    
    lim_tibTT <- NULL
    lim_tibTT <- plot_tibTT %>% dplyr::group_by( name ) %>% 
      dplyr::summarise( Min = min(value), Max=max(value) )
    
    dBp_minTT <- lim_tibTT %>% dplyr::filter(name == "dBp" ) %>% pull("Min")
    dBm_minTT <- lim_tibTT %>% dplyr::filter(name == "dBm" ) %>% pull("Min")
    mPm_minTT <- lim_tibTT %>% dplyr::filter(name == "mPm" ) %>% pull("Min")
    wSp_minTT <- lim_tibTT %>% dplyr::filter(name == "wSp" ) %>% pull("Min")
    
    plot_tibTT <- plot_tibTT %>% 
      dplyr::mutate( value = dplyr::case_when( 
        name == "dBp" ~ value - dBp_minTT, 
        name == "dBm" ~ value - dBm_minTT, 
        name == "mPm" ~ value - mPm_minTT,
        name == "wSp" ~ value - wSp_minTT,
        TRUE ~ NA_real_ ) )
    # plot_tibTT %>% dplyr::filter( is.na(value) )
    
    plot_gggLC <- NULL
    plot_gggLC <- dplyr::bind_rows( plot_tibLC ) %>%
      ggplot2::ggplot( aes(Work_Str) ) + 
      ggplot2::geom_bar( aes(weight=value, group=Probe_Type, fill=Probe_Type) ) +
      # ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
      ggplot2::facet_grid(
        rows = vars( PlotGrp,Grp,name ),
        cols = vars( Mask_Idx,Outliers ),
        scales = "free_y"
        # scales = "fixed"
      ) +
      ggplot2::labs( title = paste0("Combined Plot: LC") )
    plot(plot_gggLC)
    
    plot_gggTT <- NULL
    plot_gggTT <- dplyr::bind_rows( plot_tibTT ) %>%
      ggplot2::ggplot( aes(Work_Str) ) + 
      ggplot2::geom_bar( aes(weight=value, group=Probe_Type, fill=Probe_Type) ) +
      # ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
      ggplot2::facet_grid(
        rows = vars( PlotGrp,Grp,name ),
        cols = vars( Mask_Idx,Outliers ),
        scales = "free_y"
        # scales = "fixed"
      ) +
      ggplot2::labs( title = paste0("Combined Plot: TT") )
    plot(plot_gggTT)
    
  }
  
  #
  # Some Old Stats:: Everything improtant is above...
  #
  
  fin_sum_list <- NULL
  fin_sum_list <- fin_sum_tib %>% 
    dplyr::arrange( 
      -dBp_Dif,
      -dBm_Dif,
      -mPm_Dif,
      -wSp_Dif,
      Grp,Probe_Type ) %>%
    split(.$Probe_Type)
  
  fin_sum_list %>% print()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #.                  Post-processing:: Plotting:: Scratch
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( FALSE ) {
    
    # [TBD]: Plot each Set Separately
    wS_per_den_ggg <- NULL
    wS_per_den_ggg <- cur_dbs_tib %>% 
      # ggplot2::ggplot( aes( x=log(wS_per), color=Chop_Cut, fill=Chop_Cut) ) +
      ggplot2::ggplot( aes( x=wS_per, color=Chop_Cut, fill=Chop_Cut) ) +
      ggplot2::geom_density( alpha=0.3 ) +
      ggplot2::facet_grid(
        rows=vars(Experiment),
        # cols=vars(Set), 
        scales = "fixed"
      )
    
    wS_per_hst_ggg <- NULL
    wS_per_hst_ggg <- cur_dbs_tib %>% 
      # ggplot2::ggplot( aes( x=db_pas_per,y=wS_per), group=Chop_Cut, color=Chop_Cut, fill=Chop_Cut ) +
      ggplot2::ggplot( aes( x=wS_per,y=db_pas_per, group=Chop_Cut ), color=Chop_Cut, fill=Chop_Cut ) +
      ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
      ggplot2::facet_grid(
        rows=vars(Experiment),
        # rows=vars(Experiment),
        # cols=vars(Set), 
        scales = "fixed"
      )
    wS_per_hst_ggg
    
    cur_dbs_sum %>%
      ggplot2::ggplot( aes( x=wSp_Avg, color=Chop_Cut, fill=Chop_Cut) ) +
      ggplot2::geom_density( alpha=0.3 ) +
      ggplot2::facet_grid(
        rows=vars(Experiment),
        # cols=vars(Set), 
        scales = "fixed"
      )
    
    all_dbs_tib %>% 
      dplyr::filter( Work_Str   == "DC-O-PB" ) %>%
      dplyr::filter( Probe_Type == "cg" ) %>%
      dplyr::filter( Work_Str == "DC-O-PB" ) %>%
      dplyr::filter( Work_Str == "DC-O-PB" ) %>%
      ggplot2::ggplot( aes( x=wS_per,y=db_pas_per, group=Chop_Cut ), color=Chop_Cut, fill=Chop_Cut ) +
      ggplot2::geom_boxplot( outlier.alpha = 0.2, varwidth = TRUE ) +
      ggplot2::facet_grid(
        rows=vars(Experiment),
        # rows=vars(Experiment),
        # cols=vars(Set), 
        scales = "fixed"
      )
    
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #.                  Post-processing:: Plotting:: Scratch
  #          Post-processing:: Compare/Plot Wanding vs. Training Data
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( FALSE ) {
    top_dbs_sum <- NULL
    top_dbs_sum <- all_dbs_tib %>%
      dplyr::group_by( Work_Str,Remove_Outliers,Min_Pval,Mask_Index,Probe_Type,
                       # Experiment,
                       Set,
                       Chop_Cut ) %>%
      summarise( Tot_Cnt = n(),
                 # wSp_Min = min(wS_per, na.rm=TRUE),
                 wSp_Sds = sd(wS_per, na.rm=TRUE),
                 wSp_Avg = mean(wS_per, na.rm=TRUE),
                 wSp_Med = median(wS_per, na.rm=TRUE),
                 wSp_Mad = mad(wS_per, na.rm=TRUE),
                 # wSp_Max = max(wS_per, na.rm=TRUE), 
                 .groups = "drop" )
    
    top_dbs_sum %>% print(n=1000)
    
    top_dbs_sum_tab <- NULL
    top_dbs_sum_tab <- top_dbs_sum %>% 
      dplyr::filter( Probe_Type != "ct" ) %>%
      dplyr::select( Work_Str,
                     Remove_Outliers,
                     Min_Pval,
                     Mask_Index,
                     Probe_Type,
                     Set,
                     
                     # Experiment,
                     Chop_Cut,
                     wSp_Avg ) %>% 
      tidyr::pivot_wider( names_from = c(Chop_Cut), 
                          values_from = c(wSp_Avg) ) %>% 
      dplyr::mutate( Dif = N - Y )
    
    top_dbs_sum_tab %>% dplyr::arrange(-Dif) %>% print(n=1000)
    
    # Summary of Summary::
    top_dbs_sum_tab %>% 
      dplyr::filter( Min_Pval == 0.5 ) %>%
      dplyr::filter( !(Set == "T" & Probe_Type == "rs") ) %>%
      dplyr::filter( !(Set == "T" & Probe_Type == "ch") ) %>%
      dplyr::group_by( Probe_Type,
                       # Work_Str,
                       # Remove_Outliers,
                       Set ) %>% 
      dplyr::summarise( 
        Dif_Avg = mean( Dif, na.rm=TRUE), 
        Dif_Med = median( Dif, na.rm=TRUE), 
        .groups = "drop" ) %>%
      dplyr::arrange( -Dif_Med, Set ) %>% print(n=1000)
    
    all_dbs_tib %>% 
      dplyr::filter( Min_Pval == 0.5 ) %>%
      dplyr::filter( !(Set == "T" & Probe_Type == "rs") ) %>%
      dplyr::filter( !(Set == "T" & Probe_Type == "ch") ) %>%
      dplyr::group_by( Probe_Type,
                       # Work_Str,
                       # Remove_Outliers,
                       Set )
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
