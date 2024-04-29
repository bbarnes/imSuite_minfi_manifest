
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
par$prgm_tag <- 'stable_priym_analysis'
par$verbose  <- 3
local_paths  <- c( 
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
par$version <- 2
par$version <- 3
par$version <- 4
par$version <- 5

# par$run_name <- "EPICv1"
# par$run_name <- "Embarkv1"
# par$run_name <- "FAILv1"
# par$run_name <- "COREv1"
#
# par$run_name <- "MSAv03"
par$run_name <- "EPICv2"

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
# opt$rcpp <- 1
# opt$rcpp <- 2
# opt$rcpp <- 3
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
#
#                 Pre-processing:: Load Sample Sheet
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sam_tib <- NULL

use_new_ssh <- FALSE
use_new_ssh <- TRUE
if ( !use_new_ssh ) {
  sam_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/for_Priyam/sample_sheets/Sample_Sheet_SKY_BBarnes.csv" )
  sam_tib <- readr::read_csv( file = sam_csv, show_col_types = FALSE ) %>% 
    dplyr::select( Sentrix_Name,Sample_Name,Sample_Group, 
                   Intervention,Gender,age ) %>%
    dplyr::rename( Sample_Gender = Gender,
                   Sample_Age = age ) %>% 
    dplyr::mutate( 
      Intervention = dplyr::case_when( Sample_Name == "Post_DNA_21" ~ "SKY", 
                                       TRUE ~ Intervention ),
      Intervention = dplyr::case_when( Sentrix_Name == "207912280104_R02C01" ~ "Control",
                                       TRUE ~ Intervention ),
      Intervention = dplyr::case_when( Sentrix_Name == "208019130005_R06C01" ~ "SKY",
                                       TRUE ~ Intervention )
    )
  
} else {
  
  # Should have included the sex and gender from the corrected sample sheet...
  sam_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/for_Priyam/sample_sheets/SS_CORRECTED1.csv.gz" )
  sam_tib <- readr::read_csv( file = sam_csv, show_col_types = FALSE ) %>% 
    dplyr::mutate( Sentrix_Name = paste( Sentrix_ID,Sentrix_Position, sep="_") ) %>%
    dplyr::select( Sentrix_Name,Sample_Name,Sample_Group, 
                   Intervention ) %>%
    dplyr::mutate( 
      Intervention = dplyr::case_when( Sample_Name == "Post_DNA_21" ~ "SKY", 
                                       TRUE ~ Intervention ),
      Intervention = dplyr::case_when( Sentrix_Name == "207912280104_R02C01" ~ "Control",
                                       TRUE ~ Intervention ),
      Intervention = dplyr::case_when( Sentrix_Name == "208019130005_R06C01" ~ "SKY",
                                       TRUE ~ Intervention ),
      Intervention = dplyr::case_when( Sentrix_Name == "207910560089_R01C01" ~ "SKY",
                                       TRUE ~ Intervention ),
      Sample_Gender = 2,
      Sample_Age = 100 )
}

# Sanity Check::
sam_tib %>% dplyr::filter( Sentrix_Name == "208019130005_R06C01" ) %>% print()
sam_tib %>% dplyr::filter( Sentrix_Name == "207912280104_R02C01" ) %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 Pre-processing:: Load Swifthoof Results
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssh_list <- NULL
# ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/EPICv2_Priyam/full/swifthoof_main" )
# ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/EPICv2_Priyam_raw2/full/swifthoof_main" )
ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/join/EPICv2_Priyam_all" )
ssh_list <- file_list( path    = ssh_path, 
                       suffix  = "_EPIC_A1_AutoSampleSheet.csv.gz", 
                       pattern = "_EPIC_A1_AutoSampleSheet.csv.gz$", 
                       recursive = TRUE )

#
# ssh_list <- NULL
# ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/EPICv2_Priyam_ind/full/swifthoof_main" )
# ssh_list <- file_list( path    = ssh_path, 
#                        suffix  = "_EPIC_A1_AutoSampleSheet.csv.gz", 
#                        pattern = "_EPIC_A1_AutoSampleSheet.csv.gz$" )

ssh_tib <- NULL
ssh_tib <- ssh_list %>%
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows() %>% 
  dplyr::inner_join( sam_tib, by=c("Sentrix_Name") ) %>% 
  dplyr::mutate( 
    Sample_State = Sample_Name %>% stringr::str_remove("_.*$"),
    Sample_Index = Sample_Name %>% stringr::str_remove("^.*_") %>% as.integer(),
    Sample_Class = dplyr::case_when(
      Sample_State == "Pre"  ~ "A",
      Sample_State == "Mid"  ~ "B",
      Sample_State == "Post" ~ "C",
      TRUE ~ "z" ),
    Sample_Group = Intervention %>% stringr::str_sub(1,1),
    Sample_Name_Src = Sample_Name,
    Sample_Name = paste( Sample_Group,Sample_Class,Sample_Index, sep="_" )
  ) %>%
  dplyr::add_count( Sample_Index, name = "Index_Count" ) %>%
  dplyr::arrange( Intervention,Sample_Index,Sample_Class ) %>%
  dplyr::select( Sample_Name,Sample_Group,Sample_Class,Sample_Index,
                 Sample_Gender,Sample_Age,Index_Count,
                 Intervention,
                 dplyr::everything() )
# ssh_tib %>% head() %>% as.data.frame()

# [FIXED]: Sample Index == 21 has mixed Sample_Groups...
# ssh_tib %>% dplyr::group_by( Sample_Index,Sample_Group ) %>% dplyr::summarise( Count=n(), .groups = "drop" ) %>% dplyr::add_count( Sample_Index, name = "Cross_Cnt" ) %>% dplyr::filter( Cross_Cnt != 1 )

# sel_ssh_tib %>% dplyr::filter( Sample_Index == 19 )
# sel_ssh_tib %>% dplyr::filter( Sample_Index == 21 )

# This should be fixed above...
if ( FALSE ) {
  if ( use_new_ssh ) {
    ssh_tib <- ssh_tib %>% 
      dplyr::mutate( 
        Intervention = dplyr::case_when( Sample_Index == 21 ~ "Control", 
                                         TRUE ~ Intervention )
      ) # %>% dplyr::filter( Sample_Index == 21 )
    ssh_tib <- ssh_tib %>% 
      dplyr::mutate( 
        Intervention = dplyr::case_when( Sample_Index == 19 ~ "SKY", 
                                         TRUE ~ Intervention )
      ) # %>% dplyr::filter( Sample_Index == 19 )
  }
}

ssh_sum0 <- NULL
ssh_sum0 <- ssh_tib %>%
  dplyr::group_by( Intervention,Sample_Class,Index_Count ) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )

ssh_sum1 <- NULL
ssh_sum1 <- ssh_tib %>% 
  dplyr::group_by( Intervention,Index_Count ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) 

#
# [TBD]: Compare Sex_1 SexKaryotype_1 Ethnicity_1 AgeSkinBlood_1 
#                   Sample_Gender,Sex_1,SexKaryotype_1 ) %>%
#
sex_sum0 <- ssh_tib %>%
  dplyr::group_by( Intervention,
                   Sample_Gender,Sex_1 ) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )
print(sex_sum0)

sex_sum1 <- ssh_tib %>% dplyr::group_by( Sample_Gender ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
print(sex_sum1)
sex_sum2 <- ssh_tib %>% dplyr::group_by( Sex_1 ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
print(sex_sum2)

#
# Build Selected Sample Sheet::
#

#
# [TBD]: Fix Plot_ID using Sample_ID...
#
sel_ssh_tib <- NULL
sel_ssh_tib <- ssh_tib %>% 
  # dplyr::filter( (Sample_State=="Pre" | Sample_State=="Post") ) %>%
  # dplyr::add_count( Sample_Index, name = "Index_Count" ) %>%
  dplyr::filter( Index_Count >= 2 ) %>%
  dplyr::mutate( Sample_Base = Sample_Name ) %>%
  dplyr::group_by( Sample_Base,detect_version ) %>% 
  dplyr::mutate(
    Rep_Num = dplyr::row_number(),
    PPP_Int = cg_calls_pass_perc_1 %>% as.integer(),
    Plot_ID = paste( Sample_Name,PPP_Int, sep="_" ),
    # Plot_ID = paste( detect_version,Sample_Index,Rep_Num,PPP_Int, sep="_")
    # Plot_ID = paste( detect_version,Sample_ID,Rep_Num,PPP_Int, sep="_")
    # Plot_ID = dplyr::case_when( 
    #   opt$add_conc_name ~ paste( detect_version,Sample_ID,Concentration,Rep_Num,PPP_Int, sep="_"),
    #   TRUE ~ paste( detect_version,Sample_ID,Rep_Num,PPP_Int, sep="_") )
  ) %>% dplyr::ungroup() %>%
  dplyr::select( Sentrix_Name,Plot_ID, dplyr::everything() )

# Drop poor samples
sel_ssh_tib <- sel_ssh_tib %>% dplyr::filter( cg_calls_pass_perc_1 > 90 )
sel_ssh_tib <- sel_ssh_tib %>% dplyr::filter( Intervention != "Advanced" )

sel_ssh_sum <- NULL
sel_ssh_sum <- sel_ssh_tib %>% 
  dplyr::group_by( Intervention,Sample_State,Index_Count ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

# Failed samples with cg_calls_pass_perc_1 < 90
# "207912280101_R01C01" "207912280054_R08C01" "207912280104_R01C01" "207910560066_R08C01" "207910560089_R06C01"
failed_ssh_tib <- NULL
failed_ssh_tib <- ssh_tib %>% dplyr::anti_join( sel_ssh_tib, by=c("Sentrix_Name") )

# sel_ssh_tib %>% dplyr::arrange( Sample_Index ) %>% as.data.frame()
# sel_ssh_tib %>% dplyr::filter( Sample_Index == 3 ) %>% as.data.frame()

# ssh_tib %>% head() %>% as.data.frame()
# ssh_tib %>% dplyr::add_count( Sample_Index, name = "Index_Count" ) %>% dplyr::group_by( Index_Count ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Simple Plotting::
#      cg_pvals_PnegEcdf_pass_perc_1 vs. cg_pvals_pOOBAH_pass_perc_1
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_plot1 <- TRUE
run_plot1 <- FALSE

if ( run_plot1 ) {
  ssh_pval_pnt_ggg <- NULL
  ssh_pval_pnt_ggg <- ssh_tib %>%
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    dplyr::group_by( Sample_State ) %>%
    ggplot2::ggplot( aes( x=cg_pvals_PnegEcdf_pass_perc_1, 
                          y=cg_pvals_pOOBAH_pass_perc_1, color=Sample_State) ) +
    ggplot2::geom_point()
  # ssh_pval_pnt_ggg
  
  ssh_poob_den_ggg <- NULL
  ssh_poob_den_ggg <- ssh_tib %>%
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    dplyr::group_by( Sample_State ) %>%
    ggplot2::ggplot( aes( x=cg_pvals_pOOBAH_pass_perc_1,
                          color=Sample_State) ) +
    ggplot2::geom_density()
  # ssh_poob_den_ggg
  
  ssh_poob_den_box1_ggg <- NULL
  ssh_poob_den_box1_ggg <- ssh_tib %>% 
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    dplyr::filter( Intervention != "Advanced" ) %>%
    tidyr::pivot_longer( cols = c(cg_pvals_PnegEcdf_pass_perc_1,cg_pvals_pOOBAH_pass_perc_1), 
                         names_to = "Pval_Type", values_to = "Pval_Pass_Percent" ) %>%
    dplyr::group_by( Sample_State ) %>%
    ggplot2::ggplot( aes( x=Pval_Pass_Percent,
                          color=Sample_State) ) +
    ggplot2::geom_density() +
    ggplot2::facet_grid( rows = vars(Intervention),
                         cols = vars(Pval_Type), scales = "free" )
  # ssh_poob_den_box1_ggg
  
  ssh_poob_den_box2_ggg <- NULL
  ssh_poob_den_box2_ggg <- ssh_tib %>% 
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    dplyr::filter( Intervention != "Advanced" ) %>%
    tidyr::pivot_longer( cols = c(cg_pvals_PnegEcdf_pass_perc_1,cg_pvals_pOOBAH_pass_perc_1), 
                         names_to = "Pval_Type", values_to = "Pval_Pass_Percent" ) %>%
    dplyr::group_by( Sample_State ) %>%
    ggplot2::ggplot( aes( x=Pval_Pass_Percent,
                          color=Sample_State) ) +
    ggplot2::geom_density() +
    ggplot2::facet_grid( rows = vars(Intervention,Pval_Type),
                         scales = "free" )
  # ssh_poob_den_box2_ggg
  
  sel_poob_den_box2_ggg <- NULL
  sel_poob_den_box2_ggg <- sel_ssh_tib %>% 
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    dplyr::filter( Intervention != "Advanced" ) %>%
    tidyr::pivot_longer( cols = c(cg_pvals_PnegEcdf_pass_perc_1,cg_pvals_pOOBAH_pass_perc_1), 
                         names_to = "Pval_Type", values_to = "Pval_Pass_Percent" ) %>%
    dplyr::group_by( Sample_State ) %>%
    ggplot2::ggplot( aes( x=Pval_Pass_Percent,
                          color=Sample_State) ) +
    ggplot2::geom_density() +
    ggplot2::facet_grid( rows = vars(Intervention,Pval_Type),
                         scales = "free" )
  # sel_poob_den_box2_ggg
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Compare Sample Beta Values::
#                                   Raw
#
# Sample vs. Sample r2/dB
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

raw_list <- NULL
# raw_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/EPICv2_Priyam/full/swifthoof_main" )
raw_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/EPICv2_Priyam_raw2/full/swifthoof_main" )
raw_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/EPICv2_Priyam_all" )
raw_list <- file_list( path    = raw_path, 
                       suffix  = "_EPIC_A1_raw.call.dat.csv.gz", 
                       pattern = "_EPIC_A1_raw.call.dat.csv.gz$", 
                       recursive = TRUE )

use_ind <- FALSE
use_ind <- TRUE
if ( use_ind ) {
  
  # raw_list <- NULL
  # raw_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/EPICv2_Priyam_ind/full/swifthoof_main" )
  # raw_list <- file_list( path    = raw_path, 
  #                        suffix  = "_EPIC_A1_raw.call.dat.csv.gz", 
  #                        pattern = "_EPIC_A1_raw.call.dat.csv.gz$" )
  
  ind_list <- NULL
  ind_list <- file_list( path    = raw_path, 
                         suffix  = "_EPIC_A1_ind.call.dat.csv.gz", 
                         pattern = "_EPIC_A1_ind.call.dat.csv.gz$", 
                         recursive = TRUE )
  raw_list <- ind_list
}

# raw_beta_tib <- NULL
# raw_beta_tib <- raw_list[sel_ssh_tib$Sentrix_Name] %>%
#   lapply( readr::read_csv, show_col_types = FALSE ) %>%
#   dplyr::bind_rows( .id = "Sentrix_Name" )

sel_pair_list <- NULL
sel_pair_list <- sel_ssh_tib %>% split(.$Sample_Index)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#              Processing:: First Round Using top N hits
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load DML data::

dml_tib <- NULL
dml_txt <- file.path( opt$top_path, "Projects.new/EPIC_v2/for_Priyam/smty_sky.txt" )
dml_tib <- readr::read_delim( file = dml_txt, delim = " ", show_col_types = FALSE )

max_hit_vec <- c( 100, 1000, 10000, 100000, 600000 )
min_pval <- 0.05
str_pval <- paste0("pval-",min_pval)

for ( max_hit_cnt in max_hit_vec ) {
  cgn_dB_tib0 <- NULL
  for ( sel_idx in names(sel_pair_list) ) {
    
    cur_inv_str <- NULL
    # cur_inv_str <- sel_pair_list[[1]]$Intervention
    # cur_inv_str <- sel_pair_list[[sel_idx]]$Intervention %>% unique()
    cur_inv_str <- sel_pair_list[[sel_idx]]$Sample_Group %>% unique()
    
    if ( cur_inv_str != "S" ) next
    
    if ( length(cur_inv_str) > 1 ) {
      cat(glue::glue("{pmssg} Multiple values[{sel_idx}]: ({cur_inv_str}), Skipping...{RET}"))
      next
    }
    
    cur_idx_str <- NULL
    cur_idx_str <- sel_pair_list[[sel_idx]]$Sample_Index %>% unique()
    
    cur_sam_str <- NULL
    cur_sam_str <- paste( cur_inv_str,sel_idx, sep="_" )
    
    cur_pair_list <- NULL
    cur_pair_list <- sel_pair_list[[sel_idx]] %>% split(.$Sentrix_Name)
    
    mat_pre <- FALSE
    mat_pos <- FALSE
    
    # if ( sel_pair_list[[sel_idx]]$Sample_State[1] %>% stringr::str_detect("Pre") )
    #   mat_pre <- TRUE
    # if ( sel_pair_list[[sel_idx]]$Sample_State[1] %>% stringr::str_detect("Post") )
    #   mat_pos <- TRUE
    
    pair_data_tab <- NULL
    for ( sentrix_id in names(cur_pair_list) ) {
      cur_name_str <- NULL
      cur_name_str <- cur_pair_list[[sentrix_id]]$Sample_State
      
      cur_group_str <- NULL
      cur_group_str <- cur_pair_list[[sentrix_id]]$Sample_Group
      
      if ( cur_name_str %>% stringr::str_detect("Pre") )
        mat_pre <- TRUE
      if ( cur_name_str %>% stringr::str_detect("Post") )
        mat_pos <- TRUE
      
      unq_name_str <- NULL
      unq_name_str <- paste( cur_inv_str,cur_name_str,sel_idx, sep="_" )
      
      cur_data_tib <- NULL
      cur_data_tib <- readr::read_csv( raw_list[[sentrix_id]], 
                                       show_col_types = FALSE ) %>%
        dplyr::mutate( Sample_State = cur_name_str )
      
      pair_data_tab <- pair_data_tab %>% bind_rows(cur_data_tib)
    }
    
    if ( !mat_pre || !mat_pos ) next
    
    pair_data_tib <- NULL
    pair_data_tib <- pair_data_tab %>% 
      dplyr::filter( pvals_pOOBAH <= min_pval ) %>% 
      dplyr::select( Probe_ID,Sample_State,betas ) %>%
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sample_State), 
                          values_from = c(betas) ) %>% 
      dplyr::mutate( dB = Pre - Post,
                     Sample_Group = cur_inv_str,
                     Sample_Index = cur_idx_str ) %>%
      dplyr::arrange( -abs(dB) ) %>% head( n=max_hit_cnt )
    
    # Gather dB across samples...
    cgn_dB_tib0 <- cgn_dB_tib0 %>% 
      dplyr::bind_rows( dplyr::select( pair_data_tib, Probe_ID, Sample_Group, Sample_Index, dB ) )
    
    
    cat(glue::glue("{pmssg} r-squared({cur_sam_str}):{RET}"))
    # if ( sel_idx > 3 ) break
  }
  
  # Using both controls and SKY::
  # > cgn_dB_tib0 %>% dplyr::distinct( Probe_ID )
  # A tibble: 120,482 × 1
  #
  # Using SKY Only::
  # > cgn_dB_tib0 %>% dplyr::distinct( Probe_ID )
  # A tibble: 98,495 × 1
  # as.vector( dplyr::distinct( cgn_dB_tib0, Probe_ID ) )
  
  cgn_dB_vec0 <- NULL
  cgn_dB_vec0 <- as.vector( dplyr::distinct( cgn_dB_tib0, Probe_ID ) %>% dplyr::pull(Probe_ID) )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #              Processing:: Second Round Using top N hits
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  plot_path <- safe_mkdir( file.path( opt$out_path, "plots", paste0("top-",max_hit_cnt),str_pval ) )
  
  cgn_dB_tib <- NULL
  for ( sel_idx in names(sel_pair_list) ) {
    
    cur_inv_str <- NULL
    # cur_inv_str <- sel_pair_list[[1]]$Intervention
    # cur_inv_str <- sel_pair_list[[sel_idx]]$Intervention %>% unique()
    cur_inv_str <- sel_pair_list[[sel_idx]]$Sample_Group %>% unique()
    
    if ( length(cur_inv_str) > 1 ) {
      cat(glue::glue("{pmssg} Multiple values[{sel_idx}]: ({cur_inv_str}), Skipping...{RET}"))
      next
    }
    
    cur_idx_str <- NULL
    cur_idx_str <- sel_pair_list[[sel_idx]]$Sample_Index %>% unique()
    
    cur_sam_str <- NULL
    cur_sam_str <- paste( cur_inv_str,sel_idx, sep="_" )
    
    cur_pair_list <- NULL
    cur_pair_list <- sel_pair_list[[sel_idx]] %>% split(.$Sentrix_Name)
    
    mat_pre <- FALSE
    mat_pos <- FALSE
    
    pair_data_tab <- NULL
    for ( sentrix_id in names(cur_pair_list) ) {
      cur_name_str <- NULL
      cur_name_str <- cur_pair_list[[sentrix_id]]$Sample_State
      
      cur_group_str <- NULL
      cur_group_str <- cur_pair_list[[sentrix_id]]$Sample_Group
      
      if ( cur_name_str %>% stringr::str_detect("Pre") )
        mat_pre <- TRUE
      if ( cur_name_str %>% stringr::str_detect("Post") )
        mat_pos <- TRUE
      
      unq_name_str <- NULL
      unq_name_str <- paste( cur_inv_str,cur_name_str,sel_idx, sep="_" )
      
      cur_data_tib <- NULL
      cur_data_tib <- readr::read_csv( raw_list[[sentrix_id]], 
                                       show_col_types = FALSE ) %>%
        dplyr::mutate( Sample_State = cur_name_str )
      
      pair_data_tab <- pair_data_tab %>% bind_rows(cur_data_tib)
      
      #
      # LEFT OFF HERE: Not sure this is the right place to do this...
      # [TBD]: Join dB's across samples by CGN..
      #
      if ( FALSE ) {
        
        tmp_tib <- NULL
        tmp_tib <- dplyr::select( pair_data_tab,Probe_ID,dB ) %>%
          magrittr::set_names( c("Probe_ID","") ) %>% 
          magrittr::set_names( c("Probe_ID", unq_name_str) )
        
        if ( cgn_dB_tib %>% is.null() ) {
          cgn_dB_tib <- tmp_tib
        } else {
          cgn_dB_tib <- cgn_dB_tib %>% 
            dplyr::inner_join( tmp_tib, by=c("Probe_ID") )
        }
        
      }
    }
    
    if ( !mat_pre || !mat_pos ) next
    
    pair_data_tib <- NULL
    pair_data_tib <- pair_data_tab %>% 
      # dplyr::filter( pvals_pOOBAH < 0.05 ) %>% 
      dplyr::select( Probe_ID,Sample_State,betas ) %>%
      dplyr::filter( Probe_ID %in% cgn_dB_vec0 ) %>%
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sample_State), 
                          values_from = c(betas) ) %>% 
      dplyr::mutate( dB = Pre - Post,
                     Sample_Group = cur_inv_str,
                     Sample_Index = cur_idx_str ) %>%
      dplyr::arrange( -abs(dB) ) # %>% head( n=max_hit_cnt )
    
    pair_data_pdf <- file.path( plot_path, paste0("dB_abs.",cur_sam_str,".n-",max_hit_cnt,".",str_pval,".pdf") )
    pair_data_ggg <- NULL
    pair_data_ggg <- pair_data_tib %>%
      dplyr::arrange( -abs(dB) ) %>% head( n=max_hit_cnt ) %>%
      ggplot2::ggplot( aes( x=Pre, y=Post) ) +
      ggplot2::geom_point() + 
      ggplot2::geom_abline()
    ggplot2::ggsave( filename = pair_data_pdf, plot = pair_data_ggg, 
                     device = "pdf", width = 7, height = 7, dpi = 320)
    
    r2_mat = NULL
    r2_mat = pair_data_tib %>% 
      dplyr::select( -Sample_Group, -Sample_Index, -dB ) %>%
      tibble::column_to_rownames(var = "Probe_ID") %>% 
      as.data.frame() %>% cor( use = "pairwise.complete.obs", method = "pearson" )
    
    r2_sum <- NULL
    r2_sum <- r2_mat %>% as.data.frame() %>% 
      tibble::rownames_to_column(var = "Group" ) %>%
      tibble::as_tibble() %>% 
      tidyr::pivot_longer( cols = c(-Group), 
                           names_to = "Name", 
                           values_to = "r2" ) %>% 
      dplyr::filter( Group != Name )
    
    # [TBD]: Record r2 in tibble
    # [TBD]: Intersect Y-Probes and create average
    # [TBD]: Intersect X-Probes and create average
    # y_avg_int <- 
    
    # Gather dB across samples...
    cgn_dB_tib <- cgn_dB_tib %>% 
      dplyr::bind_rows( dplyr::select( pair_data_tib, Probe_ID, Sample_Group, Sample_Index, dB ) )
    
    cat(glue::glue("{pmssg} r-squared({cur_sam_str}):{RET}"))
    print(r2_sum)
    
    # if ( sel_idx > 3 ) break
  }
  
  cgn_dB_sum <- NULL
  cgn_dB_sum <- cgn_dB_tib %>% 
    dplyr::mutate( dB_abs = abs(dB) ) %>% 
    # dplyr::group_by( Sample_Group ) %>% 
    dplyr::group_by( Sample_Group,Sample_Index ) %>% 
    dplyr::summarise( Tot = n(), 
                      Min = min(dB_abs, na.rm = TRUE ), 
                      Max = max(dB_abs, na.rm = TRUE ), 
                      
                      Avg = mean(dB_abs, na.rm = TRUE ), 
                      Sds = sd(dB_abs, na.rm = TRUE ), 
                      Med = median(dB_abs, na.rm = TRUE ), 
                      Mad = mad(dB_abs, na.rm = TRUE ), 
                      .groups = "drop" )
  cgn_dB_sum %>% print(n=1000)
  
  # View Tibble
  # cgn_dB_tib %>% dplyr::group_by( Sample_Group ) %>% dplyr::summarise( Count = n(), .groups = "drop")  %>% print(n=1000)
  
  # Plot dB Density
  cgn_dB_ggg <- NULL
  cgn_dB_ggg <- cgn_dB_tib %>% 
    # ggplot2::ggplot( aes(x=dB, color=Sample_Group, group=Sample_Group, fill=Sample_Group) ) + 
    # ggplot2::ggplot( aes(x=abs(dB), color=Sample_Group, group=Sample_Group, fill=Sample_Group ) ) + 
    ggplot2::ggplot( aes(x=dB, color=Sample_Group, group=Sample_Group, fill=Sample_Group ) ) + 
    ggplot2::geom_density( alpha=0.2 )
  
  # Build Table
  cgn_dB_tab <- NULL
  cgn_dB_tab <- cgn_dB_tib %>% 
    tidyr::pivot_wider( id_cols = c(Probe_ID), 
                        names_from = c(Sample_Group,Sample_Index), names_sep = "_", 
                        values_from = c(dB) ) # %>% dplyr::arrange( Sky_1 )
  
}

tmp_join_dB_dml_tib <- NULL
tmp_join_dB_dml_tib <- cgn_dB_tib %>% dplyr::inner_join( dml_tib, by=c("Probe_ID"="Est_X.Intercept.") ) 

tmp_join_dB_dml_tib %>% ggplot2::ggplot( aes( x=dB, y=Est_Sample_GroupGroup_Mid, color=Sample_Group, group=Sample_Group ) ) + ggplot2::geom_point( )
tmp_join_dB_dml_tib %>% ggplot2::ggplot( aes( x=dB, y=Est_Sample_GroupGroup_Post, color=Sample_Group, group=Sample_Group ) ) + ggplot2::geom_point( )

# Sky Summary::
# cgn_dB_tib0 %>% head() %>% dplyr::add_count( Probe_ID, name="Sky_Top_Cnt" )
sky_db_cnt0 <- NULL
sky_db_cnt0 <- cgn_dB_tib0 %>% 
  dplyr::add_count( Probe_ID, name="Sky_Top_Cnt" )

sky_db_top0 <- NULL
sky_db_top0 <- sky_db_cnt0 %>% 
  dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>% dplyr::arrange( -Sky_Top_Cnt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Pre-processing:: Find Sample VCFs
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# [TBD]: Compare karyotype calls
#
# [TBD]: Fix Sample_Name,Rep_Num[0,1,2]
# [TBD]: Allow all samples
# [TBD]: Pick best mate from matrix
#
# [TBD]: Investigate missing ind data vs. raw or use/develop Rcpp version...
# [TBD]: Clean up code in functions...

# Order of Operations Next Steps:
#  1. Validate the first 8 that genotypes match Priyam's Excel
#  2. Look up dB on each pair 
#  3. Demonstrate that the Sky show the same high dB and controls don't show anything

plot_heatmaps <- TRUE
plot_heatmaps <- FALSE
if ( plot_heatmaps ) {
  raw_vcf_list <- NULL
  raw_vcf_list <- file_list( path    = raw_path, 
                             prefix  = raw_path,
                             suffix  = "_EPIC_.*_raw.snps.vcf", 
                             pattern = "_raw.snps.vcf$",
                             recursive = TRUE )
  
  snvs_raw_vcfs <- NULL
  snvs_raw_vcfs <- raw_vcf_list[sel_ssh_tib$Sentrix_Name]
  
  # opt$gts_mins <- c( 20 )
  opt$gts_mins <- c( 20, 10, 1, 0 )
  
  opt$sub_strs <- c( "All")
  # opt$sub_strs <- c( "ADD","MVP","All")
  
  for ( sub_str in opt$sub_strs ) {
    sub_vec <- NULL
    
    # [NOTE]: Left over code below from VA:MVP:EPICv2
    if ( sub_str == "ADD" ) sub_vec <- add_151_tib$IlmnID
    if ( sub_str == "MVP" ) sub_vec <- mvp_211_tib$Probe_ID
    
    for ( gts_min in opt$gts_mins ) {
      # gts_min <- 20
      gts_str <- paste0("gts",gts_min)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                           SNP Analysis:: EPICv2
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      grp_key <- paste( sub_str,"EPICv2", sep="_" )
      snp_key <- paste( grp_key, gts_str, sep="_"  )
      if ( p1 ) cat(glue::glue("{pmssg} Fingerprint Params: {snp_key}...{RET}"))
      
      # NOTE: Setting ssh_tib to null fixes the Sample_Name issue above, which
      #  should be fixed manually...
      # [TBD]: Fix Sample Name above in sample sheet...
      snvs_raw_gtc_dfs <- NULL
      snvs_raw_gtc_dfs <- analyze_snvs( 
        vcfs = snvs_raw_vcfs,
        run_tag = grp_key,
        # ssh_tib = NULL,
        ssh_tib = sel_ssh_tib,
        # ssh_tib = sel_ssh_tib %>% dplyr::mutate( Sample_Name = Sample_Index ),
        # ssh_tib = sel_ssh_tib %>% dplyr::mutate( Sample_Name = Sample_Name %>% stringr::str_remove("_DNA_") %>% stringr::str_remove("t") ),
        sub_vec = sub_vec,
        gts_min = gts_min, 
        fval = 5, 
        jcnt = 0, 
        uval = TRUE,
        outDir = file.path( opt$out_path ),
        write_out = FALSE, 
        plot_heat = TRUE
      )
      
      #
      # [TBD]: Try gts_min == 1, or some other numbers
      # [TBD]: Cluster on hit_mat
      #
      # [TBD]: Pre-Post 1-7: Compare
      #
      if ( FALSE ) {
        snvs_raw_gtc_tib <- NULL
        snvs_raw_gtc_tib <- snvs_raw_gtc_dfs %>% 
          as.data.frame() %>%
          tibble::rownames_to_column( var = "Probe_ID" ) %>%
          tibble::as_tibble()
        
        r2_tab <- NULL
        r2_tab <- snvs_raw_gtc_dfs$hit_mat %>% 
          cor( use = "pairwise.complete.obs", method = "pearson" ) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Sample_Name_A" ) %>% 
          tidyr::pivot_longer( cols = c(-Sample_Name_A), 
                               names_to = "Sample_Name_B", 
                               values_to = "r2" ) %>%
          tidyr::separate( Sample_Name_A, into=c("Sample_Group_A", "Sample_Class_A", "Sample_Index_A", "Sample_PPP_A"), sep="_", remove = FALSE ) %>%
          tidyr::separate( Sample_Name_B, into=c("Sample_Group_B", "Sample_Class_B", "Sample_Index_B", "Sample_PPP_B"), sep="_", remove = FALSE )
        
        r2_max_tib <- NULL
        r2_max_tib <- r2_tab %>%
          dplyr::filter( Sample_Name_A != Sample_Name_B ) %>%
          dplyr::arrange( -r2 ) %>% 
          dplyr::distinct( Sample_Name_A, .keep_all = TRUE )
        
        r2_max_tib %>% dplyr::filter( Sample_Index_A == Sample_Index_B )
        r2_max_tib %>% dplyr::filter( Sample_Index_A != Sample_Index_B )
        
        r2_tab %>%
          dplyr::filter( Sample_Name_A != Sample_Name_B ) %>%
          dplyr::arrange( Sample_Index_A,-r2 ) %>%
          # head() %>% 
          # dplyr::filter( Sample_Index_A != "10" ) %>%
          # dplyr::filter( Sample_Index_A != "17" ) %>%
          as.data.frame() %>%
          tibble::as_tibble()
      }
      
      # break
    }
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
