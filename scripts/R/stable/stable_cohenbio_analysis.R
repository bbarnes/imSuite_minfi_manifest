
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
par$prgm_tag <- 'stable_cohenbio_analysis'
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

ssh_list <- NULL
# ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/full/swifthoof_main" )
# ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin/cohenbio/full/swifthoof_main" )
# ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/lin" )

split_ssh_loading <- FALSE
split_ssh_loading <- TRUE

if ( !split_ssh_loading ) {
  #
  # Joined merged set of data...
  #
  ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/join" )
  ssh_list <- file_list( path    = ssh_path, 
                         suffix  = "_EPIC_A1_AutoSampleSheet.csv.gz", 
                         pattern = "_EPIC_A1_AutoSampleSheet.csv.gz$", 
                         recursive = TRUE )
  
  
  ssh_tib <- NULL
  ssh_tib <- ssh_list %>%
    lapply( readr::read_csv, show_col_types = FALSE ) %>%
    dplyr::bind_rows()
  
  ssh_sum0 <- NULL
  ssh_sum0 <- ssh_tib %>% 
    dplyr::group_by( AutoSample_R2_1_Key_2,AutoSample_dB_1_Key_2 ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
} else {
  #
  # Split the two datasets::
  #
  sam_vec <- NULL
  sam_vec <- c( "EPICv2_Priyam_all", "EPICv2_MVP", "cohenbio" )
  
  lam_vec <- NULL
  lam_vec <- c("EPICv2_0", "EPICv2_1", "EPICv2_Cohen")
  
  # dir_vec <- c(
  #   "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/join",
  #   "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/join",
  #   "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2.1/par/EPICv2_MVP/chip-207881890007"
  # )
  
  ssh_tib <- NULL
  for ( ii in c(1:length(lam_vec) ) ) {
    sam <- sam_vec[ii]
    lam <- lam_vec[ii]
    
    ssh_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/snps/docker-1.11.15.1.p.0.6.2/join", sam )
    ssh_list <- file_list( path    = ssh_path, 
                           suffix  = "_EPIC_A1_AutoSampleSheet.csv.gz", 
                           pattern = "_EPIC_A1_AutoSampleSheet.csv.gz$", 
                           recursive = TRUE )
    
    
    ssh_tib <- ssh_list %>%
      lapply( readr::read_csv, show_col_types = FALSE ) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate( Sample_Source = !!lam ) %>%
      dplyr::bind_rows( ssh_tib )
    
    ssh_sum0 <- NULL
    ssh_sum0 <- ssh_tib %>% 
      dplyr::group_by( AutoSample_R2_1_Key_2,AutoSample_dB_1_Key_2 ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
  }
}

ssh_tib <- ssh_tib %>% 
  dplyr::mutate( 
    Sample_Type = dplyr::case_when( 
      AutoSample_dB_1_Key_2 == "T99BZ" ~ "Control",
      TRUE ~ "Sample" )
  )

sam_sum0 <- NULL
sam_sum0 <- ssh_tib %>% 
  dplyr::group_by( Sample_Source,AutoSample_dB_1_Key_2 ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

sam_sum1 <- NULL
sam_sum1 <- ssh_tib %>% 
  dplyr::group_by( Sample_Source,Sample_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Simple Plotting::
#      cg_pvals_PnegEcdf_pass_perc_1 vs. cg_pvals_pOOBAH_pass_perc_1
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_plot1 <- TRUE
run_plot1 <- FALSE

if ( run_plot1 ) {
  # Groups:
  # AutoSample_R2_1_Key_2
  # AutoSample_dB_1_Key_2
  
  ssh_pval_pnt_ggg <- NULL
  ssh_pval_pnt_ggg <- ssh_tib %>%
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    dplyr::group_by( AutoSample_R2_1_Key_2 ) %>%
    ggplot2::ggplot( aes( x=cg_pvals_PnegEcdf_pass_perc_1, 
                          y=cg_pvals_pOOBAH_pass_perc_1, color=AutoSample_R2_1_Key_2) ) +
    ggplot2::geom_point()
  # ssh_pval_pnt_ggg
  
  ssh_poob_den_ggg <- NULL
  ssh_poob_den_ggg <- ssh_tib %>%
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    dplyr::group_by( AutoSample_dB_1_Key_2 ) %>%
    ggplot2::ggplot( aes( x=cg_pvals_pOOBAH_pass_perc_1,
                          color=AutoSample_dB_1_Key_2) ) +
    ggplot2::geom_density()
  # ssh_poob_den_ggg
  
  ssh_poob_den_box1_ggg <- NULL
  ssh_poob_den_box1_ggg <- ssh_tib %>% 
    dplyr::filter( cg_pvals_PnegEcdf_pass_perc_1 > 96 ) %>%
    # dplyr::filter( Intervention != "Advanced" ) %>%
    tidyr::pivot_longer( cols = c(cg_pvals_PnegEcdf_pass_perc_1,cg_pvals_pOOBAH_pass_perc_1), 
                         names_to = "Pval_Type", values_to = "Pval_Pass_Percent" ) %>%
    dplyr::group_by( AutoSample_dB_1_Key_2 ) %>%
    ggplot2::ggplot( aes( x=Pval_Pass_Percent,
                          color=AutoSample_dB_1_Key_2) ) +
    ggplot2::geom_density() +
    ggplot2::facet_grid( # rows = vars(Sex_1),
                         cols = vars(Pval_Type), scales = "free" )
  # ssh_poob_den_box1_ggg
  
}

#
# BS conversion Analysis::
#
bs_den_tab <- NULL
bs_den_tab <- ssh_tib %>% dplyr::select( 
  Sample_Source,Sample_Type,
  AutoSample_dB_1_Key_2,
  GCT_1, GCT_2, 
  BISULFITE_CONVERSION_I_2_sig_M_mean_2, 
  BISULFITE_CONVERSION_II_2_sig_M_mean_2 ) %>% 
  tidyr::pivot_longer( cols=c( 
    # AutoSample_dB_1_Key_2,
    GCT_1, GCT_2,
    BISULFITE_CONVERSION_I_2_sig_M_mean_2, 
    BISULFITE_CONVERSION_II_2_sig_M_mean_2), 
    names_to = "name", values_to = "value" )

bs_den_ggg0 <- NULL
bs_den_ggg0 <- bs_den_tab %>% 
  ggplot2::ggplot( aes(x=value, group=name, color=name) ) + 
  ggplot2::geom_density( alpha=0.2 ) +
  ggplot2::facet_grid( 
    rows = vars( Sample_Source ),
    cols = vars( AutoSample_dB_1_Key_2 ) )
# bs_den_ggg0

# GCT Score::
bs_den_ggg1 <- NULL
bs_den_ggg1 <- bs_den_tab %>% 
  dplyr::filter( name %>% stringr::str_starts("GCT") ) %>%
  ggplot2::ggplot( aes(x=value ) ) + 
  ggplot2::geom_density( alpha=0.2 ) +
  ggplot2::facet_grid( 
    rows = vars( Sample_Source,name ),
    cols = vars( AutoSample_dB_1_Key_2 ) )
# bs_den_ggg1

# BS Controls
bs_den_ggg2 <- NULL
bs_den_ggg2 <- bs_den_tab %>% 
  dplyr::filter( !name %>% stringr::str_starts("GCT") ) %>%
  ggplot2::ggplot( aes(x=value ) ) + 
  ggplot2::geom_density( alpha=0.2 ) +
  ggplot2::facet_grid( 
    rows = vars( Sample_Source,name ),
    cols = vars( AutoSample_dB_1_Key_2 ) )
# bs_den_ggg2

#
# Scatter plots
#
bs_den_tab2 <- NULL
bs_den_tab2 <- ssh_tib %>% 
  dplyr::mutate(
    cg_pvals_pOOBAH_mean_2 = ( cg_1_pvals_pOOBAH_mean_2 + cg_2_pvals_pOOBAH_mean_2 ) / 2,
    cg_pvals_pOOBAH_mean_2_int = (cg_pvals_pOOBAH_mean_2 * 100) %>% as.integer()
  ) %>%
  dplyr::select( 
    Sample_Source,Sample_Type,
    AutoSample_dB_1_Key_2,
    cg_pvals_pOOBAH_mean_2,
    cg_pvals_pOOBAH_mean_2_int,
    GCT_1, GCT_2, 
    BISULFITE_CONVERSION_I_2_sig_M_mean_2, 
    BISULFITE_CONVERSION_II_2_sig_M_mean_2 ) %>% 
  tidyr::pivot_longer( cols=c( 
    # AutoSample_dB_1_Key_2,
    # GCT_1, GCT_2,
    BISULFITE_CONVERSION_I_2_sig_M_mean_2, 
    BISULFITE_CONVERSION_II_2_sig_M_mean_2), 
    names_to = "name", values_to = "value" )

bs_den_ggg3 <- NULL
bs_den_ggg3 <- bs_den_tab2 %>% 
  # dplyr::filter( name %>% stringr::str_starts("GCT") ) %>%
  ggplot2::ggplot( aes(x=value, y=GCT_2, color=name, group=name ) ) + 
  ggplot2::geom_point() +
  ggplot2::facet_grid( 
    rows = vars( Sample_Source,Sample_Type ),
    # cols = vars( Sample_Type )
  )
# bs_den_ggg3

bs_den_pdf4 <- file.path( opt$out_path, "bs_controls_vs_gct_score.pdf" )
bs_den_ggg4 <- NULL
bs_den_ggg4 <- bs_den_tab2 %>% 
  # dplyr::filter( name %>% stringr::str_starts("GCT") ) %>%
  # ggplot2::ggplot( aes(x=value, y=GCT_2, color=cg_pvals_pOOBAH_mean_2, group=cg_pvals_pOOBAH_mean_2 ) ) + 
  ggplot2::ggplot( aes(x=value, y=GCT_2, color=factor(cg_pvals_pOOBAH_mean_2_int), group=factor(cg_pvals_pOOBAH_mean_2_int) ) ) + 
  ggplot2::geom_point() +
  ggplot2::facet_grid( 
    rows = vars( Sample_Source,Sample_Type ),
    cols = vars( name )
  )
# bs_den_ggg4

ggplot2::ggsave( filename = bs_den_pdf4, plot = bs_den_ggg4, 
                 device = "pdf", width = 7, height = 7, dpi = 320 )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
