
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
#             Pre-processing:: Load Training Data: Alpha
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

local_dat_path <- file.path( opt$top_path, "scratch/stable_imSesameCpp_MSAv1_Alpha/MSAv10-NCBI-v0" )
train_sdf_path <- file.path( local_dat_path, "sdf" )

ssh_tib <- NULL
ssh_csv <- file.path( local_dat_path, "sample_sheets/Samplesheet_48x1_preDVT_072023-SMG-Verbose.formatted.csv.gz" )
ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE )

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

# dplyr::filter( Probe_Type=="cg")
# color=Min_Pval,
# cols = vars( Cnt_Rat, dBp_Dif, dBm_Dif, mPm_Dif, wSp_Dif)
# rows = vars( Grp,Outliers,Mask_Idx,Work_Str)
# scales = "free_y"

fin_sum_tab <- NULL
fin_sum_tab <- fin_sum_tib %>% 
  # dplyr::filter( Probe_Type == "cg" ) %>%
  dplyr::select( Outliers,Min_Pval,Mask_Idx,Work_Str,Probe_Type,Set,Grp, 
                 ppp, dBp, dBm, mPm, wSp
                 # Cnt_Rat, dBp_Dif, dBm_Dif, mPm_Dif, wSp_Dif
  ) %>%
  tidyr::pivot_longer( cols = c(ppp, dBp, dBm, mPm, wSp) )
                         # c(Cnt_Rat, dBp_Dif, dBm_Dif, mPm_Dif, wSp_Dif) )

# [TBD]:: Loop over Grp's

grp_vec <- unique(fin_sum_tab$Grp)

# [TBD]:: Reduce Min_Pval, subtract min value for each group...
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


fin_sum_tib %>% 
  dplyr::group_by( Outliers, Min_Pval, Mask_Idx, Work_Str, Probe_Type, Set, Grp) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )

# print(n=nrow(fin_sum_tib))

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
  
  
  #
  # OLD::
  #
  cur_dbs_sum %>% dplyr::select( Chop_Cut,Experiment,wSp_Avg ) %>% 
    tidyr::pivot_wider( names_from = c(Chop_Cut), values_from = c(wSp_Avg) ) %>% 
    dplyr::mutate( Dif = N - Y )
  
  all_dbs_sum %>% 
    dplyr::select( Outliers,Min_Pval,Mask_Idx,Work_Str,Chop_Cut,
                   Experiment,wSp_Avg ) %>% 
    tidyr::pivot_wider( 
      names_from = c(Outliers,Min_Pval,Mask_Idx,Work_Str,Chop_Cut), 
      values_from = c(wSp_Avg) ) %>% 
    dplyr::mutate( Dif = N - Y )
  
  # Use the following dplyr code to identify duplicates.
  all_dbs_sum  %>%
    dplyr::group_by(Experiment, Chop_Cut) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n > 1L)
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
