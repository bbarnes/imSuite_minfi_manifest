
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
suppressWarnings(suppressPackageStartupMessages( 
  base::require("sesame",  quietly = TRUE) ) )

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
par$prgm_tag <- 'stable_MSA_VCF_analysis'
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
# par$version <- 1

# par$run_name <- "EPICv1"
# par$run_name <- "Embarkv1"
# par$run_name <- "FAILv1"
# par$run_name <- "COREv1"
#
# par$run_name <- "MSAv03"
# par$run_name <- "EPICv2"
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
#                    Pre-processing:: Load Sample Sheet
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

dat_dir <- file.path( opt$top_path, "Projects.new/MSA/data/MSA-Alpha-CHOP/AlphaData-07122023_PreDVT_SGrun" )

# [Note]: Everything is Semi-Auto
# [Note]: Only 50 and 250 ng inputs
# [Note]: Blood only comes in 50 ng
# [Note]: FFPE only comes in 250 ng
# [Note]: Everything else comes in 50 and 250 ng
src_ssh_tib <- NULL
src_ssh_csv <- file.path( dat_dir, "Samplesheet_48x1_preDVT_072023-SMG-Verbose.csv.gz" )
src_ssh_tib <- readr::read_csv( file = src_ssh_csv, show_col_types = FALSE, skip = 7 ) %>% 
  dplyr::rename( DNA_Input = 'DNA Input' ) %>% 
  dplyr::mutate( DNA_Input = DNA_Input %>% stringr::str_remove("ng") %>% as.integer(),
                 Sample_Group = Sample_Group %>% stringr::str_remove_all(" ") ) %>%
  tidyr::unite( Sentrix_Name, Sentrix_ID,Sentrix_Position, sep = "_" )

src_ssh_sum <- NULL
src_ssh_sum <- src_ssh_tib %>% 
  dplyr::group_by( Sample_Group,DNA_Input ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

sel_ssh_tib <- NULL
sel_ssh_tib <- src_ssh_tib %>%
  dplyr::filter( Sample_Group == "Blood" ) %>%
  dplyr::arrange( Sample ) %>%
  head( n=6 )

exp_idat_lst <- sel_ssh_tib %>% split( .$Sample_Name )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Pre-processing:: Load Manifest
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

exp_man_tib <- NULL
exp_man_csv <- file.path( opt$top_path, "Projects.new/MSA/manifests/MSA-48v1-0-Post-PQC_2B.v1.0.rccp.csv.gz" )
exp_man_tib <- readr::read_csv( file = exp_man_csv, show_col_types = FALSE )

exp_ctl_tib <- NULL
exp_ctl_tib <- exp_man_tib %>% dplyr::filter( Probe_ID %>% stringr::str_starts("ctl") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Processing:: Idats to SDF to VCF
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# file_list()
# Add file paths...
tmp_idat_lst <- NULL
tmp_idat_lst <- file_list( path = dat_dir, suffix = "_Red.idat.gz", 
                           pattern = "_Red.idat.gz$", recursive = TRUE,
                           vb=vb, vt=vt, tc=tc, tt=tt )

all_idat_lst <- NULL
all_idat_lst <- sesame::searchIDATprefixes( dir.name = dat_dir ) # %>% head()

# prefix_to_sdf()

# [TBD]: Update fields
# [TBD]: Add formatVCF() 

platform <- NULL
platform <- "MSA"

exp_out_dir <- file.path( opt$out_path )

if ( FALSE ) {
  cur_sdfs <- NULL
  if ( opt$parallel ) {
    cur_sdfs <- mclapply( all_idat_lst,
                          prefix_to_sdf, 
                          platform = platform,
                          manifest = exp_man_tib, 
                          controls = exp_ctl_tib, 
                          out_dir  = exp_out_dir,
                          run_tag  = opt$run_name,
                          reload   = opt$reload, 
                          reload_min = 10,
                          parallel = opt$parallel, 
                          vb=vb, vt=vt, tc=0, tt = tt )
  } else {
    cur_sdfs <- lapply( all_idat_lst,
                        prefix_to_sdf, 
                        platform = platform,
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



# formatVCF()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                        Processing:: SNV VCFs
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# [NOTE]: Below is example code that needs to be updated for this project
#

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
  c
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
          dplyr::filter( Sample_Index_A != "10" ) %>%
          dplyr::filter( Sample_Index_A != "17" ) %>%
          as.data.frame()
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
