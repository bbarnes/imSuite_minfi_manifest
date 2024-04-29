
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
par$prgm_tag <- 'stable_trifecta_quiz'
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
# par$run_name <- "MSAv10"
par$run_name <- "TRIv01"

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
#             Pre-processing:: Load Study Files
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

dat_dir <- file.path( opt$top_path, "Projects.new/EPIC_v2/data/trifecta/quiz" )
tri_csv <- file.path( dat_dir, "EPICv2_EAtest_trifecta.csv.gz" )
gst_tsv <- file.path( dat_dir, "EPICv2_EAtest_MethylationGenotyping.tsv" )

# A tibble: 565 × 19
tri_tib <- NULL
tri_tib <- load_genome_studio_manifest( 
  file = tri_csv, load_clean = TRUE, load_controls = FALSE ) %>% 
  dplyr::arrange( IlmnID )

# A tibble: 955 × 8
gst_tib <- NULL
gst_tib <- readr::read_tsv( file = gst_tsv, show_col_types = FALSE ) %>% 
  dplyr::arrange( Trifecta_rsID )

#
# Simple Stats::
#

tri_tib %>% head(n=3) %>% as.data.frame()
gst_tib %>% head(n=3) %>% as.data.frame()

tri_tib %>% dplyr::filter( IlmnID == "rs1000769227" ) %>% as.data.frame()
gst_tib %>% dplyr::filter( Trifecta_rsID == "rs1000769227" ) %>% as.data.frame()

#
# Disjoint::
#

gst_tib %>% dplyr::filter(  Trifecta_rsID %in% tri_tib$Name )
gst_tib %>% dplyr::filter( !Trifecta_rsID %in% tri_tib$Name )


#
# Intersections::
#

gst_tib %>% dplyr::anti_join(  tri_tib, by=c("Trifecta_rsID"="Name") ) %>% base::nrow()
gst_tib %>% dplyr::anti_join(  tri_tib, by=c("Trifecta_rsID"="Name") ) %>% head(n=3) %>% as.data.frame()
gst_tib %>% dplyr::inner_join( tri_tib, by=c("Trifecta_rsID"="Name") ) %>% base::nrow()
gst_tib %>% dplyr::inner_join( tri_tib, by=c("Trifecta_rsID"="Name") ) %>% head(n=3) %>% as.data.frame()

tri_tib %>% dplyr::anti_join(  gst_tib, by=c("Name"="Trifecta_rsID") ) %>% head(n=3) %>% as.data.frame()
tri_tib %>% dplyr::inner_join( gst_tib, by=c("Name"="Trifecta_rsID") ) %>% head(n=3) %>% as.data.frame()



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
