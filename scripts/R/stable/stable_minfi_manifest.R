
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Script for Analyzing Minfi Manifest Generation
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Minfi EPICv2 Packages:: https://support.bioconductor.org/p/9150344/#9154928
library(BiocManager)
# install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
# install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")

# Minfi Packages from old scripts::
library(limma)
library(minfi)
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# BiocManager::install("IlluminaHumanMethylation450kmanifest")
library(IlluminaHumanMethylation450kmanifest)
# BiocManager::install("IlluminaHumanMethylationEPICmanifest")
library(IlluminaHumanMethylationEPICmanifest)

library(RColorBrewer)
# library(missMethyl)
library(matrixStats)
# library(Gviz)
# library(DMRcate)
library(stringr)
library(hexbin)
library(plotly)
library(png)
library(gridExtra)
library(grid)

# Options, Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )

source("/Users/bbarnes/mylocal/Documents/tools/IlluminaHumanMethylationEPICv2manifest/inst/scripts/read.manifest.EPIC.R")

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
par$prgm_tag <- 'stable_minfi_manifest'
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
#                 Pre-processing:: Load Minfi EPICv2 Manifest
#
# Old Code: Projects.old/methylation/GravityLens/scripts/minfi.R
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ssh_dir <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted" )
epi_ssh_dir <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/idats" )

epi_ssh_dfs <- NULL
epi_ssh_dfs <- read.metharray.sheet(epi_ssh_dir, pattern = "csv$", 
                                    ignore.case = TRUE, recursive = FALSE, 
                                    verbose = TRUE) %>% head(n=2) %>%
  dplyr::mutate( Sentrix_ID = Slide, 
                 Sentrix_Position = Array )

epi_ssh_dfs$Basename <- file.path( epi_ssh_dir, epi_ssh_dfs$Sentrix_ID, 
                                   paste(epi_ssh_dfs$Sentrix_ID, 
                                         epi_ssh_dfs$Sentrix_Position, sep="_") )

epi_raw_rgSet <- NULL
epi_raw_rgSet <- read.metharray.exp( targets = epi_ssh_dfs, force = TRUE, extended = TRUE )
# rgSet

# Reproduce EPICv2
# Build MSA

# EPICv2 RDS investigation::
epicv2_str <- NULL
epicv2_rda <- file.path( opt$top_path, "tools/IlluminaHumanMethylationEPICv2manifest/data/IlluminaHumanMethylationEPICv2manifest.rda" )
epicv2_str <- load( file = epicv2_rda )

# Update Annotation
epi_ann_rgSet <- NULL
epi_ann_rgSet <- epi_raw_rgSet
epi_ann_rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
# rgSet

# annotation(epi_ann_rgSet)["array"] = "IlluminaHumanMethylationEPICv2"
# preprocessRaw(epi_ann_rgSet)

cat("---------- ---------- ---------- ---------- ---------- ----------\n\n");
cat("\t[Read]: Reading detP and RS# SNPs:\n\n");
cat("---------- ---------- ---------- ---------- ---------- ----------\n\n");

# epi_raw_detP <- detectionP(epi_raw_rgSet);
epi_ann_detP <- NULL
epi_ann_detP <- detectionP(epi_ann_rgSet);
cat("\t[Read]: EPICv2 detP:\n\n");

epi_ann_snps <- NULL
epi_ann_snps <- getSnpBeta(epi_ann_rgSet);
cat("\t[Read]: EPICv2 snps:\n\n");

epi_ann_phenoData <- NULL
epi_ann_phenoData <- pData(epi_ann_rgSet);
cat("\t[Read]: EPICv2 phenoData:\n\n");

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             MSA Development::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa_man_tib <- NULL
msa_man_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MSA-48v1-0_20102838_A1.csv.gz" )
msa_man_tib <- load_genome_studio_manifest( file = msa_man_csv, ret_data = TRUE )

msa_dat_dir <- NULL
msa_dat_dir <- file.path( opt$top_path, "data/idats/idats_MSA_v1-48x1-04072023_MSAEX01V" )

msa_ssh_dfs <- NULL
msa_ssh_dfs <- read.metharray.sheet( msa_dat_dir, pattern = "csv$", 
                                     ignore.case = TRUE, recursive = FALSE, 
                                     verbose = TRUE) %>% head(n=2) %>%
  dplyr::mutate( Sentrix_ID = Sample_ID %>% stringr::str_remove("_.*$"), 
                 Sentrix_Position = Sample_ID %>% stringr::str_remove("^.*_"),
                 Slide = Sentrix_ID )

msa_ssh_dfs$Basename <- file.path( msa_dat_dir, msa_ssh_dfs$Sentrix_ID, 
                                   paste(msa_ssh_dfs$Sentrix_ID, 
                                         msa_ssh_dfs$Sentrix_Position, sep="_") )

msa_raw_rgSet <- NULL
msa_raw_rgSet <- read.metharray.exp( targets = msa_ssh_dfs, force = TRUE, extended = TRUE )

msa_man_dat <- NULL
msa_man_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MSA-48v1-0_20102838_A1.csv" )
# msa_man_rda <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MSA-48v1-0_20102838_A1.rda" )
msa_man_dat <- read.manifest.EPIC( file = msa_man_csv )

IlluminaHumanMethylationMSAmanifest <- NULL
IlluminaHumanMethylationMSAmanifest <- IlluminaMethylationManifest(
  TypeI = msa_man_dat$manifestList$TypeI,
  TypeII = msa_man_dat$manifestList$TypeII,
  TypeControl = msa_man_dat$manifestList$TypeControl,
  TypeSnpI = msa_man_dat$manifestList$TypeSnpI,
  TypeSnpII = msa_man_dat$manifestList$TypeSnpII,
  annotation = "IlluminaHumanMethylationMSAmanifest")
# annotation = "IlluminaHumanMethylationMSA")

msa_man_rda <- file.path( opt$top_path, "data/manifests/methylation/minfi/IlluminaHumanMethylationMSAmanifest.rda" )
stopifnot(validObject(IlluminaHumanMethylationMSAmanifest))
save(IlluminaHumanMethylationMSAmanifest, compress = "xz",
     file = msa_man_rda )

msa_man_str <- NULL
msa_man_str <- load( file = msa_man_rda )

# Update Annotation
msa_ann_rgSet <- NULL
msa_ann_rgSet <- msa_raw_rgSet
msa_ann_rgSet@annotation <- c(array = "IlluminaHumanMethylationMSA", annotation = "20102838_A1.hg38")

# preprocessRaw(msa_ann_rgSet)

#
# MSA: Analysis
#
msa_ann_detP <- NULL
msa_ann_detP <- detectionP(msa_ann_rgSet);
cat("\t[Read]: MSA detP:\n\n");

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
