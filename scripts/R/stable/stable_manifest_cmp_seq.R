
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
par$prgm_tag <- 'stable_manifest_cmp_seq'
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

#
# 1. Load his/rem (historic/removal) [ 27k, 450k, EPICv1, MM280k ]
# 2. Load ref (reference) [ Horvath40k ]
# 3. Remove his vs. ref
# 4. Loop through can (candidates) [ EPICv2, MSA, Embark ]
#
# Add Control Seqs?
# +/- Mouse...

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Pre-processing:: Load & Expand Manifests
#                            Historic/Removal
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# [TBD]: Remove Mouse Array
# [TBD]: Check cg# overlap

rem_list <- list()
rem_list[["HM027"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz" )
rem_list[["HM450"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz" )
rem_list[["HM850"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz" )
# rem_list[["MM289"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MouseMethylation-12v1-0_A2.csv.gz" )

rem_man_tib <- NULL
rem_man_tib <- rem_list %>% # head(n=1) %>% tail(n=1) %>%
  lapply(function(x) {
  load_genome_studio_manifest( 
    file = x, load_clean = TRUE, load_controls = FALSE, 
    cols_convert = FALSE, write_clean = FALSE, overwrite = FALSE ) %>%
    dplyr::rename( Probe_ID = IlmnID ) %>%
      man_to_prb_stack( prb_vec = c("AlleleA_ProbeSeq","AlleleB_ProbeSeq"),
                        cgn_vec = c("Probe_ID"),
                        # pos_vec = c("MAPINFO"), 
                        expand_prb = TRUE,
                        sort_cgn = TRUE,
                        add_fwd = TRUE, 
                        add_top = FALSE, 
                        add_tbs = FALSE, 
                        add_grp = TRUE, 
                        add_ord = TRUE,
                        uc = TRUE,
                        vb=vb,vt=vt,tc=tc )
  }) %>% dplyr::bind_rows( .id = "Manifest" )



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Pre-processing:: Load & Expand Manifests
#                                Reference
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ref_man_tib <- NULL
ref_man_csv <- file.path( opt$top_path, "data/manifests/methylation/csv/HorvathMammalMethylChip40/internal/AllMammalMethyl_B1.csv.gz" )
ref_man_tib <- load_genome_studio_manifest( 
  file = ref_man_csv, load_clean = TRUE, load_controls = FALSE, 
  cols_convert = FALSE, write_clean = FALSE, overwrite = FALSE ) %>%
  dplyr::rename( Probe_ID = IlmnID )

# A tibble: 204,843 × 4
ref_exp_tib <- NULL
ref_exp_tib <- man_to_prb_stack( tib = ref_man_tib,
                                 prb_vec = c("AlleleA_ProbeSeq","AlleleB_ProbeSeq"),
                                 cgn_vec = c("Probe_ID"),
                                 # pos_vec = c("MAPINFO"), 
                                 expand_prb = TRUE,
                                 sort_cgn = TRUE,
                                 add_fwd = TRUE, 
                                 add_top = FALSE, 
                                 add_tbs = FALSE, 
                                 add_grp = TRUE, 
                                 add_ord = TRUE,
                                 uc = TRUE,
                                 vb=vb,vt=vt,tc=tc )

# ref_exp_tib %>% dplyr::distinct( Probe_ID )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Pre-processing:: Remove Historic Probes
#                             New Reference
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ref_new_tib <- NULL
ref_new_tib <- ref_exp_tib %>% dplyr::anti_join( rem_man_tib, by=c("Exp_Sequence") )

# ref_new_tib %>% dplyr::distinct( Probe_ID )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Pre-processing:: Load & Expand Manifests
#                                Candidate(s)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

can_list <- list()
can_list[["EM057"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Embark_Methylation_1-1-1_2B.csv.gz" )
can_list[["MS280"]] <- file.path( opt$top_path, "Projects.new/MSA/data/MSA-Alpha-CHOP/MSA-48v1-0-Post-PQC_3E-DoNotDistribute-Zhou.csv.gz" )
can_list[["MS290"]] <- file.path( opt$top_path, "Projects.new/MSA/manifests/MSA-48v1-0-Post-PQC_2B.csv.gz" )
can_list[["HM960"]] <- file.path( opt$top_path, "data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1.csv.gz" )

# Re-run Removal Lists as QC
can_list[["HM027"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/HumanMethylation27_270596_v.1.2.csv.gz" )
can_list[["HM450"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/HumanMethylation450_15017482_v.1.2.csv.gz" )
can_list[["HM850"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MethylationEPIC_v-1-0_B2.csv.gz" )
can_list[["MM289"]] <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/MouseMethylation-12v1-0_A2.csv.gz" )

can_list[["MM289"]] <- file.path( opt$top_path, "data/manifests/methylation/csv/HorvathMammalMethylChip40/internal/AllMammalMethyl_B1.csv.gz" )


int_all_tib <- NULL
int_all_sum <- NULL

for ( can_man_str in names(can_list) ) {
  can_man_csv <- can_list[[can_man_str]]
  
  if ( p1 ) cat(glue::glue("{pmssg} Processing manifest({can_man_str}): {can_man_csv}.{RET}"))
  
  can_man_tib <- NULL
  can_man_tib <- load_genome_studio_manifest( 
    file = can_man_csv, load_clean = TRUE, load_controls = FALSE, 
    cols_convert = FALSE, write_clean = FALSE, overwrite = FALSE ) %>%
    dplyr::rename( Probe_ID = IlmnID )
  
  can_man_cgn_unq_cnt <- NULL
  can_man_cgn_unq_cnt <- can_man_tib %>% dplyr::distinct( Probe_ID ) %>% base::nrow()
  
  # A tibble: 204,843 × 4
  can_prb_tib <- NULL
  can_prb_tib <- man_to_prb_stack( tib = can_man_tib,
                                   prb_vec = c("AlleleA_ProbeSeq","AlleleB_ProbeSeq"),
                                   cgn_vec = c("Probe_ID"),
                                   # pos_vec = c("MAPINFO"), 
                                   expand_prb = TRUE,
                                   sort_cgn = TRUE,
                                   add_fwd = TRUE, 
                                   add_top = FALSE, 
                                   add_tbs = FALSE, 
                                   add_grp = TRUE, 
                                   add_ord = TRUE,
                                   uc = TRUE,
                                   vb=vb,vt=vt,tc=tc )
  
  can_prb_cgn_unq_cnt <- NULL
  can_prb_cgn_unq_cnt <- can_prb_tib %>% dplyr::distinct( Probe_ID ) %>% base::nrow()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Processing:: Intersect Manifests: By Probe
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  int_inn_tib <- NULL
  int_inn_tib <- ref_new_tib %>% 
    dplyr::inner_join( can_prb_tib, by=c("Exp_Sequence"), 
                       suffix = c("_ref","_can"), 
                       relationship = "many-to-many" ) %>%
    dplyr::mutate( Int_Class = "Inn",
                   Can_Sample = can_man_str )

  int_ref_tib <- NULL
  int_ref_tib <- ref_new_tib %>% 
    dplyr::left_join( can_prb_tib, by=c("Exp_Sequence"), 
                       suffix = c("_ref","_can"), 
                       relationship = "many-to-many" ) %>%
    dplyr::mutate( Int_Class = "Ref",
                   Can_Sample = can_man_str )

  int_can_tib <- NULL
  int_can_tib <- ref_new_tib %>% 
    dplyr::right_join( can_prb_tib, by=c("Exp_Sequence"), 
                      suffix = c("_ref","_can"), 
                      relationship = "many-to-many" ) %>%
    dplyr::mutate( Int_Class = "Can",
                   Can_Sample = can_man_str )
  
  int_all_tib <- int_all_tib %>%
    dplyr::bind_rows( int_inn_tib ) %>% 
    dplyr::bind_rows( int_ref_tib ) %>%
    dplyr::bind_rows( int_can_tib )
  
  if ( p1 ) cat(glue::glue("{pmssg} can_man_cgn_unq_cnt={can_man_cgn_unq_cnt}.{RET}"))
  if ( p1 ) cat(glue::glue("{pmssg} can_prb_cgn_unq_cnt={can_prb_cgn_unq_cnt}.{RET2}"))

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Processing:: Summarize Intersections: By Probe
  #                                   Exp
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  int_inn_sum <- NULL
  int_inn_sum <- int_inn_tib %>% 
    dplyr::group_by( Probe_ID_ref ) %>% 
    dplyr::summarise( Tot_All_Cnt = n(), 
                      Mat_Ref_Cnt = sum ( !is.na(Fwd_Sequence_ref) ),
                      Mat_Can_Cnt = sum ( !is.na(Fwd_Sequence_can) ),
                      Mat_Ref_Per = round( (100*Mat_Ref_Cnt / Tot_All_Cnt),3 ),
                      .groups = "drop" ) %>%
    dplyr::filter( Mat_Can_Cnt != 0 ) %>%
    dplyr::filter( Mat_Ref_Cnt != 0 ) %>% 
    dplyr::mutate( Probe_ID = Probe_ID_ref %>% stringr::str_remove("_[0-9]+$"), 
                   Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   SRD_Str = Probe_ID %>% stringr::str_remove("^.*_") ) %>% 
    tidyr::separate( SRD_Str, 
                     into=c("Strand_TB","Strand_CO","Infinium_Design","Rep_Idx"), 
                     sep = c(1,2,3,4), remove = TRUE, convert = TRUE )
  
  int_exp_cnt <- NULL
  int_exp_cnt <- int_inn_sum %>% dplyr::filter( Tot_All_Cnt != Mat_Ref_Cnt ) %>%
    base::nrow()
  
  if ( int_exp_cnt ) {
    stop( glue::glue("{perrs} Inner Intersection not equal to 0 != {int_exp_cnt}.{RET2}") )
    break
  }

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Processing:: Summarize Intersections: By Probe
  #                                   Ref
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  int_ref_sum <- NULL
  int_ref_sum <- int_ref_tib %>% 
    dplyr::group_by( Probe_ID_ref ) %>% 
    dplyr::summarise( Tot_All_Cnt = n(), 
                      Mat_Ref_Cnt = sum ( !is.na(Fwd_Sequence_ref) ),
                      Mat_Can_Cnt = sum ( !is.na(Fwd_Sequence_can) ),
                      Mat_Ref_Per = round( (100*Mat_Ref_Cnt / Tot_All_Cnt),3 ),
                      .groups = "drop" ) %>%
    dplyr::filter( Mat_Can_Cnt != 0 ) %>%
    dplyr::filter( Mat_Ref_Cnt != 0 ) %>% 
    dplyr::mutate( Probe_ID = Probe_ID_ref %>% stringr::str_remove("_[0-9]+$"), 
                   Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   SRD_Str = Probe_ID %>% stringr::str_remove("^.*_") ) %>% 
    tidyr::separate( SRD_Str, 
                     into=c("Strand_TB","Strand_CO","Infinium_Design","Rep_Idx"), 
                     sep = c(1,2,3,4), remove = TRUE, convert = TRUE )
  
  int_ref_cnt <- NULL
  int_ref_cnt <- int_ref_sum %>% dplyr::filter( Tot_All_Cnt != Mat_Ref_Cnt ) %>%
    base::nrow()
  
  if ( int_ref_cnt ) {
    stop( glue::glue("{perrs} Reference Intersection not equal to 0 != {int_ref_cnt}.{RET2}") )
    break
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Processing:: Summarize Intersections: By Probe
  #                                   Can
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  int_can_sum <- NULL
  int_can_sum <- int_can_tib %>% 
    dplyr::group_by( Probe_ID_can ) %>% 
    dplyr::summarise( Tot_All_Cnt = n(), 
                      Mat_Ref_Cnt = sum ( !is.na(Fwd_Sequence_ref) ),
                      Mat_Can_Cnt = sum ( !is.na(Fwd_Sequence_can) ),
                      Mat_Ref_Per = round( (100*Mat_Ref_Cnt / Tot_All_Cnt),3 ),
                      .groups = "drop" ) %>%
    dplyr::filter( Mat_Ref_Cnt != 0 ) %>% 
    dplyr::filter( Mat_Can_Cnt != 0 ) %>%
    dplyr::mutate( Probe_ID = Probe_ID_can %>% stringr::str_remove("_[0-9]+$"), 
                   Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   SRD_Str = Probe_ID %>% stringr::str_remove("^.*_") ) %>% 
    tidyr::separate( SRD_Str, 
                     into=c("Strand_TB","Strand_CO","Infinium_Design","Rep_Idx"), 
                     sep = c(1,2,3,4), remove = TRUE, convert = TRUE )
  
  # int_can_sum %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  int_can_cnt <- NULL
  int_can_cnt <- int_can_sum %>% dplyr::filter( Tot_All_Cnt != Mat_Can_Cnt ) %>%
    base::nrow()
  
  if ( int_can_cnt ) {
    stop( glue::glue("{perrs} Candidate Intersection not equal to 0 != {int_can_cnt}.{RET2}") )
    break
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Processing:: Join All Summaries
  #                               Exp,Ref,Can
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  int_all_sum <- int_all_sum %>% 
    dplyr::bind_rows( 
      int_inn_sum %>% 
        dplyr::group_by( Probe_Type ) %>% 
        dplyr::summarise( Tot_Cnt = n(), 
                          Avg_Ref = mean(Mat_Ref_Per), 
                          Med_Ref = median(Mat_Ref_Per), 
                          .groups = "drop" ) %>%
        dplyr::mutate( Can_Man = can_man_str,
                       Can_Cmp = "Exp") ) %>%
    dplyr::bind_rows( 
      int_ref_sum %>% 
        dplyr::group_by( Probe_Type ) %>% 
        dplyr::summarise( Tot_Cnt = n(), 
                          Avg_Ref = mean(Mat_Ref_Per), 
                          Med_Ref = median(Mat_Ref_Per), 
                          .groups = "drop" ) %>%
        dplyr::mutate( Can_Man = can_man_str,
                       Can_Cmp = "Ref") ) %>%
    dplyr::bind_rows( 
      int_can_sum %>% 
        dplyr::group_by( Probe_Type ) %>% 
        dplyr::summarise( Tot_Cnt = n(), 
                          Avg_Ref = mean(Mat_Ref_Per), 
                          Med_Ref = median(Mat_Ref_Per), 
                          .groups = "drop" ) %>%
        dplyr::mutate( Can_Man = can_man_str,
                       Can_Cmp = "Can") )
  

  # if ( can_man_str == "MS280" ) break
}

# A tibble: 16 × 6
# Probe_Type Tot_Cnt Avg_Ref Med_Ref Can_Man Can_Cmp
# <chr>        <int>   <dbl>   <dbl> <chr>   <chr>  
#  1 cg           39302   100       100 EM057   Exp    
#  2 cg           39302   100       100 EM057   Ref    
#  3 cg           39552    96.8     100 EM057   Can    
#  4 cg             212   100       100 MS280   Exp    
#  5 cg             212   100       100 MS280   Ref    
#  6 cg             194    95.6     100 MS280   Can    
#  7 cg             230   100       100 MS290   Exp    
#  8 cg             230   100       100 MS290   Ref    
#  9 cg             207    95.9     100 MS290   Can    
# 10 cg             346   100       100 HM960   Exp    
# 11 cg             346   100       100 HM960   Ref    
# 12 cg             342    94.4     100 HM960   Can    
# 13 nv               1   100       100 HM960   Can    
# 14 cg             817   100       100 MM289   Exp    
# 15 cg             817   100       100 MM289   Ref    
# 16 cg             774    91.8     100 MM289   Can    

# MS290
#
# > int_prb_sum == # A tibble: 3,452 × 3
# > int_prb_sum %>% dplyr::filter( Mat_Cnt != 0 ) == # A tibble: 3,452 × 3
# > int_ref_sum %>% dplyr::filter( Mat_Cnt != Tot_Cnt ) == # A tibble: 0 × 3
# > > int_prb_sum %>% dplyr::filter( Mat_Cnt != 0 ) %>% dplyr::filter( Mat_Cnt != Tot_Cnt ) == # A tibble: 0 × 3
#
# > int_ref_sum == # A tibble: 42,585 × 3
# > int_ref_sum %>% dplyr::filter( Mat_Cnt != 0 ) == # A tibble: 3,452 × 3
# > int_ref_sum %>% dplyr::filter( Mat_Cnt != Tot_Cnt ) # A tibble: 40,174 × 3
# > int_ref_sum %>% dplyr::filter( Mat_Cnt != 0 ) %>% dplyr::filter( Mat_Cnt != Tot_Cnt ) == # A tibble: 1,041 × 3
#
# > int_can_sum == # A tibble: 3,453 × 3
# > int_can_sum %>% dplyr::filter( Mat_Cnt != 0 ) == # A tibble: 3,453 × 3
# > int_can_sum %>% dplyr::filter( Mat_Cnt != Tot_Cnt ) # A tibble: 0 × 3
# > int_can_sum %>% dplyr::filter( Mat_Cnt != 0 ) %>% dplyr::filter( Mat_Cnt != Tot_Cnt ) == # A tibble: 0 × 3
#

# A tibble: 10 × 6
# Probe_Type Tot_Cnt Avg_Ref Med_Ref Can_Man Can_Cmp
# <chr>        <int>   <dbl>   <dbl> <chr>   <chr>  
# 1 cg           39205   100       100 EM057   Ref    
# 2 cg           39205   100       100 EM057   Ref    
# 3 cg           39452    96.3     100 EM057   Can    
# 4 cg             217   100       100 MS290   Ref    
# 5 cg             217   100       100 MS290   Ref    
# 6 cg             199    95.7     100 MS290   Can    
# 7 cg             338   100       100 HM960   Ref    
# 8 cg             338   100       100 HM960   Ref    
# 9 cg             334    94.5     100 HM960   Can    
# 10 nv               1   100       100 HM960   Can    


if ( FALSE ) {
  
  # A tibble: 306,327 × 12
  int_all_sum <- NULL
  int_all_sum <- int_all_tab %>% 
    # dplyr::filter( Int_Class  == "Can" ) %>%
    # dplyr::filter( Can_Sample == "EM057" ) %>%
    # dplyr::filter( Can_Sample == "MS290" ) %>%
    # dplyr::filter( !is.na(Fwd_Sequence_can) ) %>%
    # dplyr::filter( !is.na(Fwd_Sequence_ref) ) %>%
    # head(n=10000) %>% 
    dplyr::mutate( Probe_ID = Probe_ID_can %>% stringr::str_remove("_[0-9]+$") ) %>%
    dplyr::group_by( Can_Sample,Int_Class,Probe_ID ) %>% 
    dplyr::summarise( Tot_Prb_Cnt = n(), 
                      Mat_Can_Cnt = sum( !is.na(Fwd_Sequence_can) ),
                      Mat_Ref_Cnt = sum( !is.na(Fwd_Sequence_ref) ),
                      .groups = "drop")
  
  int_all_sum2 <- NULL
  int_all_sum2 <- int_all_sum %>% 
    dplyr::filter( Mat_Ref_Cnt != 0 ) %>% 
    dplyr::mutate( Mat_Can_Per = round( 100*Mat_Can_Cnt/Tot_Prb_Cnt, 3),
                   Mat_Ref_Per = round( 100*Mat_Ref_Cnt/Tot_Prb_Cnt, 3) ) %>% 
    dplyr::group_by( Can_Sample,Int_Class,Probe_ID ) %>% 
    dplyr::summarise( 
      Avg_Can_Per = mean(Mat_Can_Per),
      Med_Can_Per = median(Mat_Can_Per),
      Avg_Ref_Per = mean(Mat_Ref_Per),
      Med_Ref_Per = median(Mat_Ref_Per),
      .groups = "drop")
  
  int_all_sum3 <- NULL
  int_all_sum3 <- int_all_sum2 %>% 
    dplyr::group_by( Can_Sample,Int_Class ) %>% 
    dplyr::summarise( Count=n(), 
                      Avg_Can_Per = mean( Avg_Can_Per ),
                      Med_Can_Per = median( Med_Can_Per ),
                      Avg_Ref_Per = mean( Avg_Ref_Per ),
                      Med_Ref_Per = median( Med_Ref_Per ),
                      
                      .groups = "drop" )
  

}

if ( FALSE ) {
  int_prb_sum %>% 
    dplyr::mutate( Probe_ID = Probe_ID_ref %>% stringr::str_remove("_[0-9]+$"), 
                   Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"), 
                   Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   Mat_Avg = round( 100 * Mat_Cnt / Tot_Cnt, 3 ) ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::reframe( Tot_Cnt = n(), 
                    Avg_Mat = round( Mat_Avg, 3 ) )
  
  int_prb_sum %>% 
    dplyr::mutate( Probe_ID = Probe_ID_ref %>% stringr::str_remove("_[0-9]+$"), 
                   Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"), 
                   Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   Mat_Avg = round( 100 * Mat_Cnt / Tot_Cnt, 3 ) ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Tot_Cnt = n(), 
                      Avg_Mat = mean( Mat_Avg ),
                      .groups = "drop" )
  
  
  int_prb_sum %>% 
    dplyr::mutate( Probe_ID = Probe_ID_ref %>% stringr::str_remove("_[0-9]+$"), 
                   Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"), 
                   Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   Mat_Avg = Mat_Cnt / Tot_Cnt ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Tot_Cnt = n(), 
                      Avg_Mat = mean( round( 100 * Mat_Avg, 3 ) ),
                      .groups = "drop" )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
