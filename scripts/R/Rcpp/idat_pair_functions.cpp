
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppCommon.h>

// My headers::
#include "idat_pair.h"

using namespace Rcpp;
using namespace idat_pair;

// Unclear why these don't allow as_tibble function to be imported, ignore for now...
Environment pkg_tidyverse = Environment::namespace_env("tidyverse");
Environment pkg_tibble = Environment::namespace_env("dplyr");
// Environment pkg_tibble = Environment::namespace_env("tibble");

Function asTibble("as_tibble");
Function typeConvert("type.convert");

/*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
 *
 *                    Rcpp to Cpp Local Conversion Function::
 * 
 * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */

void rcpp_to_cpp_vec_int( Rcpp::IntegerVector& vec_r,
                          std::vector<_i>& vec ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  for ( _I ii = 0; ii < s; ii++ ) vec[ii] = vec_r(ii);
}

void rcpp_to_cpp_vec_str( Rcpp::CharacterVector& vec_r,
                          std::vector<_S>& vec ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  for ( _I ii = 0; ii < s; ii++ ) vec[ii] = vec_r(ii);
}

void rcpp_to_cpp_vec_bool( Rcpp::LogicalVector& vec_r,
                           std::vector<_Y>& vec ) {
  _I s = vec_r.size();
  
  vec.clear();
  vec.resize( s );
  for ( _I ii = 0; ii < s; ii++ ) vec[ii] = vec_r(ii);
}

/* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
 * 
 * External Rcpp Functions::
 * 
 *                              Load Pair of Idats::
 *                         
 *                                   Goals::
 *
 *  1. Idat_Pair::
 *     - set_manifest()
 *     - set_masks()
 *        Sets Color Channel
 *        Sets SNPs:: Need to look up how formatVCF() works...
 *        { Base, DetSum,DetCut, SigSum,SigMax }
 *        
 *     - Locus Concordance 
 *        Foreach cgn with reps; foreach Mask calculate 
 *         { dB_avg,dB_sds, dB_med,dB_mad }
 *         
 *     - dyeBiasCorr()
 *     - dyeBiasCorrTypeINorm() ***
 *     - dyeBiasCorrMostBalanced()
 *     
 *     - detectionPnegEcdf()
 *     - detectionPoobEcdf()
 *     
 *     - Calculate Summary Stats::
 *     - sesameQC()
 *     - set_bacr()
 * 
 * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */

// [[Rcpp::export]]
Rcpp::DataFrame read_idat_pair_rcpp( const _S& prefix_path, 
                                     const _S& output_path,
                                     std::vector< _S >& workflow_vec,
                                     std::vector< _I >& pval_add_vec,
                                     std::vector< _I >& addU_man_vec,
                                     std::vector< _I >& addM_man_vec,
                                     
                                     std::vector< _S >& cgns_man_vec,
                                     std::vector< _S >& cols_man_vec,
                                     std::vector< _S >& keys_man_vec,
                                     std::vector< _S >& chrs_man_vec,
                                     std::vector< _S >& anns_man_vec,
                                     
                                     const _d min_pval = 0.05,
                                     const _d min_beta = 0.30,
                                     const _d max_beta = 0.70,
                                     
                                     const _d min_perO = 0.75,
                                     const _d min_perI = 0.05,
                                     
                                     const _Y read_bgz  = false,
                                     const _Y write_bgz = false,
                                     
                                     const _Y rm_pval_outliers = false,
                                     
                                     const _I return_df = 0,
                                     
                                     const _I vb = 0,
                                     const _I vt = 1,
                                     const _I tc = 0,
                                     const _S ft="read_idat_pair_rcpp" )
{
  const _S tb(tc, '\t');
  const _S _fo("["+ft+"]: "+tb);
  const _S _fe("["+ft+"]: ERROR: ");
  const _S _fw("["+ft+"]: "+tb+"Warning: ");
  
  const _Y p0 = vb > vt + 0;
  // const _Y p1 = vb > vt + 1;
  // const _Y p2 = vb > vt + 2;
  // const _Y p8 = vb > vt + 8;
  
  // Environment pkg = Environment::namespace_env("sesame");
  // Function SigDF_r("SigDF");
  // Function sesameQC_r("sesameQC");
  // Function sesameQC_calcStats_r("sesameQC_calcStats");
  
  _I pval_size = pval_add_vec.size();
  _I cgns_size = cgns_man_vec.size();
  
  if ( p0 ) Rcpp::Rcerr
  <<_fo<< "Statring...\n"
  <<_fo<<"\t  (vb,vt,tc): ("<<vb<<","<<vt<<","<<tc<<")\n"
  <<_fo<<"\t prefix_path: '"<<prefix_path<<"'\n"
  <<_fo<<"\t   pval_size: '"<<pval_size<<"'\n"
  <<_fo<<"\t   cgns_size: '"<<cgns_size<<"'\n"
  <<_fo<<"\t output_path: '"<<output_path<<"'\n"
  <<_fo<<"\t    min_perO: '"<<min_perO<<"'\n"
  <<_fo<<"\t    min_perI: '"<<min_perI<<"'\n"
  <<_fo<<"\t    min_pval: '"<<min_pval<<"'\n"
  <<_fo<<"\t    min_beta: '"<<min_beta<<"'\n"
  <<_fo<<"\t    max_beta: '"<<max_beta<<"'\n"
  <<_fo<<"\t    read_bgz: '"<<read_bgz<<"'\n"
  <<_fo<<"\t   write_bgz: '"<<write_bgz<<"'\n"
  <<_fo<<"\t   return_df: '"<<return_df<<"'\n"
  <<std::endl;
  
  Rcpp::DataFrame df;
  _Y success = true;
  
  const bool run_pair = true;
  // const bool run_pair = false;
  
  if ( run_pair ) {
    
    idat_pair::Idat_Pair idat_pair( prefix_path, output_path,
                                    pval_add_vec,
                                    addU_man_vec, addM_man_vec,
                                    cgns_man_vec, cols_man_vec,
                                    keys_man_vec, anns_man_vec,
                                    chrs_man_vec,
                                    min_pval, min_beta, max_beta,
                                    min_perO, min_perI,
                                    read_bgz, write_bgz,rm_pval_outliers,
                                    vb,vt+1,tc+1 );
    
    // if ( p0 ) Rcpp::Rcerr
    //   <<_fo<<"Idat_Pair:: to_string():\n"
    //   <<_fo<<"\n"<<idat_pair.to_string()<<"\n"
    //   <<std::endl;
    
    // Write Summary::
    //  - summarize_manifest()
    // Loop through Workflows
    //  - mutate_manifest()
    
    _S IDX_STR = "0";
    
    _S MASK_0_STR = "";
    _S MASK_N_STR = "";
    _S BETA_N_STR = "";
    _S PVAL_N_STR = "";
    _S DETP_N_STR = "";
    
    if ( return_df > 0 ) {
      // Build Manifest Data Frame::
      if ( idat_pair.get_Info_Vector( PROBE_ID ).size() != 0 )
        df.push_back( idat_pair.get_Info_Vector( PROBE_ID ), "Probe_ID" );
    }
    
    //
    // TBD:: Do we need this below if we're already going to do it for return_df > 3???
    // - Validate the values are the same and then remove code below if that's
    //.  the case... 
    //
    // if ( return_df > 1 ) {
    //   _S BETA_N_STR = BETA_STR+"_"+IDX_STR;
    //   if ( idat_pair.get_Pair_Vector( BETA_STR, IDX_STR ).size() != 0 )
    //     df.push_back( idat_pair.get_Pair_Vector( BETA_STR, IDX_STR ), BETA_N_STR );
    //   
    //   _S PVAL_N_STR = PVAL_STR+"_"+IDX_STR;
    //   if ( idat_pair.get_Pair_Vector( PVAL_STR, IDX_STR ).size() != 0 )
    //     df.push_back( idat_pair.get_Pair_Vector( PVAL_STR, IDX_STR ), PVAL_N_STR );
    // }
    
    if ( return_df > 2 ) {
      // Build Manifest Data Frame::
      
      //
      // Code below should be redundant.
      //  TBD:: Verify that with some check...
      //
      // if ( idat_pair.get_Info_Vector( PROBE_ID ).size() != 0 )
      //   df.push_back( idat_pair.get_Info_Vector( PROBE_ID ), "Probe_ID" );
      
      if ( idat_pair.get_Pair_Vector( SIG_STR, UG_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( SIG_STR, UG_STR ), UG_STR );
      if ( idat_pair.get_Pair_Vector( SIG_STR, UR_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( SIG_STR, UR_STR ), UR_STR );
      
      if ( idat_pair.get_Pair_Vector( SIG_STR, MG_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( SIG_STR, MG_STR ), MG_STR );
      if ( idat_pair.get_Pair_Vector( SIG_STR, MR_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( SIG_STR, MR_STR ), MR_STR );
      
      if ( idat_pair.get_Info_Vector( COL_STR ).size() != 0 )
        df.push_back( idat_pair.get_Info_Vector( COL_STR ), "col" );
      
      MASK_0_STR = MASK_STR;
      MASK_N_STR = MASK_STR+"_"+IDX_STR;
      if ( idat_pair.get_Bool_Vector( MASK_N_STR ).size() != 0 )
        df.push_back( idat_pair.get_Bool_Vector( MASK_N_STR ), MASK_0_STR );
      
      // Pval and Beta are not easily accesisble since they require the col 
      //   to be provided...
      //
      // - Pval => Det0 + UG_det0 UR_det0 MG_det0 MR_det0
      // - Pval => Det1 + UG_det1 UR_det1 MG_det1 MR_det1
      // - Expand negs_adds[G,R]
      // - Recalculate DetP() => set_pval()
      // - test SigDF()
      // - write_bgz()???
      //
      
      // if ( idat_pair.get_Beta_Vector().size() != 0 )
      //   df.push_back( idat_pair.get_Beta_Vector(), BETA_STR );
      // if ( idat_pair.get_Pval_Vector().size() != 0 )
      //   df.push_back( idat_pair.get_Pval_Vector(), PVAL_STR );
    }
    
    if ( return_df > 3 ) {
      
      /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
       *                              Raw Suffix 0::
       * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
      
      MASK_N_STR = MASK_STR+"_"+IDX_STR;
      if ( idat_pair.get_Bool_Vector( MASK_N_STR ).size() != 0 )
        df.push_back( idat_pair.get_Bool_Vector( MASK_N_STR ), MASK_N_STR );
      
      //
      // Redundant Reporting that is done above...
      //
      BETA_N_STR = BETA_STR+"_"+IDX_STR;
      if ( idat_pair.get_Pair_Vector( BETA_STR, IDX_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( BETA_STR, IDX_STR ), BETA_N_STR );
      
      PVAL_N_STR = PVAL_STR+"_"+IDX_STR;
      if ( idat_pair.get_Pair_Vector( PVAL_STR, IDX_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( PVAL_STR, IDX_STR ), PVAL_N_STR );
      
      DETP_N_STR = DET_STR+"_"+IDX_STR;
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, UG_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, UG_STR ), UG_STR+"_"+DETP_N_STR );
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, UR_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, UR_STR ), UR_STR+"_"+DETP_N_STR );
      
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, MG_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, MG_STR ), MG_STR+"_"+DETP_N_STR );
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, MR_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, MR_STR ), MR_STR+"_"+DETP_N_STR );
    }
    
    if ( return_df > 4 ) {
      
      /*  ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- -----
       *                           Re-calibrated Suffix 1::
       * ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- */
      
      IDX_STR = "1";
      MASK_N_STR = MASK_STR+"_"+IDX_STR;
      if ( idat_pair.get_Bool_Vector( MASK_N_STR ).size() != 0 )
        df.push_back( idat_pair.get_Bool_Vector( MASK_N_STR ), MASK_N_STR );
      
      BETA_N_STR = BETA_STR+"_"+IDX_STR;
      if ( idat_pair.get_Pair_Vector( BETA_STR, IDX_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( BETA_STR, IDX_STR ), BETA_N_STR );
      
      PVAL_N_STR = PVAL_STR+"_"+IDX_STR;
      if ( idat_pair.get_Pair_Vector( PVAL_STR, IDX_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( PVAL_STR, IDX_STR ), PVAL_N_STR );
      
      DETP_N_STR = DET_STR+"_"+IDX_STR;
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, UG_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, UG_STR ), UG_STR+"_"+DETP_N_STR );
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, UR_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, UR_STR ), UR_STR+"_"+DETP_N_STR );
      
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, MG_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, MG_STR ), MG_STR+"_"+DETP_N_STR );
      if ( idat_pair.get_Pair_Vector( DETP_N_STR, MR_STR ).size() != 0 )
        df.push_back( idat_pair.get_Pair_Vector( DETP_N_STR, MR_STR ), MR_STR+"_"+DETP_N_STR );
      
      // if ( idat_pair.get_Bool_Vector( PVAL_STR ).size() != 0 )
      //   df.push_back( idat_pair.get_Bool_Vector( PVAL_STR ), PVAL_STR );
      // if ( idat_pair.get_Bool_Vector( BETA_STR ).size() != 0 )
      //   df.push_back( idat_pair.get_Bool_Vector( BETA_STR ), BETA_STR );
      
      // if ( idat_pair.get_Man_AddM_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_AddM(), "M" );
      // 
      // if ( idat_pair.get_Man_Col1_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_Col1(), "col" );
      // if ( idat_pair.get_Man_ColS_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_ColS(), "colS" );
      // if ( idat_pair.get_Man_ColC_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_ColC(), "colC" );
      //
      // // Default Mask Vector (all unmasked)
      // Rcpp::LogicalVector mask_vec( idat_pair.get_Manifest_Loci_Count(), false );
      // df.push_back( mask_vec, "mask" );
      // 
      // if ( idat_pair.get_Man_MaskS_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_MaskS(), "maskS" );
      // if ( idat_pair.get_Man_MaskC_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_MaskC(), "maskC" );
      // 
      // if ( idat_pair.get_Man_dDGR_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_dDGR(), "dDGR" );
      // if ( idat_pair.get_Man_Code_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_Code(), "Det_Code" );
      // 
      // if ( idat_pair.get_Man_Keys_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_Keys(), "Manifest" );
      
      // if ( idat_pair.get_Man_DetUG_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_DetUG(), "UG_Det" );
      // if ( idat_pair.get_Man_DetUR_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_DetUR(), "UR_Det" );
      // if ( idat_pair.get_Man_DetMG_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_DetMG(), "MG_Det" );
      // if ( idat_pair.get_Man_DetMR_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_DetMR(), "MR_Det" );
      // 
      // if ( idat_pair.get_Man_SigUG_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_SigUG(), "UG_Sig" );
      // if ( idat_pair.get_Man_SigUR_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_SigUR(), "UR_Sig" );
      // if ( idat_pair.get_Man_SigMG_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_SigMG(), "MG_Sig" );
      // if ( idat_pair.get_Man_SigMR_Cnt() != 0 )
      //   df.push_back( idat_pair.get_Man_SigMR(), "MR_Sig" );
    }
    
    // idat_pair.to_sample_sheet();
  }
  
  /* ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** 
   *                                    DONE::
   * ***** ***** ***** ***** ***** ***** | ***** ***** ***** ***** ***** ***** */
  
  if ( p0 ) Rcpp::Rcerr<<_fo<<"Done. Return: '"<<success<<"'\n"<<std::endl;
  
  // type.convert( x, na.strings = "NA", as.is, ... )
  return( typeConvert( asTibble( df ), "", true ) );
}

/*** R

run_test <- TRUE
run_test <- FALSE

if ( run_test ) {
  suppressWarnings(suppressPackageStartupMessages( 
    base::require("tidyverse",  quietly = TRUE) ) )
  suppressWarnings(suppressPackageStartupMessages( 
    base::require("dplyr",  quietly = TRUE) ) )
  suppressWarnings(suppressPackageStartupMessages( 
    base::require("tidyr",  quietly = TRUE) ) )
  suppressWarnings(suppressPackageStartupMessages( 
    base::require("tibble",  quietly = TRUE) ) )
  suppressWarnings(suppressPackageStartupMessages( 
    base::require("sesame",  quietly = TRUE) ) )
  
  opt_rcpp <- NULL
  opt_rcpp$top_path <- "/Users/bbarnes/mylocal/Documents"
  opt_rcpp$prgm_tag <- "idat_pair"
  
  pmssg <- glue::glue("[{opt_rcpp$prgm_tag}]:")
  pwarn <- glue::glue("{RET}[{opt_rcpp$prgm_tag}]: Warning:")
  perrs <- glue::glue("{RET}[{opt_rcpp$prgm_tag}]: ERROR:")
  
  source( file.path( opt_rcpp$top_path, "tools/imSuite/scripts/R/functions/common_functions.R" ) )
  source( file.path( opt_rcpp$top_path, "tools/imSuite/scripts/R/sdf/functions/sdf_base_functions.R" ) )
  source( file.path( opt_rcpp$top_path, "tools/imSuite/scripts/R/sdf/functions/sdf_functions.R" ) )
  source( file.path( opt_rcpp$top_path, "tools/imSuite/scripts/R/functions/timeTracker.R" ) )
  
  cat(glue::glue("{pmssg} Running Test...{RET}"))
  #
  # IDATS::
  #
  pre_vec <- c(
    # EPICv1: HELA???
    file.path( opt_rcpp$top_path, "data/idats/idats_EPIC-8x1-DELTA-Core/203319730005/203319730005_R01C01" ),
    # EPICv2: HELA
    file.path( opt_rcpp$top_path, "data/idats/idats_EPIC_v2-20220912-Alpha/206891110001/206891110001_R01C01" )
    # MSA: HELA
    # file.path( opt_rcpp$top_path, "data/pre-idats/MSA/SNG.07122023_PreDVT_SGrun/07122023_PreDVT_SGrun/207705000015/207705000015_R01C01" )
  )
  
  #
  # Rcpp Manifests::
  #
  all_sub_rds  <- file.path( opt_rcpp$top_path, "data/manifests/methylation/bgz/all_manifests.sub8.rds")
  rev_sub_rds  <- file.path( opt_rcpp$top_path, "data/manifests/methylation/bgz/rev_manifests.sub8.rds")
  neg_ctl_rds  <- file.path( opt_rcpp$top_path, "data/manifests/methylation/bgz/all_negative_ctls.rds" )
  
  all_sub_tib  <- NULL
  all_sub_tib  <- readr::read_rds( all_sub_rds )
  
  neg_ctl_tib  <- NULL
  neg_ctl_tib  <- readr::read_rds( neg_ctl_rds )
  
  #
  # TBD: Update manifest/controls with MSA and EPICv2
  #
  
  #
  # Output Directory::
  #
  opt_rcpp$output_dir <- "/Users/bbarnes/mylocal/Documents/tmp/idat_pair"
  if ( !dir.exists( opt_rcpp$output_dir ) ) base::dir.create( opt_rcpp$output_dir, recursive = TRUE )
  cat( glue::glue("[idat_pair]: Prefix: 'opt_rcpp$prefix'\n\n") )
  
  #
  # Run Time Parameters::
  #
  opt_rcpp$tc <- 0
  opt_rcpp$vt <- 0
  
  # opt_rcpp$vb <- 5 # Two level deep (p3)
  # opt_rcpp$vb <- 4 # One level deep (p2)
  # opt_rcpp$vb <- 3 # Standard
  opt_rcpp$vb <- 2 # Light
  # opt_rcpp$vb <- 1 # Min
  # opt_rcpp$vb <- 0 # None
  
  success <- FALSE
  
  opt_rcpp$min_pval <- 0.05
  opt_rcpp$min_beta <- 0.3
  opt_rcpp$max_beta <- 0.7
  opt_rcpp$min_perO <- 0.75
  opt_rcpp$min_perI <- 0.05
  
  # Mask Methods:: csJjKk
  # workflow_vec = c( "dc",
  #                    "ds",
  #                    "dJ",
  #                    "dj",
  #                    "dK",
  #                    "dk" )
  workflow_vec = c( "cd",
                    "sd",
                    "Jd",
                    "jd",
                    "Kd",
                    "kd" )
  workflow_vec = c( "id" )
  
  workflow_vec = c( "i" )
  
  return_df <- 0
  return_df <- 1
  return_df <- 2
  return_df <- 3
  return_df <- 4
  
  opt_rcpp$run_pair <- FALSE
  opt_rcpp$run_pair <- TRUE
  
  opt_rcpp$run_sesame <- FALSE
  opt_rcpp$run_sesame <- TRUE
  
  for ( prefix in pre_vec ) {
    prefix <- pre_vec[1]
    # prefix <- pre_vec[2]
    # prefix <- pre_vec[3]
    
    # sesame::SigDF()
    # [TBD]: Add sesame code to run same workflow...
    ses_sdf <- NULL
    ses_beta_tib <- NULL
    if ( opt_rcpp$run_sesame ) {
      ses_sdf <- prefix_to_sdf( 
        prefix = prefix, 
        platform = "", # NULL, 
        manifest = NULL, # all_sub_tib, 
        controls = NULL, # neg_ctl_tib,
        out_dir = file.path( opt_rcpp$output_dir ),
        run_tag = "sesame",
        reload     = 0,
        reload_min = 2,
        reload_pre = NULL,
        ret_data   = FALSE,
        parallel   = FALSE,
        vb=opt_rcpp$vb,vt=opt_rcpp$vt,tc=opt_rcpp$tc ) 
      
      ses_ctl_tib <- NULL
      ses_ctl_tib <- ses_sdf %>% sesame::controls() %>% 
        # dplyr::filter( Type == "NEGATIVE" ) %>% 
        tibble::as_tibble()
      # ses_ctl_tib %>% print_sum( vec = c("Type") )
      
      # ses_sdf <- ses_sdf %>% mutate_sdf_simple( steps = "NDCDPB", negs_min = opt_rcpp$min_pval )
      
      # tmp_tib %>% head() %>% cbind() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID") %>% tibble::as_tibble() %>% magrittr::set_names( value = c("Probe_ID","Beta") )
      ses_beta_tib <- mutate_sdf_simple( 
        sdf = ses_sdf,
        # steps = "V", # With masking...
        steps = "v", # Without masking...
        negs_min = 1.0, # negs_min = 0.05,
        poob_min = 1.0, # poob_min = 0.05,
        probe_id = "Probe_ID", beta_key = "Ses_Beta",
        ret_tib = TRUE,
        vb=opt_rcpp$vb, vt=opt_rcpp$vt, tc=opt_rcpp$tc )
      
      # mutate_sdf_simple( steps = "N", negs_min = opt_rcpp$min_pval )
      # ses_sdf <- ses_sdf %>% mutate_sdf_simple( steps = workflow_vec[1], vb=opt_rcpp$vb,vt=opt_rcpp$vt,tc=opt_rcpp$tc )
      
      # Workflow Notes::
      #
      # D   = dyeBiasNL
      # S   = formatVCF
      # C/c = inferInfiniumIChannel( switchFail+/switchFail-)
      # N   = detectionPnegEcdf
      # P/p = pOOBAH (neg+/neg-)
      # B/b = noob (neg+/-neg-)
      # 
      # [TBD]: Apply methods above to both ses_sdf and cpp_tib...
      # [TBD]: Use ssh_tib -> pre_vec to choose a testable dataset...
      # [TBD]: cpp_tib add to auto sample sheet number of flips
      #
      # Workflows: ND[Cc]+[D][PpN]+[Bb]+
      
      ses_tib <- NULL
      ses_tib <- ses_sdf %>% tibble::as_tibble()
      
    }
    
    cpp_tib <- NULL
    if ( opt_rcpp$run_pair ) {
      cpp_tib <- read_idat_pair_rcpp( 
        prefix_path      = prefix, 
        output_path      = opt_rcpp$output_dir, 
        workflow_vec     = workflow_vec %>% as.vector(),
        
        pval_add_vec     = neg_ctl_tib$Address %>% as.vector(), 
        addU_man_vec     = all_sub_tib$U %>% as.vector(),
        addM_man_vec     = all_sub_tib$M %>% as.vector(),
        
        cgns_man_vec     = all_sub_tib$Probe_ID %>% as.vector(), 
        cols_man_vec     = all_sub_tib$col %>% as.vector(),
        keys_man_vec     = all_sub_tib$Manifest %>% as.vector(),
        anns_man_vec     = all_sub_tib$Annotation %>% as.vector(), 
        chrs_man_vec     = all_sub_tib$Chromosome %>% as.vector(),
        
        min_pval         = opt_rcpp$min_pval,
        min_beta         = opt_rcpp$min_beta,
        max_beta         = opt_rcpp$max_beta,
        min_perO         = opt_rcpp$min_perO, 
        min_perI         = opt_rcpp$min_perI, 
        read_bgz         = FALSE,
        write_bgz        = FALSE,
        rm_pval_outliers = FALSE,
        return_df        = return_df,
        vb=opt_rcpp$vb,vt=opt_rcpp$vt,tc=opt_rcpp$tc )
      
      cpp_ctl_tib <- NULL
      cpp_ctl_tib <- cpp_tib %>% 
        dplyr::filter( stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
        dplyr::mutate( col = as.factor(col) )
      
      cpp_ctl_tib %>% dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative") )
      
      cpp_sdf <- NULL
      cpp_sdf <- cpp_tib %>%
        dplyr::mutate( col = as.factor(col)
                       # U = dplyr::case_when( U==0 ~ NA_real_, TRUE ~ U ),
                       # M = dplyr::case_when( M==0 ~ NA_real_, TRUE ~ M )
        ) %>%
        dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") ) %>% 
        dplyr::select( Probe_ID:mask ) %>%
        as.data.frame() %>%
        sesame::SigDF( # platform = "EPIC", 
          ctl = NULL
          # cpp_ctl_tib %>% 
          #     dplyr::filter( !mask ) %>%
          #     dplyr::select( Probe_ID:mask ) %>%
          #     as.data.frame()
        ) %>% dplyr::select( Probe_ID:mask )
      
      cpp_sdf %>% dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative") ) %>% tibble::as_tibble()
      cpp_sdf %>% sesame::controls() %>% 
        dplyr::filter( Type == "NEGATIVE" ) %>% 
        tibble::as_tibble()
      
      cpp_beta_tib <- NULL
      cpp_beta_tib <- mutate_sdf_simple( 
        sdf = cpp_sdf,
        # steps = "V", # With masking...
        steps = "v", # Without masking...
        negs_min = 1.0, # negs_min = 0.05,
        poob_min = 1.0, # poob_min = 0.05,
        probe_id = "Probe_ID", beta_key = "Cpp_Beta",
        ret_tib = TRUE ) # , vb=vb, vt=vt, tc=tc )
      
      cpp_pval_tibM <- NULL
      cpp_pval_tibM <- mutate_sdf_simple( 
        sdf = cpp_sdf,
        steps = "n", # Without masking...
        negs_min = 1.0, # negs_min = 0.05,
        poob_min = 1.0, # poob_min = 0.05,
        probe_id = "Probe_ID", beta_key = "Cpp_Pval",
        ret_tib = TRUE,
        vb=vb, vt=vt, tc=tc )
      
      cpp_pval_tibm <- NULL
      cpp_pval_tibm <- mutate_sdf_simple( 
        sdf = cpp_sdf,
        steps = "n", # Without masking...
        negs_min = 1.0, # negs_min = 0.05,
        poob_min = 1.0, # poob_min = 0.05,
        probe_id = "Probe_ID", beta_key = "Cpp_Pval",
        ret_tib = TRUE,
        vb=vb, vt=vt, tc=tc )
      
      cpp_pval_tibM <- cpp_pval_tib
      
      cpp_pval_both_tib <- cpp_pval_tibm %>% dplyr::inner_join( cpp_pval_tibM, by=c("Probe_ID"), suffix=c("_m","_M") )
      
    }
    
    if ( !is.null(ses_sdf) ) {
      cat(glue::glue("{RET}{pmssg} ses_sdf::{RET}"))
      ses_tib %>% 
        dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
        dplyr::arrange(Probe_ID) %>% 
        print()
      
      # ses_sdf %>% dplyr::filter( mask ) %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) ) %>% print_sum( vec = c("Probe_Type") )
      #
      # EPIC default Controls::
      # ses_tib %>% dplyr::filter( mask ) %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) ) %>% dplyr::filter( Probe_Type == "ct" ) %>% print(n=1000)
      #
    }
    
    if ( !is.null(cpp_tib) ) {
      cat(glue::glue("{RET}{pmssg} cpp_tib::{RET}"))
      cpp_tib %>% 
        dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
        dplyr::arrange(Probe_ID) %>% 
        print()
    }
    
    
    inn_tib <- NULL
    inn_tib <- dplyr::inner_join( ses_tib,cpp_tib, by=c("Probe_ID"), suffix=c("_ses","_cpp") ) %>%
      dplyr::inner_join( ses_beta_tib, by=c("Probe_ID") )
    
    dif_tib <- inn_tib %>% 
      dplyr::mutate( 
        Beta_0 = round(Beta_0,6),
        Ses_Beta = round(Ses_Beta,6),
        dif = Beta_0 - Ses_Beta ) %>% 
      dplyr::filter( dif != 0 )
    # dplyr::arrange( dif ) %>% head() %>% as.data.frame()
    
    ant_ses_tib <- dplyr::anti_join( ses_tib,cpp_tib, by=c("Probe_ID") )
    ant_cpp_tib <- dplyr::anti_join( cpp_tib,ses_tib, by=c("Probe_ID") )
    
    # [TBD]: Compare ses_sdf vs. cpp_tib
    if ( FALSE ) {
      #
      # Test Sesame SigDF functionality::
      #
    }
  }
  cat( glue::glue("\n\nDONE(test_read_idats): success='{success}'\n\n") )
}

*/
