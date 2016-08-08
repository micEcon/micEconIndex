
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
// [[Rcpp::export]]
arma::mat fisherEKSsparse(const arma::sp_mat  Q_consol,  const arma::sp_mat  P_consol,
                             arma::mat     Q_freq,          arma::mat     P_freq,
                             arma::uvec    Q_ind,           arma::uvec    P_ind) {
                              //     arma::mat split_size) {

  // auto x = 10;
  // cout << "Got here, after data read-in " << endl;
  uword base_period = 0 ;
  double M_dbl = Q_ind.n_elem ;
  uword split_interval = ceil( as_scalar(2e+12) / (P_consol.n_rows * Q_consol.n_rows)) ;
  // 1e+08 ; 1e+6; 1e+2; 1e+12; 2e+10 2e+6; Best: split_size = 1e+13
  // The 1e+08 will likely create a bug if there are more than 100 million firms
  // 1e+08 is designed to keep the Q_split_x_P under one gig of memory
  // cout << "Got here, right before spl creation " << endl;
  uvec spl = regspace<uvec>(1, split_interval, Q_consol.n_rows - 1) - 1L  ;
  // subtracting 1 from P_consol.n_rows so that I can append it later
  // without running into a bug. Using the useful behavior of regspace here
  // cout << "Got here, after spl creation " << endl;
  spl.resize(spl.n_elem + 1) ;
  // cout << "Got here, after spl resize " << endl;
  spl(spl.n_elem - 1) = Q_consol.n_rows  ;
  // preventing off-by-one error.

  arma::vec Q_x_P_vec(Q_consol.n_rows) ;
  Q_x_P_vec.zeros() ;
  arma::rowvec P_x_Q_vec(P_consol.n_rows) ;
  P_x_Q_vec.zeros() ;
  arma::vec revenue(Q_ind.n_elem) ;
  revenue.zeros() ;

  // cout << "spl.n_rows = " << spl.n_rows << std::endl;
  // // cout << "spl values: " << endl;
  // cout << "spl element 0 = " << spl(0) << std::endl;
  // cout << "spl element 1 = " << spl(1) << std::endl;
  // cout << "spl elements = " << std::endl;
 // spl.print() ;

  // cout << "Got here, right before loop " << endl;

//   stop("") ;

  inplace_trans( Q_freq );

 // umat Q_P_ind_combined = join_rows(Q_ind, P_ind ) ;

//  inplace_trans( Q_P_ind_combined );

 // Progress p(0, false);

  for (uword spl_i=0; spl_i < (spl.n_elem - 1); spl_i++) {
    // cout << "iteration " << spl_i << " of " << (spl.n_elem - 1) << endl;
    // Note that col < split.n_elem is strictly less than, so
    // it won't run over the last element.
    // uvec split_span = span(split(split_it), split(split_it + 1)) ;
   // A(p:q, :) 	  	A.rows(p, q)

  //  if (Progress::check_abort() )
  //          stop("User abort");

     arma::wall_clock timer;
     timer.tic();

    arma::mat Q_split_x_P = log( conv_to<mat>::from(Q_consol.rows(
      spl(spl_i), spl(spl_i + 1) - 1) * P_consol.t() )) ;

    // cout << "Main matrix operation time = " << timer.toc()  << endl;

    // cout << "Got here, after Q_split_x_P  creation " << endl;

    timer.tic();

    Q_x_P_vec(span( spl(spl_i), spl(spl_i + 1) - 1 )) = Q_split_x_P * P_freq ;
    // P_freq is a column vec.
    // In a sense I am appending these values
    // cout << "Q_x_P_vec operation time = " << timer.toc()  << endl;

    // cout << "Got here, after Q_x_P_vec fill " << endl;
    // cout << "Q_split_x_P.n_rows = " << Q_split_x_P.n_rows << std::endl;
   //  // cout << "Q_split_x_P.n_cols = " << Q_split_x_P.n_cols << std::endl;

    // P_x_Q_vec(span( spl(spl_i), spl(spl_i + 1) - 1 )) =
    //   (Q_split_x_P * Q_freq).t() ;
    // I changd this, and it may be a problem
    // P_x_Q_vec += Q_freq.t() % Q_split_x_P ;
     // span( spl(spl_i), spl(spl_i + 1) - 1 )
    timer.tic();
     P_x_Q_vec += Q_freq.cols(spl(spl_i), spl(spl_i + 1) - 1 ) * Q_split_x_P ;

     // cout << "P_x_Q_vec operation time = " << timer.toc()  << endl;

     timer.tic();
    // Q_freq originally comes in as a column vector
    // P_x_Q_vec is a row vector. Add in place
    // cout << "Got here, after P_x_Q_vec addition " << endl;

   // uvec Q_split_ind = regspace<uvec>(spl(spl_i), spl(spl_i + 1) - 1)  ;
  //  // cout << "Got here, after Q_split_ind creation " << endl;

   // uvec Q_split_supra_ind = conv_to<uvec>::from(
    //  sugar_in(
    //  conv_tostd::vector<int>::from(Q_ind),
    //  conv_to<int>::from(Q_split_ind) )
    //);
   // Q_split_ind.print();
  //   arma::wall_clock timer;
   //  timer.tic();

  // uvec Q_split_supra_ind =  myInOperator(Q_ind, Q_split_ind)   ; // - 1
   // uvec Q_split_supra_ind =  Q_ind.rows(find(Q_ind >= spl(spl_i) && Q_ind <= spl(spl_i + 1) - 1)) ;

   //uvec Q_split_supra_ind =  find( ( (Q_ind >= spl(spl_i) ) +
    // (Q_ind <= (spl(spl_i + 1))) ) == 2) ;

      uvec logic_1 = Q_ind >= spl(spl_i) ;
      uvec logic_2 = Q_ind <= (spl(spl_i + 1) - 1) ;
      uvec Q_split_supra_ind = find( (logic_1 + logic_2) == 2) ;

   // Jumping through a bunch of hoops to get compound logical expressions

    // // cout << "IN operator time = " << timer.toc()  << endl;
   // Had to convert a lot since Rcpp's sugar "in"
   // is written for a certain data ype
   // cout << "Got here, after Q_split_supra_ind creation " << endl;
   // Q_split_supra_ind.print();
   // Q_ind.print() ;
   // P_ind.print() ;
   // cout << "Q_split_x_P.n_rows " << Q_split_x_P.n_rows << endl;
   // cout << "Q_split_x_P.n_cols " << Q_split_x_P.n_cols << endl;

  // uvec test = conv_to<uvec>::from(Q_ind(Q_split_supra_ind) +
  //   (conv_to<vec>::from(P_ind(Q_split_supra_ind)) - 1L) * Q_split_x_P.n_rows) ;
   // uvec test = conv_to<uvec>::from(Q_ind(Q_split_supra_ind) + (P_ind) * Q_split_x_P.n_rows) ;
   // convert to vec very important because subtracting uvec below zero makes
   // integer underflow
    umat testt = join_rows(Q_ind(Q_split_supra_ind) - spl(spl_i),
                          P_ind(Q_split_supra_ind) );
   // arma::mat testtprint = conv_to<mat>::from(testt)  ;
  // testtprint.print() ;
   //  size(Q_split_x_P).print() ;

  // umat Q_P_ind_combined_temp = Q_P_ind_combined.cols(Q_split_supra_ind) ;
  // Q_P_ind_combined_temp.row(1) -= spl(spl_i) ;

   // uvec test =  ;

    // test.print();
    // stop("STOP");
    // // cout << "max test value " << max(test) << endl;
    // cout << "Q_split_x_P.n_elem " << Q_split_x_P.n_elem << endl;
   // Q_P_ind_combined_temp.print();

    revenue(Q_split_supra_ind) =
       Q_split_x_P( arma::sub2ind( size(Q_split_x_P), testt.t() )) ;
      // Q_split_x_P( arma::sub2ind( size(Q_split_x_P), Q_P_ind_combined_temp ) ) ;

    // cout << "Got here, after revenue fill " << endl;
    // cout << "Rest of operations in the loop = " << timer.toc()  << endl;

  }

  // cout << "Got here, after loop done " << endl;

  double sum_revenue = sum(revenue) ;


   double top_row  = (-1) * M_dbl * arma::as_scalar(revenue(base_period))  +
        sum_revenue -
        Q_x_P_vec(base_period) +
        P_x_Q_vec(base_period) ;

  // cout << "Got here, after top_row calc " << endl;

  arma::vec ret = (1 / (2 * M_dbl)) * (
        top_row +
        M_dbl * revenue  -
        sum_revenue  +
        Q_x_P_vec(Q_ind) -
        P_x_Q_vec(P_ind)
    ) ;
    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "



  return(exp(ret)) ;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
// [[Rcpp::export]]
arma::mat fisherEKSdense(const arma::mat  Q_consol,  const arma::mat  P_consol,
                             arma::mat     Q_freq,          arma::mat     P_freq,
                             arma::uvec    Q_ind,           arma::uvec    P_ind) {
                              //     arma::mat split_size) {

  // auto x = 10;
  // cout << "Got here, after data read-in " << endl;
  uword base_period = 0 ;
  double M_dbl = Q_ind.n_elem ;
  uword split_interval = ceil( as_scalar(2e+12) / (P_consol.n_rows * Q_consol.n_rows)) ;
  // 1e+08 ; 1e+6; 1e+2; 1e+12; 2e+10 2e+6; Best: split_size = 1e+13
  // The 1e+08 will likely create a bug if there are more than 100 million firms
  // 1e+08 is designed to keep the Q_split_x_P under one gig of memory
  // cout << "Got here, right before spl creation " << endl;
  uvec spl = regspace<uvec>(1, split_interval, Q_consol.n_rows - 1) - 1L  ;
  // subtracting 1 from P_consol.n_rows so that I can append it later
  // without running into a bug. Using the useful behavior of regspace here
  // cout << "Got here, after spl creation " << endl;
  spl.resize(spl.n_elem + 1) ;
  // cout << "Got here, after spl resize " << endl;
  spl(spl.n_elem - 1) = Q_consol.n_rows  ;
  // preventing off-by-one error.

  arma::vec Q_x_P_vec(Q_consol.n_rows) ;
  Q_x_P_vec.zeros() ;
  arma::rowvec P_x_Q_vec(P_consol.n_rows) ;
  P_x_Q_vec.zeros() ;
  arma::vec revenue(Q_ind.n_elem) ;
  revenue.zeros() ;

  // cout << "spl.n_rows = " << spl.n_rows << std::endl;
  // // cout << "spl values: " << endl;
  // cout << "spl element 0 = " << spl(0) << std::endl;
  // cout << "spl element 1 = " << spl(1) << std::endl;
  // cout << "spl elements = " << std::endl;
 // spl.print() ;

  // cout << "Got here, right before loop " << endl;

//   stop("") ;

  inplace_trans( Q_freq );

 // umat Q_P_ind_combined = join_rows(Q_ind, P_ind ) ;

//  inplace_trans( Q_P_ind_combined );

 // Progress p(0, false);

  for (uword spl_i=0; spl_i < (spl.n_elem - 1); spl_i++) {
    // cout << "iteration " << spl_i << " of " << (spl.n_elem - 1) << endl;
    // Note that col < split.n_elem is strictly less than, so
    // it won't run over the last element.
    // uvec split_span = span(split(split_it), split(split_it + 1)) ;
   // A(p:q, :) 	  	A.rows(p, q)

  //  if (Progress::check_abort() )
  //          stop("User abort");

     arma::wall_clock timer;
     timer.tic();

    arma::mat Q_split_x_P = log( Q_consol.rows(
      spl(spl_i), spl(spl_i + 1) - 1) * P_consol.t() ) ;

    // cout << "Main matrix operation time = " << timer.toc()  << endl;

    // cout << "Got here, after Q_split_x_P  creation " << endl;

    timer.tic();

    Q_x_P_vec(span( spl(spl_i), spl(spl_i + 1) - 1 )) = Q_split_x_P * P_freq ;
    // P_freq is a column vec.
    // In a sense I am appending these values
    // cout << "Q_x_P_vec operation time = " << timer.toc()  << endl;

    // cout << "Got here, after Q_x_P_vec fill " << endl;
    // cout << "Q_split_x_P.n_rows = " << Q_split_x_P.n_rows << std::endl;
   //  // cout << "Q_split_x_P.n_cols = " << Q_split_x_P.n_cols << std::endl;

    // P_x_Q_vec(span( spl(spl_i), spl(spl_i + 1) - 1 )) =
    //   (Q_split_x_P * Q_freq).t() ;
    // I changd this, and it may be a problem
    // P_x_Q_vec += Q_freq.t() % Q_split_x_P ;
     // span( spl(spl_i), spl(spl_i + 1) - 1 )
    timer.tic();
     P_x_Q_vec += Q_freq.cols(spl(spl_i), spl(spl_i + 1) - 1 ) * Q_split_x_P ;

     // cout << "P_x_Q_vec operation time = " << timer.toc()  << endl;

     timer.tic();
    // Q_freq originally comes in as a column vector
    // P_x_Q_vec is a row vector. Add in place
    // cout << "Got here, after P_x_Q_vec addition " << endl;

   // uvec Q_split_ind = regspace<uvec>(spl(spl_i), spl(spl_i + 1) - 1)  ;
  //  // cout << "Got here, after Q_split_ind creation " << endl;

   // uvec Q_split_supra_ind = conv_to<uvec>::from(
    //  sugar_in(
    //  conv_tostd::vector<int>::from(Q_ind),
    //  conv_to<int>::from(Q_split_ind) )
    //);
   // Q_split_ind.print();
  //   arma::wall_clock timer;
   //  timer.tic();

  // uvec Q_split_supra_ind =  myInOperator(Q_ind, Q_split_ind)   ; // - 1
   // uvec Q_split_supra_ind =  Q_ind.rows(find(Q_ind >= spl(spl_i) && Q_ind <= spl(spl_i + 1) - 1)) ;

   //uvec Q_split_supra_ind =  find( ( (Q_ind >= spl(spl_i) ) +
    // (Q_ind <= (spl(spl_i + 1))) ) == 2) ;

      uvec logic_1 = Q_ind >= spl(spl_i) ;
      uvec logic_2 = Q_ind <= (spl(spl_i + 1) - 1) ;
      uvec Q_split_supra_ind = find( (logic_1 + logic_2) == 2) ;

   // Jumping through a bunch of hoops to get compound logical expressions

    // // cout << "IN operator time = " << timer.toc()  << endl;
   // Had to convert a lot since Rcpp's sugar "in"
   // is written for a certain data ype
   // cout << "Got here, after Q_split_supra_ind creation " << endl;
   // Q_split_supra_ind.print();
   // Q_ind.print() ;
   // P_ind.print() ;
   // cout << "Q_split_x_P.n_rows " << Q_split_x_P.n_rows << endl;
   // cout << "Q_split_x_P.n_cols " << Q_split_x_P.n_cols << endl;

  // uvec test = conv_to<uvec>::from(Q_ind(Q_split_supra_ind) +
  //   (conv_to<vec>::from(P_ind(Q_split_supra_ind)) - 1L) * Q_split_x_P.n_rows) ;
   // uvec test = conv_to<uvec>::from(Q_ind(Q_split_supra_ind) + (P_ind) * Q_split_x_P.n_rows) ;
   // convert to vec very important because subtracting uvec below zero makes
   // integer underflow
    umat testt = join_rows(Q_ind(Q_split_supra_ind) - spl(spl_i),
                          P_ind(Q_split_supra_ind) );
   // arma::mat testtprint = conv_to<mat>::from(testt)  ;
  // testtprint.print() ;
   //  size(Q_split_x_P).print() ;

  // umat Q_P_ind_combined_temp = Q_P_ind_combined.cols(Q_split_supra_ind) ;
  // Q_P_ind_combined_temp.row(1) -= spl(spl_i) ;

   // uvec test =  ;

    // test.print();
    // stop("STOP");
    // // cout << "max test value " << max(test) << endl;
    // cout << "Q_split_x_P.n_elem " << Q_split_x_P.n_elem << endl;
   // Q_P_ind_combined_temp.print();

    revenue(Q_split_supra_ind) =
       Q_split_x_P( arma::sub2ind( size(Q_split_x_P), testt.t() )) ;
      // Q_split_x_P( arma::sub2ind( size(Q_split_x_P), Q_P_ind_combined_temp ) ) ;

    // cout << "Got here, after revenue fill " << endl;
    // cout << "Rest of operations in the loop = " << timer.toc()  << endl;

  }

  // cout << "Got here, after loop done " << endl;

  double sum_revenue = sum(revenue) ;


   double top_row  = (-1) * M_dbl * arma::as_scalar(revenue(base_period))  +
        sum_revenue -
        Q_x_P_vec(base_period) +
        P_x_Q_vec(base_period) ;

  // cout << "Got here, after top_row calc " << endl;

  arma::vec ret = (1 / (2 * M_dbl)) * (
        top_row +
        M_dbl * revenue  -
        sum_revenue  +
        Q_x_P_vec(Q_ind) -
        P_x_Q_vec(P_ind)
    ) ;
    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "



  return(exp(ret)) ;
}







// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fisherInd (arma::mat Q, arma::mat P, int base_period) {
  // Fisher index; direct approach; yes transitivity normalization
  //int M = Q.n_rows ;
  //int K = Q.n_cols ;
  //arma::mat check_nonneg_zeros(M, K) ;
  //check_nonneg_zeros.zeros() ;

  //arma::mat Q_nonneg = Q >= check_nonneg_zeros.zeros()
  //bool all_nonneg_check_Q = arma::as_scalar( arma::all(arma::all(Q_nonneg, 0), 1) ) ;
  // bool all_nonneg_check_Q = arma::as_scalar( arma::all(arma::all(Q >= 0, 0), 1) ) ;
  // bool all_nonneg_check_P = arma::as_scalar( arma::all(arma::all(P >= 0, 0), 1) )  ;
  //if ( ! (all_nonneg_check_Q && all_nonneg_check_P) ) {
  //  stop("All elements of the quantity and price matrices must be non-negative.") ;
  //}

  // Eventually I want to check that all elemets are posi, but I failed to figure out how to do it
  // Probably better to wrap those checks within R code.

  base_period = base_period - 1 ;
  // Since C++ is zero-indexed
  int M = Q.n_rows ;
  // double M_double = M ;
  // arma::mat x(M, M) ;
  arma::vec I_row(M);
  arma::vec I_col(M);
  arma::vec ret(M);
  I_row.ones() ;
  I_col.ones() ;
  ret.ones() ;

  for (int col=0; col<M; col++) {
    double Q_ind_L = dot(P.row(base_period), Q.row(col)) / dot(P.row(base_period), Q.row(base_period)) ;
    double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(base_period)) ;
    I_row(col) = sqrt(Q_ind_L * Q_ind_P) ;
  }

  for (int col=0; col<M; col++) {

    for (int row=0; row<M; row++) {
      double Q_ind_L = dot(P.row(row), Q.row(col)) / dot(P.row(row), Q.row(row)) ;
      double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(row)) ;
      I_col(row) = sqrt(Q_ind_L * Q_ind_P) ;
    }

    arma::vec interm_vec_prod = I_row % I_col ;

      // R// cout << "The value is " << sum(interm_vec_prod) << std::endl ;
      if ( ! arma::is_finite(interm_vec_prod) ) {
        // Note that arma::is_finite checks the whole vector for any non-finite values
        stop("NaNs produced in quantity index. Check the quantity and price matrix inputs. Quantity indices must be positive, so the product of quantities and prices must be positive in all cases.") ;
        // Throws error in the case of zeros in the computed quantity index, which
        // will only occur in the case of all zeros in quantity
        // Thanks to explanation here: http://gallery.rcpp.org/articles/intro-to-exceptions/
      }
      ret(col) = exp(mean(log(interm_vec_prod))) ;
  }

  return(ret) ;
}





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
if ( T) {

library(data.table)
library(Matrix)

consol.matrix <- function(x) {
  if (!is.data.table(x)) x <- as.data.table(x)
  x.ret <- x[, .(.N), by = names(x)]
  N.ret <- matrix(x.ret$N, ncol = 1)
  x.ret[, N := NULL]
  list(mat = as.matrix(x.ret), freq = N.ret)
}



set.seed(100)
# n.col <- 100; n.row = 40000
# With these params, fastest index fn get 77 secs. Faster index fn gets 320 secs (4 times faster):
#  n.row.fact <- 100000 ; real.rows.factor = 2 ; n.col <- 400;
# n.row.fact <- 100000 ; real.rows.factor = 2 ; n.col <- 10;
# With the below, I have fastest 0.13; faster 0.185; naive 18.4 secs :
# n.row.fact <- 1000 ; real.rows.factor = 5 ; n.col <- 300;
# With below, I get fastest 0.013; faster 0.014; naive 112.533:
# n.row.fact <- 100 ; real.rows.factor = 100 ; n.col <- 100;
 n.row.fact <- 10 ; real.rows.factor = 2 ; n.col <- 4;
#n.row.fact <- 1000 ; real.rows.factor = 1 ; n.col <- 400;
n.row = real.rows.factor; n.row = n.row * n.row.fact
n.real.rows = n.row / real.rows.factor
P.mat <- matrix(runif(n.real.rows*n.col), ncol = n.col, nrow = n.row, byrow = TRUE )
P.mat <- rbind(P.mat[-1, ], P.mat[1, ])
#P.mat <- rbind(P.mat[1, ], P.mat[1, ], P.mat[2, ], P.mat[2, ])
#P.mat <- matrix(runif(n.col*n.row), nrow = n.row )
# Q.mat <- matrix(runif(n.col*n.row), ncol = n.col)
Q.mat <- matrix(runif(n.real.rows*n.col), ncol = n.col, nrow = n.row, byrow = TRUE )
Q.mat[, 4:ncol(Q.mat)] <- 0
# Making the matrix sparse


Q.mat.consol <- consol.matrix(Q.mat)
P.mat.consol <- consol.matrix(P.mat)
rm(Q.mat)
rm(P.mat)

if (F) {
  print( system.time( fisherInd.ret <- fisherInd(Q.mat, P.mat, 1) ) )
}

if (F) {
print(system.time(
fisherIndfast.ret <-
  fisherIndfast(Q = Q.mat, P = P.mat,
                Q_consol = Q.mat.consol$mat,
                P_consol = P.mat.consol$mat,
                Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq ) # t(P.mat.consol$freq ))
))
}


if (F) {
print(system.time(
fisherIndfaster.ret <- fisherIndfaster(Q_consol = Q.mat.consol$mat,
                P_consol = P.mat.consol$mat,
                Q_freq = t(Q.mat.consol$freq),
                #Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq,
                Q_ind = rep((1:n.real.rows) - 1, real.rows.factor),
                P_ind = rep((1:n.real.rows) - 1, real.rows.factor))
                # P_ind = c(rep((1:n.real.rows) - 1, real.rows.factor)[-1], 0))
))
}



if (F) {
print(system.time(
fisherIndfastest.ret <- fisherIndfastest(
              # Q_consol = Q.mat.consol$mat,
              # P_consol = P.mat.consol$mat,
                Q_consol = Matrix(Q.mat.consol$mat, sparse = TRUE),
                P_consol = Matrix(P.mat.consol$mat, sparse = TRUE),
                Q_freq = t(Q.mat.consol$freq),
                #Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq,
                Q_ind = rep((1:n.real.rows) - 1, real.rows.factor),
                P_ind = rep((1:n.real.rows) - 1, real.rows.factor))
                # P_ind = c(rep((1:n.real.rows) - 1, real.rows.factor)[-1], 0))
))
}

if (T) {
print(system.time(
fisherIndfastestfurious.ret <- fisherIndfastestfurious(
              # Q_consol = Q.mat.consol$mat,
              # P_consol = P.mat.consol$mat,
                Q_consol = Matrix(Q.mat.consol$mat, sparse = TRUE),
                P_consol = Matrix(P.mat.consol$mat, sparse = TRUE),
                Q_freq = Q.mat.consol$freq,
                #Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq,
                Q_ind = rep((1:n.real.rows) - 1, real.rows.factor),
                P_ind = rep((1:n.real.rows) - 1, real.rows.factor))
               #, split_size = 1e+13
                # P_ind = c(rep((1:n.real.rows) - 1, real.rows.factor)[-1], 0))
))
}


fisherEKSdense.ret <- micEconIndex:::fisherEKSdense(
               Q_consol = Q.mat.consol$mat,
               P_consol = P.mat.consol$mat,
                Q_freq = Q.mat.consol$freq,
                #Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq,
                Q_ind = rep((1:n.real.rows) - 1, real.rows.factor),
                P_ind = rep((1:n.real.rows) - 1, real.rows.factor))

fisherInd.ret <- micEconIndex:::fisherInd(Q.mat, P.mat, 1)

summary(fisherEKSdense.ret - fisherInd.ret )




try( print(summary(fisherInd.ret - fisherIndfast.ret)) )

try( print(summary(fisherIndfaster.ret - fisherIndfast.ret)) )

try( print(summary(fisherIndfastest.ret - fisherIndfaster.ret)) )

try( print(summary(fisherIndfastest.ret - fisherIndfastestfurious.ret)) )

try( print(summary(fisherIndfastest.ret - fisherInd.ret)) )




}
*/
