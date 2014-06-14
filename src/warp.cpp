////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////  Weak Axiom of Revealed Preferences (WARP) functions  //////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
// check WARP
RcppExport SEXP CheckWarp(SEXP x, SEXP p) { try {
  // import quantities and prices matrices (x, p)
  NumericMatrix r_x(x), r_p(p);
  unsigned n_obs= r_x.nrow();
  arma::mat mat_x(r_x.begin(), r_x.nrow(), r_x.ncol());
  arma::mat mat_p(r_p.begin(), r_p.nrow(), r_p.ncol());
  
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= i_row+1; i_col < n_obs; i_col++) {
      if (arma::dot(mat_p.row(i_row), mat_x.row(i_row)) >=
            arma::dot(mat_p.row(i_row), mat_x.row(i_col))) {
        if (arma::dot(mat_p.row(i_col), mat_x.row(i_col)) >=
              arma::dot(mat_p.row(i_col), mat_x.row(i_row))) {
          if (arma::accu(arma::abs(mat_p.row(i_col) - mat_x.row(i_row))) != 0)
          { // if p_i x_i >= p_i x_k AND p_k x_k >= p_k x_i AND x_i != x_k
            NumericVector violators(2);
            violators(0)= i_row + 1;
            violators(1)= i_col + 1;
            return List::create(Named("violation", wrap(true)), 
                                Named("violators", violators));
          }
        }
      }
    }
  }
  // if no violation found
  return List::create(Named("violation", wrap(false)));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}
