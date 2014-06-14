#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
// Compute all indirect preferences
// prefs(i, j)= 0 iff i not indirectly prefered to j
// prefs(i, j)= 1 iff i indirectly prefered to j (only equalities)
// prefs(i, j)= 2 iff i indirectly strictly prefered to j 
//                    (with at least one strict preference)
RcppExport SEXP IndirectPrefs(SEXP px) { try {
  // import prices*quantities matrix (px)
  NumericMatrix r_px(px);
  unsigned n_obs= r_px.nrow();
  arma::mat mat_px(r_px.begin(), n_obs, n_obs);
  
  // Compute direct revealed preferences from matrix p*x
  arma::mat direct_prefs(arma::zeros(n_obs, n_obs));
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= 0; i_col < n_obs; i_col++) {
      if (mat_px(i_row, i_row) > mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 2;
      } else if (mat_px(i_row, i_row) == mat_px(i_row, i_col)) { 
        direct_prefs(i_row, i_col)= 1;
      }
    }
  }
  
  // Variant of Warshall-Floyd algorithm to all find indirect preferences 
  arma::mat indirect(direct_prefs);
  for (unsigned k= 0; k < n_obs; k++) {
    for (unsigned i_row= 0; i_row < n_obs; i_row++) {
      for (unsigned i_col= 0; i_col < n_obs; i_col++) {
        if (indirect(i_row, i_col) < 2) {
          if(indirect(i_row, k) != 0 && indirect(k, i_col) != 0) {
            if(indirect(i_row, k) == 2 || indirect(k, i_col) == 2) {
              indirect(i_row, i_col)= 2;
            } else indirect(i_row, i_col)= 1;
          }
        }
      }
    }
  }
  return wrap(indirect);  

} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}
