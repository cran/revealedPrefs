////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////  Generalized Axiom of Revealed Preferences (GARP) functions  /////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*  
Copyright 2014 Julien Boelaert.

This file is part of revealedPrefs.

revealedPrefs is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

revealedPrefs is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with revealedPrefs.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
// Auxiliary functions

// Revealed Prefs from matrices x and p, for obs i and k
// afriat_par in (0,1), 1 for standard RP
// return 0 if p_i x_i < p_i x_k,
// return 1 if p_i x_i == p_i x_k,
// return 2 if p_i x_i > p_i x_k,
unsigned CheckRevPref(arma::mat *x, arma::mat *p, unsigned obs_i,
                      unsigned obs_k, double afriat_par) {
  double d_test= afriat_par * arma::dot(p[0].row(obs_i), x[0].row(obs_i)) -
    arma::dot(p[0].row(obs_i), x[0].row(obs_k));
  if (d_test < 0) return 0;
  if (d_test > 0) return 2;
  return 1;
}

// boolean check GARP from an arma::mat argument (used by CpUp)
bool ViolateGarp(arma::mat mat_px, double afriat_par) {
  unsigned n_obs= mat_px.n_rows;
  arma::mat direct_prefs(arma::zeros(n_obs, n_obs)), 
    direct_strict(arma::zeros(n_obs, n_obs));
    
  // Compute direct (and strict) revealed preferences from matrix p*x
  // exit as soon as contradiction between direct preferences (WARP violation)
  // pref(i,j)=1 iff bundle i directly prefered to bundle j
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= i_row; i_col < n_obs; i_col++) {
      // i_row prefered to i_col?
      if (afriat_par * mat_px(i_row, i_row) > mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 1;
      } else if (afriat_par * mat_px(i_row, i_row) == mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 0;
      }

      // i_col prefered to i_row?
      if (afriat_par * mat_px(i_col, i_col) > mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 1;
      } else if (afriat_par * mat_px(i_col, i_col) == mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 0;
      }
      
      // violation in direct preferences?
      if ((direct_prefs(i_row, i_col) + direct_strict(i_col, i_row) == 2) ||
            (direct_prefs(i_col, i_row) + direct_strict(i_row, i_col) == 2)) {
        return true;
      }
    }
  }

  // Variant of Warshall's algorithm to all find indirect preferences 
  // Exit as soon as a contradiction with direct strict preferences is found
  // (ie GARP violation)
  arma::mat indirect(direct_prefs);
  for (unsigned k= 0; k < n_obs; k++) {
    for (unsigned i_row= 0; i_row < n_obs; i_row++) {
      for (unsigned i_col= 0; i_col < n_obs; i_col++) {
        if (indirect(i_row, i_col) == 0) {
          if(indirect(i_row, k) != 0 && indirect(k, i_col) != 0) {
            indirect(i_row, i_col)= 1;
            if (direct_strict(i_col, i_row) != 0) {
              return true;
            }
          }
        }
      }
    }
  }
  
  // if no violation found
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// GARP - Variant of Floyd-Warshall algorithm
RcppExport SEXP CheckGarp(SEXP px, SEXP afriat) { try {
  // import prices*quantities matrix (px)
  NumericMatrix r_px(px);
  unsigned n_obs= r_px.nrow();
  arma::mat mat_px(r_px.begin(), n_obs, n_obs);
  double afriat_par= as<double>(afriat);
  
  // Compute direct (and strict) revealed preferences from matrix p*x
  // exit as soon as contradiction between direct preferences
  // pref(i,j)=1 iff bundle i directly prefered to bundle j
  arma::mat direct_prefs(arma::zeros(n_obs, n_obs)), 
    direct_strict(arma::zeros(n_obs, n_obs));
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= i_row; i_col < n_obs; i_col++) {
      // i_row prefered to i_col?
      if (afriat_par * mat_px(i_row, i_row) > mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 1;
      } else if (afriat_par * mat_px(i_row, i_row) == mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 0;
      }

      // i_col prefered to i_row?
      if (afriat_par * mat_px(i_col, i_col) > mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 1;
      } else if (afriat_par * mat_px(i_col, i_col) == mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 0;
      }
      
      // violation in direct preferences?
      if ((direct_prefs(i_row, i_col) + direct_strict(i_col, i_row) == 2) ||
            (direct_prefs(i_col, i_row) + direct_strict(i_row, i_col) == 2)) {
        NumericVector violators(2), strict(2);
        violators(0)= i_row + 1; violators(1)= i_col + 1;
        strict(0)= direct_strict(i_row, i_col); 
        strict(1)= direct_strict(i_col, i_row); 
        return List::create(Named("violation", wrap(true)),
                            Named("violators", violators), 
                            Named("strict", strict),
                            Named("direct.violation", wrap(true)));
      }
    }
  }
  
  // Variant of Warshall's algorithm to all find indirect preferences 
  // Exit as soon as a contradiction with direct strict preferences is found
  // (ie GARP violation)
  arma::mat indirect(direct_prefs);
  for (unsigned k= 0; k < n_obs; k++) {
    for (unsigned i_row= 0; i_row < n_obs; i_row++) {
      for (unsigned i_col= 0; i_col < n_obs; i_col++) {
        if (indirect(i_row, i_col) == 0) {
          if(indirect(i_row, k) != 0 && indirect(k, i_col) != 0) {
            indirect(i_row, i_col)= 1;
            if (direct_strict(i_col, i_row) != 0) {
              NumericVector violators(2), strict(2);
              violators(0)= i_row + 1; violators(1)= i_col + 1;
              strict(0)= 0; strict(1)= 1;
              return List::create(Named("violation", wrap(true)),
                                  Named("violators", violators),
                                  Named("strict", strict),
                                  Named("direct.violation", wrap(false))) ;
            }
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

////////////////////////////////////////////////////////////////////////////////
// GARP - recursive function for depth-first search with tabu list
// returns empty vector if no violation found
// returns vector of violating path if violation found
// (used by DeepGarp)
std::vector<unsigned> RecGarp(unsigned cur_obs, std::vector<bool> * tabu, 
                              unsigned *n_tabu, 
                              std::vector<unsigned> *hist_tabu,
                              std::vector<unsigned> ascendence, 
                              std::vector<bool> strict_asc, 
                              arma::mat *x, arma::mat *p, 
                              double afriat_par) {
  bool b_nonstrict_cycle= false;
  for (unsigned i_search= 0; i_search < x[0].n_rows; i_search++) {
    if (!tabu[0][i_search] && i_search!=cur_obs) {
      unsigned rev_pref= CheckRevPref(x, p, cur_obs, i_search, afriat_par);
      if (rev_pref) {
        // if cur_obs prefered to i_search, check if i_search is in ascendence
        bool b_asc= false;
        unsigned i_asc;
        for (i_asc= 0; i_asc < ascendence.size(); i_asc++) {
          if (ascendence[i_asc] == i_search) {
            b_asc= true;
            break;
          }
        }
        if (b_asc) {
          // if i_search is in ascendence, check if there is a strict ascendence
          // (start looking where it left off, at beginning of loop)
          bool b_strict= false;
          for (; i_asc < ascendence.size(); i_asc++) {
            if (strict_asc[i_asc]) {
              b_strict= true;
              break;
            }
          }
          if (b_strict) { // found strict cycle, exit
            ascendence.push_back(cur_obs);
            ascendence.push_back(i_search);
            return(ascendence);
          }
          // if non-strict cycle, don't pursue into past but don't mark tabu
          b_nonstrict_cycle= true;
        } else {
          // if not in ascendence, pursue depth-first search
          std::vector<unsigned> new_asc= ascendence;
          std::vector<bool> new_strict= strict_asc;
          new_asc.push_back(cur_obs);
          if (rev_pref == 2) {
            new_strict.push_back(true);
          } else new_strict.push_back(false);
          
          std::vector<unsigned> pursue= RecGarp(i_search, 
                                                tabu, n_tabu, hist_tabu,
                                                new_asc, new_strict, x, p, 
                                                afriat_par);
          if (pursue.size())
            return pursue; // if violation in descendence, return violation
        }
      }
    }
  }
  
  // once all the available choices have been tried,
  // if no nonstrict cycle found, tabu the current obs
  if (!b_nonstrict_cycle) {
    tabu[0][cur_obs]= true;
    (*n_tabu)++;
    hist_tabu[0].push_back(cur_obs);
  }
  return std::vector<unsigned>(0);
}

////////////////////////////////////////////////////////////////////////////////
// Check GARP using depth-first search with tabu list
RcppExport SEXP DeepGarp(SEXP quanti, SEXP prices, SEXP afriat) { try {
  // import arguments
  NumericMatrix r_quanti(quanti), r_prices(prices);
  unsigned n_obs= r_quanti.nrow();
  arma::mat mat_q(r_quanti.begin(), r_quanti.nrow(), r_quanti.ncol());
  arma::mat mat_p(r_prices.begin(), r_prices.nrow(), r_prices.ncol());
  double afriat_par= as<double>(afriat);
  
  // tabu list to be passed to the recursive search
  std::vector<bool> tabu(n_obs, false);
  unsigned n_tabu= 0;
  std::vector<unsigned> hist_tabu; // history of tabu obs, ie. utility ordering
  
  std::vector<unsigned> df_search;
  bool found_violation= false;
  for (unsigned current_obs=0; !found_violation && n_tabu < n_obs ; current_obs++) {
    if (!tabu[current_obs]) {
      // launch recursive search with a non-tabu observation
      df_search= RecGarp(current_obs, &tabu, &n_tabu, &hist_tabu,
                         std::vector<unsigned>(0), 
                         std::vector<bool>(0),
                         &mat_q, &mat_p,
                         afriat_par);
      if (df_search.size())
        found_violation= true;
    }
  }
  
  if (found_violation) {
    // establish strictness of preferences on the violating path
    std::vector<bool> strict;
    for (unsigned i_strict= 0; i_strict < df_search.size() - 1; i_strict++) {
      unsigned rev_pref= CheckRevPref(&mat_q, &mat_p, 
                                      df_search[i_strict], 
                                      df_search[i_strict + 1], 
                                      afriat_par);
      if (rev_pref == 2) {
        strict.push_back(true);
      } else strict.push_back(false);
    }
    return List::create(Named("violation", wrap(true)),
                        Named("path", wrap(df_search)),
                        Named("path.strict", wrap(strict)));
  }
  return List::create(Named("violation", wrap(false)), 
                      Named("pref.order", wrap(hist_tabu)));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}
