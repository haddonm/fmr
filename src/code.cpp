#include <Rcpp.h>

using namespace Rcpp;


//' @title calccaaC calculates the predicted catch-at-age
//'
//' @description calccaaC calculate the expected catch-at-age for a 
//'     statistical catch-at-age model. Implemented in Rcpp.
//'
//' @param pars a vector of log-transformed parameters
//' @param pnaa the predicted numbers-at-age in the population
//' @param M the natural mortality instantaneous rate
//' @param sel the selectivity at age
//' @param ages a vector identifying the ages being considered
//' @export
// [[Rcpp::export]]
NumericMatrix calccaaC(NumericVector pars, NumericMatrix pnaa, double M,
                      NumericVector sel, NumericVector ages) {
  int nage = ages.size();
  int nyr = pnaa.nrow();
  int f1 = nyr + nage - 1;
  double sF;
  NumericMatrix caa = pnaa;
  for (int yr = 0; yr < nyr; yr++) { 
    for (int age = 0; age< nage; age++) {
      sF = sel[age] * exp(pars[(yr + f1)]);
      caa(yr,age) = (sF/(M + sF)) * pnaa(yr,age) * (1 - exp(-(M + sF)));
    }
  }
  return caa;
}

//' @title calcnaaC calculates the predicted numbers-at-age and exploitable B
//'
//' @description calcnaaC calculates the expected numbers-at-age, the 
//'     exploitable-numbers-at-age, the exploitable biomass by year, and the
//'     selectivity for a statistical catch-at-age model. Implemented in Rcpp.
//'
//' @param pars a vector of log-transformed parameters
//' @param M the natural mortality instantaneous rate
//' @param yrs is a vector of the years
//' @param ages a vector identifying the ages being considered
//' @param owa is a matrix of the observed weight-at-age
//' @export
// [[Rcpp::export]]
List calcnaaC(NumericVector pars, double M, NumericVector yrs,
                       NumericVector ages, NumericMatrix owa) { 
  int nyr = yrs.size();
  int nage = ages.size();
  int f1 = nyr + nage - 2;
  NumericMatrix pnaa(nyr,nage);
  NumericMatrix enaa(nyr,nage);
  NumericVector exB(nyr);
  NumericVector sel(nage);
  double sF;
  double log19 = log(19);
  double L50 = exp(pars[26]);
  double L95 = exp(pars[27]);
  for (int age = 0; age < nage; age++) {
    sel[age] = 1/(1 + exp(-log19 * (ages[age] - L50)/(L95 - L50)));
  }
  for (int yr = 0;yr < nyr; yr++) {
    pnaa(yr,0) = exp(pars[yr]);
  }
  for (int age = 1;age< nage; age++) {
    pnaa(0,age) = exp(pars[(nyr+age-1)]);
  }
  for (int yr = 1;yr < nyr;yr++) {
    for (int age = 1;age < nage; age++) {
      sF = sel[age-1] * exp(pars[yr + f1]);
      pnaa(yr,age) = pnaa((yr-1),(age-1)) * exp(-(M + sF));
    }
  }
  for (int yr = 0; yr< nyr; yr++) {
    for (int age= 0; age < nage; age++) {
      enaa(yr,age) = pnaa(yr,age) * sel[age];
      exB[yr] += (enaa(yr,age) * owa(yr,age));
    }
    exB[yr] = exB[yr]/1000.0;
  }
  List out = List::create(
    _["pnaa"]=pnaa,
    _["enaa"]=enaa,
    _["exB"]=exB,
    _["sel"]=sel);

  return out;
}










  