#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int ruin_prob_cpp(NumericVector x, int n, int m) {
  int count = 0;
  
  for(int j = 0; j < m; ++j) {
    
    double sum = 0;
    
    for(int i = 0; i < n; ++i) {
      
      sum += x[n * j + i];
      
      if(30 + sum <= 0) {
        count++;
        break;
      }
      
    }
  }
  
  return count;
}
