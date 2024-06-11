#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List extractCovarianceMatrices(DataFrame df) {
  // Extract column names
  CharacterVector col_names = df.names();
  int num_rows = df.nrows();
  
  // Initialize an empty set to store unique variable names
  std::set<std::string> variables_set;
  
  // Iterate over column names to extract unique variable names
  for (int i = 0; i < col_names.size(); ++i) {
    std::string col_name = as<std::string>(col_names[i]);
    size_t pos = col_name.find("_");
    std::string var1 = col_name.substr(0, pos);
    std::string var2 = col_name.substr(pos + 1);
    variables_set.insert(var1);
    variables_set.insert(var2);
  }
  
  // Convert the set of variables to a vector
  std::vector<std::string> variables(variables_set.begin(), variables_set.end());
  int n = variables.size();
  
  // Create a map to quickly get the index of a variable
  std::map<std::string, int> var_index;
  for (int i = 0; i < n; ++i) {
    var_index[variables[i]] = i;
  }
  
  // Initialize a list to store the covariance matrices
  List covariance_matrices(num_rows);
  
  // Fill the list with covariance matrices for each row
  for (int row = 0; row < num_rows; ++row) {
    // Initialize an empty matrix to store covariances
    NumericMatrix covariance_matrix(n, n);
    
    // Fill the covariance matrix
    for (int i = 0; i < col_names.size(); ++i) {
      std::string col_name = Rcpp::as<std::string>(col_names[i]);
      size_t pos = col_name.find("_");
      std::string var1 = col_name.substr(0, pos);
      std::string var2 = col_name.substr(pos + 1);
      int idx1 = var_index[var1];
      int idx2 = var_index[var2];
      double value = Rcpp::as<NumericVector>(df[i])[row];
      covariance_matrix(idx1, idx2) = value;
      covariance_matrix(idx2, idx1) = value; // Because covariance matrices are symmetric
    }
    
    // Fill diagonal with variances (covariance with itself)
    for (int i = 0; i < n; ++i) {
      std::string var = variables[i];
      std::string diag_col = var + "_" + var;
      if (df.containsElementNamed(diag_col.c_str())) {
        covariance_matrix(i, i) = Rcpp::as<NumericVector>(df[diag_col])[row];
      }
    }
    // Store the covariance matrix in the list
    covariance_matrices[row] = covariance_matrix;
  }
  
  return covariance_matrices;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
