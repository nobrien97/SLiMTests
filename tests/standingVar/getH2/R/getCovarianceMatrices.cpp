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
    if (col_name.substr(0,2) == "VA") {
      // Variances
      std::string var = col_name.substr(3);
      variables_set.insert(var);
    }
    if (col_name.substr(0,3) == "CVA") {
      // Covariances
      size_t first_underscore = col_name.find("_");
      std::string vars = col_name.substr(first_underscore+1);
      size_t second_underscore = vars.find("_");
      std::string var1 = vars.substr(0, second_underscore);
      std::string var2 = vars.substr(second_underscore+1);
      variables_set.insert(var1);
      variables_set.insert(var2);
    }
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
  
  // vector to keep track of which columns have NA values
  std::vector<bool> has_na(n, false);
  
  // Fill the list with covariance matrices for each row
  for (int row = 0; row < num_rows; ++row) {
    // Initialize an empty matrix to store covariances
    NumericMatrix covariance_matrix(n, n);
    
    // Fill the covariance matrix
    for (int i = 0; i < col_names.size(); ++i) {
      std::string col_name = Rcpp::as<std::string>(col_names[i]);
      if (col_name.substr(0, 2) == "VA") {
        // Diagonal elements
        std::string var = col_name.substr(3);
        int idx = var_index[var];
        double value = Rcpp::as<NumericVector>(df[i])[row];
        covariance_matrix(idx, idx) = value;
        
        // If the variance for this variable is NA, mark it for removal later
        if (NumericVector::is_na(value)) {
          has_na[idx] = true;
        }
      } else if (col_name.substr(0, 3) == "CVA") {
        // Off-diagonal elements
        size_t first_underscore = col_name.find("_", 4);
        std::string var1 = col_name.substr(4, first_underscore - 4);
        std::string var2 = col_name.substr(first_underscore + 1);
        int idx1 = var_index[var1];
        int idx2 = var_index[var2];
        double value = Rcpp::as<NumericVector>(df[i])[row];

        covariance_matrix(idx1, idx2) = value;
        covariance_matrix(idx2, idx1) = value; // Because covariance matrices are symmetric
      }
    }
    
    // Remove NA
    std::vector<int> keep_indices;
    for (int i = 0; i < n; ++i) {
      if (!has_na[i]) {
        keep_indices.push_back(i);
      }
    }
    
    // Create a smaller matrix discarding NA values
    int m = keep_indices.size();
    NumericMatrix out_matrix(m, m);
    CharacterVector out_var_names(m);
    
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        out_matrix(i, j) = covariance_matrix(keep_indices[i], keep_indices[j]);
      }
      out_var_names[i] = variables[keep_indices[i]];
    }
    
    // Store the covariance matrix in the list
    covariance_matrices[row] = out_matrix;
    colnames(covariance_matrices[row]) = out_var_names;
    rownames(covariance_matrices[row]) = colnames(covariance_matrices[row]);
  }
  
  return covariance_matrices;
}


