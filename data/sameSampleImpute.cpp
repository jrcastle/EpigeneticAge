#include <Rcpp.h>
#include <cmath> 
#include <cstdlib>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
void sameSampleImpute(CharacterVector badCpGList, IntegerVector badCpGRows_, DataFrame MethDfIn_, int nearby_limit_, bool weighted_ = true) {
  
  CharacterVector colnames = MethDfIn_.names();
  int NCpGs = badCpGList.size();
  for(int i = 0; i < NCpGs; i++){
    
    if( i % 5000 == 0 ) Rcout << "Imputing CpG " << i+1 << "/" << NCpGs << "\t(" << (int)(((float)i+1.)/(float)NCpGs) << "% complete)" << endl;
    int bad_cpg_row = badCpGRows_[i]-1;
    
    // Get Chr and bpPos of bad CpG
    IntegerVector chr    = MethDfIn_["chr"];
    IntegerVector bppos  = MethDfIn_["bppos"];
    int bad_cpg_chr      = chr[bad_cpg_row];
    int bad_cpg_bppos    = bppos[bad_cpg_row];
    
    // Identify near-neighbor CpGs based on distance from bad cpg
    std::vector<int> neighbor_rows;
    
    // Search rows prior to the bad row to find neighbors
    int irow = bad_cpg_row - 1;
    if( irow >= 0 ){
      while( std::fabs(bppos[irow] - bad_cpg_bppos) <= nearby_limit_ ){

        // Don't consider rows not on the same chromosome 
        int irow_chr = chr[irow];
        if(irow_chr != bad_cpg_chr) break;
        
        // Drop rows located > NEARBY.LIMIT away
        int irow_bppos = bppos[irow];
        int sep        = std::fabs(irow_bppos - bad_cpg_bppos);
        if(sep > nearby_limit_) break;
        
        // All tests, passed, keep this row number
        neighbor_rows.push_back(irow);
        irow--;
        if( irow < 0 ) break;
      }
    }
    
    // Search rows after the bad row to find neighbors
    irow = bad_cpg_row + 1;
    if( irow >= 0 ){
      while( std::fabs(bppos[irow] - bad_cpg_bppos) <= nearby_limit_ ){
        
        // Don't consider rows not on the same chromosome 
        int irow_chr = chr[irow];
        if(irow_chr != bad_cpg_chr) break;
        
        // Drop rows located > NEARBY.LIMIT away
        int irow_bppos = bppos[irow];
        int sep        = std::fabs(irow_bppos - bad_cpg_bppos);
        if(sep > nearby_limit_) break;
        
        // All tests, passed, keep this row number
        neighbor_rows.push_back(irow);
        irow++;
        if( irow > chr.size()-1 ) break;
      }
    }
    
    // Check to see if any neighbors were found, skip CpG if none were found
    if( neighbor_rows.size() == 0 ){
      Rcout << "WARNING: no neighbors found for " << badCpGList[i] << " ... skipping imputation" << endl;
      continue;
    }
    
    // IMPUTE
    for(int icol = 0; icol < colnames.size(); icol++){
      
      // Get sample column
      string colname = Rcpp::as<std::string>(colnames[icol]);
      if(colname == "position" | colname == "chr" | colname == "bppos") continue;
      DoubleVector SAMPLE = MethDfIn_[colname];
      
      // If the value for this sample is not missing, don't impute it
      if( !NumericVector::is_na(SAMPLE[bad_cpg_row]) ) continue;
      
      // Mean imputation
      double imputed_value = 0.;
      double sumw          = 0.;
      for(int j = 0; j < neighbor_rows.size(); j++){
        int row = neighbor_rows[j];
        if( !NumericVector::is_na(SAMPLE[row]) ){
          if(weighted_){
            double w = 1./(double)std::fabs(bppos[row] - bad_cpg_bppos);
            imputed_value += SAMPLE[row] * w;
            sumw          += w;
          }else{
            imputed_value += SAMPLE[row];
            sumw          += 1.;
          }
        }
      }
      
      if( sumw == 0. ){
        Rcout << "WARNING: " << badCpGList[i] << " cannot be imputed for sample " << colname << ": All neigbor beta values are NA" << endl;
        continue;
      }
      
      imputed_value /= sumw;
      SAMPLE[bad_cpg_row] = imputed_value;
    }
  }
}
