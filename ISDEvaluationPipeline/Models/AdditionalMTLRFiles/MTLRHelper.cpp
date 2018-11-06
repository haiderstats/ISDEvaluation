#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

//Author: Humza Haider
//Email: hshaider@ualberta.ca

//Please see www.cs.cornell.edu/~cnyu/papers/nips11_survival.pdf for relevant information for objective value and gradient calculations.

/*mltr_objVal arguments:
 * (1) params:  The parameters for the MTLR model. These must be sorted with all the biases appearing first followed by the parameters 
 * for the first feature at all time points, second feature for all time points, etc. E.g. with m time points we have:
 * b1, b2,...,bm, theta1,1 ,....theta 1,m, theta 2,1, ... theta 2,m, ... theta p,m
 * (2) yval: A matrix with rows indicating the time points and columns being the patients status at that time point. E.g. for 5 time 
 * intervals and for a patient who died in the third interval we would have 0,0,1,1,1 as a column. yval should be a matrix of all these
 * columns. (Should be N x m where N is the number of patients and m is the number of time points).
 * (3) featureValue: These should be the matrix of feature values. A row indicates a single patient (corresponding to the first COLUMN of yval).
 * (4) C1: The regularzation parameter
 * (5) delta: A vector indicating is a patient is uncensored or censored (delta = 1 --> uncensored).
 */
//Here we will assume that the data is ordered such that all censored observations come first. We could abstract (fairly easily)
//away from this requirement but then we would be doing the same operation multiple times within the bottleneck of the code.

// [[Rcpp::export]]
double mtlr_objVal(NumericVector params, arma::mat yval, arma::mat featureValue, double C1,arma::vec delta) {
  //For the objective value see Equation 3 of the MTLR paper.
  
  //We passed params in as a NumericVector because these are easier to parse and turn into a matrix.
  double valToReturn = 0;
  int N = featureValue.n_rows;
  int m = yval.n_rows;
  Range biasIdx(0,m-1);
  
  NumericVector biasVec = params[biasIdx];
  NumericMatrix thetaMatrix(m,featureValue.n_cols);
  int thetaSize = thetaMatrix.nrow() * thetaMatrix.ncol();
  for (int i = 0; i < thetaSize; i++) {
    thetaMatrix[i] = params[m+i];
  }  
  
  arma::mat thetas = as<arma::mat>(thetaMatrix);
  arma::vec biases = as<arma::vec>(biasVec);
  
  arma::vec val(m,arma::fill::ones);
  arma::vec bigVal(m, arma::fill::ones);
  arma::vec secondPieceVal(m, arma::fill::ones);
  arma::vec firstPieceVal(m, arma::fill::ones);
  arma::uvec uncensoredIndex;
  //Censored Piece
  double NCens = sum(1-delta);
  double firstPiece;
  double biggestCensVal;
  double biggestVal;
  //For the censored individuals (see Equation 4)... and take the log to get the log-likelihood.
  for(int i=0;i<NCens;i++){
    val =  thetas * featureValue.row(i).t() + biases;
    bigVal  = reverse(cumsum(reverse(val)));
    // log trick: https://timvieira.github.io/blog/post/2014/02/11/exp-normalize-trick/
    biggestVal = bigVal.max() > 0 ? bigVal.max() : 0;
    secondPieceVal = bigVal - biggestVal;
    double secondPiece = biggestVal + log(sum(exp(secondPieceVal))+exp(-biggestVal));
    
    uncensoredIndex = find(yval.col(i) == 1);
    
    if(uncensoredIndex.size() ==0){
      firstPiece = 0;
    }else{
      biggestCensVal  = bigVal.elem(uncensoredIndex).max() > 0 ? bigVal.elem(uncensoredIndex).max() : 0 ;
      firstPieceVal = bigVal - biggestCensVal;
      //we multiply firstPieceVal by ycol to 0 out the other large values in bigVal. Note while this will make exp(0) = 1, it will again be zeroed out by multiplying by y.col(i) again.
      firstPieceVal %= yval.col(i);
      firstPiece = biggestCensVal + log(sum(yval.col(i) % exp(firstPieceVal)) + exp(-biggestCensVal));
    }

    valToReturn += firstPiece - secondPiece;
    
  }
  //For the uncensored individuals.
  for(int i=NCens;i<N;i++){
    val =  thetas * featureValue.row(i).t() + biases;
    bigVal  = reverse(cumsum(reverse(val)));
    double biggestVal  = bigVal.max();
    bigVal -= biggestVal;
    double secondPiece = biggestVal + log(sum(exp(bigVal))+exp(-biggestVal));
    double firstPiece = sum(yval.col(i).t() * val);

    valToReturn += firstPiece - secondPiece;
  }
  
  valToReturn = (C1/2)*accu(square(thetas)) - valToReturn/N;
  return valToReturn;
}


//See mtlr_objVal for argument definition.
//Here we will assume that the data is ordered such that all censored observations come first. We could abstract (fairly easily)
//away from this requirement but then we would be doing the same operation multiple times within the bottleneck of the code.
// [[Rcpp::export]]
arma::rowvec mtlr_grad(NumericVector params, arma::mat yval, arma::mat featureValue, double C1, arma::vec delta) {
  double N = featureValue.n_rows;
  double P = featureValue.n_cols;
  int m = yval.n_rows;
  Range biasIdx(0,m-1);
  
  NumericVector biasVec = params[biasIdx];
  NumericMatrix thetaMatrix(m,P);
  int thetaSize = thetaMatrix.nrow() * thetaMatrix.ncol();
  for (int i = 0; i < thetaSize; i++) {
    thetaMatrix[i] = params[m+i];
  }  
  
  arma::mat thetas = as<arma::mat>(thetaMatrix);
  arma::vec biases = as<arma::vec>(biasVec);
  arma::vec BiasGradient = arma::zeros<arma::vec>(biases.n_elem);
  arma::mat ThetaGradient = arma::zeros<arma::mat>(thetas.n_rows, thetas.n_cols);
  
  arma::rowvec Xval(P, arma::fill::ones);
  arma::mat  Numer(m,1,arma::fill::ones), ThetaSecond(m,P,arma::fill::ones), toAddWeight(m,P,arma::fill::ones);
  arma::vec toAddBias(m,arma::fill::ones), val(m, arma::fill::ones),bigVal(m, arma::fill::ones),bigValCens(m, arma::fill::ones), BiasSecond(m,arma::fill::ones);
  arma::uvec uncensoredIndex;
  double Denom;
  
  //Censored Piece
  arma::vec BiasCens(m,arma::fill::ones);
  arma::mat NumerCens(m,1,arma::fill::ones), ThetaCens(m,P,arma::fill::ones);
  double DenomCens;
  double NCens = sum(1-delta);
  double biggestCensVal;
  for(int i=0;i<NCens;i++){
    Xval = featureValue.row(i);
    val =  (thetas * Xval.t()) + biases;
    bigVal  = reverse(cumsum(reverse(val)));
    
    uncensoredIndex = find(yval.col(i) == 1);
    
    if(uncensoredIndex.size() ==0){
      bigValCens = bigVal % yval.col(i);
      biggestCensVal = 0;
    }else{
      biggestCensVal  = bigVal.elem(uncensoredIndex).max() > 0 ? bigVal.elem(uncensoredIndex).max() : 0 ;

      bigValCens = bigVal - biggestCensVal;
      //we multiply firstPieceVal by ycol to 0 out the other large values in bigVal. Note while this will make exp(0) = 1, it will again be zeroed out by multiplying by y.col(i) again.
      bigValCens %= yval.col(i);
    }
    
    NumerCens= cumsum(yval.col(i) % exp(bigValCens));
    DenomCens = NumerCens[NumerCens.n_rows-1] + exp(-biggestCensVal);

    BiasCens = NumerCens/DenomCens;
    
    NumerCens *= Xval;
    
    ThetaCens = NumerCens/DenomCens;
    
    double biggestVal;
    biggestVal = bigVal.max() > 0 ? bigVal.max() : 0;
  
    bigVal -= biggestVal;
    
    Numer= cumsum(exp(bigVal));
    Denom = Numer[Numer.n_rows-1] +exp(-biggestVal);

    BiasSecond = Numer/Denom;
    
    Numer *= Xval;
    ThetaSecond = Numer/Denom;
    
    toAddBias =  BiasCens - BiasSecond;
    toAddWeight =  ThetaCens - ThetaSecond;
    BiasGradient += toAddBias;

    ThetaGradient += toAddWeight;
  }
  
  //Uncensored Piece
  for(int i=NCens;i<N;i++){
    Xval = featureValue.row(i);
    val =  (thetas * Xval.t()) + biases;
    
    bigVal  = reverse(cumsum(reverse(val)));
    double biggestVal;
    biggestVal = bigVal.max() > 0 ? bigVal.max() : 0;
    bigVal -= biggestVal;
    
    Numer= cumsum(exp(bigVal));
    Denom = Numer[Numer.n_rows-1] +exp(-biggestVal);
    BiasSecond = Numer/Denom;
    
    Numer *= Xval;
    ThetaSecond = Numer/Denom;
    
    toAddBias =  sum(yval.col(i) - BiasSecond.each_col(), 1);
    toAddWeight =  yval.col(i) * Xval - ThetaSecond;
    BiasGradient += toAddBias;
    ThetaGradient += toAddWeight;
  }
  
  BiasGradient.transform( [N](double val) { return (val * (double)(-1/N)); } );
  ThetaGradient.transform( [N](double val) { return (val * (double)(1/N)); } );
  
  ThetaGradient = C1*thetas - ThetaGradient;
  arma::vec ThetaFlat = vectorise(ThetaGradient);
  return arma::join_cols<arma::mat>(BiasGradient, ThetaFlat).t();
}


//See mltr_objVal for argument definitions.
// [[Rcpp::export]]
arma::mat mtlr_predict(NumericVector params,arma::mat featureValue){
  double N = featureValue.n_rows;
  double P = featureValue.n_cols;
  
  
  double m = params.size()/(P+1); // +1 to account for bias terms.
  Range biasIdx(0,m-1);
  
  NumericVector biasVec = params[biasIdx];
  NumericMatrix thetaMatrix(m,featureValue.n_cols);
  int thetaSize = thetaMatrix.nrow() * thetaMatrix.ncol();
  for (int i = 0; i < thetaSize; i++) {
    thetaMatrix[i] = params[m+i];
  }  
  arma::mat thetas = as<arma::mat>(thetaMatrix);
  arma::vec biases = as<arma::vec>(biasVec);
  
  
  arma::mat resMat = thetas*featureValue.t();
  resMat.each_col() += biases;
  arma::mat B = exp(reverse(cumsum(reverse(resMat))));
  arma::mat Ones = arma::ones<arma::mat>(1,N);
  
  B.insert_rows(B.n_rows,Ones);
  arma::rowvec BDenoms = sum(B);
  arma::mat resultsMatrix = reverse(cumsum(reverse(B.each_row()/BDenoms)));
  resultsMatrix.shed_row(0);
  return resultsMatrix;
}













