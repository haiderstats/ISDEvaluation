#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

//params will be: b1, b2,...,bm, theta1,1 ,....theta 1,m, theta 2,1, ... theta 2,m, ... theta p,m
// [[Rcpp::export]]
double objValC_Cens(NumericVector params, arma::mat yval, arma::mat featureValue, double C1,arma::vec delta) {
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
  
  //Censored Piece
  double NCens = sum(1-delta);
  
  for(int i=0;i<NCens;i++){
    val =  thetas * featureValue.row(i).t() + biases;
    double secondPiece = log(sum(exp(reverse(cumsum(reverse(val)))))+1);
    double firstPiece = log(sum(yval.col(i) % exp(reverse(cumsum(reverse(val)))))+1);

    valToReturn += firstPiece - secondPiece;
  }
  
  for(int i=NCens;i<N;i++){
    val =  thetas * featureValue.row(i).t() + biases;
    double secondPiece = log(sum(exp(reverse(cumsum(reverse(val)))))+1);
    double firstPiece = sum(yval.col(i).t() * val);
    valToReturn += firstPiece - secondPiece;
  }
  
  valToReturn = (C1/2)*accu(square(thetas)) - valToReturn/N;
  return valToReturn;
}


// [[Rcpp::export]]
double objValC_Cens_LogTrick(NumericVector params, arma::mat yval, arma::mat featureValue, double C1,arma::vec delta) {
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


//Here we will assume that the data is ordered such that all censored observations come first. We could abstract (fairly easily)
//away from this requirement but then we would be doing the same operation multiple times within the bottleneck of the code.
// [[Rcpp::export]]
arma::rowvec gradC_Cens(NumericVector params, arma::mat yval, arma::mat featureValue, double C1, arma::vec delta) {
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
  arma::vec toAddBias(m,arma::fill::ones), val(m, arma::fill::ones), BiasSecond(m,arma::fill::ones);
  double Denom;
  
  //Censored Piece
  arma::vec BiasCens(m,arma::fill::ones);
  arma::mat NumerCens(m,1,arma::fill::ones), ThetaCens(m,P,arma::fill::ones);
  double DenomCens;
  double NCens = sum(1-delta);
  for(int i=0;i<NCens;i++){
    Xval = featureValue.row(i);
    val =  (thetas * Xval.t()) + biases;
    
    NumerCens= cumsum(yval.col(i) % exp(reverse(cumsum(reverse(val)))));

    DenomCens = NumerCens[NumerCens.n_rows-1] +1;
    

    BiasCens = NumerCens/DenomCens;
    
    NumerCens *= Xval;
    
    ThetaCens = NumerCens/DenomCens;
    

    Numer= cumsum(exp(reverse(cumsum(reverse(val)))));
    Denom = Numer[Numer.n_rows-1] +1;
    
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
    Numer= cumsum(exp(reverse(cumsum(reverse(val)))));
    Denom = Numer[Numer.n_rows-1] +1;
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

// [[Rcpp::export]]
arma::rowvec gradC_Cens_LogTrick(NumericVector params, arma::mat yval, arma::mat featureValue, double C1, arma::vec delta) {
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


// [[Rcpp::export]]
arma::mat predictC(NumericVector params,arma::mat featureValue){
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




/*
 // [[Rcpp::export]]
 double objValC(NumericVector params, arma::mat yval, arma::mat featureValue, double C1) {
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

arma::vec val;
for(int i=0;i<N;i++){
val =  thetas * featureValue.row(i).t() + biases;
double secondPiece = log(sum(exp(reverse(cumsum(reverse(val)))))+1);
double firstPiece = sum(yval.col(i).t() * val);
valToReturn += firstPiece - secondPiece;
}
valToReturn = (C1/2)*accu(square(thetas)) - valToReturn/N;
return valToReturn;
}


// [[Rcpp::export]]
arma::rowvec gradC(NumericVector params, arma::mat yval, arma::mat featureValue, double C1) {
double N = featureValue.n_rows;
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
arma::vec BiasGradient = arma::zeros<arma::vec>(biases.n_elem);
arma::mat ThetaGradient = arma::zeros<arma::mat>(thetas.n_rows, thetas.n_cols);

arma::rowvec Xval;
arma::mat  Numer, BiasSecond, ThetaSecond, toAddWeight;
arma::vec toAddBias, val;
double Denom;
for(int i=0;i<N;i++){
Xval = featureValue.row(i);
val =  (thetas * Xval.t()) + biases;
Numer= cumsum(exp(reverse(cumsum(reverse(val)))));
Denom = Numer[Numer.n_rows-1] +1;
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


 
 
 // [[Rcpp::export]]
 double objValC_Cens_Matrix(NumericVector params, arma::mat yval, arma::mat featureValue, double C1,arma::vec delta) {
 
 double valToReturn = 0;
 int N = featureValue.n_rows;
 int m = yval.n_rows;
 
 //Get biases and weights into armadillo vector and matrix.
 Range biasIdx(0,m-1);
 NumericVector biasVec = params[biasIdx];
 NumericMatrix thetaMatrix(m,featureValue.n_cols);
 int thetaSize = thetaMatrix.nrow() * thetaMatrix.ncol();
 for (int i = 0; i < thetaSize; i++) {
 thetaMatrix[i] = params[m+i];
 }  
 
 arma::mat thetas = as<arma::mat>(thetaMatrix);
 arma::vec biases = as<arma::vec>(biasVec);
 
 //Do a shared operation and break it into an uncensored piece and censored piece.
 double NCens = sum(1-delta);
 arma::mat valMat = thetas * featureValue.t();
 valMat.each_col() += biases;
 arma::mat CensMat = valMat.cols(0,NCens-1);
 arma::mat UncensMat = valMat.cols(NCens,N-1);
 
 //Censored Piece
 arma::mat CensIntermed = exp(reverse(cumsum(reverse(CensMat))));
 arma::mat secondMatCens = log(sum(CensIntermed)+1);
 arma::mat yvalCens = yval.cols(0,NCens-1);
 arma::mat firstMatCens = log(sum(yvalCens % CensIntermed) + 1);
 
 valToReturn += accu(firstMatCens - secondMatCens);
 
 //Uncensored Piece
 arma::mat secondMatUncens = log(sum(exp(reverse(cumsum(reverse(UncensMat)))))+1);
 arma::mat yvalUncens = yval.cols(NCens,N-1);
 arma::mat firstMatUncens = sum(yvalUncens % UncensMat);
 
 
 valToReturn += accu(firstMatUncens - secondMatUncens);
 valToReturn = (C1/2)*accu(square(thetas)) - valToReturn/N;
 return valToReturn;
 
 }
 
 // [[Rcpp::export]]
 arma::rowvec gradC_Cens_Matrix(NumericVector params, arma::mat yval, arma::mat featureValue, double C1, arma::vec delta) {
double N = featureValue.n_rows;
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
arma::vec BiasGradient = arma::zeros<arma::vec>(biases.n_elem);
arma::mat ThetaGradient = arma::zeros<arma::mat>(thetas.n_rows, thetas.n_cols);


//Censored Piece
double NCens = sum(1-delta);

arma::mat valMat = thetas * featureValue.t();
valMat.each_col() += biases;
arma::mat valMatTrans = exp(reverse(cumsum(reverse(valMat))));

arma::mat NumerMat = cumsum(valMatTrans);
arma::rowvec DenomMat = NumerMat.row(m-1) + 1;
arma::mat BiasSecondMat = NumerMat.each_row()/DenomMat;

arma::cube secondPieceWeight(m,featureValue.n_cols,N);

for(int i=0;i<secondPieceWeight.n_slices;i++){
secondPieceWeight.slice(i) = (NumerMat.col(i) * featureValue.row(i))/DenomMat.at(i); // Creates an m x p x N cube (row, col, slices)
}

//For individual operations
arma::mat CensMat = valMatTrans.cols(0,NCens-1);
arma::mat yvalCens = yval.cols(0,NCens-1);
arma::mat yvalUncens = yval.cols(NCens,N-1);
arma::mat featureValCens = featureValue.rows(0,NCens-1);
arma::mat featureValUncens = featureValue.rows(NCens,N-1);

//Censored Piece
arma::mat NumerCensMat = cumsum(yvalCens % CensMat);
arma::rowvec DenomCensMat = NumerCensMat.row(m-1) +1;
arma::mat BiasCensMat = NumerCensMat.each_row()/DenomCensMat;

arma::cube ThetaCensCube(m,featureValue.n_cols,NCens);
for(int i=0;i<ThetaCensCube.n_slices;i++){
ThetaCensCube.slice(i) = (NumerCensMat.col(i) * featureValCens.row(i))/DenomCensMat.at(i); // Creates an m x p x N cube (row, col, slices)
}

//UncensoredPiece
arma::cube ThetaUncensCube(m,featureValue.n_cols,N-NCens);
for(int i=0;i<ThetaUncensCube.n_slices;i++){
ThetaUncensCube.slice(i) = yvalUncens.col(i) * featureValUncens.row(i); // Creates an m x p x N cube (row, col, slices)
}

arma::cube ThetaCombined = join_slices(ThetaCensCube, ThetaUncensCube);
ThetaCombined -= secondPieceWeight;
ThetaGradient = mean(ThetaCombined, 2);


arma::mat BiasCombined = join_rows(BiasCensMat,yvalUncens);
BiasCombined -= BiasSecondMat;
BiasGradient =  mean(BiasCombined,1);


ThetaGradient = C1*thetas - ThetaGradient;

BiasGradient.transform( [N](double val) { return (val * (double)(-1)); } );
ThetaGradient.transform( [N](double val) { return (val * (double)(1)); } );

arma::vec ThetaFlat = vectorise(ThetaGradient);
return arma::join_cols<arma::mat>(BiasGradient, ThetaFlat).t();
}
*/















