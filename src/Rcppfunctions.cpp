#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//------------------------------------------------------------------------------
// [[Rcpp::export]]
List Kmatrix_(int N, double Pr1, double Pr2, double Pr3, double Pr4)
{
  double L,c,K_value;
  List outputL(3);
  NumericMatrix K(N,N);
  L=runif(1,Pr1,Pr2)[0];
  c=runif(1,Pr3,Pr4)[0];
  for(int j=0;j<N-1;++j){
    K(j,j)=R_PosInf;
    for(int k=j+1;k<N;++k){
      K_value=pow(c,abs(k-j)-1)*L;
      K(j,k)=K_value;
      K(k,j)=K_value;
    }
    K(N-1,N-1)=R_PosInf;
  }
  outputL[0]=K;
  outputL[1]=L;
  outputL[2]=c;
return outputL;
};

// [[Rcpp::export]]
NumericMatrix KmatrixgivenLc_(int N, double L, double c)
{
  double K_value;
  NumericMatrix K(N,N);
  for(int j=0;j<N-1;++j){
    K(j,j)=R_PosInf;
    for(int k=j+1;k<N;++k){
      K_value=pow(c,abs(k-j)-1)*L;
      K(j,k)=K_value;
      K(k,j)=K_value;
    }
    K(N-1,N-1)=R_PosInf;
  }
  return K;
};

// //------------------------------------------------------------------------------
// // [[Rcpp::export]]
// NumericVector perturb_continuous_(NumericVector theta_c_sampled, NumericMatrix sigma_kernel, NumericMatrix Pr_cont)
// {
//   NumericVector theta_c_perturbed(theta_c_sampled.size());
//   int flag=1;
//   Function f_rmvn("rmvn");
//   while(flag>0){
//   theta_c_perturbed=f_rmvn(1,theta_c_sampled,sigma_kernel);
//   flag=sum((theta_c_perturbed<=Pr_cont(_,0))|(theta_c_perturbed>=Pr_cont(_,1)));
//   }
//   return(theta_c_perturbed);
// }

//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector perturb_continuous_(NumericVector theta_c_sampled, NumericMatrix sigma_kernel)
{
  NumericVector theta_c_perturbed(theta_c_sampled.size());
  Function f_rmvn("rmvn");
    theta_c_perturbed=f_rmvn(1,theta_c_sampled,sigma_kernel);
  return(theta_c_perturbed);
}


// [[Rcpp::export]]
NumericVector perturb_sample_discrete_(int N, NumericVector hat_p_vec,double stay_prob)
{
  NumericVector Rho(N*N);
  NumericVector ber(2),ber2(2),ber3(2),perturb_prob(2);
  ber[0]=1;
  ber[1]=0;
  ber3[0]=0;
  ber3[1]=1;
  perturb_prob[0]=stay_prob;
  perturb_prob[1]=1-stay_prob;
  int perturb=0;
  int counter=0;
  for(int j=0; j<N;++j){
    for(int k=0;k<N;++k ){
      counter=counter+1;
      if(k!=j) {
        ber2[0]=hat_p_vec[counter];
        ber2[1]=1-ber2[0];
        Rho[counter]=sample(ber,1,ber2)[0];
        perturb=sample(ber3,1,perturb_prob)[0];
        Rho[counter]=abs(Rho[counter]-perturb);
        }}}
  return(Rho);}

// nJRNMM_prior_ computes the pdf of the uniform prior (if draw=0) or sample from it (if draw=1)
// Input:  - theta, needed only if draw=0 to compute the corresponding pdf
// Output: - out, pdf of the uniform prior (if draw=0) or sample from it (if draw=1)

// [[Rcpp::export]]
NumericVector nJRNMM_prior_(NumericVector theta,int draw,NumericMatrix Pr_cont){

  int N=Pr_cont.rows()-2;
  NumericVector ber(2);
  ber[0]=0;
  ber[1]=1;
  if(draw == 0){
    NumericVector out(1);
    out=1;
    for(int i=0;i<N+2;++i) {out=out*R::dunif(theta[i],Pr_cont(i,0),Pr_cont(i,1),false);}
    return(out);}
  else{
    NumericVector out(N+2+pow(N,2)); // Here NDIM=
    for(int i=0;i<N+2;++i){
      out[i]=runif(1,Pr_cont(i,0),Pr_cont(i,1))[0];}
    for(int i=N+2;i<N+2+pow(N,2);++i){
      out[i]=sample(ber,1)[0];
    }
    return(out);
  }}


// model_ Function to simulate the observed components from a system of N populations of JRNMMs
// [[Rcpp::export]]
NumericMatrix model_(int N, NumericMatrix sol)
{
  int iter=sol.cols();
  NumericMatrix Y(N,iter);
  int index_l, index_r;
  for(int j=0;j<N;j++){
    index_l=3*j+1;
    index_r=3*j+2;
    Y(j,_)=sol(index_l,_)-sol(index_r,_);
  }
  return(Y);}
