#include <Rcpp.h>
using namespace Rcpp;
#define MPI 3.141592653589793
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

int which_max(int G, NumericVector x){
  int k=0;
  for(int i=1; i<G; i++){
    if(x[i]>x[k]){
      k=i;
    }
  }
  return k;
}


double gaussian_l(int T, NumericVector x, double mu, double sigma){
  double s1=0;
  for(int i=0; i<T; i++){
    s1+=-log(2*MPI)/2-log(sigma)-(x[i]-mu)*(x[i]-mu)/(2*sigma*sigma);
  }
  return s1;
}


NumericMatrix which_A(int g, int m, NumericVector a, IntegerVector hg, NumericMatrix hA){
  int k,j,hgi;
  k=0;
  j=0;
  for(int i=0;i<m;i++){
    hgi=hg[i];
    if(hgi==g){
      k=k+1;
    }
  }
  for(int i=0;i<m;i++){
    hgi=hg[i];
    if(hgi==g){
      a[j]=i;
      j=j+1;
    }
  }
  NumericMatrix A(k,k);
  for(int u=0;u<k;u++){
    for(int v=0;v<k;v++){
      A(u,v)=hA(a[u],a[v]);
    }
  }
  return A;
}

int number(int g, int m, IntegerVector hg ){
  int k,hgi;
  k=0;
  for(int i=0;i<m;i++){
    hgi=hg[i];
    if(hgi==g){
      k=k+1;
    }
  }
  return k;
}

double mean(int k, NumericMatrix A){
  double sum,bar;
  sum=0;
  for(int i=1;i<k;i++){
    for(int j=0;j<i;j++){
      sum=sum+A(i,j);
    }
  }
  bar=sum/(k*(k-1));
  return bar;
}

double which_ci(int n,double rg,NumericMatrix A){
  double c;
  c=0;
  for(int i=1;i<n;i++){
    for(int j=0;j<i;j++){
      for(int u=1;u<n;u++){
        for(int v=1;v<u;v++){
          if(i!=u||j!=v){
            c=c+rg*A(i,j)*A(u,v);
          }
        }
      }
    }
  }
  return c;
}

// [[Rcpp::export]]
double r_g(int n,NumericMatrix A_hat){
  return (pow(std::accumulate(A_hat.begin(), A_hat.end(), 0.0), 2) - 2 * std::inner_product(A_hat.begin(), A_hat.end(), A_hat.begin(), 0.0)) / (pow(n * (n - 1), 2) - 2 * n);
}

// [[Rcpp::export]]
double which_c(int n,double rg,NumericMatrix A){
  double c;
  c=0;
  for(int i=1;i<n;i++){
    for(int j=0;j<i;j++){
      for(int u=1;u<n;u++){
        for(int v=1;v<u;v++){
          if(i!=u||j!=v){
            c=c+rg*A(i,j)*A(u,v);
          }
        }
      }
    }
  }
  return c;
}

// [[Rcpp::export]]
NumericMatrix weight_gic(int m, int T,  NumericMatrix X1, double hrho, NumericVector hmu0,
                        NumericVector hmu1, NumericVector hsigma0, NumericVector hsigma1) {
  NumericMatrix mat(m,m);
  int l=0;
  double a, b, d;
  for(int i=0;i < (m-1); i++){
    for(int j=i+1; j<m; j++){
      a=log(hrho)+gaussian_l( T, X1(l,_),hmu1[i]+hmu1[j] , hsigma1[i]*hsigma1[j]);
      b=log(1-hrho)+ gaussian_l( T, X1(l,_),hmu0[i]+hmu0[j] , hsigma0[i]*hsigma0[j]);
      d=b-a;
      if(d>30){
        mat(i,j)=0;
      }
      if(d<-30){
        mat(i,j)=1;
      }
      if(d>-30 & d<30){
        mat(i,j)=1/(1+exp(d));
      }
      mat(j,i)=mat(i,j);
      l+=1;
    }
  }
  return mat;
}



// [[Rcpp::export]]
NumericMatrix weight_gc(int m, int T, double c, IntegerVector hg, NumericMatrix X1, NumericMatrix hrho, NumericMatrix hmu0,
                        NumericMatrix hmu1, NumericMatrix hsigma0, NumericMatrix hsigma1) {
  NumericMatrix mat(m,m);
  int g1, g2, l=0;
  double a, b, d;
  for(int i=0;i < (m-1); i++){
    for(int j=i+1; j<m; j++){
      g1=hg(i);
      g2=hg(j);
      a=log(hrho(g1,g2))+gaussian_l( T, X1(l,_),hmu1(g1,g2) , hsigma1(g1,g2));
      b=log(1-hrho(g1,g2))+ gaussian_l( T, X1(l,_),hmu0(g1,g2) , hsigma0(g1,g2));
      d=b-a;
      // if(d>30){
      //   mat(i,j)=0;
      //   }
      // if(d<-30){
      //   mat(i,j)=1;
      //   }
      mat(i,j)=1/(1+exp(d)/c);
      mat(j,i)=mat(i,j);
      l+=1;
    }
    
  }
  return mat;
}


// [[Rcpp::export]]
IntegerVector hg_gc(int m, int G, NumericVector a, IntegerVector hg, NumericVector r, NumericMatrix hA, NumericMatrix N, NumericMatrix Y, NumericMatrix Y2, NumericMatrix weight,
                   NumericMatrix hrho, NumericMatrix hmu0, NumericMatrix hmu1, NumericMatrix hsigma0, NumericMatrix hsigma1) {
  int gj,n;
  double ci,rg,A_bar;
  NumericMatrix  val(m,G);
  for(int i=0;i<m;i++){
    for(int g=0;g<G;g++){
      hg[i]=g;
      ci=0;
      for(int g1=0;g1<G;g1++){
        n=number(g1,m,hg);
        if(n>1){
          NumericMatrix A(n,n);
          A=which_A(g1,m,a,hg,hA);
          A_bar=mean(n,A);
          NumericMatrix A_hat(n,n);
          A_hat=A-A_bar;
          rg=r[g1];
          ci=ci+(pow(std::accumulate(A_hat.begin(), A_hat.end(), 0.0), 2) - 2 * std::inner_product(A_hat.begin(), A_hat.end(), A_hat.begin(), 0.0)) / (pow(n * (n - 1), 2) - 2 * n);
        }
      }
      ci=ci+1;
      for(int j=0;j<m;j++){
        gj=hg[j];
        if(j!=i){
          val(i,g)=val(i,g)+(1-weight(i,j))*(-N(i,j)*log(hsigma0(g,gj)) - (Y2(i,j)-2*hmu0(g,gj)*Y(i,j)+N(i,j)*hmu0(g,gj)*hmu0(g,gj))/(2*hsigma0(g,gj)*hsigma0(g,gj)) +log(1-hrho(g,gj)))+weight(i,j)*(-N(i,j)*log(hsigma1(g,gj)) - (Y2(i,j)-2*hmu1(g,gj)*Y(i,j)+N(i,j)*hmu1(g,gj)*hmu1(g,gj))/(2*hsigma1(g,gj)*hsigma1(g,gj)) + log(hrho(g,gj)))+log(ci);
        }
      }
    }
    hg[i]=which_max(G,val(i,_));
  }
  return hg;
}


// [[Rcpp::export]]
double cal_B(int k, NumericMatrix ij_comb, NumericMatrix uv_comb,NumericVector rho_store, NumericMatrix Y_hat, NumericMatrix alpha){
  double result = 0.0;
  for (int ij = 0; ij < ij_comb.ncol(); ij++) {
    for(int uv = 0; uv < uv_comb.ncol(); uv++){
      result = result + alpha(ij_comb(0, ij)-1, k-1) * alpha(ij_comb(1, ij)-1, k-1) * alpha(ij_comb(0, uv)-1, k-1)*alpha(ij_comb(1, uv)-1, k-1)*rho_store[k-1]*Y_hat(ij_comb(0, ij)-1, ij_comb(1, ij)-1) * Y_hat(uv_comb(0, uv)-1, uv_comb(1, uv)-1);
    }
  }
  return result;
}

// [[Rcpp::export]]
double cal_A(int K, NumericMatrix Y, NumericMatrix alpha, NumericMatrix xy_comb, NumericMatrix mu){
  double log_likelihood = 0.0;
  for(int Q = 0; Q < K; Q++){
    for(int l = 0; l < K; l++){ 
      for(int xy = 0; xy < xy_comb.ncol(); xy++){
        log_likelihood = log_likelihood + alpha(xy_comb(0, xy)-1, Q) * alpha(xy_comb(1, xy)-1, l) * (Y(xy_comb(0, xy)-1, xy_comb(1, xy)-1) * log(mu(Q, l)) + (1 - Y(xy_comb(0, xy)-1, xy_comb(1, xy)-1)) * log(1 - mu(Q, l)));
      }
    }
  }
  return log_likelihood;
}


// [[Rcpp::export]]
NumericVector cal_C(NumericMatrix rho,NumericVector rho_store, NumericVector count, NumericVector z, int N,int K){
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){ 
      for(int u = 0; u < N; u++){
        for(int v = 0; v < N; v++){
          if(z[i] == z[j] && z[u] == z[v] && z[u] == z[i]){
            rho_store[z[i]-1] = rho_store[z[i]-1] + rho[j+N*i,v+N*u];
            count[z[i]-1]= count[z[i]-1] + 1;
          }
        }
      }
    }
  }
  for(int i = 0; i < K; i++){
    rho_store[i] = rho_store[i]/count[i];
  }
  return rho_store;
}



// [[Rcpp::export]]
NumericVector compute_Y_hat(NumericVector Y) {
  IntegerVector dims = Y.attr("dim");
  int N = dims[0];
  int M = dims[2];
  NumericVector Y_hat = NumericVector(Dimension(N, N, M));
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      NumericVector slice(M);
      for (int k = 0; k < M; ++k) {
        slice[k] = Y[i + j * N + k * N * N];
      }
      double mean_val = mean(slice);
      
      if (mean_val == 1.0 || mean_val == 0.0) {
        for (int k = 0; k < M; ++k) {
          Y_hat[i + j * N + k * N * N] = 0.0;
        }
      } else {
        double denominator = sqrt(mean_val * (1.0 - mean_val));
        for (int k = 0; k < M; ++k) {
          Y_hat[i + j * N + k * N * N] = (slice[k] - mean_val) / denominator;
        }
      }
    }
  }
  
  return Y_hat;
}



// [[Rcpp::export]]
NumericMatrix compute_mu(const NumericMatrix& alpha, const NumericVector& Y, int N, int M, int K) {
  NumericMatrix mu(K, K);
  for (int q = 0; q < K; ++q) {
    for (int l = 0; l < K; ++l) {
      double x1 = 0.0;
      double x2 = 0.0;
      for (int m = 0; m < M; ++m) {
        for (int i = 0; i < N; ++i) {
          for (int j = 0; j < N; ++j) {
            if (i != j) {
              int idx = i + j * N + m * N * N; 
              x1 += alpha(i, q) * alpha(j, l) * Y[idx];
              x2 += alpha(i, q) * alpha(j, l);
            }
          }
        }
      }
      mu(q, l) = x1 / x2;
    }
  }
  
  return mu;
}



// [[Rcpp::export]]
NumericMatrix update_alpha(const NumericMatrix& alpha, const NumericMatrix& result_L, int N, int K) {
  NumericMatrix alpha_new(N, K);
  
  for (int i = 0; i < N; ++i) {
    for (int q = 0; q < K; ++q) {
      double sumk = 0.0;
      for (int sum_k = 0; sum_k < K; ++sum_k) {
        double ratio = alpha(i, sum_k) / alpha(i, q);
        double exp_diff = exp(result_L(i, sum_k) - result_L(i, q));
        sumk += ratio * exp_diff;
      }
      if (std::isnan(sumk)) {
        alpha_new(i, q) = 0.0;
      } else {
        alpha_new(i, q) = 1.0 / sumk;
        if (std::isnan(alpha_new(i, q))) {
          alpha_new(i, q) = 0.0;
        }
      }
    }
  }
  
  return alpha_new;
}



// [[Rcpp::export]]
NumericMatrix compute_result_L(int N, int K, int M, NumericMatrix alpha_save, NumericVector rho_store, NumericVector Y_hat,
                               NumericMatrix z, NumericMatrix xy_comb, NumericMatrix mu, NumericVector Y) {
  NumericMatrix result_L(N, K);
  
  for (int i = 0; i < N; ++i) {
    for (int q = 0; q < K; ++q) {
      NumericMatrix alpha = clone(alpha_save);
      for (int k = 0; k < K; ++k) {
        alpha(i, k) = (k == q) ? 1 : 0;
      }
      
      double log_likelihood = 0.0;
      
      // First part: cal_A contributions
      for (int m = 0; m < M; ++m) {
        NumericMatrix Y_matrix(N, N);
        for (int row = 0; row < N; row++) {
          for (int col = 0; col < N; col++) {
            Y_matrix(row, col) = Y[row + col * N + m * N * N];
          }
        }
        double log_likelihood_a = 0.0;
        for(int Q = 0; Q < K; Q++){
          for(int l = 0; l < K; l++){ 
            for(int xy = 0; xy < xy_comb.ncol(); xy++){
              log_likelihood_a = log_likelihood_a + alpha(xy_comb(0, xy)-1, Q) * alpha(xy_comb(1, xy)-1, l) * (Y_matrix(xy_comb(0, xy)-1, xy_comb(1, xy)-1) * log(mu(Q, l)) + (1 - Y_matrix(xy_comb(0, xy)-1, xy_comb(1, xy)-1)) * log(1 - mu(Q, l)));
            }
          }
        }
        log_likelihood += log_likelihood_a;
      }
      
      // Second part: sum_result contributions
      for (int m = 0; m < M; ++m) {
        double sum_result = 0.0;
        
        for (int k = 1; k <= K; ++k) {
          std::vector<int> which_z_k;
          for (int idx = 0; idx < N; ++idx) {
            if (z[idx] == k) {
              which_z_k.push_back(idx);
            }
          }
          
          int len = which_z_k.size();
          if (len > 0) {
            NumericVector sub_alpha(len);
            NumericMatrix sub_Y_hat(len, len);
            
            for (int x = 0; x < len; ++x) {
              sub_alpha[x] = alpha(which_z_k[x], k - 1);
            }
            
            for (int x = 0; x < len; ++x) {
              for (int y = 0; y < len; ++y) {
                sub_Y_hat(x, y) = Y_hat[which_z_k[x] + which_z_k[y] * N + m * N * N];
              }
            }
            
            double inner_sum = 0.0;
            for (int i = 0; i < len; ++i) {
              for (int j = 0; j < len; ++j) {
                double outer_val = sub_alpha[i] * sub_alpha[j];
                inner_sum += outer_val * outer_val * sub_Y_hat(i, j) * sub_Y_hat(i, j);
              }
            }
            
            double result = rho_store[k - 1] * inner_sum / 2.0;
            sum_result += std::max(result, 0.0);
          }
        }
        
        log_likelihood += log(1.0 + sum_result);
      }
      result_L(i, q) = log_likelihood;
    }
  }
  
  return result_L;
}
