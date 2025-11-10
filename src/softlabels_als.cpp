#include <RcppArmadillo.h>

using arma::mat;
using arma::sp_mat;
using arma::vec;
using arma::rowvec;
using arma::speye;
using arma::norm;

static inline double sq(double x) { return x * x; }

static inline void project_nonneg(vec &v) {
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    if (v[i] < 0.0) v[i] = 0.0;
  }
}

static inline void project_nonneg_inplace(mat &M) {
  for (arma::uword j = 0; j < M.n_cols; ++j) {
    for (arma::uword i = 0; i < M.n_rows; ++i) {
      if (M(i, j) < 0.0) M(i, j) = 0.0;
    }
  }
}

static inline void compute_objective_terms(
    const arma::mat &X,
    const arma::mat &W,
    const arma::rowvec &b,
    const arma::mat &P,
    const arma::mat &DBbeta,
    double lambda_W,
    double lambda_HRF,
    double lambda_smooth,
    const arma::sp_mat &L,
    double &recon, double &regW, double &prior, double &rough) {
  mat S = X * W;
  S.each_row() += b;
  recon = arma::accu(arma::square(S - P));
  regW  = lambda_W * arma::accu(arma::square(W));
  prior = lambda_HRF * arma::accu(arma::square(P - DBbeta));
  rough = lambda_smooth * arma::accu((L * P) % P);
}

// [[Rcpp::export(name = "engine_objective")]]
Rcpp::List engine_objective_cpp(const arma::mat &X,
                                const arma::mat &W,
                                const arma::rowvec &b,
                                const arma::mat &P,
                                const arma::mat &DBbeta,
                                double lambda_W,
                                double lambda_HRF,
                                double lambda_smooth,
                                const arma::sp_mat &L) {
  double recon = 0.0, regW = 0.0, prior = 0.0, rough = 0.0;
  compute_objective_terms(X, W, b, P, DBbeta, lambda_W, lambda_HRF, lambda_smooth, L,
                          recon, regW, prior, rough);
  return Rcpp::List::create(
    Rcpp::Named("value") = recon + regW + prior + rough,
    Rcpp::Named("recon") = recon,
    Rcpp::Named("regW")  = regW,
    Rcpp::Named("prior") = prior,
    Rcpp::Named("rough") = rough
  );
}

// [[Rcpp::export]]
Rcpp::List fit_softlabels_als(const arma::mat &X,
                              const arma::mat &P0,
                              const arma::sp_mat &L,
                              const arma::mat &DBbeta,
                              double lambda_W,
                              double lambda_HRF,
                              double lambda_smooth,
                              int max_iter,
                              double tol,
                              bool nonneg,
                              bool threads) {
  (void)threads;
  const int T = X.n_rows;
  const int V = X.n_cols;
  const int K = P0.n_cols;
  if (P0.n_rows != T || DBbeta.n_rows != T || DBbeta.n_cols != K) {
    Rcpp::stop("Dimension mismatch between inputs.");
  }
  mat P = P0;
  rowvec Xbar = arma::mean(X, 0);
  mat Xc = X;
  Xc.each_row() -= Xbar;
  mat XtXc = Xc.t() * Xc;
  mat Aw = XtXc + lambda_W * arma::eye<mat>(V, V);
  sp_mat lap_op = (1.0 + lambda_HRF) * speye<sp_mat>(T, T) + lambda_smooth * L;
  mat W(V, K, arma::fill::zeros);
  rowvec b(K, arma::fill::zeros);
  mat XcW(T, K, arma::fill::zeros);
  std::vector<double> obj_val, obj_recon, obj_regW, obj_prior, obj_rough;
  obj_val.reserve(std::max(1, max_iter));
  obj_recon.reserve(std::max(1, max_iter));
  obj_regW.reserve(std::max(1, max_iter));
  obj_prior.reserve(std::max(1, max_iter));
  obj_rough.reserve(std::max(1, max_iter));
  std::vector<double> w_norm, p_norm, dW_vec, dP_vec, rel_dW_vec, rel_dP_vec;
  w_norm.reserve(std::max(1, max_iter));
  p_norm.reserve(std::max(1, max_iter));
  dW_vec.reserve(std::max(1, max_iter));
  dP_vec.reserve(std::max(1, max_iter));
  rel_dW_vec.reserve(std::max(1, max_iter));
  rel_dP_vec.reserve(std::max(1, max_iter));
  mat W_prev = W;

  double prev_norm = arma::norm(P, "fro");
  double rel_change = tol + 1.0;
  int iter = 0;

  for (iter = 0; iter < max_iter; ++iter) {
    rowvec Pbar = arma::mean(P, 0);
    mat Pc = P;
    Pc.each_row() -= Pbar;
    mat rhs = Xc.t() * Pc;
    W = arma::solve(Aw, rhs, arma::solve_opts::likely_sympd);
    b = Pbar - Xbar * W;
    double dW = arma::norm(W - W_prev, "fro");
    double rel_dW = dW / std::max(1e-12, arma::norm(W_prev, "fro"));
    W_prev = W;
    XcW = Xc * W;

    mat P_prev = P;
    mat RHS = XcW;
    RHS.each_row() += (b + Xbar * W);
    RHS += lambda_HRF * DBbeta;
    mat P_new;
    bool ok = false;
    try { ok = arma::spsolve(P_new, lap_op, RHS, "superlu"); }
    catch (std::logic_error &) { ok = false; }
    if (!ok) {
      try { ok = arma::spsolve(P_new, lap_op, RHS, "lapack"); }
      catch (std::logic_error &) { ok = false; }
    }
    if (!ok) {
      Rcpp::stop("Failed to solve smoothing system (spsolve).");
    }
    if (nonneg) {
      project_nonneg_inplace(P_new);
    }
    P = std::move(P_new);
    double diff = arma::norm(P - P_prev, "fro");
    double denom = std::max(1e-12, prev_norm);
    rel_change = diff / denom;
    prev_norm = arma::norm(P, "fro");
    double dP = diff;
    double rel_dP = rel_change;
    double Wnorm = arma::norm(W, "fro");
    double Pnorm = prev_norm;
    double recon=0.0, regW=0.0, prior=0.0, rough=0.0;
    compute_objective_terms(X, W, b, P, DBbeta, lambda_W, lambda_HRF, lambda_smooth, L,
                            recon, regW, prior, rough);
    obj_val.push_back(recon + regW + prior + rough);
    obj_recon.push_back(recon);
    obj_regW.push_back(regW);
    obj_prior.push_back(prior);
    obj_rough.push_back(rough);
    w_norm.push_back(Wnorm);
    p_norm.push_back(Pnorm);
    dW_vec.push_back(dW);
    dP_vec.push_back(dP);
    rel_dW_vec.push_back(rel_dW);
    rel_dP_vec.push_back(rel_dP);
    if (rel_change < tol) {
      break;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("W") = W,
    Rcpp::Named("P") = P,
    Rcpp::Named("b") = b,
    Rcpp::Named("iterations") = iter + 1,
    Rcpp::Named("converged") = rel_change < tol,
    Rcpp::Named("obj_trace") = Rcpp::List::create(
      Rcpp::Named("value") = obj_val,
      Rcpp::Named("recon") = obj_recon,
      Rcpp::Named("regW")  = obj_regW,
      Rcpp::Named("prior") = obj_prior,
      Rcpp::Named("rough") = obj_rough,
      Rcpp::Named("w_norm") = w_norm,
      Rcpp::Named("p_norm") = p_norm,
      Rcpp::Named("dW") = dW_vec,
      Rcpp::Named("dP") = dP_vec,
      Rcpp::Named("rel_dW") = rel_dW_vec,
      Rcpp::Named("rel_dP") = rel_dP_vec
    )
  );
}

// [[Rcpp::export]]
arma::mat predict_softlabels(const arma::mat &Xtest, const arma::mat &W, const arma::rowvec &b) {
  if (Xtest.n_cols != W.n_rows) {
    Rcpp::stop("Dimension mismatch between Xtest and W.");
  }
  mat S = Xtest * W;
  if (b.n_elem == S.n_cols) {
    S.each_row() += b;
  }
  return S;
}
