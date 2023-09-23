#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <iostream>
#include <random>
#include <vector>
#include <unordered_map>
#include <cmath>

// helper progress bar function
void updateProgressBar(double progress) {
  const int barWidth = 50;
  Rcpp::Rcout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) Rcpp::Rcout << "=";
    else if (i == pos) Rcpp::Rcout << ">";
    else Rcpp::Rcout << " ";
  }
  Rcpp::Rcout << "] " << std::setw(3) << int(progress * 100.0) << "%\r";
  Rcpp::Rcout.flush();
}

// gsdmm gibbs sampler
// [[Rcpp::export]]
Rcpp::List gsdmm_gibbs(
    std::vector<std::vector<int>> d,
    int I,
    int K,
    double alpha,
    double beta,
    int V,
    bool progress
) {
  int D = d.size();
  arma::vec Z(D, arma::fill::zeros);
  arma::vec m_z(K, arma::fill::zeros); // number of documents per cluster
  arma::vec n_z(K, arma::fill::zeros); // number of words per cluster
  std::vector<std::vector<int>> n_z_w(V, std::vector<int>(K, 0)); // number of occurrences of word w per cluster

  // generate initial values ---------------------------------------------------
  // initialize equal probabilities for each cluster
  std::vector<double> p(K, (1.0 / static_cast<double>(K)));

  // initialize multinomial
  unsigned int seed = floor((unif_rand() * 100000));
  std::mt19937 gen(seed);
  std::discrete_distribution<> distribution(p.begin(), p.end());

  std::vector<std::vector<int>> d_r(D);
  std::vector<std::unordered_map<int, int>> n_d_w(D); // word counts per document
  RcppThread::parallelFor(0, D, [&](int doc) {
    for(std::size_t w = 0; w < d[doc].size(); ++w) {
      n_d_w[doc][d[doc][w]]++;
    }
    // keep only unique words in each document
    std::vector<int> w_d = d[doc];
    std::sort(w_d.begin(), w_d.end());
    w_d.erase(std::unique(w_d.begin(), w_d.end()), w_d.end()); // keep unique words
    d_r[doc] = w_d;
  });

  // initialize values
  for(int doc = 0; doc < D; ++doc) {
    int z = distribution(gen); // sample cluster
    Z(doc) = z; // initial cluster
    m_z(z)++; // add doc to cluster
    n_z(z) += d[doc].size(); // add number of words to cluster

    // add occurences of word w to cluster
    for(std::size_t w = 0; w < d_r[doc].size(); ++w) {
      int word = d_r[doc][w];
      n_z_w[word][z] += n_d_w[doc][word];
    }
  }

  // pre compute constants
  // compute fraction 1 denominator
  double denom_1 = std::log((D - 1 + K * alpha));
  // compute fraction 1 numerator
  arma::vec num_1_vec = arma::log((m_z + alpha));

  // run gibbs sampler ---------------------------------------------------------
  int cluster_count = K;
  for(int i = 0; i < I; ++i) { // for I iterations
    if (i % 20 == 0) RcppThread::checkUserInterrupt();

    for(int doc = 0; doc < D; ++doc) { // for each document
      int z = Z(doc); // get current cluster of document
      m_z(z)--; // remove document from cluster
      n_z(z) -= d[doc].size(); // remove total document word count from cluster

      for(std::size_t w = 0; w < d_r[doc].size(); ++w) { // remove individual document word counts from cluster
        int word = d_r[doc][w];
        n_z_w[word][z] -= n_d_w[doc][word];
      }

      RcppThread::parallelFor(0, K, [&](int k) {
        // get fraction 1 numerator
        double num_1 = num_1_vec(k);
        // compute fraction 2 numerator
        double num_2 = 0.0;
        for(std::size_t w = 0; w < d_r[doc].size(); ++w) {
          int word = d_r[doc][w];
          for(int j = 0; j < n_d_w[doc][word]; ++j) {
            num_2 += std::log((n_z_w[word][k] + beta + j));
          }
        }
        // compute fraction 2 denominator
        double denom_2 = arma::sum(arma::log(n_z(k) + V * beta + arma::linspace(1, d[doc].size(), d[doc].size())));

        // compute p from eq 4
        p[k] = std::exp((num_1 - denom_1 + num_2 - denom_2));
      });

      // normalize p
      double psum = 0.0;
      for(int j = 0; j < K; ++j) {
        psum += p[j];
      }
      for(int j = 0; j < K; ++j) {
        p[j] = p[j] / psum;
      }

      // sample new cluster
      std::discrete_distribution<> distribution_d(p.begin(), p.end());
      z = distribution_d(gen);
      Z(doc) = z;

      m_z(z)++; // add doc to cluster
      n_z(z) += d[doc].size(); // add number of words to cluster
      // add occurences of word w to cluster
      for(std::size_t w = 0; w < d[doc].size(); ++w) {
        int word = d[doc][w];
        n_z_w[word][z]++;
      }
    }

    if (progress) {
      double progress = static_cast<double>(i + 1) / I;
      updateProgressBar(progress);
    }
    if (i > 25) {
      int cluster_count_new = 0;
      for(const auto& v: m_z) {
        if (v != 0) cluster_count_new++;
      }
      if (cluster_count_new == cluster_count) {
        if (progress) updateProgressBar(1.0);
        break;
      }
      cluster_count = cluster_count_new;
    }
  }
  if (progress) {
    Rcpp::Rcout << std::endl;
  }

  std::vector<int> d_doc;
  std::vector<int> d_tok;
  std::vector<int> d_val;
  for(std::size_t i = 0; i < n_z_w.size(); i++) {
    for(std::size_t j = 0; j < n_z_w[i].size(); j++) {
      if (n_z_w[i][j] != 0) {
        d_doc.push_back(i);
        d_tok.push_back(j);
        d_val.push_back(n_z_w[i][j]);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("cluster") = Rcpp::wrap(arma::conv_to<std::vector<int>>::from(Z)),
    Rcpp::Named("distribution") = Rcpp::DataFrame::create(
      Rcpp::Named("doc") = Rcpp::wrap(d_doc),
      Rcpp::Named("tok") = Rcpp::wrap(d_tok),
      Rcpp::Named("value") = Rcpp::wrap(d_val)
    )
  );
}
