#include <cmath>
#include <tuple>

class EllCalcCore {
  public:
    EllCalcCore(double n_f) {
        _n_f = n_f;
        _half_n = _n_f / 2.0;
        _n_plus_1 = _n_f + 1.0;
        double n_sq = _n_f * _n_f;
        _cst0 = 1.0 / _n_plus_1;
        _cst1 = n_sq / (n_sq - 1.0);
        _cst2 = 2.0 * _cst0;
        _cst3 = _n_f * _cst0;
    }

    std::tuple<double, double, double> calc_parallel_central_cut(double beta1, double tsq) {
        double b1sq = beta1 * beta1;
        double a1sq = b1sq / tsq;
        double temp = _half_n * a1sq;
        double mu_plus_1 = temp + std::sqrt(1.0 - a1sq + temp * temp);
        double mu_plus_2 = mu_plus_1 + 1.0;
        double rho = beta1 / mu_plus_2;
        double sigma = 2.0 / mu_plus_2;
        double temp2 = _n_f * mu_plus_1;
        double delta = temp2 / (temp2 - 1.0);
        return std::make_tuple(rho, sigma, delta);
    }

    std::tuple<double, double, double> calc_parallel_deep_cut(double beta0, double beta1,
                                                              double tsq) {
        double b0b1 = beta0 * beta1;
        double bsum = beta0 + beta1;
        double bsumsq = bsum * bsum;
        double gamma = tsq + _n_f * b0b1;
        double h = tsq + b0b1 + _half_n * bsumsq;
        double temp2 = h + std::sqrt(h * h - gamma * _n_plus_1 * bsumsq);
        double inv_mu_plus_2 = gamma / temp2;
        double inv_mu = gamma / (temp2 - 2.0 * gamma);
        double rho = bsum * inv_mu_plus_2;
        double sigma = 2.0 * inv_mu_plus_2;
        double delta = 1.0 + (-2.0 * b0b1 + bsumsq * inv_mu_plus_2) * inv_mu / tsq;
        return std::make_tuple(rho, sigma, delta);
    }

  private:
    double _n_f;
    double _half_n;
    double _n_plus_1;
    double _cst0;
    double _cst1;
    double _cst2;
    double _cst3;
};
