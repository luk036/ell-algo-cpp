#include <cassert>
#include <cmath>                      // for sqrt
#include <ellalgo/ell_calc_core.hpp>  // for EllCalcCore
#include <tuple>                      // for tuple

/**
 * The function calculates and returns three values (rho, sigma, and delta) based on the input
 * parameters (beta0, beta1, and tsq).
 * 
 * @param[in] beta0 The parameter `beta0` represents the value of beta for the first variable.
 * @param[in] beta1 The parameter `beta1` represents a value used in the calculation. 
 * @param[in] tsq tsq is a constant value of type double. It represents the square of the parameter t.
 * 
 * @return The function `calc_parallel_cut` returns a tuple containing three values: `rho`, `sigma`,
 * and `delta`.
 */
auto EllCalcCore::calc_parallel_cut(const double& beta0, const double& beta1,
                                    const double& tsq) const -> std::tuple<double, double, double> {
    auto b0b1 = beta0 * beta1;
    auto b0b1n = b0b1 / tsq;
    auto t1n = 1.0 - beta1 * (beta1 / tsq);
    auto t0n = 1.0 - beta0 * (beta0 / tsq);
    auto bsum = beta0 + beta1;
    auto bsumn = bsum / tsq;
    auto bav = bsum / 2.0;
    auto tempn = this->_half_n * bsumn * (beta1 - beta0);
    auto xi = std::sqrt(t0n * t1n + tempn * tempn);
    auto sigma = this->_cst3 + (1.0 + b0b1n - xi) / (bsumn * bav) / this->_n_plus_1;
    auto&& rho = sigma * bav;
    auto&& delta = this->_cst1 * ((t0n + t1n) / 2.0 + xi / this->_n_f);
    return {rho, sigma, delta};
}

/**
 * @brief Parallel Central Cut
 *
 * The function `calc_parallel_central_cut` calculates and returns the values of rho, sigma, and delta
 * based on the given input parameters under the parallel central cuts:
 *
 *        g' (x - xc) <= 0,
 *        g' (x - xc) + beta1 >= 0.
 * 
 * @param[in] beta1 The parameter `beta1` represents a double value. 
 * @param[in] tsq tsq is a constant value representing the square of the variable tau.
 * 
 * @return The function `calc_parallel_central_cut` returns a tuple containing three values: rho,
 * sigma, and delta.
 */
auto EllCalcCore::calc_parallel_central_cut(const double& beta1, const double& tsq) const
    -> std::tuple<double, double, double> {
    auto b1sqn = beta1 * beta1 / tsq;
    auto temp = this->_half_n * b1sqn;
    auto xi = std::sqrt(1.0 - b1sqn + temp * temp);
    auto&& delta = this->_cst1 * (1.0 - b1sqn / 2.0 + xi / this->_n_f);
    auto sigma = this->_cst3 + this->_cst2 * (1.0 - xi) / b1sqn;
    auto&& rho = sigma * beta1 / 2;
    return {rho, sigma, delta};
}

/**
 * @brief Calculate new ellipsoid under Non-central Cut
 *
 * The function `calc_bias_cut` calculates and returns the values of rho, sigma, and delta based on the
 * given beta and tau values under the bias-cut:
 *
 *        g' (x - xc) + beta \le 0
 *
 * @param[in] beta The parameter "beta" represents a value used in the calculation. It is a double value.
 * @param[in] tau The parameter "tau" represents a value used in the calculation. It is not specified in
 * the code snippet provided. You would need to provide a value for "tau" in order to use the
 * `calc_bias_cut` function.
 * 
 * @return The function `calc_bias_cut` returns a tuple containing the following values:
 */
auto EllCalcCore::calc_bias_cut(const double& beta, const double& tau) const
    -> std::tuple<double, double, double> {
    auto alpha = beta / tau;
    auto gamma = tau + this->_n_f * beta;
    auto&& sigma = this->_cst2 * gamma / (tau + beta);
    auto&& rho = gamma / this->_n_plus_1;
    auto&& delta = this->_cst1 * (1.0 - alpha * alpha);
    return {rho, sigma, delta};
}

/**
 * @brief Central Cut
 *
 * The function `_calc_deep_cut_core` calculates and returns the values of rho, sigma, and delta based
 * on the given beta, tau, and gamma values under the central-cut:
 *
 *        g' (x - xc) \le 0.
 * 
 * @param[in] tau tau is a constant value of type double. It represents the square of the variable tau.
 * 
 * @return A tuple containing the values of rho, sigma, and delta.
 */
auto EllCalcCore::calc_central_cut(const double& tau) const -> std::tuple<double, double, double> {
    auto&& sigma = this->_cst2;
    auto&& rho = tau / this->_n_plus_1;
    auto&& delta = this->_cst1;
    return {rho, sigma, delta};
}
