#pragma once

#include <cmath>
#include <tuple>

/**
 * @brief Ellipsoid Search Space
 *
 *  EllCalcCore = {x | (x - xc)' mq^-1 (x - xc) \le \kappa}
 *
 * Keep $Q$ symmetric but no promise of positive definite
 */
class EllCalcCore {
  private:
    const double _n_f;
    const double _n_plus_1;
    const double _half_n;
    const double _n_sq;
    const double _cst1;
    const double _cst2;
    const double _cst3;

    /**
     * @brief Construct a new EllCalcCore object
     *
     * @param[in] E
     */
    // auto operator=(const EllCalcCore& E) const -> EllCalcCore& = delete;

  public:
    /**
     * @brief Construct a new EllCalcCore object
     *
     * @tparam V
     * @tparam U
     * @param kappa
     * @param mq
     * @param x
     */
    EllCalcCore(size_t ndim)
        : _n_f{double(ndim)},
          _n_plus_1{_n_f + 1.0},
          _half_n{_n_f / 2.0},
          _n_sq{_n_f * _n_f},
          _cst1{_n_sq / (_n_sq - 1.0)},
          _cst2{2.0 / _n_plus_1},
          _cst3{_n_f / _n_plus_1} {}

  public:
    /**
     * @brief Construct a new EllCalcCore object
     *
     * @param[in] E (move)
     */
    EllCalcCore(EllCalcCore &&E) = default;

    /**
     * @brief Destroy the EllCalcCore object
     *
     */
    ~EllCalcCore() {}

    /**
     * @brief Construct a new EllCalcCore object
     *
     * To avoid accidentally copying, only explicit copy is allowed
     *
     * @param E
     */
    EllCalcCore(const EllCalcCore &E) = default;

    /**
     * @brief Calculate a new ellipsoid under a parallel cut
     *
     *        g' (x - xc) + beta0 \le 0
     *        g' (x - xc) + beta1 \ge 0
     *
     * @param[in] beta0
     * @param[in] beta1
     * @param[in] tsq
     * @return std::tuple<double, double, double>
     */
    auto calc_parallel_cut(const double &beta0, const double &beta1, const double &tsq) const
        -> std::tuple<double, double, double>;

    /**
     * @brief Calculate new ellipsoid under Parallel Cut, one of them is central
     *
     *        g' (x - xc) \le 0
     *        g' (x - xc) + beta1 \ge 0
     *
     * @param[in] beta1
     * @param[in] b1sq
     * @param[in] tsq
     * @return std::tuple<double, double, double>
     */
    auto calc_parallel_central_cut(const double &beta1, const double &tsq) const
        -> std::tuple<double, double, double>;

    /**
     * @brief Calculate new ellipsoid under Non-central Cut
     *
     *        g' (x - xc) + beta \le 0
     *
     * @param[in] beta
     * @param[in] tsq
     * @return std::tuple<double, double, double>
     */
    auto calc_bias_cut(const double &beta, const double &tsq) const
        -> std::tuple<double, double, double>;

    /**
     * @brief Calculate new ellipsoid under Central Cut
     *
     *        g' (x - xc) \le 0
     *
     * @param[in] tau
     * @return std::tuple<double, double, double>
     */
    auto calc_central_cut(const double &tau) const -> std::tuple<double, double, double>;
};  // } EllCalcCore
