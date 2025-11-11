/*************************************************************************************************************
 *  Class for Discrete Extended Kalman Filter
 *  The system to be estimated is defined as a discrete nonlinear dynamic dystem:
 *              x(k) = f[x(k-1), u(k-1)] + v(k)     ; x = Nx1,    u = Mx1
 *              y(k) = h[x(k)] + n(k)               ; y = Zx1
 *
 *        Where:
 *          x(k) : State Variable at time-k                          : Nx1
 *          y(k) : Measured output at time-k                         : Zx1
 *          u(k) : System input at time-k                            : Mx1
 *          v(k) : Process noise, AWGN assumed, w/ covariance Qn     : Nx1
 *          n(k) : Measurement noise, AWGN assumed, w/ covariance Rn : Nx1
 *
 *          f(..), h(..) is a nonlinear transformation of the system to be estimated.
 *
 ***************************************************************************************************
 *      Extended Kalman Filter algorithm:
 *          Initialization:
 *              x(k=0|k=0) = Expected value of x at time-0 (i.e. x(k=0)), typically set to zero.
 *              P(k=0|k=0) = Identity matrix * covariant(P(k=0)), typically initialized with some
 *                            big number.
 *              Q, R       = Covariance matrices of process & measurement. As this implementation
 *                            the noise as AWGN (and same value for every variable), this is set
 *                            to Q=diag(QInit,...,QInit) and R=diag(RInit,...,RInit).
 *
 *
 *          EKF Calculation (every sampling time):
 *              Calculate the Jacobian matrix of f (i.e. F):
 *                  F        = d(f(..))/dx |x(k-1|k-1),u(k-1)                        ...{EKF_1}
 *
 *              Predict x(k) through nonlinear function f:
 *                  x(k|k-1) = f[x(k-1|k-1), u(k-1)]                                 ...{EKF_2}
 *
 *              Predict P(k) using linearized f (i.e. F):
 *                  P(k|k-1) = F*P(k-1|k-1)*F' + Q                                   ...{EKF_3}
 *
 *              Calculate the Jacobian matrix of h (i.e. H):
 *                  H        = d(h(..))/dx |x(k|k-1)                                 ...{EKF_4}
 *
 *              Predict residual covariance using linearized h (i.e. S):
 *                  S        = H*P(k|k-1)*H' + R                                     ...{EKF_5}
 *
 *              Calculate the kalman gain:
 *                  K        = P(k|k-1)*H'*(S^-1)                                    ...{EKF_6}
 *
 *              Correct x(k) using kalman gain:
 *                  x(k|k)   = x(k|k-1) + K*[y(k) - h(x(k|k-1))]                     ...{EKF_7}
 *
 *              Correct P(k) using kalman gain:
 *                  P(k|k)   = (I - K*H)*P(k|k-1)                                    ...{EKF_8}
 *
 *
 *        *Additional Information:
 *              - Pada contoh di atas X~(k=0|k=0) = [0]. Untuk mempercepat konvergensi bisa
 *                  digunakan informasi plant-spesific. Misal pada implementasi Kalman Filter
 *                  untuk sensor IMU (Inertial measurement unit) dengan X = [quaternion], dengan
 *                  asumsi IMU awalnya menghadap ke atas tanpa rotasi: X~(k=0|k=0) = [1, 0, 0, 0]'
 *
 *
 * Based on implementation from: https://github.com/pronenewbits
 ************************************************************************************************************/
#pragma once
//
#include <concepts>
//
#include "GenericMath/GenericMatrix.h"

namespace GenericMath
{

	template<class ConcreateFilter, MatrixIdx U_LEN, MatrixIdx X_LEN, MatrixIdx Y_LEN>
	class ExtendedKalmanFilter;

	template<class ConcreateFilter, MatrixIdx U_LEN, MatrixIdx X_LEN, MatrixIdx Y_LEN>
	concept KalmanFilter = std::derived_from<ConcreateFilter, ExtendedKalmanFilter<ConcreateFilter, U_LEN, X_LEN, Y_LEN>> &&
						   requires(ConcreateFilter filter, Matrix<double, X_LEN, X_LEN> F, Matrix<double, Y_LEN, X_LEN> H,
									Vector<double, U_LEN> u, Vector<double, X_LEN> x_est, Vector<double, Y_LEN> y_est, double dt) {
							   { filter.CalcJacobianF(F, x_est, u, dt) } -> std::same_as<bool>;
							   { filter.NonlinearUpdateX(x_est, x_est, u, dt) } -> std::same_as<bool>;
							   { filter.CalcJacobianH(H, x_est, u) } -> std::same_as<bool>;
							   { filter.NonlinearUpdateY(y_est, x_est, u) } -> std::same_as<bool>;
						   };

	template<class ConcreateFilter, MatrixIdx U_LEN, MatrixIdx X_LEN, MatrixIdx Y_LEN>
	class ExtendedKalmanFilter
	{
	public:
		using InputVector = Vector<double, U_LEN>;
		using StateVector = Vector<double, X_LEN>;
		using OutputVector = Vector<double, Y_LEN>;
		using StateCovMatrix = Matrix<double, X_LEN, X_LEN>;
		using StateToOutputMatrix = Matrix<double, Y_LEN, X_LEN>;
		using OutputToStateMatrix = Matrix<double, X_LEN, Y_LEN>;
		using OutputCovMatrix = Matrix<double, Y_LEN, Y_LEN>;

	public:
		constexpr ExtendedKalmanFilter() {
			static_assert(KalmanFilter<ConcreateFilter, U_LEN, X_LEN, Y_LEN>,
						  "ExtendedKalmanFilter template parameter must implement KalmanFilter concept");
		}

		constexpr void Init(const StateVector& x_init, const StateCovMatrix& P_init) {
			/* Initialization:
			 *  x(k=0|k=0) = Expected value of x at time-0 (i.e. x(k=0)), typically set to zero.
			 *  P(k=0|k=0) = Identity matrix * covariant(P(k=0)), typically initialized with some
			 *               big number.
			 */

			y_est.SetZero();
			y_diff.SetZero();
			x_est = x_init;
			P = P_init;
		}

		constexpr void SetCovarianceMatrices(const StateCovMatrix& Q, const OutputCovMatrix& R) {
			/* Initialization:
			 *  Q, R       = Covariance matrices of process & measurement. As this implementation
			 *               the noise as AWGN (and same value for every variable), this is set
			 *               to Q=diag(QInit,...,QInit) and R=diag(RInit,...,RInit).
			 */
			this->Q = Q;
			this->R = R;
		}

		constexpr bool Update(const OutputVector& y, const InputVector& u, double dt) {
			/* Run once every sampling time */

			/* =============== Calculate the Jacobian matrix of f (i.e. F) =============== */
			/* F = d(f(..))/dx |x(k-1|k-1),u(k-1)                               ...{EKF_1} */
			StateCovMatrix F;
			if(!Self().CalcJacobianF(F, x_est, u, dt)) {
				return false;
			}

			/* =========================== Prediction of x & P =========================== */
			/* x(k|k-1) = f[x(k-1|k-1), u(k-1)]                                 ...{EKF_2} */
			if(!Self().NonlinearUpdateX(x_est, x_est, u, dt)) {
				return false;
			}

			/* P(k|k-1) = F*P(k-1|k-1)*F' + Q                                   ...{EKF_3} */
			P = F * P * F.Transpose() + Q;

			/* =============== Calculate the Jacobian matrix of h (i.e. H) =============== */
			/* H = d(h(..))/dx |x(k|k-1)                                        ...{EKF_4} */
			StateToOutputMatrix H;
			if(!Self().CalcJacobianH(H, x_est, u)) {
				return false;
			}

			/* =========================== Correction of x & P =========================== */
			/* S = H*P(k|k-1)*H' + R                                            ...{EKF_5} */
			const OutputCovMatrix S = H * P * H.Transpose() + R;

			/* K = P(k|k-1)*H'*(S^-1)                                           ...{EKF_6} */
			const OutputToStateMatrix K = P * H.Transpose() * S.Inverse();

			if(!K.IsValid()) {
				return false;
			}

			/* x(k|k) = x(k|k-1) + K*[y(k) - h(x(k|k-1))]                       ...{EKF_7} */
			if(!Self().NonlinearUpdateY(y_est, x_est, u)) {
				return false;
			}

			y_diff = y - y_est;
			x_est = x_est + K * y_diff;

			/* P(k|k) = (I - K*H)*P(k|k-1)                                      ...{EKF_8} */
			P = (StateCovMatrix::Identity() - K * H) * P;

			return true;
		}

	public:
		constexpr const StateVector& GetX() const { return x_est; }
		constexpr const OutputVector& GetY() const { return y_est; }
		constexpr const StateCovMatrix& GetP() const { return P; }
		constexpr const OutputVector& GetErrY() const { return y_diff; }

	protected:
		constexpr ConcreateFilter& Self() { return static_cast<ConcreateFilter&>(*this); }

	protected:
		StateCovMatrix Q;
		OutputCovMatrix R;

	private:
		StateVector x_est;
		OutputVector y_est;
		StateCovMatrix P;
		OutputVector y_diff;
	};
}