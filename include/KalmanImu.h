/*************************************************************************************************************
 * Class for Extended Kalman Filter to estimating orientation from IMU data
 *
 * Based on implementation from: https://github.com/pronenewbits
 ************************************************************************************************************/
#pragma once

#include "ExtendedKalmanFilter.h"
#include "GenericMath/GenericQuaternion.h"

namespace GenericMath
{
	/* State Space dimension */

	constexpr MatrixIdx U_LEN = 3;
	constexpr MatrixIdx X_LEN = 4;
	constexpr MatrixIdx Y_LEN = 6;

	class KalmanImu : public ExtendedKalmanFilter<KalmanImu, U_LEN, X_LEN, Y_LEN>
	{
	private:
		inline static constexpr float Q_init = 1e-6;
		inline static constexpr float R_initAcc = 1e-1;
		inline static constexpr float R_initMag = 1e-1;
		inline static constexpr double ACC_Z0 = 1.0;
		inline static constexpr double P_init = 1.0e-8;

	private:
		Vector3d magneticDir{1.0, 0.0, 0.0};

	public:
		constexpr KalmanImu() : ExtendedKalmanFilter() {}

		constexpr bool Init() {
			constexpr StateVector identityQ{1.0, 0.0, 0.0, 0.0};

			ExtendedKalmanFilter::Init(identityQ, StateCovMatrix::Identity() * P_init);
			SetCovarianceMatrices(Q_init, R_initAcc, R_initMag);

			return true;
		}

		constexpr void SetCovarianceMatrices(double Q_init, double R_initAcc, double R_initMag) {
			Q.SetDiagonal(Q_init);
			R.SetDiagonal({R_initAcc, R_initAcc, R_initAcc, R_initMag, R_initMag, R_initMag});
		}

		constexpr void SetEarthMagneticVector(const Vector3d& magneticDir) {
			this->magneticDir = magneticDir;
		}

		constexpr bool Update(const Vector3d& gyroRps, Vector3d accG, Vector3d mag, double yawRad, double dT) {
			InputVector u;

			/* Input 1:3 = gyroscope */

			u(0) = gyroRps.x();
			u(1) = gyroRps.y();
			u(2) = gyroRps.z();

			OutputVector y;

			/* Output 1:3 = accelerometer */

			accG.Normalize();

			y(0) = accG.x();
			y(1) = accG.y();
			y(2) = accG.z();

			/* Output 4:6 = magnetometer */

			mag.Normalize();

			y(3) = mag.x();
			y(4) = mag.y();
			y(5) = mag.z();

			/* Calculate */

			if(!ExtendedKalmanFilter::Update(y, u, dT)) {
				Init();
				return false;
			}

			return true;
		}

		constexpr Quaterniond GetQuaternion() const {
			const StateVector& quat = GetX();
			return Quaterniond(quat(0), quat(1), quat(2), quat(3));
		}

		/* Nonlinear & linearization function ------------------------------------------------------------------------------- */

		static constexpr bool NonlinearUpdateX(StateVector& x_next, const StateVector& x, const InputVector& u, double dt) {
			/* Insert the nonlinear update transformation here
			 *          x(k+1) = f[x(k), u(k)]
			 *
			 * The quaternion update function:
			 *  q0_dot = 1/2 * (  0   - u0*q1 - u1*q2 - u2*q3)
			 *  q1_dot = 1/2 * (u0*q0 +   0   + u2*q2 - u1*q3)
			 *  q2_dot = 1/2 * (u1*q0 - u2*q1 +  0    + u0*q3)
			 *  q3_dot = 1/2 * (u2*q0 + u1*q1 - u0*q2 +   0  )
			 *
			 * Euler method for integration:
			 *  q0 = q0 + q0_dot * dT;
			 *  q1 = q1 + q1_dot * dT;
			 *  q2 = q2 + q2_dot * dT;
			 *  q3 = q3 + q3_dot * dT;
			 */

			const auto& q0 = x(0);
			const auto& q1 = x(1);
			const auto& q2 = x(2);
			const auto& q3 = x(3);

			const auto& u0 = u(0);
			const auto& u1 = u(1);
			const auto& u2 = u(2);

			x_next(0) = q0 + 0.5 * (+0.00000 - u0 * q1 - u1 * q2 - u2 * q3) * dt;
			x_next(1) = q1 + 0.5 * (+u0 * q0 + 0.00000 + u2 * q2 - u1 * q3) * dt;
			x_next(2) = q2 + 0.5 * (+u1 * q0 - u2 * q1 + 0.00000 + u0 * q3) * dt;
			x_next(3) = q3 + 0.5 * (+u2 * q0 + u1 * q1 - u0 * q2 + 0.00000) * dt;

			/* ======= Additional ad-hoc quaternion normalization to make sure the quaternion is a unit vector (i.e. ||q|| = 1) ======= */
			if(!x_next.Normalize()) {
				/* System error, return false so we can reset the UKF */
				return false;
			}

			return true;
		}

		constexpr bool NonlinearUpdateY(OutputVector& y_est, const StateVector& x, const InputVector& u) {
			/* Insert the nonlinear measurement transformation here
			 *          y(k)   = h[x(k), u(k)]
			 *
			 * The measurement output is the gravitational and magnetic projection to the body:
			 *     DCM     = [(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                    2*(q1*q2+q0*q3),                    2*(q1*q3-q0*q2)]
			 *               [                   2*(q1*q2-q0*q3), (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                    2*(q2*q3+q0*q1)]
			 *               [                   2*(q1*q3+q0*q2),                    2*(q2*q3-q0*q1), (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]
			 *
			 *  G_proj_sens = DCM * [0 0 1]             --> Gravitational projection to the accelerometer sensor
			 *  M_proj_sens = DCM * [Mx My Mz]          --> (Earth) magnetic projection to the magnetometer sensor
			 */

			const auto& q0 = x(0);
			const auto& q1 = x(1);
			const auto& q2 = x(2);
			const auto& q3 = x(3);

			const auto q0_2 = q0 * q0;
			const auto q1_2 = q1 * q1;
			const auto q2_2 = q2 * q2;
			const auto q3_2 = q3 * q3;

			y_est(0) = (2.0 * (q1 * q3 - q0 * q2)) * ACC_Z0;
			y_est(1) = (2.0 * (q2 * q3 + q0 * q1)) * ACC_Z0;
			y_est(2) = (q0_2 - q1_2 - q2_2 + q3_2) * ACC_Z0;

			y_est(3) = (q0_2 + q1_2 - q2_2 - q3_2) * magneticDir(0) +
					   (2.0 * (q1 * q2 + q0 * q3)) * magneticDir(1) +
					   (2.0 * (q1 * q3 - q0 * q2)) * magneticDir(2);

			y_est(4) = (2.0 * (q1 * q2 - q0 * q3)) * magneticDir(0) +
					   (q0_2 - q1_2 + q2_2 - q3_2) * magneticDir(1) +
					   (2.0 * (q2 * q3 + q0 * q1)) * magneticDir(2);

			y_est(5) = (2.0 * (q1 * q3 + q0 * q2)) * magneticDir(0) +
					   (2.0 * (q2 * q3 - q0 * q1)) * magneticDir(1) +
					   (q0_2 - q1_2 - q2_2 + q3_2) * magneticDir(2);

			return true;
		}

		static constexpr bool CalcJacobianF(StateCovMatrix& F, const StateVector& x, const InputVector& u, double dt) {
			/* In UpdateNonlinearX():
			 *  q0 = q0 + q0_dot * dT;
			 *  q1 = q1 + q1_dot * dT;
			 *  q2 = q2 + q2_dot * dT;
			 *  q3 = q3 + q3_dot * dT;
			 */

			const auto& u0 = u(0);
			const auto& u1 = u(1);
			const auto& u2 = u(2);

			F(0, 0) = 1.000;
			F(1, 0) = +0.5 * u0 * dt;
			F(2, 0) = +0.5 * u1 * dt;
			F(3, 0) = +0.5 * u2 * dt;

			F(0, 1) = -0.5 * u0 * dt;
			F(1, 1) = 1.000;
			F(2, 1) = -0.5 * u2 * dt;
			F(3, 1) = +0.5 * u1 * dt;

			F(0, 2) = -0.5 * u1 * dt;
			F(1, 2) = +0.5 * u2 * dt;
			F(2, 2) = 1.000;
			F(3, 2) = -0.5 * u0 * dt;

			F(0, 3) = -0.5 * u2 * dt;
			F(1, 3) = -0.5 * u1 * dt;
			F(2, 3) = +0.5 * u0 * dt;
			F(3, 3) = 1.000;

			return true;
		}

		constexpr bool CalcJacobianH(StateToOutputMatrix& H, const StateVector& x, const InputVector& u) const {
			/* In UpdateNonlinearY():
			 *
			 * The measurement output is the gravitational and magnetic projection to the body:
			 *     DCM     = [(+(q0**2)+(q1**2)-(q2**2)-(q3**2)),                    2*(q1*q2+q0*q3),                    2*(q1*q3-q0*q2)]
			 *               [                   2*(q1*q2-q0*q3), (+(q0**2)-(q1**2)+(q2**2)-(q3**2)),                    2*(q2*q3+q0*q1)]
			 *               [                   2*(q1*q3+q0*q2),                    2*(q2*q3-q0*q1), (+(q0**2)-(q1**2)-(q2**2)+(q3**2))]
			 *
			 *  G_proj_sens = DCM * [0 0 -g]            --> Gravitational projection to the accelerometer sensor
			 *  M_proj_sens = DCM * [Mx My Mz]          --> (Earth) magnetic projection to the magnetometer sensor
			 */

			const auto& q0 = x(0);
			const auto& q1 = x(1);
			const auto& q2 = x(2);
			const auto& q3 = x(3);

			H(0, 0) = -2 * q2 * ACC_Z0;
			H(1, 0) = +2 * q1 * ACC_Z0;
			H(2, 0) = +2 * q0 * ACC_Z0;
			H(3, 0) = +2 * q0 * magneticDir(0) + 2 * q3 * magneticDir(1) - 2 * q2 * magneticDir(2);
			H(4, 0) = -2 * q3 * magneticDir(0) + 2 * q0 * magneticDir(1) + 2 * q1 * magneticDir(2);
			H(5, 0) = +2 * q2 * magneticDir(0) - 2 * q1 * magneticDir(1) + 2 * q0 * magneticDir(2);

			H(0, 1) = +2 * q3 * ACC_Z0;
			H(1, 1) = +2 * q0 * ACC_Z0;
			H(2, 1) = -2 * q1 * ACC_Z0;
			H(3, 1) = +2 * q1 * magneticDir(0) + 2 * q2 * magneticDir(1) + 2 * q3 * magneticDir(2);
			H(4, 1) = +2 * q2 * magneticDir(0) - 2 * q1 * magneticDir(1) + 2 * q0 * magneticDir(2);
			H(5, 1) = +2 * q3 * magneticDir(0) - 2 * q0 * magneticDir(1) - 2 * q1 * magneticDir(2);

			H(0, 2) = -2 * q0 * ACC_Z0;
			H(1, 2) = +2 * q3 * ACC_Z0;
			H(2, 2) = -2 * q2 * ACC_Z0;
			H(3, 2) = -2 * q2 * magneticDir(0) + 2 * q1 * magneticDir(1) - 2 * q0 * magneticDir(2);
			H(4, 2) = +2 * q1 * magneticDir(0) + 2 * q2 * magneticDir(1) + 2 * q3 * magneticDir(2);
			H(5, 2) = +2 * q0 * magneticDir(0) + 2 * q3 * magneticDir(1) - 2 * q2 * magneticDir(2);

			H(0, 3) = +2 * q1 * ACC_Z0;
			H(1, 3) = +2 * q2 * ACC_Z0;
			H(2, 3) = +2 * q3 * ACC_Z0;
			H(3, 3) = -2 * q3 * magneticDir(0) + 2 * q0 * magneticDir(1) + 2 * q1 * magneticDir(2);
			H(4, 3) = -2 * q0 * magneticDir(0) - 2 * q3 * magneticDir(1) + 2 * q2 * magneticDir(2);
			H(5, 3) = +2 * q1 * magneticDir(0) + 2 * q2 * magneticDir(1) + 2 * q3 * magneticDir(2);

			return true;
		}
	};
}
