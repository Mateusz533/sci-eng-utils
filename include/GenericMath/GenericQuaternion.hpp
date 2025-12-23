#pragma once
//
#include <cmath>
//
#include "GenericMatrix.hpp"

namespace GenericMath
{
	template<typename T>
	class Quaternion
	{
	public:
		constexpr Quaternion() : Quaternion{Identity()} {}
		constexpr Quaternion(T w, T x, T y, T z) : mScalar{w}, mVector{x, y, z} {}
		constexpr Quaternion(const Quaternion &other) = default;

		template<typename U>
		Quaternion(const Matrix3<U> &rotationMatrix) {
			const auto &mat = rotationMatrix;
			const double trace = mat(0, 0) + mat(1, 1) + mat(2, 2);

			if(trace > 0.0) {
				const double s = 2.0 * std::sqrt(trace + 1.0);
				w() = 0.25 * s;
				x() = (mat(2, 1) - mat(1, 2)) / s;
				y() = (mat(0, 2) - mat(2, 0)) / s;
				z() = (mat(1, 0) - mat(0, 1)) / s;
			} else if(mat(0, 0) > mat(1, 1) && mat(0, 0) > mat(2, 2)) {
				const double s = 2.0 * std::sqrt(1.0 + mat(0, 0) - mat(1, 1) - mat(2, 2));
				w() = (mat(2, 1) - mat(1, 2)) / s;
				x() = 0.25 * s;
				y() = (mat(0, 1) + mat(1, 0)) / s;
				z() = (mat(0, 2) + mat(2, 0)) / s;
			} else if(mat(1, 1) > mat(2, 2)) {
				const double s = 2.0 * std::sqrt(1.0 + mat(1, 1) - mat(0, 0) - mat(2, 2));
				w() = (mat(0, 2) - mat(2, 0)) / s;
				x() = (mat(0, 1) + mat(1, 0)) / s;
				y() = 0.25 * s;
				z() = (mat(1, 2) + mat(2, 1)) / s;
			} else {
				const double s = 2.0 * std::sqrt(1.0 + mat(2, 2) - mat(0, 0) - mat(1, 1));
				w() = (mat(1, 0) - mat(0, 1)) / s;
				x() = (mat(0, 2) + mat(2, 0)) / s;
				y() = (mat(1, 2) + mat(2, 1)) / s;
				z() = 0.25 * s;
			}
		}

	public:
		constexpr Quaternion &operator=(const Quaternion &other) = default;

		constexpr Quaternion &operator*=(const Quaternion &other) {
			return *this = *this * other;
		}
		constexpr Quaternion &operator/=(const Quaternion &other) {
			return *this = *this / other;
		}

		constexpr Quaternion operator*(const Quaternion &other) const {
			return Quaternion{
				w() * other.w() - x() * other.x() - y() * other.y() - z() * other.z(),
				w() * other.x() + x() * other.w() + y() * other.z() - z() * other.y(),
				w() * other.y() + y() * other.w() + z() * other.x() - x() * other.z(),
				w() * other.z() + z() * other.w() + x() * other.y() - y() * other.x(),
			};
		}
		constexpr Quaternion operator/(const Quaternion &other) const {
			return *this * other.Inverse();
		}

		template<typename U>
		operator Quaternion<U>() const {
			return Quaternion<U>{w(), x(), y(), z()};
		}

	public:
		constexpr T &w() { return mScalar; }
		constexpr T &x() { return mVector(0); }
		constexpr T &y() { return mVector(1); }
		constexpr T &z() { return mVector(2); }

		constexpr const T &w() const { return mScalar; }
		constexpr const T &x() const { return mVector(0); }
		constexpr const T &y() const { return mVector(1); }
		constexpr const T &z() const { return mVector(2); }

		constexpr T &Scalar() { return mScalar; }
		constexpr const T &Scalar() const { return mScalar; }
		constexpr Vector3<T> &Vector() { return mVector; }
		constexpr const Vector3<T> &Vector() const { return mVector; }

	public:
		static constexpr Quaternion<T> Identity() {
			return Quaternion<T>{1, 0, 0, 0};
		}

		constexpr Quaternion Normalize() const {
			const double normInv = 1.0 / Norm();
			return Quaternion{
				w() * normInv,
				x() * normInv,
				y() * normInv,
				z() * normInv,
			};
		}

		constexpr Quaternion Inverse() const {
			const double factor = 1.0 / SquareSum();
			return Quaternion{
				+w() * factor,
				-x() * factor,
				-y() * factor,
				-z() * factor,
			};
		}

		constexpr Matrix3<T> ToRotationMatrix() const {
			Matrix3<T> result;

			const double wSq = w() * w();
			const double xSq = x() * x();
			const double ySq = y() * y();
			const double zSq = z() * z();

			const double invSqNorm = 1.0 / (xSq + ySq + zSq + wSq);

			result(0, 0) = (+xSq - ySq - zSq + wSq);
			result(1, 1) = (-xSq + ySq - zSq + wSq);
			result(2, 2) = (-xSq - ySq + zSq + wSq);

			result(1, 0) = 2.0 * (x() * y() + z() * w());
			result(0, 1) = 2.0 * (x() * y() - z() * w());

			result(2, 0) = 2.0 * (x() * z() - y() * w());
			result(0, 2) = 2.0 * (x() * z() + y() * w());

			result(2, 1) = 2.0 * (y() * z() + x() * w());
			result(1, 2) = 2.0 * (y() * z() - x() * w());

			return result * invSqNorm;
		}

		constexpr double Norm() const {
			return std::sqrt(SquareSum());
		}

	private:
		constexpr double SquareSum() const {
			return x() * x() + y() * y() + z() * z() + w() * w();
		}

	private:
		T mScalar;
		Vector3<T> mVector;
	};

	using Quaternionf = Quaternion<float>;
	using Quaterniond = Quaternion<double>;
}