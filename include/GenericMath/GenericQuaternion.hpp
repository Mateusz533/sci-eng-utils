#pragma once
//
#include <cmath>
#include <concepts>
//
#include "GenericMatrix.hpp"

namespace GenericMath
{
	template<std::floating_point T>
	class Quaternion
	{
	public:
		constexpr Quaternion() : Quaternion{Identity()} {}
		constexpr Quaternion(T w, T x, T y, T z) : mScalar{w}, mVector{x, y, z} {}
		constexpr Quaternion(const Quaternion& other) = default;
		constexpr Quaternion(Quaternion&& other) = default;

		template<std::floating_point U>
		Quaternion(const Matrix3<U>& rotationMatrix) {
			const auto& mat = rotationMatrix;
			const double trace = mat(0, 0) + mat(1, 1) + mat(2, 2);

			if(trace > 0.0) {
				const double s = 2.0 * std::sqrt(trace + 1.0);
				W() = 0.25 * s;
				X() = (mat(2, 1) - mat(1, 2)) / s;
				Y() = (mat(0, 2) - mat(2, 0)) / s;
				Z() = (mat(1, 0) - mat(0, 1)) / s;
			} else if(mat(0, 0) > mat(1, 1) && mat(0, 0) > mat(2, 2)) {
				const double s = 2.0 * std::sqrt(1.0 + mat(0, 0) - mat(1, 1) - mat(2, 2));
				W() = (mat(2, 1) - mat(1, 2)) / s;
				X() = 0.25 * s;
				Y() = (mat(0, 1) + mat(1, 0)) / s;
				Z() = (mat(0, 2) + mat(2, 0)) / s;
			} else if(mat(1, 1) > mat(2, 2)) {
				const double s = 2.0 * std::sqrt(1.0 + mat(1, 1) - mat(0, 0) - mat(2, 2));
				W() = (mat(0, 2) - mat(2, 0)) / s;
				X() = (mat(0, 1) + mat(1, 0)) / s;
				Y() = 0.25 * s;
				Z() = (mat(1, 2) + mat(2, 1)) / s;
			} else {
				const double s = 2.0 * std::sqrt(1.0 + mat(2, 2) - mat(0, 0) - mat(1, 1));
				W() = (mat(1, 0) - mat(0, 1)) / s;
				X() = (mat(0, 2) + mat(2, 0)) / s;
				Y() = (mat(1, 2) + mat(2, 1)) / s;
				Z() = 0.25 * s;
			}
		}

	public:
		constexpr Quaternion& operator=(const Quaternion& other) = default;
		constexpr Quaternion& operator=(Quaternion&& other) = default;

		constexpr Quaternion& operator*=(const Quaternion& other) {
			return *this = *this * other;
		}
		constexpr Quaternion& operator/=(const Quaternion& other) {
			return *this = *this / other;
		}

		constexpr Quaternion operator*(const Quaternion& other) const {
			return Quaternion{
				W() * other.W() - X() * other.X() - Y() * other.Y() - Z() * other.Z(),
				W() * other.X() + X() * other.W() + Y() * other.Z() - Z() * other.Y(),
				W() * other.Y() + Y() * other.W() + Z() * other.X() - X() * other.Z(),
				W() * other.Z() + Z() * other.W() + X() * other.Y() - Y() * other.X(),
			};
		}
		constexpr Quaternion operator/(const Quaternion& other) const {
			return *this * other.Inverse();
		}

		template<std::floating_point U>
		operator Quaternion<U>() const {
			return Quaternion<U>{W(), X(), Y(), Z()};
		}

	public:
		constexpr T& W() { return mScalar; }
		constexpr T& X() { return mVector(0); }
		constexpr T& Y() { return mVector(1); }
		constexpr T& Z() { return mVector(2); }

		constexpr const T& W() const { return mScalar; }
		constexpr const T& X() const { return mVector(0); }
		constexpr const T& Y() const { return mVector(1); }
		constexpr const T& Z() const { return mVector(2); }

		constexpr T& Scalar() { return mScalar; }
		constexpr const T& Scalar() const { return mScalar; }
		constexpr Vector3<T>& Vector() { return mVector; }
		constexpr const Vector3<T>& Vector() const { return mVector; }

	public:
		static constexpr Quaternion<T> Identity() {
			return Quaternion<T>{1, 0, 0, 0};
		}

		constexpr Quaternion Normalize() const {
			const double normInv = 1.0 / Norm();
			return Quaternion{
				W() * normInv,
				X() * normInv,
				Y() * normInv,
				Z() * normInv,
			};
		}

		constexpr Quaternion Inverse() const {
			const double factor = 1.0 / SquareSum();
			return Quaternion{
				+W() * factor,
				-X() * factor,
				-Y() * factor,
				-Z() * factor,
			};
		}

		constexpr Matrix3<T> ToRotationMatrix() const {
			Matrix3<T> result;

			const double wSq = W() * W();
			const double xSq = X() * X();
			const double ySq = Y() * Y();
			const double zSq = Z() * Z();

			const double invSqNorm = 1.0 / (xSq + ySq + zSq + wSq);

			result(0, 0) = (+xSq - ySq - zSq + wSq);
			result(1, 1) = (-xSq + ySq - zSq + wSq);
			result(2, 2) = (-xSq - ySq + zSq + wSq);

			result(1, 0) = 2.0 * (X() * Y() + Z() * W());
			result(0, 1) = 2.0 * (X() * Y() - Z() * W());

			result(2, 0) = 2.0 * (X() * Z() - Y() * W());
			result(0, 2) = 2.0 * (X() * Z() + Y() * W());

			result(2, 1) = 2.0 * (Y() * Z() + X() * W());
			result(1, 2) = 2.0 * (Y() * Z() - X() * W());

			return result * invSqNorm;
		}

		constexpr double Norm() const {
			return std::sqrt(SquareSum());
		}

	private:
		constexpr double SquareSum() const {
			return X() * X() + Y() * Y() + Z() * Z() + W() * W();
		}

	private:
		T mScalar;
		Vector3<T> mVector;
	};

	using Quaternionf = Quaternion<float>;
	using Quaterniond = Quaternion<double>;
}