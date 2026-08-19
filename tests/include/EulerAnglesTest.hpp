#pragma once
//
#include <gtest/gtest.h>
//
#include <vector>
//
#include "GenericMath/GenericEulerAngles.hpp"

namespace EulerAnglesTest
{
	constexpr double EPS = 1e-6;
	constexpr long N = 18;

	using namespace GenericMath;

	inline bool AnglesEqual(const StdEulerAnglesd& l, const StdEulerAnglesd& r) {
		return (std::abs(NormalizeRad(l.Roll() - r.Roll())) < EPS &&
				std::abs(NormalizeRad(l.Pitch() - r.Pitch())) < EPS &&
				std::abs(NormalizeRad(l.Yaw() - r.Yaw())) < EPS);
	}

	inline bool AnglesOk(StdEulerAnglesd l, StdEulerAnglesd r) {
		StdEulerAnglesd lDual = l.CalcDualAngles();
		if(std::abs(std::cos(l.Pitch())) < EPS || std::abs(std::cos(r.Pitch())) < EPS) {
			l.Yaw() -= std::sin(l.Pitch()) * l.Roll();
			lDual.Yaw() -= std::sin(lDual.Pitch()) * lDual.Roll();
			l.Roll() = 0.0;
			lDual.Roll() = 0.0;
		}

		return AnglesEqual(l, r) || AnglesEqual(lDual, r);
	}

	inline bool QuaternionEqualAndNormed(const Quaternion<double>& l, const Quaternion<double>& r) {
		return ((std::abs(l.W() - r.W()) < EPS && std::abs(l.X() - r.X()) < EPS && std::abs(l.Y() - r.Y()) < EPS && std::abs(l.Z() - r.Z()) < EPS) ||
				(std::abs(l.W() + r.W()) < EPS && std::abs(l.X() + r.X()) < EPS && std::abs(l.Y() + r.Y()) < EPS && std::abs(l.Z() + r.Z()) < EPS)) &&
			   std::abs(l.W() * l.W() + l.X() * l.X() + l.Y() * l.Y() + l.Z() * l.Z() - 1.0) < EPS;
	}

	inline bool MatrixEqual(const Matrix3<double>& l, const Matrix3<double>& r) {
		for(int i = 0; i < 3; ++i) {
			for(int j = 0; j < 3; ++j) {
				if(std::abs(l(i, j) - r(i, j)) >= EPS)
					return false;
			}
		}
		return true;
	}

	inline auto GetScopeDeg(bool useSingularity = true) {
		std::vector<Vector3d> data;
		data.reserve(8 * N * N * N);
		for(int i = -N; i < N; ++i) {
			const double roll = 180. * i / N;
			for(int j = -N; j < N; ++j) {
				const double pitch = 180. * j / N;

				// singularity
				if(!useSingularity && (N == 2L * j || -N == 2L * j)) {
					continue;
				}

				for(int k = -N; k < N; ++k) {
					const double yaw = 180. * k / N;
					data.emplace_back(roll, pitch, yaw);
				}
			}
		}
		return data;
	}

	/************************************************************************************************************/

	TEST(EulerZYXTestSuite, MatrixConversionTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			const auto anglesZYX = StdEulerAnglesd::FromVectorDegRPY(vec);
			successCounter += AnglesOk(anglesZYX, StdEulerAnglesd::FromMatrix(anglesZYX.ToMatrix()));
		}
		EXPECT_EQ(successCounter, scope.size());
	}

	TEST(EulerZYXTestSuite, MatrixOrderTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			const auto rotZYX = StdEulerAnglesd::FromVectorDegRPY(vec).ToMatrix();
			const auto rotX = StdEulerAnglesd::FromVectorDegRPY({vec.X(), 0.0, 0.0}).ToMatrix();
			const auto rotY = StdEulerAnglesd::FromVectorDegRPY({0.0, vec.Y(), 0.0}).ToMatrix();
			const auto rotZ = StdEulerAnglesd::FromVectorDegRPY({0.0, 0.0, vec.Z()}).ToMatrix();
			successCounter += MatrixEqual(rotZYX, rotZ * rotY * rotX);
		}
		EXPECT_EQ(successCounter, scope.size());
	}

	TEST(EulerZYXTestSuite, QuaternionConversionTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			StdEulerAnglesd anglesZYX = StdEulerAnglesd::FromVectorDegRPY(vec);
			successCounter += AnglesOk(anglesZYX, StdEulerAnglesd::FromQuaternion(anglesZYX.ToQuaternion()));
		}
		EXPECT_EQ(successCounter, scope.size());
	}

	TEST(EulerZYXTestSuite, QuaternionOrderTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			const auto rotZYX = StdEulerAnglesd::FromVectorDegRPY(vec).ToQuaternion();
			const auto rotX = StdEulerAnglesd::FromVectorDegRPY({vec.X(), 0.0, 0.0}).ToQuaternion();
			const auto rotY = StdEulerAnglesd::FromVectorDegRPY({0.0, vec.Y(), 0.0}).ToQuaternion();
			const auto rotZ = StdEulerAnglesd::FromVectorDegRPY({0.0, 0.0, vec.Z()}).ToQuaternion();
			successCounter += QuaternionEqualAndNormed(rotZYX, rotZ * rotY * rotX);
		}
		EXPECT_EQ(successCounter, scope.size());
	}

	TEST(EulerZYXTestSuite, DualAnglesMatrixTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			const auto anglesXYZ = StdEulerAnglesd::FromVectorDegRPY(vec);
			successCounter += MatrixEqual(anglesXYZ.ToMatrix(), anglesXYZ.CalcDualAngles().ToMatrix());
		}
		EXPECT_EQ(successCounter, scope.size());
	}

	TEST(EulerZYXTestSuite, DualAnglesQuaternionTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			const auto anglesXYZ = StdEulerAnglesd::FromVectorDegRPY(vec);
			successCounter += QuaternionEqualAndNormed(anglesXYZ.ToQuaternion(), anglesXYZ.CalcDualAngles().ToQuaternion());
		}
		EXPECT_EQ(successCounter, scope.size());
	}

	TEST(EulerZYXTestSuite, SizeTest) {
		const auto anglesZYX = StdEulerAnglesd::FromVectorDegRPY({0.0, 1.0, 2.0});
		EXPECT_EQ(sizeof(anglesZYX), 3 * sizeof(double));
	}

	/************************************************************************************************************/

	TEST(EulerYZXTestSuite, MatrixOrderTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			const auto rotYZX = EulerAnglesLocalYZX<double>::FromDegrees(vec.Y(), vec.Z(), vec.X()).ToMatrix();
			const auto rotY = EulerAnglesLocalYZX<double>::FromDegrees(vec.Y(), 0.0, 0.0).ToMatrix();
			const auto rotZ = EulerAnglesLocalYZX<double>::FromDegrees(0.0, vec.Z(), 0.0).ToMatrix();
			const auto rotX = EulerAnglesLocalYZX<double>::FromDegrees(0.0, 0.0, vec.X()).ToMatrix();
			successCounter += MatrixEqual(rotYZX, rotY * rotZ * rotX);
		}
		EXPECT_EQ(successCounter, scope.size());
	}

	TEST(EulerYZXTestSuite, QuaternionOrderTest) {
		long successCounter = 0;
		const auto scope = GetScopeDeg();
		for(const auto& vec : scope) {
			const auto rotYZX = EulerAnglesLocalYZX<double>::FromDegrees(vec.Y(), vec.Z(), vec.X()).ToQuaternion();
			const auto rotY = EulerAnglesLocalYZX<double>::FromDegrees(vec.Y(), 0.0, 0.0).ToQuaternion();
			const auto rotZ = EulerAnglesLocalYZX<double>::FromDegrees(0.0, vec.Z(), 0.0).ToQuaternion();
			const auto rotX = EulerAnglesLocalYZX<double>::FromDegrees(0.0, 0.0, vec.X()).ToQuaternion();
			successCounter += QuaternionEqualAndNormed(rotYZX, rotY * rotZ * rotX);
		}
		EXPECT_EQ(successCounter, scope.size());
	}
}