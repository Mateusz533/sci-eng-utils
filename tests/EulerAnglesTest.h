#pragma once

#include <gtest/gtest.h>

#include <vector>

#include "GenericMath/GenericEulerAngles.h"

namespace EulerAnglesTest
{
	constexpr double EPS = 1e-6;
	constexpr int N = 18;

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
		return ((std::abs(l.w() - r.w()) < EPS && std::abs(l.x() - r.x()) < EPS && std::abs(l.y() - r.y()) < EPS && std::abs(l.z() - r.z()) < EPS) ||
				(std::abs(l.w() + r.w()) < EPS && std::abs(l.x() + r.x()) < EPS && std::abs(l.y() + r.y()) < EPS && std::abs(l.z() + r.z()) < EPS)) &&
			   std::abs(l.w() * l.w() + l.x() * l.x() + l.y() * l.y() + l.z() * l.z() - 1.0) < EPS;
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
				if(!useSingularity && (N == 2 * j || -N == 2 * j)) {
					continue;
				}

				for(int k = -N; k < N; ++k) {
					const double yaw = 180. * k / N;
					data.push_back(Vector3d{std::array{roll, pitch, yaw}});
				}
			}
		}
		return std::move(data);
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
			const auto rotX = StdEulerAnglesd::FromVectorDegRPY(Vector3d{std::array{vec.x(), 0.0, 0.0}}).ToMatrix();
			const auto rotY = StdEulerAnglesd::FromVectorDegRPY(Vector3d{std::array{0.0, vec.y(), 0.0}}).ToMatrix();
			const auto rotZ = StdEulerAnglesd::FromVectorDegRPY(Vector3d{std::array{0.0, 0.0, vec.z()}}).ToMatrix();
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
			const auto rotX = StdEulerAnglesd::FromVectorDegRPY(Vector3d{std::array{vec.x(), 0.0, 0.0}}).ToQuaternion();
			const auto rotY = StdEulerAnglesd::FromVectorDegRPY(Vector3d{std::array{0.0, vec.y(), 0.0}}).ToQuaternion();
			const auto rotZ = StdEulerAnglesd::FromVectorDegRPY(Vector3d{std::array{0.0, 0.0, vec.z()}}).ToQuaternion();
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
		const auto anglesZYX = StdEulerAnglesd::FromVectorDegRPY(Vector3d{std::array{0.0, 1.0, 2.0}});
		EXPECT_EQ(sizeof(anglesZYX), 3 * sizeof(double));
	}
}