#pragma once

#include <gtest/gtest.h>

#include <array>

#include "GenericMath/GenericMatrix.h"

namespace MatrixTest
{
	constexpr double EPS = 1e-6;
	constexpr int N = 18;

	using namespace GenericMath;

	/************************************************************************************************************/

	TEST(MatrixTestSuite, MatrixInversionQR) {
		const std::array<std::array<double, 3>, 3> matArray = {
			54, 97, -21,
			-96, 353, 2,
			-39, 5, 24};

		const Matrix3d mat = matArray;

		EXPECT_EQ(mat * mat.Inverse(), Matrix3d::Identity());
	}

	TEST(MatrixTestSuite, MatrixInversionLU) {
		const std::array<std::array<double, 3>, 3> matArray = {
			54, 97, -21,
			-96, 353, 2,
			-39, 5, 24};

		const Matrix3d mat = matArray;

		EXPECT_EQ(mat * mat.InverseWithLu(), Matrix3d::Identity());
	}
}