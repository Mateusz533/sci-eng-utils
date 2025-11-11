#pragma once
//
#include <gtest/gtest.h>
//
#include <array>
//
#include "GenericMath/GenericMatrix.h"

namespace MatrixTest
{
	constexpr double EPS = 1e-6;
	constexpr int N = 18;

	using namespace GenericMath;

	/************************************************************************************************************/

	TEST(MatrixTestSuite, MatrixInversionQRTest) {
		const std::array<std::array<double, 3>, 3> matArray = {
			54, 97, -21,
			-96, 353, 2,
			-39, 5, 24};

		const Matrix3d mat = matArray;

		EXPECT_EQ(mat * mat.Inverse(), Matrix3d::Identity());
	}

	TEST(MatrixTestSuite, MatrixInversionLUTest) {
		const std::array<std::array<double, 3>, 3> matArray = {
			54, 97, -21,
			-96, 353, 2,
			-39, 5, 24};

		const Matrix3d mat = matArray;

		EXPECT_EQ(mat * mat.InverseWithLu(), Matrix3d::Identity());
	}

	TEST(MatrixTestSuite, DynamicMatrixBasicsTest) {
		const MatrixXd mat(5, 7, 10);
		const MatrixXd matInv(7, 5, 0.1);

		EXPECT_EQ(mat.GetRows(), 5);
		EXPECT_EQ(mat.GetCols(), 7);
		EXPECT_EQ(mat(1, 3), 10);

		EXPECT_TRUE((mat + mat).IsValid());
		EXPECT_FALSE((mat * mat).IsValid());
		EXPECT_TRUE((mat * matInv).IsValid());
	}

	TEST(MatrixTestSuite, DynamicMatrixOperatorsTest) {
		const MatrixXd mat(5, 7, 10);
		const MatrixXd matInv(7, 5, 0.1);

		EXPECT_EQ((mat + mat)(4, 2), 20);
		EXPECT_EQ((mat - mat)(3, 1), 0);
		EXPECT_EQ((mat / 5)(2, 1), 2);
		EXPECT_EQ((7 * mat)(1, 5), 70);
		EXPECT_EQ((mat * matInv)(3, 3), 7);
		EXPECT_EQ((mat * matInv)(4, 1), 7);

		EXPECT_FALSE((mat * mat).IsValid());
	}
}