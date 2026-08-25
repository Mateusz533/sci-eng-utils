#pragma once
//
#include <gtest/gtest.h>
//
#include "GenericMath/GenericTensor.hpp"

namespace TensorTest
{
	constexpr double EPS = 1e-6;
	constexpr int DIM = 3;

	using namespace GenericMath;

	/************************************************************************************************************/

	TEST(TensorTestSuite, TensorTest) {
		const std::array<std::array<double, 3>, 3> matArray{
			std::array<double, 3>{54, 97, -21},
			std::array<double, 3>{-96, 353, 2},
			std::array<double, 3>{-39, 5, 24},
		};

		const Tensor<double, 2, DIM> mat{matArray};

		const auto matCopy = mat;
		const auto outer = mat.OuterProduct(matCopy);
		const auto fasterInner = mat.InnerProduct(matCopy, 1, 0);
		const auto slowerInner = outer.Contraction(1, 2);

		EXPECT_EQ(mat.Contraction(0, 1), 431);

		for(TensorIdx i = 0; i < DIM; ++i) {
			for(TensorIdx j = 0; j < DIM; ++j) {
				EXPECT_EQ(fasterInner(i, j), slowerInner(i, j));
			}
		}
	}
}