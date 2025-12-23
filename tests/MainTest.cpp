#include "EulerAnglesTest.hpp"
#include "MatrixTest.hpp"
#include "TensorTest.hpp"

int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}