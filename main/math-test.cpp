#include <iostream>

#include "GenericMatrix.h"

int main() {
	/* Quick test for matrix operations */

	GenericMath::Matrix<float, 3, 1> v{};
	std::cout << "Random vector 3D:  ";
	v.Print();

	GenericMath::Matrix<float, 2, 3> m{};
	std::cout << "Random matrix 2x3: ";
	m.Print();

	GenericMath::Matrix3d a;
	std::cout << "Random matrix 3x3: ";
	a.Print();

	GenericMath::Matrix3d b = {1.0};
	std::cout << "Ones:              ";
	b.Print();
	std::cout << "Threes:            ";
	(b * b).Print();

	GenericMath::Matrix3d c = b * GenericMath::Matrix3d::Zero();
	std::cout << "Zeros:             ";
	c.Print();

	b.SetInvalid();
	std::cout << "Nans:              ";
	b.Print();
}
