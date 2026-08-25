#include <algorithm>
#include <iostream>
#include <iterator>
//
#include "GenericMath/GenericMatrix.hpp"

template<typename T>
	requires(std::is_trivially_destructible_v<T>)
inline bool HasEconomicalDefaultInitialization() noexcept {
	constexpr auto RANDOM_BYTE = std::byte{0x33};

	alignas(T) std::byte scope[sizeof(T)];
	std::fill(std::begin(scope), std::end(scope), RANDOM_BYTE);

	::new(static_cast<void*>(scope)) T;
	return std::all_of(std::begin(scope), std::end(scope), [](const std::byte& byte) {
		// std::cout << int(byte) << ' ';
		return byte == RANDOM_BYTE;
	});
}

int main() {
	/* Quick test for matrix operations */

	std::cout << "Economical default-initialization vector 3D:  " << HasEconomicalDefaultInitialization<GenericMath::Matrix<float, 3, 1>>() << '\n';
	std::cout << "Economical default-initialization matrix 2x3: " << HasEconomicalDefaultInitialization<GenericMath::Matrix<float, 2, 3>>() << '\n';
	std::cout << "Economical default-initialization matrix 3x3: " << HasEconomicalDefaultInitialization<GenericMath::Matrix3d>() << '\n';
	GenericMath::Matrix3d a{};
	std::cout << "Zeros:   ";
	a.Print();

	GenericMath::Matrix3d b = {1.0};
	std::cout << "Ones:    ";
	b.Print();
	std::cout << "Threes:  ";
	(b * b).Print();

	GenericMath::Matrix3d c = b * GenericMath::Matrix3d::Zero();
	std::cout << "Zeros:   ";
	c.Print();

	b.SetInvalid();
	std::cout << "Invalid: ";
	b.Print();
}