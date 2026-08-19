#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <type_traits>
//
#include "GenericMath/GenericMatrix.hpp"

constexpr std::size_t N = 3;
constexpr std::size_t MAX_SAMPLES = 300'000;

template<class Matrix3d>
inline Matrix3d MakeMatrix(const std::array<double, N * N>& values) {
	Matrix3d res;
	for(std::size_t i = 0; i < N; ++i)
		for(std::size_t j = 0; j < N; ++j)
			res(i, j) = values[i * N + j];

	return res;
}

template<class Matrix3d>
inline auto InverseMatrix(const Matrix3d& mat) {
	if constexpr(std::is_same_v<Matrix3d, Eigen::Matrix3d>)
		return Matrix3d(mat.inverse());
	else if constexpr(std::is_same_v<Matrix3d, GenericMath::Matrix3d>)
		return Matrix3d(mat.Inverse());
	return mat;
}

template<class Matrix3d>
inline bool ValidateMatrix(const Matrix3d& mat) {
	if constexpr(std::is_same_v<Matrix3d, Eigen::Matrix3d>)
		return true;
	else if constexpr(std::is_same_v<Matrix3d, GenericMath::Matrix3d>)
		return mat.IsValid();
	return false;
}

template<class Matrix3d>
inline double MaxDiff(const Matrix3d& mat) {
	double max = 0.0;
	for(std::size_t i = 0; i < N; ++i) {
		for(std::size_t j = 0; j < N; ++j) {
			const double val = std::abs(mat(i, j));
			if(val > max) max = val;
		}
	}
	return max;
}

template<class Matrix3d>
inline void PrintMatrix(const Matrix3d& mat) {
	for(std::size_t i = 0; i < N; ++i) {
		for(std::size_t j = 0; j < N; ++j)
			printf("% 13.6f\t", mat(i, j));
		std::cout << "\n";
	}
}

inline bool IsSingular(const std::array<double, N * N>& matrixArray) {
	bool singular = false;
	for(std::size_t i = 0; i < N; ++i) {
		bool anyNonZeroInRow = false;
		bool anyNonZeroInCol = false;
		for(std::size_t j = 0; j < N; ++j) {
			anyNonZeroInRow = anyNonZeroInRow || matrixArray[i + N * j] != 0;
			anyNonZeroInCol = anyNonZeroInCol || matrixArray[j + N * i] != 0;
		}
		if(!anyNonZeroInRow || !anyNonZeroInCol) {
			singular = true;
			break;
		}
	}
	return singular;
}

template<class Matrix3d>
inline double TestInversionForMatrix(const std::array<double, N * N>& data, int count) {
	try {
		const Matrix3d mat = MakeMatrix<Matrix3d>(data);
		// Matrix3d matInv = mat.InverseWithLu();
		const Matrix3d matInv = InverseMatrix(mat);
		const Matrix3d matIdent = mat * matInv;
		const double diff = MaxDiff(matIdent - Matrix3d::Identity());
		const bool isValid = ValidateMatrix(matInv);

		if(!isValid || diff > 1e-11) {
			std::cout << "Permutation #" << count << ":\n";
			std::cout << "Matrix:\n";
			PrintMatrix(mat);
			std::cout << "Inverse:\n";
			PrintMatrix(matInv);
			std::cout << "Ident:\n";
			PrintMatrix(matIdent);
			std::cout << "Max diff: " << diff << "\n";
			std::cout << "-----------------------------------------------------------\n";
		}
		return isValid ? diff : 0.0;

	} catch(const std::exception& ex) {
		std::cout << "Permutation #" << count << ": singular matrix.\n";
		const GenericMath::Matrix3d mat = MakeMatrix<GenericMath::Matrix3d>(data);
		PrintMatrix(mat);
		std::cout << "Determinant: " << mat.Determinant() << "\n";
		std::cout << "-----------------------------------------------------------\n";
		return 0.0;
	}
}

template<class Matrix3d>
void MatrixGeneralTest() {
	std::array<double, N * N> randomNumbers = {12321, 73434, 5424,
											   8787, 35346, 9037,
											   5824, 9573, 2973};
	std::sort(randomNumbers.begin(), randomNumbers.end());

	const auto startTime = std::chrono::high_resolution_clock::now();

	double maxDiff = 0.0;
	std::size_t count = 0;
	do {
		if(count >= MAX_SAMPLES) break;
		maxDiff = std::max(maxDiff, TestInversionForMatrix<Matrix3d>(randomNumbers, ++count));
	} while(std::next_permutation(randomNumbers.begin(), randomNumbers.end()));

	const auto endTime = std::chrono::high_resolution_clock::now();
	const auto timeMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

	std::cout << "General test:\n";
	std::cout << "\tTotal permutations processed: " << count << "\n";
	std::cout << "\tMax absolute cell difference: " << maxDiff << "\n";
	std::cout << "\tTotal elapsed time [ms]: " << timeMs << "\n";
}

template<class Matrix3d>
void MatrixSpecialTest() {
	std::array<double, N * N> specialNumbers = {1, 0, 0,
												0, 2, 0,
												0, 0, 3};
	std::sort(specialNumbers.begin(), specialNumbers.end());

	double maxDiff = 0.0;
	std::size_t count = 0;
	do {
		if(!IsSingular(specialNumbers)) {
			maxDiff = std::max(maxDiff, TestInversionForMatrix<Matrix3d>(specialNumbers, ++count));
		}
	} while(std::next_permutation(specialNumbers.begin(), specialNumbers.end()));

	std::cout << "Special test:\n";
	std::cout << "\tTotal permutations processed: " << count << "\n";
	std::cout << "\tMax absolute cell difference: " << maxDiff << "\n";
}

int main() {
	std::cout << "-----------------------------------------------------------\n";
	std::cout << "-------------------------- Eigen --------------------------\n";
	std::cout << "-----------------------------------------------------------\n";
	MatrixSpecialTest<Eigen::Matrix3d>();
	std::cout << "-----------------------------------------------------------\n";
	MatrixGeneralTest<Eigen::Matrix3d>();
	std::cout << "-----------------------------------------------------------\n";
	std::cout << "------------------------- Generic -------------------------\n";
	std::cout << "-----------------------------------------------------------\n";
	MatrixSpecialTest<GenericMath::Matrix3d>();
	std::cout << "-----------------------------------------------------------\n";
	MatrixGeneralTest<GenericMath::Matrix3d>();

	return 0;
}