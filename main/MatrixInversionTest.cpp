#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <type_traits>
//
#include "GenericMath/GenericMatrix.hpp"

constexpr int N = 3;
constexpr int MAX_SAMPLES = 300'000;

template<class Matrix3d>
inline Matrix3d makeMatrix(const std::array<double, N * N>& values) {
	Matrix3d res;
	for(size_t i = 0; i < N; ++i)
		for(size_t j = 0; j < N; ++j)
			res(i, j) = values[i * N + j];

	return res;
}

template<class Matrix3d>
inline auto inverseMatrix(const Matrix3d& mat) {
	if constexpr(std::is_same_v<Matrix3d, Eigen::Matrix3d>)
		return Matrix3d(mat.inverse());
	else if constexpr(std::is_same_v<Matrix3d, GenericMath::Matrix3d>)
		return Matrix3d(mat.Inverse());
	return mat;
}

template<class Matrix3d>
inline bool validateMatrix(const Matrix3d& mat) {
	if constexpr(std::is_same_v<Matrix3d, Eigen::Matrix3d>)
		return true;
	else if constexpr(std::is_same_v<Matrix3d, GenericMath::Matrix3d>)
		return mat.IsValid();
	return false;
}

template<class Matrix3d>
inline double maxDiff(const Matrix3d& mat) {
	double max = 0.0;
	for(size_t i = 0; i < N; ++i) {
		for(size_t j = 0; j < N; ++j) {
			const double val = std::abs(mat(i, j));
			if(val > max) max = val;
		}
	}
	return max;
}

template<class Matrix3d>
inline void printMatrix(const Matrix3d& mat) {
	for(size_t i = 0; i < N; ++i) {
		for(size_t j = 0; j < N; ++j)
			printf("% 13.6f\t", mat(i, j));
		std::cout << "\n";
	}
}

inline bool isSingular(const std::array<double, N * N>& matrixArray) {
	bool singular = false;
	for(int i = 0; i < N; ++i) {
		bool anyNonZeroInRow = false;
		bool anyNonZeroInCol = false;
		for(int j = 0; j < N; ++j) {
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
inline double testInversionForMatrix(const std::array<double, N * N>& data, int count) {
	try {
		const Matrix3d mat = makeMatrix<Matrix3d>(data);
		// Matrix3d matInv = mat.InverseWithLu();
		const Matrix3d matInv = inverseMatrix(mat);
		const Matrix3d matIdent = mat * matInv;
		const double diff = maxDiff(matIdent - Matrix3d::Identity());
		const bool isValid = validateMatrix(matInv);

		if(!isValid || diff > 1e-11) {
			std::cout << "Permutation #" << count << ":\n";
			std::cout << "Matrix:\n";
			printMatrix(mat);
			std::cout << "Inverse:\n";
			printMatrix(matInv);
			std::cout << "Ident:\n";
			printMatrix(matIdent);
			std::cout << "Max diff: " << diff << "\n";
			std::cout << "-----------------------------------------------------------\n";
		}
		return isValid ? diff : 0.0;

	} catch(const std::exception& ex) {
		std::cout << "Permutation #" << count << ": singular matrix.\n";
		const GenericMath::Matrix3d mat = makeMatrix<GenericMath::Matrix3d>(data);
		printMatrix(mat);
		std::cout << "Determinant: " << mat.Determinant() << "\n";
		std::cout << "-----------------------------------------------------------\n";
		return 0.0;
	}
}

template<class Matrix3d>
void matrixGeneralTest() {
	std::array<double, N * N> randomNumbers = {12321, 73434, 5424,
											   8787, 35346, 9037,
											   5824, 9573, 2973};
	std::sort(randomNumbers.begin(), randomNumbers.end());

	const auto startTime = std::chrono::high_resolution_clock::now();

	double maxDiff = 0.0;
	int count = 0;
	do {
		if(count >= MAX_SAMPLES) break;
		maxDiff = std::max(maxDiff, testInversionForMatrix<Matrix3d>(randomNumbers, ++count));
	} while(std::next_permutation(randomNumbers.begin(), randomNumbers.end()));

	const auto endTime = std::chrono::high_resolution_clock::now();
	const auto timeMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

	std::cout << "General test:\n";
	std::cout << "\tTotal permutations processed: " << count << "\n";
	std::cout << "\tMax absolute cell difference: " << maxDiff << "\n";
	std::cout << "\tTotal elapsed time [ms]: " << timeMs << "\n";
}

template<class Matrix3d>
void matrixSpecialTest() {
	std::array<double, N * N> specialNumbers = {1, 0, 0,
												0, 2, 0,
												0, 0, 3};
	std::sort(specialNumbers.begin(), specialNumbers.end());

	double maxDiff = 0.0;
	int count = 0;
	do {
		if(!isSingular(specialNumbers)) {
			maxDiff = std::max(maxDiff, testInversionForMatrix<Matrix3d>(specialNumbers, ++count));
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
	matrixSpecialTest<Eigen::Matrix3d>();
	std::cout << "-----------------------------------------------------------\n";
	matrixGeneralTest<Eigen::Matrix3d>();
	std::cout << "-----------------------------------------------------------\n";
	std::cout << "------------------------- Generic -------------------------\n";
	std::cout << "-----------------------------------------------------------\n";
	matrixSpecialTest<GenericMath::Matrix3d>();
	std::cout << "-----------------------------------------------------------\n";
	matrixGeneralTest<GenericMath::Matrix3d>();

	return 0;
}