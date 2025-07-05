#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <tuple>
#include <type_traits>
//
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#define EIGEN_AVAILABLE

namespace GenericMath
{
	using MatrixIdx = int16_t;

	namespace _
	{
		inline constexpr MatrixIdx MAX_MATRIX_ALLOCATION_SIZE = 6;

		template<typename>
		struct SiblingTrait;
	}

	template<class DataClass>
	class AbstractMatrix;

	template<class DataClass>
	class AbstractSquareMatrix;

	template<class DataClass>
	class AbstractVector;

	template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
	class Matrix;

	/* ================================================================================================== */

	template<class DataClass>
	class AbstractMatrix
	{
	protected:
		constexpr AbstractMatrix() {
			static_assert(std::is_base_of_v<AbstractMatrix, DataClass>,
						  "AbstractMatrix template parameter must derive from AbstractMatrix");
		}
		constexpr ~AbstractMatrix() = default;

		using Matrix = DataClass;
		using Idx = MatrixIdx;
		using T = typename _::SiblingTrait<DataClass>::RawType;
		inline static constexpr Idx ROWS = _::SiblingTrait<DataClass>::R;
		inline static constexpr Idx COLS = _::SiblingTrait<DataClass>::C;

		template<typename U, Idx R = ROWS, Idx C = COLS>
		using Sibling = typename _::SiblingTrait<DataClass>::template Type<U, R, C>;

		constexpr Matrix &Self() { return static_cast<Matrix &>(*this); }
		constexpr const Matrix &Self() const { return static_cast<const Matrix &>(*this); }
		constexpr T &Data(Idx row, Idx col) { return Self()(row, col); }
		constexpr const T &Data(Idx row, Idx col) const { return Self()(row, col); }

		static constexpr bool IsCloseToZeroEco(T value) {
			return (std::abs(value) < T(1e-8));
		}

		template<typename Arg>
		static constexpr auto ExtractIfMatrixLike(const Arg &v, auto r, auto c) {
			if constexpr(std::is_arithmetic_v<Arg>) {
				return v;
			} else {
				return v(r, c);
			}
		}
		template<auto Function, typename... Args>
		constexpr void ForEachElementCall(const Args &...args) const {
			for(Idx row = 0; row < Self().GetRows(); ++row)
				for(Idx col = 0; col < Self().GetCols(); ++col)
					Function(ExtractIfMatrixLike(args, row, col)...);
		}
		template<auto Function, typename... Args>
		constexpr void ForEachElementAssign(const Args &...args) {
			for(Idx row = 0; row < Self().GetRows(); ++row)
				for(Idx col = 0; col < Self().GetCols(); ++col)
					Data(row, col) = Function(ExtractIfMatrixLike(args, row, col)...);
		}
		template<auto Function, typename U, typename... Args>
		constexpr Matrix &CalcFrom(const Sibling<U> &sizeDescriptor, const Args &...args) {
			Self().ResizeIfDynamic(sizeDescriptor);
			ForEachElementAssign<Function>(args...);
			return Self();
		}

	protected:
		inline static constexpr T EPSILON = T(1e-13);

	public:
		constexpr Matrix &operator=(const Matrix &mat) {
			return (this == &mat) ? Self() : CalcFrom<std::identity{}>(mat, mat);
		}

		template<typename U>
		constexpr operator Sibling<U>() const { return Sibling<U>().template CalcFrom<std::identity{}>(Self(), Self()); }

		constexpr bool operator==(const Matrix &compare) const {
			if constexpr(Matrix::IsDynamic()) {
				if((Self().GetRows() != compare.GetRows()) || (Self().GetCols() != compare.GetCols())) {
					return false;
				}
			}

			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < Self().GetCols(); ++j) {
					if(!IsCloseToZeroEco(Data(i, j) - compare(i, j))) {
						return false;
					}
				}
			}
			return true;
		}

		constexpr bool operator!=(const Matrix &compare) const { return (!(Self() == compare)); }
		constexpr Matrix &operator+=(const T scalar) { return Self() = Self() + scalar; }
		constexpr Matrix &operator-=(const T scalar) { return Self() = Self() - scalar; }
		constexpr Matrix &operator*=(const T scalar) { return Self() = Self() * scalar; }
		constexpr Matrix &operator/=(const T scalar) { return Self() = Self() / scalar; }
		constexpr Matrix &operator+=(const Matrix &mat) { return Self() = Self() + mat; }
		constexpr Matrix &operator-=(const Matrix &mat) { return Self() = Self() - mat; }
		constexpr Matrix &operator*=(const Matrix &mat) { return Self() = Self() * mat; }

		constexpr Matrix operator-() const { return Matrix().template CalcFrom<std::negate<>{}>(Self(), Self()); }
		constexpr Matrix operator+(const T scalar) const { return Matrix().template CalcFrom<std::plus<>{}>(Self(), Self(), scalar); }
		constexpr Matrix operator-(const T scalar) const { return Matrix().template CalcFrom<std::minus<>{}>(Self(), Self(), scalar); }
		constexpr Matrix operator*(const T scalar) const { return Matrix().template CalcFrom<std::multiplies<>{}>(Self(), Self(), scalar); }
		constexpr Matrix operator/(const T scalar) const {
			return (scalar == T(0)) ? CreateInvalid() : Matrix().template CalcFrom<std::divides<>{}>(Self(), Self(), scalar);
		}

		friend constexpr Matrix operator+(const T scalar, const Matrix &mat) { return Matrix().template CalcFrom<std::plus<>{}>(mat, scalar, mat); }
		friend constexpr Matrix operator-(const T scalar, const Matrix &mat) { return Matrix().template CalcFrom<std::minus<>{}>(mat, scalar, mat); }
		friend constexpr Matrix operator*(const T scalar, const Matrix &mat) { return Matrix().template CalcFrom<std::multiplies<>{}>(mat, scalar, mat); }

		constexpr Matrix operator+(const Matrix &mat) const {
			if constexpr(Matrix::IsDynamic()) {
				if(Self().GetRows() != mat.GetRows() || Self().GetCols() != mat.GetCols()) {
					return CreateInvalid();
				}
			}
			return Matrix().template CalcFrom<std::plus<>{}>(Self(), Self(), mat);
		}

		constexpr Matrix operator-(const Matrix &mat) const {
			if constexpr(Matrix::IsDynamic()) {
				if((Self().GetRows() != mat.GetRows()) || (Self().GetCols() != mat.GetCols())) {
					return CreateInvalid();
				}
			}
			return Matrix().template CalcFrom<std::minus<>{}>(Self(), Self(), mat);
		}

		template<typename U, Idx R, Idx C>
		constexpr Sibling<decltype(T() * U()), ROWS, C> operator*(const Sibling<U, R, C> &mat) const
			requires(COLS == 0 || R == 0 || COLS == R)
		{
			if constexpr(Matrix::IsDynamic() || Sibling<U, R, C>::IsDynamic()) {
				if(Self().GetCols() != mat.GetRows()) {
					return CreateInvalid();
				}
			}
			auto res = Sibling<decltype(T() * U()), ROWS, C>::AllocateIfDynamic(Self().GetRows(), mat.GetCols());

			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < mat.GetCols(); ++j) {
					res(i, j) = T(0);
					for(Idx k = 0; k < Self().GetCols(); ++k) {
						res(i, j) += Data(i, k) * mat(k, j);
					}
				}
			}
			return res;
		}

	public:
		/* Debug */
		constexpr void Print() const {
			ForEachElementCall<[](const auto &v) { std::cout << v << ", "; }>(Self());
			std::cout << '\n';
		}

		static constexpr Matrix CreateInvalid() {
			Matrix res;
			res.SetInvalid();
			return res;
		}

		static constexpr Matrix Zero()
			requires(Matrix::IsStatic())
		{
			Matrix res;
			res.SetZero();
			return res;
		}

		static constexpr Matrix Zero(Idx rows, Idx cols)
			requires(Matrix::IsDynamic())
		{
			Matrix res = Matrix::AllocateIfDynamic(rows, cols);
			res.SetZero();
			return res;
		}

		constexpr bool IsValid() const
			requires(Matrix::IsStatic())
		{
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < Self().GetCols(); ++j) {
					if(std::isnan(Data(i, j)))
						return false;
				}
			}
			return true;
		}

		constexpr bool IsValid() const
			requires(Matrix::IsDynamic())
		{
			return ((Self().GetRows() > 0) && (Self().GetRows() <= _::MAX_MATRIX_ALLOCATION_SIZE) &&
					(Self().GetCols() > 0) && (Self().GetCols() <= _::MAX_MATRIX_ALLOCATION_SIZE));
		}

		constexpr void SetInvalid() {
			if constexpr(Matrix::IsDynamic()) {
				Self().ResizeIfDynamic(0, 0);
			} else {
				Data(0, 0) = std::numeric_limits<double>::quiet_NaN();
			}
		}

		constexpr void SetHomogen(T val) {
			ForEachElementAssign<std::identity{}>(val);
		}

		constexpr void SetZero() {
			SetHomogen(T(0));
		}

		constexpr void SetRandom(int32_t minRand, int32_t maxRand) {
			const int32_t scope = (maxRand - minRand + 1);
			ForEachElementAssign<[](auto minRand, auto scope) { return T((rand() % scope) + minRand); }>(minRand, scope);
		}

		constexpr auto Transpose() const {
			auto res = Sibling<T, COLS, ROWS>::AllocateIfDynamic(Self().GetCols(), Self().GetRows());
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < Self().GetCols(); ++j) {
					res(j, i) = Data(i, j);
				}
			}
			return res;
		}

	public:
#ifdef EIGEN_AVAILABLE
		using SiblingEigenMatrix = Eigen::Matrix<T, (ROWS > 0) ? ROWS : Eigen::Dynamic, (COLS > 0) ? COLS : Eigen::Dynamic>;

		constexpr operator SiblingEigenMatrix() const {
			SiblingEigenMatrix eigenMatrix(Self().GetRows(), Self().GetCols());
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < Self().GetCols(); ++j) {
					eigenMatrix(i, j) = Data(i, j);
				}
			}
			return eigenMatrix;
		}
#endif
	};

	/* ================================================================================================== */

	template<class DataClass>
	class AbstractSquareMatrix : public AbstractMatrix<DataClass>
	{
	protected:
		static_assert(_::SiblingTrait<DataClass>::R == _::SiblingTrait<DataClass>::C,
					  "AbstractSquareMatrix must have same static rows and columns number or both danamic!");
		constexpr AbstractSquareMatrix(){};
		constexpr ~AbstractSquareMatrix() = default;

		using Matrix = DataClass;
		using Base = AbstractMatrix<DataClass>;
		using typename Base::Idx, typename Base::T;
		using Base::Self, Base::Data, Base::CreateInvalid;
		using Base::EPSILON;

	private:
		// Container needed for some internal calculations
		class Vector : private Base::template Sibling<T, Base::ROWS, 1>
		{
		private:
			using VectorBase = Base::template Sibling<T, Base::ROWS, 1>;

		public:
			constexpr Vector() = delete;
			constexpr Vector(Idx rows) { VectorBase::ResizeIfDynamic(rows, Idx(1)); }
			constexpr Vector(Idx rows, T initialValue) : Vector(rows) { VectorBase::SetHomogen(initialValue); }
			constexpr ~Vector() = default;

			constexpr T &operator()(Idx index) { return VectorBase::Self()(index, Idx(0)); }
			constexpr const T &operator()(Idx index) const { return VectorBase::Self()(index, Idx(0)); }
		};

	public:
		static consteval Matrix Identity()
			requires(Matrix::IsStatic())
		{
			Matrix res;
			res.SetIdentity();
			return res;
		}

		static constexpr Matrix Identity(Idx rows, Idx cols)
			requires(Matrix::IsDynamic())
		{
			Matrix res = Matrix::AllocateIfDynamic(rows, cols);
			res.SetIdentity();
			return res;
		}

		constexpr void SetDiagonal(const T val) {
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < Self().GetCols(); ++j) {
					Data(i, j) = (i == j) ? val : T(0);
				}
			}
		}

		constexpr void SetIdentity() {
			SetDiagonal(T(1));
		}

		constexpr Matrix Inverse() const {
			const auto &[Q, R] = QrDecomposition();

			// Check if R is invertible (no zeros on the diagonal)
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				if(std::abs(R(i, i)) < EPSILON) {
					return CreateInvalid();
				}
			}

			// Solve Rx = Q^T b for each column of the identity
			Matrix inv;

			for(Idx col = 0; col < Self().GetRows(); ++col) {
				Vector x{Self().GetRows()};

				// Backward substitution
				for(int i = Self().GetRows() - 1; i >= 0; --i) {
					const auto &b_i = Q(col, i);  // b(i) = Q^T(i,col)

					T sum = T(0);
					for(Idx j = i + 1; j < Self().GetRows(); ++j) {
						sum += R(i, j) * x(j);
					}
					x(i) = (b_i - sum) / R(i, i);
				}

				for(Idx i = 0; i < Self().GetRows(); ++i) {
					inv(i, col) = x(i);
				}
			}

			return inv;
		}

		constexpr std::tuple<Matrix, Matrix> QrDecomposition() const {
			Matrix R = Self();
			Matrix Q = Identity();

			for(Idx k = 0; k < Self().GetRows() - 1; ++k) {
				// Create Householder vector 'v' from a part of column 'k'
				Vector v{Self().GetRows()};

				T normSqCol = T(0);
				for(Idx i = k; i < Self().GetRows(); ++i) {
					const T colVal = R(i, k);
					v(i) = colVal;
					normSqCol += colVal * colVal;
				}

				if(normSqCol < EPSILON * EPSILON) {
					continue;  // Column near to zeros
				}

				const T normCol = std::sqrt(normSqCol);
				v(k) += (v(k) >= T(0) ? normCol : -normCol);

				// Normalization of 'v' moved to next loop
				const T halfNormSqV = normCol * std::abs(v(k));
				if(halfNormSqV < T(0.5) * EPSILON * EPSILON) {
					continue;
				}

				// Update matrix R: R = H * R
				for(Idx j = k; j < Self().GetRows(); ++j) {
					T dot = T(0);
					for(Idx i = k; i < Self().GetRows(); ++i) {
						dot += v(i) * R(i, j);
					}
					dot /= halfNormSqV;
					for(Idx i = k; i < Self().GetRows(); ++i) {
						R(i, j) -= v(i) * dot;
					}
				}

				// Update matrix Q: Q = Q * H
				for(Idx i = 0; i < Self().GetRows(); ++i) {
					T dot = T(0);
					for(Idx j = k; j < Self().GetRows(); ++j) {
						dot += v(j) * Q(i, j);
					}
					dot /= halfNormSqV;
					for(Idx j = k; j < Self().GetRows(); ++j) {
						Q(i, j) -= v(j) * dot;
					}
				}
			}

			return std::make_tuple(Q, R);
		}

		constexpr Matrix InverseWithLu() const {
			const auto &[LU, P] = LuDecomposition();
			if(!LU.IsValid()) {
				return CreateInvalid();
			}

			Matrix inv;

			for(Idx bInd = 0; bInd < Self().GetRows(); ++bInd) {
				Vector y{Self().GetRows(), T(0)};
				const Idx pbInd = P(bInd);
				y(bInd) = 1.0;

				// Solve L * y = Pb
				for(Idx i = bInd + 1; i < Self().GetRows(); ++i) {
					for(Idx j = bInd; j < i; ++j) {
						y(i) -= LU(i, j) * y(j);
					}
				}

				auto &x = y;  // Cost reduction instead of vector copying

				// Solve U * x = y
				for(Idx i = Self().GetRows() - 1; i >= 0; --i) {
					for(Idx j = i + 1; j < Self().GetRows(); ++j) {
						x(i) -= LU(i, j) * x(j);
					}

					if(std::abs(LU(i, i)) < EPSILON) {
						return CreateInvalid();
					}
					x(i) /= LU(i, i);

					inv(i, pbInd) = x(i);
				}
			}

			return inv;
		}

		constexpr std::tuple<Matrix, Vector> LuDecomposition() const {
			Matrix LU = Self();
			Vector P{Self().GetRows()};

			for(Idx i = 0; i < Self().GetRows(); ++i) {
				P(i) = i;
			}

			for(Idx k = 0; k < Self().GetRows(); ++k) {
				Idx pivot = k;
				T maxVal = std::abs(LU(k, k));

				for(Idx i = k + 1; i < Self().GetRows(); ++i) {
					const T val = std::abs(LU(i, k));
					if(val > maxVal) {
						maxVal = val;
						pivot = i;
					}
				}

				if(maxVal < EPSILON) {
					return std::make_tuple(CreateInvalid(), P);
				}

				std::swap(P(k), P(pivot));
				for(Idx i = 0; i < Self().GetRows(); ++i) {
					std::swap(LU(k, i), LU(pivot, i));
				}

				for(Idx i = k + 1; i < Self().GetRows(); ++i) {
					const T factor = LU(i, k) / LU(k, k);
					LU(i, k) = factor;
					for(Idx j = k + 1; j < Self().GetRows(); ++j) {
						LU(i, j) -= factor * LU(k, j);
					}
				}
			}

			return std::make_tuple(LU, P);
		}

		constexpr T Determinant() const {
			const auto &[Q, R] = QrDecomposition();

			T detR = T(1);
			T detQ = T(1);

			// Determinant of an orthogonal matrix Q should be either +1 or -1
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				if(Q(i, i) < T(0)) {
					detQ *= -T(1);
				}
				detR *= R(i, i);
			}

			return detQ * detR;
		}
	};

	/* ================================================================================================== */

	template<class DataClass>
	class AbstractVector : public AbstractMatrix<DataClass>
	{
	protected:
		static_assert(_::SiblingTrait<DataClass>::C == 1,
					  "AbstractVector must have only one and static column!");
		constexpr AbstractVector(){};
		constexpr ~AbstractVector() = default;

		using Matrix = DataClass;
		using Base = AbstractMatrix<DataClass>;
		using Base::Self, Base::CreateInvalid, Base::EPSILON, Base::ROWS;
		using typename Base::Idx, typename Base::T;

		constexpr T &Data(Idx index) { return Base::Data(index, 0); }
		constexpr const T &Data(Idx index) const { return Base::Data(index, 0); }

	public:
		constexpr T DotProduct(const Matrix &other) const {
			T res = T(0);
			for(Idx i = 0; i < other.GetRows(); ++i) {
				res += Data(i) * other.Data(i);
			}
			return res;
		}

		constexpr T Norm() const {
			return std::sqrt(DotProduct(Self()));
		}

		constexpr bool Normalize() {
			const T norm = Norm();
			if(norm < EPSILON) {
				return false;
			}

			for(Idx i = 0; i < Self().GetRows(); ++i) {
				Data(i) /= norm;
			}
			return true;
		}

		constexpr const T &x() const
			requires(ROWS >= 1)
		{ return Data(0); }

		constexpr const T &y() const
			requires(ROWS >= 2)
		{ return Data(1); }

		constexpr const T &z() const
			requires(ROWS >= 3)
		{ return Data(2); }
	};

	/* ================================================================================================== */

	namespace _
	{
		template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
		struct SiblingTrait<Matrix<T, ROWS, COLS>> {
			template<typename U, MatrixIdx R = ROWS, MatrixIdx C = COLS>
			using Type = Matrix<U, R, C>;
			using RawType = T;
			inline static constexpr MatrixIdx R = ROWS;
			inline static constexpr MatrixIdx C = COLS;
		};

		template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
			requires(ROWS >= 0 && ROWS <= MAX_MATRIX_ALLOCATION_SIZE && COLS >= 0 && COLS <= MAX_MATRIX_ALLOCATION_SIZE)
		class CompactMatrixAllocator
		{
		private:
			inline static constexpr bool IS_STATIC = ROWS > 0 && COLS > 0;
			inline static constexpr MatrixIdx ROWS_ALLOC = (ROWS > 0) ? ROWS : MAX_MATRIX_ALLOCATION_SIZE;

			using Row = std::array<T, (COLS > 0) ? COLS : MAX_MATRIX_ALLOCATION_SIZE>;
			using MatrixData = std::conditional_t<IS_STATIC, std::array<Row, ROWS_ALLOC>,
												  std::pair<std::array<MatrixIdx, ROWS == COLS ? 2 : 1>, std::array<Row, ROWS_ALLOC>>>;

		public:
			template<typename U>
			constexpr CompactMatrixAllocator(const CompactMatrixAllocator<U, ROWS, COLS> &other) {
				if constexpr(!IS_STATIC)
					matrixData.first = other.matrixData.first;
			}

			constexpr Row &operator[](MatrixIdx row) {
				if constexpr(IS_STATIC)
					return matrixData[row];
				else
					return matrixData.second[row];
			}
			constexpr const Row &operator[](MatrixIdx row) const {
				if constexpr(IS_STATIC)
					return matrixData[row];
				else
					return matrixData.second[row];
			}

			static constexpr MatrixIdx GetRows()
				requires(ROWS > 0)
			{ return ROWS; }
			static constexpr MatrixIdx GetCols()
				requires(COLS > 0)
			{ return COLS; }
			constexpr MatrixIdx GetRows() const
				requires(ROWS == 0)
			{ return matrixData.first[0]; }
			constexpr MatrixIdx GetCols() const
				requires(COLS == 0)
			{ return matrixData.first[matrixData.first.size() - 1]; }

			template<typename U>
			constexpr void ResizeIfDynamic(const CompactMatrixAllocator<U, ROWS, COLS> &other) {
				if constexpr(!IS_STATIC)
					matrixData.first = other.matrixData.first;
			}
			constexpr void ResizeIfDynamic(MatrixIdx rows, MatrixIdx cols) {
				if constexpr(ROWS == 0)
					matrixData.first[0] = rows;
				if constexpr(COLS == 0)
					matrixData.first[matrixData.first.size() - 1] = cols;
			}

		protected:
			constexpr CompactMatrixAllocator() {
				if constexpr(!IS_STATIC)
					matrixData.first = decltype(matrixData.first){};
			}

		private:
			MatrixData matrixData;
		};
	}

	/* ================================================================================================== */

	template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
	class Matrix : public _::CompactMatrixAllocator<T, ROWS, COLS>,
				   public std::conditional_t<ROWS == COLS, AbstractSquareMatrix<Matrix<T, ROWS, COLS>>,
											 std::conditional_t<COLS == 1, AbstractVector<Matrix<T, ROWS, COLS>>,
																AbstractMatrix<Matrix<T, ROWS, COLS>>>>
	{
	public:
		using Allocator = _::CompactMatrixAllocator<T, ROWS, COLS>;
		using Base = std::conditional_t<ROWS == COLS, AbstractSquareMatrix<Matrix<T, ROWS, COLS>>,
										std::conditional_t<COLS == 1, AbstractVector<Matrix<T, ROWS, COLS>>,
														   AbstractMatrix<Matrix<T, ROWS, COLS>>>>;

		using Allocator::GetRows, Allocator::GetCols, Allocator::ResizeIfDynamic;
		using Allocator::operator[];
		using Base::SetHomogen, Base::CalcFrom, Base::ForEachElementAssign;
		using typename Base::Idx;

		static consteval bool IsRowDynamic() { return ROWS == 0; }
		static consteval bool IsColDynamic() { return COLS == 0; }
		static consteval bool IsRowStatic() { return !IsRowDynamic(); }
		static consteval bool IsColStatic() { return !IsColDynamic(); }
		static consteval bool IsStatic() { return IsRowStatic() && IsColStatic(); }
		static consteval bool IsDynamic() { return !IsStatic(); }
		static consteval bool IsSquareAllocated() { return ROWS == COLS; }

		static consteval bool IsSquare()
			requires(IsStatic())
		{ return (GetRows() == GetCols()); }
		constexpr bool IsSquare() const
			requires(IsDynamic())
		{ return (GetRows() == GetCols()); }

		static constexpr Matrix AllocateIfDynamic(Idx rows, Idx cols) {
			Matrix res;
			res.ResizeIfDynamic(rows, cols);
			return res;
		}
		template<typename U, Idx R, Idx C>
		static constexpr Matrix AllocateIfDynamic(const Matrix<U, R, C> &other) {
			Matrix res;
			res.ResizeIfDynamic(other);
			return res;
		}

		constexpr T &At(Idx row, Idx col) { return (*this)(row, col); }
		constexpr const T &At(Idx row, Idx col) const { return (*this)(row, col); }

		constexpr T &operator()(Idx row, Idx col) { return (*this)[row][col]; }
		constexpr const T &operator()(Idx row, Idx col) const { return (*this)[row][col]; }

		constexpr T &operator()(Idx index)
			requires(COLS == 1)
		{ return (*this)(index, 0); }
		constexpr const T &operator()(Idx index) const
			requires(COLS == 1)
		{ return (*this)(index, 0); }

	public:
		constexpr Matrix() {}

		constexpr Matrix(Idx rows, Idx cols)
			requires(IsDynamic())
		{ ResizeIfDynamic(rows, cols); }

		constexpr Matrix(Idx rows, Idx cols, T initialValue)
			requires(IsDynamic())
			: Matrix(rows, cols) { SetHomogen(initialValue); }

		constexpr Matrix(T initialValue)
			requires(IsStatic())
		{ SetHomogen(initialValue); }

		template<typename U>
		constexpr Matrix(const Matrix<U, ROWS, COLS> &other) {
			CalcFrom<std::identity{}>(other, other);
		}

#ifdef EIGEN_AVAILABLE
		using SiblingEigenMatrix = Eigen::Matrix<T, (ROWS > 0) ? ROWS : Eigen::Dynamic, (COLS > 0) ? COLS : Eigen::Dynamic>;

		constexpr Matrix(const SiblingEigenMatrix &eigenMatrix) {
			ResizeIfDynamic(eigenMatrix.rows(), eigenMatrix.cols());
			ForEachElementAssign<std::identity{}>(eigenMatrix);
		}
#endif
	};

	/* ================================================================================================== */

	using Matrix3f = Matrix<float, 3, 3>;
	using Matrix4f = Matrix<float, 4, 4>;
	using MatrixXf = Matrix<float, 0, 0>;
	using Matrix3d = Matrix<double, 3, 3>;
	using Matrix4d = Matrix<double, 4, 4>;
	using MatrixXd = Matrix<double, 0, 0>;

	template<typename T>
	using Matrix3 = Matrix<T, 3, 3>;

	template<typename T>
	using Matrix4 = Matrix<T, 4, 4>;

	template<typename T>
	using MatrixX = Matrix<T, 0, 0>;

	/* ================================================================================================== */

	using Vector3f = Matrix<float, 3, 1>;
	using Vector4f = Matrix<float, 4, 1>;
	using VectorXf = Matrix<float, 0, 1>;
	using Vector3d = Matrix<double, 3, 1>;
	using Vector4d = Matrix<double, 4, 1>;
	using VectorXd = Matrix<double, 0, 1>;

	template<typename T>
	using Vector3 = Matrix<T, 3, 1>;

	template<typename T>
	using Vector4 = Matrix<T, 4, 1>;

	template<typename T>
	using VectorX = Matrix<T, 0, 1>;

	template<typename T, MatrixIdx LENGTH>
	using Vector = Matrix<T, LENGTH, 1>;
}
