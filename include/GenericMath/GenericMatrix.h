#pragma once
//
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>
//
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
//
#define EIGEN_AVAILABLE

namespace GenericMath
{
	using MatrixIdx = int16_t;

	namespace _
	{
		static inline constexpr MatrixIdx DYNAMIC_SIZE = 0;
		static inline constexpr MatrixIdx MAX_MATRIX_ALLOCATION_SIZE = 6;

		template<typename>
		struct SiblingTrait;

		template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
			requires((0 <= ROWS && ROWS <= MAX_MATRIX_ALLOCATION_SIZE) && (0 <= COLS && COLS <= MAX_MATRIX_ALLOCATION_SIZE))
		class CompactMatrixAllocator;

		template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
			requires(ROWS >= 0 && COLS >= 0)
		class StdMatrixAllocator;
	}

	template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
	using DefaultMatrixAllocator = _::StdMatrixAllocator<T, ROWS, COLS>;

	template<class DataClass>
	class AbstractMatrix;

	template<class DataClass>
	class AbstractSquareMatrix;

	template<class DataClass>
	class AbstractVector;

	template<typename T, MatrixIdx ROWS, MatrixIdx COLS, template<typename, MatrixIdx, MatrixIdx> typename Allocator = DefaultMatrixAllocator>
	class Matrix;

	/* ================================================================================================== */

	template<class DataClass>
	class AbstractMatrix
	{
	protected:
		using Matrix = DataClass;
		using Idx = MatrixIdx;

		using T = typename _::SiblingTrait<DataClass>::RawType;
		static inline constexpr Idx ROWS = _::SiblingTrait<DataClass>::R;
		static inline constexpr Idx COLS = _::SiblingTrait<DataClass>::C;

		template<typename U, Idx R = ROWS, Idx C = COLS>
		using Sibling = typename _::SiblingTrait<DataClass>::template Type<U, R, C>;

		static inline constexpr T EPSILON = T(1e-13);

	protected:
		constexpr AbstractMatrix() {
			static_assert(std::is_base_of_v<AbstractMatrix, DataClass>,
						  "AbstractMatrix template parameter must derive from AbstractMatrix");
		}
		constexpr ~AbstractMatrix() = default;

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
		constexpr bool HasDifferentSize(const Matrix &compare) const {
			if constexpr(Matrix::IsRowDynamic()) {
				if(Self().GetRows() != compare.GetRows()) {
					return true;
				}
			}
			if constexpr(Matrix::IsColDynamic()) {
				if(Self().GetCols() != compare.GetCols()) {
					return true;
				}
			}
			return false;
		}

	public:
		constexpr Matrix &operator=(const Matrix &mat) {
			return (this == &mat) ? Self() : CalcFrom<std::identity{}>(mat, mat);
		}

		template<typename U>
		constexpr operator Sibling<U>() const { return Sibling<U>().template CalcFrom<std::identity{}>(Self(), Self()); }

		constexpr bool operator==(const Matrix &compare) const {
			if(HasDifferentSize(compare)) {
				return false;
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

		constexpr Matrix operator-() const { return std::move(Matrix().template CalcFrom<std::negate<>{}>(Self(), Self())); }
		constexpr Matrix operator+(const T scalar) const { return std::move(Matrix().template CalcFrom<std::plus<>{}>(Self(), Self(), scalar)); }
		constexpr Matrix operator-(const T scalar) const { return std::move(Matrix().template CalcFrom<std::minus<>{}>(Self(), Self(), scalar)); }
		constexpr Matrix operator*(const T scalar) const { return std::move(Matrix().template CalcFrom<std::multiplies<>{}>(Self(), Self(), scalar)); }
		constexpr Matrix operator/(const T scalar) const {
			return (scalar == T(0)) ? CreateInvalid() : std::move(Matrix().template CalcFrom<std::divides<>{}>(Self(), Self(), scalar));
		}

		friend constexpr Matrix operator+(const T scalar, const Matrix &mat) { return std::move(Matrix().template CalcFrom<std::plus<>{}>(mat, scalar, mat)); }
		friend constexpr Matrix operator-(const T scalar, const Matrix &mat) { return std::move(Matrix().template CalcFrom<std::minus<>{}>(mat, scalar, mat)); }
		friend constexpr Matrix operator*(const T scalar, const Matrix &mat) { return std::move(Matrix().template CalcFrom<std::multiplies<>{}>(mat, scalar, mat)); }

		constexpr Matrix operator+(const Matrix &mat) const {
			if(HasDifferentSize(mat)) {
				return CreateInvalid();
			}
			return std::move(Matrix().template CalcFrom<std::plus<>{}>(Self(), Self(), mat));
		}

		constexpr Matrix operator-(const Matrix &mat) const {
			if(HasDifferentSize(mat)) {
				return CreateInvalid();
			}
			return std::move(Matrix().template CalcFrom<std::minus<>{}>(Self(), Self(), mat));
		}

		template<typename U, Idx R, Idx C>
		constexpr Sibling<decltype(T() * U()), ROWS, C> operator*(const Sibling<U, R, C> &mat) const
			requires(COLS == R || Matrix::IsColDynamic() || Sibling<U, R, C>::IsRowDynamic())
		{
			using ResultType = Sibling<decltype(T() * U()), ROWS, C>;

			if constexpr(Matrix::IsColDynamic() || Sibling<U, R, C>::IsRowDynamic()) {
				if(Self().GetCols() != mat.GetRows()) {
					return AbstractMatrix<ResultType>::CreateInvalid();
				}
			}
			ResultType res = ResultType::AllocateIfDynamic(Self().GetRows(), mat.GetCols());

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
			return ((0 < Self().GetRows() && Self().GetRows() <= Self().GetRowsLimit()) &&
					(0 < Self().GetCols() && Self().GetCols() <= Self().GetColsLimit()));
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
	};

	/* ================================================================================================== */

	template<class DataClass>
	class AbstractSquareMatrix : public AbstractMatrix<DataClass>
	{
		static_assert(_::SiblingTrait<DataClass>::R == _::SiblingTrait<DataClass>::C,
					  "AbstractSquareMatrix must have same static rows and columns number or both danamic!");

		using Matrix = DataClass;
		using Base = AbstractMatrix<DataClass>;

	protected:
		using typename Base::Idx, typename Base::T;
		using Base::Self, Base::Data, Base::CreateInvalid;
		using Base::EPSILON;

		constexpr AbstractSquareMatrix(){};
		constexpr ~AbstractSquareMatrix() = default;

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

		constexpr void SetDiagonal(const Base::template Sibling<T, Base::ROWS, 1> &vec)
			requires(Matrix::IsStatic())
		{
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < Self().GetCols(); ++j) {
					Data(i, j) = (i == j) ? vec(i) : T(0);
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
				for(Idx i = Self().GetRows() - 1; i >= 0; --i) {
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
				v(k) += (v(k) >= T(0)) ? normCol : -normCol;

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

		constexpr T Trace() const {
			T trace = T(0);

			for(Idx i = 0; i < Self().GetRows(); ++i) {
				trace += Data(i, i);
			}

			return trace;
		}
	};

	/* ================================================================================================== */

	template<class DataClass>
	class AbstractVector : public AbstractMatrix<DataClass>
	{
		static_assert(_::SiblingTrait<DataClass>::C == 1,
					  "AbstractVector must have only one and static column!");

		using Matrix = DataClass;
		using Base = AbstractMatrix<DataClass>;

	protected:
		using typename Base::Idx, typename Base::T;
		using Base::Self, Base::CreateInvalid;
		using Base::EPSILON, Base::ROWS;

		constexpr AbstractVector(){};
		constexpr ~AbstractVector() = default;

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
			requires((1 <= ROWS && ROWS <= 3))
		{ return Data(0); }

		constexpr const T &y() const
			requires((2 <= ROWS && ROWS <= 3))
		{ return Data(1); }

		constexpr const T &z() const
			requires(ROWS == 3)
		{ return Data(2); }

		template<Idx OFFSET, typename... Args>
			requires(OFFSET >= 0 && sizeof...(Args) == ROWS - OFFSET - 1)
		constexpr void Fill(T v, Args... args) {
			Data(OFFSET) = v;
			if constexpr(sizeof...(Args) > 0)
				Fill<OFFSET + 1>(args...);
		}
	};

	/* ================================================================================================== */

	namespace _
	{
		template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
		struct SiblingTrait<Matrix<T, ROWS, COLS>> {
			template<typename U, MatrixIdx R = ROWS, MatrixIdx C = COLS>
			using Type = Matrix<U, R, C>;
			using RawType = T;
			static inline constexpr MatrixIdx R = ROWS;
			static inline constexpr MatrixIdx C = COLS;
		};

		template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
			requires((0 <= ROWS && ROWS <= MAX_MATRIX_ALLOCATION_SIZE) && (0 <= COLS && COLS <= MAX_MATRIX_ALLOCATION_SIZE))
		class CompactMatrixAllocator
		{
		private:
			static inline constexpr bool IS_DYNAMIC = (ROWS == DYNAMIC_SIZE || COLS == DYNAMIC_SIZE);
			static inline constexpr MatrixIdx ROWS_ALLOC = (ROWS == DYNAMIC_SIZE) ? MAX_MATRIX_ALLOCATION_SIZE : ROWS;
			static inline constexpr MatrixIdx COLS_ALLOC = (COLS == DYNAMIC_SIZE) ? MAX_MATRIX_ALLOCATION_SIZE : COLS;

			using Row = std::array<T, COLS_ALLOC>;
			using MatrixData = std::conditional_t<!IS_DYNAMIC, std::array<Row, ROWS_ALLOC>,
												  std::pair<std::array<MatrixIdx, (ROWS == COLS) ? 2 : 1>, std::array<Row, ROWS_ALLOC>>>;

		public:
			template<typename U>
			constexpr CompactMatrixAllocator(const CompactMatrixAllocator<U, ROWS, COLS> &other) {
				ResizeIfDynamic(other);
			}

			constexpr Row &operator[](MatrixIdx row) {
				if constexpr(IS_DYNAMIC)
					return matrixData.second[row];
				else
					return matrixData[row];
			}
			constexpr const Row &operator[](MatrixIdx row) const {
				if constexpr(IS_DYNAMIC)
					return matrixData.second[row];
				else
					return matrixData[row];
			}

			static constexpr MatrixIdx GetRows()
				requires(ROWS != DYNAMIC_SIZE)
			{ return ROWS; }
			static constexpr MatrixIdx GetCols()
				requires(COLS != DYNAMIC_SIZE)
			{ return COLS; }
			constexpr MatrixIdx GetRows() const
				requires(ROWS == DYNAMIC_SIZE)
			{ return DynamicSize()[0]; }
			constexpr MatrixIdx GetCols() const
				requires(COLS == DYNAMIC_SIZE)
			{ return DynamicSize()[ROWS == DYNAMIC_SIZE ? 1 : 0]; }

			static constexpr MatrixIdx GetRowsLimit() { return ROWS_ALLOC; }
			static constexpr MatrixIdx GetColsLimit() { return COLS_ALLOC; }

			template<typename U>
			constexpr void ResizeIfDynamic(const CompactMatrixAllocator<U, ROWS, COLS> &other) {
				if constexpr(IS_DYNAMIC)
					DynamicSize() = other.DynamicSize();
			}
			constexpr void ResizeIfDynamic(MatrixIdx rows, MatrixIdx cols) {
				if constexpr(ROWS == DYNAMIC_SIZE)
					DynamicSize()[0] = rows;
				if constexpr(COLS == DYNAMIC_SIZE)
					DynamicSize()[ROWS == DYNAMIC_SIZE ? 1 : 0] = cols;
			}

		protected:
			constexpr CompactMatrixAllocator() {
				if constexpr(IS_DYNAMIC)
					DynamicSize() = std::remove_reference_t<decltype(DynamicSize())>{};
			}

			constexpr auto &DynamicSize()
				requires(IS_DYNAMIC)
			{ return matrixData.first; }
			constexpr const auto &DynamicSize() const
				requires(IS_DYNAMIC)
			{ return matrixData.first; }

		private:
			MatrixData matrixData;
		};

		namespace
		{
			static_assert(sizeof(CompactMatrixAllocator<double, 3, 5>) == 3 * 5 * sizeof(double));
			static_assert(sizeof(CompactMatrixAllocator<float, 0, 0>) == 2 * sizeof(MatrixIdx) + (MAX_MATRIX_ALLOCATION_SIZE * MAX_MATRIX_ALLOCATION_SIZE) * sizeof(float));
			static_assert(sizeof(CompactMatrixAllocator<char, 0, 4>) == sizeof(MatrixIdx) + (MAX_MATRIX_ALLOCATION_SIZE * 4) * sizeof(char));
			static_assert(sizeof(CompactMatrixAllocator<short, 2, 0>) == sizeof(MatrixIdx) + (2 * MAX_MATRIX_ALLOCATION_SIZE) * sizeof(short));
		}

		template<typename T, MatrixIdx ROWS, MatrixIdx COLS>
			requires(ROWS >= 0 && COLS >= 0)
		class StdMatrixAllocator
		{
		private:
			static inline constexpr bool IS_DYNAMIC = (ROWS == DYNAMIC_SIZE || COLS == DYNAMIC_SIZE);

			using Row = std::conditional_t<(COLS == DYNAMIC_SIZE), T[], std::array<T, COLS>>;
			using MatrixData = std::conditional_t<!IS_DYNAMIC, std::array<Row, ROWS>,
												  std::pair<std::array<MatrixIdx, (ROWS == COLS) ? 2 : 1>, std::unique_ptr<T[]>>>;

		public:
			template<typename U>
			constexpr StdMatrixAllocator(const StdMatrixAllocator<U, ROWS, COLS> &other) {
				ResizeIfDynamic(other);
			}

			constexpr Row &operator[](MatrixIdx row) {
				if constexpr(IS_DYNAMIC)
					return (T(&)[])matrixData.second[row * GetCols()];
				else
					return matrixData[row];
			}
			constexpr const Row &operator[](MatrixIdx row) const {
				if constexpr(IS_DYNAMIC)
					return (const T(&)[])matrixData.second[row * GetCols()];
				else
					return matrixData[row];
			}

			static constexpr MatrixIdx GetRows()
				requires(ROWS != DYNAMIC_SIZE)
			{ return ROWS; }
			static constexpr MatrixIdx GetCols()
				requires(COLS != DYNAMIC_SIZE)
			{ return COLS; }
			constexpr MatrixIdx GetRows() const
				requires(ROWS == DYNAMIC_SIZE)
			{ return DynamicSize()[0]; }
			constexpr MatrixIdx GetCols() const
				requires(COLS == DYNAMIC_SIZE)
			{ return DynamicSize()[ROWS == DYNAMIC_SIZE ? 1 : 0]; }

			static constexpr MatrixIdx GetRowsLimit() { return (ROWS == DYNAMIC_SIZE) ? std::numeric_limits<MatrixIdx>::max() : ROWS; }
			static constexpr MatrixIdx GetColsLimit() { return (COLS == DYNAMIC_SIZE) ? std::numeric_limits<MatrixIdx>::max() : COLS; }

			template<typename U>
			constexpr void ResizeIfDynamic(const StdMatrixAllocator<U, ROWS, COLS> &other) {
				if constexpr(IS_DYNAMIC) {
					DynamicSize() = other.DynamicSize();
					matrixData.second = std::make_unique<T[]>(GetRows() * GetCols());
				}
			}
			constexpr void ResizeIfDynamic(MatrixIdx rows, MatrixIdx cols) {
				if constexpr(ROWS == DYNAMIC_SIZE)
					DynamicSize()[0] = rows;
				if constexpr(COLS == DYNAMIC_SIZE)
					DynamicSize()[ROWS == DYNAMIC_SIZE ? 1 : 0] = cols;
				if constexpr(IS_DYNAMIC)
					matrixData.second = std::make_unique<T[]>(GetRows() * GetCols());
			}

		protected:
			constexpr StdMatrixAllocator() {
				if constexpr(IS_DYNAMIC)
					DynamicSize() = std::remove_reference_t<decltype(DynamicSize())>{};
			}

			constexpr auto &DynamicSize()
				requires(IS_DYNAMIC)
			{ return matrixData.first; }
			constexpr const auto &DynamicSize() const
				requires(IS_DYNAMIC)
			{ return matrixData.first; }

		private:
			MatrixData matrixData;
		};

		namespace
		{
			static_assert(sizeof(StdMatrixAllocator<double, 3, 5>) == 3 * 5 * sizeof(double));
			static_assert(sizeof(StdMatrixAllocator<float, 0, 0>) == (1 + 1) * sizeof(void *));
			static_assert(sizeof(StdMatrixAllocator<char, 0, 4>) == (1 + 1) * sizeof(void *));
			static_assert(sizeof(StdMatrixAllocator<short, 2, 0>) == (1 + 1) * sizeof(void *));
		}

		template<class DataClass>
		using AbstractUnifiedMatrix = std::conditional_t<SiblingTrait<DataClass>::R == SiblingTrait<DataClass>::C, AbstractSquareMatrix<DataClass>,
														 std::conditional_t<SiblingTrait<DataClass>::C == 1, AbstractVector<DataClass>,
																			AbstractMatrix<DataClass>>>;
	}

	/* ================================================================================================== */

	template<typename T, MatrixIdx ROWS, MatrixIdx COLS, template<typename, MatrixIdx, MatrixIdx> typename AllocatorT>
	class Matrix : public AllocatorT<T, ROWS, COLS>,
				   public _::AbstractUnifiedMatrix<Matrix<T, ROWS, COLS>>
	{
	public:
		using Allocator = AllocatorT<T, ROWS, COLS>;
		using Base = _::AbstractUnifiedMatrix<Matrix<T, ROWS, COLS>>;

		using Allocator::GetRows, Allocator::GetCols, Allocator::GetRowsLimit, Allocator::GetColsLimit, Allocator::ResizeIfDynamic;
		using Allocator::operator[];
		using Base::SetHomogen, Base::CalcFrom, Base::ForEachElementAssign;
		using typename Base::Idx;

		static consteval bool IsRowDynamic() { return ROWS == _::DYNAMIC_SIZE; }
		static consteval bool IsColDynamic() { return COLS == _::DYNAMIC_SIZE; }
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

		constexpr Matrix(const std::array<T, ROWS> &data)
			requires(IsStatic() && COLS == 1)
		{
			for(Idx i = 0; i < ROWS; ++i)
				Base::Data(i) = data[i];
		}

		template<typename... Args>
		constexpr Matrix(const Args... args)
			requires(IsStatic() && COLS == 1 && (std::is_same_v<Args, T> && ...) && sizeof...(Args) == ROWS)
		{ Base::template Fill<0>(args...); }

		constexpr Matrix(const std::array<std::array<T, COLS>, ROWS> &data)
			requires(IsStatic() && COLS > 1)
		{
			for(Idx i = 0; i < ROWS; ++i)
				for(Idx j = 0; j < COLS; ++j)
					Base::Data(i, j) = data[i][j];
		}

		template<typename U>
		constexpr Matrix(const Matrix<U, ROWS, COLS> &other) {
			CalcFrom<std::identity{}>(other, other);
		}

#ifdef EIGEN_AVAILABLE
		using SiblingEigenMatrix = Eigen::Matrix<T, IsRowDynamic() ? Eigen::Dynamic : ROWS, IsColDynamic() ? Eigen::Dynamic : COLS>;

		constexpr Matrix(const SiblingEigenMatrix &eigenMatrix) {
			ResizeIfDynamic(eigenMatrix.rows(), eigenMatrix.cols());
			ForEachElementAssign<std::identity{}>(eigenMatrix);
		}

		constexpr operator SiblingEigenMatrix() const {
			SiblingEigenMatrix eigenMatrix(GetRows(), GetCols());
			for(Idx i = 0; i < GetRows(); ++i) {
				for(Idx j = 0; j < GetCols(); ++j) {
					eigenMatrix(i, j) = Base::Data(i, j);
				}
			}
			return eigenMatrix;
		}
#endif
	};

	/* ================================================================================================== */

	template<typename T>
	using MatrixX = Matrix<T, _::DYNAMIC_SIZE, _::DYNAMIC_SIZE>;

	template<typename T>
	using Matrix3 = Matrix<T, 3, 3>;

	template<typename T>
	using Matrix4 = Matrix<T, 4, 4>;

	using MatrixXf = MatrixX<float>;
	using Matrix3f = Matrix3<float>;
	using Matrix4f = Matrix4<float>;
	using MatrixXd = MatrixX<double>;
	using Matrix3d = Matrix3<double>;
	using Matrix4d = Matrix4<double>;

	/* ================================================================================================== */

	template<typename T, MatrixIdx LENGTH>
	using Vector = Matrix<T, LENGTH, 1>;

	template<typename T>
	using VectorX = Vector<T, _::DYNAMIC_SIZE>;

	template<typename T>
	using Vector3 = Vector<T, 3>;

	template<typename T>
	using Vector4 = Vector<T, 4>;

	using VectorXf = VectorX<float>;
	using Vector3f = Vector3<float>;
	using Vector4f = Vector4<float>;
	using VectorXd = VectorX<double>;
	using Vector3d = Vector3<double>;
	using Vector4d = Vector4<double>;
}