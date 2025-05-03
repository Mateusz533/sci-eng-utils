#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
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
		constexpr T &Data(MatrixIdx row, MatrixIdx col) { return Self()(row, col); }
		constexpr const T &Data(MatrixIdx row, MatrixIdx col) const { return Self()(row, col); }

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
		constexpr Matrix &operator=(const Matrix &obj) {
			if(this == &obj) { return Self(); }
			Self().ResizeIfDynamic(obj);

			const T *src = &obj.Data(0, 0);
			T *dst = &Data(0, 0);

			if constexpr(Matrix::IsColStatic()) {
				std::memcpy(dst, src, sizeof(T) * static_cast<size_t>(Self().GetRows() * Self().GetCols()));
			} else {
				for(Idx i = 0; i < Self().GetRows(); ++i) {
					std::memcpy(dst, src, sizeof(T) * static_cast<size_t>(Self().GetCols()));
					src += _::MAX_MATRIX_ALLOCATION_SIZE;
					dst += _::MAX_MATRIX_ALLOCATION_SIZE;
				}
			}

			return Self();
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
		constexpr Matrix operator*(const Sibling<U, R, C> &mat) const
			requires(!(Matrix::IsStatic() && Sibling<U, R, C>::IsStatic() && ROWS != C))
		{
			if constexpr(Matrix::IsDynamic() || Sibling<U, R, C>::IsDynamic()) {
				if(Self().GetCols() != mat.GetRows()) {
					return CreateInvalid();
				}
			}
			Matrix res = Matrix::AllocateIfDynamic(Self().GetRows(), mat.GetCols());

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

		static consteval Matrix CreateInvalid() {
			Matrix res;
			res.SetInvalid();
			return res;
		}

		static consteval Matrix Zero()
			requires(Matrix::IsStatic())
		{
			Matrix res;
			res.SetToZero();
			return res;
		}

		static constexpr Matrix Zero(Idx rows, Idx cols)
			requires(Matrix::IsDynamic())
		{
			Matrix res = Matrix::AllocateIfDynamic(rows, cols);
			res.SetToZero();
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
				SetHomogen(std::numeric_limits<double>::quiet_NaN());
			}
		}

		constexpr void SetHomogen(const T val) {
			ForEachElementAssign<std::identity{}>(val);
		}

		constexpr void SetToZero() {
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

		constexpr void SetDiag(const T val) {
			for(Idx i = 0; i < Self().GetRows(); ++i) {
				for(Idx j = 0; j < Self().GetCols(); ++j) {
					Data(i, j) = (i == j) ? val : T(0);
				}
			}
		}

		constexpr void SetIdentity() {
			SetDiag(T(1));
		}

		constexpr void Inverse() {
			Self() = GetInverse();
		}

		constexpr Matrix GetInverse() const {}	// TODO: Implement QR algorithm
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
		using Base::Self, Base::Data, Base::CreateInvalid, Base::EPSILON;
		using typename Base::Idx, typename Base::T;

		constexpr T &Data(Idx index) { return Data(index, 0); }
		constexpr const T &Data(Idx index) const { return Data(index, 0); }

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
			using Data = std::conditional_t<IS_STATIC, std::array<Row, ROWS_ALLOC>,
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
			Data matrixData;
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
}
