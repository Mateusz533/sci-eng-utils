#pragma once
//
#include <array>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

namespace GenericMath
{
	using TensorIdx = int32_t;

	namespace Private
	{
		template<typename Data, TensorIdx DIM, TensorIdx RANK>
			requires(RANK >= 0)
		struct MultiArrayComposer {
			using ArrayType = std::array<typename MultiArrayComposer<Data, DIM, RANK - 1>::ArrayType, DIM>;
		};

		template<typename Data, TensorIdx DIM>
		struct MultiArrayComposer<Data, DIM, 0> {
			using ArrayType = Data;
		};
	}

	template<typename Data, TensorIdx RANK, TensorIdx DIM>
		requires(std::is_floating_point_v<Data> && DIM > 0 && RANK >= 0)
	class Tensor
	{
		template<TensorIdx _RANK>
		using ArrayType = typename Private::MultiArrayComposer<Data, DIM, _RANK>::ArrayType;

	public:
		template<typename... Indices>
			requires(sizeof...(Indices) == RANK && (std::is_same_v<Indices, TensorIdx> && ...))
		constexpr Data &operator()(Indices... indices) noexcept {
			return At(std::array<TensorIdx, RANK>{indices...});
		}

		template<typename... Indices>
			requires(sizeof...(Indices) == RANK && (std::is_same_v<Indices, TensorIdx> && ...))
		constexpr const Data &operator()(Indices... indices) const noexcept {
			return At(std::array<TensorIdx, RANK>{indices...});
		}

		constexpr operator Data() const noexcept
			requires(RANK == 0)
		{
			return mData;
		}

		constexpr Data &At(std::span<const TensorIdx, RANK> indices) noexcept {
			return At<RANK>(mData, indices);
		}
		constexpr const Data &At(std::span<const TensorIdx, RANK> indices) const noexcept {
			return At<RANK>(mData, indices);
		}

		template<TensorIdx _RANK>
		[[nodiscard]] constexpr Tensor<Data, RANK + _RANK, DIM> OuterProduct(const Tensor<Data, _RANK, DIM> &other) const noexcept {
			Tensor<Data, RANK + _RANK, DIM> result{};

			ForEachIndex<RANK + _RANK>([&](const auto &resultIndices) {
				const auto rIdxView = std::span<const TensorIdx, RANK + _RANK>(resultIndices);

				const Data a = this->At(rIdxView.template subspan<0, RANK>());
				const Data b = other.At(rIdxView.template subspan<RANK, _RANK>());
				result.At(resultIndices) = a * b;
			});

			return result;
		}

		template<TensorIdx I, TensorIdx J, TensorIdx _RANK>
			requires((0 <= I && I < RANK && 0 <= J && J < _RANK))
		[[nodiscard]] constexpr Tensor<Data, RANK + _RANK - 2, DIM> InnerProduct(const Tensor<Data, _RANK, DIM> &other) const noexcept {
			return InnerProduct(other, I, J);
		}

		template<TensorIdx _RANK>
			requires(RANK >= 1 && _RANK >= 1)
		[[nodiscard]] constexpr Tensor<Data, RANK + _RANK - 2, DIM>
			InnerProduct(const Tensor<Data, _RANK, DIM> &other, TensorIdx i, TensorIdx j) const noexcept {
			if(!(0 <= i && i < RANK && 0 <= j && j < _RANK)) {
				return Tensor<Data, RANK + _RANK - 2, DIM>{};
			}

			// mapping from (R+S-2) result indices to (R) and (S) indices excluding contracted positions
			std::array<TensorIdx, RANK> lhsIndicesMap;
			std::array<TensorIdx, _RANK> rhsIndicesMap;

			TensorIdx currentIndex = 0;
			for(TensorIdx p = 0; p < RANK; ++p) {
				if(p == i) {
					lhsIndicesMap[p] = TensorIdx(0);
					continue;
				}
				lhsIndicesMap[p] = currentIndex++;
			}
			for(TensorIdx p = 0; p < _RANK; ++p) {
				if(p == j) {
					rhsIndicesMap[p] = TensorIdx(0);
					continue;
				}
				rhsIndicesMap[p] = currentIndex++;
			}

			Tensor<Data, RANK + _RANK - 2, DIM> result;

			ForEachIndex<RANK + _RANK - 2>([&](const auto &resultIndices) {
				std::array<TensorIdx, RANK> lhsIndices;
				if constexpr(RANK + _RANK - 2 > 0) {
					for(TensorIdx p = 0; p < RANK; ++p)
						lhsIndices[p] = resultIndices[lhsIndicesMap[p]];
				}

				std::array<TensorIdx, _RANK> rhsIndices;
				if constexpr(RANK + _RANK - 2 > 0) {
					for(TensorIdx p = 0; p < _RANK; ++p)
						rhsIndices[p] = resultIndices[rhsIndicesMap[p]];
				}

				Data diagonalSum{0};
				for(TensorIdx k = 0; k < DIM; ++k) {
					lhsIndices[i] = k;
					rhsIndices[j] = k;
					diagonalSum += this->At(lhsIndices) * other.At(rhsIndices);
				}
				result.At(resultIndices) = diagonalSum;
			});

			return result;
		}

		template<TensorIdx I, TensorIdx J>
		[[nodiscard]] constexpr auto Contraction() const noexcept
			requires((0 <= I && I < J && J < RANK))
		{
			return Contraction(I, J);
		}

		[[nodiscard]] constexpr auto Contraction(TensorIdx i, TensorIdx j) const noexcept
			requires(RANK >= 2)
		{
			// ensure i < j for simpler mapping
			if(j < i) std::swap(i, j);

			if(!(0 <= i && i < j && j < RANK))
				return Tensor<Data, RANK - 2, DIM>{};

			// map positions excluding i and j into resultIndices
			std::array<TensorIdx, RANK> indicesMap;

			TensorIdx currentIndex = 0;
			for(TensorIdx p = 0; p < RANK; ++p) {
				if(p == i || p == j) {
					indicesMap[p] = TensorIdx(0);
					continue;
				}
				indicesMap[p] = currentIndex++;
			}

			Tensor<Data, RANK - 2, DIM> result;

			ForEachIndex<RANK - 2>([&](const auto &resultIndices) {
				std::array<TensorIdx, RANK> sourceIndices;
				if constexpr(RANK - 2 > 0) {
					for(TensorIdx p = 0; p < RANK; ++p) {
						sourceIndices[p] = resultIndices[indicesMap[p]];
					}
				}

				Data diagonalSum = Data{0};
				for(TensorIdx k = 0; k < DIM; ++k) {
					sourceIndices[i] = k;
					sourceIndices[j] = k;
					diagonalSum += this->At(sourceIndices);
				}
				result.At(resultIndices) = diagonalSum;
			});

			return result;
		}

	private:
		template<TensorIdx N>
		static constexpr Data &At(ArrayType<N> &array, std::span<const TensorIdx, N> indices) noexcept {
			if constexpr(N == 0)
				return array;
			else
				return At<N - 1>(array[indices[0]], indices.template subspan<1, N - 1>());
		}

		template<TensorIdx N>
		static constexpr const Data &At(const ArrayType<N> &array, std::span<const TensorIdx, N> indices) noexcept {
			if constexpr(N == 0)
				return array;
			else
				return At<N - 1>(array[indices[0]], indices.template subspan<1, N - 1>());
		}

		// Generic N-dimensional index iteration: calls f(indices) for all indices in [0..DIM)^N
		template<TensorIdx N, typename F>
			requires(std::is_invocable_v<F, const std::array<TensorIdx, N> &>)
		static constexpr void ForEachIndex(F &&f) {
			std::array<TensorIdx, N> indices{};
			ForEachIndex<N, F>(indices, std::forward<F>(f));
		}

		// Generic N-dimensional index iteration: calls f(indices) for all indices in [0..DIM)^N
		template<TensorIdx N, typename F, TensorIdx POS = 0>
			requires((0 <= POS && POS <= N) && std::is_invocable_v<F, const std::array<TensorIdx, N> &>)
		static constexpr void ForEachIndex(std::array<TensorIdx, N> &indices, F &&f) {
			if constexpr(POS == N) {
				f(indices);
			} else {
				for(TensorIdx i = 0; i < DIM; ++i) {
					indices[POS] = i;
					ForEachIndex<N, F, POS + 1>(indices, std::forward<F>(f));
				}
			}
		}

	private:
		ArrayType<RANK> mData{};
	};
}