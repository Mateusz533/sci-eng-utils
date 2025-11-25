#pragma once
//
#include <array>
#include <cstdint>
#include <type_traits>

namespace GenericMath
{
	using TensorIdx = int32_t;

	template<typename Data, TensorIdx RANK, TensorIdx DIM>
		requires(std::is_floating_point_v<Data> && DIM > 0 && RANK >= 0)
	class Tensor;

	template<typename Data, TensorIdx RANK, TensorIdx DIM>
		requires(std::is_floating_point_v<Data> && DIM > 0 && RANK >= 0)
	class Tensor
	{
		template<typename _Data, TensorIdx _RANK, TensorIdx _DIM>
			requires(std::is_floating_point_v<_Data> && _DIM > 0 && _RANK >= 0)
		friend class Tensor;

		template<TensorIdx _RANK>
			requires(_RANK >= 0)
		struct Helper {
			using ArrayType = std::array<typename Helper<_RANK - 1>::ArrayType, DIM>;
		};

		template<>
		struct Helper<0> {
			using ArrayType = Data;
		};

		template<TensorIdx _RANK>
		using ArrayType = typename Helper<_RANK>::ArrayType;

	public:
		template<typename... Indices>
			requires(sizeof...(Indices) == RANK && (std::is_same_v<Indices, TensorIdx> && ...))
		constexpr Data &operator()(Indices... indices) {
			return At<RANK>(mData, {indices...});
		}

		template<typename... Indices>
			requires(sizeof...(Indices) == RANK && (std::is_same_v<Indices, TensorIdx> && ...))
		constexpr const Data &operator()(Indices... indices) const {
			return At<RANK>(mData, {indices...});
		}

		constexpr operator Data() const
			requires(RANK == 0)
		{
			return mData;
		}

		// Outer product: (R) x (S) -> (R+S)
		template<TensorIdx _RANK>
		constexpr Tensor<Data, RANK + _RANK, DIM> OuterProduct(const Tensor<Data, _RANK, DIM> &other) const {
			Tensor<Data, RANK + _RANK, DIM> result{};
			std::array<TensorIdx, RANK + _RANK> rIdx{};

			ForEachIndex<RANK + _RANK>(rIdx, [&](const auto &ri) {
				const auto &ai = (const std::array<TensorIdx, RANK> &)ri;
				const auto &bi = (const std::array<TensorIdx, _RANK> &)ri[RANK];

				const Data a = At<RANK>(mData, ai);
				const Data b = At<_RANK>(other.mData, bi);
				At<RANK + _RANK>(result.mData, rIdx) = a * b;
			});

			return result;
		}

		// Inner product (single-index contraction between two tensors):
		// contract index i (this) with index j (other) -> (R+S-2)
		template<TensorIdx _RANK>
			requires(RANK >= 1 && _RANK >= 1)
		constexpr Tensor<Data, RANK + _RANK - 2, DIM> InnerProduct(const Tensor<Data, _RANK, DIM> &other, TensorIdx i, TensorIdx j) const {
			if(!(0 <= i && i < RANK && 0 <= j && j < _RANK)) {
				return Tensor<Data, RANK + _RANK - 2, DIM>{};
			}

			Tensor<Data, RANK + _RANK - 2, DIM> result{};
			std::array<TensorIdx, RANK> aIdx{};
			std::array<TensorIdx, _RANK> bIdx{};

			// mapping from (R+S-2) result indices to (R) and (S) indices excluding contracted positions
			std::array<int, RANK> mapA{};
			std::array<int, _RANK> mapB{};
			int cur = 0;
			for(TensorIdx p = 0; p < RANK; ++p) {
				if(p == i) {
					mapA[p] = -1;
					continue;
				}
				mapA[p] = cur++;
			}
			for(TensorIdx p = 0; p < _RANK; ++p) {
				if(p == j) {
					mapB[p] = -1;
					continue;
				}
				mapB[p] = cur++;
			}

			std::array<TensorIdx, RANK + _RANK - 2> rIdx{};

			// iterate over all result indices and sum over the contracted dimension
			ForEachIndex<RANK + _RANK - 2>(rIdx, [&](const auto &ri) {
				Data acc = Data{0};
				for(TensorIdx v = 0; v < DIM; ++v) {
					// build aIdx/bIdx from ri and contracted value v
					for(TensorIdx p = 0; p < RANK; ++p)
						aIdx[p] = (p == i) ? v : ri[mapA[p]];
					for(TensorIdx p = 0; p < _RANK; ++p)
						bIdx[p] = (p == j) ? v : ri[mapB[p]];

					acc += At<RANK>(mData, aIdx) * At<_RANK>(other.mData, bIdx);
				}
				At<RANK + _RANK - 2>(result.mData, ri) = acc;
			});

			return result;
		}

		// Contraction within a single tensor on two indices i and j: (R) -> (R-2)
		constexpr auto Contraction(TensorIdx i, TensorIdx j) const
			requires(RANK >= 2)
		{
			if(i == j || !(0 <= i && i < RANK && 0 <= j && j < RANK))
				return Tensor<Data, RANK - 2, DIM>{};

			// ensure i < j for simpler mapping
			if(j < i) std::swap(i, j);

			Tensor<Data, RANK - 2, DIM> result{};

			std::array<TensorIdx, RANK> srcIdx{};
			std::array<TensorIdx, RANK - 2> dstIdx{};

			// map positions excluding i and j into dstIdx
			std::array<TensorIdx, RANK> map{};
			TensorIdx cur = 0;
			for(TensorIdx p = 0; p < RANK; ++p) {
				if(p == i || p == j) {
					map[p] = -1;
					continue;
				}
				map[p] = cur++;
			}

			// iterate over all destination indices and sum over the diagonal (k on i and j)
			ForEachIndex<RANK - 2>(dstIdx, [&](const auto &di) {
				Data acc = Data{0};
				for(TensorIdx k = 0; k < DIM; ++k) {
					for(TensorIdx p = 0; p < RANK; ++p) {
						if(p == i || p == j)
							srcIdx[p] = k;
						else
							srcIdx[p] = di[map[p]];
					}
					acc += At<RANK>(mData, srcIdx);
				}
				At<RANK - 2>(result.mData, di) = acc;
			});

			return result;
		}

	private:
		// Access by runtime array of indices (works for any rank)
		template<TensorIdx N>
		static constexpr Data &At(ArrayType<N> &array, const std::array<TensorIdx, N> &indices) {
			if constexpr(N == 0)
				return array;
			else
				return At<N - 1>(array[indices[0]], (const std::array<TensorIdx, N - 1> &)indices[1]);
		}

		template<TensorIdx N>
		static constexpr const Data &At(const ArrayType<N> &array, const std::array<TensorIdx, N> &indices) {
			if constexpr(N == 0)
				return array;
			else
				return At<N - 1>(array[indices[0]], (const std::array<TensorIdx, N - 1> &)indices[1]);
		}

		// Generic N-dimensional index iteration: calls f(indices) for all indices in [0..DIM)^N
		template<TensorIdx N, typename F, TensorIdx POS = 0>
			requires((0 <= POS && POS <= N))
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