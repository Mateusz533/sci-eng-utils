#pragma once
//
#include <algorithm>
#include <cmath>
#include <concepts>
#include <ratio>
#include <type_traits>
//
#include "Utils/Types.hpp"

namespace Utils::CompileTime
{
	using namespace Utils::Types;

	template<Arithmetic T, u8 N>
	static constexpr T Pow(T base) noexcept {
		if constexpr(N > 0) {
			return base * Pow<T, N - 1U>(base);
		} else {
			return T{1};
		}
	}

	static consteval i64 Exp10(u8 exp) noexcept {
		return (exp > 0) ? 10L * Exp10(exp - 1) : 1L;
	}

	template<std::floating_point T, u8 N>
		requires(N > 0)
	static constexpr T RootNewtonRaphson(T x, T curr, T prev) noexcept {
		return (curr == prev) ? curr : RootNewtonRaphson<T, N>(x, (T{N - 1U} * curr + x / Pow<T, N - 1U>(curr)) / T{N}, curr);
	}

	template<std::floating_point T>
	static constexpr T SqrtNewtonRaphson(T x, T curr, T prev) noexcept {
		return RootNewtonRaphson<T, 2>(x, curr, prev);
	}

	template<std::floating_point T>
	static constexpr T CbrtNewtonRaphson(T x, T curr, T prev) noexcept {
		return RootNewtonRaphson<T, 3>(x, curr, prev);
	}

	consteval i64 CeiledExp2(i64 value) noexcept {
		i64 result = 1;
		while(--value >= 0) {
			result *= 2L;
		}
		return result;
	}
	consteval i64 FlooredLog2(f64 value) noexcept {
		if(value <= 0) return 0;
		i64 log = 0;
		while(value > 1.0) {
			value /= 2.0;
			++log;
		}
		while(value < 1.0) {
			value *= 2.0;
			--log;
		}
		return log;
	}
	consteval i64 ExtractOddPart(i64 value) noexcept {
		while(value != 0 && value % 2 == 0) {
			value /= 2;
		}
		return value;
	}
	consteval i64 ExtractExp2(i64 value) noexcept {
		i64 log = 0;
		while(value != 0 && value % 2 == 0) {
			value /= 2;
			++log;
		}
		return log;
	}

	template<i64 NUMERATOR = 0, i64 DENOMINATOR = 1, i64 EXPONENT = 0>
		requires(DENOMINATOR > 0)
	class Fraction;

	template<typename T>
	concept NormedFraction = std::is_same_v<T, Fraction<T::NUM, T::DEN, T::EXP>>;

	template<typename Ratio, i64 EXP = 0>
	using FractionFromRatio = Fraction<Ratio::num, Ratio::den, EXP>::Norm;

	template<i8 EXP>
	using FractionFromExp10 = Fraction<Utils::CompileTime::Exp10(EXP > 0 ? +EXP : 0), Utils::CompileTime::Exp10(EXP < 0 ? -EXP : 0)>::Norm;

	template<i64 NUMERATOR, i64 DENOMINATOR, i64 EXPONENT>
		requires(DENOMINATOR > 0)
	class Fraction
	{
	public:
		static constexpr i64 NUM = std::ratio<ExtractOddPart(NUMERATOR), ExtractOddPart(DENOMINATOR)>::num;
		static constexpr i64 DEN = std::ratio<ExtractOddPart(NUMERATOR), ExtractOddPart(DENOMINATOR)>::den;
		static constexpr i64 EXP = (NUM == 0L) ? 0L : (EXPONENT + ExtractExp2(NUMERATOR) - ExtractExp2(DENOMINATOR));

		static consteval bool IsNormalized() noexcept {
			return NUM == NUMERATOR && DEN == DENOMINATOR && EXP == EXPONENT;
		}
		static consteval bool IsIdentity() noexcept {
			return NUMERATOR == 1 && DENOMINATOR == 1 && EXPONENT == 0;
		}
		static consteval bool IsZero() noexcept {
			return NUMERATOR == 0;
		}
		static consteval bool IsPositive() noexcept {
			return NUMERATOR > 0;
		}
		template<Arithmetic T = f64>
		static consteval T ToDecimal() noexcept {
			return static_cast<T>(NUM) * static_cast<T>(CeiledExp2(EXP)) /
				   (static_cast<T>(DEN) * static_cast<T>(CeiledExp2(-EXP)));
		}

		using Norm = Fraction<NUM, DEN, EXP>;
		using Opposite = Fraction<-NUM, DEN, EXP>;

		template<NormedFraction Other>
		using Sum = FractionFromRatio<std::ratio_add<std::ratio<NUM * CeiledExp2(EXP - Other::EXP), DEN>,
													 std::ratio<Other::NUM * CeiledExp2(Other::EXP - EXP), Other::DEN>>,
									  std::min(EXP, Other::EXP)>;
		template<NormedFraction Other>
		using Diff = FractionFromRatio<std::ratio_subtract<std::ratio<NUM * CeiledExp2(EXP - Other::EXP), DEN>,
														   std::ratio<Other::NUM * CeiledExp2(Other::EXP - EXP), Other::DEN>>,
									   std::min(EXP, Other::EXP)>;
		template<NormedFraction Other>
		using Product = FractionFromRatio<std::ratio_multiply<std::ratio<NUM, DEN>, std::ratio<Other::NUM, Other::DEN>>,
										  EXP + Other::EXP>::Norm;
		template<NormedFraction Other>
		using Quotient = FractionFromRatio<std::ratio_divide<std::ratio<NUM, DEN>, std::ratio<Other::NUM, Other::DEN>>,
										   EXP - Other::EXP>::Norm;
	};

	template<f64 VALUE>
		requires((VALUE == 0) || std::isnormal(VALUE))
	class FloatDecomposer
	{
	public:
		static consteval i64 Exponent() noexcept {
			return FlooredLog2(std::abs(VALUE)) - 52L;
		}
		static consteval i64 Numerator() noexcept {
			f64 value = VALUE;
			i64 exp = -Exponent();
			while(exp < 0) {
				value /= 2.0;
				++exp;
			}
			while(exp > 0) {
				value *= 2.0;
				--exp;
			}
			return static_cast<i64>(value);
		}
	};
}