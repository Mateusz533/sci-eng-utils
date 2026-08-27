#pragma once
//
#include <cmath>
#include <cstdlib>
#include <numbers>
#include <ratio>
#include <type_traits>
//
#include "UnitsSI.hpp"

namespace Physics::Units::NStd
{
	namespace Detail
	{
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
				return static_cast<T>(NUM) / static_cast<T>(DEN) *
					   static_cast<T>(CeiledExp2(EXP)) / static_cast<T>(CeiledExp2(-EXP));	// TODO: Consider accuracy
			}

			using Norm = Fraction<NUM, DEN, EXP>;
			using Opposite = Fraction<-NUM, DEN, EXP>;

			template<NormedFraction Other>
			using Sum = Fraction<NUM * Other::DEN * CeiledExp2(EXP - Other::EXP) + DEN * Other::NUM * CeiledExp2(Other::EXP - EXP),
								 DEN * Other::DEN, std::min(EXP, Other::EXP)>::Norm;
			template<NormedFraction Other>
			using Diff = Fraction<NUM * Other::DEN * CeiledExp2(EXP - Other::EXP) - DEN * Other::NUM * CeiledExp2(Other::EXP - EXP),
								  DEN * Other::DEN, std::min(EXP, Other::EXP)>::Norm;
			template<NormedFraction Other>
			using Product = Fraction<std::ratio_multiply<std::ratio<NUM, DEN>, std::ratio<Other::NUM, Other::DEN>>::num,
									 std::ratio_multiply<std::ratio<NUM, DEN>, std::ratio<Other::NUM, Other::DEN>>::den,
									 EXP + Other::EXP>::Norm;
			template<NormedFraction Other>
			using Quotient = Fraction<std::ratio_divide<std::ratio<NUM, DEN>, std::ratio<Other::NUM, Other::DEN>>::num,
									  std::ratio_divide<std::ratio<NUM, DEN>, std::ratio<Other::NUM, Other::DEN>>::den,
									  EXP - Other::EXP>::Norm;
		};

		class SymbolicBase
		{};

		template<typename T>
		concept NonStandardUnit = std::is_base_of_v<SymbolicBase, T>;
	}

	template<typename StandardUnit, Detail::NormedFraction Scale, Detail::NormedFraction Offset, Arithmetic AccuracyType = f64>
		requires(!Detail::NonStandardUnit<StandardUnit> && StandardUnit::HasNoPrefix() && Scale::IsPositive() && !(Offset::IsZero() && Scale::IsIdentity()))
	class GenerativeUnit : Detail::SymbolicBase
	{
	private:
		using Type = decltype(std::declval<StandardUnit>().ToRaw());
		template<Arithmetic OtherType = Type>
		using SiblingStdUnit = decltype(std::declval<StandardUnit>().template Cast<OtherType>());

		using Self = GenerativeUnit;

		template<Arithmetic OtherType = Type>
		using Sibling = GenerativeUnit<SiblingStdUnit<OtherType>, Scale, Offset, AccuracyType>;

		template<typename OtherStandardUnit = StandardUnit, Detail::NormedFraction OtherScale = Scale,
				 Detail::NormedFraction OtherOffset = Offset, Arithmetic OtherAccuracyType = AccuracyType>
		using ExtendedSibling = GenerativeUnit<OtherStandardUnit, OtherScale, OtherOffset, OtherAccuracyType>;

		static consteval bool HasNoOffset() noexcept {
			return Offset::IsZero();
		}
		static consteval bool HasIdentScale() noexcept {
			return Scale::IsIdentity();
		}
		template<Detail::NormedFraction OtherScale>
		static consteval bool HasSameScale() noexcept {
			return std::is_same_v<Scale, OtherScale>;
		}
		template<Detail::NonStandardUnit OtherUnit>
		static consteval bool HasSameStdUnitBase() noexcept {
			using OtherType = decltype(std::declval<OtherUnit>().ToRaw());
			return std::is_same_v<decltype(std::declval<OtherUnit>().ToStandardUnit() / SiblingStdUnit<OtherType>()), SI::Scale<OtherType>>;
		}

	public:
		using OffsetType = Offset;
		using ScaleType = Scale;

		/* Constructors */;

		constexpr GenerativeUnit() = default;
		constexpr GenerativeUnit(const Self& value) = default;
		constexpr GenerativeUnit(Self&& value) = default;
		constexpr GenerativeUnit(Type data) noexcept : mData{data} {}
		template<Arithmetic OtherType = Type>
		constexpr GenerativeUnit(const SiblingStdUnit<OtherType>& value) noexcept : mData{FromStandardUnit(value)} {}
		template<Arithmetic OtherType = Type>
		constexpr GenerativeUnit(const Sibling<OtherType>& value) noexcept : mData{value.ToRaw()} {}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>())
		constexpr GenerativeUnit(const OtherUnit& value) noexcept : mData{FromOtherNStd(value)} {};

		/* Assignment operators */;

		constexpr Self& operator=(const Self& value) = default;
		constexpr Self& operator=(Self&& value) = default;
		template<Arithmetic OtherType = Type>
		constexpr Self& operator=(const SiblingStdUnit<OtherType>& value) noexcept {
			mData = FromStandardUnit(value);
			return *this;
		}
		template<Arithmetic OtherType = Type>
		constexpr Self& operator=(const Sibling<OtherType>& value) noexcept {
			mData = value.ToRaw();
			return *this;
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>())
		constexpr Self& operator=(const OtherUnit& value) noexcept {
			mData = FromOtherNStd(value);
			return *this;
		}

		template<Arithmetic OtherType = Type>
		constexpr Self& operator+=(const SiblingStdUnit<OtherType>& value) noexcept {
			mData += value.ToRaw() / SCALE;
			return *this;
		}
		template<Arithmetic OtherType = Type>
			requires(HasNoOffset())
		constexpr Self& operator+=(const Sibling<OtherType>& value) noexcept {
			mData += value.ToRaw();
			return *this;
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>() && OtherUnit::HasNoOffset())
		constexpr Self& operator+=(const OtherUnit& value) noexcept {
			mData += FromOtherNStd(value);
			return *this;
		}
		template<Arithmetic OtherType = Type>
		constexpr Self& operator-=(const SiblingStdUnit<OtherType>& value) noexcept {
			mData -= value.ToRaw() / SCALE;
			return *this;
		}
		template<Arithmetic OtherType = Type>
			requires(HasNoOffset())
		constexpr Self& operator-=(const Sibling<OtherType>& value) noexcept {
			mData -= value.ToRaw();
			return *this;
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>() && OtherUnit::HasNoOffset())
		constexpr Self& operator-=(const OtherUnit& value) noexcept {
			mData -= FromOtherNStd(value);
			return *this;
		}
		template<Arithmetic OtherType = Type>
			requires(HasNoOffset())
		constexpr Self& operator*=(const SI::Scale<OtherType>& value) noexcept {
			mData *= value.ToRaw();
			return *this;
		}
		template<Arithmetic OtherType = Type>
			requires(HasNoOffset())
		constexpr Self& operator/=(const SI::Scale<OtherType>& value) noexcept {
			mData /= value.ToRaw();
			return *this;
		}

		/* Comparison operators */;

		template<Arithmetic OtherType = Type>
		constexpr bool operator<(const Sibling<OtherType>& value) const noexcept {
			return mData < value.ToRaw();
		}
		template<Arithmetic OtherType = Type>
		constexpr bool operator>(const Sibling<OtherType>& value) const noexcept {
			return mData > value.ToRaw();
		}
		template<Arithmetic OtherType = Type>
		constexpr bool operator<=(const Sibling<OtherType>& value) const noexcept {
			return mData <= value.ToRaw();
		}
		template<Arithmetic OtherType = Type>
		constexpr bool operator>=(const Sibling<OtherType>& value) const noexcept {
			return mData >= value.ToRaw();
		}
		template<Arithmetic OtherType = Type>
		constexpr bool operator==(const Sibling<OtherType>& value) const noexcept {
			return mData == value.ToRaw();
		}
		template<Arithmetic OtherType = Type>
		constexpr bool operator!=(const Sibling<OtherType>& value) const noexcept {
			return mData != value.ToRaw();
		}
		template<Arithmetic OtherType = Type>
		constexpr auto operator<=>(const Sibling<OtherType>& value) const noexcept {
			return mData <=> value.ToRaw();
		}

		/* Logical operators */;

		constexpr bool operator!() const noexcept {
			return !mData;
		}

		/* Arithmetic operators */;

		constexpr auto operator-() const noexcept {
			using NewUnit = ExtendedSibling<SiblingStdUnit<decltype(-std::declval<Type>())>, Scale, typename Offset::Opposite>;
			return NewUnit{-mData};
		}

		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameScale<typename OtherUnit::ScaleType>() && HasSameStdUnitBase<OtherUnit>())
		constexpr auto operator+(const OtherUnit& value) const noexcept {
			using RawType = decltype(std::declval<Type>() + std::declval<OtherUnit>().ToRaw());
			using ResultOffset = Offset::template Sum<typename OtherUnit::OffsetType>;
			using AccType = decltype(std::declval<AccuracyType>() * OtherUnit::GetScale());

			if constexpr(Scale::IsIdentity() && ResultOffset::IsZero()) {
				return SiblingStdUnit<RawType>{mData + value.ToRaw()};
			} else {
				using NewUnit = ExtendedSibling<SiblingStdUnit<RawType>, Scale, ResultOffset, AccType>;
				return NewUnit{mData + value.ToRaw()};
			}
		}
		template<Arithmetic OtherType = Type>
			requires(HasIdentScale())
		constexpr auto operator+(const SiblingStdUnit<OtherType>& value) const noexcept {
			using NewUnit = Sibling<decltype(std::declval<Type>() + std::declval<OtherType>())>;
			return NewUnit{mData + value.ToRaw()};
		}
		template<Arithmetic OtherType = Type>
			requires(HasIdentScale())
		friend constexpr auto operator+(const SiblingStdUnit<OtherType>& value, const Self& self) noexcept {
			using NewUnit = Sibling<decltype(std::declval<OtherType>() + std::declval<Type>())>;
			return NewUnit{value.ToRaw() + self.ToRaw()};
		}

		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameScale<typename OtherUnit::ScaleType>() && HasSameStdUnitBase<OtherUnit>())
		constexpr auto operator-(const OtherUnit& value) const noexcept {
			using RawType = decltype(std::declval<Type>() - std::declval<OtherUnit>().ToRaw());
			using ResultOffset = Offset::template Diff<typename OtherUnit::OffsetType>;
			using AccType = decltype(std::declval<AccuracyType>() * OtherUnit::GetScale());

			if constexpr(Scale::IsIdentity() && ResultOffset::IsZero()) {
				return SiblingStdUnit<RawType>{mData - value.ToRaw()};
			} else {
				using NewUnit = ExtendedSibling<SiblingStdUnit<RawType>, Scale, ResultOffset, AccType>;
				return NewUnit{mData - value.ToRaw()};
			}
		}
		template<Arithmetic OtherType = Type>
			requires(HasIdentScale())
		constexpr auto operator-(const SiblingStdUnit<OtherType>& value) const noexcept {
			using NewUnit = Sibling<decltype(std::declval<Type>() - std::declval<OtherType>())>;
			return NewUnit{mData - value.ToRaw()};
		}
		template<Arithmetic OtherType = Type>
			requires(HasIdentScale())
		friend constexpr auto operator-(const SiblingStdUnit<OtherType>& value, const Self& self) noexcept {
			using NewUnit = decltype(-std::declval<Sibling<decltype(std::declval<OtherType>() - std::declval<Type>())>>());
			return NewUnit{value.ToRaw() - self.ToRaw()};
		}

		template<typename OtherStandardUnit, Detail::NormedFraction OtherScale, Detail::NormedFraction OtherOffset, Arithmetic OtherAccuracyType = f64>
			requires(HasNoOffset() && OtherOffset::IsZero() && !Detail::NonStandardUnit<OtherStandardUnit>)
		constexpr auto operator*(const ExtendedSibling<OtherStandardUnit, OtherScale, OtherOffset, OtherAccuracyType>& value) const noexcept {
			using NewStdUnit = decltype(std::declval<StandardUnit>() * std::declval<OtherStandardUnit>());
			using AccType = decltype(std::declval<AccuracyType>() * std::declval<OtherAccuracyType>());
			using ResultScale = Scale::template Product<OtherScale>;

			if constexpr(ResultScale::IsIdentity()) {
				return NewStdUnit{mData * value.ToRaw()};
			} else {
				using NewUnit = ExtendedSibling<NewStdUnit, ResultScale, Detail::Fraction<0, 1>, AccType>;
				return NewUnit{mData * value.ToRaw()};
			}
		}
		template<typename OtherStandardUnit = StandardUnit>
			requires(HasNoOffset() && !Detail::NonStandardUnit<OtherStandardUnit>)
		constexpr auto operator*(const OtherStandardUnit& value) const noexcept {
			using NewUnit = ExtendedSibling<decltype(std::declval<StandardUnit>() * std::declval<OtherStandardUnit>())>;
			return NewUnit{mData * value.ToRaw()};
		}
		template<typename OtherStandardUnit = StandardUnit>
			requires(HasNoOffset() && !Detail::NonStandardUnit<OtherStandardUnit>)
		friend constexpr auto operator*(const OtherStandardUnit& value, const Self& self) noexcept {
			using NewUnit = ExtendedSibling<decltype(std::declval<OtherStandardUnit>() * std::declval<StandardUnit>())>;
			return NewUnit{value.ToRaw() * self.ToRaw()};
		}

		template<typename OtherStandardUnit, Detail::NormedFraction OtherScale, Detail::NormedFraction OtherOffset, Arithmetic OtherAccuracyType = f64>
			requires(HasNoOffset() && OtherOffset::IsZero() && !Detail::NonStandardUnit<OtherStandardUnit>)
		constexpr auto operator/(const ExtendedSibling<OtherStandardUnit, OtherScale, OtherOffset, OtherAccuracyType> value) const noexcept {
			using NewStdUnit = decltype(std::declval<StandardUnit>() / std::declval<OtherStandardUnit>());
			using AccType = decltype(std::declval<AccuracyType>() / std::declval<OtherAccuracyType>());
			using ResultScale = Scale::template Quotient<OtherScale>;

			if constexpr(ResultScale::IsIdentity()) {
				return NewStdUnit{mData / value.ToRaw()};
			} else {
				using NewUnit = ExtendedSibling<NewStdUnit, ResultScale, Detail::Fraction<0, 1>, AccType>;
				return NewUnit{mData / value.ToRaw()};
			}
		}
		template<typename OtherStandardUnit = StandardUnit>
			requires(HasNoOffset() && !Detail::NonStandardUnit<OtherStandardUnit>)
		constexpr auto operator/(const OtherStandardUnit& value) const noexcept {
			using NewUnit = ExtendedSibling<decltype(std::declval<StandardUnit>() / std::declval<OtherStandardUnit>())>;
			return NewUnit{mData / value.ToRaw()};
		}
		template<typename OtherStandardUnit = StandardUnit>
			requires(HasNoOffset() && !Detail::NonStandardUnit<OtherStandardUnit>)
		friend constexpr auto operator/(const OtherStandardUnit& value, const Self& self) noexcept {
			using NewStdUnit = decltype(std::declval<OtherStandardUnit>() / std::declval<StandardUnit>());
			using NewUnit = ExtendedSibling<NewStdUnit, Detail::Fraction<1, 1>::Quotient<Scale>>;
			return NewUnit{value.ToRaw() / self.ToRaw()};
		}

		friend std::ostream& operator<<(std::ostream& os, const Self& obj) {
			os << obj.mData;
			return os;
		}

		/* Other methods */;

		template<Arithmetic OtherType = Type>
		constexpr Sibling<OtherType> Cast() const noexcept {
			return mData;
		}

		constexpr Type ToRaw() const noexcept {
			return mData;
		}

		template<Arithmetic OtherType = Type>
		constexpr SiblingStdUnit<OtherType> ToStandardUnit() const noexcept {
			if constexpr(SCALE == 1) {
				return mData + OFFSET;
			} else if constexpr(OFFSET == 0) {
				return mData * SCALE;
			} else {
				return mData * SCALE + OFFSET;
			}
		}

		static constexpr AccuracyType GetOffset() noexcept {
			return OFFSET;
		}

		static constexpr AccuracyType GetScale() noexcept {
			return SCALE;
		}

	private:
		template<Arithmetic OtherType = Type>
		static constexpr Type FromStandardUnit(const SiblingStdUnit<OtherType>& value) noexcept {
			if constexpr(SCALE == 1) {
				return value.ToRaw() - OFFSET;
			} else if constexpr(OFFSET == 0) {
				return value.ToRaw() / SCALE;
			} else {
				return (value.ToRaw() - OFFSET) / SCALE;
			}
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>())
		static constexpr Type FromOtherNStd(const OtherUnit& value) noexcept {
			if constexpr(std::is_same_v<Sibling<decltype(std::declval<OtherUnit>().ToRaw())>, OtherUnit>) {
				return value;
			} else {
				constexpr auto COMP_SCALE = OtherUnit::GetScale() * SCALE;
				constexpr auto COMP_OFFSET = (OtherUnit::GetOffset() - OFFSET) * SCALE;

				return value.ToRaw() * COMP_SCALE + COMP_OFFSET;
			}
		}

	private:
		/* TODO: Consider how to handle `AccuracyType` */;
		static constexpr AccuracyType SCALE = Scale::template ToDecimal<AccuracyType>();
		static constexpr AccuracyType OFFSET = Offset::template ToDecimal<AccuracyType>();
		Type mData;
	};

	namespace Detail
	{
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

		template<Arithmetic Type, template<Arithmetic> typename StdU, i64 ScNum, i64 ScDenom, i64 ScExp = 0, i64 OffNum = 0, i64 OffDenom = 1>
		using Simplifier = GenerativeUnit<StdU<Type>, typename Fraction<ScNum, ScDenom, ScExp>::Norm, typename Fraction<OffNum, OffDenom>::Norm>;
	}

#define GENERATE_NSTD_FROM_DOUBLE(Name, StdType, ScFloat, ScDenom)               \
	template<Arithmetic T = f64> /* NOLINTNEXTLINE(bugprone-macro-parentheses)*/ \
	using Name = Detail::Simplifier<T, StdType, Detail::FloatDecomposer<ScFloat>::Numerator(), ScDenom, Detail::FloatDecomposer<ScFloat>::Exponent()>;

#define GENERATE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom)               \
	template<Arithmetic T = f64> /* NOLINTNEXTLINE(bugprone-macro-parentheses)*/ \
	using Name = Detail::Simplifier<T, StdType, ScNum, ScDenom>;

#define GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom, OffNum, OffDenom) \
	template<Arithmetic T = f64> /* NOLINTNEXTLINE(bugprone-macro-parentheses)*/                 \
	using Name = Detail::Simplifier<T, StdType, ScNum, ScDenom, 0, OffNum, OffDenom>;

	/* --- GENERATE NEEDED UNITS HERE --- */

	/* Basic units */;

	GENERATE_NSTD_FROM_FRACTION(Percents, SI::Scale, 1, 100)
	GENERATE_NSTD_FROM_FRACTION(Permille, SI::Scale, 1, 1'000)
	GENERATE_NSTD_FROM_DOUBLE(Degrees, SI::Radians, std::numbers::pi, 180L)
	GENERATE_NSTD_FROM_DOUBLE(Arcminutes, SI::Radians, std::numbers::pi, 180L * 60)
	GENERATE_NSTD_FROM_DOUBLE(Arcseconds, SI::Radians, std::numbers::pi, 180L * 60 * 60)

	/* Distance */;

	GENERATE_NSTD_FROM_FRACTION(Mils, SI::Meters, 9'144L, 10'000L * 3 * 12 * 1'000)
	GENERATE_NSTD_FROM_FRACTION(Inches, SI::Meters, 9'144L, 10'000L * 3 * 12)
	GENERATE_NSTD_FROM_FRACTION(Hands, SI::Meters, 9'144L * 4, 10'000L * 3 * 12)
	GENERATE_NSTD_FROM_FRACTION(Links, SI::Meters, 9'144L * 22, 10'000L * 100)
	GENERATE_NSTD_FROM_FRACTION(Feet, SI::Meters, 9'144L, 10'000L * 3)
	GENERATE_NSTD_FROM_FRACTION(Yards, SI::Meters, 9'144L, 10'000L)
	GENERATE_NSTD_FROM_FRACTION(Fathoms, SI::Meters, 9'144L * 2, 10'000L)
	GENERATE_NSTD_FROM_FRACTION(Rods, SI::Meters, 9'144L * 55, 10'000L * 10)
	GENERATE_NSTD_FROM_FRACTION(Chains, SI::Meters, 9'144L * 22, 10'000L)
	GENERATE_NSTD_FROM_FRACTION(Furlongs, SI::Meters, 9'144L * 22 * 10, 10'000L)
	GENERATE_NSTD_FROM_FRACTION(StatuteMiles, SI::Meters, 9'144L * 22 * 10 * 8, 10'000L)
	GENERATE_NSTD_FROM_FRACTION(NauticalMiles, SI::Meters, 1'852, 1)

	/* Surface */;

	GENERATE_NSTD_FROM_FRACTION(SquareInches, SI::MetersSquare, 9'144L * 9'144, 10'000L * 3 * 12 * 10'000 * 3 * 12)
	GENERATE_NSTD_FROM_FRACTION(SquareFeet, SI::MetersSquare, 9'144L * 9'144, 10'000L * 3 * 10'000 * 3)
	GENERATE_NSTD_FROM_FRACTION(SquareYards, SI::MetersSquare, 9'144L * 9'144, 10'000L * 10'000)
	GENERATE_NSTD_FROM_FRACTION(SquareRods, SI::MetersSquare, 9'144L * 55 * 9'144 * 55, 10'000L * 10 * 10'000 * 10)
	GENERATE_NSTD_FROM_FRACTION(SquareChains, SI::MetersSquare, 9'144L * 22 * 9'144 * 22, 10'000L * 10'000)
	GENERATE_NSTD_FROM_FRACTION(SquareStatuteMiles, SI::MetersSquare, 9'144L * 22 * 10 * 8 * 9'144 * 22 * 10 * 8, 10'000L * 10'000)
	GENERATE_NSTD_FROM_FRACTION(Roods, SI::MetersSquare, 9'144L * 9'144 * 4840, 10'000L * 10'000 * 4)
	GENERATE_NSTD_FROM_FRACTION(Acres, SI::MetersSquare, 9'144L * 9'144 * 4840, 10'000L * 10'000)

	/* Volume */;

	GENERATE_NSTD_FROM_FRACTION(CubicInches, SI::MetersCube, 9'144L * 9'144 * 9'144, 10'000L * 3 * 12 * 10'000 * 3 * 12 * 10'000 * 3 * 12)
	GENERATE_NSTD_FROM_FRACTION(CubicFeet, SI::MetersCube, 9'144L * 9'144 * 9'144, 10'000L * 3 * 10'000 * 3 * 10'000 * 3)
	GENERATE_NSTD_FROM_FRACTION(CubicYards, SI::MetersCube, 9'144L * 9'144 * 9'144, 10'000L * 10'000 * 10'000)
	GENERATE_NSTD_FROM_FRACTION(AcreFeet, SI::MetersCube, 9'144L * 9'144 * 4840 * 9'144, 10'000L * 10'000 * 10'000L * 3)
	GENERATE_NSTD_FROM_FRACTION(BoardFeet, SI::MetersCube, 9'144L * 9'144 * 9'144, 10'000L * 3 * 10'000 * 3 * 10'000 * 3 * 12)
	GENERATE_NSTD_FROM_FRACTION(Cords, SI::MetersCube, 9'144L * 9'144 * 9'144 * 128, 10'000L * 3 * 10'000 * 3 * 10'000 * 3)

	GENERATE_NSTD_FROM_FRACTION(UsMinims, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 128 * 8 * 60)
	GENERATE_NSTD_FROM_FRACTION(UsFluidDrams, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 128 * 8)
	GENERATE_NSTD_FROM_FRACTION(UsFluidOunces, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 128)
	GENERATE_NSTD_FROM_FRACTION(UsGills, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 32)
	GENERATE_NSTD_FROM_FRACTION(UsCups, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 16)
	GENERATE_NSTD_FROM_FRACTION(UsMeasures, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 16)
	GENERATE_NSTD_FROM_FRACTION(UsLiquidPints, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 8)
	GENERATE_NSTD_FROM_FRACTION(UsLiquidQuarts, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L * 4)
	GENERATE_NSTD_FROM_FRACTION(UsGallons, SI::MetersCube, 3'785'411'784L, 1'000'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(PetroleumBarrels, SI::MetersCube, 3'785'411'784L * 42, 1'000'000'000'000L)

	GENERATE_NSTD_FROM_FRACTION(UsDryPints, SI::MetersCube, 2'150'420L * 16'387'064, 1'000L * 1'000'000'000'000 * 64)
	GENERATE_NSTD_FROM_FRACTION(UsDryQuarts, SI::MetersCube, 2'150'420L * 16'387'064, 1'000L * 1'000'000'000'000 * 32)
	GENERATE_NSTD_FROM_FRACTION(UsPecks, SI::MetersCube, 2'150'420L * 16'387'064, 1'000L * 1'000'000'000'000 * 4)
	GENERATE_NSTD_FROM_FRACTION(UsBushels, SI::MetersCube, 2'150'420L * 16'387'064, 1'000L * 1'000'000'000'000)

	GENERATE_NSTD_FROM_FRACTION(ImperialFluidMinims, SI::MetersCube, 4'546'090L, 1'000'000'000L * 8 * 20 * 8 * 60)
	GENERATE_NSTD_FROM_FRACTION(ImperialFluidDrachms, SI::MetersCube, 4'546'090L, 1'000'000'000L * 8 * 20 * 8)
	GENERATE_NSTD_FROM_FRACTION(ImperialFluidOunces, SI::MetersCube, 4'546'090L, 1'000'000'000L * 8 * 20)
	GENERATE_NSTD_FROM_FRACTION(ImperialGills, SI::MetersCube, 4'546'090L, 1'000'000'000L * 32)
	GENERATE_NSTD_FROM_FRACTION(ImperialPints, SI::MetersCube, 4'546'090L, 1'000'000'000L * 8)
	GENERATE_NSTD_FROM_FRACTION(ImperialQuarts, SI::MetersCube, 4'546'090L, 1'000'000'000L * 4)
	GENERATE_NSTD_FROM_FRACTION(ImperialGallons, SI::MetersCube, 4'546'090L, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(ImperialPecks, SI::MetersCube, 4'546'090L * 2, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(ImperialBushels, SI::MetersCube, 4'546'090L * 8, 1'000'000'000L)

	/* Time */;

	GENERATE_NSTD_FROM_FRACTION(Minutes, SI::Seconds, 60L, 1L)
	GENERATE_NSTD_FROM_FRACTION(Hours, SI::Seconds, 60L * 60, 1L)
	GENERATE_NSTD_FROM_FRACTION(Days, SI::Seconds, 60L * 60 * 24, 1L)
	GENERATE_NSTD_FROM_FRACTION(Weeks, SI::Seconds, 60L * 60 * 24 * 7, 1L)

	/* Velocity */;

	GENERATE_NSTD_FROM_FRACTION(Knots, SI::MetersPerSecond, 1'852L, 60L * 60)
	GENERATE_NSTD_FROM_FRACTION(MilesPerHour, SI::MetersPerSecond, 9'144L * 22 * 10 * 8, 10'000L * 60 * 60)
	GENERATE_NSTD_FROM_FRACTION(FeetPerMinute, SI::MetersPerSecond, 9'144L, 10'000L * 3 * 60)
	GENERATE_NSTD_FROM_FRACTION(FeetPerSecond, SI::MetersPerSecond, 9'144L, 10'000L * 3)

	/* Mass & Force */;

	GENERATE_NSTD_FROM_FRACTION(Grains, SI::KiloGrams, 453'592'370L, 1'000'000'000L * 16 * 16 * 7'000)
	GENERATE_NSTD_FROM_FRACTION(Pennyweights, SI::KiloGrams, 453'592'370L * 24, 1'000'000'000L * 16 * 16 * 7'000)
	GENERATE_NSTD_FROM_FRACTION(TroyOunces, SI::KiloGrams, 453'592'370L * 24 * 20, 1'000'000'000L * 16 * 16 * 7'000)
	GENERATE_NSTD_FROM_FRACTION(TroyPounds, SI::KiloGrams, 453'592'370L * 24 * 20 * 12, 1'000'000'000L * 16 * 16 * 7'000)
	GENERATE_NSTD_FROM_FRACTION(Drams, SI::KiloGrams, 453'592'370L, 1'000'000'000L * 16 * 16)
	GENERATE_NSTD_FROM_FRACTION(Ounces, SI::KiloGrams, 453'592'370L, 1'000'000'000L * 16)
	GENERATE_NSTD_FROM_FRACTION(PoundsMass, SI::KiloGrams, 453'592'370L, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(Stones, SI::KiloGrams, 453'592'370L * 14, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(Quarters, SI::KiloGrams, 453'592'370L * 28, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(Centals, SI::KiloGrams, 453'592'370L * 100, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(ShortHundredweights, SI::KiloGrams, 453'592'370L * 100, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(LongHundredweights, SI::KiloGrams, 453'592'370L * 112, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(ShortTons, SI::KiloGrams, 453'592'370L * 2'000, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(LongTons, SI::KiloGrams, 453'592'370L * 2'240, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(Slug, SI::KiloGrams, 453'592'370L * 9'806'650, 3'048L * 100'000 * 1'000'000)

	GENERATE_NSTD_FROM_FRACTION(StandardGravity, SI::MetersPerSecondSquare, 9'806'650L, 1'000'000L)
	GENERATE_NSTD_FROM_FRACTION(OuncesForce, SI::Newtons, 453'592'370L * 9'806'650, 1'000'000'000L * 1'000'000 * 16)
	GENERATE_NSTD_FROM_FRACTION(PoundsForce, SI::Newtons, 453'592'370L * 9'806'650, 1'000'000'000L * 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(Kips, SI::Newtons, 453'592'370L * 9'806'650 * 1'000, 1'000'000'000L * 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(Poundals, SI::Newtons, 453'592'370L * 9'144, 1'000'000'000L * 10'000 * 3)

	/* Torque, Energy, Power */;

	GENERATE_NSTD_FROM_FRACTION(PoundForceFeet, SI::NewtonMeters, 13'558'179'483'314'004L, 10'000'000'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(PoundForceInches, SI::NewtonMeters, 13'558'179'483'314'004L, 10'000'000'000'000'000L * 12)
	GENERATE_NSTD_FROM_FRACTION(FootPoundsForce, SI::Joules, 13'558'179'483'314'004L, 10'000'000'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(UsTherms, SI::Joules, 105'480'400L, 1L)
	GENERATE_NSTD_FROM_FRACTION(BtuIt, SI::Joules, 1'055'055'852'620L, 1'000'000'000L)
	GENERATE_NSTD_FROM_FRACTION(BtuTh, SI::Joules, 4'184L, 1'000L)
	GENERATE_NSTD_FROM_FRACTION(TonsOfRefrigeration, SI::Watts, 1'055'055'852'620L * 12'000, 1'000'000'000L * 60 * 60)
	GENERATE_NSTD_FROM_FRACTION(Horsepower, SI::Watts, 13'558'179'483'314'004L * 550, 10'000'000'000'000'000L)

	/* Pressure */;

	GENERATE_NSTD_FROM_FRACTION(InchesOfMercury, SI::Pascals, 3'386'389L, 1'000L)
	GENERATE_NSTD_FROM_FRACTION(Bars, SI::Pascals, 100'000L, 1L)
	GENERATE_NSTD_FROM_FRACTION(MilliBars, SI::Pascals, 100'000L, 1'000L)
	GENERATE_NSTD_FROM_FRACTION(Atmospheres, SI::Pascals, 101'325L, 1L)
	GENERATE_NSTD_FROM_FRACTION(PoundsPerSquareFoot, SI::Pascals, 453'592'370L * 9'806'650, 25'400L * 25'400 * 1'000 * 12 * 12)
	GENERATE_NSTD_FROM_FRACTION(PoundsPerSquareInch, SI::Pascals, 453'592'370L * 9'806'650, 25'400L * 25'400 * 1'000)
	GENERATE_NSTD_FROM_FRACTION(KipsPerSquareInch, SI::Pascals, 453'592'370L * 9'806'650 * 1'000, 25'400L * 25'400 * 1'000)

	/* Temperature */;

	GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION(DegreesCelsius, SI::Kelvins, 1, 1, 273'150, 1'000)
	GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION(DegreesFahrenheit, SI::Kelvins, 5, 9, 9 * 273'150L - 4 * 40'000L, 9 * 1'000L)
	GENERATE_NSTD_FROM_FRACTION(DegreesRankin, SI::Kelvins, 5, 9)

#undef GENERATE_NSTD_FROM_DOUBLE
#undef GENERATE_NSTD_FROM_FRACTION
#undef GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION
}