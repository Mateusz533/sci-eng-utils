#pragma once
//
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
		static consteval i64 CeiledExp2(i64 value) noexcept {
			i64 result = 1;
			while(--value >= 0) {
				result *= 2L;
			}
			return result;
		}
		static consteval i64 FlooredLog2(f64 value) noexcept {
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
		static consteval i64 ExtractOddPart(i64 value) noexcept {
			while(value != 0 && value % 2 == 0) {
				value /= 2;
			}
			return value;
		}
		static consteval i64 ExtractExp2(i64 value) noexcept {
			i64 log = 0;
			while(value != 0 && value % 2 == 0) {
				value /= 2;
				++log;
			}
			return log;
		}

		template<i64 Numerator = 0, i64 Denominator = 1, i64 Exponent = 0>
			requires(Denominator > 0)
		class Fraction;

		template<class T>
		concept NormedFraction = std::is_same_v<T, Fraction<T::Num, T::Denom, T::Exp>>;

		template<i64 Numerator, i64 Denominator, i64 Exponent>
			requires(Denominator > 0)
		class Fraction
		{
		public:
			static constexpr i64 Num = std::ratio<ExtractOddPart(Numerator), ExtractOddPart(Denominator)>::num;
			static constexpr i64 Denom = std::ratio<ExtractOddPart(Numerator), ExtractOddPart(Denominator)>::den;
			static constexpr i64 Exp = (Num == 0L) ? 0L : (Exponent + ExtractExp2(Numerator) - ExtractExp2(Denominator));

			static consteval bool IsNormalized() noexcept {
				return Num == Numerator && Denom == Denominator && Exp == Exponent;
			}
			static consteval bool IsIdentity() noexcept {
				return Numerator == 1 && Denominator == 1 && Exponent == 0;
			}
			static consteval bool IsZero() noexcept {
				return Numerator == 0;
			}
			static consteval bool IsPositive() noexcept {
				return Numerator > 0;
			}
			template<Arithmetic T = f64>
			static consteval T ToDecimal() noexcept {
				return static_cast<T>(Num) / static_cast<T>(Denom) *
					   static_cast<T>(CeiledExp2(Exp)) / static_cast<T>(CeiledExp2(-Exp));	// TODO: Consider accuracy
			}

			using Norm = Fraction<Num, Denom, Exp>;
			using Opposite = Fraction<-Num, Denom, Exp>;

			template<NormedFraction Other>
			using Sum = Fraction<Num * Other::Denom * CeiledExp2(Exp - Other::Exp) + Denom * Other::Num * CeiledExp2(Other::Exp - Exp),
								 Denom * Other::Denom, std::min(Exp, Other::Exp)>::Norm;
			template<NormedFraction Other>
			using Diff = Fraction<Num * Other::Denom * CeiledExp2(Exp - Other::Exp) - Denom * Other::Num * CeiledExp2(Other::Exp - Exp),
								  Denom * Other::Denom, std::min(Exp, Other::Exp)>::Norm;
			template<NormedFraction Other>
			using Product = Fraction<Num * Other::Num, Denom * Other::Denom, Exp + Other::Exp>::Norm;
			template<NormedFraction Other>
			using Quotient = Fraction<Num * Other::Denom, Denom * Other::Num, Exp - Other::Exp>::Norm;
		};

		class SymbolicBase
		{};

		template<class T>
		concept NonStandardUnit = std::is_base_of_v<SymbolicBase, T>;
	}

	template<Arithmetic Type, template<Arithmetic> class StandardUnit, Detail::NormedFraction Scale, Detail::NormedFraction Offset, Arithmetic AccuracyType = f64>
		requires(!Detail::NonStandardUnit<StandardUnit<i8>> && StandardUnit<i8>::HasNoPrefix() && Scale::IsPositive() && !(Offset::IsZero() && Scale::IsIdentity()))
	class GenerativeUnit : Detail::SymbolicBase
	{
	private:
		static consteval bool HasNoOffset() noexcept {
			return Offset::IsZero();
		}
		static consteval bool HasIdentScale() noexcept {
			return Scale::IsIdentity();
		}
		template<Detail::NormedFraction _Scale>
		static consteval bool HasSameScale() noexcept {
			return std::is_same_v<Scale, _Scale>;
		}

		template<Arithmetic _Type = Type, template<Arithmetic> class _StandardUnit = StandardUnit>
		using Sibling = GenerativeUnit<_Type, _StandardUnit, Scale, Offset, AccuracyType>;

		using Self = GenerativeUnit;

		template<Arithmetic _Type = Type, template<Arithmetic> class _StandardUnit = StandardUnit, Detail::NormedFraction _Scale = Scale,
				 Detail::NormedFraction _Offset = Offset, Arithmetic _AccuracyType = AccuracyType>
		using AdaptiveSibling = std::conditional_t<_Scale::IsIdentity() && _Offset::IsZero(), _StandardUnit<_Type>,
												   GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType>>;

		template<class Base>
		struct Decomposer {
			template<Arithmetic T>
			using OuterType = decltype(Base{}.template Cast<T>());
			using InnerType = decltype(Base{}.ToRaw());
		};
		template<class Complex>
		struct StdToSibling {
			using NewUnit = Sibling<decltype(Complex{}.ToRaw()), Decomposer<Complex>::template OuterType>;
		};
		template<Detail::NonStandardUnit OtherUnit>
		static consteval bool HasSameStdUnitBase() noexcept {
			return std::is_same_v<decltype(OtherUnit{}.ToStandardUnit() / StandardUnit<decltype(OtherUnit{}.ToRaw())>()),
								  SI::Scale<decltype(OtherUnit{}.ToRaw())>>;
		}

	public:
		/* Constructors */;

		constexpr GenerativeUnit() = default;
		constexpr GenerativeUnit(const Self& value) = default;
		constexpr GenerativeUnit(Self&& value) = default;
		constexpr GenerativeUnit(Type data) noexcept : mData(data) {}
		template<Arithmetic _Type = Type>
		constexpr GenerativeUnit(const StandardUnit<_Type>& value) noexcept : mData(FromStandardUnit(value)) {}
		template<Arithmetic _Type = Type>
		constexpr GenerativeUnit(const Sibling<_Type>& value) noexcept : mData(value.ToRaw()) {}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>())
		constexpr GenerativeUnit(const OtherUnit& value) noexcept : mData(FromOtherNStd(value)){};

		/* Assignment operators */;

		constexpr Self& operator=(const Self& value) = default;
		constexpr Self& operator=(Self&& value) = default;
		template<Arithmetic _Type = Type>
		constexpr Self& operator=(const StandardUnit<_Type>& value) noexcept {
			mData = FromStandardUnit(value);
			return *this;
		}
		template<Arithmetic _Type = Type>
		constexpr Self& operator=(const Sibling<_Type>& value) noexcept {
			mData = value.ToRaw();
			return *this;
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>())
		constexpr Self& operator=(const OtherUnit& value) noexcept {
			mData = FromOtherNStd(value);
			return *this;
		}

		template<Arithmetic _Type = Type>
		constexpr Self& operator+=(const StandardUnit<_Type>& value) noexcept {
			mData += value.ToRaw() / SCALE;
			return *this;
		}
		template<Arithmetic _Type = Type>
			requires(HasNoOffset())
		constexpr Self& operator+=(const Sibling<_Type>& value) noexcept {
			mData += value.ToRaw();
			return *this;
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>() && OtherUnit::HasNoOffset())
		constexpr Self& operator+=(const OtherUnit& value) noexcept {
			mData += FromOtherNStd(value);
			return *this;
		}
		template<Arithmetic _Type = Type>
		constexpr Self& operator-=(const StandardUnit<_Type>& value) noexcept {
			mData -= value.ToRaw() / SCALE;
			return *this;
		}
		template<Arithmetic _Type = Type>
			requires(HasNoOffset())
		constexpr Self& operator-=(const Sibling<_Type>& value) noexcept {
			mData -= value.ToRaw();
			return *this;
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(HasSameStdUnitBase<OtherUnit>() && OtherUnit::HasNoOffset())
		constexpr Self& operator-=(const OtherUnit& value) noexcept {
			mData -= FromOtherNStd(value);
			return *this;
		}
		template<Arithmetic _Type = Type>
			requires(HasNoOffset())
		constexpr Self& operator*=(const SI::Scale<_Type>& value) noexcept {
			mData *= value.ToRaw();
			return *this;
		}
		template<Arithmetic _Type = Type>
			requires(HasNoOffset())
		constexpr Self& operator/=(const SI::Scale<_Type>& value) noexcept {
			mData /= value.ToRaw();
			return *this;
		}

		/* Comparison operators */;

		template<Arithmetic _Type = Type>
		constexpr bool operator<(const Sibling<_Type>& value) const noexcept {
			return mData < value.ToRaw();
		}
		template<Arithmetic _Type = Type>
		constexpr bool operator>(const Sibling<_Type>& value) const noexcept {
			return mData > value.ToRaw();
		}
		template<Arithmetic _Type = Type>
		constexpr bool operator<=(const Sibling<_Type>& value) const noexcept {
			return mData <= value.ToRaw();
		}
		template<Arithmetic _Type = Type>
		constexpr bool operator>=(const Sibling<_Type>& value) const noexcept {
			return mData >= value.ToRaw();
		}
		template<Arithmetic _Type = Type>
		constexpr bool operator==(const Sibling<_Type>& value) const noexcept {
			return mData == value.ToRaw();
		}
		template<Arithmetic _Type = Type>
		constexpr bool operator!=(const Sibling<_Type>& value) const noexcept {
			return mData != value.ToRaw();
		}
		template<Arithmetic _Type = Type>
		constexpr auto operator<=>(const Sibling<_Type>& value) const noexcept {
			return mData <=> value.ToRaw();
		}

		/* Logical operators */;

		constexpr bool operator!() const noexcept {
			return !mData;
		}

		/* Arithmetic operators */;

		constexpr AdaptiveSibling<decltype(-Type{}), StandardUnit, Scale, typename Offset::Opposite> operator-() const noexcept {
			return -mData;
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(OtherUnit::template HasSameScale<Scale>() && HasSameStdUnitBase<OtherUnit>())
		constexpr auto operator+(const OtherUnit& value) const noexcept {
			using RawType = decltype(Type{} + OtherUnit{}.ToRaw());
			using ResultOffset = Offset::template Sum<typename OtherUnit::OffsetT>;
			using AccType = decltype(AccuracyType{} * OtherUnit::SCALE);
			using NewUnit = AdaptiveSibling<RawType, StandardUnit, Scale, ResultOffset, AccType>;
			return NewUnit{mData + value.ToRaw()};
		}
		template<Arithmetic _Type = Type>
			requires(HasIdentScale())
		constexpr Sibling<decltype(Type{} + _Type{})> operator+(const StandardUnit<_Type>& value) const noexcept {
			return mData + value.ToRaw();
		}
		template<Arithmetic _Type = Type>
			requires(HasIdentScale())
		friend constexpr Sibling<decltype(_Type{} + Type{})> operator+(const StandardUnit<_Type>& value, const Self& self) noexcept {
			return value.ToRaw() + self.ToRaw();
		}
		template<Detail::NonStandardUnit OtherUnit = Self>
			requires(OtherUnit::template HasSameScale<Scale>() && HasSameStdUnitBase<OtherUnit>())
		constexpr auto operator-(const OtherUnit& value) const noexcept {
			using RawType = decltype(Type{} - OtherUnit{}.ToRaw());
			using ResultOffset = Offset::template Diff<typename OtherUnit::OffsetT>;
			using AccType = decltype(AccuracyType{} * OtherUnit::SCALE);
			using NewUnit = AdaptiveSibling<RawType, StandardUnit, Scale, ResultOffset, AccType>;
			return NewUnit{mData - value.ToRaw()};
		}
		template<Arithmetic _Type = Type>
			requires(HasIdentScale())
		constexpr Sibling<decltype(Type{} - _Type{})> operator-(const StandardUnit<_Type>& value) const noexcept {
			return mData - value.ToRaw();
		}
		template<Arithmetic _Type = Type>
			requires(HasIdentScale())
		friend constexpr auto operator-(const StandardUnit<_Type>& value, const Self& self) noexcept {
			using NewUnit = decltype(-Sibling<decltype(_Type{} - Type{})>{});
			return NewUnit{value.ToRaw() - self.ToRaw()};
		}
		template<Arithmetic _Type, template<Arithmetic> class _StandardUnit, Detail::NormedFraction _Scale, Detail::NormedFraction _Offset, Arithmetic _AccuracyType = f64>
			requires(HasNoOffset() && _Offset::IsZero())
		constexpr auto operator*(const GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType>& value) const noexcept {
			using RawType = decltype(Type{} * _Type{});
			using ComposeType = decltype(StandardUnit<Type>{} * _StandardUnit<_Type>{});
			using AccType = decltype(AccuracyType{} * _AccuracyType{});
			using ResultScale = Scale::template Product<_Scale>;
			using NewUnit = AdaptiveSibling<RawType, Decomposer<ComposeType>::template OuterType, ResultScale, Detail::Fraction<0, 1>, AccType>;
			return NewUnit{mData * value.ToRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		constexpr auto operator*(const _StdUnitT& value) const noexcept {
			using NewUnit = StdToSibling<decltype(StandardUnit<Type>{} * _StdUnitT{})>::NewUnit;
			return NewUnit{mData * value.ToRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		friend constexpr auto operator*(const _StdUnitT& value, const Self& self) noexcept {
			using NewUnit = StdToSibling<decltype(_StdUnitT{} * StandardUnit<Type>{})>::NewUnit;
			return NewUnit{value.ToRaw() * self.ToRaw()};
		}
		template<Arithmetic _Type, template<Arithmetic> class _StandardUnit, Detail::NormedFraction _Scale, Detail::NormedFraction _Offset, Arithmetic _AccuracyType = f64>
			requires(HasNoOffset() && _Offset::IsZero())
		constexpr auto operator/(const GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType> value) const noexcept {
			using RawType = decltype(Type{} / _Type{});
			using ComposeType = decltype(StandardUnit<Type>{} / _StandardUnit<_Type>{});
			using AccType = decltype(AccuracyType{} / _AccuracyType{});
			using ResultScale = Scale::template Quotient<_Scale>;
			using NewUnit = AdaptiveSibling<RawType, Decomposer<ComposeType>::template OuterType, ResultScale, Detail::Fraction<0, 1>, AccType>;
			return NewUnit{mData / value.ToRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		constexpr auto operator/(const _StdUnitT& value) const noexcept {
			using NewUnit = StdToSibling<decltype(StandardUnit<Type>{} / _StdUnitT{})>::NewUnit;
			return NewUnit{mData / value.ToRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		friend constexpr auto operator/(const _StdUnitT& value, const Self& self) noexcept {
			using NewUnit = StdToSibling<decltype(_StdUnitT{} / StandardUnit<Type>{})>::NewUnit;
			return NewUnit{value.ToRaw() / self.ToRaw()};
		}

		friend std::ostream& operator<<(std::ostream& os, const Self& obj) {
			os << obj.mData;
			return os;
		}

		/* Other methods */;

		template<Arithmetic _Type = Type>
		constexpr Sibling<_Type> Cast() const noexcept {
			return mData;
		}

		constexpr Type ToRaw() const noexcept {
			return mData;
		}

		template<Arithmetic _Type = Type>
		constexpr StandardUnit<_Type> ToStandardUnit() const noexcept {
			if constexpr(SCALE == 1) {
				return mData + OFFSET;
			} else if constexpr(OFFSET == 0) {
				return mData * SCALE;
			} else {
				return mData * SCALE + OFFSET;
			}
		}

		using OffsetT = Offset;
		using ScaleT = Scale;

	private:
		template<Arithmetic _Type = Type>
		static constexpr Type FromStandardUnit(const StandardUnit<_Type>& value) noexcept {
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
			if constexpr(std::is_same_v<Sibling<decltype(OtherUnit{}.ToRaw())>, OtherUnit>) {
				return value;
			} else {
				constexpr auto COMP_SCALE = OtherUnit::SCALE * SCALE;
				constexpr auto COMP_OFFSET = (OtherUnit::OFFSET - OFFSET) * SCALE;

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
		class FractionParser
		{
		public:
			static consteval i64 Exponent() noexcept {
				return FlooredLog2(std::abs(VALUE)) - 52L;
			}
			static consteval i64 Denominator() noexcept {
				return 1L;
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

		template<Arithmetic Type, template<Arithmetic> class StdU, i64 ScNum, i64 ScDenom, i64 ScExp = 0, i64 OffNum = 0, i64 OffDenom = 1>
		using Simplifier = GenerativeUnit<Type, StdU, typename Fraction<ScNum, ScDenom, ScExp>::Norm, typename Fraction<OffNum, OffDenom>::Norm>;

		inline constexpr f64 DEG2RAD = std::numbers::pi / 180.0;
	}

#define GENERATE_NSTD_FROM_DOUBLE(Name, StdType, Scale) \
	template<Arithmetic T = f64>                        \
	using Name = Detail::Simplifier<T, StdType, Detail::FractionParser<Scale>::Numerator(), Detail::FractionParser<Scale>::Denominator(), Detail::FractionParser<Scale>::Exponent()>;

#define GENERATE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom) \
	template<Arithmetic T = f64>                                   \
	using Name = Detail::Simplifier<T, StdType, ScNum, ScDenom>;

#define GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom, OffNum, OffDenom) \
	template<Arithmetic T = f64>                                                                 \
	using Name = Detail::Simplifier<T, StdType, ScNum, ScDenom, 0, OffNum, OffDenom>;

	/* --- GENERATE NEEDED UNITS HERE --- */

	/* Basic units */;

	GENERATE_NSTD_FROM_FRACTION(Percents, SI::Scale, 1, 100)
	GENERATE_NSTD_FROM_FRACTION(Permille, SI::Scale, 1, 1'000)
	GENERATE_NSTD_FROM_DOUBLE(Degrees, SI::Radians, Detail::DEG2RAD)
	// GENERATE_NSTD_FROM_DOUBLE(Arcminutes, SI::Radians, Detail::DEG2RAD / 60.0)
	// GENERATE_NSTD_FROM_DOUBLE(Arcseconds, SI::Radians, Detail::DEG2RAD / 3'600.0)

	/* Distance */;

	GENERATE_NSTD_FROM_FRACTION(Inches, SI::Meters, 25'400, 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(Feet, SI::Meters, 304'800, 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(Yards, SI::Meters, 914'400, 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(Chains, SI::Meters, 20'116'800, 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(Furlongs, SI::Meters, 201'168, 1'000)
	GENERATE_NSTD_FROM_FRACTION(StatuteMiles, SI::Meters, 1'609'344, 1'000)
	GENERATE_NSTD_FROM_FRACTION(NauticalMiles, SI::Meters, 1'852, 1)

	/* Time */;

	GENERATE_NSTD_FROM_FRACTION(Minutes, SI::Seconds, 60, 1)
	GENERATE_NSTD_FROM_FRACTION(Hours, SI::Seconds, 3'600, 1)
	GENERATE_NSTD_FROM_FRACTION(Days, SI::Seconds, 86'400, 1)
	GENERATE_NSTD_FROM_FRACTION(Weeks, SI::Seconds, 604'800, 1)

	/* Velocity */;

	GENERATE_NSTD_FROM_FRACTION(Knots, SI::MetersPerSecond, 1'852, 3'600)
	GENERATE_NSTD_FROM_FRACTION(MilesPerHour, SI::MetersPerSecond, 1'609'344, 3'600'000)
	GENERATE_NSTD_FROM_FRACTION(FeetPerMinute, SI::MetersPerSecond, 5'080, 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(FeetPerSecond, SI::MetersPerSecond, 304'800, 1'000'000)

	/* Pressure */;

	GENERATE_NSTD_FROM_FRACTION(InchesOfMercury, SI::Pascals, 3'386'389, 1'000)
	GENERATE_NSTD_FROM_FRACTION(Bars, SI::Pascals, 100'000, 1)
	GENERATE_NSTD_FROM_FRACTION(MilliBars, SI::Pascals, 100, 1)
	GENERATE_NSTD_FROM_FRACTION(Atmospheres, SI::Pascals, 101'325, 1)
	GENERATE_NSTD_FROM_FRACTION(PoundsPerSquareInch, SI::Pascals, 4'448'221'615'260'500, 645'160'000'000)

	/* Mass & Force */;

	GENERATE_NSTD_FROM_FRACTION(Ounces, SI::KiloGrams, 28'349'523'125, 1'000'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(PoundsMass, SI::KiloGrams, 453'592'370, 1'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(ShortTons, SI::KiloGrams, 907'184'740, 1'000'000)
	GENERATE_NSTD_FROM_FRACTION(Slug, SI::KiloGrams, 4'448'221'615'260'500, 304'800'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(PoundsForce, SI::Newtons, 4'448'221'615'260'500, 1'000'000'000'000'000)

	/* Volume */;

	GENERATE_NSTD_FROM_FRACTION(UsFluidOunces, SI::MetersCube, 3'785'411'784, 128'000'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(UsPints, SI::MetersCube, 3'785'411'784, 8'000'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(UsQuarts, SI::MetersCube, 3'785'411'784, 4'000'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(UsGallons, SI::MetersCube, 3'785'411'784, 1'000'000'000'000)

	GENERATE_NSTD_FROM_FRACTION(ImperialFluidOunces, SI::MetersCube, 4'546'090, 160'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(ImperialPints, SI::MetersCube, 4'546'090, 8'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(ImperialQuarts, SI::MetersCube, 4'546'090, 4'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(ImperialGallons, SI::MetersCube, 4'546'090, 1'000'000'000)

	/* Power, Energy, Torque */;

	GENERATE_NSTD_FROM_FRACTION(Horsepower, SI::Watts, 550 * 13'558'179'483'314'004, 10'000'000'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(FootPoundsForce, SI::Joules, 13'558'179'483'314'004, 10'000'000'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(Btu, SI::Joules, 1055, 1)
	GENERATE_NSTD_FROM_FRACTION(PoundForceFeet, SI::NewtonMeters, 13'558'179'483'314'004, 10'000'000'000'000'000)
	GENERATE_NSTD_FROM_FRACTION(PoundForceInches, SI::NewtonMeters, 13'558'179'483'314'004, 1'200'000'000'000'000'000)

	/* Temperature */;

	GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION(DegreesCelsius, SI::Kelvins, 1, 1, 273'150, 1'000)
	GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION(DegreesFahrenheit, SI::Kelvins, 5, 9, 9 * 273'150 - 4 * 40'000, 9 * 1'000)

#undef GENERATE_NSTD_FROM_DOUBLE
#undef GENERATE_NSTD_FROM_FRACTION
#undef GENERATE_OFFSETTABLE_NSTD_FROM_FRACTION
}