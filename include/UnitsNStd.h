#pragma once

#include "UnitsSI.h"

namespace Physics::Units::NStd
{
	namespace _
	{
		template<i64 Numerator = 0, i64 Denominator = 1>
			requires(Denominator > 0)
		class Fraction
		{
		private:
			template<i64 A, i64 B>
				requires(B > 0)
			static consteval i64 greatestCommonDivisor() {
				i64 a = A < 0 ? -A : A;
				i64 b = B;
				while(b != 0) {
					i64 temp = b;
					b = a % b;
					a = temp;
				}
				return a;
			}

		public:
			static constexpr i64 Num = Numerator / greatestCommonDivisor<Numerator, Denominator>();
			static constexpr i64 Denom = Denominator / greatestCommonDivisor<Numerator, Denominator>();
			static consteval bool isNormalized() {
				return greatestCommonDivisor<Numerator, Denominator>() == 1;
			}
			static consteval bool isIdentity() {
				return Numerator == 1 && Denominator == 1;
			}
			static consteval bool isZero() {
				return Numerator == 0;
			}
			static consteval bool isPositive() {
				return Numerator > 0;
			}
			template<typename T = f64>
				requires(std::is_arithmetic_v<T>)
			static consteval T toDecimal() {
				return static_cast<T>(Num) / static_cast<T>(Denom);
			}

			using Norm = Fraction<Num, Denom>;
			using Opposite = Fraction<-Num, Denom>;

			template<class _Fraction>
			using Sum = Fraction<Num * _Fraction::Denom + Denom * _Fraction::Num, Denom * _Fraction::Denom>::Norm;

			template<class _Fraction>
			using Diff = Fraction<Num * _Fraction::Denom - Denom * _Fraction::Num, Denom * _Fraction::Denom>::Norm;

			template<class _Fraction>
			using Product = Fraction<Num * _Fraction::Num, Denom * _Fraction::Denom>::Norm;

			template<class _Fraction>
			using Quotient = Fraction<Num * _Fraction::Denom, Denom * _Fraction::Num>::Norm;
		};

		class Base
		{};
	}

	template<typename Type, template<typename> class StandardUnit, class Scale, class Offset, typename AccuracyType = f64>
		requires(!std::is_base_of_v<_::Base, StandardUnit<u8>> && std::is_arithmetic_v<Type> && StandardUnit<i8>::hasNoPrefix() &&
				 Scale::isNormalized() && Offset::isNormalized() && Scale::isPositive() && !(Offset::isZero() && Scale::isIdentity()))
	class GenerativeUnit : _::Base
	{
	private:
		static consteval bool hasNoOffset() {
			return Offset::isZero();
		}
		static consteval bool hasIdentScale() {
			return Scale::isIdentity();
		}
		template<class _Fraction>
		static consteval bool hasSameScale() {
			return std::is_same_v<Scale, _Fraction>;
		}

		template<typename _Type = Type, template<typename> class _StandardUnit = StandardUnit>
		using Self = GenerativeUnit<_Type, _StandardUnit, Scale, Offset, AccuracyType>;

		template<typename _Type = Type, template<typename> class _StandardUnit = StandardUnit, class _Scale = Scale,
				 class _Offset = Offset, typename _AccuracyType = AccuracyType>
		using AdaptiveSelf = std::conditional_t<_Scale::isIdentity() && _Offset::isZero(), _StandardUnit<_Type>,
												GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType>>;

		template<class Base>
		struct Decomposer {
			template<typename T>
			using OuterType = decltype(Base().template cast<T>());
			using InnerType = decltype(Base().toRaw());
		};
		template<class Complex>
		struct StdToSelf {
			using NewUnit = Self<decltype(Complex().toRaw()), Decomposer<Complex>::template OuterType>;
		};
		template<typename SiblingUnit>
			requires(std::is_base_of_v<_::Base, SiblingUnit>)
		static consteval bool hasSameStdUnitBase() {
			return std::is_same_v<decltype(SiblingUnit().toStandardUnit() / StandardUnit<decltype(SiblingUnit().toRaw())>()),
								  SI::Scale<decltype(SiblingUnit().toRaw())>>;
		}

	public:
		/* Constructors */;

		constexpr GenerativeUnit(const Self<> &value) : mData(value.mData) {}
		constexpr GenerativeUnit(const Self<> &&value) : mData(value.mData) {}
		constexpr GenerativeUnit() : mData(0) {}
		constexpr GenerativeUnit(const Type data) : mData(data) {}
		template<typename _Type = Type>
		constexpr GenerativeUnit(const StandardUnit<_Type> value) : mData(fromStandardUnit(value)) {}
		template<typename _Type = Type>
		constexpr GenerativeUnit(const Self<_Type> value) : mData(value.toRaw()) {}
		template<typename SiblingUnit = Self<>>
			requires(hasSameStdUnitBase<SiblingUnit>())
		constexpr GenerativeUnit(const SiblingUnit value) : mData(fromOtherNStd(value)){};

		/* Assignment operators */;

		constexpr Self<> &operator=(const Self<> &value) {
			mData = value.mData;
			return *this;
		}
		constexpr Self<> &operator=(const Self<> &&value) {
			mData = value.mData;
			return *this;
		}
		template<typename _Type = Type>
		constexpr Self<> &operator=(const StandardUnit<_Type> value) {
			mData = fromStandardUnit(value);
			return *this;
		}
		template<typename _Type = Type>
		constexpr Self<> &operator=(const Self<_Type> value) {
			mData = value.toRaw();
			return *this;
		}
		template<typename SiblingUnit = Self<>>
			requires(hasSameStdUnitBase<SiblingUnit>())
		constexpr Self<> &operator=(const SiblingUnit value) {
			mData = fromOtherNStd(value);
			return *this;
		}

		template<typename _Type = Type>
		constexpr Self<> &operator+=(const StandardUnit<_Type> value) {
			mData += value.toRaw() / SCALE;
			return *this;
		}
		template<typename _Type = Type>
			requires(hasNoOffset())
		constexpr Self<> &operator+=(const Self<_Type> value) {
			mData += value.toRaw();
			return *this;
		}
		template<typename SiblingUnit = Self<>>
			requires(hasSameStdUnitBase<SiblingUnit>() && SiblingUnit::hasNoOffset())
		constexpr Self<> &operator+=(const SiblingUnit value) {
			mData += fromOtherNStd(value);
			return *this;
		}
		template<typename _Type = Type>
		constexpr Self<> &operator-=(const StandardUnit<_Type> value) {
			mData -= value.toRaw() / SCALE;
			return *this;
		}
		template<typename _Type = Type>
			requires(hasNoOffset())
		constexpr Self<> &operator-=(const Self<_Type> value) {
			mData -= value.toRaw();
			return *this;
		}
		template<typename SiblingUnit = Self<>>
			requires(hasSameStdUnitBase<SiblingUnit>() && SiblingUnit::hasNoOffset())
		constexpr Self<> &operator-=(const SiblingUnit value) {
			mData -= fromOtherNStd(value);
			return *this;
		}
		template<typename _Type = Type>
			requires(hasNoOffset())
		constexpr Self<> &operator*=(const SI::Scale<_Type> value) {
			mData *= value.toRaw();
			return *this;
		}
		template<typename _Type = Type>
			requires(hasNoOffset())
		constexpr Self<> &operator/=(const SI::Scale<_Type> value) {
			mData /= value.toRaw();
			return *this;
		}

		/* Comparison operators */;

		template<typename _Type = Type>
		constexpr bool operator<(const Self<_Type> value) const {
			return mData < value.toRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator>(const Self<_Type> value) const {
			return mData > value.toRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator<=(const Self<_Type> value) const {
			return mData <= value.toRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator>=(const Self<_Type> value) const {
			return mData >= value.toRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator==(const Self<_Type> value) const {
			return mData == value.toRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator!=(const Self<_Type> value) const {
			return mData != value.toRaw();
		}
		template<typename _Type = Type>
		constexpr auto operator<=>(const Self<_Type> value) const {
			return mData <=> value.toRaw();
		}

		/* Logical operators */;

		constexpr bool operator!() const {
			return !mData;
		}

		/* Arithmetic operators */;

		constexpr AdaptiveSelf<decltype(-Type()), StandardUnit, Scale, typename Offset::Opposite>
			operator-() const {
			return -mData;
		}
		template<typename SiblingUnit = Self<>>
			requires(std::is_base_of_v<_::Base, SiblingUnit> && SiblingUnit::template hasSameScale<Scale>() && hasSameStdUnitBase<SiblingUnit>())
		constexpr auto operator+(const SiblingUnit value) const {
			using RawType = decltype(Type() + SiblingUnit().toRaw());
			using ResultOffset = Offset::template Sum<_::Fraction<SiblingUnit::getOffsetNumerator(), SiblingUnit::getOffsetDenominator()>>;
			using AccType = decltype(AccuracyType() * SiblingUnit::SCALE);
			using NewUnit = AdaptiveSelf<RawType, StandardUnit, Scale, ResultOffset, AccType>;
			return NewUnit{mData + value.toRaw()};
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		constexpr Self<decltype(Type() + _Type())> operator+(const StandardUnit<_Type> value) const {
			return mData + value.toRaw();
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		friend constexpr Self<decltype(_Type() + Type())> operator+(const StandardUnit<_Type> &value, const Self<> &self) {
			return value.toRaw() + self.toRaw();
		}
		template<typename SiblingUnit = Self<>>
			requires(std::is_base_of_v<_::Base, SiblingUnit> && SiblingUnit::template hasSameScale<Scale>() && hasSameStdUnitBase<SiblingUnit>())
		constexpr auto operator-(const SiblingUnit value) const {
			using RawType = decltype(Type() - SiblingUnit().toRaw());
			using ResultOffset = Offset::template Diff<_::Fraction<SiblingUnit::getOffsetNumerator(), SiblingUnit::getOffsetDenominator()>>;
			using AccType = decltype(AccuracyType() * SiblingUnit::SCALE);
			using NewUnit = AdaptiveSelf<RawType, StandardUnit, Scale, ResultOffset, AccType>;
			return NewUnit{mData - value.toRaw()};
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		constexpr Self<decltype(Type() - _Type())> operator-(const StandardUnit<_Type> value) const {
			return mData - value.toRaw();
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		friend constexpr auto operator-(const StandardUnit<_Type> &value, const Self<> &self) {
			using NewUnit = decltype(-Self<decltype(_Type() - Type())>());
			return NewUnit{value.toRaw() - self.toRaw()};
		}
		template<typename _Type, template<typename> class _StandardUnit, class _Scale, class _Offset, typename _AccuracyType = f64>
			requires(hasNoOffset() && _Offset::isZero())
		constexpr auto operator*(const GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType> value) const {
			using RawType = decltype(Type() * _Type());
			using ComposeType = decltype(StandardUnit<Type>() * _StandardUnit<_Type>());
			using AccType = decltype(AccuracyType() * _AccuracyType());
			using ResultScale = Scale::template Product<_Scale>;
			using NewUnit = AdaptiveSelf<RawType, Decomposer<ComposeType>::template OuterType, ResultScale, _::Fraction<0, 1>, AccType>;
			return NewUnit{mData * value.toRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(hasNoOffset())
		constexpr typename StdToSelf<decltype(StandardUnit<Type>() * _StdUnitT())>::NewUnit
			operator*(const _StdUnitT value) const {
			return mData * value.toRaw();
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(hasNoOffset())
		friend constexpr typename StdToSelf<decltype(_StdUnitT() * StandardUnit<Type>())>::NewUnit
			operator*(const _StdUnitT &value, const Self<> &self) {
			return value.toRaw() * self.toRaw();
		}
		template<typename _Type, template<typename> class _StandardUnit, class _Scale, class _Offset, typename _AccuracyType = f64>
			requires(hasNoOffset() && _Offset::isZero())
		constexpr auto operator/(const GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType> value) const {
			using RawType = decltype(Type() / _Type());
			using ComposeType = decltype(StandardUnit<Type>() / _StandardUnit<_Type>());
			using AccType = decltype(AccuracyType() / _AccuracyType());
			using ResultScale = Scale::template Quotient<_Scale>;
			using NewUnit = AdaptiveSelf<RawType, Decomposer<ComposeType>::template OuterType, ResultScale, _::Fraction<0, 1>, AccType>;
			return NewUnit{mData / value.toRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(hasNoOffset())
		constexpr typename StdToSelf<decltype(StandardUnit<Type>() / _StdUnitT())>::NewUnit
			operator/(const _StdUnitT value) const {
			return mData / value.toRaw();
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(hasNoOffset())
		friend constexpr typename StdToSelf<decltype(_StdUnitT() / StandardUnit<Type>())>::NewUnit
			operator/(const _StdUnitT &value, const Self<> &self) {
			return value.toRaw() / self.toRaw();
		}

		friend std::ostream &operator<<(std::ostream &os, const Self<> &obj) {
			os << obj.mData;
			return os;
		}

		/* Other methods */;

		template<typename _Type = Type>
		constexpr Self<_Type> cast() const {
			return mData;
		}

		constexpr Type toRaw() const {
			return mData;
		}

		template<typename _Type = Type>
		constexpr StandardUnit<_Type> toStandardUnit() const {
			if constexpr(SCALE == 1)
				return mData + OFFSET;
			else if constexpr(OFFSET == 0)
				return mData * SCALE;
			else
				return mData * SCALE + OFFSET;
		}

		static consteval auto getScaleNumerator() {
			return Scale::Num;
		}
		static consteval auto getScaleDenominator() {
			return Scale::Denom;
		}
		static consteval auto getOffsetNumerator() {
			return Offset::Num;
		}
		static consteval auto getOffsetDenominator() {
			return Offset::Denom;
		}

	private:
		template<typename _Type = Type>
		static constexpr Type fromStandardUnit(const StandardUnit<_Type> value) {
			if constexpr(SCALE == 1)
				return value.toRaw() - OFFSET;
			else if constexpr(OFFSET == 0)
				return value.toRaw() / SCALE;
			else
				return (value.toRaw() - OFFSET) / SCALE;
		}
		template<typename SiblingUnit = Self<>>
			requires(hasSameStdUnitBase<SiblingUnit>())
		static constexpr Type fromOtherNStd(const SiblingUnit value) {
			if constexpr(std::is_same_v<Self<decltype(SiblingUnit().toRaw())>, SiblingUnit>) {
				return value;
			} else {
				constexpr auto compScale = SiblingUnit::SCALE * SCALE;
				constexpr auto compOffset = (SiblingUnit::OFFSET - OFFSET) * SCALE;

				return value.toRaw() * compScale + compOffset;
			}
		}

		/* TODO: Consider how to handle `AccuracyType` */;
		inline static constexpr AccuracyType SCALE = Scale::template toDecimal<AccuracyType>();
		inline static constexpr AccuracyType OFFSET = Offset::template toDecimal<AccuracyType>();
		Type mData;
	};

	namespace _
	{
		class FractionParser
		{
		private:
			static consteval i64 exp2(i64 value) noexcept {
				if(value < 0) return 0;
				i64 result = 1;
				while(--value >= 0) {
					result *= 2L;
				}
				return result;
			}
			static consteval i64 flooredLog2(f64 value) noexcept {
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
			static consteval i64 calcExponent(const f64 value) noexcept {
				return flooredLog2(std::abs(value));
			}

		public:
			static consteval i64 denominator(const f64 value) noexcept {
				const i64 exp = calcExponent(value);
				return exp2(exp < 0L ? 62L : 62L - exp);
			}
			static consteval i64 numerator(const f64 value) noexcept {
				return value * denominator(value);
			}
			static consteval bool isConvertible(const f64 value) noexcept {
				const i64 exp = calcExponent(value);
				return exp > -11 && exp < 63;
			}
		};

		template<typename Type, template<typename> class StdU, i64 ScNum, i64 ScDenom, i64 OffNum = 0, i64 OffDenom = 1>
		using Simplifier = GenerativeUnit<Type, StdU, typename Fraction<ScNum, ScDenom>::Norm, typename Fraction<OffNum, OffDenom>::Norm>;

		inline static constexpr f64 DEG2RAD = M_PI / 180;
	}

#define GENERATE_NSTD_FROM_DOUBLE(Name, StdType, Scale)                                                   \
	static_assert(_::FractionParser::isConvertible(Scale), "Value '" #Scale "' out of supported range!"); \
	template<typename T = f64>                                                                            \
	using Name = _::Simplifier<T, StdType, _::FractionParser::numerator(Scale), _::FractionParser::denominator(Scale)>;

#define GENERATE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom) \
	template<typename T = f64>                                     \
	using Name = _::Simplifier<T, StdType, ScNum, ScDenom>;

#define GENERATE_OFFSETABLE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom, OffNum, OffDenom) \
	template<typename T = f64>                                                                  \
	using Name = _::Simplifier<T, StdType, ScNum, ScDenom, OffNum, OffDenom>;

	/* --- GENERATE NEEDED UNITS HERE --- */

	GENERATE_NSTD_FROM_FRACTION(Percent, SI::Scale, 1, 100)
	GENERATE_NSTD_FROM_FRACTION(Permille, SI::Scale, 1, 1'000)
	GENERATE_NSTD_FROM_FRACTION(Inch, SI::Meters, 25'400, 1'000'000)
	GENERATE_NSTD_FROM_DOUBLE(Degree, SI::Radians, _::DEG2RAD)
	GENERATE_OFFSETABLE_NSTD_FROM_FRACTION(DegreeCelsius, SI::Kelvins, 1, 1, 273'150, 1'000)
	GENERATE_OFFSETABLE_NSTD_FROM_FRACTION(DegreeFahrenheit, SI::Kelvins, 5, 9, 9 * 273'150 - 4 * 40'000, 9 * 1'000)
}