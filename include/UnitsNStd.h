#pragma once

#include "UnitsSI.h"

namespace Physics::Units::NStd
{
	template<i64 Numerator = 0, i64 Denominator = 1>
		requires(Denominator > 0)
	class FractionNorm
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
		static constexpr bool isUnchanged() {
			return greatestCommonDivisor<Numerator, Denominator>() == 1;
		}
	};
	class Fraction
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
			// TODO: verify if it is correct in all cases
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

	class _Base
	{};

	template<typename Type, template<typename> class StandardUnit, i64 ScaleNumerator, i64 ScaleDenominator,
			 i64 OffsetNumerator, i64 OffsetDenominator, typename AccuracyType = f64>
		requires(!std::is_base_of_v<_Base, StandardUnit<u8>> && std::is_arithmetic_v<Type> &&
				 ScaleDenominator > 0 && OffsetDenominator > 0 && ScaleNumerator > 0 && StandardUnit<i8>::hasNoPrefix() &&
				 !(OffsetNumerator == 0 && ScaleNumerator == 1 && ScaleDenominator == 1))
	class GenerativeUnit : _Base
	{
	private:
		static constexpr bool hasNoOffset() {
			return OffsetNumerator == 0;
		}
		static constexpr bool hasIdentScale() {
			return ScaleNumerator == 1 && ScaleDenominator == 1;
		}
		static constexpr bool hasSameScale(const i64 _ScaleNumerator, const i64 _ScaleDenominator) {
			return ScaleNumerator == _ScaleNumerator && ScaleDenominator == _ScaleDenominator;
		}

		template<typename _Type = Type, template<typename> class _StandardUnit = StandardUnit>
		using Self = GenerativeUnit<_Type, _StandardUnit, ScaleNumerator, ScaleDenominator,
									OffsetNumerator, OffsetDenominator, AccuracyType>;

		template<typename _Type = Type, template<typename> class _StandardUnit = StandardUnit, i64 _ScaleNumerator = ScaleNumerator,
				 i64 _ScaleDenominator = ScaleDenominator, i64 _OffsetNumerator = OffsetNumerator,
				 i64 _OffsetDenominator = OffsetDenominator, typename _AccuracyType = AccuracyType>
		using AdaptiveSelf = std::conditional_t<_ScaleNumerator == 1 && _ScaleDenominator == 1 && _OffsetNumerator == 0, _StandardUnit<_Type>,
												GenerativeUnit<_Type, _StandardUnit, _ScaleNumerator, _ScaleDenominator, _OffsetNumerator, _OffsetDenominator, _AccuracyType>>;

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
			requires(std::is_base_of_v<_Base, SiblingUnit>)
		static consteval bool hasSameStdUnitBase() {
			return std::is_same_v<decltype(SiblingUnit().toStandardUnit() / StandardUnit<decltype(SiblingUnit().toRaw())>()),
								  SI::Scale<decltype(SiblingUnit().toRaw())>>;
		}

	public:
		// Constructors

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
		constexpr GenerativeUnit(const SiblingUnit value) : mData(fromOtherNStd(value)) {}

		// Assignment operators

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

		// Comparison operators

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

		// Logical operators

		constexpr bool operator!() const {
			return !mData;
		}

		// Arithmetic operators

		constexpr AdaptiveSelf<decltype(-Type()), StandardUnit, ScaleNumerator, ScaleDenominator, -OffsetNumerator>
			operator-() const {
			return -mData;
		}
		template<typename SiblingUnit = Self<>>
			requires(std::is_base_of_v<_Base, SiblingUnit> && SiblingUnit::hasSameScale(ScaleNumerator, ScaleDenominator) &&
					 hasSameStdUnitBase<SiblingUnit>())
		constexpr auto operator+(const SiblingUnit value) const {
			using RawType = decltype(SiblingUnit().toRaw() * Type());
			using FrNm = FractionNorm<OffsetNumerator * SiblingUnit::getOffsetDenominator() +
										  OffsetDenominator * SiblingUnit::getOffsetNumerator(),
									  OffsetDenominator * SiblingUnit::getOffsetDenominator()>;
			constexpr i64 OffNum = FrNm::Num;
			constexpr i64 OffDenom = FrNm::Denom;
			using AccType = decltype(AccuracyType() * SiblingUnit::SCALE);
			using NewUnit = AdaptiveSelf<RawType, StandardUnit, ScaleNumerator, ScaleDenominator,
										 OffNum, OffDenom, AccType>;
			return NewUnit{mData + value.toRaw()};
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		constexpr Self<decltype(Type() * _Type())> operator+(const StandardUnit<_Type> value) const {
			return mData + value.toRaw();
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		friend constexpr Self<decltype(Type() * _Type())> operator+(const StandardUnit<_Type> &value, const Self<> &self) {
			return value.toRaw() + self.toRaw();
		}
		template<typename SiblingUnit = Self<>>
			requires(std::is_base_of_v<_Base, SiblingUnit> && SiblingUnit::hasSameScale(ScaleNumerator, ScaleDenominator) &&
					 hasSameStdUnitBase<SiblingUnit>())
		constexpr auto operator-(const SiblingUnit value) const {
			using RawType = decltype(SiblingUnit().toRaw() * Type());
			using FrNm = FractionNorm<OffsetNumerator * SiblingUnit::getOffsetDenominator() -
										  OffsetDenominator * SiblingUnit::getOffsetNumerator(),
									  OffsetDenominator * SiblingUnit::getOffsetDenominator()>;
			constexpr i64 OffNum = FrNm::Num;
			constexpr i64 OffDenom = FrNm::Denom;
			using AccType = decltype(AccuracyType() * SiblingUnit::SCALE);
			using NewUnit = AdaptiveSelf<RawType, StandardUnit, ScaleNumerator, ScaleDenominator,
										 OffNum, OffDenom, AccType>;
			return NewUnit{mData - value.toRaw()};
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		constexpr Self<decltype(Type() * _Type())> operator-(const StandardUnit<_Type> value) const {
			return mData - value.toRaw();
		}
		template<typename _Type = Type>
			requires(hasIdentScale())
		friend constexpr auto operator-(const StandardUnit<_Type> &value, const Self<> &self) {
			using NewUnit = decltype(-Self<decltype(Type() * _Type())>());
			return NewUnit{value.toRaw() - self.toRaw()};
		}
		template<typename _Type, template<typename> class _StandardUnit, i64 _ScaleNumerator, i64 _ScaleDenominator,
				 i64 _OffsetNumerator, i64 _OffsetDenominator, typename _AccuracyType = f64>
			requires(hasNoOffset() && _OffsetNumerator == 0)
		constexpr auto operator*(const GenerativeUnit<_Type, _StandardUnit, _ScaleNumerator, _ScaleDenominator,
													  _OffsetNumerator, _OffsetDenominator, _AccuracyType>
									 value) const {
			using RawType = decltype(_Type() * Type());
			using ComposeType = decltype(StandardUnit<Type>() * _StandardUnit<_Type>());
			using AccType = decltype(AccuracyType() * _AccuracyType());
			using FrNm = FractionNorm<ScaleNumerator * _ScaleNumerator, ScaleDenominator * _ScaleDenominator>;
			constexpr i64 ScNum = FrNm::Num;
			constexpr i64 ScDenom = FrNm::Denom;
			using NewUnit = AdaptiveSelf<RawType, Decomposer<ComposeType>::template OuterType, ScNum, ScDenom, 0, 1, AccType>;
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
		template<typename _Type, template<typename> class _StandardUnit, i64 _ScaleNumerator, i64 _ScaleDenominator,
				 i64 _OffsetNumerator, i64 _OffsetDenominator, typename _AccuracyType = f64>
			requires(hasNoOffset() && _OffsetNumerator == 0)
		constexpr auto operator/(const GenerativeUnit<_Type, _StandardUnit, _ScaleNumerator, _ScaleDenominator,
													  _OffsetNumerator, _OffsetDenominator, _AccuracyType>
									 value) const {
			using RawType = decltype(_Type() * Type());
			using ComposeType = decltype(StandardUnit<Type>() / _StandardUnit<_Type>());
			using AccType = decltype(AccuracyType() * _AccuracyType());
			constexpr i64 ScNum = FractionNorm<ScaleNumerator * _ScaleDenominator, ScaleDenominator * _ScaleNumerator>::Num;
			constexpr i64 ScDenom = FractionNorm<ScaleNumerator * _ScaleDenominator, ScaleDenominator * _ScaleNumerator>::Denom;
			using NewUnit = AdaptiveSelf<RawType, Decomposer<ComposeType>::template OuterType, ScNum, ScDenom, 0, 1, AccType>;
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

		// Other methods

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
			return ScaleNumerator;
		}
		static consteval auto getScaleDenominator() {
			return ScaleDenominator;
		}
		static consteval auto getOffsetNumerator() {
			return OffsetNumerator;
		}
		static consteval auto getOffsetDenominator() {
			return OffsetDenominator;
		}

	private:
		template<typename _Type = Type>
		static constexpr auto fromStandardUnit(const StandardUnit<_Type> value) {
			if constexpr(SCALE == 1)
				return value.toRaw() - OFFSET;
			else if constexpr(OFFSET == 0)
				return value.toRaw() / SCALE;
			else
				return (value.toRaw() - OFFSET) / SCALE;
		}
		template<typename SiblingUnit = Self<>>
			requires(hasSameStdUnitBase<SiblingUnit>())
		static constexpr auto fromOtherNStd(const SiblingUnit value) {
			// TODO: optimize calculation
			return fromStandardUnit(value.toStandardUnit());
		}

		static constexpr AccuracyType SCALE = AccuracyType(ScaleNumerator) / ScaleDenominator;
		static constexpr AccuracyType OFFSET = AccuracyType(OffsetNumerator) / OffsetDenominator;
		Type mData;
	};

	template<typename Type, template<typename> class StdU, i64 ScNum, i64 ScDenom, i64 OffNum = 0, i64 OffDenom = 1>
		requires(ScNum > 0 && ScDenom > 0 && OffDenom > 0)
	using Simplifier = GenerativeUnit<Type, StdU, FractionNorm<ScNum, ScDenom>::Num, FractionNorm<ScNum, ScDenom>::Denom,
									  FractionNorm<OffNum, OffDenom>::Num, FractionNorm<OffNum, OffDenom>::Denom>;

#define GENERATE_NSTD_FROM_DOUBLE(Name, StdType, Scale)                                          \
	static_assert(Fraction::isConvertible(Scale), "Value '" #Scale "' out of supported range!"); \
	template<typename T = f64>                                                                   \
	using Name = Simplifier<T, StdType, Fraction::numerator(Scale), Fraction::denominator(Scale)>;

#define GENERATE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom) \
	template<typename T = f64>                                     \
	using Name = Simplifier<T, StdType, ScNum, ScDenom>;

#define GENERATE_OFFSETABLE_NSTD_FROM_FRACTION(Name, StdType, ScNum, ScDenom, OffNum, OffDenom) \
	template<typename T = f64>                                                                  \
	using Name = Simplifier<T, StdType, ScNum, ScDenom, OffNum, OffDenom>;

	/* --- GENERATE NEEDED UNITS HERE --- */

	GENERATE_NSTD_FROM_FRACTION(Percent, SI::Scale, 1, 100)
	GENERATE_NSTD_FROM_FRACTION(Permille, SI::Scale, 1, 1'000)
	GENERATE_NSTD_FROM_FRACTION(Inch, SI::Meters, 25'400, 1'000'000)
	constexpr f64 DEG2RAD = M_PI / 180;
	GENERATE_NSTD_FROM_DOUBLE(Degree, SI::Radians, DEG2RAD)
	GENERATE_OFFSETABLE_NSTD_FROM_FRACTION(DegreeCelsius, SI::Kelvins, 1, 1, 273'150, 1'000)
	GENERATE_OFFSETABLE_NSTD_FROM_FRACTION(DegreeFahrenheit, SI::Kelvins, 5, 9, 9 * 273'150 - 4 * 40'000, 9 * 1'000)
}