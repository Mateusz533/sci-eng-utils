#pragma once
//
#include "UnitsSI.hpp"

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
			static consteval i64 GreatestCommonDivisor() {
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
			static constexpr i64 Num = Numerator / GreatestCommonDivisor<Numerator, Denominator>();
			static constexpr i64 Denom = Denominator / GreatestCommonDivisor<Numerator, Denominator>();

			static consteval bool IsNormalized() {
				return GreatestCommonDivisor<Numerator, Denominator>() == 1;
			}
			static consteval bool IsIdentity() {
				return Numerator == 1 && Denominator == 1;
			}
			static consteval bool IsZero() {
				return Numerator == 0;
			}
			static consteval bool IsPositive() {
				return Numerator > 0;
			}
			template<typename T = f64>
				requires(std::is_arithmetic_v<T>)
			static consteval T ToDecimal() {
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
		requires(!std::is_base_of_v<_::Base, StandardUnit<u8>> && std::is_arithmetic_v<Type> && StandardUnit<i8>::HasNoPrefix() &&
				 Scale::IsNormalized() && Offset::IsNormalized() && Scale::IsPositive() && !(Offset::IsZero() && Scale::IsIdentity()))
	class GenerativeUnit : _::Base
	{
	private:
		static consteval bool HasNoOffset() {
			return Offset::IsZero();
		}
		static consteval bool HasIdentScale() {
			return Scale::IsIdentity();
		}
		template<class _Fraction>
		static consteval bool HasSameScale() {
			return std::is_same_v<Scale, _Fraction>;
		}

		template<typename _Type = Type, template<typename> class _StandardUnit = StandardUnit>
		using Self = GenerativeUnit<_Type, _StandardUnit, Scale, Offset, AccuracyType>;

		template<typename _Type = Type, template<typename> class _StandardUnit = StandardUnit, class _Scale = Scale,
				 class _Offset = Offset, typename _AccuracyType = AccuracyType>
		using AdaptiveSelf = std::conditional_t<_Scale::IsIdentity() && _Offset::IsZero(), _StandardUnit<_Type>,
												GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType>>;

		template<class Base>
		struct Decomposer {
			template<typename T>
			using OuterType = decltype(Base().template Cast<T>());
			using InnerType = decltype(Base().ToRaw());
		};
		template<class Complex>
		struct StdToSelf {
			using NewUnit = Self<decltype(Complex().ToRaw()), Decomposer<Complex>::template OuterType>;
		};
		template<typename SiblingUnit>
			requires(std::is_base_of_v<_::Base, SiblingUnit>)
		static consteval bool HasSameStdUnitBase() {
			return std::is_same_v<decltype(SiblingUnit().toStandardUnit() / StandardUnit<decltype(SiblingUnit().ToRaw())>()),
								  SI::Scale<decltype(SiblingUnit().ToRaw())>>;
		}

	public:
		/* Constructors */;

		constexpr GenerativeUnit(const Self<> &value) : mData(value.mData) {}
		constexpr GenerativeUnit(const Self<> &&value) : mData(value.mData) {}
		constexpr GenerativeUnit() : mData(0) {}
		constexpr GenerativeUnit(const Type data) : mData(data) {}
		template<typename _Type = Type>
		constexpr GenerativeUnit(const StandardUnit<_Type> value) : mData(FromStandardUnit(value)) {}
		template<typename _Type = Type>
		constexpr GenerativeUnit(const Self<_Type> value) : mData(value.ToRaw()) {}
		template<typename SiblingUnit = Self<>>
			requires(HasSameStdUnitBase<SiblingUnit>())
		constexpr GenerativeUnit(const SiblingUnit value) : mData(FromOtherNStd(value)){};

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
			mData = FromStandardUnit(value);
			return *this;
		}
		template<typename _Type = Type>
		constexpr Self<> &operator=(const Self<_Type> value) {
			mData = value.ToRaw();
			return *this;
		}
		template<typename SiblingUnit = Self<>>
			requires(HasSameStdUnitBase<SiblingUnit>())
		constexpr Self<> &operator=(const SiblingUnit value) {
			mData = FromOtherNStd(value);
			return *this;
		}

		template<typename _Type = Type>
		constexpr Self<> &operator+=(const StandardUnit<_Type> value) {
			mData += value.ToRaw() / SCALE;
			return *this;
		}
		template<typename _Type = Type>
			requires(HasNoOffset())
		constexpr Self<> &operator+=(const Self<_Type> value) {
			mData += value.ToRaw();
			return *this;
		}
		template<typename SiblingUnit = Self<>>
			requires(HasSameStdUnitBase<SiblingUnit>() && SiblingUnit::HasNoOffset())
		constexpr Self<> &operator+=(const SiblingUnit value) {
			mData += FromOtherNStd(value);
			return *this;
		}
		template<typename _Type = Type>
		constexpr Self<> &operator-=(const StandardUnit<_Type> value) {
			mData -= value.ToRaw() / SCALE;
			return *this;
		}
		template<typename _Type = Type>
			requires(HasNoOffset())
		constexpr Self<> &operator-=(const Self<_Type> value) {
			mData -= value.ToRaw();
			return *this;
		}
		template<typename SiblingUnit = Self<>>
			requires(HasSameStdUnitBase<SiblingUnit>() && SiblingUnit::HasNoOffset())
		constexpr Self<> &operator-=(const SiblingUnit value) {
			mData -= FromOtherNStd(value);
			return *this;
		}
		template<typename _Type = Type>
			requires(HasNoOffset())
		constexpr Self<> &operator*=(const SI::Scale<_Type> value) {
			mData *= value.ToRaw();
			return *this;
		}
		template<typename _Type = Type>
			requires(HasNoOffset())
		constexpr Self<> &operator/=(const SI::Scale<_Type> value) {
			mData /= value.ToRaw();
			return *this;
		}

		/* Comparison operators */;

		template<typename _Type = Type>
		constexpr bool operator<(const Self<_Type> value) const {
			return mData < value.ToRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator>(const Self<_Type> value) const {
			return mData > value.ToRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator<=(const Self<_Type> value) const {
			return mData <= value.ToRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator>=(const Self<_Type> value) const {
			return mData >= value.ToRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator==(const Self<_Type> value) const {
			return mData == value.ToRaw();
		}
		template<typename _Type = Type>
		constexpr bool operator!=(const Self<_Type> value) const {
			return mData != value.ToRaw();
		}
		template<typename _Type = Type>
		constexpr auto operator<=>(const Self<_Type> value) const {
			return mData <=> value.ToRaw();
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
			requires(std::is_base_of_v<_::Base, SiblingUnit> && SiblingUnit::template hasSameScale<Scale>() && HasSameStdUnitBase<SiblingUnit>())
		constexpr auto operator+(const SiblingUnit value) const {
			using RawType = decltype(Type() + SiblingUnit().ToRaw());
			using ResultOffset = Offset::template Sum<_::Fraction<SiblingUnit::GetOffsetNumerator(), SiblingUnit::GetOffsetDenominator()>>;
			using AccType = decltype(AccuracyType() * SiblingUnit::SCALE);
			using NewUnit = AdaptiveSelf<RawType, StandardUnit, Scale, ResultOffset, AccType>;
			return NewUnit{mData + value.ToRaw()};
		}
		template<typename _Type = Type>
			requires(HasIdentScale())
		constexpr Self<decltype(Type() + _Type())> operator+(const StandardUnit<_Type> value) const {
			return mData + value.ToRaw();
		}
		template<typename _Type = Type>
			requires(HasIdentScale())
		friend constexpr Self<decltype(_Type() + Type())> operator+(const StandardUnit<_Type> &value, const Self<> &self) {
			return value.ToRaw() + self.ToRaw();
		}
		template<typename SiblingUnit = Self<>>
			requires(std::is_base_of_v<_::Base, SiblingUnit> && SiblingUnit::template hasSameScale<Scale>() && HasSameStdUnitBase<SiblingUnit>())
		constexpr auto operator-(const SiblingUnit value) const {
			using RawType = decltype(Type() - SiblingUnit().ToRaw());
			using ResultOffset = Offset::template Diff<_::Fraction<SiblingUnit::GetOffsetNumerator(), SiblingUnit::GetOffsetDenominator()>>;
			using AccType = decltype(AccuracyType() * SiblingUnit::SCALE);
			using NewUnit = AdaptiveSelf<RawType, StandardUnit, Scale, ResultOffset, AccType>;
			return NewUnit{mData - value.ToRaw()};
		}
		template<typename _Type = Type>
			requires(HasIdentScale())
		constexpr Self<decltype(Type() - _Type())> operator-(const StandardUnit<_Type> value) const {
			return mData - value.ToRaw();
		}
		template<typename _Type = Type>
			requires(HasIdentScale())
		friend constexpr auto operator-(const StandardUnit<_Type> &value, const Self<> &self) {
			using NewUnit = decltype(-Self<decltype(_Type() - Type())>());
			return NewUnit{value.ToRaw() - self.ToRaw()};
		}
		template<typename _Type, template<typename> class _StandardUnit, class _Scale, class _Offset, typename _AccuracyType = f64>
			requires(HasNoOffset() && _Offset::IsZero())
		constexpr auto operator*(const GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType> value) const {
			using RawType = decltype(Type() * _Type());
			using ComposeType = decltype(StandardUnit<Type>() * _StandardUnit<_Type>());
			using AccType = decltype(AccuracyType() * _AccuracyType());
			using ResultScale = Scale::template Product<_Scale>;
			using NewUnit = AdaptiveSelf<RawType, Decomposer<ComposeType>::template OuterType, ResultScale, _::Fraction<0, 1>, AccType>;
			return NewUnit{mData * value.ToRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		constexpr typename StdToSelf<decltype(StandardUnit<Type>() * _StdUnitT())>::NewUnit
			operator*(const _StdUnitT value) const {
			return mData * value.ToRaw();
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		friend constexpr typename StdToSelf<decltype(_StdUnitT() * StandardUnit<Type>())>::NewUnit
			operator*(const _StdUnitT &value, const Self<> &self) {
			return value.ToRaw() * self.ToRaw();
		}
		template<typename _Type, template<typename> class _StandardUnit, class _Scale, class _Offset, typename _AccuracyType = f64>
			requires(HasNoOffset() && _Offset::IsZero())
		constexpr auto operator/(const GenerativeUnit<_Type, _StandardUnit, _Scale, _Offset, _AccuracyType> value) const {
			using RawType = decltype(Type() / _Type());
			using ComposeType = decltype(StandardUnit<Type>() / _StandardUnit<_Type>());
			using AccType = decltype(AccuracyType() / _AccuracyType());
			using ResultScale = Scale::template Quotient<_Scale>;
			using NewUnit = AdaptiveSelf<RawType, Decomposer<ComposeType>::template OuterType, ResultScale, _::Fraction<0, 1>, AccType>;
			return NewUnit{mData / value.ToRaw()};
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		constexpr typename StdToSelf<decltype(StandardUnit<Type>() / _StdUnitT())>::NewUnit
			operator/(const _StdUnitT value) const {
			return mData / value.ToRaw();
		}
		template<typename _StdUnitT = StandardUnit<Type>>
			requires(HasNoOffset())
		friend constexpr typename StdToSelf<decltype(_StdUnitT() / StandardUnit<Type>())>::NewUnit
			operator/(const _StdUnitT &value, const Self<> &self) {
			return value.ToRaw() / self.ToRaw();
		}

		friend std::ostream &operator<<(std::ostream &os, const Self<> &obj) {
			os << obj.mData;
			return os;
		}

		/* Other methods */;

		template<typename _Type = Type>
		constexpr Self<_Type> Cast() const {
			return mData;
		}

		constexpr Type ToRaw() const {
			return mData;
		}

		template<typename _Type = Type>
		constexpr StandardUnit<_Type> ToStandardUnit() const {
			if constexpr(SCALE == 1)
				return mData + OFFSET;
			else if constexpr(OFFSET == 0)
				return mData * SCALE;
			else
				return mData * SCALE + OFFSET;
		}

		static consteval auto GetScaleNumerator() {
			return Scale::Num;
		}
		static consteval auto GetScaleDenominator() {
			return Scale::Denom;
		}
		static consteval auto GetOffsetNumerator() {
			return Offset::Num;
		}
		static consteval auto GetOffsetDenominator() {
			return Offset::Denom;
		}

	private:
		template<typename _Type = Type>
		static constexpr Type FromStandardUnit(const StandardUnit<_Type> value) {
			if constexpr(SCALE == 1)
				return value.ToRaw() - OFFSET;
			else if constexpr(OFFSET == 0)
				return value.ToRaw() / SCALE;
			else
				return (value.ToRaw() - OFFSET) / SCALE;
		}
		template<typename SiblingUnit = Self<>>
			requires(HasSameStdUnitBase<SiblingUnit>())
		static constexpr Type FromOtherNStd(const SiblingUnit value) {
			if constexpr(std::is_same_v<Self<decltype(SiblingUnit().ToRaw())>, SiblingUnit>) {
				return value;
			} else {
				constexpr auto compScale = SiblingUnit::SCALE * SCALE;
				constexpr auto compOffset = (SiblingUnit::OFFSET - OFFSET) * SCALE;

				return value.ToRaw() * compScale + compOffset;
			}
		}

	private:
		/* TODO: Consider how to handle `AccuracyType` */;
		static constexpr AccuracyType SCALE = Scale::template ToDecimal<AccuracyType>();
		static constexpr AccuracyType OFFSET = Offset::template ToDecimal<AccuracyType>();
		Type mData;
	};

	namespace _
	{
		class FractionParser
		{
		private:
			static consteval i64 Exp2(i64 value) noexcept {
				if(value < 0) return 0;
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
			static consteval i64 CalcExponent(f64 value) noexcept {
				return FlooredLog2(std::abs(value));
			}

		public:
			static consteval i64 Denominator( f64 value) noexcept {
				const i64 exp = CalcExponent(value);
				return Exp2(exp < 0L ? 62L : 62L - exp);
			}
			static consteval i64 Numerator( f64 value) noexcept {
				return value * Denominator(value);
			}
			static consteval bool IsConvertible( f64 value) noexcept {
				const i64 exp = CalcExponent(value);
				return exp > -11 && exp < 63;
			}
		};

		template<typename Type, template<typename> class StdU, i64 ScNum, i64 ScDenom, i64 OffNum = 0, i64 OffDenom = 1>
		using Simplifier = GenerativeUnit<Type, StdU, typename Fraction<ScNum, ScDenom>::Norm, typename Fraction<OffNum, OffDenom>::Norm>;

		static constexpr f64 DEG2RAD = M_PI / 180;
	}

#define GENERATE_NSTD_FROM_DOUBLE(Name, StdType, Scale)                                                   \
	static_assert(_::FractionParser::IsConvertible(Scale), "Value '" #Scale "' out of supported range!"); \
	template<typename T = f64>                                                                            \
	using Name = _::Simplifier<T, StdType, _::FractionParser::Numerator(Scale), _::FractionParser::Denominator(Scale)>;

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

#undef GENERATE_NSTD_FROM_DOUBLE
#undef GENERATE_NSTD_FROM_FRACTION
#undef GENERATE_OFFSETABLE_NSTD_FROM_FRACTION
}