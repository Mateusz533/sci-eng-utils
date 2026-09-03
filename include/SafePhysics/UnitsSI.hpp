#pragma once
//
#include <cmath>
#include <limits>
#include <ostream>
#include <string_view>
#include <type_traits>
#include <utility>
//
#include "Utils/CompileTime.hpp"
#include "Utils/Types.hpp"

namespace Physics
{
	using namespace Utils::Types;
}

namespace Physics::Units::SI
{
	template<Arithmetic Type, i8 M, i8 S, i8 KG, i8 A, i8 K, i8 MOL, i8 CD, i8 RAD, i8 SR, i8 PREFIX>
		requires(PREFIX == 0 || (M != 0 || S != 0 || KG != 0 || A != 0 || K != 0 || MOL != 0 || CD != 0 || RAD != 0 || SR != 0))
	class GenerativeUnit;

	template<Arithmetic T = f64>
	using Scale = GenerativeUnit<T, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>;

	template<Arithmetic Type, i8 M, i8 S, i8 KG, i8 A, i8 K, i8 MOL, i8 CD, i8 RAD, i8 SR, i8 PREFIX>
		requires(PREFIX == 0 || (M != 0 || S != 0 || KG != 0 || A != 0 || K != 0 || MOL != 0 || CD != 0 || RAD != 0 || SR != 0))
	class GenerativeUnit
	{
	private:
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		using Sibling = GenerativeUnit<OtherType, M, S, KG, A, K, MOL, CD, RAD, SR, OTHER_PREFIX>;

		using Self = GenerativeUnit;

	public:
		/* Constructors */;

		constexpr GenerativeUnit() = default;
		constexpr GenerativeUnit(const Self& value) = default;
		constexpr GenerativeUnit(Self&& value) = default;
		constexpr GenerativeUnit(Type data) noexcept : mData{data} {}
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		constexpr GenerativeUnit(const Sibling<OtherType, OTHER_PREFIX>& value) noexcept : mData{ScaleUnit(value)} {};

		/* Assignment operators */;

		constexpr Self& operator=(const Self& value) = default;
		constexpr Self& operator=(Self&& value) = default;
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		constexpr Self& operator=(const Sibling<OtherType, OTHER_PREFIX>& value) noexcept {
			mData = ScaleUnit(value);
			return *this;
		}
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		constexpr Self& operator+=(const Sibling<OtherType, OTHER_PREFIX>& value) noexcept {
			mData += ScaleUnit(value);
			return *this;
		}
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		constexpr Self& operator-=(const Sibling<OtherType, OTHER_PREFIX>& value) noexcept {
			mData -= ScaleUnit(value);
			return *this;
		}
		template<Arithmetic OtherType = Type>
		constexpr Self& operator*=(const Scale<OtherType>& value) noexcept {
			mData *= value.ToRaw();
			return *this;
		}
		template<Arithmetic OtherType = Type>
		constexpr Self& operator/=(const Scale<OtherType>& value) noexcept {
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
			using NewType = decltype(-std::declval<Type>());
			return Sibling<NewType>{-mData};
		}
		template<Arithmetic OtherType = Type>
		constexpr auto operator+(const Sibling<OtherType>& value) const noexcept {
			using NewType = decltype(std::declval<Type>() + std::declval<OtherType>());
			return Sibling<NewType>{mData + value.ToRaw()};
		}
		template<Arithmetic OtherType = Type>
		constexpr auto operator-(const Sibling<OtherType>& value) const noexcept {
			using NewType = decltype(std::declval<Type>() - std::declval<OtherType>());
			return Sibling<NewType>{mData - value.ToRaw()};
		}
		template<Arithmetic OtherType, i8 N_M, i8 N_S, i8 N_KG, i8 N_A, i8 N_K, i8 N_MOL, i8 N_CD, i8 N_RAD, i8 N_SR, i8 OTHER_PREFIX>
		constexpr auto operator*(const GenerativeUnit<OtherType, N_M, N_S, N_KG, N_A, N_K, N_MOL, N_CD, N_RAD, N_SR, OTHER_PREFIX>& value) const noexcept {
			using RawType = decltype(std::declval<Type>() * std::declval<OtherType>());
			using NewType = GenerativeUnit<RawType, M + N_M, S + N_S, KG + N_KG, A + N_A, K + N_K,
										   MOL + N_MOL, CD + N_CD, RAD + N_RAD, SR + N_SR, PREFIX + OTHER_PREFIX>;
			return NewType{mData * value.ToRaw()};
		}
		template<Arithmetic OtherType, i8 N_M, i8 N_S, i8 N_KG, i8 N_A, i8 N_K, i8 N_MOL, i8 N_CD, i8 N_RAD, i8 N_SR, i8 OTHER_PREFIX>
		constexpr auto operator/(const GenerativeUnit<OtherType, N_M, N_S, N_KG, N_A, N_K, N_MOL, N_CD, N_RAD, N_SR, OTHER_PREFIX>& value) const noexcept {
			using RawType = decltype(std::declval<Type>() / std::declval<OtherType>());
			using NewType = GenerativeUnit<RawType, M - N_M, S - N_S, KG - N_KG, A - N_A, K - N_K,
										   MOL - N_MOL, CD - N_CD, RAD - N_RAD, SR - N_SR, PREFIX - OTHER_PREFIX>;
			return NewType{mData / value.ToRaw()};
		}

		friend std::ostream& operator<<(std::ostream& os, const Self& obj) {
			os << obj.mData << ' ';
			if constexpr(PREFIX != 0) {
				os << "* ";
			}
			os << Self::GetTypeView();
			return os;
		}

		/* Other methods */;

		template<Arithmetic OtherType = Type>
		constexpr Sibling<OtherType> Cast() const noexcept {
			return static_cast<OtherType>(mData);
		}

		constexpr Type ToRaw() const noexcept {
			return mData;
		}

		constexpr GenerativeUnit<Type, M, S, KG, A, K, MOL, CD, RAD, SR, 0> ToCoherentUnit() const noexcept {
			if constexpr(PREFIX > 0) {
				constexpr Type SCALE_FACTOR = Utils::CompileTime::Exp10(+PREFIX);
				return mData.ToRaw() * SCALE_FACTOR;
			} else if constexpr(PREFIX < 0) {
				constexpr Type SCALE_FACTOR = Utils::CompileTime::Exp10(-PREFIX);
				return mData.ToRaw() / SCALE_FACTOR;
			} else {
				return mData.ToRaw();
			}
		}

		template<i8 EXPONENT>
		constexpr auto Power() const noexcept {
			if constexpr(EXPONENT == 0) {
				return Scale<Type>{1};
			} else if constexpr(EXPONENT < 0) {
				return Scale<Type>{1} / Power<-EXPONENT>();
			} else {
				return *this * Power<EXPONENT - 1>();
			}
		}

		template<u8 N>
			requires((N > 0) && (M % N == 0) && (S % N == 0) && (KG % N == 0) && (A % N == 0) && (K % N == 0) &&
					 (MOL % N == 0) && (CD % N == 0) && (RAD % N == 0) && (SR % N == 0) && (PREFIX % N == 0))
		constexpr auto Root() const noexcept {
			using RawType = std::conditional_t<std::is_same_v<Type, f128>, f128, f64>;
			using NewType = GenerativeUnit<RawType, M / N, S / N, KG / N, A / N, K / N, MOL / N, CD / N, RAD / N, SR / N, PREFIX / N>;
			return NewType{RawType{0} <= mData && mData < std::numeric_limits<RawType>::infinity()
							   ? Utils::CompileTime::RootNewtonRaphson<RawType, N>(mData, mData, 0)
							   : std::numeric_limits<f64>::quiet_NaN()};
		}
		constexpr auto Sqrt() const noexcept {
			return Root<2>();
		}
		constexpr auto Cbrt() const noexcept {
			return Root<3>();
		}

		static consteval bool HasNoPrefix() noexcept {
			return PREFIX == 0;
		}

		static consteval i8 GetPrefix() noexcept {
			return PREFIX;
		}

		static consteval std::string_view GetTypeView() noexcept {
			return TYPE_TEXT.View();
		}

		static consteval std::string_view GetTypeCString() noexcept {
			return TYPE_TEXT.CString();
		}

	private:
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		static constexpr auto ScaleUnit(const Sibling<OtherType, OTHER_PREFIX>& value) noexcept {
			using NewType = decltype(std::declval<Type>() * std::declval<OtherType>());

			if constexpr(OTHER_PREFIX - PREFIX > 0) {
				constexpr NewType SCALE_FACTOR = Utils::CompileTime::Exp10(OTHER_PREFIX - PREFIX);
				return value.ToRaw() * SCALE_FACTOR;
			} else if constexpr(OTHER_PREFIX - PREFIX < 0) {
				constexpr NewType SCALE_FACTOR = Utils::CompileTime::Exp10(PREFIX - OTHER_PREFIX);
				return value.ToRaw() / SCALE_FACTOR;
			} else {
				return value.ToRaw();
			}
		}

		static consteval auto GetTypeText() noexcept {
			constexpr auto TEXT = Utils::CompileTime::String("").Join(
				DimensionToString<PREFIX>("10"),
				DimensionToString<M>("m"),
				DimensionToString<S>("s"),
				DimensionToString<KG>("kg"),
				DimensionToString<A>("A"),
				DimensionToString<K>("K"),
				DimensionToString<MOL>("mol"),
				DimensionToString<CD>("cd"),
				DimensionToString<RAD>("rad"),
				DimensionToString<SR>("sr"));

			if constexpr(TEXT.Size() != 0) {
				constexpr auto PREFIX_SIZE = Utils::CompileTime::String(" * ").Size();
				return TEXT.template SubString<PREFIX_SIZE, TEXT.Size() - PREFIX_SIZE>();
			} else {
				return Utils::CompileTime::String("(dimensionless)");
			}
		}

		template<i8 Dimension, usize N>
		static consteval auto DimensionToString(const char (&baseText)[N]) noexcept {
			if constexpr(Dimension == 0) {
				return Utils::CompileTime::String("");
			} else {
				return Utils::CompileTime::String(" * ").Join(
					Utils::CompileTime::String(baseText),
					Utils::CompileTime::String("^"),
					Utils::CompileTime::IntegerToString<Dimension>());
			}
		}

	private:
		static constexpr auto TYPE_TEXT = GetTypeText();
		Type mData;
	};

#define GENERATE_SI_UNIT(Name, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift, ScaleFactor, Prefix)                      \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Nano##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) - 9 * (ScaleFactor)>;  \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Micro##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) - 6 * (ScaleFactor)>; \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Milli##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) - 3 * (ScaleFactor)>; \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Centi##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) - 2 * (ScaleFactor)>; \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Deci##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) - 1 * (ScaleFactor)>;  \
	template<Arithmetic T = f64> /* NOLINTNEXTLINE(bugprone-macro-parentheses) */                                     \
	using Prefix##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift>;                              \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Deca##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) + 1 * (ScaleFactor)>;  \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Hecto##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) + 2 * (ScaleFactor)>; \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Kilo##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) + 3 * (ScaleFactor)>;  \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Mega##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) + 6 * (ScaleFactor)>;  \
	template<Arithmetic T = f64>                                                                                      \
	using Prefix##Giga##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, (BaseShift) + 9 * (ScaleFactor)>;

	/* --- GENERATE NEEDED UNITS HERE --- */

	/* Basic units */;

	GENERATE_SI_UNIT(Meters, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Seconds, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Grams, 0, 0, 1, 0, 0, 0, 0, 0, 0, -3, 1, );
	GENERATE_SI_UNIT(Amperes, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Kelvins, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Moles, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Candelas, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Radians, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, );
	GENERATE_SI_UNIT(Steradians, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, );

	/* Kinematics */;

	GENERATE_SI_UNIT(MetersPerSecond, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(RadiansPerSecond, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 1, );
	GENERATE_SI_UNIT(MetersPerSecondSquared, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(RadiansPerSecondSquared, 0, -2, 0, 0, 0, 0, 0, 1, 0, 0, 1, );
	GENERATE_SI_UNIT(MetersPerRadian, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, );
	GENERATE_SI_UNIT(RadialMeters, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, );

	/* Dynamics */;

	GENERATE_SI_UNIT(GramMeters, 1, 0, 1, 0, 0, 0, 0, -2, 0, -3, 1, );
	GENERATE_SI_UNIT(GramMetersSquared, 2, 0, 1, 0, 0, 0, 0, -2, 0, -3, 1, );
	GENERATE_SI_UNIT(NewtonSeconds, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(NewtonMeterSeconds, 2, -1, 1, 0, 0, 0, 0, -1, 0, 0, 1, );
	GENERATE_SI_UNIT(Newtons, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(NewtonMeters, 2, -2, 1, 0, 0, 0, 0, -1, 0, 0, 1, );
	GENERATE_SI_UNIT(Joules, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Watts, 2, -3, 1, 0, 0, 0, 0, 0, 0, 0, 1, );

	/* Electrics */;

	GENERATE_SI_UNIT(Coulombs, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(NewtonsPerCoulomb, 1, -3, 1, -1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(VoltsPerMeter, 1, -3, 1, -1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Volts, 2, -3, 1, -1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Ohms, 2, -3, 1, -2, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Siemens, -2, 3, -1, 2, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Farads, -2, 4, -1, 2, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(FaradsPerMeter, -3, 4, -1, 2, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(AmperesPerSecond, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 1, );

	/* Magnetism */;

	GENERATE_SI_UNIT(AmperesPerMeter, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Teslas, 0, -2, 1, -1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Webers, 2, -2, 1, -1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(TeslaMeters, 1, -2, 1, -1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Henries, 2, -2, 1, -2, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(HenriesPerMeter, 1, -2, 1, -2, 0, 0, 0, 0, 0, 0, 1, );

	/* Materials */;

	GENERATE_SI_UNIT(Meters, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, Square);
	GENERATE_SI_UNIT(Meters, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, Cubic);
	GENERATE_SI_UNIT(GramsPerSquareMeter, -2, 0, 1, 0, 0, 0, 0, 0, 0, -3, 1, );
	GENERATE_SI_UNIT(GramsPerCubicMeter, -3, 0, 1, 0, 0, 0, 0, 0, 0, -3, 1, );
	GENERATE_SI_UNIT(MetersPerKilogram, 3, 0, -1, 0, 0, 0, 0, 0, 0, 0, 3, Cubic);
	GENERATE_SI_UNIT(Pascals, -1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(PascalSeconds, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(NewtonsPerMeter, 0, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(JoulesPerCubicMeter, -1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(CoulombsPerSquareMeter, -2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(CoulombsPerCubicMeter, -3, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(AmperesPerSquareMeter, -2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(GramsPerMole, 0, 0, 1, 0, 0, -1, 0, 0, 0, -3, 1, );
	GENERATE_SI_UNIT(PartsPerMole, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(MolesPerCubicMeter, -3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Katals, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(KatalsPerCubicMeter, -3, -1, 0, 0, 0, 1, 0, 0, 0, 0, 1, );

	/* Vibrations and waves */;

	GENERATE_SI_UNIT(Hertzes, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(RadiansPerMeter, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, );

	/* Thermodynamics */;

	GENERATE_SI_UNIT(WattsPerSquareMeter, 0, -3, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(JoulesPerKelvin, 2, -2, 1, 0, -1, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(JoulesPerKilogramKelvin, 2, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(JoulesPerKilogram, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(WattsPerMeterKelvin, 1, -3, 1, 0, -1, 0, 0, 0, 0, 0, 1, );

	/* Optics */;

	GENERATE_SI_UNIT(Lumens, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, );
	GENERATE_SI_UNIT(Luxes, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, );
	GENERATE_SI_UNIT(CandelasPerSquareMeter, -2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(WattsPerSteradian, 2, -3, 1, 0, 0, 0, 0, 0, -1, 0, 1, );
	GENERATE_SI_UNIT(WattsPerSquareMeterSteradian, 0, -3, 1, 0, 0, 0, 0, 0, -1, 0, 1, );

	/* Quantum mechanics */;

	GENERATE_SI_UNIT(JouleSeconds, 2, -1, 1, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(JouleSecondsPerRadian, 2, -1, 1, 0, 0, 0, 0, -1, 0, 0, 1, );

	/* Nuclear physics */;

	GENERATE_SI_UNIT(Becquerels, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Grays, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(Sieverts, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(GraysPerSecond, 2, -3, 0, 0, 0, 0, 0, 0, 0, 0, 1, );
	GENERATE_SI_UNIT(CoulombsPerKilogram, 0, 1, -1, 1, 0, 0, 0, 0, 0, 0, 1, );

#undef GENERATE_SI_UNIT
}