#pragma once
//
#include <cmath>
#include <iostream>
#include <string_view>
//
#include "SafePhysics.hpp"
#include "Utils/CompileTime.hpp"

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
		constexpr GenerativeUnit(Type data) noexcept : mData(data) {}
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		constexpr GenerativeUnit(const Sibling<OtherType, OTHER_PREFIX>& value) noexcept : mData(ScaleUnit(value)){};

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
			using NewType = decltype(-Type{});
			return Sibling<NewType>{-mData};
		}
		template<Arithmetic OtherType = Type>
		constexpr auto operator+(const Sibling<OtherType>& value) const noexcept {
			using NewType = decltype(Type{} + OtherType{});
			return Sibling<NewType>{mData + value.ToRaw()};
		}
		template<Arithmetic OtherType = Type>
		constexpr auto operator-(const Sibling<OtherType>& value) const noexcept {
			using NewType = decltype(Type{} - OtherType{});
			return Sibling<NewType>{mData - value.ToRaw()};
		}
		template<Arithmetic OtherType, i8 N_M, i8 N_S, i8 N_KG, i8 N_A, i8 N_K, i8 N_MOL, i8 N_CD, i8 N_RAD, i8 N_SR, i8 OTHER_PREFIX>
		constexpr auto operator*(const GenerativeUnit<OtherType, N_M, N_S, N_KG, N_A, N_K, N_MOL, N_CD, N_RAD, N_SR, OTHER_PREFIX> value) const noexcept {
			using NewType = GenerativeUnit<decltype(Type{} * OtherType{}), M + N_M, S + N_S, KG + N_KG, A + N_A, K + N_K,
										   MOL + N_MOL, CD + N_CD, RAD + N_RAD, SR + N_SR, PREFIX + OTHER_PREFIX>;
			return NewType{mData * value.ToRaw()};
		}
		template<Arithmetic OtherType, i8 N_M, i8 N_S, i8 N_KG, i8 N_A, i8 N_K, i8 N_MOL, i8 N_CD, i8 N_RAD, i8 N_SR, i8 OTHER_PREFIX>
		constexpr auto operator/(const GenerativeUnit<OtherType, N_M, N_S, N_KG, N_A, N_K, N_MOL, N_CD, N_RAD, N_SR, OTHER_PREFIX> value) const noexcept {
			using NewType = GenerativeUnit<decltype(Type{} / OtherType{}), M - N_M, S - N_S, KG - N_KG, A - N_A, K - N_K,
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

		template<u8 EXPONENT>
		constexpr auto Power() const noexcept {
			if constexpr(EXPONENT == 0) {
				return Scale<Type>{1};
			} else {
				return *this * Power<EXPONENT - 1>();
			}
		}

		constexpr auto Sqrt() const noexcept {
			static_assert(M % 2 == 0);
			static_assert(S % 2 == 0);
			static_assert(KG % 2 == 0);
			static_assert(A % 2 == 0);
			static_assert(K % 2 == 0);
			static_assert(MOL % 2 == 0);
			static_assert(CD % 2 == 0);
			static_assert(RAD % 2 == 0);
			static_assert(SR % 2 == 0);
			static_assert(PREFIX % 2 == 0);
			using RawType = std::conditional_t<std::is_same_v<Type, f128>, f128, f64>;
			using NewType = GenerativeUnit<RawType, M / 2, S / 2, KG / 2, A / 2, K / 2, MOL / 2, CD / 2, RAD / 2, SR / 2, PREFIX / 2>;
			return NewType{mData >= 0 && mData < std::numeric_limits<double>::infinity()
							   ? SqrtNewtonRaphson(mData, mData, 0)
							   : std::numeric_limits<f64>::quiet_NaN()};
		}

		static consteval bool HasNoPrefix() noexcept {
			return PREFIX == 0;
		}

		static consteval std::string_view GetTypeView() noexcept {
			return sTypeText.View();
		}

		static consteval std::string_view GetTypeCString() noexcept {
			return sTypeText.CString();
		}

	private:
		template<Arithmetic OtherType = Type, i8 OTHER_PREFIX = PREFIX>
		static constexpr auto ScaleUnit(const Sibling<OtherType, OTHER_PREFIX> value) noexcept {
			using NewType = decltype(Type{} * OtherType{});

			if constexpr(OTHER_PREFIX - PREFIX > 0) {
				constexpr NewType scaleFactor = Pow10(OTHER_PREFIX - PREFIX);
				return value.ToRaw() * scaleFactor;
			} else if constexpr(OTHER_PREFIX - PREFIX < 0) {
				constexpr NewType scaleFactor = Pow10(PREFIX - OTHER_PREFIX);
				return value.ToRaw() / scaleFactor;
			} else {
				return value.ToRaw();
			}
		}

		static consteval i64 Pow10(const u8 pow) noexcept {
			if(pow > 0) {
				return 10L * Pow10(pow - 1);
			}
			return 1;
		}

		static constexpr f64 SqrtNewtonRaphson(f64 x, f64 curr, f64 prev) noexcept {
			return curr == prev ? curr : SqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
		}

		static consteval auto GetTypeText() noexcept {
			constexpr auto text = Utils::CompileTime::String("").Join(
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

			if constexpr(text.Size() != 0) {
				constexpr auto prefixSize = Utils::CompileTime::String(" * ").Size();
				return text.template SubString<prefixSize, text.Size() - prefixSize>();
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
		static constexpr auto sTypeText = GetTypeText();
		Type mData;
	};

#define GENERATE_SI_UNIT(Name, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift)                 \
	template<Arithmetic T = f64>                                                            \
	using Nano##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 9>;  \
	template<Arithmetic T = f64>                                                            \
	using Micro##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 6>; \
	template<Arithmetic T = f64>                                                            \
	using Milli##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 3>; \
	template<Arithmetic T = f64>                                                            \
	using Centi##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 2>; \
	template<Arithmetic T = f64>                                                            \
	using Deci##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 1>;  \
	template<Arithmetic T = f64>                                                            \
	using Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift>;            \
	template<Arithmetic T = f64>                                                            \
	using Deca##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 1>;  \
	template<Arithmetic T = f64>                                                            \
	using Hecto##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 2>; \
	template<Arithmetic T = f64>                                                            \
	using Kilo##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 3>;  \
	template<Arithmetic T = f64>                                                            \
	using Mega##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 6>;  \
	template<Arithmetic T = f64>                                                            \
	using Giga##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 9>;  \
	template<Arithmetic T = f64, i8 Power = 0>                                              \
	using Any##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + Power>;

	/* --- GENERATE NEEDED UNITS HERE --- */

	/* Basic units */;

	GENERATE_SI_UNIT(Meters, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Seconds, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Grams, 0, 0, 1, 0, 0, 0, 0, 0, 0, -3);
	GENERATE_SI_UNIT(Amperes, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Kelvins, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Moles, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Candelas, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0);
	GENERATE_SI_UNIT(Radians, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0);
	GENERATE_SI_UNIT(Steradians, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0);

	/* Kinematics */;

	GENERATE_SI_UNIT(MetersPerSecond, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(RadiansPerSecond, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0);
	GENERATE_SI_UNIT(MetersPerSecondSquare, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(RadiansPerSecondSquare, 0, -2, 0, 0, 0, 0, 0, 1, 0, 0);
	GENERATE_SI_UNIT(MetersPerRadian, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0);
	GENERATE_SI_UNIT(RadialMeters, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0);

	/* Dynamics */;

	GENERATE_SI_UNIT(GramMeters, 1, 0, 1, 0, 0, 0, 0, -2, 0, -3);
	GENERATE_SI_UNIT(GramMetersSquare, 2, 0, 1, 0, 0, 0, 0, -2, 0, -3);
	GENERATE_SI_UNIT(NewtonSeconds, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(NewtonMeterSeconds, 2, -1, 1, 0, 0, 0, 0, -1, 0, 0);
	GENERATE_SI_UNIT(Newtons, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(NewtonMeters, 2, -2, 1, 0, 0, 0, 0, -1, 0, 0);
	GENERATE_SI_UNIT(Joules, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Watts, 2, -3, 1, 0, 0, 0, 0, 0, 0, 0);

	/* Electrics */;

	GENERATE_SI_UNIT(Coulombs, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(NewtonsPerCoulomb, 1, -3, 1, -1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(VoltsPerMeter, 1, -3, 1, -1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Volts, 2, -3, 1, -1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Ohms, 2, -3, 1, -2, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Siemens, -2, 3, -1, 2, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Farads, -2, 4, -1, 2, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(FaradsPerMeter, -3, 4, -1, 2, 0, 0, 0, 0, 0, 0);

	/* Magnetism */;

	GENERATE_SI_UNIT(AmperesPerMeter, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Teslas, 0, -2, 1, -1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Webers, 2, -2, 1, -1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(TeslaMeters, 1, -2, 1, -1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(Henries, 2, -2, 1, -2, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(HenriesPerMeter, 1, -2, 1, -2, 0, 0, 0, 0, 0, 0);

	/* Materials */;

	GENERATE_SI_UNIT(MetersSquare, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(MetersCube, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(GramsPerMeterSquare, -2, 0, 1, 0, 0, 0, 0, 0, 0, -3);
	GENERATE_SI_UNIT(GramsPerMeterCube, -3, 0, 1, 0, 0, 0, 0, 0, 0, -3);
	GENERATE_SI_UNIT(Pascals, -1, -2, 1, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(CoulombsPerMeterCube, -3, 1, 0, 1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(AmperesPerMeterSquare, -2, 0, 0, 1, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(GramsPerMole, 0, 0, 1, 0, 0, -1, 0, 0, 0, -3);
	GENERATE_SI_UNIT(PartsPerMole, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0);

	/* Vibrations and waves */;

	GENERATE_SI_UNIT(Hertzes, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(RadiansPerMeter, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0);

	/* Thermodynamics */;

	GENERATE_SI_UNIT(WattsPerMeterSquare, 0, -3, 1, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(JoulesPerKelvin, 2, -2, 1, 0, -1, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(JoulesPerKilogramKelvin, 2, -2, 0, 0, -1, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(WattsPerMeterKelvin, 1, -3, 1, 0, -1, 0, 0, 0, 0, 0);

	/* Optics */;

	GENERATE_SI_UNIT(Lumens, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0);
	GENERATE_SI_UNIT(Luxes, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0);

	/* Quantum mechanics */;

	GENERATE_SI_UNIT(JouleSeconds, 2, -1, 1, 0, 0, 0, 0, 0, 0, 0);
	GENERATE_SI_UNIT(JouleSecondsPerRadian, 2, -1, 1, 0, 0, 0, 0, -1, 0, 0);

#undef GENERATE_SI_UNIT
}