#pragma once
//
#include <array>
#include <iostream>
#include <string_view>
//
#include "SafePhysics.h"

namespace Physics::Units::SI
{
	template<typename Type, i8 M, i8 S, i8 Kg, i8 A, i8 K, i8 Mol, i8 Cd, i8 Rad, i8 Sr, i8 Prefix>
		requires(std::is_arithmetic_v<Type> && !!(Prefix == 0 || (M || S || Kg || A || K || Mol || Cd || Rad || Sr)))
	class GenerativeUnit
	{
	private:
		template<typename _Type = Type, i8 _Prefix = Prefix>
		using Self = GenerativeUnit<_Type, M, S, Kg, A, K, Mol, Cd, Rad, Sr, _Prefix>;

	public:
		/* Constructors */;

		constexpr GenerativeUnit(const Self<> &value) : mData(value.mData) {}
		constexpr GenerativeUnit(const Self<> &&value) : mData(value.mData) {}
		constexpr GenerativeUnit() : mData(0) {}
		constexpr GenerativeUnit(const Type data) : mData(data) {}
		template<typename _Type = Type, i8 _Prefix = Prefix>
		constexpr GenerativeUnit(const Self<_Type, _Prefix> value) : mData(scale(value)){};

		/* Assignment operators */;

		constexpr Self<> &operator=(const Self<> &value) {
			mData = value.mData;
			return *this;
		}
		constexpr Self<> &operator=(const Self<> &&value) {
			mData = value.mData;
			return *this;
		}
		template<typename _Type = Type, i8 _Prefix = Prefix>
		constexpr Self<> &operator=(const Self<_Type, _Prefix> value) {
			mData = scale(value);
			return *this;
		}
		template<typename _Type = Type, i8 _Prefix = Prefix>
		constexpr Self<> &operator+=(const Self<_Type, _Prefix> value) {
			mData += scale(value);
			return *this;
		}
		template<typename _Type = Type, i8 _Prefix = Prefix>
		constexpr Self<> &operator-=(const Self<_Type, _Prefix> value) {
			mData -= scale(value);
			return *this;
		}
		template<typename _Type = Type>
		constexpr Self<> &operator*=(const GenerativeUnit<_Type, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> value) {
			mData *= value.toRaw();
			return *this;
		}
		template<typename _Type = Type>
		constexpr Self<> &operator/=(const GenerativeUnit<_Type, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0> value) {
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

		constexpr auto operator-() const {
			return Self<decltype(-Type())>{-mData};
		}
		template<typename _Type = Type>
		constexpr auto operator+(const Self<_Type> value) const {
			using NewType = decltype(Type() + _Type());
			return Self<NewType>{mData + value.toRaw()};
		}
		template<typename _Type = Type>
		constexpr auto operator-(const Self<_Type> value) const {
			using NewType = decltype(Type() - _Type());
			return Self<NewType>{mData - value.toRaw()};
		}
		template<typename _Type, i8 _M, i8 _S, i8 _Kg, i8 _A, i8 _K, i8 _Mol, i8 _Cd, i8 _Rad, i8 _Sr, i8 _Prefix>
		constexpr auto operator*(const GenerativeUnit<_Type, _M, _S, _Kg, _A, _K, _Mol, _Cd, _Rad, _Sr, _Prefix> value) const {
			return GenerativeUnit<decltype(Type() * _Type()), M + _M, S + _S, Kg + _Kg, A + _A, K + _K, Mol + _Mol, Cd + _Cd,
								  Rad + _Rad, Sr + _Sr, Prefix + _Prefix>{mData * value.toRaw()};
		}
		template<typename _Type, i8 _M, i8 _S, i8 _Kg, i8 _A, i8 _K, i8 _Mol, i8 _Cd, i8 _Rad, i8 _Sr, i8 _Prefix>
		constexpr auto operator/(const GenerativeUnit<_Type, _M, _S, _Kg, _A, _K, _Mol, _Cd, _Rad, _Sr, _Prefix> value) const {
			return GenerativeUnit<decltype(Type() / _Type()), M - _M, S - _S, Kg - _Kg, A - _A, K - _K, Mol - _Mol, Cd - _Cd,
								  Rad - _Rad, Sr - _Sr, Prefix - _Prefix>{mData / value.toRaw()};
		}

		friend std::ostream &operator<<(std::ostream &os, const Self<> &obj) {
			os << obj.mData << ' ';
			if constexpr(Prefix != 0)
				os << "* ";
			os << obj.getType();
			return os;
		}

		/* Other methods */;

		template<typename _Type = Type>
		constexpr Self<_Type> cast() const {
			return static_cast<_Type>(mData);
		}

		constexpr Type toRaw() const {
			return mData;
		}

		template<u8 EXPONENT>
		constexpr auto power() const {
			if constexpr(EXPONENT == 0)
				return GenerativeUnit<Type, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>{1};
			else
				return *this * power<EXPONENT - 1>();
		}

		constexpr auto sqrt() const {
			static_assert(M % 2 == 0);
			static_assert(S % 2 == 0);
			static_assert(Kg % 2 == 0);
			static_assert(A % 2 == 0);
			static_assert(K % 2 == 0);
			static_assert(Mol % 2 == 0);
			static_assert(Cd % 2 == 0);
			static_assert(Rad % 2 == 0);
			static_assert(Sr % 2 == 0);
			static_assert(Prefix % 2 == 0);
			using NewType = GenerativeUnit<std::conditional_t<std::is_same_v<Type, f128>, f128, f64>, M / 2, S / 2,
										   Kg / 2, A / 2, K / 2, Mol / 2, Cd / 2, Rad / 2, Sr / 2, Prefix / 2>;
			return NewType{mData >= 0 && mData < std::numeric_limits<double>::infinity()
							   ? sqrtNewtonRaphson(mData, mData, 0)
							   : std::numeric_limits<f64>::quiet_NaN()};
		}

		static consteval bool hasNoPrefix() {
			return Prefix == 0;
		}

		static consteval std::string_view getType() {
			return std::string_view(sTypeText.begin(), sTypeText.size());
		}

	private:
		template<typename _Type = Type, i8 _Prefix = Prefix>
		static constexpr auto scale(const Self<_Type, _Prefix> value) {
			using NewType = decltype(Type() * _Type());

			if constexpr(_Prefix - Prefix > 0) {
				constexpr NewType scaleFactor = pow10(_Prefix - Prefix);
				return value.toRaw() * scaleFactor;
			} else if constexpr(_Prefix - Prefix < 0) {
				constexpr NewType scaleFactor = pow10(Prefix - _Prefix);
				return value.toRaw() / scaleFactor;
			} else
				return value.toRaw();
		}

		static consteval i64 pow10(const u8 pow) {
			if(pow > 0) {
				return 10L * pow10(pow - 1);
			}
			return 1;
		}

		static constexpr f64 sqrtNewtonRaphson(f64 x, f64 curr, f64 prev) {
			return curr == prev ? curr : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
		}

		static consteval auto getTypeArray() {
			constexpr auto withExp = addDimention<Prefix>(toArray(""), toArray(" * 10^"));
			constexpr auto withMet = addDimention<M>(withExp, toArray(" * m^"));
			constexpr auto withSec = addDimention<S>(withMet, toArray(" * s^"));
			constexpr auto withKgm = addDimention<Kg>(withSec, toArray(" * kg^"));
			constexpr auto withAmp = addDimention<A>(withKgm, toArray(" * A^"));
			constexpr auto withKel = addDimention<K>(withAmp, toArray(" * K^"));
			constexpr auto withMol = addDimention<Mol>(withKel, toArray(" * mol^"));
			constexpr auto withCnd = addDimention<Cd>(withMol, toArray(" * cd^"));
			constexpr auto withRad = addDimention<Rad>(withCnd, toArray(" * rad^"));
			constexpr auto withAll = addDimention<Sr>(withRad, toArray(" * sr^"));

			if constexpr(withAll.size() != 0) {
				std::array<char, withAll.size() - 3> result = {};
				for(usize i = 0; i < result.size(); ++i)
					result[i] = withAll[i + 3];
				return result;
			} else {
				constexpr auto result = toArray("dimensionless");
				return result;
			}
		}

		template<i8 Dimention, usize N, usize N1>
		static consteval auto addDimention(const std::array<char, N> &text, const std::array<char, N1> &dimentionText) {
			if constexpr(Dimention > 0) {
				constexpr std::array<char, 2> num = {'0' + Dimention};
				return concatenate(concatenate(text, dimentionText), num);
			} else if constexpr(Dimention < 0) {
				constexpr std::array<char, 2> num = {'-', '0' - Dimention};
				return concatenate(concatenate(text, dimentionText), num);
			} else
				return text;
		}

		template<usize N, usize N1>
		static consteval auto concatenate(const std::array<char, N> &arr1, const std::array<char, N1> &arr2) {
			std::array<char, N + N1> result = {};
			for(usize i = 0; i < N; ++i)
				result[i] = arr1[i];
			for(usize i = 0; i < N1; ++i)
				result[N + i] = arr2[i];

			return result;
		}

		template<usize N>
		static consteval auto toArray(const char (&arr)[N]) {
			std::array<char, N - 1> result = {};
			for(usize i = 0; i < N - 1; ++i)
				result[i] = arr[i];

			return result;
		}

	private:
		inline static constexpr std::array sTypeText = getTypeArray();
		Type mData;
	};

#define GENERATE_SI_UNIT(Name, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift)                 \
	template<typename T = f64>                                                              \
	using Nano##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 9>;  \
	template<typename T = f64>                                                              \
	using Micro##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 6>; \
	template<typename T = f64>                                                              \
	using Milli##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 3>; \
	template<typename T = f64>                                                              \
	using Centi##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 2>; \
	template<typename T = f64>                                                              \
	using Deci##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift - 1>;  \
	template<typename T = f64>                                                              \
	using Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift>;            \
	template<typename T = f64>                                                              \
	using Deca##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 1>;  \
	template<typename T = f64>                                                              \
	using Hecto##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 2>; \
	template<typename T = f64>                                                              \
	using Kilo##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 3>;  \
	template<typename T = f64>                                                              \
	using Mega##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 6>;  \
	template<typename T = f64>                                                              \
	using Giga##Name = GenerativeUnit<T, M, S, Kg, A, K, Mol, Cd, Rad, Sr, BaseShift + 9>;  \
	template<typename T = f64, i8 Power = 0>                                                \
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
	GENERATE_SI_UNIT(Siemenses, -2, 3, -1, 2, 0, 0, 0, 0, 0, 0);
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

	template<typename T = f64>
	using Scale = GenerativeUnit<T, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>;
}