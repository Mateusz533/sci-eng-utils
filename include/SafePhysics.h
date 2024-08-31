#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <type_traits>

namespace Physics
{
	using usize = size_t;
	using u8 = u_int8_t;
	using u16 = u_int16_t;
	using u32 = u_int32_t;
	using u64 = u_int64_t;
	using i8 = int8_t;
	using i16 = int16_t;
	using i32 = int32_t;
	using i64 = int64_t;
	using f32 = float;
	using f64 = double;
	using f128 = long double;

	namespace Units
	{
		namespace SI
		{
			template<typename Type, i8 M, i8 S, i8 Kg, i8 A, i8 K, i8 Mol, i8 Cd, i8 Rad, i8 Sr, i8 Prefix>
				requires(std::is_arithmetic_v<Type> && !!(Prefix == 0 || (M || S || Kg || A || K || Mol || Cd || Rad || Sr)))
			class GenerativeUnit
			{
			public:
				template<typename _Type = Type, i8 _Prefix = Prefix>
				using Self = GenerativeUnit<_Type, M, S, Kg, A, K, Mol, Cd, Rad, Sr, _Prefix>;

				// Constructors

				constexpr GenerativeUnit(const Self<> &value) : mData(value.mData) {}
				constexpr GenerativeUnit(const Self<> &&value) : mData(value.mData) {}
				constexpr GenerativeUnit() : mData(0) {}
				constexpr GenerativeUnit(const Type data) : mData(data) {}
				template<typename _Type = Type, i8 _Prefix = Prefix>
				constexpr GenerativeUnit(const Self<_Type, _Prefix> value) : mData(scale(value)) {}

				// Assignment operators

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

				constexpr auto operator-() const {
					return Self<decltype(-Type())>{-mData};
				}
				template<typename _Type = Type>
				constexpr auto operator+(const Self<_Type> value) const {
					using NewType = decltype(Type() * _Type());
					return Self<NewType>{mData + value.toRaw()};
				}
				template<typename _Type = Type>
				constexpr auto operator-(const Self<_Type> value) const {
					using NewType = decltype(Type() * _Type());
					return Self<NewType>{mData - value.toRaw()};
				}
				template<typename _Type, i8 _M, i8 _S, i8 _Kg, i8 _A, i8 _K, i8 _Mol, i8 _Cd, i8 _Rad, i8 _Sr, i8 _Prefix>
				constexpr auto operator*(const GenerativeUnit<_Type, _M, _S, _Kg, _A, _K, _Mol, _Cd, _Rad, _Sr, _Prefix> value) const {
					return GenerativeUnit<decltype(Type() * _Type()), M + _M, S + _S, Kg + _Kg, A + _A, K + _K, Mol + _Mol, Cd + _Cd,
										  Rad + _Rad, Sr + _Sr, Prefix + _Prefix>{mData * value.toRaw()};
				}
				template<typename _Type, i8 _M, i8 _S, i8 _Kg, i8 _A, i8 _K, i8 _Mol, i8 _Cd, i8 _Rad, i8 _Sr, i8 _Prefix>
				constexpr auto operator/(const GenerativeUnit<_Type, _M, _S, _Kg, _A, _K, _Mol, _Cd, _Rad, _Sr, _Prefix> value) const {
					return GenerativeUnit<decltype(Type() * _Type()), M - _M, S - _S, Kg - _Kg, A - _A, K - _K, Mol - _Mol, Cd - _Cd,
										  Rad - _Rad, Sr - _Sr, Prefix - _Prefix>{mData / value.toRaw()};
				}

				friend std::ostream &operator<<(std::ostream &os, const Self<> &obj) {
					os << obj.mData << ' ';
					if constexpr(Prefix != 0)
						os << "* ";
					os << obj.getType();
					return os;
				}

				// Other methods

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

				static std::string getType() {
					constexpr std::array type = getTypeArray();
					constexpr auto size = type.size();
					return std::string(type.begin(), size);
				}

				static consteval bool hasNoPrefix() {
					return Prefix == 0;
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
				static constexpr i64 pow10(const u8 pow) {
					if(pow > 0) {
						return 10L * pow10(pow - 1);
					}
					return 1;
				}
				static constexpr f64 sqrtNewtonRaphson(f64 x, f64 curr, f64 prev) {
					return curr == prev ? curr : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
				}

				template<usize N>
				static consteval auto toArray(const char (&arr)[N]) {
					std::array<char, N - 1> result = {};
					for(usize i = 0; i < N - 1; ++i)
						result[i] = arr[i];

					return result;
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

				static constexpr auto getTypeArray() {
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

			// Basic units
			GENERATE_SI_UNIT(Meters, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Seconds, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Grams, 0, 0, 1, 0, 0, 0, 0, 0, 0, -3);
			GENERATE_SI_UNIT(Amperes, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Kelvins, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Moles, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Candelas, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0);
			GENERATE_SI_UNIT(Radians, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0);
			GENERATE_SI_UNIT(Steradians, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0);
			// Kinematics
			GENERATE_SI_UNIT(MetersPerSecond, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(RadiansPerSecond, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0);
			GENERATE_SI_UNIT(MetersPerSecondSquare, 1, -2, 0, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(RadiansPerSecondSquare, 0, -2, 0, 0, 0, 0, 0, 1, 0, 0);
			GENERATE_SI_UNIT(MetersPerRadian, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0);
			GENERATE_SI_UNIT(RadialMeters, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0);
			// Dynamics
			GENERATE_SI_UNIT(GramMeters, 1, 0, 1, 0, 0, 0, 0, -2, 0, -3);
			GENERATE_SI_UNIT(GramMetersSquare, 2, 0, 1, 0, 0, 0, 0, -2, 0, -3);
			GENERATE_SI_UNIT(NewtonSeconds, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(NewtonMeterSeconds, 2, -1, 1, 0, 0, 0, 0, -1, 0, 0);
			GENERATE_SI_UNIT(Newtons, 1, -2, 1, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(NewtonMeters, 2, -2, 1, 0, 0, 0, 0, -1, 0, 0);
			GENERATE_SI_UNIT(Joules, 2, -2, 1, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Watts, 2, -3, 1, 0, 0, 0, 0, 0, 0, 0);
			// Electrics
			GENERATE_SI_UNIT(Coulombs, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(NewtonsPerCoulomb, 1, -3, 1, -1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(VoltsPerMeter, 1, -3, 1, -1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Volts, 2, -3, 1, -1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Ohms, 2, -3, 1, -2, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Siemenses, -2, 3, -1, 2, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Farads, -2, 4, -1, 2, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(FaradsPerMeter, -3, 4, -1, 2, 0, 0, 0, 0, 0, 0);
			// Magnetism
			GENERATE_SI_UNIT(AmperesPerMeter, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Teslas, 0, -2, 1, -1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Webers, 2, -2, 1, -1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(TeslaMeters, 1, -2, 1, -1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(Henries, 2, -2, 1, -2, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(HenriesPerMeter, 1, -2, 1, -2, 0, 0, 0, 0, 0, 0);
			// Materials
			GENERATE_SI_UNIT(MetersSquare, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(MetersCube, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(GramsPerMeterSquare, -2, 0, 1, 0, 0, 0, 0, 0, 0, -3);
			GENERATE_SI_UNIT(GramsPerMeterCube, -3, 0, 1, 0, 0, 0, 0, 0, 0, -3);
			GENERATE_SI_UNIT(Pascals, -1, -2, 1, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(CoulombsPerMeterCube, -3, 1, 0, 1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(AmperesPerMeterSquare, -2, 0, 0, 1, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(GramsPerMole, 0, 0, 1, 0, 0, -1, 0, 0, 0, -3);
			GENERATE_SI_UNIT(PartsPerMole, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0);
			// Vibrations and waves
			GENERATE_SI_UNIT(Hertzes, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(RadiansPerMeter, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0);
			// Thermodynamics
			GENERATE_SI_UNIT(WattsPerMeterSquare, 0, -3, 1, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(JoulesPerKelvin, 2, -2, 1, 0, -1, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(JoulesPerKilogramKelvin, 2, -2, 0, 0, -1, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(WattsPerMeterKelvin, 1, -3, 1, 0, -1, 0, 0, 0, 0, 0);
			// Optics
			GENERATE_SI_UNIT(Lumens, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0);
			GENERATE_SI_UNIT(Luxes, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0);
			// Quantum mechanics
			GENERATE_SI_UNIT(JouleSeconds, 2, -1, 1, 0, 0, 0, 0, 0, 0, 0);
			GENERATE_SI_UNIT(JouleSecondsPerRadian, 2, -1, 1, 0, 0, 0, 0, -1, 0, 0);

			template<typename T = f64>
			using Scale = GenerativeUnit<T, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0>;
		}

		namespace NStd
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
				template<typename SiblingUnit>
					requires(std::is_base_of_v<_Base, SiblingUnit>)
				static consteval bool hasSameStdUnitBase() {
					return std::is_same_v<decltype(SiblingUnit().toStandardUnit() / StandardUnit<decltype(SiblingUnit().toRaw())>()),
										  SI::Scale<decltype(SiblingUnit().toRaw())>>;
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
	}

	namespace Constants
	{
		using namespace Units::SI;
		constexpr inline Radians<> RADIAN_PER_REVOLUTION{2 * M_PI};
		constexpr inline MetersPerSecond<> SPEED_OF_LIGHT{2.99792458e8};
		constexpr inline JouleSeconds<> PLANCK_CONSTANT{6.62607015e-34};
		constexpr inline Coulombs<> ELEMENTARY_CHARGE{1.602176634e-19};
		constexpr inline JoulesPerKelvin<> BOLTZMANN_CONSTANT{1.380649e-23};
		constexpr inline PartsPerMole<> AVOGADRO_CONSTANT{6.02214076e23};

		constexpr inline JouleSecondsPerRadian<> DIRAC_CONSTANT{PLANCK_CONSTANT / RADIAN_PER_REVOLUTION};
		constexpr inline HenriesPerMeter<> MAGNETIC_CONSTANT{4e-7 * M_PI};
		constexpr inline FaradsPerMeter<> ELECTRIC_CONSTANT = Scale<>{1} / (MAGNETIC_CONSTANT * SPEED_OF_LIGHT.power<2>());
	}

	namespace Calculate
	{
		using namespace Units::SI;
		// mechanics
		template<typename T>
		constexpr inline MetersPerSecond<T> tangentialVelocity(RadiansPerSecond<T> angularVelocity, RadialMeters<T> radius) {
			return angularVelocity * radius;
		}
		template<typename T>
		constexpr inline MetersPerSecondSquare<T> tangentialAcceleration(RadiansPerSecondSquare<T> angularAcceleration, RadialMeters<T> radius) {
			return angularAcceleration * radius;
		}
		template<typename T>
		constexpr inline NewtonMeterSeconds<T> angularMomentum(NewtonSeconds<T> momentum, RadialMeters<T> radius) {
			return momentum * radius;
		}
		template<typename T>
		constexpr inline Joules<T> energy(KiloGrams<T> mass, MetersPerSecond<T> velocity) {
			return mass * velocity * velocity / Scale<T>(2);
		}
		template<typename T>
		constexpr inline Joules<T> energy(KiloGramMetersSquare<T> momentOfInertia, RadiansPerSecond<T> angularVelocity) {
			return momentOfInertia * angularVelocity * angularVelocity / Scale<T>(2);
		}
		template<typename T>
		constexpr inline Watts<T> power(Newtons<T> force, MetersPerSecond<T> velocity) {
			return force * velocity;
		}
		template<typename T>
		constexpr inline Watts<T> power(NewtonMeters<T> torque, RadiansPerSecond<T> angularVelocity) {
			return torque * angularVelocity;
		}
		template<typename T>
		constexpr inline KiloGramMetersSquare<T> inertia(KiloGrams<T> mass, RadialMeters<T> radius, Scale<T> factor) {
			return factor * mass * radius * radius;
		}
		// electronics
		template<typename T>
		constexpr inline Amperes<T> current(Coulombs<T> charge, Seconds<T> time) {
			return charge / time;
		}
		template<typename T>
		constexpr inline Coulombs<T> charge(Amperes<T> current, Seconds<T> time) {
			return current * time;
		}
		template<typename T>
		constexpr inline Ohms<T> resistance(Volts<T> voltage, Amperes<T> current) {
			return voltage / current;
		}
		template<typename T>
		constexpr inline Coulombs<T> charge(Farads<T> capacity, Volts<T> voltage) {
			return capacity * voltage;
		}
		template<typename T>
		constexpr inline Webers<T> magneticFlux(Henries<T> inductance, Amperes<T> current) {
			return inductance * current;
		}
		// oscilations and waves
		template<typename T>
		constexpr inline RadiansPerMeter<T> waveNumber(Meters<T> waveLength) {
			return Constants::RADIAN_PER_REVOLUTION / waveLength;
		}
		template<typename T>
		constexpr inline Seconds<T> period(Hertzes<T> frequency) {
			return Scale<T>(1) / frequency;
		}
		// relativistic mechanics
		template<typename T>
		constexpr inline Scale<T> timeDilation(MetersPerSecond<T> velocity) {
			return Scale<T>(1) / (Scale<T>(1) - (velocity / Constants::SPEED_OF_LIGHT).template power<2>()).sqrt();
		}

		// TODO: Add more patterns
	}
}