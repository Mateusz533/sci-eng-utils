#include "SafePhysics/Patterns.hpp"
#include "SafePhysics/UnitsNStd.hpp"

using namespace Physics;
using namespace Units;

inline SI::KiloGramMetersSquare<f64> cylinderInertia(SI::RadialMeters<f64> radius, SI::KiloGrams<f64> mass) {
	return SI::Scale<f64>(0.5) * mass * radius * radius;
}

/// @return Moment of inertia in kg * m^2
inline double cylinderInertia(double radiusInMeters, double massInKilograms) {
	return 0.5 * massInKilograms * radiusInMeters * radiusInMeters;
}

inline double cylinderInertia2(double radius, double mass) {
	return 0.5 * mass * radius * radius;
}

inline SI::Meters<f64> distance(SI::MetersPerSecondSquare<f64> acceleration, SI::Seconds<f64> time, SI::MetersPerSecond<f64> startVelocity = 0) {
	return (startVelocity + SI::Scale<f64>(0.5) * acceleration * time) * time;
}

int main() {
	/*----------------------------------------------------------------------*/
	/*------------------------ Standard Units tests ------------------------*/
	/*----------------------------------------------------------------------*/
	std::cout << "\n--- Standard units tests ---\n\n";

	constexpr SI::Meters<> distance(1);
	constexpr SI::Seconds<> time(0);
	static_assert(!distance == false);
	static_assert(!time);
	constexpr SI::KiloGrams<u16> mass(1);
	static_assert((-mass).toRaw() < 0);

	constexpr auto ingMM = SI::MilliMeters<>(1);
	constexpr auto ingM = SI::Meters<>(1);
	auto sumM = SI::Meters<>(1);
	auto sumMM = SI::MilliMeters<>(1);
	// 1 m += 1 mm = 1.001 m
	sumM += ingMM;
	// 1 mm += 1 m = 1001 mm
	sumMM += ingM;
	// 1 m + 1 mm = ERROR
	// ingM + ingMM;
	std::cout << "1.001 m: " << sumM << std::endl;
	std::cout << "1001 mm: " << sumMM << std::endl;

	static_assert(SI::Meters<int>{5} / SI::Seconds<int>{2} == SI::MetersPerSecond<int>{2});
	static_assert(SI::NewtonSeconds<double>{5} / SI::Seconds<int>{2} == SI::Newtons<double>{2.5});
	static_assert(SI::MetersSquare<>{4}.sqrt() == SI::Meters<>{2});
	static_assert(std::isnan(SI::MetersSquare<>{-1}.sqrt().toRaw()));
	// 5 m / 0 s = ERROR
	// SI::Meters<int>{5} / SI::Seconds<int>{0}

	std::cout << "Scale type: " << SI::Scale<>().getType() << std::endl;
	std::cout << "Meters type: " << SI::Meters<>().getType() << std::endl;
	std::cout << "NanoWatts type: " << SI::NanoWatts<>().getType() << std::endl;
	std::cout << "Electric constant value: " << Constants::ELECTRIC_CONSTANT << std::endl;

	// Patterns
	constexpr SI::MetersPerSecond<f128> _99percentSpeedOfLight{SI::Scale<>{0.99} * Constants::SPEED_OF_LIGHT};
	std::cout << "Gamma for 99 \% speed of light: " << Calculate::timeDilation(_99percentSpeedOfLight) << std::endl;

	/*----------------------------------------------------------------------*/
	/*-------------------------- NStd Units tests --------------------------*/
	/*----------------------------------------------------------------------*/
	std::cout << "\n--- NStd units tests ---\n\n";

	NStd::Inch<double> resultInch = ingM;
	std::cout << "1 meter in inches: " << resultInch << std::endl;
	std::cout << "1 meter after double conversion: " << resultInch.toStandardUnit<>() << std::endl;

	NStd::Degree<> degree;
	std::cout << "Offset Numerator:   " << degree.getOffsetNumerator() << std::endl;
	std::cout << "Offset Denominator: " << degree.getOffsetDenominator() << std::endl;
	std::cout << "Scale Numerator:    " << degree.getScaleNumerator() << std::endl;
	std::cout << "Scale Denominator:  " << degree.getScaleDenominator() << std::endl;
	std::cout << "Degree per radian:  " << 1.0 * degree.getScaleNumerator() / degree.getScaleDenominator() << std::endl;
	std::cout << "Radian per degree:  " << 1.0 * degree.getScaleDenominator() / degree.getScaleNumerator() << std::endl;

	std::cout << "  0 *F = " << NStd::DegreeFahrenheit<>{0}.toStandardUnit() << std::endl;
	std::cout << "100 *F = " << NStd::DegreeFahrenheit<>{100}.toStandardUnit() << std::endl;
	std::cout << "  0 *C = " << NStd::DegreeCelsius<>{0}.toStandardUnit() << std::endl;
	std::cout << "100 *C = " << NStd::DegreeCelsius<>{100}.toStandardUnit() << std::endl;

	std::cout << "\n--- The end of the tests ---\n\n";

	return 0;
}