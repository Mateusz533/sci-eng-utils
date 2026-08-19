#include "SafePhysics/Patterns.hpp"
#include "SafePhysics/UnitsNStd.hpp"

using namespace Physics;
using namespace Units;

constexpr SI::KiloGramMetersSquare<f64> CylinderInertia(SI::RadialMeters<f64> radius, SI::KiloGrams<f64> mass) noexcept {
	return SI::Scale<f64>{0.5} * mass * radius * radius;
}

constexpr SI::Meters<f64> Distance(SI::MetersPerSecondSquare<f64> acceleration, SI::Seconds<f64> time, SI::MetersPerSecond<f64> startVelocity = 0) noexcept {
	return (startVelocity + SI::Scale<f64>{0.5} * acceleration * time) * time;
}

int main() {
	/*----------------------------------------------------------------------*/
	/*------------------------ Standard Units tests ------------------------*/
	/*----------------------------------------------------------------------*/
	std::cout << "\n--- Standard units tests ---\n\n";

	constexpr SI::Meters<> DISTANCE{1};
	constexpr SI::Seconds<> TIME{0};
	static_assert(!DISTANCE == false);
	static_assert(!TIME);
	constexpr SI::KiloGrams<u16> MASS{1};
	static_assert((-MASS).ToRaw() < 0);

	constexpr auto ING_MM = SI::MilliMeters<>{1};
	constexpr auto ING_M = SI::Meters<>{1};
	auto sumM = SI::Meters<>{1};
	auto sumMM = SI::MilliMeters<>{1};
	// 1 m += 1 mm = 1.001 m
	sumM += ING_MM;
	// 1 mm += 1 m = 1001 mm
	sumMM += ING_M;
	// 1 m + 1 mm = ERROR
	// ingM + ingMM;
	std::cout << "1.001 m: " << sumM << std::endl;
	std::cout << "1001 mm: " << sumMM << std::endl;

	static_assert(SI::Meters<int>{5} / SI::Seconds<int>{2} == SI::MetersPerSecond<int>{2});
	static_assert(SI::NewtonSeconds<double>{5} / SI::Seconds<int>{2} == SI::Newtons<double>{2.5});
	static_assert(SI::MetersSquare<>{4}.Sqrt() == SI::Meters<>{2});
	static_assert(std::isnan(SI::MetersSquare<>{-1}.Sqrt().ToRaw()));
	// 5 m / 0 s = ERROR
	// SI::Meters<int>{5} / SI::Seconds<int>{0}

	std::cout << "Scale type: " << SI::Scale<>().GetTypeView() << std::endl;
	std::cout << "Meters type: " << SI::Meters<>().GetTypeView() << std::endl;
	std::cout << "NanoWatts type: " << SI::NanoWatts<>().GetTypeView() << std::endl;
	std::cout << "Electric constant value: " << Constants::ELECTRIC_CONSTANT << std::endl;

	// Patterns
	constexpr SI::MetersPerSecond<f128> _99percentSpeedOfLight{SI::Scale<>{0.99} * Constants::SPEED_OF_LIGHT};
	std::cout << "Gamma for 99 \% speed of light: " << Calculate::TimeDilation(_99percentSpeedOfLight) << std::endl;

	/*----------------------------------------------------------------------*/
	/*-------------------------- NStd Units tests --------------------------*/
	/*----------------------------------------------------------------------*/
	std::cout << "\n--- NStd units tests ---\n\n";

	NStd::Inches<double> resultInch = ING_M;
	std::cout << "1 meter in inches: " << resultInch << std::endl;
	std::cout << "1 meter after double conversion: " << resultInch.ToStandardUnit<>() << std::endl;

	NStd::Degrees<> degree{1};
	std::cout << "Degree per radian:  " << degree.ToStandardUnit().ToRaw() / degree.ToRaw() << std::endl;
	std::cout << "Radian per degree:  " << degree.ToRaw() / degree.ToStandardUnit().ToRaw() << std::endl;

	std::cout << "  0 *F = " << NStd::DegreesFahrenheit<>{0}.ToStandardUnit() << std::endl;
	std::cout << "100 *F = " << NStd::DegreesFahrenheit<>{100}.ToStandardUnit() << std::endl;
	std::cout << "  0 *C = " << NStd::DegreesCelsius<>{0}.ToStandardUnit() << std::endl;
	std::cout << "100 *C = " << NStd::DegreesCelsius<>{100}.ToStandardUnit() << std::endl;

	std::cout << "\n--- The end of the tests ---\n\n";

	return 0;
}