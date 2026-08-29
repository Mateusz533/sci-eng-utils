#pragma once
//
#include <gtest/gtest.h>
//
#include "SafePhysics/UnitsNStd.hpp"
#include "SafePhysics/UnitsSI.hpp"

namespace PhysicsTest
{
	using namespace Physics;

	constexpr double EPS = 1e-12;

	/************************************************************************************************************/

	TEST(PhysicsTestSuite, OffsettableNStdUnitsTest) {
		using namespace Units;

		static_assert(SI::Kelvins<>{1} + NStd::DegreesCelsius<>{2} == NStd::DegreesCelsius<>{1} + SI::Kelvins<>{2});
		static_assert(NStd::DegreesRankin<>{1} + NStd::DegreesFahrenheit<>{2} == NStd::DegreesFahrenheit<>{1} + NStd::DegreesRankin<>{2});

		static_assert(NStd::DegreesCelsius<>{2} - NStd::DegreesCelsius<>{1} == SI::Kelvins<>{1});
		static_assert(NStd::DegreesFahrenheit<>{2} - NStd::DegreesFahrenheit<>{1} == NStd::DegreesRankin<>{1});

		constexpr auto FAHRENHEIT_NEGATIVE_40 = NStd::DegreesFahrenheit<>{-40}.ToStandardUnit();
		constexpr auto CELSIUS_NEGATIVE_40 = NStd::DegreesCelsius<>{-40}.ToStandardUnit();

		static_assert(FAHRENHEIT_NEGATIVE_40 - CELSIUS_NEGATIVE_40 < SI::Kelvins<>{EPS});
		static_assert(CELSIUS_NEGATIVE_40 - FAHRENHEIT_NEGATIVE_40 < SI::Kelvins<>{EPS});

		static_assert(-(-NStd::DegreesCelsius<>{1}) == NStd::DegreesCelsius<>{1});
	}

	TEST(PhysicsTestSuite, OrdinaryNStdUnitsTest) {
		using namespace Units;

		static_assert(NStd::Inches<>{4} * SI::Meters<>{1} == SI::Meters<>{2} * NStd::Inches<>{2});
		static_assert(NStd::Inches<>{1} * SI::Meters<>{4} == SI::Meters<>{2} * NStd::Inches<>{2});

		constexpr auto LEFT = NStd::Inches<>{1'000} / SI::Meters<>{25.4};
		constexpr auto RIGHT = SI::Meters<>{25.4} / NStd::Inches<>{1'000};
		static_assert(LEFT * RIGHT == SI::Scale<>{1});
		static_assert(LEFT / RIGHT == SI::Scale<>{1});
	}

	TEST(PhysicsTestSuite, NStdUnitRelationsTest) {
		using namespace Units;

		/* Basic units */;

		static_assert(NStd::Percents<>{100}.ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::Permille<>{1000}.ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Degrees<>{180} / SI::Radians<>{std::numbers::pi}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Arcminutes<>{60} / NStd::Degrees<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Arcseconds<>{60} / NStd::Arcminutes<>{1}).ToStandardUnit() == SI::Scale<>{1});

		/* Distance */;

		static_assert((NStd::Yards<>{1} / SI::Meters<>{0.9144}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::NauticalMiles<>{1} / SI::Meters<>{1'852}).ToStandardUnit() == SI::Scale<>{1});

		static_assert((NStd::Inches<>{1} / NStd::Mils<>{1000}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Hands<>{1} / NStd::Inches<>{4}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Links<>{100} / NStd::Chains<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Feet<>{1} / NStd::Inches<>{12}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Yards<>{1} / NStd::Feet<>{3}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Fathoms<>{1} / NStd::Yards<>{2}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Rods<>{4} / NStd::Chains<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Chains<>{1} / NStd::Fathoms<>{11}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Furlongs<>{1} / NStd::Chains<>{10}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::StatuteMiles<>{1} / NStd::Furlongs<>{8}).ToStandardUnit() == SI::Scale<>{1});

		static_assert(NStd::AstronomicalUnits<>{1}.ToStandardUnit() == SI::Meters<>{149'597'870'700});

		/* Surface */;

		static_assert(NStd::SquareInches<>{1} == NStd::Inches<>{1} * NStd::Inches<>{1});
		static_assert(NStd::SquareFeet<>{1} == NStd::Feet<>{1} * NStd::Feet<>{1});
		static_assert(NStd::SquareYards<>{1} == NStd::Yards<>{1} * NStd::Yards<>{1});
		static_assert(NStd::SquareRods<>{1} == NStd::Rods<>{1} * NStd::Rods<>{1});
		static_assert(NStd::SquareChains<>{1} == NStd::Chains<>{1} * NStd::Chains<>{1});
		static_assert(NStd::SquareStatuteMiles<>{1} == NStd::StatuteMiles<>{1} * NStd::StatuteMiles<>{1});
		static_assert((NStd::Roods<>{4} / NStd::Acres<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Acres<>{1} / NStd::SquareYards<>{4'840}).ToStandardUnit() == SI::Scale<>{1});

		static_assert((NStd::Ares<>{1} / SI::SquareMeters<>{100}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Hectares<>{1} / NStd::Ares<>{100}).ToStandardUnit() == SI::Scale<>{1});

		/* Volume */;

		static_assert(NStd::CubicInches<>{1} == NStd::Inches<>{1} * NStd::Inches<>{1} * NStd::Inches<>{1});
		static_assert(NStd::CubicFeet<>{1} == NStd::Feet<>{1} * NStd::Feet<>{1} * NStd::Feet<>{1});
		static_assert(NStd::CubicYards<>{1} == NStd::Yards<>{1} * NStd::Yards<>{1} * NStd::Yards<>{1});
		static_assert(NStd::AcreFeet<>{1} == NStd::Acres<>{1} * NStd::Feet<>{1});
		static_assert(NStd::BoardFeet<>{1} == NStd::SquareFeet<>{1} * NStd::Inches<>{1});
		static_assert((NStd::Cords<>{1} / NStd::CubicFeet<>{128}).ToStandardUnit() == SI::Scale<>{1});

		static_assert((NStd::UsMinims<>{60} / NStd::UsFluidDrams<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsFluidDrams<>{8} / NStd::UsFluidOunces<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsFluidOunces<>{4} / NStd::UsGills<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsGills<>{2} / NStd::UsMeasures<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::UsCups<>{1} == NStd::UsMeasures<>{1});
		static_assert((NStd::UsMeasures<>{2} / NStd::UsLiquidPints<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsLiquidPints<>{2} / NStd::UsLiquidQuarts<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsLiquidQuarts<>{4} / NStd::UsGallons<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsGallons<>{1} / NStd::Liters<>{3.785'411'784}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::PetroleumBarrels<>{1} / NStd::UsGallons<>{42}).ToStandardUnit() == SI::Scale<>{1});

		static_assert((NStd::UsDryPints<>{2} / NStd::UsDryQuarts<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsDryQuarts<>{8} / NStd::UsPecks<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsPecks<>{4} / NStd::UsBushels<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::UsBushels<>{1} / NStd::CubicInches<>{2'150.42}).ToStandardUnit() == SI::Scale<>{1});

		static_assert((NStd::ImperialFluidMinims<>{60} / NStd::ImperialFluidDrachms<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialFluidDrachms<>{8} / NStd::ImperialFluidOunces<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialFluidOunces<>{5} / NStd::ImperialGills<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialGills<>{4} / NStd::ImperialPints<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialPints<>{2} / NStd::ImperialQuarts<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialQuarts<>{4} / NStd::ImperialGallons<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialGallons<>{1} / NStd::Liters<>{4.546'090}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialPecks<>{1} / NStd::ImperialGallons<>{2}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ImperialBushels<>{1} / NStd::ImperialPecks<>{4}).ToStandardUnit() == SI::Scale<>{1});

		static_assert((NStd::Liters<>{1'000} / SI::CubicMeters<>{1}).ToStandardUnit() == SI::Scale<>{1});

		/* Time */;

		static_assert((NStd::Minutes<>{1} / SI::Seconds<>{60}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Hours<>{1} / NStd::Minutes<>{60}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Days<>{1} / NStd::Hours<>{24}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Weeks<>{1} / NStd::Days<>{7}).ToStandardUnit() == SI::Scale<>{1});

		/* Velocity */;

		static_assert(NStd::NauticalMiles<>{1} / NStd::Hours<>{1} == NStd::Knots<>{1});
		static_assert(NStd::StatuteMiles<>{1} / NStd::Hours<>{1} == NStd::MilesPerHour<>{1});
		static_assert(NStd::Feet<>{1} / NStd::Minutes<>{1} == NStd::FeetPerMinute<>{1});
		static_assert(NStd::Feet<>{1} / SI::Seconds<>{1} == NStd::FeetPerSecond<>{1});

		/* Mass & Force */;

		static_assert((NStd::Grains<>{7'000} / NStd::Drams<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Pennyweights<>{1} / NStd::Grains<>{24}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::TroyOunces<>{1} / NStd::Pennyweights<>{20}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::TroyPounds<>{1} / NStd::TroyOunces<>{12}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Drams<>{16} / NStd::Ounces<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Ounces<>{16} / NStd::PoundsMass<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::PoundsMass<>{1} / SI::KiloGrams<>{0.453'592'370}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Stones<>{1} / NStd::PoundsMass<>{14}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Quarters<>{1} / NStd::Stones<>{2}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Centals<>{1} / NStd::PoundsMass<>{100}).ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::ShortHundredweights<>{1} == NStd::Centals<>{1});
		static_assert((NStd::LongHundredweights<>{1} / NStd::Quarters<>{4}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ShortTons<>{1} / NStd::ShortHundredweights<>{20}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::LongTons<>{1} / NStd::LongHundredweights<>{20}).ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::Slug<>{1} == NStd::PoundsForce<>{1} / (NStd::FeetPerSecond<>{1} / SI::Seconds<>{1}));

		static_assert((NStd::Tons<>{1} / SI::KiloGrams<>{1'000}).ToStandardUnit() == SI::Scale<>{1});

		static_assert((NStd::StandardGravity<>{1} / SI::MetersPerSecondSquared<>{9.806'650}).ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::OuncesForce<>{1} == NStd::StandardGravity<>{1} * NStd::Ounces<>{1});
		static_assert(NStd::PoundsForce<>{1} == NStd::StandardGravity<>{1} * NStd::PoundsMass<>{1});
		static_assert((NStd::Kips<>{1} / NStd::PoundsForce<>{1'000}).ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::Poundals<>{1} == NStd::FeetPerSecond<>{1} / SI::Seconds<>{1} * NStd::PoundsMass<>{1});

		/* Torque, Energy, Power */;

		static_assert(NStd::PoundForceFeet<>{1} == NStd::PoundsForce<>{1} * NStd::Feet<>{1} / SI::Radians<>{1});
		static_assert(NStd::PoundForceInches<>{1} == NStd::PoundsForce<>{1} * NStd::Inches<>{1} / SI::Radians<>{1});
		static_assert(NStd::FootPoundsForce<>{1} == NStd::Feet<>{1} * NStd::PoundsForce<>{1});
		static_assert((NStd::ThCalories<>{1} / SI::Joules<>{4.184}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::ItCalories<>{1} / SI::Joules<>{4.1868}).ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::ThBtu<>{1} == NStd::ThCalories<>{1} * (NStd::DegreesRankin<>{1} * NStd::PoundsMass<>{1}) /
											  (SI::Kelvins<>{1} * SI::KiloGrams<>{0.001}));
		static_assert(NStd::ItBtu<>{1} == NStd::ItCalories<>{1} * (NStd::DegreesRankin<>{1} * NStd::PoundsMass<>{1}) /
											  (SI::Kelvins<>{1} * SI::KiloGrams<>{0.001}));
		static_assert((NStd::UsTherms<>{1} / SI::Joules<>{105'480'400}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::TonsOfRefrigeration<>{1} / (NStd::ItBtu<>{12'000} / NStd::Hours<>{1})).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Horsepower<>{1} / (NStd::PoundsForce<>{550} * NStd::FeetPerSecond<>{1})).ToStandardUnit() == SI::Scale<>{1});

		/* Pressure */;

		static_assert((NStd::InchesOfMercury<>{1'000} / SI::Pascals<>{3'386'389}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Bars<>{1} / SI::Pascals<>{100'000}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::MilliBars<>{1'000} / NStd::Bars<>{1}).ToStandardUnit() == SI::Scale<>{1});
		static_assert((NStd::Atmospheres<>{1} / SI::Pascals<>{101'325}).ToStandardUnit() == SI::Scale<>{1});
		static_assert(NStd::PoundsPerSquareFoot<>{1} == NStd::PoundsForce<>{1} / NStd::SquareFeet<>{1});
		static_assert(NStd::PoundsPerSquareInch<>{1} == NStd::PoundsForce<>{1} / NStd::SquareInches<>{1});
		static_assert(NStd::KipsPerSquareInch<>{1} == NStd::Kips<>{1} / NStd::SquareInches<>{1});

		/* Temperature */;

		static_assert((NStd::DegreesRankin<>{9} / SI::Kelvins<>{5}).ToStandardUnit() == SI::Scale<>{1});
	}
}