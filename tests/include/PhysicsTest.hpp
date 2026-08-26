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

		static_assert(NStd::DegreesCelsius<>{2} - NStd::DegreesCelsius<>{1} == SI::Kelvins<>{1});
		static_assert(NStd::DegreesFahrenheit<>{2} - NStd::DegreesFahrenheit<>{1} == NStd::DegreesRankin<>{1});

		constexpr auto FAHRENHEIT_NEGATIVE_40 = NStd::DegreesFahrenheit<>{-40}.ToStandardUnit();
		constexpr auto CELSIUS_NEGATIVE_40 = NStd::DegreesCelsius<>{-40}.ToStandardUnit();

		static_assert(FAHRENHEIT_NEGATIVE_40 - CELSIUS_NEGATIVE_40 < SI::Kelvins<>{EPS});
		static_assert(CELSIUS_NEGATIVE_40 - FAHRENHEIT_NEGATIVE_40 < SI::Kelvins<>{EPS});

		static_assert(-(-NStd::DegreesCelsius<>{1}) == NStd::DegreesCelsius<>{1});
	}
}