#pragma once
//
#include <gtest/gtest.h>
//
#include "SafePhysics/UnitsNStd.hpp"
#include "SafePhysics/UnitsSI.hpp"

namespace PhysicsTest
{
	using namespace Physics;

	/************************************************************************************************************/

	TEST(PhysicsTestSuite, PhysicsTest) {
		using namespace Units;
		static_assert(NStd::DegreesCelsius<>{2} - NStd::DegreesCelsius<>{1} == SI::Kelvins<>{1}, "");
		static_assert(NStd::DegreesFahrenheit<>{-40}.ToStandardUnit() - NStd::DegreesCelsius<>{-40}.ToStandardUnit() < SI::Kelvins<>{1e-13}, "");
	}
}