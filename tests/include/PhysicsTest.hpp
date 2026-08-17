#pragma once
//
#include <gtest/gtest.h>
//
#include "SafePhysics/SafePhysics.hpp"
#include "SafePhysics/UnitsNStd.hpp"
#include "SafePhysics/UnitsSI.hpp"

namespace PhysicsTest
{
	using namespace Physics;

	/************************************************************************************************************/

	TEST(PhysicsTestSuite, PhysicsTest) {
		using namespace Units;
		static_assert(NStd::DegreesCelsius<>{2} - NStd::DegreesCelsius<>{1} == SI::Kelvins<>{1}, "");
	}
}