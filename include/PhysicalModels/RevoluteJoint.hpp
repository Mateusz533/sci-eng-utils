#pragma once
//
#include <algorithm>
//
#include "SafePhysics/UnitsSI.hpp"

namespace PhysicalModels
{
	using namespace Physics::Units;

	class RevoluteJoint
	{
	public:
		constexpr RevoluteJoint() = default;
		constexpr RevoluteJoint(SI::NewtonMeters<> staticFriction) : staticFriction{staticFriction} {};

		constexpr void SetStaticFriction(SI::NewtonMeters<> newStaticFriction) {
			staticFriction = newStaticFriction;
		}
		constexpr SI::Radians<> GetAngle() const {
			return angle;
		}
		constexpr SI::RadiansPerSecond<> GetRate() const {
			return angularRate;
		}
		constexpr void Step(SI::NewtonMeters<> activeTorque, SI::RadiansPerSecondSquared<> baseLinkAcc,
							SI::KiloGramMetersSquared<> drivenLinkInertia, SI::Seconds<> timeStep) {
			const SI::RadiansPerSecondSquared<> externalAcc = activeTorque / drivenLinkInertia - baseLinkAcc;
			const SI::RadiansPerSecondSquared<> maxFrictionAcc = staticFriction / drivenLinkInertia;

			if(angularRate != SI::RadiansPerSecond<>{0.0}) {
				const SI::RadiansPerSecondSquared<> totalAcc = externalAcc - Sign(angularRate) * maxFrictionAcc;
				const bool mayChangeDirection = (Sign(angularRate) != Sign(angularRate + totalAcc * timeStep));
				const SI::Seconds<> timeStep0 = mayChangeDirection ? -angularRate / totalAcc : timeStep;

				angle += (angularRate + totalAcc * timeStep0 / SI::Scale<>{2}) * timeStep0;
				angularRate += totalAcc * timeStep0;
				timeStep -= timeStep0;
			}

			const SI::RadiansPerSecondSquared<> totalAcc = externalAcc - std::clamp(externalAcc, -maxFrictionAcc, maxFrictionAcc);

			angle += totalAcc * timeStep * timeStep / SI::Scale<>{2};
			angularRate += totalAcc * timeStep;
		}

		template<typename T>
		static constexpr SI::Scale<> Sign(T value) {
			return (value < T{0}) ? -1 : 1;
		}

	private:
		SI::NewtonMeters<> staticFriction = 0.0;
		SI::Radians<> angle = 0.0;
		SI::RadiansPerSecond<> angularRate = 0.0;
	};
}