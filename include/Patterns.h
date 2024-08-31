#pragma once

#include "UnitsSI.h"

namespace Physics
{
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