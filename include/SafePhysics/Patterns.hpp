#pragma once
//
#include <numbers>
//
#include "UnitsSI.hpp"
#include "Utils/CompileTime/Math.hpp"
#include "Utils/Types.hpp"

namespace Physics
{
	namespace Constants
	{
		using namespace Units::SI;

		inline constexpr Radians<> RADIAN_PER_REVOLUTION{2.0 * std::numbers::pi};
		inline constexpr Hertzes<> CESIUM_133_HYPERFINE_TRANSITION_FREQUENCY{9.19263177e9};
		inline constexpr MetersPerSecond<> SPEED_OF_LIGHT{2.99792458e8};
		inline constexpr JouleSeconds<> PLANCK_CONSTANT{6.62607015e-34};
		inline constexpr Coulombs<> ELEMENTARY_CHARGE{1.602176634e-19};
		inline constexpr JoulesPerKelvin<> BOLTZMANN_CONSTANT{1.380649e-23};
		inline constexpr PartsPerMole<> AVOGADRO_CONSTANT{6.02214076e23};
		inline constexpr LumensPerWatt<> LUMINOUS_EFFICACY_FOR_540_THZ{683};

		inline constexpr JouleSecondsPerRadian<> DIRAC_CONSTANT = PLANCK_CONSTANT / RADIAN_PER_REVOLUTION;
		inline constexpr JoulesPerMoleKelvin<> UNIVERSAL_GAS_CONSTANT = AVOGADRO_CONSTANT * BOLTZMANN_CONSTANT;
		inline constexpr CoulombsPerMole<> FARADAY_CONSTANT = AVOGADRO_CONSTANT * ELEMENTARY_CHARGE;
		inline constexpr auto STEFAN_BOLTZMANN_CONSTANT = Scale<>{2.0 * Utils::CompileTime::Pow<f64, 5>(std::numbers::pi) / 15.0} *
														  BOLTZMANN_CONSTANT.Power<4>() / (PLANCK_CONSTANT.Power<3>() * SPEED_OF_LIGHT.Power<2>());
		inline constexpr HertzesPerVolt<> JOSEPHSON_CONSTANT = Scale<>{2} * ELEMENTARY_CHARGE / PLANCK_CONSTANT;
		inline constexpr Ohms<> VON_KLITZING_CONSTANT = PLANCK_CONSTANT / ELEMENTARY_CHARGE.Power<2>();

		inline constexpr CubicMetersPerKilogramSecondSquared<> GRAVITATIONAL_CONSTANT{6.6743e-11};
		inline constexpr Scale<> FINE_STRUCTURE_CONSTANT{7.2973525643e-3};
		inline constexpr HenriesPerMeter<> MAGNETIC_CONSTANT = Scale<>{2} * PLANCK_CONSTANT * FINE_STRUCTURE_CONSTANT /
															   (ELEMENTARY_CHARGE.Power<2>() * SPEED_OF_LIGHT);
		inline constexpr FaradsPerMeter<> ELECTRIC_CONSTANT = Scale<>{1} / (MAGNETIC_CONSTANT * SPEED_OF_LIGHT.Power<2>());
	}

	namespace Calculate
	{
		using namespace Units::SI;

		/* Mechanics */;

		template<Arithmetic T>
		constexpr MetersPerSecond<T> TangentialVelocity(RadiansPerSecond<T> angularVelocity, RadialMeters<T> radius) noexcept {
			return angularVelocity * radius;
		}
		template<Arithmetic T>
		constexpr MetersPerSecondSquared<T> TangentialAcceleration(RadiansPerSecondSquared<T> angularAcceleration, RadialMeters<T> radius) noexcept {
			return angularAcceleration * radius;
		}
		template<Arithmetic T>
		constexpr KiloGramMetersSquared<T> MomentOfInertia(KiloGrams<T> mass, RadialMeters<T> radius, Scale<T> factor) noexcept {
			return factor * mass * radius * radius;
		}
		template<Arithmetic T>
		constexpr NewtonSeconds<T> Momentum(KiloGrams<T> mass, MetersPerSecond<T> velocity) noexcept {
			return mass * velocity;
		}
		template<Arithmetic T>
		constexpr NewtonMeterSeconds<T> AngularMomentum(NewtonSeconds<T> momentum, RadialMeters<T> radius) noexcept {
			return momentum * radius;
		}
		template<Arithmetic T>
		constexpr NewtonMeterSeconds<T> AngularMomentum(KiloGramMetersSquared<T> momentOfInertia, RadiansPerSecond<T> angularVelocity) noexcept {
			return momentOfInertia * angularVelocity;
		}
		template<Arithmetic T>
		constexpr Newtons<T> Force(KiloGrams<T> mass, MetersPerSecondSquared<T> acceleration) noexcept {
			return mass * acceleration;
		}
		template<Arithmetic T>
		constexpr NewtonMeters<T> Torque(Newtons<T> force, RadialMeters<T> radius) noexcept {
			return force * radius;
		}
		template<Arithmetic T>
		constexpr NewtonMeters<T> Torque(KiloGramMetersSquared<T> momentOfInertia, MetersPerSecondSquared<T> angularAcceleration) noexcept {
			return momentOfInertia * angularAcceleration;
		}
		template<Arithmetic T>
		constexpr Joules<T> Work(Newtons<T> force, Meters<T> shift) noexcept {
			return force * shift;
		}
		template<Arithmetic T>
		constexpr Joules<T> Work(NewtonMeters<T> torque, Radians<T> angularShift) noexcept {
			return torque * angularShift;
		}
		template<Arithmetic T>
		constexpr Joules<T> KineticEnergy(KiloGrams<T> mass, MetersPerSecond<T> velocity) noexcept {
			return mass * velocity * velocity / Scale<T>{2};
		}
		template<Arithmetic T>
		constexpr Joules<T> KineticEnergy(KiloGramMetersSquared<T> momentOfInertia, RadiansPerSecond<T> angularVelocity) noexcept {
			return momentOfInertia * angularVelocity * angularVelocity / Scale<T>{2};
		}
		template<Arithmetic T>
		constexpr Watts<T> Power(Newtons<T> force, MetersPerSecond<T> velocity) noexcept {
			return force * velocity;
		}
		template<Arithmetic T>
		constexpr Watts<T> Power(NewtonMeters<T> torque, RadiansPerSecond<T> angularVelocity) noexcept {
			return torque * angularVelocity;
		}

		/* Electromagnetism */;

		template<Arithmetic T>
		constexpr Amperes<T> Current(Coulombs<T> charge, Seconds<T> time) noexcept {
			return charge / time;
		}
		template<Arithmetic T>
		constexpr Coulombs<T> Charge(Amperes<T> current, Seconds<T> time) noexcept {
			return current * time;
		}
		template<Arithmetic T>
		constexpr Ohms<T> Resistance(Volts<T> voltage, Amperes<T> current) noexcept {
			return voltage / current;
		}
		template<Arithmetic T>
		constexpr Volts<T> Voltage(Ohms<T> resistance, Amperes<T> current) noexcept {
			return resistance * current;
		}
		template<Arithmetic T>
		constexpr Coulombs<T> Charge(Farads<T> capacity, Volts<T> voltage) noexcept {
			return capacity * voltage;
		}
		template<Arithmetic T>
		constexpr Farads<T> Capacity(Coulombs<T> charge, Volts<T> voltage) noexcept {
			return charge / voltage;
		}
		template<Arithmetic T>
		constexpr Webers<T> MagneticFlux(Henries<T> inductance, Amperes<T> current) noexcept {
			return inductance * current;
		}
		template<Arithmetic T>
		constexpr Henries<T> Inductance(Webers<T> magneticFlux, Amperes<T> current) noexcept {
			return magneticFlux / current;
		}
		template<Arithmetic T>
		constexpr Volts<T> ElectromotiveForce(Henries<T> inductance, AmperesPerSecond<T> currentRate) noexcept {
			return inductance * currentRate;
		}

		/* Vibrations and waves */;

		template<Arithmetic T>
		constexpr RadiansPerMeter<T> WaveNumber(Meters<T> waveLength) noexcept {
			return Constants::RADIAN_PER_REVOLUTION / waveLength;
		}
		template<Arithmetic T>
		constexpr Meters<T> WaveLength(RadiansPerMeter<T> waveNumber) noexcept {
			return Constants::RADIAN_PER_REVOLUTION / waveNumber;
		}
		template<Arithmetic T>
		constexpr Seconds<T> Period(Hertzes<T> frequency) noexcept {
			return Scale<T>{1} / frequency;
		}
		template<Arithmetic T>
		constexpr Hertzes<T> Frequency(Seconds<T> period) noexcept {
			return Scale<T>{1} / period;
		}

		/* Relativistic mechanics */;

		template<Arithmetic T>
		constexpr Scale<T> LorentzFactor(MetersPerSecond<T> velocity) noexcept {
			return Scale<T>{1} / (Scale<T>{1} - (velocity / Constants::SPEED_OF_LIGHT).template Power<2>()).Sqrt();
		}
		template<Arithmetic T>
		constexpr NewtonSeconds<T> RelativisticMomentum(KiloGrams<T> mass, MetersPerSecond<T> velocity) noexcept {
			return mass * velocity * LorentzFactor(velocity);
		}
		template<Arithmetic T>
		constexpr Joules<T> RelativisticEnergy(KiloGrams<T> mass, MetersPerSecond<T> velocity) noexcept {
			return mass * Constants::SPEED_OF_LIGHT.Power<2>() * LorentzFactor(velocity);
		}

		/* Quantum mechanics */;

		template<Arithmetic T>
		constexpr NewtonSeconds<T> ParticleMomentum(Meters<T> waveLength) noexcept {
			return Constants::PLANCK_CONSTANT / waveLength;
		}
		template<Arithmetic T>
		constexpr Joules<T> ParticleEnergy(Hertzes<T> frequency) noexcept {
			return Constants::PLANCK_CONSTANT * frequency;
		}
	}
}