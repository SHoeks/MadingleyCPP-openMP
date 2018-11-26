#include "MetabolismEndotherm.h"
#include "Environment.h"

MetabolismEndotherm::MetabolismEndotherm( std::string globalModelTimeStepUnit ) {
    mTimeUnitImplementation = "day";
    // Parameters from fitting to Nagy 1999 Field Metabolic Rates for mammals and birds, and assuming that these endotherms are metabolising with a body temperature of 310K (37C)
    mMetabolismMassExponent = 0.7;
    mNormalizationConstant = 9.0809083973E+11;
    mActivationEnergy = 0.69; // includes endotherms in hibernation and torpor
    mBoltzmannConstant = 8.617e-5;
    // Currently a very rough estimate based on calorific values of fat, protein and carbohydrate - assumes organism is metabolising mass of 1/4 protein, 1/4 carbohydrate and 1/2 fat 
    mEnergyScalar = 1 / 27.25;
    // Set the constant to convert temperature in degrees Celsius to Kelvin
    mTemperatureUnitsConvert = 273.0;
    // Assume all endotherms have a constant body temperature of 37degC
    mEndothermBodyTemperature = 37.0 + mTemperatureUnitsConvert;
    // Calculate the scalar to convert from the time step units used by this implementation of metabolism to the global  model time step units
    mDeltaT = mUtilities.ConvertTimeUnits( globalModelTimeStepUnit, mTimeUnitImplementation );
}

void MetabolismEndotherm::Run( Cohort* actingCohort, unsigned currentTimestep, unsigned currentMonth ) {
    // Calculate metabolic loss for an individual and add the value to the delta biomass for metabolism
    Cohort::mMassAccounting[ "biomass" ][ "metabolism" ] = -CalculateIndividualMetabolicRate( actingCohort->mIndividualBodyMass, Environment::Get( "Temperature", actingCohort->GetCurrentCell( ) ) + mTemperatureUnitsConvert ) * mDeltaT;

    // If metabolic loss is greater than individual body mass after herbivory and predation, then set equal to individual body mass
    Cohort::mMassAccounting[ "biomass" ][ "metabolism" ] = std::max( Cohort::mMassAccounting[ "biomass" ][ "metabolism" ], -( actingCohort->mIndividualBodyMass + Cohort::mMassAccounting[ "biomass" ][ "predation" ] + Cohort::mMassAccounting[ "biomass" ][ "herbivory" ] ) );

    // Add total metabolic loss for all individuals in the cohort to delta biomass for metabolism in the respiratory CO2 pool
    Cohort::mMassAccounting[ "respiratoryCO2pool" ][ "metabolism" ] = -Cohort::mMassAccounting[ "biomass" ][ "metabolism" ] * actingCohort->mCohortAbundance;

}

double MetabolismEndotherm::CalculateIndividualMetabolicRate( double individualBodyMass, double temperature ) {
    // Calculate metabolic loss in kJ
    double metabolicLosskJ = mNormalizationConstant * pow( individualBodyMass, mMetabolismMassExponent ) * exp( -( mActivationEnergy / ( mBoltzmannConstant * mEndothermBodyTemperature ) ) );

    // Return metabolic loss in grams
    return metabolicLosskJ * mEnergyScalar;
}
