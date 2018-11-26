#include "MortalityBackground.h"

MortalityBackground::MortalityBackground( std::string globalModelTimeStepUnit ) {
    mTimeUnitImplementation = "Day";
    mMortailtyRate = 0.001;
    // Calculate the scalar to convert from the time step units used by this implementation of mortality to the global model time step units
    mDeltaT = mUtilities.ConvertTimeUnits( globalModelTimeStepUnit, mTimeUnitImplementation );
}

double MortalityBackground::CalculateMortalityRate( Cohort* actingCohort, double bodyMassIncludingChangeThisTimeStep, unsigned currentTimestep ) {
    // Convert from mortality rate per mortality formulation time step to mortality rate per model time step
    return mMortailtyRate * mDeltaT;
}
