#include "MortalitySenescence.h"

/** \brief Constructor for senscence mortality: assigns all parameter values */
MortalitySenescence::MortalitySenescence( std::string globalModelTimeStepUnit ) {
    mTimeUnitImplementation = "Day";
    mMortalityRate = 0.003;
    // Calculate the scalar to convert from the time step units used by this implementation of mortality to the global model time step units
    mDeltaT = mUtilities.ConvertTimeUnits( globalModelTimeStepUnit, mTimeUnitImplementation );
}

double MortalitySenescence::CalculateMortalityRate( Cohort* actingCohort, double bodyMassIncludingChangeThisTimeStep, unsigned currentTimestep ) {
    
    double TimeToMaturity, AgePostMaturity;

    // Check for model restart, cohort gained maturity reached restarted from past model state (in current simulation)
    if( currentTimestep < actingCohort->mBirthTimeStep && currentTimestep > actingCohort->mMaturityTimeStep && Parameters::Get( )->GetApplyModelSpinup( ) == 1 ){
        
        // Calculate the age (in model time steps) that the cohort reached maturity
        TimeToMaturity = actingCohort->mMaturityTimeStep + ( Parameters::Get( )->GetRestartedFromTimeStep( ) -  actingCohort->mBirthTimeStep);

        // Calculate how many model time steps since the cohort reached maturity
        AgePostMaturity = currentTimestep - actingCohort->mMaturityTimeStep;
    
    // Check for model restart, cohort gained maturity before restarted from past model state
    }else if( currentTimestep < actingCohort->mBirthTimeStep && currentTimestep < actingCohort->mMaturityTimeStep && Parameters::Get( )->GetApplyModelSpinup( ) == 1 ){
    
        // Calculate the age (in model time steps) that the cohort reached maturity
        TimeToMaturity = actingCohort->mMaturityTimeStep - actingCohort->mBirthTimeStep;

        // Calculate how many model time steps since the cohort reached maturity
        AgePostMaturity = currentTimestep + ( Parameters::Get( )->GetRestartedFromTimeStep( ) - actingCohort->mMaturityTimeStep);

    }else{
        
        // Calculate the age (in model time steps) that the cohort reached maturity
        TimeToMaturity = actingCohort->mMaturityTimeStep - actingCohort->mBirthTimeStep;

        // Calculate how many model time steps since the cohort reached maturity
        AgePostMaturity = currentTimestep - actingCohort->mMaturityTimeStep;
    }
    
    // Calculate the time since maturity as a fraction of the time that it took the cohort to reach maturity
    double FractionalAgePostMaturity = AgePostMaturity / ( TimeToMaturity + 1 );

    // Calculate the mortality rate per mortality formulation time step as a function of the exponential of the previous fraction
    double AgeRelatedMortalityRate = mMortalityRate * exp( FractionalAgePostMaturity );

    // Convert the mortality rate from formulation time step units to model time step units
    return AgeRelatedMortalityRate * mDeltaT;
}
