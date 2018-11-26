#include <Cohort.h>
#include <limits.h>
#include <GridCell.h>
#include <MadingleyInitialisation.h>
#include <DispersalSet.h>

#include "Parameters.h"

unsigned Cohort::mNextID = 0;
Types::Double2DMap Cohort::mMassAccounting;

Cohort::Cohort( GridCell& gcl, unsigned functionalGroupIndex, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, 
                double initialAbundance, double optimalPreyBodySizeRatio, unsigned short birthTimeStep, double proportionTimeActive, 
                long long &nextCohortID, double trophicindex, double individualReproductivePotentialMass, unsigned maturityTimeStep ) {
    mFunctionalGroupIndex = functionalGroupIndex;
    mJuvenileMass = juvenileBodyMass;
    mAdultMass = adultBodyMass;
    mIndividualBodyMass = initialBodyMass;
    mCohortAbundance = initialAbundance;
    mBirthTimeStep = birthTimeStep;
    mMaturityTimeStep = maturityTimeStep; // std::numeric_limits<unsigned>::max( );
    mLogOptimalPreyBodySizeRatio = log( optimalPreyBodySizeRatio );
    mMaximumAchievedBodyMass = juvenileBodyMass;
    mMerged = false;
    mProportionTimeActive = proportionTimeActive;
    mCurrentCell = &gcl;
    mDestinationCell = mCurrentCell;
    mCurrentLocation.SetIndices( gcl.GetLatitudeIndex( ), gcl.GetLongitudeIndex( ) );
    mDestinationLocation = mCurrentLocation;

    mIndividualReproductivePotentialMass = individualReproductivePotentialMass;

    mID = mNextID; //MB added to track this object.

    mNextID++;
    nextCohortID++; // FIX - Is this increment required?

    // track cohort consumption
    //##############
    mConsumed_Autotroph = 0.0;
    mConsumed_Herbivore = 0.0;
    mConsumed_Omnivore  = 0.0;
    mConsumed_Carnivore = 0.0;
    mTrophicIndex = trophicindex;
    //##############

    // track cohort bodymass change
    //##############
    mOriginal_AdultMass = adultBodyMass;
    //##############

}

Cohort::Cohort( Cohort* actingCohort, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, double orignalAdultMass, unsigned birthTimeStep, long long& nextCohortID, double trophicindex ) {
    mFunctionalGroupIndex = actingCohort->mFunctionalGroupIndex;
    mJuvenileMass = juvenileBodyMass;
    mAdultMass = adultBodyMass;
    mIndividualBodyMass = initialBodyMass;
    mCohortAbundance = initialAbundance;
    mBirthTimeStep = birthTimeStep;
    mMaturityTimeStep = std::numeric_limits<unsigned>::max( );
    mLogOptimalPreyBodySizeRatio = actingCohort->mLogOptimalPreyBodySizeRatio;
    mMaximumAchievedBodyMass = juvenileBodyMass;
    mMerged = false;
    mProportionTimeActive = actingCohort->mProportionTimeActive;
    mCurrentCell = actingCohort->mCurrentCell;
    mDestinationCell = mCurrentCell;
    mCurrentLocation = actingCohort->mCurrentLocation;
    mDestinationLocation = mCurrentLocation;
    mIndividualReproductivePotentialMass = 0;
    mID = mNextID; //MB added to track this object.
    mNextID++;
    nextCohortID++; // FIX - Is this increment required?

    // track cohort consumption
    //##############
    mConsumed_Autotroph = 0.0;
    mConsumed_Herbivore = 0.0;
    mConsumed_Omnivore  = 0.0;
    mConsumed_Carnivore = 0.0;
    mConsumed_Total = 0.0;
    mTrophicIndex = trophicindex;
    //##############

    // track cohort bodymass change
    //##############
    mOriginal_AdultMass = orignalAdultMass;
    //##############

}

bool Cohort::IsMature( ) {
    return ( mMaturityTimeStep < std::numeric_limits<unsigned>::max( ) );
}

void Cohort::ResetMassFluxes( ) {
    // Initialize delta abundance sorted list with appropriate processes
    mMassAccounting["abundance"]["mortality"] = 0.0;

    // Initialize delta biomass sorted list with appropriate processes
    mMassAccounting["biomass"]["metabolism"] = 0.0;
    mMassAccounting["biomass"]["predation"] = 0.0;
    mMassAccounting["biomass"]["herbivory"] = 0.0;
    mMassAccounting["biomass"]["reproduction"] = 0.0;

    // Initialize delta reproductive biomass vector with appropriate processes

    mMassAccounting["reproductivebiomass"]["reproduction"] = 0.0;

    // Initialize organic pool delta vector with appropriate processes
    mMassAccounting["organicpool"]["herbivory"] = 0.0;
    mMassAccounting["organicpool"]["predation"] = 0.0;
    mMassAccounting["organicpool"]["mortality"] = 0.0;

    // Initialize respiratory CO2 pool delta vector with appropriate processes
    mMassAccounting["respiratoryCO2pool"]["metabolism"] = 0.0;
}

double Cohort::GetRealm( ) {
    return mCurrentCell->GetRealm( );
}

void Cohort::TryLivingAt( Types::GridCellPointer destination, Location& L ) {
    if( destination != 0 && destination->GetRealm( ) == GetRealm( ) ) {
        mDestinationCell = destination;
        mDestinationLocation = L;
    }
}

GridCell& Cohort::GetDestination( ) {
    return *mDestinationCell;
}

GridCell& Cohort::GetCurrentCell( ) {
    return *mCurrentCell;
}

void Cohort::SetCurrentCell( Types::GridCellPointer gclp ) {
    mCurrentCell = gclp;
}

bool Cohort::IsMoving( ) {
    return mCurrentCell != mDestinationCell;
}

void Cohort::Move( ) {
    mCurrentCell->MoveCohort( this );
}

bool Cohort::IsMarine( ) {
    return mCurrentCell->IsMarine( );
}

bool Cohort::IsPlanktonic( MadingleyInitialisation& params ) {
    return ( IsMarine( ) && ( ( mIndividualBodyMass <= Parameters::Get( )->GetPlanktonSizeThreshold( ) ) || ( params.mCohortFunctionalGroupDefinitions.GetTraitNames( "Mobility", mFunctionalGroupIndex ) == "planktonic" ) ) );
}

std::string Cohort::GetDispersalType( MadingleyInitialisation& params ) {
    std::string dispersalName;
    if( IsPlanktonic( params ) ) {
        // Advective dispersal
        dispersalName = "advective";
    }// Otherwise, if mature do responsive dispersal
    else if( IsMature( ) ) {
        dispersalName = "responsive";
    }// If the cohort is immature, run diffusive dispersal
    else {
        dispersalName = "diffusive";
    }
    return dispersalName;
}
