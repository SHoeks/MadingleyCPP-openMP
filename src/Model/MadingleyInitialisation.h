#ifndef MADINGLEYINITIALISATION
#define MADINGLEYINITIALISATION

#include "Types.h"
#include "MassBinsHandler.h"
#include "FunctionalGroupDefinitions.h"
#include "UtilityFunctions.h"
#include "Grid.h"
#include "Stopwatch.h"
#include "Constants.h"
#include "Parameters.h"
#include "NonStaticSimpleRNG.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <list>
#include <vector>

/** \ brief* Initialization information for Madingley model simulations */
class MadingleyInitialisation {
public:
    MadingleyInitialisation( );

    /** \brief Reads the initialization file to get information for the set of simulations to be run
     @param initialisationFile The name of the initialization file with information on the simulations to be run
     @param outputPath The path to folder in which outputs will be stored */
    MadingleyInitialisation( long long&, double&, double&, Grid&, bool );

    /** \brief The functional group definitions of cohorts in the model */
    FunctionalGroupDefinitions mCohortFunctionalGroupDefinitions;
    /** \brief The functional group definitions of stocks in the model */
    FunctionalGroupDefinitions mStockFunctionalGroupDefinitions;



private:
    /** \brief Reads in all initialisation files */
    // FIX - Need to adjust this file to deal with incorrect inputs, extra columns etc by throwing an error.
    // FIX - Also need to strip leading spaces.
    void ReadInitialisationFiles( );
    /** \brief  Seed grid cells with cohorts, as specified in the model input files
    @param g A reference to a grid cell */
    long SeedGridCellCohorts( GridCell& );
    /** \brief    Seed grid cell with stocks, as specified in the model input files
    @param gcl The grid cell  */
    long SeedGridCellStocks( GridCell& );

    //# Initialise model using Spin-up outputs
    long SeedCohortsApplySpinUpFast( Grid&, std::vector< std::vector<std::string> > );
    long SeedStocksApplySpinUp( GridCell&, std::vector< std::vector<std::string> > );
    std::vector< std::vector<std::string> > readSpinUpStateCohort( std::string );
    std::vector< std::vector<std::string> > readSpinUpStateStock( std::string );
    bool UseMadingleySpinUp;

    std::vector< std::vector<std::string> > CohortData; // Cohort data
    std::vector< std::vector<std::string> > StockData; // Stock data
    std::vector< std::string > DataProperties; // vector with cell ids
    //#

    /** \brief An instance of the mass bin handler for the current model run */
    MassBinsHandler mModelMassBins;
    /** Instance of Utilities for timestep conversions */
    UtilityFunctions mUtilities;
    Stopwatch mInitializationTimer;
    NonStaticSimpleRNG mRandomNumber;
    /** \brief Information from the initialization file  */
    Types::StringMap mInitialisationFileStrings;
    /** \brief The paths and filenames for the diagnostics for the ecological processes */
    Types::StringMap mProcessTrackingOutputs;
    /** The full path for the output files for a set of simulations */
    std::string mOutputPath;
    /** track cohort ID number*/
    long long mNextCohortID;
};

#endif
