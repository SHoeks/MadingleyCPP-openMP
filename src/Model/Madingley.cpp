#include "Madingley.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <list>

//# new
#include <dirent.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include "FileWriter.h"
#include "GridDatum.h"
#include "BasicDatum.h"
//#include "WriteModelState.h"
//# end new

//extern std::string NumberOfThreads;

Madingley::Madingley( ) {
	
    // Set up list of global diagnostics
    SetUpGlobalDiagnosticsList( );

    std::cout << "Spatial Selection: " << Parameters::Get()->GetSpatialSelection( ) << std::endl;
    std::cout << "Maximum Cohorts: " << Parameters::Get()->GetMaximumNumberOfCohorts( ) << std::endl;
    std::cout << "Minimum Longitude: " << Parameters::Get()->GetUserMinimumLongitude( ) << std::endl;
    std::cout << "Maximum Longitude: " << Parameters::Get()->GetUserMaximumLongitude( ) << std::endl;
    std::cout << "Minimum Latitude: " << Parameters::Get()->GetUserMinimumLatitude( ) << std::endl;
    std::cout << "Maximum Latitude: " << Parameters::Get()->GetUserMaximumLatitude( ) << std::endl;
    std::cout << std::endl;

    std::cout << "Outputs" << std::endl;
    std::cout << "Write Cohort Specifics: " << Parameters::Get()->GetWriteCohortSpecifics( ) << std::endl;
    std::cout << "Write Cohort Consumption: " << Parameters::Get()->GetWriteCohortConsumption( ) << std::endl;
    std::cout << "Write Consumption Summary: " << Parameters::Get()->GetWriteConsumptionSummary( ) << std::endl;
    
    //# Data writing during model run
    FileWriter fileWriter;

    //############## new output for csv files
    mOutputDirectory = fileWriter.GetOutputDirectory( );
    std::cout << "Cohort outputs: " << mOutputDirectory << std::endl;
    std::cout << std::endl;
    //############## end new

    // Initialise the cohort ID to zero
    mIniTimer.Start( );
    mNextCohortID = 0;
    bool UseMadingleySpinUp = Parameters::Get( )->GetApplyModelSpinup( ); // read yes(1)/no(0) from input csv 
    mParams = MadingleyInitialisation( mNextCohortID, mGlobalDiagnosticVariables["NumberOfCohortsInModel"],
      mGlobalDiagnosticVariables["NumberOfStocksInModel"], mModelGrid, UseMadingleySpinUp );
    mIniTimer.Stop( );
    std::cout << "Initialisation took: " << mIniTimer.GetElapsedTimeSecs( ) << std::endl;

    mDispersalSet = new DispersalSet( );

    mStockLeafStrategy = mParams.mStockFunctionalGroupDefinitions.mTraitLookupFromIndex[ "leaf strategy" ];
    mCohortNutritionSource = mParams.mCohortFunctionalGroupDefinitions.mTraitLookupFromIndex[ "nutrition source" ];
    mCohortThermoregulation = mParams.mCohortFunctionalGroupDefinitions.mTraitLookupFromIndex[ "endo/ectotherm" ];
    mCohortReproductiveStrategy = mParams.mCohortFunctionalGroupDefinitions.mTraitLookupFromIndex[ "reproductive strategy" ];

}

void Madingley::Run( ) {

    //# cohorts ini properties
    WriteModelState writemodelstate;
    writemodelstate.CohortSpinUpOutputCheck( mModelGrid, mOutputDirectory );
    writemodelstate.StockSpinUpOutputCheck( mModelGrid, mOutputDirectory );

    //# Console colors
    std::string default_console = "\033[0m"; //# default
    std::string clr1 = "\033[0;31m"; //# red
    std::string clr2 = "\033[0;35m"; //# magenta
    std::string clr3 = "\033[0;33m"; //# yellow
    std::string clr4 = "\033[0;32m"; //# green

    std::cout << clr4 << std::endl;
    std::cout << "Output model state by computation time: " << Parameters::Get( )->GetWriteStateInModelTime( ) << std::endl;
    std::cout << "Overwrite model state: " << Parameters::Get( )->GetWriteStateOverwrite( ) << std::endl;
    std::cout << "Output model state every: " << Parameters::Get( )->GetWriteStateInterval( ) << " " << Parameters::Get( )->GetWriteStateIntervalUnit( ) << std::endl;

    // Write out model run details to the console
    std::cout << clr4 << std::endl;
    std::cout << "Starting model run" << std::endl;
    unsigned RunParallel = 0;
    if(Parameters::Get( )->GetRunParallel( )==1){
        std::cout << "Running in parallel, using " << Parameters::Get( )->GetThreadNumber( ) << " threads" << std::endl;
        RunParallel = 1;
    }else{
        std::cout << "Running in serial" << std::endl;
        RunParallel = 0;
    }
    std::cout<<default_console<<std::endl;

    // Store EcoTime
    double EcoTime = 0;
    double TotalTime = 0;
    double TimeToWrite = 0;
    bool DataWritten = false;
    mDispersals = 0;
    int SimulationInMonths_print = Parameters::Get( )->GetLengthOfSimulationInMonths( );



    //# store biomass diagnostics
    std::vector< std::vector<double> > StableBiomass( SimulationInMonths_print ,std::vector<double>( 5 ));
    std::vector< std::vector<double> > StableBiomassYearly( (SimulationInMonths_print/12) ,std::vector<double>( 5 ));
    int yearcounter = 0;
    std::vector<double> WindowMeanPercChange;
    std::vector<double> YearChangePerc;

    /// Run the model
    for( unsigned timeStep = 0; timeStep < Parameters::Get( )->GetLengthOfSimulationInMonths( ); timeStep += 1 ) {

        // Set current month of the year
        TimeStep::Get( )->SetMonthly( timeStep );

        // Output timestep info
        std::cout << clr2<< "Running time step " << timeStep + 1 << " / " << SimulationInMonths_print << default_console << std::endl;

        // Get current time step and month
        mCurrentTimeStep = timeStep;
        mCurrentMonth = mUtilities.GetCurrentMonth( timeStep );
        
        // Start ecology timer
        mEcologyTimer.Start( );
        
        // Reset consumption by cohorts
        mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
            gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                c->mConsumed_Autotroph = 0.0;
                c->mConsumed_Herbivore = 0.0;
                c->mConsumed_Omnivore  = 0.0;
                c->mConsumed_Carnivore = 0.0;            
            } );
        } );
        
        // Update environment according to month
        Environment::Update( mCurrentMonth );

        // Run in gridcell ecology
        if( RunParallel == 1 ) {
            RunWithinCellsInParallel(  );
        } else {
            RunWithinCells(  );
        }

        // Calculate cohort mTrophicIndex
        if( timeStep > 0 ) CalculateTrophicIndex( );

        // Stop ecology timer
        mEcologyTimer.Stop( );
        std::cout << "Within grid ecology took: " << mEcologyTimer.GetElapsedTimeSecs( ) << std::endl;
        EcoTime += mEcologyTimer.GetElapsedTimeSecs( );
        
        // Start between gridcell ecology timer
        mDispersalTimer.Start( );

        // Output Dispersal
        std::vector<double> OutputDispersalVector;
        OutputDispersalVector = RunCrossGridCellEcology( mDispersals );

        // Stop between gridcell ecology timer
        mDispersalTimer.Stop( );
        std::cout << "Across grid ecology took: " << mDispersalTimer.GetElapsedTimeSecs( ) << std::endl;

        // Start output timer
        mOutputTimer.Start( );
        
        // Save basic and grid outputs and save biomass diagnostics
        StableBiomass[ timeStep ] = Output( timeStep, OutputDispersalVector );      
        PrintStableBiomass( StableBiomass, timeStep );

        // Calculate yearly average biomass
        if( (timeStep + 1) % 12 == 0 && timeStep != 0 ) {
            for( int v = 0; v < StableBiomass[ 0 ].size(); v++) {
                double TempYearly = 0;
                for( int jj = timeStep-11; jj < timeStep+1; jj++ ) {
                    TempYearly += StableBiomass[ jj ][ v ];
                }
                StableBiomassYearly[ yearcounter ][ v ] =  TempYearly/12;
            }
            yearcounter++;
        }

        // Write the results of dispersal to the console
        std::cout << "Total Cohorts remaining: " << mGlobalDiagnosticVariables["NumberOfCohortsInModel"] << std::endl;
        
        // Check model runtime
        TotalTime += mEcologyTimer.GetElapsedTimeSecs( ) + mOutputTimer.GetElapsedTimeSecs( ) + mDispersalTimer.GetElapsedTimeSecs( );
        TimeToWrite += mEcologyTimer.GetElapsedTimeSecs( ) + mOutputTimer.GetElapsedTimeSecs( ) + mDispersalTimer.GetElapsedTimeSecs( );
        std::cout << " " << std::endl;
        std::cout << "Current total model runtime: " << TotalTime << " seconds" << std::endl;
        std::cout << "Current total model runtime: " << TotalTime/60/60 << " hours" << std::endl;
        
        // Write model state
        if ( writemodelstate.WriteState( mModelGrid, mOutputDirectory, mCurrentMonth, timeStep, 
                                         TimeToWrite, StableBiomassYearly, StableBiomass ) == 0 ) {
            
			TimeToWrite = 0;
            OutputConsumptionCSV( timeStep );
            OutputConsumptionSummarizedCSV( timeStep );

        }

        // Stop output timer
        mOutputTimer.Stop( );
        std::cout << default_console << "Global outputs took: " << mOutputTimer.GetElapsedTimeSecs( ) << std::endl;
        std::cout << " " << std::endl;
    }

    std::cout << " " << std::endl;
    std::cout << "Total model runtime: " << TotalTime << std::endl;
    std::cout << "Total in gridcell computation time: " << EcoTime << std::endl;
    std::cout << "Mean EcoTime: " << EcoTime/SimulationInMonths_print << std::endl;
    
    writemodelstate.CohortSpinUpOutput( mModelGrid, mOutputDirectory, 99999 );
    writemodelstate.StockSpinUpOutput( mModelGrid, mOutputDirectory, 99999 );


}

void Madingley::RunWithinCells( ) {
    // Instantiate a class to hold thread locked global diagnostic variables
    ThreadVariables singleThreadDiagnostics( 0, 0, 0, mNextCohortID );

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gcl ) {

        RunWithinCellStockEcology( gcl );
        RunWithinCellCohortEcology( gcl, singleThreadDiagnostics );

    } );
    // Update the variable tracking cohort unique IDs
    mNextCohortID = singleThreadDiagnostics.mNextCohortID;

    //std::cout << mNextCohortID << std::endl;

    // Take the results from the thread local variables and apply to the global diagnostic variables
    mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] = singleThreadDiagnostics.mExtinctions - singleThreadDiagnostics.mCombinations;
    mGlobalDiagnosticVariables["NumberOfCohortsProduced"] = singleThreadDiagnostics.mProductions;
    mGlobalDiagnosticVariables["NumberOfCohortsInModel"] = mGlobalDiagnosticVariables["NumberOfCohortsInModel"] + singleThreadDiagnostics.mProductions - singleThreadDiagnostics.mExtinctions;
    mGlobalDiagnosticVariables["NumberOfCohortsCombined"] = singleThreadDiagnostics.mCombinations;
}

void Madingley::CalculateTrophicIndex( ) {
    // Calculate trophic index based on cohort consumption
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
            c->mConsumed_Total = c->mConsumed_Autotroph + c->mConsumed_Herbivore + c->mConsumed_Omnivore + c->mConsumed_Carnivore;
            if( c->mConsumed_Total > 0.0 ) {
                c->mTrophicIndex = 1 + 
                c->mConsumed_Autotroph/c->mConsumed_Total + 
                c->mConsumed_Herbivore/c->mConsumed_Total * 2 + 
                c->mConsumed_Omnivore/c->mConsumed_Total * 2.5 + 
                c->mConsumed_Carnivore/c->mConsumed_Total * 3; 
            }        
        } );
    } );
}

void Madingley::RunWithinCellsInParallel( ) {
    // Instantiate a class to hold thread locked global diagnostic variables

    #ifdef _OPENMP
    std::cout<<"Running RunWithinCellsInParallel..."<<endl;
    double startTimeTest = omp_get_wtime( );
    #endif

    list<ThreadVariables> partialsDiagnostics;
    unsigned gridCellIndexS;

    #pragma omp parallel num_threads(Parameters::Get( )->GetThreadNumber( )) shared(partialsDiagnostics)
    {
        ThreadVariables singleThreadDiagnostics( 0, 0, 0, mNextCohortID );

        #pragma omp for schedule(dynamic)
        for( unsigned gridCellIndex = 0; gridCellIndex < Parameters::Get( )->GetNumberOfGridCells( ); gridCellIndex++ )
        {   
            RunWithinCellStockEcology( mModelGrid.GetACell( gridCellIndex) );
            RunWithinCellCohortEcology( mModelGrid.GetACell( gridCellIndex ), singleThreadDiagnostics ); 
        }
        partialsDiagnostics.push_back(singleThreadDiagnostics);
    }//end parallel

    ThreadVariables globalDiagnostics( 0, 0, 0, mNextCohortID);
    for (list<ThreadVariables>::iterator it=partialsDiagnostics.begin(); it != partialsDiagnostics.end(); it++)
    {
        ThreadVariables tmp=*it;
        globalDiagnostics.mProductions+=tmp.mProductions;
        globalDiagnostics.mExtinctions+=tmp.mExtinctions;
        globalDiagnostics.mCombinations+=tmp.mCombinations;
    }

    // Update the variable tracking cohort unique IDs
    mNextCohortID = globalDiagnostics.mNextCohortID;

    // Take the results from the thread local variables and apply to the global diagnostic variables
    mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] = globalDiagnostics.mExtinctions - globalDiagnostics.mCombinations;
    mGlobalDiagnosticVariables["NumberOfCohortsProduced"] = globalDiagnostics.mProductions;
    mGlobalDiagnosticVariables["NumberOfCohortsInModel"] = mGlobalDiagnosticVariables["NumberOfCohortsInModel"] + globalDiagnostics.mProductions - globalDiagnostics.mExtinctions;
    mGlobalDiagnosticVariables["NumberOfCohortsCombined"] = globalDiagnostics.mCombinations;

    #ifdef _OPENMP
    double endTimeTest = omp_get_wtime( );
    //std::cout << "RunWithinCellsInParallel( ) took: " << endTimeTest - startTimeTest << endl;
    #endif

}

void Madingley::RunWithinCellStockEcology( GridCell& gcl ) {

    // Create a local instance of the stock ecology class
    EcologyStock MadingleyEcologyStock;
    // Get the list of functional group indices for autotroph stocks
    std::vector<int> AutotrophStockFunctionalGroups = mParams.mStockFunctionalGroupDefinitions.GetFunctionalGroupIndex( "Heterotroph/Autotroph", "Autotroph", false );
    // Loop over autotroph functional groups
    for( unsigned FunctionalGroup: AutotrophStockFunctionalGroups ) {
        for( auto& ActingStock: gcl.mStocks[FunctionalGroup] ) {

            // Run stock ecology
            MadingleyEcologyStock.RunWithinCellEcology( gcl, ActingStock, mCurrentTimeStep, mCurrentMonth, mParams );
        }
    }
 
}

void Madingley::RunWithinCellCohortEcology( GridCell& gcl, ThreadVariables& partial ) {
    // Local instances of classes
    // Initialize ecology for stocks and cohorts - needed fresh every timestep?

	EcologyCohort mEcologyCohort;
	mEcologyCohort.InitialiseEating( gcl, mParams );
	Activity CohortActivity;

	// Loop over randomly ordered gridCellCohorts to implement biological functions
	gcl.ApplyFunctionToAllCohortsWithStaticRandomness( [&]( Cohort* c ) {
		// Perform all biological functions except dispersal (which is cross grid cell)

		if( gcl.mCohorts[c->mFunctionalGroupIndex].size( ) != 0 && c->mCohortAbundance > Parameters::Get( )->GetExtinctionThreshold( ) ) {

			CohortActivity.AssignProportionTimeActive( gcl, c, mCurrentTimeStep, mCurrentMonth, mParams );

			// Run ecology
			mEcologyCohort.RunWithinCellEcology( gcl, c, mCurrentTimeStep, partial, mCurrentMonth, mParams );
			// Update the properties of the acting cohort
			mEcologyCohort.UpdateEcology( gcl, c, mCurrentTimeStep );
			Cohort::ResetMassFluxes( );
			// Check that the mass of individuals in this cohort is still >= 0 after running ecology
			assert( c->mIndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort" );
		}

		// Check that the mass of individuals in this cohort is still >= 0 after running ecology
		if( gcl.mCohorts[c->mFunctionalGroupIndex].size( ) > 0 ) assert( c->mIndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort" );

	}, mCurrentTimeStep );



	for( auto c: GridCell::mNewCohorts ) {
		gcl.InsertCohort( c );
		if( c->mDestinationCell != &gcl ) std::cout << "whut? wrong cell?" << std::endl;
	}
	partial.mProductions += GridCell::mNewCohorts.size( );
	GridCell::mNewCohorts.clear( );

	RunExtinction( gcl, partial );

	// Merge cohorts, if necessary
	if( gcl.GetNumberOfCohorts( ) > Parameters::Get( )->GetMaximumNumberOfCohorts( ) ) {
		mCohortMerger.ResetRandom( );
		partial.mCombinations += mCohortMerger.MergeToReachThresholdFast( gcl, mParams );

		//Run extinction a second time to remove those cohorts that have been set to zero abundance when merging
		RunExtinction( gcl, partial );
	}

}

void Madingley::RunExtinction( GridCell& gcl, ThreadVariables& partial ) {

    // Loop over cohorts and remove any whose abundance is below the extinction threshold
    std::vector<Cohort*> CohortsToRemove;
    gcl.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
        if( c->mCohortAbundance - Parameters::Get( )->GetExtinctionThreshold( ) <= 0 || c->mIndividualBodyMass <= 0 ) {
            CohortsToRemove.push_back( c );
            partial.mExtinctions += 1;
        }
    } );

    // Code to add the biomass to the biomass pool and dispose of the cohort
    for( auto c: CohortsToRemove ) {

        // Add biomass of the extinct cohort to the organic matter pool
        double deadMatter = ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance;
        if( deadMatter < 0 ) std::cout << "Dead " << deadMatter << std::endl;
        Environment::Get( "Organic Pool", c->GetCurrentCell( ) ) += deadMatter;
        assert( Environment::Get( "Organic Pool", c->GetCurrentCell( ) ) >= 0 && "Organic pool < 0" );

        // Remove the extinct cohort from the list of cohorts
        gcl.RemoveCohort( c );
    }
    for( auto c: CohortsToRemove ) {delete(c);}
}

std::vector<double> Madingley::RunCrossGridCellEcology( unsigned& dispersals ) {
    // Loop through each grid cell, and run dispersal for each.
    // In the original model a new dispersal object is made every timestep - this resets the random number generators
    mDispersalSet->ResetRandoms( );
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & c ) {
        mDispersalSet->RunCrossGridCellEcologicalProcess( c, mModelGrid, mParams, mCurrentMonth );
    } );

    int cellCounter = 0;
    std::vector<double> DispersalVector(Parameters::Get( )->GetNumberOfGridCells( ));
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        DispersalVector[cellCounter] = 0;
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
            if( c->IsMoving( ) ) {
                //std::cout << "cohort dispersed" << std::endl;
                DispersalVector[cellCounter]+=1;
            }
        } );
        cellCounter++;
    } );

    // Apply the changes from dispersal
    mDispersalSet->UpdateCrossGridCellEcology( dispersals );
 
    return DispersalVector;
}

void Madingley::SetUpGlobalDiagnosticsList( ) {
    // Add global diagnostic variables
    mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfCohortsProduced"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfCohortsCombined"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfCohortsInModel"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfStocksInModel"] = 0.0;
    mGlobalDiagnosticVariables["ExtAdultMass"] = 999999.0; //#
    mGlobalDiagnosticVariables["NumberOfCohortForcedToExtinction"] = 0.0; //#
}


//## Write all cohort properties to csv
void Madingley::OutputCSV( unsigned step ) {
    std::string csvName = mOutputDirectory + std::to_string(step) + "_month.csv";
    ofstream CSV;
    CSV.open (csvName);
    CSV << "c.mFunctionalGroupIndex" << "," << "c.birthtimestep" << "," << "c.mAdultMass" << "," << "c.mOriginal_AdultMass" << "," << "c.mIndividualBodyMass" << "," <<
        "c.mCohortAbundance" << "," << "c.mLogOptimalPreyBodySizeRatio" << "," << "c.TrophicIndex" << "," << "c.Diet" << "," << 
        "c.ReproductiveStrategy" << "," << "c.EndoEctotherm" << "," << "step" << "," <<  "gridCell" << "," << "lat" << "," << "long" << std::endl;

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                //std::cout << exp(c->mLogOptimalPreyBodySizeRatio) << std::endl;
                CSV << c->mFunctionalGroupIndex << "," << c->mBirthTimeStep << "," << c->mAdultMass << "," << c->mOriginal_AdultMass << "," << 
                c->mIndividualBodyMass << "," << c->mCohortAbundance << "," <<
                c->mLogOptimalPreyBodySizeRatio << "," << c->mTrophicIndex << "," << mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Diet", c->mFunctionalGroupIndex ) << "," << 
                mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "reproductive strategy", c->mFunctionalGroupIndex ) << "," << 
                mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Endo/Ectotherm", c->mFunctionalGroupIndex ) << "," << 
                step << "," <<  gridCell.GetIndex( ) << "," <<
                gridCell.GetLatitudeIndex( ) << "," << gridCell.GetLongitudeIndex( ) << std::endl;

        } );
    } );
    CSV.close();
}
//##

//## Write consupmtion to csv
void Madingley::OutputConsumptionCSV( unsigned step ) {
    std::string csvName = mOutputDirectory + std::to_string(step) + "_consumption.csv";
    ofstream CSV;
    CSV.open (csvName);
    CSV << "gridCell" << "," << "c.mFunctionalGroupIndex" << "," << 
        "Consumed_Autotroph" << "," << 
        "Consumed_Herbivore" << "," << 
        "Consumed_Omnivore" << "," <<
        "Consumed_Carnivore" << std::endl;

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
            CSV << gridCell.GetIndex( ) << "," << c->mFunctionalGroupIndex << "," <<
            c->mConsumed_Autotroph << "," << c->mConsumed_Herbivore << "," << 
            c->mConsumed_Omnivore  << "," << c->mConsumed_Carnivore << std::endl;           
        } );
    } );
    CSV.close();
}
//##

//## Write summerized consupmtion to csv
void Madingley::OutputConsumptionSummarizedCSV( unsigned step ) {
    int N_gridcell = Parameters::Get( )->GetNumberOfGridCells( );
    std::string csvName = mOutputDirectory + std::to_string(step) + "_consumptionSUM.csv";
    ofstream CSV;
    CSV.open (csvName);
    CSV << "gridCell" << "," << "c.mFunctionalGroupIndex" << "," << 
        "Consumed_Autotroph" << "," << 
        "Consumed_Herbivore" << "," << 
        "Consumed_Omnivore" << "," <<
        "Consumed_Carnivore" << std::endl;

    std::vector< std::vector<double> > v_Consumed_Autotroph(20,std::vector<double>(N_gridcell));
    std::vector< std::vector<double> > v_Consumed_Herbivore(20,std::vector<double>(N_gridcell));
    std::vector< std::vector<double> > v_Consumed_Omnivore(20,std::vector<double>(N_gridcell));
    std::vector< std::vector<double> > v_Consumed_Carnivore(20,std::vector<double>(N_gridcell));
    
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
            v_Consumed_Autotroph[ c->mFunctionalGroupIndex ][ gridCell.GetIndex( ) ] += c->mConsumed_Autotroph;
            v_Consumed_Herbivore[ c->mFunctionalGroupIndex ][ gridCell.GetIndex( ) ] += c->mConsumed_Herbivore;
            v_Consumed_Omnivore[  c->mFunctionalGroupIndex ][ gridCell.GetIndex( ) ] += c->mConsumed_Omnivore;
            v_Consumed_Carnivore[ c->mFunctionalGroupIndex ][ gridCell.GetIndex( ) ] += c->mConsumed_Carnivore;
        } ); 
    } );

    for( int gc = 0; gc < N_gridcell; gc++ ) {
        for( int ii = 0; ii < 20; ii++ ) {
            if(v_Consumed_Autotroph[ii][gc] !=0 || v_Consumed_Herbivore[ii][gc] != 0 || v_Consumed_Omnivore[ii][gc] != 0 || v_Consumed_Carnivore[ii][gc] != 0) {
                CSV << gc << "," << ii << "," <<
                v_Consumed_Autotroph[ ii ][ gc ] << "," << 
                v_Consumed_Herbivore[ ii ][ gc ] << "," << 
                v_Consumed_Omnivore[  ii ][ gc ] << "," << 
                v_Consumed_Carnivore[ ii ][ gc ] << std::endl;           
            }
    }}
    CSV.close();
}
//##

//##
void Madingley::OutputGridCSV( unsigned step ) {
    std::string csvName = mOutputDirectory + std::to_string(step) + "_grid.csv";
    ofstream CSV;
    CSV.open (csvName);
    CSV << "Index" << "," << "Lat_Index" << "," << "Long_Index" << "," << "Cell_Area" << "," <<
        "Cell_Width" << "," << "Cell_Height" << std::endl;

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
                CSV << gridCell.GetIndex( ) << "," <<
                gridCell.GetLatitudeIndex( ) << "," << 
                gridCell.GetLongitudeIndex( ) << "," << 
                gridCell.GetCellArea( ) << "," << 
                gridCell.GetCellWidth( ) << "," << 
                gridCell.GetCellHeight( ) << std::endl;

    } );
   CSV.close();
}
//##

std::vector<double> Madingley::Output( unsigned step, std::vector<double> OutputDispersalVector ) {
    double totalLivingBiomass = 0;
    double totalBiomass = 0;

    double organicMatterPool = 0;
    double respiratoryPool = 0;

    double totalStockBiomass = 0;

    double totalCohortBiomass = 0;
    long totalCohorts = 0;
    double totalCohortAbundance = 0;

    //#####################################
    double totalBiomassCarnivores = 0;
    double totalBiomassHerbivores = 0;
    double totalBiomassOmnivores = 0;
    double totalBiomassAutotrophs = 0;
    //#####################################

    int cellCounter = 0;


    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        double organicMatterThisCell = Environment::Get( "Organic Pool", gridCell ) / 1000.;
        double respirationThisCell = Environment::Get( "Respiratory CO2 Pool", gridCell ) / 1000.;

                double cohortBiomassThisCell = 0;
                double stockBiomassThisCell = 0;

                double phytoplanktonBiomassThisCell = 0;
                double evergreenBiomassThisCell = 0;
                double deciduousBiomassThisCell = 0;

                double cohortAbundanceThisCell = 0;
                double herbivoreBiomassThisCell = 0;
                double herbivoreAbundanceThisCell = 0;
                double omnivoreBiomassThisCell = 0;
                double omnivoreAbundanceThisCell = 0;
                double carnivoreBiomassThisCell = 0;
                double carnivoreAbundanceThisCell = 0;

                double ectothermBiomassThisCell = 0;
                double ectothermAbundanceThisCell = 0;
                double endothermBiomassThisCell = 0;
                double endothermAbundanceThisCell = 0;

                double iteroparousBiomassThisCell = 0;
                double iteroparousAbundanceThisCell = 0;
                double semelparousBiomassThisCell = 0;
                double semelparousAbundanceThisCell = 0;

                //###### custom madingley outputs (turn on/off in OutputControlParameters.csv)

                // split carnivores (endo/ecto)
                double carnivoreEndothermBiomassThisCell = 0;
                double carnivoreEndothermAbundanceThisCell = 0;
                double carnivoreEctothermBiomassThisCell = 0;
                double carnivoreEctothermAbundanceThisCell = 0;

                // split omnivores (endo/ecto)
                double omnivoreEndothermBiomassThisCell = 0;
                double omnivoreEndothermAbundanceThisCell = 0;
                double omnivoreEctothermBiomassThisCell = 0;
                double omnivoreEctothermAbundanceThisCell = 0;

                // split herbivores (endo/ecto)
                double herbivoreEndothermBiomassThisCell = 0;
                double herbivoreEndothermAbundanceThisCell = 0;
                double herbivoreEctothermBiomassThisCell = 0;
                double herbivoreEctothermAbundanceThisCell = 0;

                organicMatterPool += organicMatterThisCell;
                respiratoryPool += respirationThisCell;

                gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                    totalCohorts += 1;
                    totalCohortAbundance += c->mCohortAbundance;

                    double cohortBiomass = ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                    cohortBiomassThisCell += cohortBiomass;
                    totalCohortBiomass += cohortBiomass;
                    cohortAbundanceThisCell += c->mCohortAbundance;
                    

                    
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" ) totalBiomassCarnivores += cohortBiomass;
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" ) totalBiomassHerbivores += cohortBiomass;
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" ) totalBiomassOmnivores += cohortBiomass;

                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" ) {
                        herbivoreBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        herbivoreAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" ) {
                        omnivoreBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        omnivoreAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" ) {
                        carnivoreBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        carnivoreAbundanceThisCell += c->mCohortAbundance;
                    }

                    if( mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        ectothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        ectothermAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        endothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        endothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    if( mCohortReproductiveStrategy[ c->mFunctionalGroupIndex ] == "iteroparity" ) {
                        iteroparousBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        iteroparousAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortReproductiveStrategy[ c->mFunctionalGroupIndex ] == "semelparity" ) {
                        semelparousBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        semelparousAbundanceThisCell += c->mCohortAbundance;
                    }

                    //###### custom madingley outputs (turn on/off in OutputControlParameters.csv)
                    // carni + endo
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        carnivoreEndothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        carnivoreEndothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // carni + ecto
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        carnivoreEctothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        carnivoreEctothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // omni + endo
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        omnivoreEndothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        omnivoreEndothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // omni + ecto
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        omnivoreEctothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        omnivoreEctothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // herbi + endo
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        herbivoreEndothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        herbivoreEndothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // herbi + ecto
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        herbivoreEctothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        herbivoreEctothermAbundanceThisCell += c->mCohortAbundance;
                    }
                    //###### end new code



                } );
                
        gridCell.ApplyFunctionToAllStocks( [&]( Stock & s ) {
            double thisStockBiomass = s.mTotalBiomass / 1000.; //convert from g to kg
            stockBiomassThisCell += thisStockBiomass;
            totalStockBiomass += thisStockBiomass;

            if( mStockLeafStrategy[ s.mFunctionalGroupIndex ] == "na" ) phytoplanktonBiomassThisCell += thisStockBiomass;
            else if( mStockLeafStrategy[ s.mFunctionalGroupIndex ] == "deciduous" ) deciduousBiomassThisCell += thisStockBiomass;
            else if( mStockLeafStrategy[ s.mFunctionalGroupIndex ] == "evergreen" ) evergreenBiomassThisCell += thisStockBiomass;
            } );

        double biomassThisCell = cohortBiomassThisCell + stockBiomassThisCell + respirationThisCell + organicMatterThisCell;

                //# outputs included in MonthlyGridOutputs.nc
                DataRecorder::Get( )->SetDataOn( "BiomassDensity", gridCell.GetIndex( ), biomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "AbundanceDensity", gridCell.GetIndex( ), cohortAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "AutotrophBiomassDensity", gridCell.GetIndex( ), stockBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HeterotrophBiomassDensity", gridCell.GetIndex( ), cohortBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "PhytoplanktonBiomassDensity", gridCell.GetIndex( ), phytoplanktonBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "DeciduousBiomassDensity", gridCell.GetIndex( ), deciduousBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EvergreenBiomassDensity", gridCell.GetIndex( ), evergreenBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreBiomassDensity", gridCell.GetIndex( ), herbivoreBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreAbundanceDensity", gridCell.GetIndex( ), herbivoreAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreBiomassDensity", gridCell.GetIndex( ), omnivoreBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreAbundanceDensity", gridCell.GetIndex( ), omnivoreAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreBiomassDensity", gridCell.GetIndex( ), carnivoreBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreAbundanceDensity", gridCell.GetIndex( ), carnivoreAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EctothermBiomassDensity", gridCell.GetIndex( ), ectothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EctothermAbundanceDensity", gridCell.GetIndex( ), ectothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EndothermBiomassDensity", gridCell.GetIndex( ), endothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EndothermAbundanceDensity", gridCell.GetIndex( ), endothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "IteroparityBiomassDensity", gridCell.GetIndex( ), iteroparousBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "IteroparityAbundanceDensity", gridCell.GetIndex( ), iteroparousAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "SemelparityBiomassDensity", gridCell.GetIndex( ), semelparousBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "SemelparityAbundanceDensity", gridCell.GetIndex( ), semelparousAbundanceThisCell / gridCell.GetCellArea( ) );


                //###### new madingley grid outputs (turn on/off in OutputControlParameters.csv)
                DataRecorder::Get( )->SetDataOn( "CarnivoreEndothermBiomass", gridCell.GetIndex( ), carnivoreEndothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreEndothermAbundance", gridCell.GetIndex( ), carnivoreEndothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreEctothermBiomass", gridCell.GetIndex( ), carnivoreEctothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreEctothermAbundance", gridCell.GetIndex( ), carnivoreEctothermAbundanceThisCell / gridCell.GetCellArea( ) );

                DataRecorder::Get( )->SetDataOn( "OmnivoreEndothermBiomass", gridCell.GetIndex( ), omnivoreEndothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreEndothermAbundance", gridCell.GetIndex( ), omnivoreEndothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreEctothermBiomass", gridCell.GetIndex( ), omnivoreEctothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreEctothermAbundance", gridCell.GetIndex( ), omnivoreEctothermAbundanceThisCell / gridCell.GetCellArea( ) );

                DataRecorder::Get( )->SetDataOn( "HerbivoreEndothermBiomass", gridCell.GetIndex( ), herbivoreEndothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreEndothermAbundance", gridCell.GetIndex( ), herbivoreEndothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreEctothermBiomass", gridCell.GetIndex( ), herbivoreEctothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreEctothermAbundance", gridCell.GetIndex( ), herbivoreEctothermAbundanceThisCell / gridCell.GetCellArea( ) );

                //# output gridcell and environment properties
                DataRecorder::Get( )->SetDataOn( "GridCellIndex", gridCell.GetIndex( ), gridCell.GetIndex( ) ); // gridcell index
                DataRecorder::Get( )->SetDataOn( "GridCellArea", gridCell.GetIndex( ), gridCell.GetCellArea( ) ); // gridcell area km2
                DataRecorder::Get( )->SetDataOn( "GridCellWidth", gridCell.GetIndex( ), gridCell.GetCellWidth( ) ); // gridcell width km
                DataRecorder::Get( )->SetDataOn( "GridCellHeight", gridCell.GetIndex( ), gridCell.GetCellHeight( ) ); // gridcell height km

                DataRecorder::Get( )->SetDataOn( "GridCellTemperature", gridCell.GetIndex( ), Environment::Get( "Temperature", gridCell.GetIndex( ) ) ); // gridcell temperature
                DataRecorder::Get( )->SetDataOn( "GridCelluVel", gridCell.GetIndex( ), Environment::Get( "uVel", gridCell.GetIndex( ) ) ); // gridcell uVel
                DataRecorder::Get( )->SetDataOn( "GridCellvVel", gridCell.GetIndex( ), Environment::Get( "vVel", gridCell.GetIndex( ) ) ); // gridcell vVel
                DataRecorder::Get( )->SetDataOn( "GridCellDiurnalTemperatureRange", gridCell.GetIndex( ), Environment::Get( "DiurnalTemperatureRange", gridCell.GetIndex( ) ) ); // DiurnalTemperatureRange
                DataRecorder::Get( )->SetDataOn( "GridCellTotalPrecip", gridCell.GetIndex( ), Environment::Get( "TotalPrecip", gridCell.GetIndex( ) ) ); // gridcell total precipitation
                DataRecorder::Get( )->SetDataOn( "GridCellPrecipitation", gridCell.GetIndex( ), Environment::Get( "Precipitation", gridCell.GetIndex( ) ) ); // gridcell precipitation
                DataRecorder::Get( )->SetDataOn( "GridCellNPP", gridCell.GetIndex( ), Environment::Get( "NPP", gridCell.GetIndex( ) ) ); // gridcell NPP
                //###### end new code

                //#####################################
                DataRecorder::Get( )->SetDataOn( "NumberOfDispersalsThisCell", gridCell.GetIndex( ), OutputDispersalVector[cellCounter] );
                //#####################################


            cellCounter++;
    } );

    //# outputs included in MonthlyBasicOutputs.nc
    totalLivingBiomass = totalCohortBiomass + totalStockBiomass;
    totalBiomass = totalCohortBiomass + totalStockBiomass + respiratoryPool + organicMatterPool;
    totalBiomassAutotrophs = totalStockBiomass;

    std::vector<double> StableBiomass(5, 0);
    StableBiomass[0] = totalBiomassCarnivores; 
    StableBiomass[1] = totalBiomassHerbivores;
    StableBiomass[2] = totalBiomassOmnivores;
    StableBiomass[3] = totalCohortBiomass;
    StableBiomass[4] = totalBiomassAutotrophs;

    DataRecorder::Get( )->SetDataOn( "InCellTime", mEcologyTimer.GetElapsedTimeSecs( ) );
    DataRecorder::Get( )->SetDataOn( "DispersalTime", mDispersalTimer.GetElapsedTimeSecs( ) );

    DataRecorder::Get( )->SetDataOn( "TotalBiomass", totalBiomass );
    DataRecorder::Get( )->SetDataOn( "TotalLivingBiomass", totalLivingBiomass );
    DataRecorder::Get( )->SetDataOn( "TotalStockBiomass", totalStockBiomass );
    DataRecorder::Get( )->SetDataOn( "TotalCohortBiomass", totalCohortBiomass );
    DataRecorder::Get( )->SetDataOn( "OrganicMatterPool", organicMatterPool );
    DataRecorder::Get( )->SetDataOn( "RespiratoryCO2Pool", respiratoryPool );

    DataRecorder::Get( )->SetDataOn( "NumberOfStocks", mGlobalDiagnosticVariables["NumberOfStocksInModel"] );
    DataRecorder::Get( )->SetDataOn( "NumberOfCohorts", mGlobalDiagnosticVariables["NumberOfCohortsInModel"] );

    DataRecorder::Get( )->SetDataOn( "CohortsProduced", mGlobalDiagnosticVariables["NumberOfCohortsProduced"] );
    DataRecorder::Get( )->SetDataOn( "CohortsExtinct", mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] );
    DataRecorder::Get( )->SetDataOn( "CohortsCombined", mGlobalDiagnosticVariables["NumberOfCohortsCombined"] );
    DataRecorder::Get( )->SetDataOn( "CohortsDispersed", mDispersals );
    DataRecorder::Get( )->SetDataOn( "CohortAbundance", totalCohortAbundance );

    return StableBiomass;
}

void Madingley::PrintStableBiomass( std::vector< std::vector<double> >  StableBiomass, unsigned timeStep ) {
    std::cout << "[ TotalCarnivoreBiomass ] " << StableBiomass[ timeStep ][0] << std::endl;
    std::cout << "[ TotalHerbivoreBiomass ] " << StableBiomass[ timeStep ][1] << std::endl;
    std::cout << "[ TotalOmnivoreBiomass ] "  << StableBiomass[ timeStep ][2] << std::endl;
    std::cout << "[ TotalCohortBiomass ] "    << StableBiomass[ timeStep ][3] << std::endl;
    std::cout << "[ TotalAutotrophBiomass ] " << StableBiomass[ timeStep ][4] << std::endl;
}
