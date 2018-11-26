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
    
    //# Data wrting during model run
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

    //############## new create vector with relevent grid indices (terrestrial only)
    std::vector<int> TerrestrialGridcellIndices; 
    unsigned cellCounter = 0;
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gcl ) {
        if ( !gcl.IsMarine( ) ) {
            TerrestrialGridcellIndices.push_back( gcl.GetIndex( ) );
            cellCounter++;
        } 
    } );

    std::cout << "Running terrestrial only ( "<< cellCounter << " of the " << 
        Parameters::Get( )->GetNumberOfGridCells( ) << " gridcells )"<< std::endl;
    std::cout<< " " <<std::endl;
    
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
            RunWithinCellsInParallel( cellCounter, TerrestrialGridcellIndices );
        } else {
            RunWithinCells( cellCounter, TerrestrialGridcellIndices );
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
            // write pools
            std::string csvNamepools = mOutputDirectory + std::to_string(timeStep) +"_Pools.csv";
            ofstream CSVpools;
            CSVpools.open (csvNamepools);
            CSVpools << "GridcellIndex" << "," << "OrgPool" << "," << "C02Pool" << std::endl;
            for( unsigned gridCellIndex = 0; gridCellIndex < Parameters::Get( )->GetNumberOfGridCells( ); gridCellIndex++ ) {
                CSVpools << gridCellIndex << "," <<
                Environment::Get( "Organic Pool", gridCellIndex ) << "," <<
                Environment::Get( "Respiratory CO2 Pool", gridCellIndex ) << std::endl; }
            CSVpools.close();

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

    // write pools
    std::string csvNamepools = mOutputDirectory + std::to_string(99999) +"_Pools.csv";
    ofstream CSVpools;
    CSVpools.open (csvNamepools);
    CSVpools << "GridcellIndex" << "," << "OrgPool" << "," << "C02Pool" << std::endl;
    for( unsigned gridCellIndex = 0; gridCellIndex < Parameters::Get( )->GetNumberOfGridCells( ); gridCellIndex++ ) {
        CSVpools << gridCellIndex << "," <<
        Environment::Get( "Organic Pool", gridCellIndex ) << "," <<
        Environment::Get( "Respiratory CO2 Pool", gridCellIndex ) << std::endl; }
    CSVpools.close();

}

void Madingley::RunWithinCells( unsigned cellCounter, std::vector<int> TerrestrialGridcellIndices ) {
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

void Madingley::RunWithinCellsInParallel( unsigned cellCounter, std::vector<int> TerrestrialGridcellIndices ) {
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

    if ( !gcl.IsMarine( ) ) {
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

}

void Madingley::RunWithinCellCohortEcology( GridCell& gcl, ThreadVariables& partial ) {
    // Local instances of classes
    // Initialize ecology for stocks and cohorts - needed fresh every timestep?
    if ( !gcl.IsMarine( ) ) {

        EcologyCohort mEcologyCohort;
        mEcologyCohort.InitialiseEating( gcl, mParams );
        Activity CohortActivity;

        //clock_t begin = clock();
        std::vector< std::vector<int> > SortedCohortIndices;
        SortedCohortIndices = DetSortIndicesCohorts( gcl, false );
        //clock_t end = clock();
        //std::cout << double(end - begin) / CLOCKS_PER_SEC << std::endl;
        //std::cout << "##############" << std::endl;

        // std::vector< std::vector<int> > SortedCohortIndices;
        // SortedCohortIndices = mGridCell.GetSortedCohortIndices();
        // for(unsigned i = 0; i < SortedCohortIndices.size(); i++) std::cout << SortedCohortIndices[i].size() << std::endl;

        // Loop over randomly ordered gridCellCohorts to implement biological functions
        gcl.ApplyFunctionToAllCohortsWithStaticRandomness( [&]( Cohort* c ) {
            // Perform all biological functions except dispersal (which is cross grid cell)

            if( gcl.mCohorts[c->mFunctionalGroupIndex].size( ) != 0 && c->mCohortAbundance > Parameters::Get( )->GetExtinctionThreshold( ) ) {

                CohortActivity.AssignProportionTimeActive( gcl, c, mCurrentTimeStep, mCurrentMonth, mParams );

                // Run ecology
                mEcologyCohort.RunWithinCellEcology( gcl, c, mCurrentTimeStep, partial, mCurrentMonth, mParams, SortedCohortIndices);
                // Update the properties of the acting cohort
                mEcologyCohort.UpdateEcology( gcl, c, mCurrentTimeStep );
                Cohort::ResetMassFluxes( );
                // Check that the mass of individuals in this cohort is still >= 0 after running ecology
                assert( c->mIndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort" );
            }

            // Check that the mass of individuals in this cohort is still >= 0 after running ecology
            if( gcl.mCohorts[c->mFunctionalGroupIndex].size( ) > 0 ) assert( c->mIndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort" );

        }, mCurrentTimeStep );

        //for(unsigned i = 0; i < SortedCohortIndices[10].size(); i++) std::cout << SortedCohortIndices[10][i] << "  ";
        //SortedCohortIndices.clear();
        //std::cout << "##" << std::endl;

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

//##
std::vector< std::vector<int> > Madingley::DetSortIndicesCohorts( GridCell& gcl, bool print ) {

    // Vectors of mIndividualBodyMass, stored per function group
    // ("columns" are functional group index, "rows" are filled with the corrosponding mIndividualBodyMass)
    std::vector< std::vector<double> > x(20,std::vector<double>(0));
    gcl.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
        x[c->mFunctionalGroupIndex].push_back(c->mIndividualBodyMass);
    } );

    // Empty vector for storing sorted indices
    std::vector< std::vector<int> > y2(20,std::vector<int>(0));

    // Determine and store sorted indeces per functional group
    for( int ii = 0; ii < 20; ii++ ) {
      if( x[ii].size()>0 ){
        std::vector<double> x2 = x[ii];

        std::vector<int> y(x2.size());

        std::size_t n(0);
        std::generate(std::begin(y), std::end(y), [&]{ return n++; });
        std::sort( std::begin(y), std::end(y), [&](int i1, int i2) { return x2[i1] < x2[i2]; } );

        for( unsigned yy = 0; yy < y.size(); yy++) y2[ii].push_back(y[yy]);
        y.clear();
      }
    }

    x.clear();

    // Print y2 to console
    if (print == true){
      for ( unsigned yy = 0; yy < y2.size(); yy++ ) {
        for ( unsigned tt = 0; tt < y2[yy].size(); tt++ ) std::cout << y2[yy][tt] << ' ';
        std::cout << "|||" <<std::endl;
      }
    }
    return y2;
}
//##

//##
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

//##
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

//##
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

                //##### end new code


                //######  additional megafauna specific cohorts statistics (turn on/off in OutputControlParameters.csv)

                // Mega carnivores endotherm
                double MCarnivoreEndothermBMThisCell = 0;
                double MCarnivoreEndothermAThisCell = 0;
                // Mega carnivores ectotherm
                double MCarnivoreEctothermBMThisCell = 0;
                double MCarnivoreEctothermAThisCell = 0;

                // Mega Herbivores endotherm
                double MHerbivoreEndothermBMThisCell = 0;
                double MHerbivoreEndothermAThisCell = 0;
                // Mega Herbivores ectotherm
                double MHerbivoreEctothermBMThisCell = 0;
                double MHerbivoreEctothermAThisCell = 0;
                
                // Mega Omnivores endotherm
                double MOmnivoreEndothermBMThisCell = 0;
                double MOmnivoreEndothermAThisCell = 0;
                // Mega Omnivores ectotherm
                double MOmnivoreEctothermBMThisCell = 0;
                double MOmnivoreEctothermAThisCell = 0;

                ////////////////////////////////////////
                // Store largest Carnivores
                double MaxCarnivores = 0;
                // Store largest Omnivores
                double MaxOmnivores = 0;
                // Store largest Herbivores
                double MaxHerbivores = 0;
                /////////////////////////////////////

                //######  end additional megafauna specific cohorts statistics

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


                    //######  additional megafauna specific cohorts statistics (turn on/off in OutputControlParameters.csv)

                    // Mega carnivores endotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" && 
                        c->mIndividualBodyMass > 21000 ) {
                        MCarnivoreEndothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MCarnivoreEndothermAThisCell  += c->mCohortAbundance;
                    }
                    // Mega carnivores ectotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" && 
                        c->mIndividualBodyMass > 21000) {
                        MCarnivoreEctothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MCarnivoreEctothermAThisCell  += c->mCohortAbundance;
                    }

                    // Mega Herbivores endotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MHerbivoreEndothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MHerbivoreEndothermAThisCell  += c->mCohortAbundance;
                    }
                    // Mega Herbivores ectotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MHerbivoreEctothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MHerbivoreEctothermAThisCell  += c->mCohortAbundance;
                    }
                    
                    // Mega Omnivores endotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MOmnivoreEndothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MOmnivoreEndothermAThisCell  += c->mCohortAbundance;
                    }
                    // Mega Omnivores ectotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MOmnivoreEctothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MOmnivoreEctothermAThisCell  += c->mCohortAbundance;
                    }
                    
                    ////////////////////////////////////////
                    // Store largest Carnivores
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" ) {
                        if( c->mAdultMass > MaxCarnivores ) MaxCarnivores = c->mAdultMass;
                    }
                    // Store largest Omnivores
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" ) {
                        if( c->mAdultMass > MaxOmnivores ) MaxOmnivores = c->mAdultMass;
                    }
                    // Store largest Herbivores
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" ) {
                        if( c->mAdultMass > MaxHerbivores ) MaxHerbivores = c->mAdultMass;
                    }
                    /////////////////////////////////////

                    //######  end additional megafauna specific cohorts statistics
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

                //######  additional megafauna specific cohorts statistics (turn on/off in OutputControlParameters.csv)

                // Mega carnivores endotherm
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEndothermBiomassDensity", gridCell.GetIndex( ), MCarnivoreEndothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEndothermAbundance", gridCell.GetIndex( ), MCarnivoreEndothermAThisCell  / gridCell.GetIndex( ) ); 
                // Mega carnivores ectotherm
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEctothermBiomassDensity", gridCell.GetIndex( ), MCarnivoreEctothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEctothermAbundance", gridCell.GetIndex( ), MCarnivoreEctothermAThisCell  / gridCell.GetIndex( ) ); 

                // Mega Herbivores endotherm
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEndothermBiomassDensity", gridCell.GetIndex( ), MHerbivoreEndothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEndothermAbundance", gridCell.GetIndex( ), MHerbivoreEndothermAThisCell  / gridCell.GetIndex( ) ); 
                // Mega Herbivores ectotherm
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEctothermBiomassDensity", gridCell.GetIndex( ), MHerbivoreEctothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEctothermAbundance", gridCell.GetIndex( ), MHerbivoreEctothermAThisCell  / gridCell.GetIndex( ) ); 
                    
                // Mega Omnivores endotherm
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEndothermBiomassDensity", gridCell.GetIndex( ), MOmnivoreEndothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEndothermAbundance", gridCell.GetIndex( ), MOmnivoreEndothermAThisCell  / gridCell.GetIndex( ) );   
                // Mega Omnivores ectotherm
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEctothermBiomassDensity", gridCell.GetIndex( ), MOmnivoreEctothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEctothermAbundance", gridCell.GetIndex( ), MOmnivoreEctothermAThisCell  / gridCell.GetIndex( ) ); 
                
                ////////////////////////////////////////
                // Store largest Carnivores
                DataRecorder::Get( )->SetDataOn( "MaxBodyMassCarnivores", gridCell.GetIndex( ), MaxCarnivores );
                // Store largest Omnivores
                DataRecorder::Get( )->SetDataOn( "MaxBodyMassOmnivores", gridCell.GetIndex( ), MaxOmnivores );
                // Store largest Herbivores
                DataRecorder::Get( )->SetDataOn( "MaxBodyMassHerbivores", gridCell.GetIndex( ), MaxHerbivores );
                /////////////////////////////////////

                //#####################################
                DataRecorder::Get( )->SetDataOn( "NumberOfDispersalsThisCell", gridCell.GetIndex( ), OutputDispersalVector[cellCounter] );
                //#####################################

                //######  end additional megafauna specific cohorts statistics

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