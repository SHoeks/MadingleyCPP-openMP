#include "MadingleyInitialisation.h"
#include <random>
#include <iostream>

MadingleyInitialisation::MadingleyInitialisation( ) {
    mOutputPath = "";
    mNextCohortID = 0;
}

MadingleyInitialisation::MadingleyInitialisation( long long& nextCohortID, double& numberOfCohorts, double& numberOfStocks, Grid& modelGrid, bool UseMadingleySpinUp ) {

    //Write to console
    if(UseMadingleySpinUp==1) {
        std::cout << "\033[0;33m" << "Applying spin-up state.." << "\033[0m" << std::endl;
        // Read state model spin-up
        //std::cout << "- " << Parameters::Get( )->GetApplyModelSpinup( ) << std::endl;
        std::cout << "Importing cohort spin-up properties from: " << Parameters::Get( )->GetCohortCSVLocation( ) << std::endl;
        CohortData = readSpinUpStateCohort( Parameters::Get( )->GetCohortCSVLocation( ) ); // Cohort data
        // CohortData[0] ==> cell id
        // CohortData[1] ==> functional index
        // CohortData[2] ==> JuvenileMass
        // CohortData[3] ==> AdultMass
        // CohortData[4] ==> IndividualBodyMass
        // CohortData[5] ==> CohortAbundance
        // CohortData[6] ==> LogOptimalPreyBodySizeRatio
        // CohortData[7] ==> BirthTimeStep
        // CohortData[8] ==> ProportionTimeActive

        std::cout << "Importing stock spin-up properties from: " << Parameters::Get( )->GetStockCSVLocation( ) << std::endl;
        StockData = readSpinUpStateStock( Parameters::Get( )->GetStockCSVLocation( ) ); // Stock data
        // StockData[0] ==> cell id
        // StockData[1] ==> functional index
        // StockData[2] ==> total biomass

    } else {
        std::cout << "\033[0;33m" << "Initialising model..." <<  "\033[0m" << std::endl;
    }

    //read and store environmental layers
    Environment::Get( );
    mNextCohortID = 0;

    ReadInitialisationFiles( );
    modelGrid.SetUpGrid( );

    // Set up the cohorts and stocks
    mInitializationTimer.Start( );

    long totalCohorts = 0;
    long totalStocks = 0;

    unsigned int NumberGridCellsToInitiate = Parameters::Get( )->GetNumberOfGridCells( );
    unsigned int GridCellCounter = 1;
    std::cout << "\033[0;32m";
    
    // Initialize cohorts and stock in gridcells
    modelGrid.ApplyFunctionToAllCells( [&]( GridCell & c ) {

        // Insert cohorts and stocks
        if(UseMadingleySpinUp==1) {      
				// Initialise stocks using spin-up data
                totalStocks += SeedStocksApplySpinUp( c, StockData );
        } else {
                // Create cohorts and stocks using user input (cohort and stock defs)
                totalCohorts += SeedGridCellCohorts( c );  
                totalStocks += SeedGridCellStocks( c );
                
        }
        // Output ini progress
        std::cout << "\rProgress: " << NumberGridCellsToInitiate << " / " << GridCellCounter << flush; 
        GridCellCounter++;

    } );
	
	//# use spin-ups csv for seeding of cohorts 
    if(UseMadingleySpinUp==1) totalCohorts = SeedCohortsApplySpinUpFast( modelGrid, CohortData );
    CohortData.clear();
    StockData.clear();

    std::cout << "\033[0m" << std::endl;
    std::cout << "Total cohorts initialised: " << totalCohorts << std::endl;
    std::cout << "Total stocks created " << totalStocks << std::endl << std::endl;

    nextCohortID = mNextCohortID;
    numberOfCohorts = totalCohorts;
    numberOfStocks = totalStocks;
    mInitializationTimer.Stop( );
    Cohort::ResetMassFluxes( );
    //std::cout << "Time required: " << mInitializationTimer.GetElapsedTimeSecs( ) << std::endl;
}

void MadingleyInitialisation::ReadInitialisationFiles( ) {

    std::cout << "Reading functional group definitions..." << std::endl;
    mInitialisationFileStrings["CohortFunctional"] = Constants::cCohortDefinitionsFileName;
    mCohortFunctionalGroupDefinitions = FunctionalGroupDefinitions( Constants::cCohortDefinitionsFileName );
    mInitialisationFileStrings["StockFunctional"] = Constants::cStockDefinitionsFileName;
    mStockFunctionalGroupDefinitions = FunctionalGroupDefinitions( Constants::cStockDefinitionsFileName );
    mModelMassBins.SetUpMassBins( Constants::cMassBinDefinitionsFileName );

    //assert( CellRarefaction >= 1 && "Cell rarefaction cannot be less than 1" ); // FIX - Does this need to be uncommented?
}

long MadingleyInitialisation::SeedGridCellCohorts( GridCell& gcl ) {
    long totalCohorts = 0;
    unsigned numCohortsThisCell = 0;
    // Define local variables
    double cohortJuvenileMass;
    double cohortAdultMassRatio;
    double cohortAdultMass;
    double expectedLnAdultMassRatio;
    double totalNewBiomass = 0.0;
    double optimalPreyBodySizeRatio;

    gcl.SetCohortSize( mCohortFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex.size( ) );
    for( int FunctionalGroup: mCohortFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex ) {
        int N = mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "Initial number of GridCellCohorts", FunctionalGroup );
        if( ( mCohortFunctionalGroupDefinitions.GetTraitNames( "Realm", FunctionalGroup ) == "terrestrial" && !gcl.IsMarine( ) ) ||
                ( mCohortFunctionalGroupDefinitions.GetTraitNames( "Realm", FunctionalGroup ) == "marine" && gcl.IsMarine( ) ) ) {

            numCohortsThisCell += N;
        }
    }
    if( numCohortsThisCell > 0 );
    {
        //Loop over all functional groups in the model
        for( int functionalGroup: mCohortFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex ) {
            // If it is a functional group that corresponds to the current realm, then seed cohorts
            if( ( mCohortFunctionalGroupDefinitions.GetTraitNames( "Realm", functionalGroup ) == "terrestrial" && !gcl.IsMarine( ) ) ||
                    ( mCohortFunctionalGroupDefinitions.GetTraitNames( "Realm", functionalGroup ) == "marine" && gcl.IsMarine( ) ) ) {
                // Get the minimum and maximum possible body masses for organisms in each functional group
                double massMinimum = mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "minimum mass", functionalGroup );
                double massMaximum = mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "maximum mass", functionalGroup );

                double proportionTimeActive = mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "proportion suitable time active", functionalGroup );

                // Loop over the initial number of cohorts
                unsigned numberOfCohortsInThisFunctionalGroup = 1;

                numberOfCohortsInThisFunctionalGroup = mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "initial number of gridcellcohorts", functionalGroup );

                for( unsigned jj = 0; jj < numberOfCohortsInThisFunctionalGroup; jj++ ) {
                    mRandomNumber.SetSeed( ( unsigned int )( jj + 1 ), ( unsigned int )( ( jj + 1 ) * 3 ) );

                    // setup uniform sample
                    std::random_device rd;  //Will be used to obtain a seed for the random number engine
                    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
                    std::uniform_real_distribution<> dis(0.0, 1.0);

                    // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
                    // within the bounds of the minimum and maximum body masses for the functional group
                    cohortAdultMass = pow( 10, ( dis(gen) * ( log10( massMaximum ) - log10( 50 * massMinimum ) ) + log10( 50 * massMinimum ) ) );

                    //# update to improve e.g. predation process of large carnivores (only applies to terrestrial fauna)
                    if ( cohortAdultMass < 21000 ||
                        mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", functionalGroup ) == "herbivore" ||
                        mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", functionalGroup ) == "omnivore" ||
                        mCohortFunctionalGroupDefinitions.GetTraitNames( "Realm", functionalGroup ) == "marine"
                        ) {
                        optimalPreyBodySizeRatio = std::max( 0.01, mRandomNumber.GetNormal( 0.1, 0.02 ) ); // use regular opt pred-prey value
                    } else {
                        optimalPreyBodySizeRatio = std::max( 0.01, mRandomNumber.GetNormal( 1.0, 0.02 ) ); // use updated opt pred-prey value
                    }

                    if( !gcl.IsMarine( ) ) {
                        do {
                            expectedLnAdultMassRatio = 2.24 + 0.13 * log( cohortAdultMass );

                            cohortAdultMassRatio = 1.0 + mRandomNumber.GetLogNormal( expectedLnAdultMassRatio, 0.5 );
                            cohortJuvenileMass = cohortAdultMass * 1.0 / cohortAdultMassRatio;
                        } while( cohortAdultMass <= cohortJuvenileMass || cohortJuvenileMass < massMinimum );
                    } else {
                        do {
                            expectedLnAdultMassRatio = 2.24 + 0.13 * log( cohortAdultMass );

                            cohortAdultMassRatio = 1.0 + 10 * mRandomNumber.GetLogNormal( expectedLnAdultMassRatio, 0.5 );
                            ;
                            cohortJuvenileMass = cohortAdultMass * 1.0 / cohortAdultMassRatio;
                        } while( cohortAdultMass <= cohortJuvenileMass || cohortJuvenileMass < massMinimum );
                    }

                    double NewBiomass = ( 3300. / numCohortsThisCell ) * 100 * 3000 *
                            pow( 0.6, ( log10( cohortJuvenileMass ) ) ) * ( gcl.GetCellArea( ) );

                    totalNewBiomass += NewBiomass;
                    double NewAbund = 0.0;

                    NewAbund = NewBiomass / cohortJuvenileMass;

                    //#### ini trophic index
                    double trophicindex = 0.0;
                    if( mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", functionalGroup ) == "herbivore" ) {
                        trophicindex = 2.0;
                    }
                    if( mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", functionalGroup ) == "omnivore" ) {
                        trophicindex = 2.5;
                    }
                    if( mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", functionalGroup ) == "carnivore" ) {
                        trophicindex = 3.0;
                    }
                    //####

                    double individualReproductivePotentialMass = 0;
                    unsigned maturityTimeStep = std::numeric_limits<unsigned>::max( );

                    // Initialise the new cohort with the relevant properties
                    Cohort* NewCohort=new Cohort( gcl, functionalGroup, cohortJuvenileMass, cohortAdultMass, cohortJuvenileMass, NewAbund,
                            optimalPreyBodySizeRatio, 0, proportionTimeActive, mNextCohortID, trophicindex, individualReproductivePotentialMass, maturityTimeStep );

                    // Add the new cohort to the list of grid cell cohorts
                    gcl.mCohorts[functionalGroup].push_back( NewCohort );
                    // Increment the variable tracking the total number of cohorts in the model
                    totalCohorts++;

                }
            }
        }
    }
    return totalCohorts;
}

long MadingleyInitialisation::SeedGridCellStocks( GridCell& gcl ) {

    long totalStocks = 0;

    // Loop over all stock functional groups in the model
    for( int functionalGroupIndex: mStockFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex ) {
        // Initialise the new stock with the relevant properties
        bool success;
        Stock NewStock( mStockFunctionalGroupDefinitions, functionalGroupIndex, gcl, success );
        // Add the new stock to the list of grid cell stocks
        if( success && !gcl.IsMarine( ) ) { //# test no marine stock ini
            gcl.mStocks[functionalGroupIndex].push_back( NewStock );

            totalStocks++;
        }
    }
    return totalStocks;
}

/////////////////////////////////////////////////////////////
// Functions to seed model using Spin-up outputs           //
/////////////////////////////////////////////////////////////

long MadingleyInitialisation::SeedStocksApplySpinUp( GridCell& gcl, std::vector< std::vector<std::string> > StockData ) {

    long totalStocks = 0;
    unsigned numStocksThisCell = 0;
    int gridCellIndex = gcl.GetIndex( );

    //std::cout << gridCellIndex << std::endl;
    std::vector<int> GCindics;
    std::transform(StockData[0].begin(), StockData[0].end(), std::back_inserter(GCindics),[](const std::string& str) { return std::stoi(str); });

    int firstIndex;
    bool CellPopulated;
    if (std::binary_search (GCindics.begin(), GCindics.end(), gridCellIndex)) {
      CellPopulated = true;
      firstIndex = find(GCindics.begin(), GCindics.end(), gridCellIndex) - GCindics.begin();
    } else CellPopulated = false;

    //std::cout << "  stocks?  ==>" <<  CellPopulated << std::endl;

    // Loop over all stock functional groups in the model
    for( int functionalGroupIndex: mStockFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex ) {
        // Initialise the new stock with the relevant properties
        bool success;
        Stock NewStock( mStockFunctionalGroupDefinitions, functionalGroupIndex, gcl, success );
        // Add the new stock to the list of grid cell stocks
        if( success ) {
            gcl.mStocks[functionalGroupIndex].push_back( NewStock );
        }
    }
    //std::cout << " " << std::endl;

    // reassing mTotalBiomass using spinup stock data
    if ( CellPopulated > 0 ) {
        //std::cout << "\033[0;35m" << numStocksThisCell << " " << totalStocks << "\033[0m" << std::endl;
        int StockCounter = firstIndex;
        gcl.ApplyFunctionToAllStocks( [&]( Stock & s ) {
            //cout << "ini mass: " << s.mTotalBiomass << " spinup mass: " << std::stod(StockData[2][StockCounter]) << std::endl; 
            s.mTotalBiomass = std::stod(StockData[2][StockCounter]);
            StockCounter++;
            totalStocks++;
        } );
    }
    //std::cout << "  " << std::endl;
    return totalStocks;
}

std::vector< std::vector<std::string> > MadingleyInitialisation::readSpinUpStateCohort( std::string filedir ) {

    std::cout << filedir << std::endl;
    std::vector< std::vector<std::string> > v(12,std::vector<std::string>(0));
    std::ifstream file ( filedir );
    unsigned col = 0;
    unsigned row = 0;
    std::string value;

    while ( file.good() )
    {
        getline ( file, value, ',' );
        int pos;
        if ((pos=value.find('\n')) != std::string::npos) {
          std::string end = value; end.erase(pos);
          std::string start = value; start.erase(0,pos+1);
          // last value current row
          if(end.length() != 0) { //cout << " v = " << end << " | row = " << row << " | col = " << col << endl;
            v[col].push_back(end);
          }
          row++; col = 0; // update row and col index
          // first value next row
          if(start.length() != 0) { //cout << " v = " << start << " | row = " << row << " | col = " << col << endl;
            v[col].push_back(start);
          }
        } else { //cout << " v = " << value << " | row = " << row << " | col = " << col << endl;
          v[col].push_back(value);
        }
        col++;
    }
    return v;
}

std::vector< std::vector<std::string> > MadingleyInitialisation::readSpinUpStateStock( std::string filedir ) {

    std::cout << filedir << std::endl;
    std::vector< std::vector<std::string> > v(3,std::vector<std::string>(0));
    std::ifstream file ( filedir );
    unsigned col = 0;
    unsigned row = 0;
    std::string value;

    while ( file.good() )
    {
        getline ( file, value, ',' );
        int pos;
        if ((pos=value.find('\n')) != std::string::npos) {
          std::string end = value; end.erase(pos);
          std::string start = value; start.erase(0,pos+1);
          // last value current row
          if(end.length() != 0) { //cout << " v = " << end << " | row = " << row << " | col = " << col << endl;
            v[col].push_back(end);
          }
          row++; col = 0; // update row and col index
          // first value next row
          if(start.length() != 0) { //cout << " v = " << start << " | row = " << row << " | col = " << col << endl;
            v[col].push_back(start);
          }
        } else { //cout << " v = " << value << " | row = " << row << " | col = " << col << endl;
          v[col].push_back(value);
        }
        col++;
    }
    return v;
}

/////////////////// ini using spin-up fast ///////////////////
long MadingleyInitialisation::SeedCohortsApplySpinUpFast( Grid& modelGrid, std::vector< std::vector<std::string> > CohortData ) {

    // number of cohorts in total
    long nrows = CohortData[0].size();
    //std::cout << "n _ cohorts to ini: " << nrows << std::endl;
    long totalCohorts = 0;

    for( unsigned jj = 0; jj < nrows; jj++ ) {

        // CohortData[0] ==> gridcell index
        unsigned gridCellIndex = std::stoul(CohortData[0][jj]);

        // CohortData[1] ==> functional index
        int functionalGroup = std::stoi(CohortData[1][jj]);

        // CohortData[2] ==> JuvenileMass
        double cohortJuvenileMass = std::stod(CohortData[2][jj]);

        // CohortData[3] ==> AdultMass
        double cohortAdultMass = std::stod(CohortData[3][jj]);

        // CohortData[4] ==> IndividualBodyMass
        double IndividualBodyMass = std::stod(CohortData[4][jj]);

        // CohortData[5] ==> CohortAbundance
        double NewAbund = std::stod(CohortData[5][jj]);

        // CohortData[6] ==> LogOptimalPreyBodySizeRatio
        double optimalPreyBodySizeRatio = exp(std::stod(CohortData[6][jj]));

        // CohortData[7] ==> BirthTimeStep
        unsigned BirthTimeStep = std::stod(CohortData[7][jj]);;

        // CohortData[8] ==> ProportionTimeActive
        double ProportionTimeActive = std::stod(CohortData[8][jj]);

        // CohortData[9] ==> TrophicIndex
        double trophicindex = std::stod(CohortData[9][jj]);

        // CohortData[10] ==> individualReproductivePotentialMass
        //double individualReproductivePotentialMass = 0;
        double individualReproductivePotentialMass = std::stod(CohortData[10][jj]);
        
        // CohortData[11] ==> maturityTimeStep
        //unsigned maturityTimeStep = std::numeric_limits<unsigned>::max( );
        unsigned maturityTimeStep = std::stod(CohortData[11][jj]);

        // std::cout << 
        // gridCellIndex << ", "<< 
        // functionalGroup << ", "<< 
        // cohortJuvenileMass << ", "<< 
        // cohortAdultMass << ", "<< 
        // IndividualBodyMass << ", "<< 
        // NewAbund << ", "<< 
        // optimalPreyBodySizeRatio << ", "<< 
        // BirthTimeStep << ", "<< 
        // ProportionTimeActive << ", "<< 
        // trophicindex << ", "<< 
        // individualReproductivePotentialMass << ", "<< 
        // maturityTimeStep << std::endl;

        modelGrid.GetACell( gridCellIndex ).SetCohortSize( mCohortFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex.size( ) );

        Cohort* NewCohort=new Cohort( modelGrid.GetACell( gridCellIndex ), functionalGroup, cohortJuvenileMass, cohortAdultMass, IndividualBodyMass, NewAbund,
                optimalPreyBodySizeRatio, BirthTimeStep, ProportionTimeActive, mNextCohortID, trophicindex, individualReproductivePotentialMass, maturityTimeStep );
        
        totalCohorts++; //  total number of cohorts in the model

        //std::cout << "iterator jj: " << jj << " gridcellindex: " << gridCellIndex << " totalCohorts: " << totalCohorts << endl;
        //std::cout << " " << endl;

        modelGrid.GetACell( gridCellIndex ).mCohorts[functionalGroup].push_back( NewCohort );
    }

    ////////////////////////////////////////////////////////////////////////////////////
    std::vector<double> v;
    for( unsigned jj = 0; jj < nrows; jj++ ) {
        v.push_back(std::stoul(CohortData[0][jj]));
    }

    //# Insert cohort in empty cell (will die immediately), not possible to have empty grid cells?
    modelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {      

        double x = gridCell.GetIndex();
        if(std::find(v.begin(), v.end(), x) == v.end()) {
            //std::cout << "empty cell "<< x << std::endl;
            modelGrid.GetACell( x ).SetCohortSize( mCohortFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex.size( ) );
            Cohort* NewCohort=new Cohort( modelGrid.GetACell( x ), 10, 1, 1, 1, 1,1, 0, 0, mNextCohortID, 2, 0, std::numeric_limits<unsigned>::max( ) );
            modelGrid.GetACell( x ).mCohorts[10].push_back( NewCohort );
            totalCohorts++;
        }

    } );
    return totalCohorts;
}
