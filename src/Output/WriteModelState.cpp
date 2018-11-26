#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <dirent.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include "GridDatum.h"
#include "BasicDatum.h"
#include "Stopwatch.h"
#include "ThreadVariables.h"
#include "Environment.h"
#include "TimeStep.h"
#include "Maths.h"
#include "Parameters.h"
#include "DataRecorder.h"
#include "Types.h"
#include "WriteModelState.h"
#include <fstream>
#include "Grid.h"
#include "GridCell.h"

double WriteModelState::WriteState( Grid& madingleyGrid, std::string output_directory, unsigned current_month, unsigned step, double time_to_write, 
                              std::vector< std::vector<double> > StableBiomassYearly, 
                              std::vector< std::vector<double> > StableBiomass ) {

    //std::cout << "T " << time_to_write << std::endl;

    if ( Parameters::Get( )->GetWriteStateIntervalUnit( ) == "s") {
        std::cout << "Time to next data write: " << 
        Parameters::Get( )->GetWriteStateInterval( ) - time_to_write << " seconds, time counter: " << 
        time_to_write << " seconds" << std::endl;
        
        if( time_to_write > Parameters::Get( )->GetWriteStateInterval( ) && current_month == 11 ){
            std::cout << "\033[0;32m" << std::endl;
            std::cout << "Writing model state..." << std::endl;

            if (Parameters::Get( )->GetWriteStateOverwrite( ) == 0) {
                CohortSpinUpOutput( madingleyGrid, output_directory, step + 1 );
                StockSpinUpOutput( madingleyGrid, output_directory, step + 1 );
            } else {
                CohortSpinUpOutput( madingleyGrid, output_directory, 0 );
                StockSpinUpOutput( madingleyGrid, output_directory, 0 );
            }
            //WriteGridOutputs( output_directory ); WriteBasicOutputs( output_directory );
            WriteTimeLineBiomass( output_directory, step + 1, StableBiomassYearly, StableBiomass);
            time_to_write = 0;
        }
            
    } 
    else if( Parameters::Get( )->GetWriteStateIntervalUnit( ) == "m") {
        std::cout << "Time to next data write: " << 
        Parameters::Get( )->GetWriteStateInterval( ) - time_to_write/60 << " minutes, time counter: " << 
        time_to_write << " seconds" << std::endl;

        if( time_to_write/60 > Parameters::Get( )->GetWriteStateInterval( ) && current_month == 11 ){
            std::cout << "\033[0;32m" << std::endl;
            std::cout << "Writing model state..." << std::endl;
            
            if (Parameters::Get( )->GetWriteStateOverwrite( ) == 0) {
                CohortSpinUpOutput( madingleyGrid, output_directory, step + 1 );
                StockSpinUpOutput( madingleyGrid, output_directory, step + 1 );
            } else {
                CohortSpinUpOutput( madingleyGrid, output_directory, 0 );
                StockSpinUpOutput( madingleyGrid, output_directory, 0 );
            }
            //WriteGridOutputs( output_directory ); WriteBasicOutputs( output_directory );
            WriteTimeLineBiomass( output_directory, step + 1, StableBiomassYearly, StableBiomass);
            time_to_write = 0;
        }
            
    } 
    else if( Parameters::Get( )->GetWriteStateIntervalUnit( ) == "h") {
        std::cout << "Time to next data write: " << 
        Parameters::Get( )->GetWriteStateInterval( ) - time_to_write/60/60 << " hours, time counter: " << 
        time_to_write << " seconds" << std::endl;
        
        if( time_to_write/60/60 > Parameters::Get( )->GetWriteStateInterval( ) && current_month == 11 ){
            std::cout << "\033[0;32m" << std::endl;
            std::cout << "Writing model state..." << std::endl;

            if (Parameters::Get( )->GetWriteStateOverwrite( ) == 0) {
                CohortSpinUpOutput( madingleyGrid, output_directory, step + 1 );
                StockSpinUpOutput( madingleyGrid, output_directory, step + 1 );
            } else {
                CohortSpinUpOutput( madingleyGrid, output_directory, 0 );
                StockSpinUpOutput( madingleyGrid, output_directory, 0 );
            }
            //WriteGridOutputs( output_directory ); WriteBasicOutputs( output_directory );
            WriteTimeLineBiomass( output_directory, step + 1, StableBiomassYearly, StableBiomass);
            time_to_write = 0;
        }
            
    }

    //std::cout << "T " << time_to_write << std::endl;

    return time_to_write;
}

void WriteModelState::WriteTimeLineBiomass( std::string output_directory, unsigned SimulationInMonths_print, 
                                        std::vector< std::vector<double> > StableBiomassYearly, 
                                        std::vector< std::vector<double> > StableBiomass ) {
    
    // write yearly biomass
    std::string csvName1 = output_directory +  "Timeline_YearlyBiomass.csv";
    ofstream CSV1;
    CSV1.open (csvName1);
    CSV1 << "Year" << "," << "TotalCarnivoreBiomass" << "," << "TotalHerbivoreBiomass" << "," << 
        "TotalOmnivoreBiomass" << "," << "TotalCohortBiomass" << "," << "TotalAutotrophBiomass" << std::endl;

    for(int ii = 0; ii < (SimulationInMonths_print/12); ii++){
        CSV1 << ii+1 << ", " <<
        StableBiomassYearly[ii][0] << ", " <<
        StableBiomassYearly[ii][1] << ", " <<
        StableBiomassYearly[ii][2] << ", " <<
        StableBiomassYearly[ii][3] << ", " <<
        StableBiomassYearly[ii][4] << std::endl;
    }
    CSV1.close();

    // write monthly biomass
    std::string csvName2 = output_directory +  "Timeline_MontlyBiomass.csv";
    ofstream CSV2;
    CSV2.open (csvName2);
    CSV2 << "Month" << "," << "Year" << "," << "TotalCarnivoreBiomass" << "," << "TotalHerbivoreBiomass" << "," << 
        "TotalOmnivoreBiomass" << "," << "TotalCohortBiomass" << "," << "TotalAutotrophBiomass" << std::endl;

    int yearii = 1;
    for(int ii = 0; ii < (SimulationInMonths_print); ii++){
        CSV2 << ii+1 << ", " << yearii << ", " <<
        StableBiomass[ii][0] << ", " <<
        StableBiomass[ii][1] << ", " <<
        StableBiomass[ii][2] << ", " <<
        StableBiomass[ii][3] << ", " <<
        StableBiomass[ii][4] << std::endl;
        if((ii+1) % 12 == 0) yearii++;
    }
    CSV2.close();

}

//##
void WriteModelState::CohortSpinUpOutput( Grid& madingleyGrid, std::string output_directory, unsigned step ) {
    std::string csvName = "C.csv";
    if(step == 0) {
        csvName = output_directory + "CohortState.csv";
    } else {
        csvName = output_directory + "CohortState_" + std::to_string(step) + ".csv";
    }
    ofstream CSV;
    CSV.open (csvName);
    madingleyGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                CSV <<
                gridCell.GetIndex( ) << "," <<
                c->mFunctionalGroupIndex << "," <<
                c->mJuvenileMass << "," <<
                c->mAdultMass << "," <<
                c->mIndividualBodyMass << "," <<
                c->mCohortAbundance << "," <<
                c->mLogOptimalPreyBodySizeRatio << "," <<
                c->mBirthTimeStep << "," <<
                c->mProportionTimeActive << "," <<
                c->mTrophicIndex << "," <<
                c->mIndividualReproductivePotentialMass << "," <<
                c->mMaturityTimeStep <<
                std::endl;
        } );
    } );
   CSV.close();
}
//##

//##
void WriteModelState::CohortSpinUpOutputCheck( Grid& madingleyGrid, std::string output_directory ) {
    std::string csvName = output_directory + "CohortState_AfterIni" + ".csv";
    ofstream CSV;
    CSV.open (csvName);
    madingleyGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                CSV <<
                gridCell.GetIndex( ) << "," <<
                c->mFunctionalGroupIndex << "," <<
                c->mJuvenileMass << "," <<
                c->mAdultMass << "," <<
                c->mIndividualBodyMass << "," <<
                c->mCohortAbundance << "," <<
                c->mLogOptimalPreyBodySizeRatio << "," <<
                c->mBirthTimeStep << "," <<
                c->mProportionTimeActive << "," <<
                c->mTrophicIndex << "," <<
                c->mIndividualReproductivePotentialMass << "," <<
                c->mMaturityTimeStep <<
                std::endl;
        } );
    } );
   CSV.close();
}
//##

//##
void WriteModelState::StockSpinUpOutput( Grid& madingleyGrid, std::string output_directory, unsigned step ) {
    std::string csvName2 = "S.csv";
    if(step == 0) {
        csvName2 = output_directory + "StockState.csv";
    } else {
        csvName2 = output_directory + "StockState_" + std::to_string(step) + ".csv";
    }
    ofstream CSV2;
    CSV2.open (csvName2);
    madingleyGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllStocks( [&]( Stock & s ) {
              CSV2 <<
              gridCell.GetIndex( ) << "," <<
              s.mFunctionalGroupIndex << "," <<
              s.mTotalBiomass <<
              std::endl;
        } );
    } );
}
//##

//##
void WriteModelState::StockSpinUpOutputCheck( Grid& madingleyGrid, std::string output_directory ) {
    std::string csvName2 = output_directory + "StockState_AfterIni" + ".csv";
    ofstream CSV2;
    CSV2.open (csvName2);
    madingleyGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllStocks( [&]( Stock & s ) {
              CSV2 <<
              gridCell.GetIndex( ) << "," <<
              s.mFunctionalGroupIndex << "," <<
              s.mTotalBiomass <<
              std::endl;
        } );
    } );
}
//##

//##
void WriteModelState::WriteGridOutputs( std::string output_directory ) {

    Types::GridDatumMap gridDatumMap = DataRecorder::Get( )->GetGridDatumMap( );

    // Separate monthly and annual grid datums
    if( gridDatumMap.size( ) > 0 ) {
        Types::GridDatumVector annualGridDatumVector;
        Types::GridDatumVector monthlyGridDatumVector;

        //int counter = 0;
        for( Types::GridDatumMap::iterator iter = gridDatumMap.begin( ); iter != gridDatumMap.end( ); ++iter ) {
            //std::cout << counter << " ";
            if( iter->second->GetTimeUnit( ) == Constants::cAnnualTimeUnitName ) {
                annualGridDatumVector.push_back( iter->second );
            } else if( iter->second->GetTimeUnit( ) == Constants::cMonthlyTimeUnitName ) {
                monthlyGridDatumVector.push_back( iter->second );
            }
            //counter++;
        }

        // Write annual grid datums
        if( annualGridDatumVector.size( ) > 0 ) {
            std::string filePath = output_directory + Constants::cAnnualGridOutputsFileName_State;

            try {
                netCDF::NcFile annualGridOutputsNcFile( filePath, netCDF::NcFile::replace ); // Creates file

                netCDF::NcDim longitudeDim = annualGridOutputsNcFile.addDim( Constants::cLongitudeVariableNames[ 0 ], Parameters::Get( )->GetLengthUserLongitudeArray( ) );
                netCDF::NcVar longitudeNcVar = annualGridOutputsNcFile.addVar( Constants::cLongitudeVariableNames[ 0 ], netCDF::ncFloat, longitudeDim );
                longitudeNcVar.putVar( Parameters::Get( )->GetUserLongitudeArray( ) );
                longitudeNcVar.putAtt( Constants::cUnitsString, Constants::cLongitudeVariableUnit );

                netCDF::NcDim latitudeDim = annualGridOutputsNcFile.addDim( Constants::cLatitudeVariableNames[ 0 ], Parameters::Get( )->GetLengthUserLatitudeArray( ) );
                netCDF::NcVar latitudeNcVar = annualGridOutputsNcFile.addVar( Constants::cLatitudeVariableNames[ 0 ], netCDF::ncFloat, latitudeDim );
                latitudeNcVar.putVar( Parameters::Get( )->GetUserLatitudeArray( ) );
                latitudeNcVar.putAtt( Constants::cUnitsString, Constants::cLatitudeVariableUnit );

                netCDF::NcDim annualTimeNcDim = annualGridOutputsNcFile.addDim( Constants::cTimeVariableNames[ 0 ], Parameters::Get( )->GetLengthOfSimulationInYears( ) );
                netCDF::NcVar annualTimeNcVar = annualGridOutputsNcFile.addVar( Constants::cTimeVariableNames[ 0 ], netCDF::ncFloat, annualTimeNcDim );
                annualTimeNcVar.putVar( Parameters::Get( )->GetAnnualTimeStepArray( ) );
                annualTimeNcVar.putAtt( Constants::cUnitsString, Constants::cAnnualTimeUnitName );

                Types::NcDimVector dataDimensions;
                dataDimensions.push_back( annualTimeNcDim );
                dataDimensions.push_back( latitudeDim );
                dataDimensions.push_back( longitudeDim );

                for( Types::GridDatumVector::iterator iter = annualGridDatumVector.begin( ); iter != annualGridDatumVector.end( ); ++iter ) {
                    netCDF::NcVar annualGridDatumNcVar = annualGridOutputsNcFile.addVar( ( *iter )->GetName( ), netCDF::ncFloat, dataDimensions );
                    annualGridDatumNcVar.putAtt( Constants::cUnitsString, ( *iter )->GetDataUnit( ) );
                    annualGridDatumNcVar.putVar( ( *iter )->GetData( ) );
                }

            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }

        //std::cout << monthlyGridDatumVector.size( ) << std::endl;
        // Write monthly grid datums
        if( monthlyGridDatumVector.size( ) > 0 ) {
            std::string filePath = output_directory + Constants::cMonthlyGridOutputsFileName_State;

            try {
                netCDF::NcFile monthlyGridOutputsNcFile( filePath, netCDF::NcFile::replace ); // Creates file

                netCDF::NcDim longitudeDim = monthlyGridOutputsNcFile.addDim( Constants::cLongitudeVariableNames[ 0 ], Parameters::Get( )->GetLengthUserLongitudeArray( ) );
                netCDF::NcVar longitudeNcVar = monthlyGridOutputsNcFile.addVar( Constants::cLongitudeVariableNames[ 0 ], netCDF::ncFloat, longitudeDim );
                longitudeNcVar.putVar( Parameters::Get( )->GetUserLongitudeArray( ) );
                longitudeNcVar.putAtt( Constants::cUnitsString, Constants::cLongitudeVariableUnit );

                netCDF::NcDim latitudeDim = monthlyGridOutputsNcFile.addDim( Constants::cLatitudeVariableNames[ 0 ], Parameters::Get( )->GetLengthUserLatitudeArray( ) );
                netCDF::NcVar latitudeNcVar = monthlyGridOutputsNcFile.addVar( Constants::cLatitudeVariableNames[ 0 ], netCDF::ncFloat, latitudeDim );
                latitudeNcVar.putVar( Parameters::Get( )->GetUserLatitudeArray( ) );
                latitudeNcVar.putAtt( Constants::cUnitsString, Constants::cLatitudeVariableUnit );

                netCDF::NcDim monthlyTimeNcDim = monthlyGridOutputsNcFile.addDim( Constants::cTimeVariableNames[ 0 ], Parameters::Get( )->GetLengthOfSimulationInMonths( ) );
                netCDF::NcVar monthlyTimeNcVar = monthlyGridOutputsNcFile.addVar( Constants::cTimeVariableNames[ 0 ], netCDF::ncFloat, monthlyTimeNcDim );
                monthlyTimeNcVar.putVar( Parameters::Get( )->GetMonthlyTimeStepArray( ) );
                monthlyTimeNcVar.putAtt( Constants::cUnitsString, Constants::cMonthlyTimeUnitName );

                Types::NcDimVector dataDimensions;
                dataDimensions.push_back( monthlyTimeNcDim );
                dataDimensions.push_back( latitudeDim );
                dataDimensions.push_back( longitudeDim );

                for( Types::GridDatumVector::iterator iter = monthlyGridDatumVector.begin( ); iter != monthlyGridDatumVector.end( ); ++iter ) {
                    netCDF::NcVar monthlyGridDatumNcVar = monthlyGridOutputsNcFile.addVar( ( *iter )->GetName( ), netCDF::ncFloat, dataDimensions );
                    monthlyGridDatumNcVar.putAtt( Constants::cUnitsString, ( *iter )->GetDataUnit( ) );
                    monthlyGridDatumNcVar.putVar( ( *iter )->GetData( ) );
                }


            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }
    } 
}
//##

//##
void WriteModelState::WriteBasicOutputs( std::string output_directory ) {

    Types::BasicDatumMap basicDatumMap = DataRecorder::Get( )->GetBasicDatumMap( );

    // Separate monthly and annual basic datums
    if( basicDatumMap.size( ) > 0 ) {
        Types::BasicDatumVector annualBasicDatumVector;
        Types::BasicDatumVector monthlyBasicDatumVector;

        for( Types::BasicDatumMap::iterator iter = basicDatumMap.begin( ); iter != basicDatumMap.end( ); ++iter ) {
            if( iter->second->GetTimeUnit( ) == Constants::cAnnualTimeUnitName ) {
                annualBasicDatumVector.push_back( iter->second );
            } else if( iter->second->GetTimeUnit( ) == Constants::cMonthlyTimeUnitName ) {
                monthlyBasicDatumVector.push_back( iter->second );
            }
        }

        // Write annual basic datums
        if( annualBasicDatumVector.size( ) > 0 ) {
            std::string filePath = output_directory + Constants::cAnnualBasicOutputsFileName_State;

            try {
                netCDF::NcFile annualBasicOutputsNcFile( filePath, netCDF::NcFile::replace ); // Creates file
                netCDF::NcDim annualTimeNcDim = annualBasicOutputsNcFile.addDim( Constants::cTimeVariableNames[ 0 ], Parameters::Get( )->GetLengthOfSimulationInYears( ) ); // Creates dimension
                netCDF::NcVar annualTimeNcVar = annualBasicOutputsNcFile.addVar( Constants::cTimeVariableNames[ 0 ], netCDF::ncFloat, annualTimeNcDim ); // Creates variable
                annualTimeNcVar.putVar( Parameters::Get( )->GetAnnualTimeStepArray( ) );
                annualTimeNcVar.putAtt( Constants::cUnitsString, Constants::cAnnualTimeUnitName );

                Types::NcDimVector dataDimensions;
                dataDimensions.push_back( annualTimeNcDim );

                for( Types::BasicDatumVector::iterator iter = annualBasicDatumVector.begin( ); iter != annualBasicDatumVector.end( ); ++iter ) {
                    netCDF::NcVar annualBasicDatumNcVar = annualBasicOutputsNcFile.addVar( ( *iter )->GetName( ), netCDF::ncFloat, dataDimensions );
                    annualBasicDatumNcVar.putAtt( Constants::cUnitsString, ( *iter )->GetDataUnit( ) );
                    annualBasicDatumNcVar.putVar( ( *iter )->GetData( ) );
                }

            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }

        // Write monthly basic datums
        if( monthlyBasicDatumVector.size( ) > 0 ) {
            std::string filePath = output_directory + Constants::cMonthlyBasicOutputsFileName_State;

            try {
                netCDF::NcFile monthlyBasicOutputsNcFile( filePath, netCDF::NcFile::replace ); // Creates file
                netCDF::NcDim monthlyTimeNcDim = monthlyBasicOutputsNcFile.addDim( Constants::cTimeVariableNames[ 0 ], Parameters::Get( )->GetLengthOfSimulationInMonths( ) ); // Creates dimension
                netCDF::NcVar monthlyTimeNcVar = monthlyBasicOutputsNcFile.addVar( Constants::cTimeVariableNames[ 0 ], netCDF::ncFloat, monthlyTimeNcDim ); // Creates variable
                monthlyTimeNcVar.putVar( Parameters::Get( )->GetMonthlyTimeStepArray( ) );
                monthlyTimeNcVar.putAtt( Constants::cUnitsString, Constants::cMonthlyTimeUnitName );

                Types::NcDimVector dataDimensions;
                dataDimensions.push_back( monthlyTimeNcDim );

                for( Types::BasicDatumVector::iterator iter = monthlyBasicDatumVector.begin( ); iter != monthlyBasicDatumVector.end( ); ++iter ) {
                    netCDF::NcVar monthlyBasicDatumNcVar = monthlyBasicOutputsNcFile.addVar( ( *iter )->GetName( ), netCDF::ncFloat, dataDimensions );
                    monthlyBasicDatumNcVar.putAtt( Constants::cUnitsString, ( *iter )->GetDataUnit( ) );
                    monthlyBasicDatumNcVar.putVar( ( *iter )->GetData( ) );
                }

            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }
    }
}
//##