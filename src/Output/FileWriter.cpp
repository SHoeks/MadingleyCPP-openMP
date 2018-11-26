#include "FileWriter.h"

#include "Constants.h"
#include "Date.h"
#include "Convertor.h"
#include "Parameters.h"
#include "DataRecorder.h"
#include "BasicDatum.h"
#include "GridDatum.h"

FileWriter::~FileWriter( ) {
}

FileWriter::FileWriter( ) {
    InitialiseOutputDirectory( );
    WriteInputFiles( );
}

bool FileWriter::WriteFiles( ) {
    bool completedSuccessfully = WriteBasicOutputs( );

    if( completedSuccessfully == true ) {
        completedSuccessfully = WriteGridOutputs( );
    }

    return completedSuccessfully;
}

std::string& FileWriter::GetOutputDirectory( ) {
    return mOutputDirectory;
}

void FileWriter::InitialiseOutputDirectory( ) {
    mOutputDirectory = Constants::cOutputBaseDirectory;
    mkdir( mOutputDirectory.c_str( ), Constants::cOutputFolderPermissions );
    mOutputDirectory.append( Date::GetDateAndTimeString( ) );
    mOutputDirectory.append( "_" );
    mOutputDirectory.append( Parameters::Get( )->GetSpatialSelection( ) );
    mkdir( mOutputDirectory.c_str( ), Constants::cOutputFolderPermissions );
    mOutputDirectory.append( Convertor::Get( )->ToString( Constants::cFolderDelimiter ) );
}

void FileWriter::WriteInputFiles( ) {
    Types::StringVector inputFilePaths = DataRecorder::Get( )->GetInputFilePaths( );

    for( unsigned stringIndex = 0; stringIndex < inputFilePaths.size( ); ++stringIndex ) {

        std::ifstream sourceFileStream( inputFilePaths[ stringIndex ].c_str( ), std::ios::in );

        std::string outputFilePath = mOutputDirectory;
        Types::StringVector inputFilePathComponents = Convertor::Get( )->StringToWords( inputFilePaths[ stringIndex ], Constants::cFolderDelimiter );

        std::string fileName = inputFilePathComponents[ inputFilePathComponents.size( ) - 1 ];
        outputFilePath.append( fileName );

        std::ofstream destinationFileStream( outputFilePath.c_str( ), std::ios::out );

        destinationFileStream << sourceFileStream.rdbuf( );

        sourceFileStream.close( );
        destinationFileStream.close( );
    }
}

bool FileWriter::WriteBasicOutputs( ) const {

    bool completedSuccessfully = false;

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
            std::string filePath = mOutputDirectory + Constants::cAnnualBasicOutputsFileName;

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

                completedSuccessfully = true;

            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }

        // Write monthly basic datums
        if( monthlyBasicDatumVector.size( ) > 0 ) {
            std::string filePath = mOutputDirectory + Constants::cMonthlyBasicOutputsFileName;

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

                completedSuccessfully = true;

            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }

    } else {
        completedSuccessfully = true;
    }

    return completedSuccessfully;
}

bool FileWriter::WriteGridOutputs( ) const {

    bool completedSuccessfully = false;

    Types::GridDatumMap gridDatumMap = DataRecorder::Get( )->GetGridDatumMap( );

    // Separate monthly and annual grid datums
    if( gridDatumMap.size( ) > 0 ) {
        Types::GridDatumVector annualGridDatumVector;
        Types::GridDatumVector monthlyGridDatumVector;

        for( Types::GridDatumMap::iterator iter = gridDatumMap.begin( ); iter != gridDatumMap.end( ); ++iter ) {
            if( iter->second->GetTimeUnit( ) == Constants::cAnnualTimeUnitName ) {
                annualGridDatumVector.push_back( iter->second );
            } else if( iter->second->GetTimeUnit( ) == Constants::cMonthlyTimeUnitName ) {
                monthlyGridDatumVector.push_back( iter->second );
            }
        }

        // Write annual grid datums
        if( annualGridDatumVector.size( ) > 0 ) {
            std::string filePath = mOutputDirectory + Constants::cAnnualGridOutputsFileName;

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

                completedSuccessfully = true;

            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }

        // Write monthly grid datums
        if( monthlyGridDatumVector.size( ) > 0 ) {
            std::string filePath = mOutputDirectory + Constants::cMonthlyGridOutputsFileName;

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

                completedSuccessfully = true;

            } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
            }
        }
    } else {
        completedSuccessfully = true;
    }

    return completedSuccessfully;
}
