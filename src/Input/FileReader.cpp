#include "FileReader.h"

#include "Convertor.h"
#include "Constants.h"
#include "InputData.h"
#include "DataLayerSet.h"
#include "DataRecorder.h"
#include "Parameters.h"
#include "DataCoords.h"

#include <netcdf>

FileReader::FileReader( ) {

}

FileReader::~FileReader( ) {

}

bool FileReader::ReadFiles( ) {

    bool success = false;

    success = ReadInputParameters( );

    if( success == true )
        success = SetUpOutputVariables( );

    if( success == true )
        success = ReadInputDataFiles( );
    
    return success;
}

bool FileReader::ReadTextFile( const std::string& filePath ) {

    bool success = false;

    ClearMetadata( );

    std::cout << "Reading text file \"" << filePath << "\"..." << std::endl;
    std::ifstream fileStream( filePath.c_str( ), std::ios::in );

    if( fileStream.is_open( ) ) {
        std::string readLine;
        unsigned lineCount = 0;

        while( std::getline( fileStream, readLine ) ) {
            if( lineCount > 0 && readLine[ 0 ] != Constants::cTextFileCommentCharacter ) {
                mMetadata.push_back( Convertor::Get( )->StringToWords( readLine, Constants::cDataDelimiterValue ) );
            } else if( lineCount == 0 ) {
                mMetadataHeadings = Convertor::Get( )->StringToWords( readLine, Constants::cDataDelimiterValue );
            }
            ++lineCount;
        }
        success = true;
        fileStream.close( );
        DataRecorder::Get( )->AddInputFilePath( filePath );
    } else {
        std::cout << "File path \"" << filePath << "\" is invalid." << std::endl;
    }

    return success;
}

bool FileReader::ReadInputParameters( ) {

    bool success = ReadTextFile( Constants::cConfigurationDirectory + Constants::cInputParametersFileName );

    if( success == true )
        success = Parameters::Get( )->Initialise( mMetadata );

    return success;
}

bool FileReader::SetUpOutputVariables( ) {

    bool success = ReadTextFile( Constants::cConfigurationDirectory + Constants::cOutputVariablesFileName );

    if( success == true )
        success = DataRecorder::Get( )->Initialise( mMetadata );

    return success;
}

bool FileReader::ReadInputDataFiles( ) {

    bool success = ReadTextFile( Constants::cConfigurationDirectory + Constants::cInputDataFileName );

    if( success == true ) {
        if( mMetadata.size( ) > 0 ) {
            Types::InputDataPointer initialInputData = new InputData( );

            for( unsigned environmentalDataFileIndex = 0; environmentalDataFileIndex < mMetadata.size( ); ++environmentalDataFileIndex ) {

                std::string filePath = Parameters::Get( )->GetRootDataDirectory( );
                filePath.append( Convertor::Get( )->ToString( Parameters::Get( )->GetGridCellSize( ) ) );

                if(Parameters::Get( )->GetSpatialSelection( )!="australia" && Parameters::Get( )->GetSpatialSelection( )!="eurasiaafrica" &&
                   Parameters::Get( )->GetSpatialSelection( )!="eurasia" && Parameters::Get( )->GetSpatialSelection( )!="africa" &&
                   Parameters::Get( )->GetSpatialSelection( )!="eurasia1" && Parameters::Get( )->GetSpatialSelection( )!="eurasia2" &&
                   Parameters::Get( )->GetSpatialSelection( )!="northamerica" && Parameters::Get( )->GetSpatialSelection( )!="southamerica" &&
                   Parameters::Get( )->GetSpatialSelection( )!="northamerica" && Parameters::Get( )->GetSpatialSelection( )!="asia" &&
                   Parameters::Get( )->GetSpatialSelection( )!="northamerica" && Parameters::Get( )->GetSpatialSelection( )!="europe" ) {
                    filePath.append( "deg/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="australia") {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( Parameters::Get( )->GetSpatialSelection( ) );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="eurasiaafrica") {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( Parameters::Get( )->GetSpatialSelection( ) );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="africa") {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( Parameters::Get( )->GetSpatialSelection( ) );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="eurasia" || 
                   Parameters::Get( )->GetSpatialSelection( )=="europe" ||
                   Parameters::Get( )->GetSpatialSelection( )=="asia" ) {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( "eurasia" );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="eurasia1") {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( Parameters::Get( )->GetSpatialSelection( ) );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="eurasia2") {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( Parameters::Get( )->GetSpatialSelection( ) );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="northamerica") {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( Parameters::Get( )->GetSpatialSelection( ) );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }
                if(Parameters::Get( )->GetSpatialSelection( )=="southamerica") {
                    //std::cout << "file reader: "<< Parameters::Get( )->GetSpatialSelection( ) << std::endl;
                    filePath.append( "deg_" );
                    filePath.append( Parameters::Get( )->GetSpatialSelection( ) );
                    filePath.append( "/" );
                    filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                }

                std::cout << "Reading NetCDF file \"" << filePath << "\"..." << std::endl;

                try {
                    netCDF::NcFile inputNcFile( filePath, netCDF::NcFile::read ); // Open the file for read access
                    std::multimap< std::string, netCDF::NcVar > multiMap = inputNcFile.getVars( );

                    // Outer variable loop
                    for( std::multimap<std::string, netCDF::NcVar>::iterator it = multiMap.begin( ); it != multiMap.end( ); ++it ) {
                        std::string variableName = ( *it ).first;
                        netCDF::NcVar variableNcVar = ( *it ).second;
                        std::vector< netCDF::NcDim > varDims = variableNcVar.getDims( );

                        Types::UnsignedVector variableDimensions;
                        unsigned variableSize = 1;

                        // Inner variable dimension loop
                        for( unsigned dimIndex = 0; dimIndex < varDims.size( ); ++dimIndex ) {
                            variableDimensions.push_back( varDims[ dimIndex ].getSize( ) );
                            variableSize *= varDims[ dimIndex ].getSize( );
                        }

                        float* variableData = new float[ variableSize ];
                        variableNcVar.getVar( variableData );
                        bool isDefault = variableName == Convertor::Get( )->ToLowercase( mMetadata[ environmentalDataFileIndex ][ Constants::eDefaultVariableName ] );
                        initialInputData->AddVariableToDatum( mMetadata[ environmentalDataFileIndex ][ Constants::eInternalName ], variableName, variableDimensions, variableSize, variableData, isDefault );
                    }

                } catch( netCDF::exceptions::NcException& e ) {
                    e.what( );
                    std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
                }
            }

            DataLayerSet::Get( )->SetDataLayers( initialInputData );
            success = true;
        } else {
            success = false;
        }
    }

    return success;
}

void FileReader::ClearMetadata( ) {
    for( unsigned rowIndex = 0; rowIndex < mMetadata.size( ); ++rowIndex ) {
        mMetadata[ rowIndex ].clear( );
    }
    mMetadata.clear( );
    mMetadataHeadings.clear( );
}
