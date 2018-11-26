#include "Environment.h"
#include "GridCell.h"
#include "ClimateVariablesCalculator.h"
#include "Types.h"
#include "FileReader.h"
#include "Parameters.h"
#include "DataLayerSet.h"
#include "DataIndices.h"
#include "TimeStep.h"
#include "Layer2D.h"
#include "Layer3D.h"

#include <string.h>

Types::EnvironmentPointer Environment::mThis = NULL;
Types::LayerMap Environment::mLayers;

Environment::Environment( ) {
    AddLayer2D( "Realm" );
    SetRealm( );

    AddLayer2D( "TerrestrialHANPP" );
    SetHANPP( );

    AddLayer3D( "uVel" );
    SetUVel( );
    AddLayer3D( "vVel" );
    SetVVel( );

    AddLayer3D( "Temperature" );
    SetTemperature( );
    AddLayer3D( "DiurnalTemperatureRange" );

    SetDiurnalTemperatureRange( );

    AddLayer3D( "Precipitation" );

    AddLayer2D( "TotalPrecip" );
    SetPrecipitation( );

    AddLayer3D( "NPP" );
    SetNPP( );

    AddLayer3D( "Seasonality" );
    SetNPPSeasonality( );
    AddLayer3D( "Breeding Season" );
    SetBreeding( );


    AddLayer2D( "Organic Pool" );
    SetOrganicPool( );

    AddLayer2D( "Respiratory CO2 Pool" );
    SetRespiratoryCO2Pool( );

    AddLayer2D( "AnnualTemperature" );
    AddLayer2D( "SDTemperature" );
    AddLayer3D( "ExpTDevWeight" );
    SetAVGSDTemp( );

    AddLayer3D( "AET" );
    AddLayer2D( "TotalAET" );
    AddLayer3D( "Fraction Month Frost" );
    AddLayer2D( "Fraction Year Frost" );
    AddLayer2D( "Fraction Year Fire" );
    SetFrostandFire( );

    //make sure all time dependent fields set to the start
    Update( 0 );
}

void Environment::Update( unsigned currentMonth ) {
    for( auto& L: mLayers ) L.second->SetTime( currentMonth );
}

void Environment::AddLayer2D( std::string name ) {
    mLayers[ name ] = new Layer2D( Parameters::Get( )->GetNumberOfGridCells( ) );
}

void Environment::AddLayer3D( std::string name ) {
    mLayers[ name ] = new Layer3D( 12, Parameters::Get( )->GetNumberOfGridCells( ) );
}

Environment* Environment::Get( ) {
    if( mThis == NULL ) {
        mThis = new Environment( );
    }
    return mThis;
}

double& Environment::Get( std::string name, unsigned cellIndex ) {
    if( mThis == NULL ) mThis = new Environment( );
    if( mLayers.count( name ) == 0 ) {
        std::cout << "Invalid Layer Request in Environment:: " << name << std::endl;
        exit( 0 );
    }
    return ( *mLayers[ name ] )[ cellIndex ];
}

double& Environment::Get( std::string name, GridCell& gcl ) {
    if( mThis == 0 )mThis = new Environment( );
    if( mLayers.count( name ) == 0 ) {
        std::cout << "Invalid Layer Request in Environment:: " << name << std::endl;
        exit( 0 );
    }
    return (*mLayers[ name ] )[ gcl.GetIndex( ) ];
}

void Environment::SetTemperature( ) {
    for( unsigned timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        TimeStep::Get( )->SetMonthly( timeIndex );
        mLayers[ "Temperature" ]->SetTime( timeIndex );

        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            double d = Constants::cMissingValue;
            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 1 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "MarineTemp", cellIndex );
            } else if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialTemp", cellIndex );
            }

            if( d == Constants::cMissingValue ) {
                d = 0;
                if(Parameters::Get( )->GetSpatialSelection( )=="longlat" || Parameters::Get( )->GetSpatialSelection( )=="global"){
                    //std::cout << "Warning Environment::setTemperature- missing values in temperature field!!" << std::endl;
                }
                
            }
            ( *mLayers[ "Temperature" ] )[ cellIndex ] = d;
        }
    }
}

void Environment::SetUVel( ) {
    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        TimeStep::Get( )->SetMonthly( timeIndex );
        mLayers[ "uVel" ]->SetTime( timeIndex );
        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            double d = Constants::cMissingValue;

            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 1 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "MarineEastVel", cellIndex );
            }
            ( *mLayers[ "uVel" ] )[ cellIndex ] = d;
        }
    }
}

void Environment::SetVVel( ) {
    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        TimeStep::Get( )->SetMonthly( timeIndex );
        mLayers[ "vVel" ]->SetTime( timeIndex );

        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            double d = Constants::cMissingValue;

            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 1 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "MarineNorthVel", cellIndex );
            }
            ( *mLayers[ "vVel" ] )[ cellIndex ] = d;
        }
    }
}

void Environment::SetDiurnalTemperatureRange( ) {
    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        TimeStep::Get( )->SetMonthly( timeIndex );
        mLayers[ "DiurnalTemperatureRange" ]->SetTime( timeIndex );
        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {

            double d = Constants::cMissingValue;
            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialDTR", cellIndex );
            }

            ( *mLayers[ "DiurnalTemperatureRange" ] )[ cellIndex ] = d;
        }
    }
}

void Environment::SetPrecipitation( ) {
    for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
        ( *mLayers[ "TotalPrecip" ] )[cellIndex] = 0;
    }
    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        TimeStep::Get( )->SetMonthly( timeIndex );
        mLayers[ "Precipitation" ]->SetTime( timeIndex );

        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            double d = 0; // There are missing values here. No marine precipitation data.

            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialPre", cellIndex );
            }

            if( d == Constants::cMissingValue ) {
                d = 0;
                if(Parameters::Get( )->GetSpatialSelection( )=="longlat" || Parameters::Get( )->GetSpatialSelection( )=="global"){
                    //std::cout << "Warning Environment::setPrecipitation- missing values in precipitation field!!" << std::endl;
                }
            }
            ( *mLayers[ "Precipitation" ] )[ cellIndex ] = d;
            ( *mLayers[ "TotalPrecip" ] )[ cellIndex ] += d;
        }
    }
}

void Environment::SetNPP( ) {
    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        TimeStep::Get( )->SetMonthly( timeIndex );
        mLayers[ "NPP" ]->SetTime( timeIndex );

        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            double d = Constants::cMissingValue;
            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 1 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "MarineNPP", cellIndex );
            } else if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialNPP", cellIndex );
            }

            if( d == Constants::cMissingValue ) {
                d = 0;
                if(Parameters::Get( )->GetSpatialSelection( )=="longlat" || Parameters::Get( )->GetSpatialSelection( )=="global"){
                    //std::cout << "Warning Environment::setNPP- missing values in NPP field!!" << std::endl;
                }
            }
            ( *mLayers[ "NPP" ] )[ cellIndex ] = d;
        }
    }
}
//------------------------------------------------------------------------------

void Environment::SetRealm( ) {
    for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
        float d = 0;
        if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 1 ) {
            d = 2.0;
        } else if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
            d = 1.0;
        }
        ( *mLayers[ "Realm" ] )[ cellIndex ] = d;
    }
}
//------------------------------------------------------------------------------

void Environment::SetOrganicPool( ) {

    if(Parameters::Get( )->GetApplyModelSpinup()==0){
        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            ( *mLayers[ "Organic Pool" ] )[ cellIndex ] = 0;
        }
    }
    if(Parameters::Get( )->GetApplyModelSpinup()==1){

        std::cout << "opening: "<< Parameters::Get( )->GetPoolCSVLocation( ) << std::endl;
        std::vector< std::vector<std::string> > v(3,std::vector<std::string>(0));
        std::ifstream file ( Parameters::Get( )->GetPoolCSVLocation( ) );
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
        std::cout << v.size() << std::endl; 
        std::cout << v[0].size() << std::endl; 
        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            //std::cout << v[0][cellIndex+1] << " | " << v[1][cellIndex+1] << " | " << v[2][cellIndex+1] << std::endl;  
            ( *mLayers[ "Organic Pool" ] )[ cellIndex ] = std::stod(v[1][cellIndex+1]);
        }
    }


 
}
//------------------------------------------------------------------------------

void Environment::SetRespiratoryCO2Pool( ) {
    if(Parameters::Get( )->GetApplyModelSpinup()==0){
        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            ( *mLayers[ "Respiratory CO2 Pool" ] )[ cellIndex ] = 0;
        }
    }
    if(Parameters::Get( )->GetApplyModelSpinup()==1){
        std::cout << "opening: "<< Parameters::Get( )->GetPoolCSVLocation( ) << std::endl;
        std::vector< std::vector<std::string> > v(3,std::vector<std::string>(0));
        std::ifstream file ( Parameters::Get( )->GetPoolCSVLocation( ) );
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
        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            ( *mLayers[ "Respiratory CO2 Pool" ] )[ cellIndex ] = std::stod(v[2][cellIndex+1]);
        }

    }
}
//------------------------------------------------------------------------------

void Environment::SetAVGSDTemp( ) {
    for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
        double avg = 0;
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {

            double d = Constants::cMissingValue;

            TimeStep::Get( )->SetMonthly( timeIndex );
            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 1 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "MarineTemp", cellIndex );
            } else if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialTemp", cellIndex );
            }

            if( d == Constants::cMissingValue ) d = 0;
            avg += d;
        }
        avg = avg / 12;
        double sota = 0, sumExp = 0;

        std::vector< double > exptdev( 12 );
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
            double d = Constants::cMissingValue;
            TimeStep::Get( )->SetMonthly( timeIndex );

            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 1 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "MarineTemp", cellIndex );
            } else if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialTemp", cellIndex );
            }

            if( d == Constants::cMissingValue )d = 0;
            sota += ( d - avg )*( d - avg );
            exptdev[ timeIndex ] = exp( -( d - avg ) / 3 );
            sumExp += exptdev[ timeIndex ];
        }
        for( int tm = 0; tm < 12; tm++ ) {
            mLayers[ "ExpTDevWeight" ]->SetTime( tm );
            ( *mLayers[ "ExpTDevWeight" ] )[ cellIndex ] = exptdev[tm] / sumExp;
        }

        ( *mLayers[ "AnnualTemperature" ] )[cellIndex] = avg;
        ( *mLayers[ "SDTemperature" ] )[ cellIndex ] = sqrt( sota / 12 );
    }
}
//----------------------------------------------------------------------------------------------

/** \brief Calculate monthly seasonality values of Net Primary Production - ignores missing values. If there is no NPP data (ie all zero or missing values)
then assign 1/12 for each month.
 */
void Environment::SetNPPSeasonality( ) {
    for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
        // Loop over months and calculate total annual NPP
        double totalNPP = 0.0;
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
            mLayers["NPP"]->SetTime( timeIndex );
            double N = ( *mLayers["NPP"] )[ cellIndex ];
            if( N != Constants::cMissingValue && N > 0 ) totalNPP += N;
        }
        if( totalNPP == 0 ) {
            // Loop over months and calculate seasonality
            // If there is no NPP value then assign a uniform flat seasonality
            for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
                mLayers["Seasonality"]->SetTime( timeIndex );
                ( *mLayers["Seasonality"] )[ cellIndex ] = 1.0 / 12.0;
            }

        } else {
            // Some NPP data exists for this grid cell so use that to infer the NPP seasonality
            // Loop over months and calculate seasonality
            for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
                mLayers["NPP"]->SetTime( timeIndex );
                mLayers["Seasonality"]->SetTime( timeIndex );
                double N = ( *mLayers["NPP"] )[ cellIndex ];
                if( N != Constants::cMissingValue && N > 0 ) {
                    ( *mLayers["Seasonality"] )[ cellIndex ] = N / totalNPP;
                } else {
                    ( *mLayers["Seasonality"] )[ cellIndex ] = 0.0;
                }
            }
        }

    }
}
//----------------------------------------------------------------------------------------------

void Environment::SetFrostandFire( ) {
    // Calculate other climate variables from temperature and precipitation
    // Declare an instance of the climate variables calculator
    ClimateVariablesCalculator CVC;

    for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
        // Calculate the fraction of the year that experiences frost
        std::vector< double > FrostDays( 12 ), Temperature( 12 ), Precipitation( 12 );
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {

            TimeStep::Get( )->SetMonthly( timeIndex );
            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                FrostDays[timeIndex] = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialFrost", cellIndex );
                Precipitation[timeIndex] = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialPre", cellIndex );
                Temperature[timeIndex] = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialTemp", cellIndex );
            }
        }
        ( *mLayers[ "Fraction Year Frost" ] )[ cellIndex ] = CVC.GetNDF( FrostDays, Temperature, Constants::cMissingValue );

        std::vector< double > MonthDays = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
            mLayers[ "Fraction Month Frost" ]->SetTime( timeIndex );
            ( *mLayers[ "Fraction Month Frost" ] )[ cellIndex ] = std::min( FrostDays[timeIndex] / MonthDays[timeIndex], ( double )1.0 );
        }
        double AWC = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialAWC", cellIndex );

        std::tuple< std::vector< double >, double, double > TempTuple = CVC.MonthlyActualEvapotranspirationSoilMoisture( AWC, Precipitation, Temperature );
        ( *mLayers[ "TotalAET" ] )[cellIndex ] = 0;
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
            mLayers[ "AET" ]->SetTime( timeIndex );
            ( *mLayers[ "AET" ] )[ cellIndex ] = std::get< 0 >( TempTuple )[ timeIndex ];
            ( *mLayers[ "TotalAET" ] )[ cellIndex ] += std::get< 0 >( TempTuple )[ timeIndex ];
        }
        ( *mLayers["Fraction Year Fire"] )[ cellIndex ] = ( std::get< 2 > ( TempTuple ) / 360 );
    }
}
//----------------------------------------------------------------------------------------------

void Environment::SetBreeding( ) {
    // Designate a breeding season for this grid cell, where a month is considered to be part of the breeding season if its NPP is at
    // least 80% of the maximum NPP throughout the whole year
    for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
        double maxSeason = -1;
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
            mLayers[ "Seasonality" ]->SetTime( timeIndex );
            maxSeason = std::max( maxSeason, ( *mLayers[ "Seasonality" ] )[ cellIndex ] );
        }
        for( int i = 0; i < 12; i++ ) {
            mLayers[ "Seasonality" ]->SetTime( i );
            mLayers[ "Breeding Season" ]->SetTime( i );

            if( ( *mLayers[ "Seasonality" ] )[ cellIndex ] / maxSeason > 0.5 ) {
                ( *mLayers[ "Breeding Season" ] )[ cellIndex ] = 1.0;
            } else {
                ( *mLayers[ "Breeding Season" ] )[ cellIndex ] = 0.0;
            }
        }
    }
}
//------------------------------------------------------------------------------

void Environment::SetHANPP( ) {
    for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
        ( *mLayers[ "TerrestrialHANPP" ] )[ cellIndex ] = 0;
    }

    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        TimeStep::Get( )->SetMonthly( timeIndex );
        mLayers[ "TerrestrialHANPP" ]->SetTime( timeIndex );

        for( unsigned cellIndex = 0; cellIndex < Parameters::Get( )->GetNumberOfGridCells( ); cellIndex++ ) {
            double d = 0; // There are missing values here. No marine precipitation data.

            if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", cellIndex ) == 2 ) {
                d = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialHANPP", cellIndex );
            }

            if( d == Constants::cMissingValue ) {
                d = 0;
                //std::cout << "Warning Environment::setHANPP- missing values in TerrestrialHANPP field!!" << std::endl;
            }
            ( *mLayers[ "TerrestrialHANPP" ] )[ cellIndex ] = d;
        }
    }
}