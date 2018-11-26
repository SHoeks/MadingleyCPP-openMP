#ifndef GRIDCELL
#define GRIDCELL

#include "UtilityFunctions.h"
#include "Stock.h"
#include "Cohort.h"
#include "Environment.h"
#include "Parameters.h"
#include "DataIndices.h"
#include "Convertor.h"

class GridCell {
public:
    GridCell( );

    void SetCellCoords( unsigned );

    double GetRealm( );
    bool IsMarine( );

    void InsertCohort( Cohort* );
    void RemoveCohort( Cohort* );
    void MoveCohort( Cohort* );
    void CheckGone( Cohort* );
    unsigned GetIndex( ) const;
    unsigned GetLatitudeIndex( ) const;
    unsigned GetLongitudeIndex( ) const;
    double GetCellArea( ) const;
    double GetCellHeight( ) const;
    double GetCellWidth( ) const;
    void SetCohortSize( unsigned );

    /** \brief Gets the number of cohorts in this grid cell */
    unsigned GetNumberOfCohorts( );

    static std::vector<Cohort*> mNewCohorts;
    #pragma omp threadprivate(mNewCohorts)

    template <typename F>
    void ApplyFunctionToAllCohorts( F f ) {
        for( int index = 0; index < mCohorts.size( ); index++ ) {
            // Work through the list of cohorts
            for( Cohort* c : mCohorts[ index ] ) {
                f( c );
            }
        }
    }

    template <typename F>
    void ApplyFunctionToAllCohortsWithStaticRandomness( F f, int CurrentTimeStep ) {
        std::vector<unsigned > RandomCohortOrder;
        std::vector<std::pair<int, int> > indexedList;
        unsigned TotalCohorts = 0;
        for( int functionalTypeIndex = 0; functionalTypeIndex < mCohorts.size( ); functionalTypeIndex++ ) {
            // Work through the list of cohorts and create a list of pairs so as to be able to lookup cohorts easily
            for( int cohortNum = 0; cohortNum < mCohorts[ functionalTypeIndex ].size( ); cohortNum++ ) {
                indexedList.push_back( std::make_pair( functionalTypeIndex, cohortNum ) );
                TotalCohorts++;
            }
        }
        //despite the name of this utility, it actually returns a random list, but with a given fixed seed
        RandomCohortOrder = mUtilities.NonRandomlyOrderedCohorts( TotalCohorts, CurrentTimeStep );

        for( int i = 0; i < RandomCohortOrder.size( ); i++ ) {
            Cohort* c = mCohorts[indexedList[RandomCohortOrder[i]].first][indexedList[RandomCohortOrder[i]].second];
            f( c );
        }
    }

    template <typename F>
    void ApplyFunctionToAllStocks( F f ) {
        for( int index = 0; index < mStocks.size( ); index++ ) {
            // Work through the list of cohorts
            for( Stock& s : mStocks[ index ] ) {
                f( s );
            }
        }
    }

    vector< vector<Cohort*> > mCohorts;

    Types::StocksMap mStocks;



private:
    void RandomizeCohorts( );

    UtilityFunctions mUtilities;
    unsigned mIndex;
    unsigned mLatitudeIndex;
    unsigned mLongitudeIndex;
    double mCellArea;
    double mCellHeightKm;
    double mCellWidthKm;

};
#endif
