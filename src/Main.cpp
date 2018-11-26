#include <fenv.h>

#include <Madingley.h>

#include "FileReader.h"
#include "FileWriter.h"
#include "Date.h"

std::string SpatialSelectionFlag;
int NumberOfThreads = 0;
int NumberFlags;

int main( int argc, char* argv[] ) {
    //std::cout << argc << std::endl;
    NumberFlags = argc;
    if( argc > 1 ) SpatialSelectionFlag = argv[1];
    if( argc > 2 ) {NumberOfThreads = std::stoi(argv[2]);
        std::cout << "Threads: " << NumberOfThreads << std::endl;
    }

    //this line enables the gdb debugger to catch Nan or floating point problems
    feraiseexcept( FE_INVALID | FE_OVERFLOW );
    // Write out model details to the console
    std::cout << "Madingley model C++ v. 0." << std::endl;
    std::cout << "Model Run started at " << Date::GetDateAndTimeString( Constants::cCompleteDateFormat ) << std::endl;

    FileReader fileReader;

    if( fileReader.ReadFiles( ) == true ) {
        std::cout << "Files read successfully..." << std::endl << std::endl;
        FileWriter fileWriter;

        // Initialise the model
        // Declare an instance of the class that runs a Madingley model simulation
        Madingley madingleyModel;

        // Declare and start a timer
        Stopwatch s;
        s.Start( );

        // Run the simulation
        madingleyModel.Run( );

        if( fileWriter.WriteFiles( ) == true )
            std::cout << "Files written to \"" << fileWriter.GetOutputDirectory( ) << "\" successfully..." << std::endl;

        // Stop the timer and write out the time taken to run this simulation
        s.Stop( );
        std::cout << "Model run finished" << std::endl;
        std::cout << "Total elapsed time was " + Convertor::Get( )->ToString( s.GetElapsedTimeSecs( ) ) + " seconds." << std::endl;
    } else {
        std::cout << "ERROR> File reading failed. System exiting..." << std::endl;
    }
    return 0;
}
