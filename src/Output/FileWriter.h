#ifndef FILEWRITER
#define	FILEWRITER

#include "Types.h"

class FileWriter {
public:
    FileWriter( );
    ~FileWriter( );
    
    bool WriteFiles( );
    
    std::string& GetOutputDirectory( );
    
private:
    void InitialiseOutputDirectory( );
    void WriteInputFiles( );
    bool WriteBasicOutputs( ) const;
    bool WriteGridOutputs( ) const;
    
    std::string mOutputDirectory;
};

#endif

