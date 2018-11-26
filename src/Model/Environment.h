#ifndef ENVIRONMENT
#define	ENVIRONMENT

#include "Types.h"

class Environment {
public:
    static Environment* Get( );
    static double& Get( std::string, GridCell& );
    static double& Get( std::string, unsigned );
    static void Update( unsigned );
    
    Environment( );

    void AddLayer2D( std::string );
    void AddLayer3D( std::string );
    void SetUVel( );
    void SetVVel( );
    void SetTemperature( );
    void SetDiurnalTemperatureRange( );
    void SetPrecipitation( );
    void SetNPP( );
    void SetRealm( );
    void SetOrganicPool( );
    void SetRespiratoryCO2Pool( );
    void SetAVGSDTemp( );
    void SetNPPSeasonality( );
    void SetBreeding( );
    void SetFrostandFire( );
    void SetHANPP( );
    
    static Types::EnvironmentPointer mThis;
    static Types::LayerMap mLayers;
};

#endif

