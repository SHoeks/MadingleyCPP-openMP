#ifndef WMS_H 
#define WMS_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "Grid.h"
#include "GridCell.h"

#include "Environment.h"
#include "TimeStep.h"
#include "Maths.h"
#include "Parameters.h"
#include "DataRecorder.h"
#include "Types.h"


class WriteModelState {
public:
    void CohortSpinUpOutput( Grid&, std::string, unsigned );
    void StockSpinUpOutput( Grid&, std::string, unsigned );
    void CohortSpinUpOutputCheck( Grid&, std::string );
    void StockSpinUpOutputCheck( Grid&, std::string );
    double WriteState( Grid&, std::string, unsigned, unsigned, double, std::vector< std::vector<double> >, std::vector< std::vector<double> > );
    void WriteTimeLineBiomass( std::string, unsigned, std::vector< std::vector<double> >, std::vector< std::vector<double> > );
    void WriteGridOutputs( std::string );
    void WriteBasicOutputs( std::string );
private:
    unsigned temp;
};
#endif 