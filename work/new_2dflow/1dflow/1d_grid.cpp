#include "vischydro.hpp"
#include <iostream>
#include <fstream>
#include "json/json.h"
#include <string>

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    // Local scope for PETSc objects
    {
        // 1) Read inputs.json
        Json::Value input;
        std::ifstream in("inputs.json");
        Json::CharReaderBuilder reader;
        std::string errs;
        if (!Json::parseFromStream(reader, in, &input, &errs)) {
            std::cerr << "JSON parse error: " << errs << std::endl;
            PetscFinalize();
            return 1;
        }

        // 2) Initialize EOS and Vischydro
        EOS eos;
        Vischydro visch(input, eos);

        // 3) Configure time-stepping monitor
        TSMonitorSet(visch.stepper, VischydroMonitor, &visch, NULL);
        PetscViewerHDF5PushTimestepping(visch.H5viewer);

        // 4) Run time evolution
        TSSolve(visch.stepper, visch.solution);

        // 5) Write final state
        PetscViewerHDF5PopTimestepping(visch.H5viewer);
        PetscObjectSetName((PetscObject)visch.solution, "solution");
        VecView(visch.solution, visch.H5viewer);
    }
    
    PetscFinalize();
    return 0;
}
