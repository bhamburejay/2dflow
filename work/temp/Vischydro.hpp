#ifndef VISCHYDRO_HPP
#define VISCHYDRO_HPP

#include "EOS.hpp"
#include "VischydroNode.hpp"
#include "limitter.hpp"
#include <petscts.h>
#include <petscdmda.h>
#include "json/json.h"
#include <memory>

class Vischydro {
public:
    Vischydro(const Json::Value &inputs, const EOS &eos);
    ~Vischydro();

    // PETSc callback wrappers
    static PetscErrorCode EulerRHSFunction(TS, PetscReal, Vec, Vec, void*);
    static PetscErrorCode LHSIFunction(TS, PetscReal, Vec, Vec, Vec, void*);
    static PetscErrorCode LHSIJacobian(TS, PetscReal, Vec, Vec, PetscReal, Mat, Mat, void*);
    static PetscErrorCode PostStepInversion(TS);

    // Utility
    Json::Value get_inputs(const std::string& key) const;

    // Members (public or private as appropriate)
    const Json::Value &inputs;
    const EOS &eos;
    DM domain;
    Vec solution, solution_local, solution_last, Residual;
    Mat Jacobian;
    TS stepper;
    PetscViewer H5viewer;
    double xmin, xmax, dx;
};

#endif // VISCHYDRO_HPP
