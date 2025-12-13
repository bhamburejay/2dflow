#include "nlohmann/json.hpp"
#include <fstream>
#include <iostream>
#include <petsc.h>

int main(int argc, char **args) {

  PetscInitialize(&argc, &args, (char *)0, (char *)0);

  // This is an example of how to use the PetscOptions to get command line
  // options help. Calling -help produces a list of all options.
  double lambda = 1.0;
  PetscOptionsBegin(PETSC_COMM_WORLD, "myprefix_", "options", __FILE__);
  PetscCall(PetscOptionsReal("-lambda", "my parameter lambda", "", lambda,
                             &lambda, NULL));
  PetscOptionsEnd();

  // Example of reading an options file.
  // char inputfile[PETSC_MAX_PATH_LEN] = "test_petsc_inputs.json";
  // PetscBool found(PETSC_TRUE);
  // PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, inputfile, found);

  // Now the value of lambda is set from the command line if provided
  // -myprefix_lambda <value>
  std::cout << "Value of lambda: " << lambda << std::endl;

  // Even though an option was not "defined" above, we can still retrieve it if
  // it was provided on the command line or in an options file. The options file
  // can be specified with -options_file <filename>
  int input_parameter = 0;
  PetscOptionsGetInt(NULL, "", "-input_parameter", &input_parameter, NULL);
  std::cout << "Value of input_parameter: " << input_parameter << std::endl;

  PetscFinalize();
  return 0;
}