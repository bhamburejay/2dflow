#include <iostream>
#include <fstream>
#include "Vischydro.hpp"


void RunCode() {
  // Json setup
  Json::Value config ;
  std::ifstream in("2dflow_input.json"); 
  in >> config  ; 
  std::cout << config ;

  std::unique_ptr<EOS> eos ; 
  Vischydro hydro(config, eos.get()) ;

}
int main(int argc, char **argv) 
{
  std::cout << "Hello World" << std::endl;
  int ierr ;
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  RunCode() ;
  ierr = PetscFinalize(); CHKERRQ(ierr);

}

