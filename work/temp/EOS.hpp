#ifndef EOS_HPP
#define EOS_HPP
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class EOS {
public:
   double Nc;
   double Nf;
   
   EOS(double Nc_in=3, double Nf_in=0) : Nc(Nc_in),  Nf(Nf_in) {;}
   ~EOS() {}
    
   void   initialize_eos();
   double get_cs2        (double e) const {return(1./3.);}
   double get_pressure   (double e) const {return(1./3.*e);}
   double get_temperature(double e) const {return pow(90.0/M_PI/M_PI*(e/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf), .25);}
};
#endif