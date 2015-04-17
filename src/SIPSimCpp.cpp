#include <iostream>
#include <math.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
 
using namespace std;
using namespace boost::python;
 


double rand_norm(double mean, double stdev){
  boost::mt19937 rng(rand()); 

  boost::normal_distribution<> nd(mean, stdev);

  boost::variate_generator<boost::mt19937&, 
                           boost::normal_distribution<> > var_nor(rng, nd);
  
  double d = var_nor();

  return d;
}

double rand_norm_range(double mean, double stdev, double min, double max){
  int max_tries = 1000;
  int tries = 0;

  while (1){
    tries++;

    double x = rand_norm(mean, stdev);
    if ((x >= min) & (x <= max)){
      return x;
    }
    else if (tries >= max_tries){
      cout << "ERROR: exceeded max tries (n=" 
	   << max_tries 
	   << ") to find a random variable" << endl;
      exit(1);
    }    
  }
}


double add_diffusion_old(double frag_gc, double frag_len){
  double diff_coef = 44500;    
  double stdev = sqrt(diff_coef / frag_len);
  return frag_gc + rand_norm(0, stdev);
}


double add_diffusion(double frag_GC, double frag_len, 
		     double T = 298, double B = 1.195e9,
		     double G = 7.87e-10, int M = 882){
  /*
    Calculating diffusion in standard deviation of %G+C equivalents.
    Adding diffusion G+C (drawn from normal distribution with s.d.=calculated s.d.) 
    to input G+C value.
    Args:
    frag_BD = rho (buoyant density)
    frag_len = fragment length (bp)
    T = absolute temperature
    B = beta
    G = G coefficient (see Clay et al., 2003)
    M = molecular weight per base pair of dry cesium DNA
   */

  double frag_BD = frag_GC / 100 * 0.098 + 1.66;

  double R = 8.3145e7;    
  double GC_var = pow(100 / 0.098, 2) * ((frag_BD*R*T)/(pow(B,2)*G*M*frag_len));
  double GC_sd = sqrt(GC_var);

  return frag_GC + rand_norm(0, GC_sd);
  
}


double GC2BD(double frag_GC){  
  return frag_GC / 100.0 * 0.098 + 1.66;
}

 
double addIncorpBD(double frag_BD, double incorp_perc, double isoMaxBD){
  return incorp_perc / 100 * isoMaxBD + frag_BD;
}
 

BOOST_PYTHON_MODULE(SIPSimCpp)
{
  def("rand_norm", rand_norm);
  def("rand_norm_range", rand_norm_range);
  def("add_diffusion", add_diffusion);
  def("GC2BD", GC2BD);
  def("addIncorpBD", addIncorpBD);
}
