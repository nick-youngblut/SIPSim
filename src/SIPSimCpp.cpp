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


double calc_diffusion_simple(double frag_gc, double frag_len){
  /*
    Simple calculation of diffusion based on fragment length
  */
  double diff_coef = 44500;    
  double stdev = sqrt(diff_coef / frag_len);
  return frag_gc + rand_norm(0, stdev);
}


double calc_diffusion_GC(double frag_GC, double frag_len,
		         double T, double B, double G, int M){
  /*
    Calculating diffusion in standard deviation of %G+C equivalents.
    Adding diffusion G+C (drawn from normal distribution with s.d.=calculated s.d.) 
    to input G+C value.
    Args:
    frag_GC = G+C content of DNA fragment
    frag_len = DNA fragment length (bp)
    T = absolute temperature
    B = beta
    G = G coefficient (see Clay et al., 2003)
    M = molecular weight per base pair of dry cesium DNA
   */

  double frag_BD = frag_GC / 100 * 0.098 + 1.66;
  double R = 8.3145e7;    
  double GC_sd = sqrt(pow(100 / 0.098, 2) * ((frag_BD*R*T)/(pow(B,2)*G*M*frag_len)));

  return rand_norm(0, GC_sd);  
}

double calc_diffusion_BD(double frag_BD, double frag_len,
		         double T, double B, double G, int M){
  /*
    Calculating diffusion in standard deviation of buoyant_density equivalents (rho).
    Args:
    frag_BD = rho (buoyant density)
    frag_len = fragment length (bp)
    T = absolute temperature
    B = beta
    G = G coefficient (see Clay et al., 2003)
    M = molecular weight per base pair of dry cesium DNA
    Return:
    BD error due to diffusion value drawn from a normal distribution with a 
    standard deviation determined by calculated diffusion
   */

  double R = 8.3145e7;    
  double sd_BD = sqrt((frag_BD*R*T)/(pow(B,2)*G*M*frag_len));

  return rand_norm(0, sd_BD);  
}


double GC2BD(double GC){  
  /*
    Calaculate buoyant density from G+C.
    Args:
    GC = % GC of DNA fragment
  */
  return GC / 100.0 * 0.098 + 1.66;
}

 
double addIncorpBD(double frag_BD, double incorp_perc, double isoMaxBD){
  return incorp_perc / 100 * isoMaxBD + frag_BD;
}
 

BOOST_PYTHON_MODULE(SIPSimCpp)
{
  def("rand_norm", rand_norm);
  def("rand_norm_range", rand_norm_range);
  def("calc_diffusion_GC", calc_diffusion_GC);
  def("calc_diffusion_BD", calc_diffusion_BD);
  def("GC2BD", GC2BD);
  def("addIncorpBD", addIncorpBD);
}
