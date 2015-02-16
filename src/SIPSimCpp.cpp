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
      cout << "ERROR: exceeded max tries (n=" << max_tries << ") to find a random variable" << endl;
      exit(1);
    }    
  }
}


double add_diffusion(double frag_gc, double frag_len){
  double diff_coef = 44500;    
  double stdev = sqrt(diff_coef / frag_len);
  return frag_gc + rand_norm(0, stdev);
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
