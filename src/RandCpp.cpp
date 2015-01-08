#include <iostream>
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
 

BOOST_PYTHON_MODULE(RandCpp)
{
  def("rand_norm", rand_norm);
  def("rand_norm_range", rand_norm_range);
}
