#include "population.h"

void population::initialise(const trait_set& init_ts, const int discrete_threshold, const int init_n)
{
  if(init_n < discrete_threshold){
    n = init_n;
    m = 0.0;
    discrete = true;
  }else{
    n = 0;
    m = static_cast<double>(init_n);
    discrete = false;
  }
  traits = init_ts;
  // reactions: birth death migration mutation
  for(int i = 0; i < num_reactions; i++){
    propensities.at(i) = 0.0;
  }
  threshold = discrete_threshold;
  active = true;
}

double population::population_propensities(const double patch_total_pop)
{
  if(discrete){
    double x = static_cast<double>(n);
    if((traits.b - traits.d) > 0.0){
      propensities.at(0) =
	(1.0 - traits.g)*x*(traits.b - (traits.b - traits.t)*patch_total_pop/traits.k);
      propensities.at(1) = x*(traits.d + (traits.t - traits.d)*patch_total_pop/traits.k);
      propensities.at(2) = traits.m*x;
      propensities.at(3) = traits.g*x*(traits.b - (traits.b - traits.t)*patch_total_pop/traits.k);
    }else{
      propensities.at(0) = (1.0 - traits.g)*x*traits.b;
      propensities.at(1) = x*traits.d;
      propensities.at(2) = traits.m*x;
      propensities.at(3) = traits.g*x*traits.b;
    }
  }else{ // uses m
    if((traits.b - traits.d) > 0.0){
      propensities.at(0) = 0.0;
      propensities.at(1) = 0.0;
      propensities.at(2) = traits.m*m;
      propensities.at(3) = traits.g*m*(traits.b - (traits.b - traits.t)*patch_total_pop/traits.k);
    }else{
      propensities.at(0) = 0.0;
      propensities.at(1) = 0.0;
      propensities.at(2) = traits.m*m;
      propensities.at(3) = traits.g*m*traits.b;
    }
  }
  return (propensities.at(0) + propensities.at(1) + propensities.at(2) + propensities.at(3));
}

void population::advance(const double dt, const double patch_total_pop){
  double birthrate;
  double deathrate;
  if((traits.b - traits.d) > 0.0){
    birthrate = (1.0 - traits.g)*(traits.b - (traits.b - traits.t)*patch_total_pop/traits.k);
    deathrate = (traits.d + (traits.t - traits.d)*patch_total_pop/traits.k);
  }else{
    birthrate = (1.0 - traits.g)*traits.b;
    deathrate = traits.d;
  }
  m = m + dt*(birthrate*m - deathrate*m);
}

// TODO: this should be called from simulation::simulate
// via patch::population_extinction
bool population::deactivate()
{
  bool ok = true;
  if(!discrete) ok = false;
  if(n != 0) ok = false;
  active = false;
  m = 0.0;
  for(auto &i : propensities) i = 0.0; 
  return ok;
}

// bool if change disc -> cont
void population::plus()
{
  if(discrete){
    n += 1;
    if(n > threshold){
      m = static_cast<double>(n);
      n = 0;
      discrete = false;
    }
  }else{
    m += 1.0;
  }
}

// TODO: test this
// TODO: should not call deactivate()
bool population::minus()
{
  if(discrete){
    n -= 1;
    if(n < 1){
      if(!deactivate()){
	std::cout << " error : population.deactivate() returned false" << std::endl;
      }
      return true;
    }
  }else{
    m -= 1.0;
    if(m < threshold){
      n = static_cast<int>(m);
      m = 0.0;
      discrete = true;
    }
  }
  return false;
}

double population::read_n() const
{
  if(discrete){
    return static_cast<double>(n);
  }else{
    return m;
  }
}
  
void population::print_propensities() const
{
  std::cout << "pop props: " <<
    propensities.at(0) << " " <<
    propensities.at(1) << " " <<
    propensities.at(2) << " " <<
    propensities.at(3) << std::endl;
}
