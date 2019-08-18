#include "patch.h"

void patch::initialise(const trait_set& init_ts, const int disc_thr, const int init_n)
{
  populations.push_back(population());
  propensities.push_back(0.0);
  patch_prop = 0.0;
  active = true;
  // initialise the first population
  populations.at(0).initialise(init_ts,disc_thr,init_n);
}

double patch::patch_propensities()
{
  double total_pop = count_total_population();
  patch_prop = 0.0; // member
  for(size_t i = 0; i < populations.size(); i++){ // auto?
    if(populations.at(i).active){
      propensities.at(i) = populations.at(i).population_propensities(total_pop);
      patch_prop += propensities.at(i);
    }
  }
  return patch_prop;
}

void patch::advance_all_continuous(const double dt)
{
  // check if any are !discrete
  double total_pop = count_total_population();
  for(auto &pop : populations){
    if(!pop.discrete){
      pop.advance(dt,total_pop);
    }
  }
}

bool patch::deactivate(){
  bool ok = true;
  for(auto &pop : populations){
    if(pop.active) ok = false;
  }
  for(auto &i : propensities) i = 0.0;
  patch_prop = 0.0;
  active = false;
  return ok;
}

void patch::add_population(const trait_set& ts, const int disc_thr)
{
  populations.push_back(population());
  propensities.push_back(0.0);
  populations.at(populations.size()-1).initialise(ts, disc_thr, 1);
}

bool patch::check_continuous() const
{
  bool any_cont = false;
  for(auto &pop : populations){
    if(!pop.discrete) any_cont = true;
  }
  return any_cont;
}

void patch::print_propensities() const
{
  for(auto &pop : populations){
    pop.print_propensities();
  }
}

double patch::count_total_population() const
{
  // std::accumulate ?
  double totpop = 0.0;
  for(auto &pop : populations){
    if(pop.active){
      if(pop.discrete){
	totpop += static_cast<double>(pop.n);
      }else{
	totpop += pop.m;
      }
    }
  }
  return totpop;
}
