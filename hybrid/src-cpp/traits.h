#ifndef TRAITS_H_
#define TRAITS_H_

#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include <assert.h>
#include "json/json.h"

/**
 * contains simulation specifications common to all runs
 */
struct par
{
public:
  
  int num_sims;     //!< number of simulations (with different seeds) to spawn
  int threshold;    //!< pop size threshold to switch discrete -> continuous
  double init_n;    //!< first patch initial population size
  double target;    //!< stop condition: population size reached
  double max_time;  //!< stop condition: time elapsed
  int max_events;   //!< stop condition: total number of discrete events elapsed
  int res;          //!< writing time-series: write only at every res step
  double dt;        //!< largest possible dt step for continuous populations
  int printmode;    //!< 0: full time series 1: ntot time series, 2: end state time, 3: end time
  std::string pre;  //!< prefix for output filez
  std::string traitfile;    //!< trait file name
  std::vector<int> seeds;   //!< array of seed values, one for each sim

  bool read_sim_json(const std::string sim_json_file); //!< populate from json
};

struct base_traits
{
public:

  // members
  double b_primary;    //!< primary population birth rate
  double d_primary;    //!< primary population death rate
  double b_metastatic; //!< metastasis birth rate
  double d_metastatic; //!< metastasis death rate
  double b_mutant;     //!< mutant birth rate
  double d_mutant;     //!< mutant death rate
  double t;            //!< turnover, note b > t > d
  double k_primary;    //!< primary population carrying capacity
  double k_metastatic; //!< metastasis carrying capacity
  double k_mutant;     //!< mutant carrying capacity
  double g;            //!< mutation probability
  double m;            //!< migration rate

  // constructors and destructors
  base_traits() {
    // std::cout << "trait_set::trait_set()" << std::endl;
  }

  ~base_traits() {
    // std::cout << "trait_set::~trait_set()" << std::endl;
  }

  // copy
  base_traits(const base_traits& bt2){
    b_primary    = bt2.b_primary;    d_primary    = bt2.d_primary;
    b_metastatic = bt2.b_metastatic; d_metastatic = bt2.d_metastatic;
    b_mutant     = bt2.b_mutant;     d_mutant     = bt2.d_mutant;
    t = bt2.t;
    k_primary    = bt2.k_primary;
    k_metastatic = bt2.k_metastatic;
    k_mutant     = bt2.k_mutant;
    g = bt2.g; m = bt2.m;
    // std::cout << "trait_set::trait_set(const trait_set& bt2)" << std::endl;
  }

  bool read_traits_from_file(const std::string trait_file); //!< populate from json
}; 

/**
 * contains parameter values for a single simulation run
 */
struct trait_set
{ 
public:

  // members
  double b;    //!< primary population birth rate
  double d;    //!< primary population death rate
  double t;            //!< turnover, note b > t > d
  double k;    //!< primary population carrying capacity
  double g;            //!< mutation probability
  double m;            //!< migration rate

  // constructors and destructors
  trait_set() {
    // std::cout << "trait_set::trait_set()" << std::endl;
  }

  ~trait_set() {
    // std::cout << "trait_set::~trait_set()" << std::endl;
  }

  // copy
  trait_set(const trait_set& ts2){
    b = ts2.b;
    d = ts2.d;
    t = ts2.t;
    k = ts2.k;
    g = ts2.g;
    m = ts2.m;
    // std::cout << "trait_set::trait_set(const trait_set& ts2)" << std::endl;
  }

  //bool read_traits_from_file(const std::string trait_file); //!< populate from json
  
};

trait_set trait_set_from_base(const base_traits& base, const bool meta, const bool muta);

#endif
