#ifndef POPULATION_H_
#define POPULATION_H_

#include <iostream>
#include <memory>
#include <array>

#include "traits.h"

const int num_reactions = 4; //! birth death migration mutation

class population
{

 public:

  // members
  
  int n;         //! holds population size if discrete
  double m;      //! holds population size if continuous
  int threshold; //! threshold for discrete/continuous

  trait_set traits;
  std::array<double,num_reactions> propensities;

  bool discrete; //! small population size -> discrete reactions
  bool active;   //! n > 0 -> active

  // constructors and destructors (these are never copied)
  
  population() {
    // std::cout << "population::population()" << std::endl;
  }

  ~population() {
    // std::cout << "population::~population()" << std::endl;
  }

  // member functions

  void initialise(const trait_set& init_ts, const int discrete_threshold, const int init_n);
  double population_propensities(const double patch_total_pop);
  void advance(const double dt, const double patch_total_pop);
  bool deactivate();
  // discrete -> continuous
  // add one (birth event) disc vs. cont
  // remove one (death, migr, mut event) disc vs. cont
  void plus();  //! adds 1 or 1.0 to population
  bool minus(); //! reduces 1 or 1.0 from population

  double read_n() const;
  void print_propensities() const;

};


#endif
