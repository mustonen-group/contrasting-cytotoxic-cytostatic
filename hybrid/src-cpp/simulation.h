#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <iostream>
#include <random>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "patch.h"
#include "traits.h"
#include "pcg_random.hpp"

class simulation
{
  
 public:

  // members
  int index;                        //!< simulation index for output writing
  double time;                      //!< keeps track of elapsed time
  std::vector<patch> patches;       //!< container for patch* objects
  std::vector<double> propensities; //!< container for patch propensities
  int num_active;                   //!< keeps count of active patches
  int num_events;                   //!< keeps count of elapsed events
  int num_mutations;                //!< keeps count of occurred mutations
  base_traits base_trait_set;       //!< new
  pcg32 rng;                        //!< a simulation object has its own rng, seed from main
  std::string outfilename;          //!< a simulation writes to outfile
  par pars;                         //!< pars common to all sims (threshold, stop conditions, etc.)

  // constructors and destructors
  simulation() {
    // std::cout << "simulation::simulation" << std::endl;
  }
  
  ~simulation() {
    // std::cout << "simulation::~simulation" << std::endl;
  }

  // member functions
  void initialise(base_traits base_ts,
		  const std::string outfile, par& p, const int ind);
  int simulate(); //!< a simulation object handles one simulation and writes to outfilename, return stop condition
  void advance_all_continuous(const double dt);
  double calculate_propensities();                  //!< returns total propensity
  void add_patch(trait_set ntraits); //, const int patch_choice, const int population_choice); // use base_traits
  //trait_set mutate(const trait_set &ts); // change this to use base traits

  template<typename T>
  int linear_sampler(const T& v); //!< selects and index from a probability vector or array

  bool check_continuous() const;  //!< returns true if non-discrete populations exist
  void print_propensities() const;
  void print_patch_count() const;
  double count_total_population() const;
  
  void print_state() const;
  void write_end_state(std::ofstream &out) const;
  void write_end_time(std::ofstream &out) const;
  void write_summary(std::ofstream &o, const int id) const; //!< printmode 1
  void write_all(std::ofstream &o, const int id) const; //!< used for printmode 0, id for cont, disc
  
};

#endif
