#ifndef PATCH_H_
#define PATCH_H_

#include <iostream>
#include <memory>
#include <vector>

#include "traits.h"
#include "population.h"

class patch
{
  
 public:

  // members
  std::vector<population> populations;
  std::vector<double> propensities;
  double patch_prop;
  bool active;

  // constructors and destructors
  patch() {
    // std::cout << "patch::patch()" << std::endl;
  }
  
  ~patch() {
    // std::cout << "patch::~patch()" << std::endl;
  }

  /*
  patch(const patch &p2){
    std::cout << "copy patch" << std::endl;
  }
  */

  // member functions
  void initialise(const trait_set& init_ts, const int disc_thr, const int init_n); //! initialises the first
  double patch_propensities();
  void advance_all_continuous(const double dt);
  bool deactivate();
  void add_population(const trait_set& ts, const int disc_thr);

  bool check_continuous() const;
  void print_propensities() const;
  double count_total_population() const;
  
};

#endif
