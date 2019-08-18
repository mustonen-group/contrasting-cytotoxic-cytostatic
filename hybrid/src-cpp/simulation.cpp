#include "simulation.h"

void simulation::initialise(base_traits base_ts, const std::string outfile,
			    par &p, const int ind)
{
  index = ind;
  pars = p; /// calls copy constructor
  time = 0.0;
  patches.push_back(patch());
  propensities.push_back(0.0);
  num_active = 1;
  num_events = 0;
  num_mutations = 0;
  base_trait_set = base_ts; /// calls copy constructor
  trait_set init_ts = trait_set_from_base(base_ts,false,false);
  patches.at(0).initialise(init_ts,pars.threshold,pars.init_n); /// initialise the first patch
  rng.seed(pars.seeds.at(index));
  outfilename = outfile;
}

/**
 * should be run once per simulation
 */
int simulation::simulate()
{
  static std::uniform_real_distribution<double> unif(0.0,1.0);
  std::ofstream out;

  // printmodes 0, 1 print from the simulation loop
  if((pars.printmode == 0) || (pars.printmode == 1)){
    out.open(outfilename);
  }
  double dt;                  //!< should initialise from pars, can be smaller
  double integrand;           //!< for ode advancing until discrete reaction
  double rnum;                //!< from unit exp dist
  double total_propensity;
  double total_population;
  int patch_choice;           //!< randomised with linear_sampler
  int population_choice;      //!< randomised with linear_sampler
  int reaction_choice;        //!< randomised with linear_sampler
  int ode_step_counter = 1;   //!< used for writing at pars.res steps
  int pm0counter = 0;         //!< hack counter for printmode 0 output :/
  int stop_code = -1;
  bool population_extinction; //!< returned from pop.minus()
  bool continuous_exist;      //!< for checking if ode systems are needed at all
  bool ode_stop_condition;
  //std::cout << "simulating..." << std::endl;
  bool cont = true; // TODO: rename this

  /// this is the main sim loop
  while(cont){
    
    /// check if there are continuous reactions
    continuous_exist = check_continuous();
    
    if(continuous_exist){
      
      /// numerically advance continuous states
      rnum = -log(unif(rng)); // unit exp dist
      integrand = 0.0;
      ode_stop_condition = false;
      while(!ode_stop_condition){
	dt = pars.dt;
	total_propensity = calculate_propensities();
	
	// check if integral is reached with smaller than dt step
	// now needs to be calculated twice per step...
	if((integrand + dt*total_propensity) < rnum){
	  integrand += dt*total_propensity;
	  advance_all_continuous(dt); // dt from par or smaller
	  time += dt;
	}else{ // TODO: step back ?
	  double difference = rnum - integrand;
	  double new_dt = difference/total_propensity;
	  integrand += new_dt*total_propensity;
	  advance_all_continuous(new_dt);
	  time += new_dt;
	  ode_stop_condition = true;
	}
	
	/// stop at time reached
	if(time > pars.max_time){
	  ode_stop_condition = true;
	}
	
	/// stop at a) target reached b) ntot < 1.0
	/// is this necessary?
	total_population = count_total_population();
	if(total_population > pars.target){
	  ode_stop_condition = true;
	  cont = false;
	  stop_code = 0;
	}else if(total_population < 1.0){
	  ode_stop_condition = true;
	  cont = false;
	  stop_code = 1;
	}

	/// write everything at pars.res steps if printmode == 0
	if((pars.printmode == 0) && (ode_step_counter % pars.res == 0)){
	  write_all(out,pm0counter);
	  pm0counter++;
	}
	
	/// write total population at pars.res steps if printmode == 1
	if((pars.printmode == 1) && (ode_step_counter % pars.res == 0)){
	  out << time << " " << num_events << " " << count_total_population() << std::endl;
	}

	ode_step_counter++;			     
      }
    }else{
      total_propensity = calculate_propensities();
      time += -log(unif(rng))/total_propensity;
    }
    patch_choice = linear_sampler(propensities);
    population_choice = linear_sampler(patches.at(patch_choice).propensities);
    reaction_choice = 
      linear_sampler(patches.at(patch_choice).populations.at(population_choice).propensities);
    if(!patches.at(patch_choice).active){
      std::cout << " error : !active patch_choice" << std::endl;
    }
    if(!patches.at(patch_choice).populations.at(population_choice).active){
      std::cout << " error : !active population_choice" << std::endl;
    }
    
    /// resolve discrete
    population_extinction = false;
    if(reaction_choice == 0){
      patches.at(patch_choice).populations.at(population_choice).plus();
    }else if(reaction_choice == 1){
      population_extinction = patches.at(patch_choice).populations.at(population_choice).minus();
    }else if(reaction_choice == 2){
      population_extinction = patches.at(patch_choice).populations.at(population_choice).minus();
      trait_set new_traits = trait_set_from_base(base_trait_set,true,false); // todo : only if not already muta
      add_patch(new_traits); // ,patch_choice,population_choice);
    }else if(reaction_choice == 3){
      trait_set new_traits = trait_set_from_base(base_trait_set,false,true); // todo : only if not already meta
      patches.at(patch_choice).add_population(new_traits,pars.threshold);
      num_mutations += 1;
    }else{
      std::cout << " error : resolver received inappropriate reaction_choice" << std::endl;
    }
    if(population_extinction){
      patches.at(patch_choice).propensities.at(population_choice) = 0.0; // this is the wrong way
      // should call simulation::population_deactivation from here
      if(patches.at(patch_choice).count_total_population() < 1){
	propensities.at(patch_choice) = 0.0;
	// should call simulation::patch_deactivation from here
	if(!patches.at(patch_choice).deactivate()){
	  std::cout << " error : patch.deactivate() returned false" << std::endl;
	}
      }
    }
    
    //print_propensities();
    num_events += 1;
    total_population = count_total_population();

    /// file output if printmode == 0
    if(pars.printmode == 0){
      write_all(out,pm0counter);
      pm0counter++;
    }
    
    /// file output if printmode == 1
    if(pars.printmode == 1){
      out << time << " " << num_events << " " << total_population << std::endl;;
    }
    
    /// stop conditions
    if(total_population < 1){
      cont = false;
      stop_code = 1;
      //std::cout << "stopping: extinction-condition " << std::endl;
    }
    if(num_events > pars.max_events){
      cont = false;
      stop_code = 2;
      //std::cout << "stopping: max-events-condition " << std::endl;
    }
    if(time > pars.max_time){
      cont = false;
      stop_code = 3;
      //std::cout << "stopping: max-time-condition " << std::endl;
    }
    if(total_population >= pars.target){
      cont = false;
      stop_code = 0;
      //std::cout << "stopping: target-n-condition " << std::endl;
    }
  }

  if((pars.printmode == 0) || (pars.printmode == 1)){
    out.close();
  }

  return stop_code;
}

// TODO: should thede be abs() ?
/*
trait_set simulation::mutate(const trait_set &old_ts)
{
  static std::normal_distribution<double> gaussian(0.0,1.0);
  trait_set new_ts;
  if(pars.mutationmode == 1){
    new_ts.b = old_ts.b + mutation_effects.b*abs(gaussian(rng));
    new_ts.d = old_ts.d - abs(mutation_effects.d*abs(gaussian(rng))); // stay positive :)
    new_ts.t = old_ts.t; // now should not mutate! perhaps in the future?
    new_ts.k = old_ts.k + mutation_effects.k*abs(gaussian(rng));
    new_ts.g = old_ts.g + mutation_effects.g*abs(gaussian(rng));
    new_ts.m = old_ts.m + mutation_effects.m*abs(gaussian(rng));
    new_ts.s = old_ts.s + mutation_effects.s*abs(gaussian(rng));
  }else if(pars.mutationmode == 2){
    new_ts.b = mutation_effects.b;
    new_ts.d = mutation_effects.d;
    new_ts.t = mutation_effects.t;
    new_ts.k = mutation_effects.k;
    new_ts.g = mutation_effects.g;
    new_ts.m = mutation_effects.m;
    new_ts.s = mutation_effects.s;
  }else{
    std::cout << " error : mutation mode not 1 or 2 " << std::endl;
  }
  return new_ts;
}
*/

void simulation::add_patch(trait_set ntraits) //, const int patch_choice, const int population_choice)
{
  patches.push_back(patch());
  propensities.push_back(0.0);
  num_active += 1;
  //trait_set traits_copy = patches.at(patch_choice).populations.at(population_choice).traits;
  patches.at(patches.size()-1).initialise(ntraits, pars.threshold, 1); // initialise the new patch
}

void simulation::advance_all_continuous(const double dt)
{
  for(auto &pt : patches){
    pt.advance_all_continuous(dt);
  }
}

double simulation::calculate_propensities()
{
  double total_propensity = 0.0;
  // only active patches
  for(size_t i = 0; i < patches.size(); ++i){
    if(patches.at(i).active){ // TODO: check that !active are 0.0
      propensities.at(i) = patches.at(i).patch_propensities();
      total_propensity += propensities.at(i);
    }
  }
  return total_propensity;
}

template<typename T>
int simulation::linear_sampler(const T& v)
{
  static std::uniform_real_distribution<double> unif(0.0,1.0);
  int m = 0;
  int imax = v.size();
  double cprob = v[0];
  double rmax = 0.0;
  for(auto j = 0; j < imax; j++){ // int
    rmax += v[j];
  }
  double rand = unif(rng);
  double target = rmax*rand;
  while( target > cprob ){
    m++;
    cprob += v[m]; // should try this: cprob += v[++m];
    if(m >= imax){
      std::cout << "sampler: m " << m << " imax " << imax << std::endl;
      m = imax - 1;
      break;
    }
  }
  return m;
}

// TODO: check that this works!
bool simulation::check_continuous() const
{
  bool any_cont = false;
  for(auto &pt : patches){
    if(pt.check_continuous()) any_cont = true;
  }
  return any_cont;
}

void simulation::print_propensities() const
{
  for(auto &pt : patches){
    pt.print_propensities();
  }
}

void simulation::print_patch_count() const
{
  std::cout << "num patches: " << patches.size() << std::endl;
}

double simulation::count_total_population() const
{
  double totpop = 0.0;
  for(auto &pt : patches){
    totpop += pt.count_total_population();
  }
  return totpop;
}

void simulation::print_state() const
{
  std::cout << "simulation end state: " << std::endl;
  std::cout << " num patches: " << patches.size() << " prop size " << propensities.size() << std::endl;
  for(size_t i = 0; i < patches.size(); ++i){
    std::cout << " patch " << i << ": " << std::endl;
    std::cout << " num populations: " << patches.at(i).populations.size() << std::endl;
    for(size_t j = 0; j < patches.at(i).populations.size(); ++j){
      std::cout << "  population " << j << " n " << patches.at(i).populations.at(j).read_n();
      std::cout << " b " << patches.at(i).populations.at(j).traits.b;
      std::cout << " d " << patches.at(i).populations.at(j).traits.d;
      std::cout << " t " << patches.at(i).populations.at(j).traits.t;
      std::cout << " k " << patches.at(i).populations.at(j).traits.k;
      std::cout << " g " << patches.at(i).populations.at(j).traits.g;
      std::cout << " m " << patches.at(i).populations.at(j).traits.m;
      std::cout << std::endl;
    }
  }
}

// rework this!
// printmode 0
void simulation::write_all(std::ofstream &o, const int id) const
{
  //static int counter = 0;
  // o << counter << " " << 0 << " " << time << " " << 0 << " " << 0 << " " << 0 << std::endl;
  // todo : modernise
  for(auto i = 0; i < patches.size(); i++){
    for(auto j = 0; j < patches.at(i).populations.size(); j++){
      o << id << " " << time << " " << i << " " << j << " ";
      o << patches.at(i).populations.at(j).read_n() << " ";
      o << patches.at(i).populations.at(j).traits.k;
      o << std::endl;
    }
  }
    
  //o << counter << " " << 1 << " " << 0 << " " << 0 << " " << patches.at(0).populations.at(0).read_n() << std::endl;
  //counter++;
}

// printmode 1
void simulation::write_summary(std::ofstream &o, const int id) const
{
  int cont_exist = 0;
  if(check_continuous()) cont_exist = 1;
  o << time << " " << cont_exist << " " << id << " " << num_events << " "
    << num_active << " "  << num_mutations << " " << count_total_population()
    << std::endl;
}

// printmode 2
void simulation::write_end_state(std::ofstream &out) const
{
  for(size_t i = 0; i < patches.size(); i++){
    for(size_t j = 0; j < patches.at(i).populations.size(); j++){
      out << index << " " << i << " " << j << " " << time << " ";
      out << patches.at(i).populations.at(j).read_n() << " ";
      out << patches.at(i).populations.at(j).traits.k;
      out << std::endl;
    }
  }
}

// printmode 3
void simulation::write_end_time(std::ofstream &out) const
{
  out << index << " " << time << std::endl;
}

