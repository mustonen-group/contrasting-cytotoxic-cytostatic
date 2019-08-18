#include <algorithm>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <memory>
#include <string>

#include "traits.h"     /// par and trait_set structs
#include "simulation.h" /// simulation object
#include "json/json.h"

/**
 *  used for '-f sim.json' sim specifications infile
 */
char* get_cmd_option(char **begin, char **end, const std::string &option);

/**
 * debugging only, not used currently
 */
void print_parameters(const par &p);

int main(int argc, char *argv[])
{
  char *filename = get_cmd_option(argv, argv + argc, "-f");
  std::string sim_json = filename;
  
  if(!filename){
    std::cout << "usage: ./sim -f sim.json" << std::endl;
  }else{
    
    std::cout << "reading main parameters from: " << filename << std::endl;
    par p; //!< holds sim specifications common to all simulations (+ seed values)
    bool sim_json_ok = p.read_sim_json(sim_json);
    
    if(!sim_json_ok){
      std::cout << " error : failed reading sim json" << std::endl;
    }else{
    
      // print_parameters(p);
    
      std::vector<std::unique_ptr<simulation>> sims; //!< container: individual simulations
      std::vector<int> stop_conditions;
      for(int i = 0; i < p.num_sims; i++){
	sims.push_back(std::make_unique<simulation>());
	stop_conditions.push_back(-1);
      }

      /// read trait jsons
      base_traits init_traits;
      bool read_traits_successful;
      //std::string trait_infile = p.traitfile;
      //std::string mutation_infile = p.mutationfile;
      // filenames should be in sim json
      // should check that file exists
      read_traits_successful = init_traits.read_traits_from_file(p.traitfile);
      
      if(!read_traits_successful){
	std::cout << " error : trait/mutation set reading unsuccessful" << std::endl;
      }else{

	// generate outfile names for printmodes 0 and 1
	// not necessary otherwise, but generated anyway
	std::vector<std::string> outfilenames;
	outfilenames.resize(p.num_sims);
	for(int i = 0; i < p.num_sims; i++){
	  std::stringstream out_filestream;
	  std::string num = std::to_string(i);
	  out_filestream << "out" << p.pre << num << ".txt";
	  outfilenames.at(i) = out_filestream.str();
	}

	if((p.printmode == 0) || (p.printmode == 1)){
	  //std::cout << "writing time-series to: " << std::endl;
	  for(auto &of : outfilenames) std::cout << of << std::endl;
	}
	
	// initialise simulations
	for(int i = 0; i < p.num_sims; i++){
	  sims.at(i)->initialise(init_traits,outfilenames.at(i),p,i);
	}

	// simulate and write
	// should restrict number of cores
        #pragma omp parallel for
	for(int i = 0; i < p.num_sims; i++){
	  stop_conditions.at(i) = sims.at(i)->simulate();
	}

	// see where we got
	/*
	for(auto &s : sims){
	  s->print_state();
	}
	*/
	
	// write endfile if printmode 2
	if((p.printmode == 2) || (p.printmode == 3)){
	  std::stringstream coll_filestream;
	  coll_filestream << "out" << p.pre << ".txt";
	  std::string coll_filename = coll_filestream.str();
	  std::cout << "printing end states to: " << coll_filename << std::endl;
	  std::ofstream coll;
	  coll.open(coll_filename);
	  if(p.printmode == 2){
	    // sim patch pop n b d t k g m s
	    for(auto &s : sims){
	      s->write_end_state(coll);
	    }
	  }else if(p.printmode == 3){
	    for(auto &s : sims){
	      s->write_end_time(coll);
	    }
	  }
	  coll.close();
	}

	// stop conditions to stdout
	for(auto i = 0; i < p.num_sims; i++){
	  std::cout << i << " " << stop_conditions.at(i) << std::endl;
	}
      }
    }
  }
  
  std::cout << "all done! exiting..." << std::endl;
  return 0;
}

char* get_cmd_option(char **begin, char **end, const std::string &option)
{
  char **itr = std::find(begin, end, option);
  if(itr != end && ++itr != end){
    return *itr;
  }else{
    return 0;
  }
}

void print_parameters(const par &p)
{
  std::cout << "num_sims:   " << p.num_sims   << std::endl;
  std::cout << "threshold:  " << p.threshold  << std::endl;
  std::cout << "init_n:     " << p.init_n     << std::endl;
  std::cout << "target:     " << p.target     << std::endl;
  std::cout << "max_time:   " << p.max_time   << std::endl;
  std::cout << "max_events: " << p.max_events << std::endl;
  std::cout << "res:        " << p.res        << std::endl;
  std::cout << "dt:         " << p.dt         << std::endl;
  std::cout << "printmode:  " << p.printmode  << std::endl;
  std::cout << "pre:        " << p.pre        << std::endl;
}

// c++17 tulossa:
/*
  std::vector<int> foo;
  foo.push_back(1);
  foo.push_back(2);
  std::for_each(
  		std::execution::par_unseq,
		foo.begin(),
		foo.end(),
		[](auto&& item)
		{
		  std::cout << item << std::endl;
		});		
*/
