#include "traits.h"

// TODO: should check that all fields exist
bool par::read_sim_json(const std::string sim_json_file)
{
  std::ifstream ifs(sim_json_file);
  std::string sim_str((std::istreambuf_iterator<char>(ifs)),
		       std::istreambuf_iterator<char>());
  Json::CharReaderBuilder builder;
  Json::CharReader *reader = builder.newCharReader(); // remember to delete
  Json::Value obj;
  std::string errors;
  bool parsing_successful = reader->parse(sim_str.c_str(),
                                          sim_str.c_str() + sim_str.size(),
                                          &obj, &errors);
  num_sims   = obj["num_sims"].asInt();
  threshold  = obj["threshold"].asInt();
  init_n     = obj["init_n"].asDouble();
  target     = obj["target"].asDouble();
  max_time   = obj["max_time"].asDouble();
  max_events = obj["max_events"].asInt();
  res        = obj["res"].asInt();
  dt         = obj["dt"].asDouble();
  printmode  = obj["printmode"].asInt();
  pre        = obj["pre"].asString();
  traitfile  = obj["traitfile"].asString();
  seeds.resize(num_sims);
  for(int i = 0; i < num_sims; ++i){
    seeds.at(i) = obj["seeds"][i].asInt();
  }
  delete reader;
  if(!parsing_successful){
    std::cout << " error : parsing sim json failed" << std::endl;
    std::cout << errors << std::endl;
    return false;
  }else{
    return true;
  }
}

// TODO: should check that all fields exist
bool base_traits::read_traits_from_file(const std::string trait_file)
{
  std::ifstream ifs(trait_file);
  std::string traits_str((std::istreambuf_iterator<char>(ifs)),
			  std::istreambuf_iterator<char>());
  Json::CharReaderBuilder builder;
  Json::CharReader *reader = builder.newCharReader(); // remember to delete
  Json::Value obj;
  std::string errors;
  bool parsing_successful = reader->parse(traits_str.c_str(),
					  traits_str.c_str() + traits_str.size(),
					  &obj, &errors);
  b_primary    = obj["b_primary"].asDouble();
  b_metastatic = obj["b_metastatic"].asDouble();
  b_mutant     = obj["b_mutant"].asDouble();
  d_primary    = obj["d_primary"].asDouble();
  d_metastatic = obj["d_metastatic"].asDouble();
  d_mutant     = obj["d_mutant"].asDouble();
  t = obj["t"].asDouble();
  k_primary    = obj["k_primary"].asDouble();
  k_metastatic = obj["k_metastatic"].asDouble();
  k_mutant     = obj["k_mutant"].asDouble();
  m = obj["m"].asDouble();
  g = obj["g"].asDouble();
  delete reader;
  if(!parsing_successful){
    std::cout << " error : parsing traits json failed" << std::endl;
    return false;
  }else{
    return true;
  }
}

trait_set trait_set_from_base(const base_traits& base, const bool meta, const bool muta)
{
  trait_set ts;
  ts.t = base.t;
  ts.m = base.m;
  ts.g = base.g;
  if(muta){
    ts.b = base.b_mutant;
    ts.d = base.d_mutant;
    ts.k = base.k_mutant;
  }else if(meta){
    ts.b = base.b_metastatic;
    ts.d = base.d_metastatic;
    ts.k = base.k_metastatic;
  }else{
    ts.b = base.b_primary;
    ts.d = base.d_primary;
    ts.k = base.k_primary;
  }
  return ts;
}
