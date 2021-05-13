#include "header.h"

void Artic_sea::output_results(){
    //Print general information of the program
    if (config.master)
      cout  << config.ode_solver_type <<"\t"<< config.order<<"\t"<<config.refinements <<"\t"
	    << total_time <<"\t"<< total_error/vis_impressions <<"\t"<< iteration <<"\n";
}
