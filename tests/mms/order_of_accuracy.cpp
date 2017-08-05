#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "geom.h"
#include "mesh.h"
#include "cassel_mms.h"
#include "utils.h"

void compute_soln(input_data<double>& data, const std::string& suffix) {

  // Setting up geometry
  geom<double> f(data.r,data.alpha);

  // Grid inputs
  double ymax_hat = 2.0/M_PI*atan(data.ymax/data.b);
  std::cout << "ymax_hat = " << ymax_hat << std::endl;

  // Creating grid in the computational plane
  mesh<double> grid(data.imax,data.jmax,data.a,data.b);
  grid.find_spacings(-1.0, 1.0, ymax_hat);
  generate_grid(grid);

  // Solving interaction problem
  soln<double> s = solve(f,grid,data.dt,data.tol,data.max_iter,data.surface_file);

  // Writing results
  s.write(grid,data.tecplot_file);
  save_object(f, std::string("geom")+suffix+std::string(".td"));
  save_object(grid, std::string("mesh")+suffix+std::string(".td"));
  save_object(s, std::string("soln")+suffix+std::string(".td"));

}

int main() {

  // Reading token-based input file
  input_data<double> data_coarse = read_input_file<double>("input_coarse.inp");
  input_data<double> data_medium = read_input_file<double>("input_medium.inp");
  input_data<double> data_fine   = read_input_file<double>("input_fine.inp");

  // Solving
  compute_soln(data_coarse, "_coarse");
  compute_soln(data_medium, "_medium");
  compute_soln(data_fine, "_fine");


  return 0;

}
