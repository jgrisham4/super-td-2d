/*
 * This file is part of super-td-2d.
 * 
 * super-td-2d is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * super-td-2d is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with super-td-2d.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include "geom.h"
#include "mesh.h"
#include "solver.h"
#include "utils.h"

int main() {

  // Reading token-based input file
  input_data<double> data = read_input_file<double>("input.inp");

  // Setting up geometry
  geom<double> f(data.r,data.alpha);
  
  // Grid inputs
  double ymax_hat = 2.0/M_PI*atan(data.ymax/data.b);
  std::cout << "ymax_hat = " << ymax_hat << std::endl;
  
  // Creating grid in the computational plane
  mesh<double> grid(data.imax,data.jmax,data.a,data.b);
  grid.find_spacings(-1.0, 1.0, ymax_hat);
  generate_grid(grid);

  // Checking for restart
  if (data.restart) {

    // Instantiating objects
    soln<double> initial_guess;

    // Reading serialized object from files
    std::cout << "Reading restart data from " << data.restart_soln_file << std::endl;
    read_object<soln<double> >(initial_guess,data.restart_soln_file);

    // Making sure that the number of points is the same
    if (initial_guess.tau.size()!=data.imax*data.jmax) {
      std::cerr << "\nError: number of grid points in restart file: ";
      std::cerr << initial_guess.tau.size() << std::endl;
      std::cerr << "Number of points in mesh: " << data.imax*data.jmax << std::endl;
      std::cerr << "\nExiting.\n\n";
      exit(1);
    }

    // Solving interaction problem
    soln<double> s = solve(f,grid,data.dt,data.tol,data.max_iter,data.surface_file,initial_guess);
  
    // Writing results
    s.write(grid,data.tecplot_file);
    save_object(f, "geom.td");
    save_object(grid, "mesh.td");
    save_object(s, "soln.td");

  }
  else {

    // Solving interaction problem
    soln<double> s = solve(f,grid,data.dt,data.tol,data.max_iter,data.surface_file);
  
    // Writing results
    s.write(grid,data.tecplot_file);
    save_object(f, "geom.td");
    save_object(grid, "mesh.td");
    save_object(s, "soln.td");
    
  }

  return 0;

}
