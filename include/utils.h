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

#ifndef UTILS_HEADER_DEF
#define UTILS_HEADER_DEF

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>

/**
 * This header contains utility functions for reading data from
 * an input file, and serialization.
 */

template <typename T>
struct input_data {

  std::string input_file;          ///< Input file from which data is read.
  std::string surface_file;        ///< Output file which contains data at the wall.
  std::string tecplot_file;        ///< Output Tecplot file.
  std::string restart_soln_file;   ///< Solution restart file.
  T tol;                           ///< Relative tolerance used for stopping.
  T dt;                            ///< Time step.
  T alpha;                         ///< Scaled ramp angle.
  T r;                             ///< Corner smoothing factor.
  T a;                             ///< Parameter which control grid clustering.
  T b;                             ///< Parameter which control grid clustering.
  T ymax;                          ///< Used to set y-extent of domain.
  int imax;                        ///< Number of points in the x-direction.
  int jmax;                        ///< Number of points in the y-direction.
  int max_iter;                    ///< Max number of time steps allowed.
  bool restart;                    ///< Determines whether to restart or not.

};

/**
 * Function used to read data from input file.
 *
 * @param[in] input_file std::string which contains input file name.
 * @return input_data struct which contains input data.
 */
template <typename T>
input_data<T> read_input_file(const std::string& input_file) {

  input_data<T> tmp;
  tmp.restart = false;
  std::string token,line;
  std::istringstream iss;
  std::ifstream infile(input_file.c_str());
  bool found_soln_restart = false;
  if (!infile.is_open()) {
    std::cerr << "\nError: Can't open input.inp.\nExiting.\n\n";
    exit(1);
  }
  while (std::getline(infile,line)) {

    // Making sure current line isn't a comment
    if (line.find("#")==std::string::npos) {

      if (line.find("tol")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.tol;
        iss.clear();
        std::cout << "Tolerance: " << tmp.tol << std::endl;
      }
      if (line.find("dt")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.dt;
        iss.clear();
        std::cout << "Time step: " << tmp.dt << std::endl;
      }
      if (line.find("alpha")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.alpha;
        iss.clear();
        std::cout << "Scaled ramp angle: " << tmp.alpha << std::endl;
      }
      if (line.find("R")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.r;
        iss.clear();
        std::cout << "Smoothing factor: " << tmp.r << std::endl;
      }
      if (line.find("imax")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.imax;
        iss.clear();
        std::cout << "imax: " << tmp.imax << std::endl;
      }
      if (line.find("jmax")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.jmax;
        iss.clear();
        std::cout << "jmax: " << tmp.jmax << std::endl;
      }
      if (line.find("A")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.a;
        iss.clear();
        std::cout << "a: " << tmp.a << std::endl;
      }
      if (line.find("B")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.b;
        iss.clear();
        std::cout << "b: " << tmp.b << std::endl;
      }
      if (line.find("ymax")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.ymax;
        iss.clear();
        std::cout << "ymax: " << tmp.ymax << std::endl;
      }
      if (line.find("max_iter")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.max_iter;
        iss.clear();
        std::cout << "Max iterations: " << tmp.max_iter << std::endl;
      }
      if (line.find("surface_file")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.surface_file;
        iss.clear();
        std::cout << "Surface output file: " << tmp.surface_file << std::endl;
      }
      if (line.find("tecplot_file")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.tecplot_file;
        iss.clear();
        std::cout << "Tecplot output file: " << tmp.tecplot_file << std::endl;
      }
      if (line.find("restart_soln")!=std::string::npos) {
        iss.str(line);
        iss >> token >> tmp.restart_soln_file;
        iss.clear();
        std::cout << "Solution restart file: " << tmp.restart_soln_file << std::endl;
        found_soln_restart = true;
      }
    }
  }
  infile.close();

  // Checking if restart is necessary
  if (found_soln_restart) {
    tmp.restart = true;
  }

  return tmp;

}

/**
 * Function for writing object to file.
 *
 * @param[in] obj object of generic type which will be serialized.
 * @param[in] file_name std::string which contains file name.
 */
template <typename A>
void save_object(const A& obj, const std::string file_name) {
  std::ofstream ofs(file_name.c_str());
  boost::archive::text_oarchive oa(ofs);
  oa << obj;
}

/**
 * Function for reading an object from file.
 *
 * @param[in] obj object of generic type which will be deserialized.
 * @param[in] file_name std::string which contains file name.
 */
template <typename A>
void read_object(A& obj, const std::string file_name) {
  std::ifstream ifs(file_name.c_str());
  boost::archive::text_iarchive ia(ifs);
  ia >> obj;
}


#endif
