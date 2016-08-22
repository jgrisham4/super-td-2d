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

#ifndef LINALG_HEADER_DEF
#define LINALG_HEADER_DEF

#include <vector>

/**
 * Function for solving tridiagonal linear system using Thomas algorithm.
 * 
 * Meant for tridiagonal problems of the form
 *
 * b(j) u(j-1) + d(j) u(j) + a(j) u(j+1) = c(j)
 *
 */
template <typename T>
std::vector<T> thomas(std::vector<T>& b, std::vector<T>& d, std::vector<T>& a, std::vector<T>& c) {

  // Making sure that vectors have been initialized
  int nj = c.size();
  if (nj==0) {
    std::cerr << "\nError: a.size() == 0.\n\nExiting.\n\n";
    exit(1);
  }
  std::vector<T> x(nj);

  // Converting tridiagonal system to upper triangular
  for (int j=1; j<nj; ++j) {
    d[j] = d[j] - b[j]/d[j-1]*a[j-1];
    c[j] = c[j] - b[j]/d[j-1]*c[j-1];
  }

  // Using back substitution to find the solution
  x[nj-1] = c[nj-1]/d[nj-1];
  for (int j=nj-2; j>=0; --j) {
    x[j] = (c[j]-a[j]*x[j+1])/d[j];
  }

  return x;

}

#endif
