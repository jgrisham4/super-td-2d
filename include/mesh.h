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

#ifndef MESH_HEADER_DEF
#define MESH_HEADER_DEF

#include <vector>
#include <boost/serialization/base_object.hpp>

/**
 * Simple container which holds mesh information.
 */
template <typename T>
struct mesh {

  // Below is required for serialization.
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive& ar, const unsigned int version) {
    ar & I & J & dx & dy & x1 & a & b & x & y;
  }

  // CTORs
  mesh() : I{0}, J{0} {};
  mesh(const unsigned int imax, const unsigned int jmax, const T A, const T B) : I{imax}, J{jmax}, a{A}, b{B} {};

  // Members
  T dx,dy,x1;
  T a,b;
  unsigned int I;
  unsigned int J;
  std::vector<T> x,y;

  // Methods
  void find_spacings(const T xmin, const T xmax, const T ymax);
};

/**
 * Method for finding dx and dy given I, J, along with the 
 * dimensions of the domain.  This method assumes that I and
 * J have already been set.
 */
template <typename T>
void mesh<T>::find_spacings(const T xmin, const T xmax, const T ymax) {
  x1 = xmin;
  dx = (xmax - xmin)/T(I-1);
  dy = ymax/T(J-1);
}

/**
 * Function for mesh generation.
 *
 * @param[in] m mesh object. 
 */
template <typename T>
void generate_grid(mesh<T>& m) {

  // Making sure that inputs were set in the mesh object
  if ((m.I==0)||(m.J==0)) {
    std::cerr << "\nError: inputs weren't set in the mesh object.\n\nExiting\n\n.";
    exit(1);
  }

  // Making sure that space for x and y has been allocated
  m.x.resize(m.I);
  m.y.resize(m.J);

  // Generating the mesh
  for (unsigned int i=0; i<m.I; ++i) {
    m.x[i] = T(i) * m.dx + m.x1;
  }
  for (unsigned int j=0; j<m.J; ++j) {
    m.y[j] = T(j) * m.dy;
  }

}

#endif
