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

#ifndef SCALING_HEADER_DEF
#define SCALING_HEADER_DEF

#include <cmath>
#include <algorithm>
#include "geom.h"
#include "mesh.h"
#include "solver.h"

/**
 * The code in this header effectively reverses the triple-deck scalings.
 *
 * NOT COMPLETE.
 *
 * \author James Grisham
 * \date 06/24/2016
 */

/**
 * This function is for unscaling the triple-deck results so that they are in
 * dimensionless form.
 */
template <typename T>
void unscale(const T rho_w, const T mu_w, const T lambda, const T beta, const T Re_0, geom<T>& f, mesh<T>& grid, soln<T>& s) {

  // Computing actual ramp angle in radians
  T alpha_prime = sqrt(lambda*mu_w*beta)*pow(Re_0, -0.25)*f.alpha;

  // Filling in x and y vectors
  std::vector<T> xv(grid.I*grid.J),yv(grid.I*grid.J);
  for (int j=0; j<grid.J; ++j) {
    for (int i=0; i<grid.I; ++i) {
      xv[j*grid.I+i] = grid.x[i];
      yv[j*grid.I+i] = grid.y[j];
    }
  }

  // Overwriting old x- and y-vectors
  grid.x = xv;
  grid.y = yv;

  // Setting up lambda functions to be applied to vectors
  auto xhat_to_x = [&] (T xhat)   { return grid.a*atan(M_PI*xhat/2.0); };
  auto yhat_to_y = [&] (T yhat)   { return grid.b*atan(M_PI*yhat/2.0); };
  auto convert_x = [&] (T x)      { return pow(rho_w,-0.5)*pow(mu_w,-0.25)*pow(lambda,-1.25)*pow(beta,-0.75)*pow(Re_0,-3.0/8.0)*x; };
  auto convert_y = [&] (T x, T y) { return pow(rho_w,-0.5)*pow(mu_w,0.25)*pow(lambda,-0.75)*pow(beta,-0.25)*pow(Re_0,-5.0/8.0)*(y - f(x)); };
  auto convert_u = [&] (T u)      { return pow(rho_w,-0.5)*pow(mu_w,0.25)*pow(lambda,0.25)*pow(beta,-0.25)*pow(Re_0,-1.0/8.0)*u; };
  auto convert_v = [&] (T x, T v) { return pow(rho_w,-0.5)*pow(mu_w,0.75)*pow(lambda,0.75)*pow(beta,0.25)*pow(Re_0,-3.0/8.0)*(v - f.d(x)); };
  auto convert_p = [&] (T p)      { return sqrt(lambda*mu_w/beta)*pow(Re_0,-0.25)*p; };

  // Using std::transform to apply above lambda functions to vectors.
  std::transform(grid.x.begin(),grid.x.end(),grid.x.begin(),xhat_to_x);
  std::transform(grid.y.begin(),grid.y.end(),grid.y.begin(),yhat_to_y);
  std::transform(grid.x.begin(),grid.x.end(),grid.x.begin(),convert_x);
  std::transform(grid.x.begin(),grid.x.end(),grid.y.begin(),grid.y.begin(),convert_y);
  std::transform(s.u.begin(),s.u.end(),s.u.begin(),convert_u);
  std::transform(grid.x.begin(),grid.x.end(),s.v.begin(),s.v.begin(),convert_v);
  std::transform(s.p.begin(),s.p.end(),s.p.begin(),convert_p);

}

#endif
