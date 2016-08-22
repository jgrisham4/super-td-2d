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

#ifndef GEOM_HEADER_DEF
#define GEOM_HEADER_DEF

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/serialization/base_object.hpp>
#include "viscosity.h"

#define LAMBDA 0.332

/**
 * Functor for geometry.
 */
template <typename T>
class geom {

  // Below is required for serialization.
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive& ar, const unsigned int version) {
    ar &r;
    ar &alpha;
  }


  public:
    geom() {};
    geom(const T R, const T A) : r{R}, alpha{A} {};
    T operator() (const T x) const;
    T d(const T x) const;
    T d2(const T x) const;
    T r;
    T alpha;
};

/**
 * Defining the surface shape.
 */
template <typename T>
T geom<T>::operator() (const T x) const {
  return 0.5*alpha*(x + sqrt(x*x + r*r));
}

/**
 * First derivative of the shape.
 */
template <typename T>
T geom<T>::d(const T x) const {
  return 0.5*alpha*(1.0 + x/sqrt(r*r + x*x));
}

/**
 * Second derivative of the shape.
 */
template <typename T>
T geom<T>::d2(const T x) const {
  return r*r*alpha/(2.0*pow(r*r + x*x,3.0/2.0));
}

/**
 * Function for writing surface to file.
 */
template <typename T>
void write_geometry(const T xmin, const T xmax, const unsigned int npts, const std::string& fname, const geom<T>& g) {

  // Opening stream
  std::ofstream outfile(fname.c_str());
  if (!outfile.is_open()) {
    std::cerr << "\nCan't open " << fname << " for writing." << std::endl;
    std::cerr << "Exiting." << std::endl;
    exit(-1);
  }

  // Constructing x-vector
  std::vector<T> x(npts);
  T dx = (xmax - xmin)/(T(npts)-1.0);
  for (unsigned int i=0; i<npts; ++i) {
    x[i] = T(i)*dx + xmin;
  }

  // Writing data to file
  outfile.precision(10);
  outfile.setf(std::ios_base::scientific);
  for (unsigned int i=0; i<npts; ++i) {
    outfile << x[i] << " " << g(x[i]) << "\n";
  }

  // Closing file stream
  outfile.close();

}

/**
 * Alternative functor for my geometry.
 */
template <typename T>
class circ_ramp {

  public:

    circ_ramp(const T radius, const T alpha_p, const T Me, const T ReL, const T T0, const T L, const T p0);
    T operator() (const T x) const;
    T d(const T x) const;
    T d2(const T x) const;
    T r;
    T alpha;
    T x1,x2,y1,y2;
};

/**
 * CTOR.
 * 
 * @param[in] radius radius of the circular arc in dimensional form.
 * @param[in] alpha_p deflection angle of the ramp (radians).
 * @param[in] Me Mach number at the boundary layer edge.
 * @param[in] ReL Reynolds number based on some reference length L.
 * @param[in] Te Temperature at the boundary layer edge.
 * @param[in] L reference length (dimensional).
 * @param[in] pe pressure at the edge of the boundary layer (Pa).
 */
template <typename T>
circ_ramp<T>::circ_ramp(const T radius, const T alpha_p, const T Me, const T ReL, const T T0, const T L, const T p0) {

  // Declaring variables
  std::cout << "\nConstructing scaled ramp geometry...\n" << std::endl;
  T Taw, muw, mu0, recovery_factor, Ue;
  T gamma = 1.4;
  T R = 287.0;
  T rhow,Re0,Te,pe;

  // Defining recovery factor
  recovery_factor = sqrt(0.715);

  // Computing Te and pe
  Te = T0/(1.0 + (gamma-1.0)/2.0*Me*Me);
  pe = p0/pow(T0/Te,gamma/(gamma-1.0));
  std::cout << "Edge temperature: " << Te << std::endl;
  std::cout << "Edge pressure: " << pe << std::endl;

  // Finding adiabatic wall temperature
  Taw = Te*(1.0 + recovery_factor*(gamma -1.0)/2.0*Me*Me);
  rhow = pe/(R*Taw);  // because p = const in the normal direction and Tw = Taw

  // Finding muw (assuming adiabatic wall so that Tw = Taw)
  muw = mu(Taw);
  std::cout << "muw = " << muw << std::endl;

  // Finding mu0 (assuming CPG again)
  Ue = Me*sqrt(gamma*R*Te);
  mu0 = mu(Ue*Ue*(gamma-1.0)/(gamma*R));
  std::cout << "mu0 = " << mu0 << std::endl;

  // Finding Re based on mu_0
  Re0 = ReL*mu(Te)/mu0;
  std::cout << "Re0 = " << Re0 << std::endl;

  // Nondimensionalizing muw and rhow by the freestream values
  rhow /= pe/(R*Te);
  //muw /= mu(Te);  // This works for some reason, but it is wrong.
  muw /= mu0;
  std::cout << "Scaled wall viscosity: " << muw << std::endl;


  // Finding the scaled ramp angle
  T beta = sqrt(Me*Me - 1.0);
  alpha = alpha_p/(sqrt(LAMBDA*muw*beta)*pow(Re0, -0.25));
  std::cout << "Scaled ramp angle: " << alpha << std::endl;

  // Defining lambda functions for use in scaling x- and y-coordinates
  // x and y values should be normalized by L before calling this function.
  auto scale_x = [&] (const T x_norm) { return (x_norm - 1.0)/(pow(rhow,-0.5)*pow(muw,-0.25)*pow(LAMBDA,-1.25)*pow(beta,-0.75)*pow(Re0,-3.0/8.0)); };
  auto scale_y = [&] (const T x_scaled, const T y_norm) { return y_norm/(pow(rhow,-0.5)*pow(muw,0.25)*pow(LAMBDA,-0.75)*pow(beta,-0.25)*pow(Re0,-5.0/8.0)) + (*this)(x_scaled); };

  // vvvvvv NOT SURE ABOUT THE BELOW vvvvvv

  // Computing some geometric quantities 
  x1 = scale_x((radius*(1.0 - cos(alpha_p))/tan(alpha_p) - radius*sin(alpha_p))/L);
  y1 = 0.0;
  x2 = scale_x(radius*(1.0 - cos(alpha_p))/tan(alpha_p)/L);
  y2 = scale_y(x2,radius*(1.0 - cos(alpha_p))/L);

  // ********* Computing new r in scaled space *********
  r = - (x1 - x2)*sqrt(1.0 - alpha*alpha)/alpha;
  std::cout << "Physical radius: " << radius << std::endl;
  std::cout << "Scaled radius: " << r << "\n\n";

}

/**
 * Defining the surface shape.  This shape is defined in terms of the scaled
 * x- and y-values
 */
template <typename T>
T circ_ramp<T>::operator() (const T x) const {
  auto f = [&] (T x) { };
  return x < x1 ? 0.0 : x >= x1 && x <= x2 ? sqrt(r*r-pow(x-x1,2))+r : alpha*x;
}

/**
 * First derivative of the shape.
 */
template <typename T>
T circ_ramp<T>::d(const T x) const {
  return x < x1 ? 0.0 : x >= x1 && x <= x2 ? (x-x1)/sqrt(r*r-pow(x-x1,2)) : alpha;
}

/**
 * Second derivative of the shape.
 */
template <typename T>
T circ_ramp<T>::d2(const T x) const {
  return x >= x1 && x <= x2 ? 1.0/sqrt(r*r-pow(x-x1,2))+(x-x1)*(x-x1)/pow(r*r-pow(x-x1,2),3.0/2.0) : 0.0;
}


#endif
