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

#ifndef CASSEL_HEADER_DEF
#define CASSEL_HEADER_DEF

#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/serialization/base_object.hpp>
#include "integration.h"
#include "linalg.h"


/**
 * Function used in coordinate transformation.
 */
template <typename T>
constexpr T Gamma(const T zeta) {
  return 1.0/M_PI*(1.0 + cos(M_PI*zeta));
}

/**
 * Derivative of the above.
 */
template <typename T>
constexpr T dGamma(const T zeta) {
  return -sin(M_PI*zeta);
}

/**
 * Simple container for the solution.
 */
template <typename T>
struct soln {

  // Below is required for serialization.
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive& ar, const unsigned int version) {
    ar & I & J & tau & u & v & psi & p;
  }

  // Declaring variables
  int I,J;
  std::vector<T> tau, u, v, psi, p;

  // CTORs
  soln() : I{0}, J{0} {};
  soln(const int imax, const int jmax);

  // Methods
  void find_u   (const mesh<T>& grid);
  void find_v   (const mesh<T>& grid);
  void find_psi (const mesh<T>& grid);
  void write    (const mesh<T>& grid, const std::string& fname) const;

};

/**
 * Function for finding the x-component of velocity when the
 * shear stress, tau, is known.
 */
template <typename T>
void soln<T>::find_u(const mesh<T>& grid) {

  // Declaring variables
  T sum, dy, ga, gb;
  dy = grid.dy;

  // Making sure that u_{i,0} = 0
  for (int i=0; i<I; ++i) {
    u[i] = 0.0;
    psi[i] = 0.0;
  }

  // Integrating
  for (int i=1; i<I; ++i) {

    sum = T(0.0);
    for (int j=1; j<J; ++j) {

      ga = Gamma<T>(grid.y[j]);
      gb = Gamma<T>(grid.y[j-1]);
      sum += grid.b*integrate_trap_disc<T>(dy,tau[j*I+i]/ga,tau[(j-1)*I+i]/gb);
      u[j*I+i] = sum;

    }

  }

}

/**
 * Function for finding the y-component of velocity.  The stream
 * function, psi, must be determined prior to calling this function.
 * This function uses the fact that v = -d/dx ( psi ).
 */
template <typename T>
void soln<T>::find_v(const mesh<T>& grid) {

  // Declaring variables
  T dx = grid.dx;
  T val1,valI;

  // Assigning v on the bottom
  for (int i=0; i<I; ++i) {
    v[i] = 0.0;
  }

  // Computing v at the inlet and exit
  val1 = -Gamma<T>(grid.x[0])/grid.a;
  valI = -Gamma<T>(grid.x[I-1])/grid.a;
  for (int j=1; j<J; ++j) {
    v[j*I] = val1*(-3.0*psi[j*I] + 4.0*psi[j*I+1] - psi[j*I+2])/(2.0*dx);
    v[j*I+I-1] = valI*(3.0*psi[j*I+I-1] - 4.0*psi[j*I+I-2] + psi[j*I+I-3])/(2.0*dx);
  }

  // Computing v on the interior
  for (int j=1; j<J; ++j) {
    for (int i=1; i<I-1; ++i) {
      v[j*I+i] = -Gamma<T>(grid.x[i])/grid.a*(psi[j*I+i+1] - psi[j*I+i-1])/(2.0*dx);
    }
  }

}

/**
 * Function for computing the stream function by integration.  The
 * x-component of velocity must already be known.
 */
template <typename T>
void soln<T>::find_psi(const mesh<T>& grid) {

  // Declaring variables
  T sum, dy, ga, gb;
  dy = grid.dy;

  // Integrating
  for (int i=0; i<I; ++i) {

    sum = T(0.0);
    for (int j=1; j<J; ++j) {

      ga = Gamma<T>(grid.y[j]);
      gb = Gamma<T>(grid.y[j-1]);
      sum += integrate_trap_disc<T>(dy,grid.b*u[j*I+i]/ga,grid.b*u[(j-1)*I+i]/gb);
      psi[j*I+i] = sum;

    }

  }

}

/**
 * CTOR.
 */
template <typename T>
soln<T>::soln(const int imax, const int jmax) {
  I = imax;
  J = jmax;
  tau.resize(I*J);
  u.resize(I*J);
  v.resize(I*J);
  psi.resize(I*J);
  p.resize(I);
}

/**
 * Method for writing data to file.
 */
template <typename T>
void soln<T>::write(const mesh<T>& grid, const std::string& fname) const {

  std::cout << "\nWriting solution to file named " << fname << std::endl;

  std::ofstream tecfile(fname.c_str());
  if (!tecfile.is_open()) {
    std::cerr << "Error: Can't write data to file named " << fname << "\n";
    std::cerr << "Exiting.\n\n";
    exit(-1);
  }
  tecfile.precision(16);
  tecfile.setf(std::ios_base::scientific);
  tecfile << "variables=\"x\",\"y\",\"tau\",\"u\",\"v\",\"psi\"\n";
  tecfile << "zone i=" << I << ", j=" << J << ", f=point\n";
  for (int j=0; j<J; ++j) {
    for (int i=0; i<I; ++i) {
      tecfile << grid.x[i] << " " << grid.y[j] << " " << tau[j*I+i] << " " << u[j*I+i]
        << " " << v[j*I+i] << " " << psi[j*I+i] << "\n";
    }
  }
  tecfile.close();

}

/**
 * Function used to find L2 norm between two solutions.  This
 * function actually finds the norm between wall shear stresses.
 */
template <typename T>
T norm(const soln<T>& s1, const soln<T>& s2) {

  // Declaring variables
  int imax = s1.I;
  T sum = T(0.0);
  T value;

  // Computing norm
  for (int i=0; i<imax; ++i) {
    value = s1.tau[i] - s2.tau[i];
    sum += value*value;
  }

  return sqrt(sum);

}

/**
 * Function for computing the convective term.  If reversed flow exists,
 * a downstream stencil is used instead of an upstream stencil.
 *
 * @param[in] s solution object from previous timestep.
 */
template <typename T>
std::vector<T> compute_convective_term(const mesh<T>& grid, soln<T>& s) {

  // Declaring variables
  T dx = grid.dx;
  std::vector<T> udtaudx(s.I*s.J);

  // Assigning conditions although these values should never be used.
  // This can be inferred by looking at the difference equations.
  for (int i=0; i<s.I; ++i) {
    udtaudx[i] = 0.0;
    udtaudx[(s.J-1)*s.I+i] = 0.0;
  }
  for (int j=0; j<s.J; ++j) {
    udtaudx[j*s.I+0] = 0.0;
    udtaudx[j*s.I+s.I-1] = 0.0;
  }

  // Computing u*d(tau)/dx
  for (int i=1; i<s.I-1; ++i) {
    for (int j=1; j<s.J; ++j) {

      // Checking for reversed flow. Also using knowledge of the boundary
      // conditions to avoid having to use ghost nodes outside the domain.
      // Specifically, tau -> 1 on the inlet and outlet.
      if (s.u[j*s.I+i]>=T(0.0)) {
        if (i!=1) {
          udtaudx[j*s.I+i] = s.u[j*s.I+i]*(3.0*s.tau[j*s.I+i]-4.0*s.tau[j*s.I+i-1]+s.tau[j*s.I+i-2])/(2.0*dx);
        }
        else {
          udtaudx[j*s.I+i] = s.u[j*s.I+i]*(3.0*s.tau[j*s.I+i]-4.0*s.tau[j*s.I+i-1]+T(1.0))/(2.0*dx);
        }
      }
      else {
        if (i!=(s.I-1)) {
          udtaudx[j*s.I+i] = -s.u[j*s.I+i]*(3.0*s.tau[j*s.I+i]-4.0*s.tau[j*s.I+i+1]+s.tau[j*s.I+i+2])/(2.0*dx);
        }
        else {
          udtaudx[j*s.I+i] = -s.u[j*s.I+i]*(3.0*s.tau[j*s.I+i]-4.0*s.tau[j*s.I+i+1]+T(1.0))/(2.0*dx);
        }
      }
    }
  }

  return udtaudx;
}

/**
 * Function for advancing one step in time.
 */
template <typename T, typename G>
void update(const G& surf, const mesh<T>& grid, const T dt, soln<T>& s_n, soln<T>& s_np1, const bool compute_pressure=false) {

  // Declaring variables
  static int call_num = 0;
  int imax,jmax;
  T dx,dy,a,b;
  std::vector<T> R,Q,C,B;

  // Getting grid info
  imax = grid.I;
  jmax = grid.J;
  dx = grid.dx;
  dy = grid.dy;
  a = grid.a;
  b = grid.b;

  // Allocating space for vectors
  R.resize(imax*jmax);
  Q.resize(imax*jmax);
  C.resize(imax*jmax);
  B.resize(imax*jmax);

  // Indices used in accumulate below
  std::vector<int> jind(jmax-1);
  std::vector<T> Nv(jind.size());
  std::vector<T> Mv(jind.size());
  for (int j=1; j<jmax; ++j) {
    jind[j-1] = j;
  }

  // Applying boundary conditions
  for (int i=0; i<imax; ++i) {
    s_np1.u[i] = T(0);
    s_np1.v[i] = T(0);
    s_np1.psi[i] = T(0);
    s_np1.tau[(jmax-1)*imax+i] = T(1);
  }
  for (int j=0; j<jmax; ++j) {
    s_np1.tau[j*imax] = T(1);
    s_np1.tau[j*imax+imax-1] = T(1);
    s_np1.u[j*imax] = grid.y[j];     // Not sure if this should be yhat or y.
  }

  // Computing convective term from previous time step
  std::vector<T> udtaudx = compute_convective_term(grid,s_n);

  // Computing A once and for all
  static std::vector<T> A;
  A.resize(jmax);
  if (call_num==0) {
    for (int j=0; j<jmax; ++j) {
      A[j] =  Gamma<T>(grid.y[j])*Gamma<T>(grid.y[j])/((b*dy)*(b*dy));
    }
  }

  // Lambda functions used in finite difference scheme
  auto c          = [&] (const int i, const int j) { return -2.0*A[j]-1.0/dt; };
  auto c_minus    = [&] (const int i, const int j) { return A[j] - Gamma<T>(grid.y[j])/(2.0*b*dy)*(dGamma<T>(grid.y[j])/b - s_n.v[j*imax+i]); };
  auto c_plus     = [&] (const int i, const int j) { return 2.0*A[j] - c_minus(i,j); };
  auto d          = [&] (const int i, const int j) { return Gamma<T>(grid.x[i])/a*udtaudx[j*imax+i]-s_n.tau[j*imax+i]/dt; };

  // Filling in vectors rather than using recursive lambdas.
  for (int i=1; i<imax-1; ++i) {
    R[(jmax-1)*imax+i] = T(0);
    Q[(jmax-1)*imax+i] = T(1);
    C[i] = T(1);
    B[i] = T(0);
    for (int j=jmax-2; j>0; --j) {
      R[j*imax+i] = -c_minus(i,j)/(c(i,j)+c_plus(i,j)*R[(j+1)*imax+i]);
      Q[j*imax+i] = -(c_plus(i,j)*Q[(j+1)*imax+i] - d(i,j))/(c(i,j) + c_plus(i,j)*R[(j+1)*imax+i]);
    }
    for (int j=1; j<jmax; ++j) {
      C[j*imax+i] = R[j*imax+i]*C[(j-1)*imax+i];
      B[j*imax+i] = R[j*imax+i]*B[(j-1)*imax+i] + Q[j*imax+i];
    }
  }

  // These are the terms inside the sums in (4.13a,b).
  auto Nterm      = [&] (const int i, const int j) { return dy/2.0*(b/Gamma<T>(grid.y[j])*C[j*imax+i] + b/Gamma<T>(grid.y[j-1])*C[(j-1)*imax+i]); };
  auto Mterm      = [&] (const int i, const int j) { return dy/2.0*(b/Gamma<T>(grid.y[j])*(B[j*imax+i]-1.0)+b/Gamma<T>(grid.y[j-1])*(B[(j-1)*imax+i]-1.0)); };

  // These are the sums from (4.13a,b) using functional concepts
  auto N          = [&] (const int i) {
    auto dummy1 = [&] (const int jv) { return Nterm(i,jv); };   // Binding argument of Nterm
    std::transform(jind.begin(),jind.end(),Nv.begin(),dummy1);
    return std::accumulate(Nv.begin(),Nv.end(),T(0));
  };
  auto M          = [&] (const int i) {
    auto dummy2 = [&] (const int jv) { return Mterm(i,jv); };
    std::transform(jind.begin(),jind.end(),Mv.begin(),dummy2);
    return std::accumulate(Mv.begin(),Mv.end(),T(0));
  };

  // Coefficients in tridiagonal problem given in (4.15)
  auto cbar       = [&] (const int i) { return pow(Gamma<T>(grid.x[i])/a,2)*2.0*N(i)/(dx*dx)-Gamma<T>(0)/b*(4.0*C[1*imax+i]-3.0-C[2*imax+i])/(2.0*dy); };
  auto cbar_minus = [&] (const int i) { return -pow(Gamma<T>(grid.x[i])/a,2)*N(i-1)/(dx*dx)+Gamma<T>(grid.x[i])*dGamma<T>(grid.x[i])/(a*a)*N(i-1)/(2.0*dx); };
  auto cbar_plus  = [&] (const int i) { return -pow(Gamma<T>(grid.x[i])/a,2)*N(i+1)/(dx*dx)-Gamma<T>(grid.x[i])*dGamma<T>(grid.x[i])/(a*a)*N(i+1)/(2.0*dx); };
  auto dbar       = [&] (const int i) { return pow(Gamma<T>(grid.x[i])/a,2)*(M(i+1)-2.0*M(i)+M(i-1))/(dx*dx)+Gamma<T>(grid.x[i])*dGamma<T>(grid.x[i])/(a*a)*(M(i+1) - M(i-1))/(2.0*dx) - surf.d2(a*tan(grid.x[i]*M_PI/2.0)) + Gamma<T>(0)/b*(4.0*B[1*imax+i] - B[2*imax+i])/(2.0*dy); };

  // Finding solution for tau_w at current time step
  std::vector<T> tau_w(imax), a_w(imax), b_w(imax), c_w(imax), d_w(imax);
  d_w[0] = T(1);
  d_w[imax-1] = T(1);
  c_w[0] = T(0);
  b_w[0] = T(1);
  b_w[imax-1] = T(1);
  a_w[imax-1] = T(0);
  for (int i=1; i<imax-1; ++i) {
    a_w[i] = cbar_minus(i);
    b_w[i] = cbar(i);
    c_w[i] = cbar_plus(i);
    d_w[i] = dbar(i);
  }
  tau_w = thomas<T>(a_w,b_w,c_w,d_w);
  for (int i=0; i<imax; ++i) {
    s_np1.tau[i] = tau_w[i];
  }

  // Solving for tau
  for (int i=1; i<imax-1; ++i) {
    for (int j=0; j<jmax; ++j) {
      s_np1.tau[j*imax+i] = C[j*imax+i]*s_np1.tau[i] + B[j*imax+i];
    }
  }

  // Finding x-component of velocity
  s_np1.find_u(grid);

  // Finding stream function
  s_np1.find_psi(grid);

  // Finding y-component of velocity
  s_np1.find_v(grid);

  // Computing pressure
  if (compute_pressure) {
    std::cout << "Computing pressure..." << std::endl;
    s_np1.p[0] = T(0.0);
    s_np1.p[imax-1] = T(surf.alpha);
    for (int i=1; i<imax-1; ++i) {
      s_np1.p[i] = surf.d(a*tan(M_PI*grid.x[i]/2.0))-Gamma<T>(grid.x[i])/a*((N(i+1)*s_np1.tau[i] + M(i+1) - N(i-1)*s_np1.tau[i-1] - M(i-1))/(2.0*dx));
    }

  }

  ++call_num;

}

/**
 * Function for advancing until the steady solution has been achieved.
 */
template <typename T, typename G>
soln<T> solve(const G& f, const mesh<T>& grid, const T dt, const T tol, const int max_iter, const std::string& surf_file, const soln<T>& initial_guess=soln<T>()) {

  // Creating solution objects
  soln<T> s_n(grid.I, grid.J), s_np1(grid.I, grid.J);
  std::vector<T> residuals;
  int imax, jmax, iter;
  imax = grid.I;
  jmax = grid.J;
  T eps = 100.0;
  T t = 0.0;
  T initial_norm = 1.0;
  iter = 0;

  // Checking for restart
  if ((initial_guess.I==0)&&(initial_guess.J==0)) {

    // Assigning initial guess
    for (int i=0; i<s_n.I; ++i) {
      for (int j=0; j<s_n.J; ++j) {
        s_n.tau[j*imax+i] = 1.0;
      }
    }
    s_n.find_u(grid);
    s_n.find_psi(grid);
    s_n.find_v(grid);

  }
  else {
    s_n = initial_guess;
  }

  // Marching in time until the steady solution has been found
  while ((eps/initial_norm > tol)&&(iter < max_iter)) {
  //while ((eps > tol)&&(iter < max_iter)) {
    update(f,grid,dt,s_n,s_np1);
    eps = norm<T>(s_n,s_np1);
    residuals.push_back(eps/initial_norm);
    if (iter==1) {
      //initial_norm = eps;
      initial_norm = 1.0;
    }
    if (iter==100) {
      auto result = std::max_element(s_np1.u.begin(),s_np1.u.end());
      std::cout << "required time step: " << M_PI*grid.a*grid.dx/(4.0*fabs(*result)) << std::endl;;
    }
    //if (iter%20==0) {
    //  std::cout << std::setw(10) << "iteration" << std::setw(10) << "time " << std::setw(10) << "residual " << std::endl;
    //  std::cout << "------------------------------------------------------------" << std::endl;
    //}
    std::cout << std::setw(10) << iter << std::left << std::setw(10) << t << std::setw(10) << eps/initial_norm << std::endl;
    s_n = soln<T>(s_np1);
    t += dt;
    ++iter;
  }

  // Computing pressure
  update(f,grid,dt,s_n,s_np1,true);

  // Writing surface profiles
  T xpt;
  std::ofstream sfile(surf_file.c_str());
  sfile.precision(16);
  sfile.setf(std::ios_base::scientific);
  for (int i=0; i<imax; ++i) {

    // Transforming back to the physical plane
    xpt = grid.a*tan(M_PI*grid.x[i]/2.0);

    // Only writing out data for -20 <= x <= 20
    if ((xpt>=-20.0)&&(xpt<=20.0)) {
      sfile << xpt << " " << s_np1.tau[i] << " " << s_np1.p[i] << "\n";
    }

  }
  sfile.close();

  // Writing residuals to file
  std::ofstream rfile("residuals.dat");
  rfile.precision(16);
  rfile.setf(std::ios_base::scientific);
  for (T r : residuals) {
    rfile << r << "\n";
  }
  rfile.close();

  return s_np1;

}

#endif
