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

#ifndef INTEGRATION_HEADER_DEF
#define INTEGRATION_HEADER_DEF

/**
 * Function for integrating a function of one variable
 * using the trapezoidal rule.  Must pass something that
 * is callable.  For example, a functor would work.
 * The operator overloading must be const so as to 
 * avoid discarding qualifiers in the below function
 * call.  This function is meant for continuous data.
 */
template <typename T, typename F> 
T integrate_trap_cont(const T a, const T b, const F& f) {
  return (b - a)/2.0*(f(a) + f(b));
}

/**
 * Function for computing the integral using trapezoidal rule
 * which is meant for discrete data.
 */
template <typename T> 
T integrate_trap_disc(const T dx, const T fa, const T fb) {
  return dx/2.0*(fa + fb);
}

#endif
