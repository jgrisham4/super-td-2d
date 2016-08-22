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

#ifndef VISCOSITYHEADER
#define VISCOSITYHEADER

// Function for Sutherland's law
template <typename U>
U mu(U T) {
  U T0,S,mu0;
  T0 = (U) 273.1;
  S = (U) 110.6;
  mu0 = (U) 1.716e-5;
  return (U) mu0*(pow(T/T0,1.5)*(T0+S)/(T+S));
}

#endif
