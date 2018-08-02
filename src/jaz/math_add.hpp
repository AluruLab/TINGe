/***
 *  $Id: math_add.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: math_add.hpp
 *  Developed: Mar 11, 2005
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2008 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_MATH_ADD_HPP
#define JAZ_MATH_ADD_HPP

#include "plain_array.hpp"
#include <functional>
#include <math.h>
#include <stddef.h>


namespace jaz {

  /** Constant Pi.
   */
  template <typename T> inline T pi() {
      static const T PI = static_cast<T>(3.141592653589793238462643383279502884197169399375105820974944);
      return PI;
  }; // struct PI


  /** Simple, iterative function returning factorial of x.
   */
  inline unsigned long int fac(unsigned short int x) {
      unsigned long int f = 1;
      for (int i = 2; i <= x; ++i) f *= i;
      return f;
  } // fac


  /** Simple functor to compute binomial coefficient B(n,k).
   *  It requires O(n) memory.
   */
  template <typename T> class Bin : public std::binary_function<T, T, T> {
  public:

      T operator()(T n, T k) const {
	  ++n;
	  bi_.resize(n);

	  bi_[0] = 1;
	  for (T i = 1; i < n; ++i) {
	      bi_[i] = 1;
	      for (T j = i - 1; j > 0; --j) bi_[j] += bi_[j - 1];
	  }

	  return bi_[k];
      } // operator()

  private:
      mutable plain_array<T> bi_;

  }; // class Bin


  /** Class that generates all possible B(n,k) combinations
   *  of the set [0,n). It requires O(k) memory.
   */
  class Combination {
  public:
      /**
       */
      typedef unsigned int value_type;

      /** Create Combination object and set it to store first
       *  k-element combination over n-element set.
       */
      explicit Combination(unsigned int n = 0, unsigned int k = 0) { reset(n, k); }

      /** Reset to first k-element combination over n-element set.
       */
      void reset(unsigned int n, unsigned int k) {
	  n_ = n;
	  k_ = k;
	  k1_ = k - 1;
	  nk_ = n - k;
	  data_.resize(k_);
	  for (unsigned int i = 0; i < k_; ++i) data_[i] = i;
      } // reset

      /** Generate next combination (in the lexicographical order).
       *  @return true if the combination is available, false otherwise.
       */
      bool next() {
	  unsigned int i = k1_;
	  while ((i > 0) && (data_[i] == nk_ + i)) --i;
	  if ((i == 0) && (data_[i] == nk_)) return false;
	  ++data_[i];
	  for (; i < k1_; ++i) data_[i + 1] = data_[i] + 1;
	  return true;
      } // next

      /** @return pointer to k-element array with generated combination.
       *  This pointer is changed only by reset() function.
       */
      const value_type* data() const { return data_.begin(); }

      /** @ return i-th element of the current combination.
       */
      value_type operator[](unsigned int i) const { return data_[i]; }


  private:
      unsigned int n_;
      unsigned int k_;
      unsigned int k1_;
      unsigned int nk_;
      jaz::plain_array<value_type> data_;

  }; // class Combination

} // namespace jaz

#endif // JAZ_MATH_ADD_HPP
