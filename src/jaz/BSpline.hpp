/***
 *  $Id: BSpline.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: BSpline.hpp
 *  Created: May 24, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2009 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_B_SPLINE_HPP
#define JAZ_B_SPLINE_HPP

#include "plain_array.hpp"
#include <cstddef>
#include <limits>


namespace jaz {

  /** This class allows to evaluate B-spline functions.
   *  @param Float must be a buil-it floating point type.
   */
  template <typename Float> class BSpline {
  public:
      /**
       */
      typedef unsigned short int index_type;

      /**
       */
      typedef Float value_type;


      /** Construct B-spline.
       *  @a k is a spline order. If @a k < 1 behavior is undefined.
       *  [@a first, @a last) is sequence of knot points with requirement
       *  that *(@a first + i) <= *(@a first + i + 1). If this condition
       *  is not satisfied result is undefined. This class requires
       *  O(@a last - @a first) memory.
       */
      template <typename Iter>
      BSpline(index_type k, Iter first, Iter last)
	  : k_(k), t_(first, last), B_(t_.size() + 1) { }


      /** @return order of the B-spline.
       */
      index_type k() const { return k_; }


      /** Compute k non-zero basis functions for @a x using de Boor algorithm.
       *  @param x is a point for which basis functions should be computed,
       *  and t[first + k - 1] <= x <= t[last - k], where t stores knot points.
       *  If x is not in this  range result is undefined.
       *  @return index to first value of k base functions.
       *  This index can be passed to operator[] to retrieve the value.
       */
      index_type operator()(value_type x) {
	  B_.zero();

	  unsigned int n = t_.size();
	  n -= k_;

	  // case: k = 1
	  index_type idx = 0;

          // assure that right most element is detected correctly
	  if ((x - t_[n] <= std::numeric_limits<float>::epsilon())
	      && (x - t_[n]) >= -std::numeric_limits<float>::epsilon()) {
	      idx = n - 1;
	      B_[idx] = 1;
	  } else {
	      for (idx = k_ - 1; idx < n; ++idx) {
		  if ((t_[idx] <= x) && (x < t_[idx + 1])) {
		      B_[idx] = 1;
		      break;
		  }
	      }
	  }

	  value_type t_a, t_b, t_den[2];

	  // case: k > 1
	  for (index_type k = 2; k <= k_; ++k) {
	      for (index_type i = idx - k + 1; i <= idx; ++i) {
		  t_den[0] = (t_[i + k - 1] - t_[i]);
		  t_den[1] = (t_[i + k] - t_[i + 1]);

		  t_a = (t_den[0] != 0) ? B_[i] * (x - t_[i]) / t_den[0] : 0;
		  t_b = (t_den[1] != 0) ? B_[i + 1] * (t_[i + k] - x) / t_den[1] : 0;

		  B_[i] = t_a + t_b;
	      }
	  }

	  return idx - k_ + 1;
      } // operator()


      /** @return value of the i-th base function computed by operator().
       *  If operator() was not called result is undefined.
       */
      const value_type& operator[](index_type i) const { return B_[i]; }


  private:
      index_type k_;
      jaz::plain_array<value_type> t_;
      jaz::plain_array<value_type> B_;

  }; // class BSpline

} // namespace jaz

#endif // JAZ_B_SPLINE_HPP
