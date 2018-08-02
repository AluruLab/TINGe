/***
 *  $Id: stat_util.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: stat_util.hpp
 *  Created: Apr 30, 2005
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2007-2008 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying file LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_STAT_UTIL_HPP
#define JAZ_STAT_UTIL_HPP

#include "plain_array.hpp"

#include <algorithm>
#include <iterator>
#include <vector>

#include <math.h>


namespace jaz {

  /** Compute arithmetic mean of the sample [@a first, @a last).
   *  We use the most basic algorithm. In some cases it can be
   *  numerically unstable!
   */
  template <typename T, typename Iter> inline T mean(Iter first, Iter last) {
      T S = 0;
      unsigned int n = std::distance(first, last);
      for (; first != last; ++first) S += *first;
      return (S / n);
  } // mean


  /** Compute unbiased variance of the sample [@a first, @a last).
   *  We use the most basic algorithm. In some cases it can be
   *  numerically unstable!
   *  @param mean must be the mean of the sample [@a first, @a last).
   */
  template <typename T, typename Iter>
  inline T var(Iter first, Iter last, T mean) {
      T S = 0;
      unsigned int n = std::distance(first, last);
      for (; first != last; ++first) S += (*first - mean) * (*first - mean);
      return (S / (n - 1));
  } // var


  template <typename T, typename IterX, typename IterY>
  T cor(IterX first_x, IterX last_x, IterY first_y, IterY last_y) {
      unsigned int i = 1;

      T sum_sq_x = 0.0;
      T sum_sq_y = 0.0;

      T sum_cp = 0.0;

      T mean_x = *(first_x++);
      T mean_y = *(first_y++);

      for (; first_x != last_x; ++first_x, ++first_y, ++i) {
	  T sweep = static_cast<T>(i - 1) / i;

	  T delta_x = *first_x - mean_x;
	  T delta_y = *first_y - mean_y;

	  sum_sq_x += delta_x * delta_x * sweep;
	  sum_sq_y += delta_y * delta_y * sweep;
	  sum_cp += delta_x * delta_y * sweep;

	  mean_x += delta_x / i;
	  mean_y += delta_y / i;
      }

      T sd_x = sqrt(sum_sq_x / i);
      T sd_y = sqrt(sum_sq_y / i);
      T cov_xy = sum_cp / i;

      return cov_xy / (sd_x * sd_y);
  } // cor


  /** Function that standardizes sequence [@a first, @a last).
   *  @param IterIn must be a model of Forward Iterator.
   *  @param IterOut must be a model of Forward Iterator.
   *  @param out must be Output Iterator to a sequence that
   *  can store (@a last - @a first) standardized elements.
   */
  template <typename IterIn, typename IterOut>
  void standardize(IterIn first, IterIn last, IterOut out) {
      typedef typename std::iterator_traits<IterOut>::value_type T;
      T avg = mean<T>(first, last);
      T sd = sqrt(var<T>(first, last, avg));
      for (; first != last; ++first, ++out) *out = (*first - avg) / sd;
  } // standardize


  template <typename T> struct rank_transform_sorter__ {
      unsigned int id;
      T val;

      bool operator<(const rank_transform_sorter__& obj) const {
	  if (val < obj.val) return true;
	  if (obj.val < val) return false;
	  return (id < obj.id);
      } // operator()

  }; // struct rank_transform_sorter__

  /** Functor that converts values in the range [0,n) into (0, 1).
   */
  class zero_to_one : public std::unary_function<unsigned int, double> {
  public:
      explicit zero_to_one(unsigned int n) : n_(n) { }

      double operator()(unsigned int i) const { return (i + 0.5) / n_; }

  private:
      unsigned int n_;

  }; // class zero_to_one

  /** Identity functor.
   */
  struct pos_to_rank : public std::unary_function<unsigned int, unsigned int> {

      unsigned int operator()(unsigned int i) const { return i; }

  }; // pos_to_rank

  /** Performs rank transformation of the sequence [@a first, @a last)
   *  and applies unary function @a op to transformed data.
   *  @param IterIn must be a model of Random Access iterator.
   *  @param IterOut must be a model of Random Access iterator.
   *  @param UnaryFun must be a model of Unary Function.
   *  @param out must be Output Iterator to a sequence
   *  that can store (@a last - @a first) elements.
   */
  template <typename IterIn, typename IterOut, typename UnaryFun>
  void rank_transform(IterIn first, IterIn last, IterOut out, UnaryFun op) {
      typedef typename std::iterator_traits<IterIn>::value_type value_type;

      unsigned int n = std::distance(first, last);
      std::vector<rank_transform_sorter__<value_type> > v(n);

      for (unsigned int i = 0; i < n; ++i) {
	  v[i].id = i;
	  v[i].val = first[i];
      }

      std::sort(v.begin(), v.end());

      for (unsigned int i = 0; i < n; ++i) {
	  *(out + v[i].id) = op(i);
      }
  } // rank_transform

  template <typename IterIn, typename IterOut>
  void rank_transform(IterIn first, IterIn last, IterOut out) {
      rank_transform(first, last, out, pos_to_rank());
  } // rank_transform

} // namespace jaz

#endif // JAZ_STAT_UTIL_HPP
