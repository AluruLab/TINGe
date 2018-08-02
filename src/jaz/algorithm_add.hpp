/***
 *  $Id: algorithm_add.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: algorithm_add.hpp
 *  Created: May 22, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2008 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_ALGORITHM_ADD_HPP
#define JAZ_ALGORITHM_ADD_HPP

#include <cstddef>
#include <iterator>
#include <ostream>
#include <string>
#include <utility>
#include <stdlib.h>


namespace jaz {

  /** Function to find smallest and largest element
   *  of the sequence [@a first, @a last).
   *  @param Iter must be a model of forward iterator.
   */
  template <typename Iter>
  std::pair<Iter, Iter> min_max_element(Iter first, Iter last) {
      std::pair<Iter, Iter> p(first, first);
      ++first;

      for (; first != last; ++first) {
	  if (*first < *p.first) p.first = first;
	  else if (*p.second < *first) p.second = first;
      }

      return p;
  } // min_max_element

  template <typename Iter, typename BinPred>
  std::pair<Iter, Iter> min_max_element(Iter first, Iter last, BinPred comp) {
      std::pair<Iter, Iter> p(first, first);
      ++first;

      for (; first != last; ++first) {
	  if (comp(*first, *p.first)) p.first = first;
	  else if (comp(*p.second, *first)) p.second = first;
      }

      return p;
  } // min_max_element


  /** Function to check if the sequence [@a first, @a last) is sorted.
   *  @param Iter must be a model of forward iterator.
   */
  template <typename Iter> bool is_sorted(Iter first, Iter last) {
      Iter next = first;
      ++next;
      for (; next != last; ++first, ++next) if (!(*first < *next)) return false;
      return true;
  } // is_sorted

  template <typename Iter, typename BinPred>
  bool is_sorted(Iter first, Iter last, BinPred comp) {
      Iter next = first;
      ++next;

      for (; next != last; ++first, ++next) {
	  if (comp(*first, *next) == false) return false;
      }

      return true;
  } // is_sorted


  /** Modulo n implementation of RandomNumberGenerator.
   *  @param T must be an integer type.
   *  @param n must be a non-negative value.
   *  @return number in the range [0, n).
   */
  template <typename T> T random(T n) { return rand() % n; }


  /**
   */
  template <typename InputIter, typename OutputIter>
  OutputIter random_sample_n(InputIter first, InputIter last, OutputIter out,
			     std::size_t n) {
      std::size_t r = std::distance(first, last);
      std::size_t m = std::min(n, r);

      while (m > 0) {
	  if (random(r) < m) {
	      *out = *first;
	      ++out;
	      --m;
	  }
	  --r;
	  ++first;
      }

      return out;
  } // random_sample_n

} // namespace jaz

#endif // JAZ_ALGORITHM_ADD_HPP
