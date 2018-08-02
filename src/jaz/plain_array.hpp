/***
 *  $Id: plain_array.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: plain_array.hpp
 *  Created: May 26, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2009 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_PLAIN_ARRAY_HPP
#define JAZ_PLAIN_ARRAY_HPP

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <new>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


namespace jaz {

  template <typename T> class plain_array {
  public:
      /**
       */
      typedef std::size_t size_type;

      /**
       */
      typedef T value_type;

      /**
       */
      typedef T& reference;

      /**
       */
      typedef const T& const_reference;

      /**
       */
      typedef T* iterator;

      /**
       */
      typedef const T* const_iterator;


      // *** CONSTRUCTORS + CO. ***

      explicit plain_array(size_type n = 0) : n_(n), q_(n), tab_(0) {
	  if (n != 0) {
	      tab_ = static_cast<T*>(malloc(n * sizeof(T)));
	      if (tab_ == 0) throw std::bad_alloc();
	  }
      } // plain_array

      plain_array(size_type n, T t) : n_(n), q_(n), tab_(0) {
	  if (n != 0) {
	      tab_ = static_cast<T*>(malloc(n * sizeof(T)));
	      if (tab_ == 0) throw std::bad_alloc();
	      for (size_type i = 0; i < n_; ++i) tab_[i] = t;
	  }
      } // plain_array

      template <typename Iter> plain_array(Iter first, Iter last) : tab_(0) {
	  n_ = std::distance(first, last);
	  q_ = n_;
	  if (n_ != 0) {
	      tab_ = static_cast<T*>(malloc(n_ * sizeof(T)));
	      if (tab_ == 0) throw std::bad_alloc();
	      std::copy(first, last, tab_);
	  }
      } // plain_array

      plain_array(const plain_array& a) : tab_(0) {
	  n_ = a.n_;
	  q_ = n_;
	  if (n_ != 0) {
	      tab_ = static_cast<T*>(malloc(n_ * sizeof(T)));
	      if (tab_ == 0) throw std::bad_alloc();
	      memcpy(tab_, a.tab_, n_ * sizeof(T));
	  }
      } // plain_array

      // move constructor
      plain_array(plain_array& a, int) {
	  n_ = a.n_;
	  q_ = a.q_;
	  tab_ = a.tab_;
	  a.n_ = 0;
	  a.q_ = 0;
	  a.tab_ = 0;
      } // plain_array_

      ~plain_array() { free(tab_); }

      void operator=(const plain_array& a) {
	  if (this != &a) {
	      if (a.n_ != 0) {
		  value_type* ntab = static_cast<T*>(realloc(tab_, a.n_ * sizeof(T)));
		  if (ntab == 0) throw std::bad_alloc();
		  n_ = a.n_;
		  q_ = n_;
		  tab_ = ntab;
		  memcpy(tab_, a.tab_, n_ * sizeof(T));
	      } else {
		  clear();
	      }
	  }
      } // operator=


      // *** STORAGE MANAGEMENT ***

      void zero() { memset(tab_, 0, n_ * sizeof(T)); }

      void clear() {
	  n_ = 0;
	  q_ = 0;
	  free(tab_);
	  tab_ = 0;
      } // clear

      void reserve(size_type n) {
	  if (n_ < n) {
	      value_type* ntab = static_cast<T*>(realloc(tab_, n * sizeof(T)));
	      if (ntab == 0) throw std::bad_alloc();
	      q_ = n;
	      tab_ = ntab;
	  }
      } // reserve

      void resize(size_type n) {
	  if (n_ != n) {
	      if (n == 0) clear();
	      else if (n + BLCK_SZ < q_) {
		  value_type* ntab = static_cast<T*>(realloc(tab_, n * sizeof(T)));
		  if (ntab == 0) throw std::bad_alloc();
		  q_ = n;
		  tab_ = ntab;
	      }
	      else if (q_ < n) {
		  unsigned int nq = n + std::min<size_type>(BLCK_SZ, static_cast<size_type>(1.25 * q_));
		  value_type* ntab = static_cast<T*>(realloc(tab_, nq * sizeof(T)));
		  if (ntab == 0) throw std::bad_alloc();
		  q_ = nq;
		  tab_ = ntab;
	      }
              n_ = n;
	  }
      } // resize

      void resize(size_type n, T t) {
	  size_type p = n_;
	  resize(n);
	  for (size_type i = p; i < n_; ++i) { tab_[i] = t; }
      } // resize

      void push_back(const T* a, size_type an) {
	  size_type n = n_;
	  resize(n + an);
	  memcpy(tab_ + n, a, an * sizeof(T));
      } // push_back

      void push_back(const plain_array& a) {
	  size_type n = n_;
	  resize(n + a.n_);
	  memcpy(tab_ + n, a.tab_, a.n_ * sizeof(T));
      } // push_back

      void push_back(T t) {
	  size_type n = n_;
	  if (q_ < n + 1) {
	      unsigned int nq = n + std::min<size_type>(BLCK_SZ, static_cast<size_type>(1.25 * q_)) + 1;
	      value_type* ntab = static_cast<T*>(realloc(tab_, nq * sizeof(T)));
	      if (ntab == 0) throw std::bad_alloc();
	      q_ = nq;
	      tab_ = ntab;
	  }
	  n_++;
	  tab_[n] = t;
      } // push_back

      iterator insert(iterator pos, T t) {
	  size_type p = pos - tab_;
	  resize(n_ + 1);
	  if (p < n_) memmove(tab_ + p + 1, tab_ + p, (n_ - p - 1) * sizeof(T));
	  tab_[p] = t;
	  return tab_ + p;
      } // insert

      void insert(iterator pos, iterator first, iterator last) {
	  size_type n = n_;
	  size_type p = pos - tab_;
	  size_type m = last - first;
	  if (pos == (tab_ + n_)) push_back(first, m);
	  else {
	      resize(n_ + m);
	      memmove(tab_ + p + m, tab_ + p, (n - p) * sizeof(T));
	      memcpy(tab_ + p, first, m * sizeof(T));
	  }
      } // insert

      iterator erase(iterator pos) {
	  size_type p = pos - tab_;
	  if (p < (n_ - 1)) memmove(pos, tab_ + p + 1, (n_ - p - 1) * sizeof(T));
	  resize(n_ - 1);
	  return tab_ + p;
      } // erase

      iterator erase(iterator first, iterator last) {
	  size_type m = n_ - (last - tab_);
	  if (last < tab_ + n_) memmove(first, last, m * sizeof(T));
	  resize(n_ - (last - first));
	  return first;
      } // erase

      // *** ACCESS INTERFACE ***

      iterator begin() { return tab_; }

      iterator end() { return tab_ + n_; }

      const_iterator begin() const { return tab_; }

      const_iterator end() const { return tab_ + n_; }

      bool empty() const { return (n_ == 0); }

      size_type size() const { return n_; }

      const T& operator[](size_type n) const { return tab_[n]; }

      T& operator[](size_type n) { return tab_[n]; }


  private:
      enum { BLCK_SZ = 4096 };

      size_type n_;
      size_type q_;
      value_type* tab_;

  }; // class plain_array

} // namespace jaz

#endif // JAZ_PLAIN_ARRAY_HPP
