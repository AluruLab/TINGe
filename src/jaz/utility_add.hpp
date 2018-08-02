/***
 *  $Id: utility_add.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: utility_add.hpp
 *  Created: Jun 18, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2009 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying file LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef UTILITY_ADD_HPP
#define UTILITY_ADD_HPP

#include <iostream>
#include <utility>


namespace jaz {

  /**
   */
  template <typename T> class print {
  public:
      explicit print(const std::string& sep) : os_(std::cout), sep_(sep) { }

      explicit print(std::ostream& os = std::cout, const std::string& sep = "\n")
	  : os_(os), sep_(sep) { }

      void operator()(const T& t) const { os_ << t << sep_; }

  private:
      std::ostream& os_;
      std::string sep_;

  }; // class print


  template <typename T1, typename T2>
  class tie_proxy__ {
  public:
      tie_proxy__(T1* t1, T2* t2) : t1_(t1), t2_(t2) { }

      template <typename U1, typename U2>
      void operator=(const std::pair<U1, U2>& p) {
	  *t1_ = p.first;
	  *t2_ = p.second;
      } // operator

  private:
      T1* t1_;
      T2* t2_;

  }; // class tie_proxy__

  /** This function mimics std::tr1::tie for std::pair.
   */
  template <typename T1, typename T2>
  inline tie_proxy__<T1, T2> tie(T1& t1, T2& t2) {
      return tie_proxy__<T1, T2>(&t1, &t2);
  } // tie

} // namespace jaz

#endif // UTILITY_ADD_HPP
