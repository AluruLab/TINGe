/***
 *  $Id: string_add.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: string_add.hpp
 *  Created: Jun 03, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2004-2008 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_STRING_ADD_HPP
#define JAZ_STRING_ADD_HPP

#include <cstddef>
#include <iostream>
#include <functional>
#include <locale>
#include <sstream>
#include <string>


namespace jaz {

  template <typename charT, typename traits, typename Alloc, typename Iter>
  void split(charT pat,
	     const std::basic_string<charT, traits, Alloc>& s, Iter out) {
      unsigned int pos = 0;

      for (unsigned int i = 0; i < s.size(); ++i) {
	  if (s[i] == pat) {
	      if (i - pos > 0) {
		  *(out++) =
		      std::basic_string<charT, traits, Alloc>(s, pos, i - pos);
	      }
	      pos = i + 1;
	  }
      }

      *(out++) =
	  std::basic_string<charT, traits, Alloc>(s, pos, s.size() - pos);
  } // split


  template <typename charT, typename traits, typename Alloc>
  std::basic_string<charT, traits, Alloc>
  remove_space(const std::basic_string<charT, traits, Alloc>& s,
	       const std::locale& loc = std::locale::classic()) {
      std::basic_string<charT, traits, Alloc> ns;
      for (std::size_t i = 0; i < s.size(); ++i) {
	  if (std::isspace(s[i], loc) == false) ns.push_back(s[i]);
      }
      return ns;
  } // remove_space


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_uc_compare :
	public std::binary_function<std::basic_string<charT, traits, Alloc>,
				    std::basic_string<charT, traits, Alloc>,
				    int> {
  public:
      explicit basic_uc_compare(const std::locale& loc = std::locale::classic())
	  : ct_(std::use_facet<std::ctype<char> >(loc)) { }

      int operator()(const std::basic_string<charT, traits, Alloc>& s1,
		     const std::basic_string<charT, traits, Alloc>& s2) const {

	  std::size_t l1 = s1.size();
	  std::size_t l2 = s2.size();

	  if (l1 < l2) return -1; else if (l2 < l1) return 1;

	  charT c1, c2;

	  for (std::size_t i = 0; i < l1; ++i) {
	      c1 = ct_.toupper(s1[i]);
	      c2 = ct_.toupper(s2[i]);
	      if (c1 < c2) return -1; else if (c2 < c1) return 1;
	  }

	  return 0;
      } // operator()

  private:
      const std::ctype<charT>& ct_;

  }; // class basic_uc_compare

  typedef basic_uc_compare<char> uc_compare;


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_to_upper :
	public std::unary_function<std::basic_string<charT, traits, Alloc>,
				   std::basic_string<charT, traits, Alloc> > {
  public:
      explicit basic_to_upper(const std::locale& L = std::locale::classic())
	  : ct_(std::use_facet<std::ctype<char> >(L)), s_() { }

      const std::basic_string<charT, traits, Alloc>&
      operator()(const std::basic_string<charT, traits, Alloc>& s) const {
	  s_ = s;
	  for (std::size_t i = 0; i < s.size(); ++i) s_[i] = ct_.toupper(s[i]);
	  return s_;
      } // operator

  private:
      const std::ctype<charT>& ct_;
      mutable std::basic_string<charT, traits, Alloc> s_;

  }; // class basic_to_upper

  typedef basic_to_upper<char> to_upper;


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_to_string {
  public:

      template <typename T>
      std::basic_string<charT, traits, Alloc> operator()(const T& t) {
	  os_.clear();
	  os_.str("");
	  os_ << t;
	  return os_.str();
      } // operator()

  private:
      std::basic_ostringstream<charT, traits, Alloc> os_;

  }; // class basic_to_string

  typedef basic_to_string<char> to_string;


  template <typename charT,
	    typename traits = std::char_traits<charT>,
	    typename Alloc = std::allocator<charT> >
  class basic_string_to {
  public:

      template <typename T>
      bool operator()(const std::basic_string<charT, traits, Alloc>& s, T& t) {
	  is_.clear();
	  is_.str(s + '\n');
	  is_ >> t;
	  return is_.good();
      } // operator()

  private:
      std::basic_istringstream<charT, traits, Alloc> is_;

  }; // class basic_string_to

  typedef basic_string_to<char> string_to;

} // namespace jaz

#endif // JAZ_STRING_ADD_HPP
