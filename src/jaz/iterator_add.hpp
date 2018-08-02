/***
 *  $Id: iterator_add.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: iterator_add.hpp
 *  Created: Apr 28, 2010
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2010 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef JAZ_ITERATOR_ADD_HPP
#define JAZ_ITERATOR_ADD_HPP

#include <iterator>
#include <string>


namespace jaz {

  template <typename charT = char,
	    typename traits = std::char_traits<charT>,
	    typename dist = std::ptrdiff_t>
  class getline_iterator : public std::iterator<std::input_iterator_tag, std::basic_string<charT, traits>, dist> {
  public:
      typedef charT char_type;
      typedef traits traits_type;
      typedef std::basic_string<charT, traits> value_type;
      typedef std::basic_istream<char_type, traits_type> istream_type;

      getline_iterator() : delim_(), value_(), state_(false), is_(0) { }

      getline_iterator(istream_type& is) : is_(&is) {
	  delim_ = std::use_facet<std::ctype<char_type> >(is_->getloc()).widen('\n');
	  read__();
      } // getline_iterator

      getline_iterator(istream_type& is, char_type delim)
	  : delim_(delim), is_(&is) { read__(); }

      getline_iterator(const getline_iterator& gi)
	  : delim_(gi.delim_), value_(gi.value_), state_(gi.state_), is_(gi.is_) { }


      const value_type& operator*() const { return value_; }

      const value_type* operator->() const { return &(operator*()); }


      getline_iterator& operator++() {
	  read__();
	  return *this;
      } // operator++

      getline_iterator operator++(int) {
	  getline_iterator tmp = *this;
	  read__();
	  return tmp;
      } // operator++


  private:
      void read__() {
	  state_ = (is_ && *is_) ? true : false;
	  if (is_) {
	      std::getline(*is_, value_, delim_);
	      state_ = *is_ ? true : false;
	  }
      } // read__

      char_type delim_;
      value_type value_;

      bool state_;
      istream_type* is_;


      friend bool operator==(const getline_iterator& lhs, const getline_iterator& rhs) {
	  return ((lhs.state_ == rhs.state_) && (!lhs.state_ || (lhs.is_ == rhs.is_)));
      } // operator==

      friend bool operator!=(const getline_iterator& lhs, const getline_iterator& rhs) {
	  return !(lhs == rhs);
      } // operator!=

  }; // class getline_iterator

} // namespace jaz

#endif // JAZ_ITERATOR_ADD_HPP
