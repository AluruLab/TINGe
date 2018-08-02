/***
 *  $Id: numeric_table.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: numeric_table.hpp
 *  Created: Dec 12, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2007-2008 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef NUMERIC_TABLE_HPP
#define NUMERIC_TABLE_HPP

#include "plain_array.hpp"
#include <iostream>
#include <string>
#include <vector>


namespace jaz {

  template <typename T> class numeric_table {
  public:
      /**
       */
      typedef std::size_t size_type;

      /**
       */
      typedef T value_type;

      /**
       */
      typedef T* row_iterator;

      /**
       */
      typedef const T* const_row_iterator;


      numeric_table() : nrow_(0), ncol_(0) { }

      numeric_table(size_type nr, size_type nc) : nrow_(nr), ncol_(nc) {
	  data_.resize(nrow_ * ncol_);
	  row_names_.resize(nrow_);
	  col_names_.resize(ncol_);
      } // numeric_table


      // TODO:
      // insert_row
      // insert_col
      // col_iterator


      template <typename Iter>
      bool add_row(Iter first, Iter last, const std::string& name = "") {
	  if ((ncol_ != 0) && (last - first) != ncol_) return false;
	  if ((last - first) == 0) return false;

	  row_names_.push_back(name);

	  data_.resize(nrow_ * ncol_ + (last - first));
	  for (size_type i = 0; i < last - first; ++i) {
	      data_[nrow_ * ncol_ + i] = first[i];
	  }

	  ++nrow_;

	  ncol_ = (last - first);
	  col_names_.resize(ncol_, "");

	  return true;
      } // add_row

      template <typename Iter>
      bool add_col(Iter first, Iter last, const std::string& name = "") {
	  if ((nrow_ != 0) && (last - first) != nrow_) return false;
	  if ((last - first) == 0) return false;

	  col_names_.push_back(name);

	  data_.resize(nrow_ * ncol_ + (last - first));
	  data_.insert(data_.begin() + ncol_, *first);
	  ++ncol_;

	  for (size_type i = 1; i < last - first; ++i) {
	      data_.insert(data_.begin() + (i + 1) * ncol_ - 1, first[i]);
	  }

	  nrow_ = (last - first);
	  row_names_.resize(nrow_, "");

	  return true;
      } // add_col

      void row_name(size_type n, const std::string& s) { row_names_[n] = s; }

      void col_name(size_type n, const std::string& s) { col_names_[n] = s; }

      template <typename Iter> bool row_names(Iter first, Iter last) {
	  if (static_cast<unsigned int>(last - first) != nrow_) return false;
	  std::copy(first, last, row_names_.begin());
	  return true;
      } // row_names

      template <typename Iter> bool col_names(Iter first, Iter last) {
	  if (static_cast<unsigned int>(last - first) != ncol_) return false;
	  std::copy(first, last, col_names_.begin());
	  return true;
      } // col_names

      // *** ACCESS INTERFACE ***

      size_type nrow() const { return nrow_; }

      size_type ncol() const { return ncol_; }

      std::string row_name(size_type n) const { return row_names_[n]; }

      std::string col_name(size_type n) const { return col_names_[n]; }

      const_row_iterator row_begin(size_type n) const {
	  return data_.begin() + n * ncol_;
      } // row_begin

      const_row_iterator row_end(size_type n) const {
	  return data_.begin() + (n + 1) * ncol_;
      } // row_end

      row_iterator row_begin(size_type n) { return data_.begin() + n * ncol_; }

      row_iterator row_end(size_type n) { return data_.begin() + (n + 1) * ncol_; }


      row_iterator operator[](size_type n) { return data_.begin() + n * ncol_; }

      const_row_iterator operator[](size_type n) const {
	  return data_.begin() + n * ncol_;
      } // operator[]


  private:
      size_type nrow_;
      size_type ncol_;

      plain_array<value_type> data_;

      std::vector<std::string> row_names_;
      std::vector<std::string> col_names_;

  }; // class numeric_table


  template <typename T>
  inline std::ostream& operator<<(std::ostream& os, const numeric_table<T>& t) {
      for (unsigned int i = 0; i < t.ncol(); ++i) {
	  os << "\t" << t.col_name(i);
      }

      for (unsigned int i = 0; i < t.nrow(); ++i) {
	  os << "\n";
	  os << t.row_name(i) << std::flush;
	  typename numeric_table<T>::const_row_iterator iter(t.row_begin(i));
	  for (unsigned int j = 0; j < t.ncol(); ++j) {
	      os << "\t" << iter[j];
	  }
      }
      return os;
  } // operator <<

} // namespace jaz

#endif // NUMERIC_TABLE_HPP
