/***
 *  $Id: mutual_info.hpp 596 2010-06-17 16:19:55Z zola $
 **
 *  File: mutual_info.hpp
 *  Created: Apr 10, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2009 Jaroslaw Zola
 *  Distributed under the Boost Software License.
 *  See accompanying file LICENSE.
 *
 *  This file is part of jaz.
 */

#ifndef MUTUAL_INFO_HPP
#define MUTUAL_INFO_HPP

#include "BSpline.hpp"
#include "math_add.hpp"
#include <math.h>


namespace jaz {

  /**
   */
  template <typename Float> class bs_mutual_info {
  private:
      typedef BSpline<Float> spline_type;
      typedef typename spline_type::index_type index_type;


  public:
      /**
       */
      typedef Float value_type;


      /**
       */
      bs_mutual_info(unsigned int m, index_type b, index_type k)
	  : m_(m), b_(b), k_(k), xi_(0), shash_idx_(m), shash_val_(m * k), P_(b * b) {
	  bsqr_ = b * b;
	  logm_ = log2(static_cast<value_type>(m));
	  hash_spline_();
	  H2_ = 2 * marginal_();
      } // bs_mutual_info


      /**
       */
      value_type operator()(const unsigned int* yi) {
	  return this->operator()(xi_, yi);
      } // operator()

      /**
       */
      value_type operator()(const unsigned int* xi, const unsigned int* yi) {
	  P_.zero();

	  const value_type* B = shash_val_.begin();

	  const value_type* bx = 0;
	  const value_type* by = 0;

	  for (unsigned int l = 0; l < m_; ++l) {
	      unsigned int pos = xi[l];

	      bx = B + (pos * k_);
	      by = B + (yi[l] * k_);

	      index_type off = shash_idx_[yi[l]];
	      value_type* p = P_.begin() + shash_idx_[pos] * b_;

	      for (index_type i = 0; i < k_; ++i) {
		  for (index_type j = 0; j < k_; ++j) p[off + j] += bx[i] * by[j];
		  p += b_;
	      } // i

	  } // l

	  xi_ = xi;

	  return (H2_ - entropy_(bsqr_)) / m_;
      } // operator()


  private:
      void hash_spline_() {
	  shash_idx_.zero();
	  shash_val_.zero();

	  jaz::plain_array<index_type> t(b_ + k_);

	  for (index_type i = 0; i < k_; ++i) t[i] = 0;
	  for (index_type i = k_; i < b_; ++i) t[i] = i - k_ + 1;
	  for (index_type i = b_; i < b_ + k_; ++i) t[i] = b_ - k_ + 1;

	  spline_type bs(k_, t.begin(), t.end());

	  value_type mk_norm = (b_ - k_ + 1.0) / (m_ - 1);
	  value_type* sh = shash_val_.begin();

	  for (unsigned int i = 0; i < m_; ++i, sh += k_) {
	      value_type x = i * mk_norm;
	      index_type idx = shash_idx_[i] = bs(x);
	      memcpy(sh, &bs[idx], k_ * sizeof(value_type));
	  }
      } // hash_spline_

      value_type entropy_(index_type b) {
	  value_type p = 0.0;
	  value_type H = 0.0;

	  for (index_type i = 0; i < b; ++i) {
	      p = P_[i];
	      if (std::numeric_limits<float>::epsilon() < p) {
		  H += p * (log2(p) - logm_);
	      }
	  }

	  // *** we do not divide by m_ here ***
	  // you must do it later when H is used

	  return -H;
      } // entropy_

      value_type marginal_() {
	  P_.zero();

	  for (unsigned int i = 0; i < m_; ++i) {
	      value_type* p = P_.begin() + shash_idx_[i];
	      value_type* B = shash_val_.begin() + i * k_;
	      for (index_type j = 0; j < k_; ++j) p[j] += B[j];
	  }

	  return entropy_(b_);
      } // marginal_

      unsigned int m_;
      index_type b_;
      index_type k_;

      index_type bsqr_;
      value_type logm_;
      value_type H2_;

      const unsigned int* xi_;

      jaz::plain_array<index_type> shash_idx_;
      jaz::plain_array<value_type> shash_val_;
      jaz::plain_array<value_type> P_;

  }; // class bs_mutual_info


  /**
   */
  template <typename Float> class gk_mutual_info {
  public:
      /**
       */
      typedef Float value_type;


      /**
       */
      gk_mutual_info(unsigned int m, value_type h_1, value_type h_2)
	  : m_(m), xi_(0), ehash_(m) {
	  hash_exp_(h_1);
	  H2_ = 2 * marginal_(h_1);
	  hash_exp_(h_2);
	  Kxy_ = m_ * log2(1.0 / (m_ * 2 * pi<value_type>() * h_2 * h_2));
      } // gs_mutual_info


      /**
       */
      value_type operator()(const unsigned int* yi) {
	  return this->operator()(xi_, yi);
      } // operator()

      /**
       */
      value_type operator()(const unsigned int* xi, const unsigned int* yi) {
	  value_type H = 0;

	  for (unsigned int i = 0; i < m_; ++i) {
	      unsigned int x = xi[i];
	      unsigned int y = yi[i];

	      value_type S = 0;

	      for (unsigned int j = 0; j < m_; ++j) {
		  unsigned int dx = abs(static_cast<int>(x - xi[j]));
		  unsigned int dy = abs(static_cast<int>(y - yi[j]));
		  S += (ehash_[dx] * ehash_[dy]);
	      }

	      H += log2(S);
	  } // for i

	  H += Kxy_;
	  xi_ = xi;

	  return (H - H2_) / m_;
      } // operator()


  private:
      void hash_exp_(value_type h) {
	  value_type h2_sqr = 1.0 / (2 * h * h);
	  for (unsigned int i = 0; i < m_; ++i) {
	      // value_type x = (i + 0.5) / m_;
	      // ehash_[i] = 1.0 / exp(x * x * h2_sqr);
	      ehash_[i] = 1.0 / exp(i * i * h2_sqr);
	  }
      } // hash_exp_

      value_type marginal_(value_type h) {
	  value_type H = 0;
	  value_type K = m_ * log2(1.0 / (m_ * sqrt(2 * pi<value_type>()) * h));

	  for (unsigned int i = 0; i < m_; ++i) {
	      value_type S = 0;
	      for (unsigned int j = 0; j < m_; ++j) {
		  S += ehash_[abs(static_cast<int>(i - j))];
	      }
	      H += log2(S);
	  }

	  H += K;

	  // *** we do not divide by m_ here ***
	  // you must do it later when H is used

	  return H;
      } // marginal_

      unsigned int m_;

      value_type Kxy_;
      value_type H2_;

      const unsigned int* xi_;

      jaz::plain_array<value_type> ehash_;

  }; // class gk_mutual_info


  /**
   */
  template <typename Float, int D> Float h_silverman(unsigned int m) {
      Float sd = sqrt(static_cast<Float>(m * m + m) / 12);
      return pow(static_cast<Float>(4.0 / ((D + 2) * m)), 1.0 / (D + 4)) * sd;
  } // h_silverman

} // namespace jaz

#endif // MUTUAL_INFO_HPP
