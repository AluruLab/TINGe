/***
 *  $Id: AppConfig.hpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: AppConfig.hpp
 *  Created: Mar 28, 2008
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *
 *  Copyright 2007-2010 Jaroslaw Zola
 *
 *  This file is part of TINGe.
 *
 *  TINGe is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TINGe is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TINGe. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef APP_CONFIG_HPP
#define APP_CONFIG_HPP

#include <jaz/string_add.hpp>
#include <jaz/sys_tools.hpp>
#include "config.h"
#include "iomanip.hpp"
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <limits>
#include <unistd.h>


struct AppConfig {
    enum convert_t { NONE, CORRELATION };
    enum estimator_t { BSPLINE, GAUSSIAN };

    std::string convert2str(convert_t t) const {
	switch (t) {
	  case CORRELATION : return "correlation";
	}
	return "none";
    } // convert2str

    std::string estimator2str(estimator_t t) const {
	switch (t) {
	  case BSPLINE : return "B-spline";
	  case GAUSSIAN : return "Gaussian kernel";
	}
	return "unknown";
    } // convert2str

    AppConfig() {
	has_pubsetbuf = jaz::iss_pubsetbuf_test();

	input = "";
	output = "";
	output_mi = "";
	input_mi = "";
	input_tf = "";

	output_lab = "";

	mi_estim = BSPLINE;
	mi_b = 10;
	mi_k = 4;
	mi_pval = 1.0;
	mi_thr = 0.0;
	dpi_tol = 1.0;
	mi_conv = NONE;

	boot_size = 0;
	sample_size = 0;
	rng_seed = -1;

	ptest_size = 10;

	mem_report = false;
	verbose = false;
    } // AppConfig


    void print_help(std::ostream& os) {
	os << std::endl;
	os << "Available options:\n\n";
	os << "\t-i <file> \t read expression data from file\n";
	os << "\t-o <file> \t output final network to file\n";
	os << "\t-j <file> \t read MI relevance network from file\n";
	os << "\t-w <file> \t write MI relevance network to file\n";
	os << "\t-l <file> \t get list of TFs from file\n";
	os << "\t-a <estimator> \t MI estimator type\n";
	os << "\t-b <bins> \t number of bins for MI estimation\n";
	os << "\t-k <order>\t B-spline order for MI estimation\n";
	os << "\t-p <pval>\t MI p-value\n";
	os << "\t-t <threshold>\t MI threshold value\n";
	os << "\t-e <tolerance>\t DPI tolerance\n";
	os << "\t-C <method>\t convert MI\n";
	os << "\t-r <size>\t generate size bootstrap networks\n";
	os << "\t-y <size>\t use sample of size columns\n";
	os << "\t-x <seed>\t random number generator seed\n";
	os << "\t-m        \t turn on memory reports\n";
	os << "\t-v        \t show progress\n";
	os << "\t-h        \t print this help\n";
	os << std::endl;
    } // print_help

    void print_help_ext(std::ostream& os) {
	print_help(os);
	os << "Additional options:\n\n";
	os << "\t-W <file> \t write MCL label output to file\n";
	os << "\t-q <size> \t permutation test size\n";
	os << std::endl;
    } // print_help_ext

    void print_config(std::ostream& os) const {
	const int N = 32;
	os.flags(std::ios_base::left);
	os << std::setfill('.');
	os << std::setw(N) << "> Size of double: " << ' ' << sizeof(double) << "\n";
	os << std::setw(N) << "> Input file: " << ' ' << input << "\n";
	os << std::setw(N) << "> Output file: " << ' ' << output << "\n";
	os << std::setw(N) << "> MI network input: " << ' ';
	if (input_mi.empty() == true) os << "none"; else os << input_mi;
	os << "\n";
	os << std::setw(N) << "> MI network output: " << ' ';
	if (output_mi.empty() == true) os << "none"; os << output_mi;
	os << "\n";
	os << std::setw(N) << "> TFs list: " << ' ';
	if (input_tf.empty() == true) os << "none"; else os << input_tf;
	os << "\n";
	os << std::setw(N) << "> MI estimator: " << ' ' << estimator2str(mi_estim) << "\n";
	os << std::setw(N) << "> MI bins number: " << ' ' << mi_b << "\n";
	os << std::setw(N) << "> MI B-spline order: " << ' ' << mi_k << "\n";
	os << std::setw(N) << "> MI p-value: " << ' ' << mi_pval << "\n";
	os << std::setw(N) << "> MI threshold: " << ' ' << mi_thr << "\n";
	os << std::setw(N) << "> DPI tolerance: " << ' ' << dpi_tol << "\n";
	os << std::setw(N) << "> MI conversion: " << ' ' << convert2str(mi_conv) << "\n";
	os << std::setw(N) << "> Bootstrap size: " << ' ' << boot_size << "\n";
	os << std::setw(N) << "> Column sample: " << ' ' << sample_size << "\n";
	os << std::setw(N) << "> RNG seed: " << ' ' << rng_seed << "\n";
    } // print_config


    bool parse(int aargc, char* aargv[], std::ostream& os, int print_err = 1) {
	argc = aargc;
	argv = aargv;

	int c = 0;

	extern char* optarg;
	extern int opterr;

	opterr = print_err;

	jaz::uc_compare ucmp;

	bool has_output_mi = false;
	bool has_output_lab = false;
	bool has_mi_pval = false;
	bool has_mi_thr = false;

	// we use gnu extension (double colon)
	const char* params = "i:o:j:w:l:W:a:b:k:p:t:e:C:y:x:r:q:mvhH";

	while ((c = getopt(argc, argv, params)) != -1) {
	    switch (c) {
	      case 'i':
		  input = optarg;
		  break;

	      case 'o':
		  output = optarg;
		  break;

	      case 'j':
		  input_mi = optarg;
		  break;

	      case 'w':
		  output_mi = optarg;
		  has_output_mi = true;
		  break;

	      case 'l':
		  input_tf = optarg;
		  break;

	      case 'W':
		  output_lab = optarg;
		  has_output_lab = true;
		  break;

	      case 'a':
		  switch (optarg[0]) {
		    case 'B':
		    case 'b':
			mi_estim = BSPLINE;
		    break;

		    case 'G':
		    case 'g':
			mi_estim = GAUSSIAN;
		    break;

		    default:
			os << error << "unknown MI estimator\n";
			return false;
		  }
		  break;

	      case 'b':
		  mi_b = std::atoi(optarg);
		  break;

	      case 'k':
		  mi_k = std::atoi(optarg);
		  break;

	      case 'p':
		  mi_pval = std::atof(optarg);
		  has_mi_pval = true;
		  break;

	      case 't':
		  mi_thr = std::atof(optarg);
		  has_mi_thr = true;
		  break;

	      case 'e':
		  dpi_tol = std::atof(optarg);
		  break;

	      case 'C':
		  switch (optarg[0]) {
		    case 'C':
		    case 'c':
			mi_conv = CORRELATION;
		    break;

		    default:
			os << error << "unknown MI conversion mode\n";
			return false;
		  }
		  break;

	      case 'r':
		  boot_size = std::atoi(optarg);
		  break;

	      case 'y':
		  sample_size = std::atoi(optarg);
		  break;

	      case 'x':
		  rng_seed = std::atoi(optarg);
		  break;

	      case 'q':
		  ptest_size = std::atoi(optarg);
		  break;

	      case 'm':
		  mem_report = true;
		  break;

	      case 'v':
		  verbose = true;
		  break;

	      case 'h' :
		  if (print_err == 1) print_help(os);
		  return false;

	      case 'H' :
		  if (print_err == 1) print_help_ext(os);
		  return false;

	      default:
		  os << error << "wrong parameters\n";
		  return false;
	    } // switch
	} // while

	// check input file names
	if (input.empty() == true) {
	    os << error << "no input file specified\n";
	    return false;
	}

	if (output.empty() == true) {
	    std::string::size_type pos = input.rfind('.');
	    if (pos != std::string::npos) {
		output = std::string(input.begin(), input.begin() + pos);
	    } else output = input;
	    output += ".adj";
	}

	if (has_output_mi == true) {
	    if ((output_mi.empty() == true) || (output_mi == "-")) {
		std::string::size_type pos = input.rfind('.');
		if (pos != std::string::npos) {
		    output_mi = std::string(input.begin(), input.begin() + pos);
		} else output_mi = input;
		output_mi += ".mi.adj";
	    }
	}

	if (has_output_lab == true) {
	    if ((output_lab.empty() == true) || (output_lab == "-")) {
		std::string::size_type pos = input.rfind('.');
		if (pos != std::string::npos) {
		    output_lab = std::string(input.begin(), input.begin() + pos);
		} else output_lab = input;
		output_lab += ".lab";
	    }
	}

	// check estimator

	// check MI estimation parameters
	if (mi_estim == BSPLINE) {
	    if ((mi_b <= mi_k) || (mi_k < 2) || (mi_b < 4)) {
		os << error << "wrong parameters of MI estimator\n";
		return false;
	    }
	}

	if (has_mi_pval && has_mi_thr) {
	    os << error << "-p and -t are exclusive, choose one\n";
	    return false;
	}

	if ((input_mi.empty() == false) && (has_mi_pval == true)) {
	    os << error << "MI threshold must be provided\n";
	    return false;
	}

	if ((mi_pval < std::numeric_limits<double>::epsilon()) || (mi_pval > 1.0)) {
	    os << error << "MI p-value must be in range (0,1]\n";
	    return false;
	}

	if (mi_thr < 0) {
	    os << error << "MI threshold less than 0\n";
	    return false;
	}

	if ((dpi_tol < std::numeric_limits<double>::epsilon()) || (dpi_tol > 1)) {
	    os << error << "DPI tolerance must be in range (0,1]\n";
	    return false;
	}

	// check sample
	if ((sample_size != 0) && (sample_size < mi_b)) {
	    os << error << "column sample too small\n";
	    return false;
	}

	if ((boot_size > 0) && (sample_size > 0)) {
	    os << error << "-r and -y are exclusive, choose one\n";
	    return false;
	}

	// check seed
	if (rng_seed < 0) rng_seed = std::time(0) + getpid();

	// ptest size
	if (ptest_size < 1) {
	    os << error << "permutation test size less than 1\n";
	    return false;
	}

	return true;
    } // parse


    int argc;
    char** argv;

    // STL check, true if STL exploits pubsetbuf
    bool has_pubsetbuf;

    // name of the input file (with genes expression data)
    std::string input;

    // name of the adj file to create
    std::string output;

    // name of the adj file to store MI relevance network
    std::string output_mi;

    // name of adj file from which adj should be taken
    std::string input_mi;

    // name of the file with TF list
    std::string input_tf;

    // name of the lab file for MCL
    std::string output_lab;

    // MI estimator type
    estimator_t mi_estim;

    // number of bins of MI estimator
    unsigned int mi_b;

    // number of basis functions for MI estimator
    unsigned int mi_k;

    // p-value for MI
    double mi_pval;

    // MI threshold
    double mi_thr;

    // DPI tolerance
    double dpi_tol;

    // convert MI
    convert_t mi_conv;

    // size of bootstrap
    unsigned short int boot_size;

    // column sample size
    unsigned int sample_size;

    // seed for random number generator
    long int rng_seed;

    // permutation test size
    long int ptest_size;

    // memory reports on/off
    bool mem_report;

    // verbosity in main loop
    bool verbose;

}; // struct AppConfig

#endif // APP_CONFIG_HPP
