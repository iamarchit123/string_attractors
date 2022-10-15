/**
 * @file    construction.cpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
 *
 * Copyright (C) 2021
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>
#include <unistd.h>

#include "../include/utils.hpp"
#include "../include/compute_sa.hpp"
#include "../include/compute_lz77.hpp"
#include "../include/uint40.hpp"
#include "../include/compute_st_att.hpp"

//=============================================================================
// Run tests for text up to a given length on a given number of cases.
//=============================================================================
template<typename char_type, typename text_offset_type>
double test(
    const std::uint64_t max_text_length,
    const std::uint64_t testcases) {
  double tot_time = 0.0;
  // Print initial message.
  fprintf(stderr, "TEST, max_length = %lu, testtimes = %lu\n",
      max_text_length, testcases);

  
  // Allocate text and SA.
  char_type * const text = (char_type *)calloc(max_text_length,sizeof(char_type));
  if(!text){
    fprintf(stderr,"Archit couldn't allocate memory for %ld\n",max_text_length);
    exit(-1);
  }
  else
    fprintf(stderr,"Archit allocated memory for %ld\n",max_text_length);
  // Run tests.
  for (std::uint64_t testid = 0; testid < testcases; ++testid) {

    // Print progress message.
    if (testid % 10 == 0)
      fprintf(stderr, "%.2Lf%%\r", (100.L * testid) / testcases);

    // Generate the text.
    const std::uint64_t text_length = max_text_length;
    //  utils::random_int<std::uint64_t>(
    //      (std::uint64_t)1,
    //      (std::uint64_t)max_text_length);
    fprintf(stderr,"Archit trying to work with text length %ld\n",text_length);
    for (std::uint64_t i = 0; i < text_length; ++i)
      text[i] = 'a' + utils::random_int<std::uint64_t>(0UL, 4);
    fprintf(stderr,"Archit wrote letters with text length %ld\n",text_length);
    /* Compute string attractor structure*/
    double t1 = utils::wclock();
    st_att<> * st_att_file = new st_att<>(2, text,  text_length);
    tot_time+=(utils::wclock()-t1);
    fprintf(stderr,"Archit string attractor work complete: %ld\n",text_length);
    //fprintf(stderr,"Tot time is %f",tot_time);
    delete(st_att_file);
  }
  delete[] text;
  return tot_time/testcases;
}

int main() {

  // Init random number generator.
  srand(time(0) + getpid());
  vector<double> time_to_construct(21,0);
  static const std::uint64_t text_length_limit = (1 << 17);
  static const std::uint64_t n_tests = 5;

  typedef std::uint8_t char_type;
  typedef long unsigned int text_offset_type;

  // Run tests.
  int ind=0;
  for (std::uint64_t max_text_length = 1;
      max_text_length <= text_length_limit; max_text_length *= 2)
    time_to_construct[ind++] = test<char_type, text_offset_type>(max_text_length, n_tests);

  // Print summary.
  for(long int i=0; i<=20;i++)
    fprintf(stderr,"Size of string %ld Time to construct : %f\n",(long)1 << i , time_to_construct[i]);
}


