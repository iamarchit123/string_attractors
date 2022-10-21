/**
 * @file    main.cpp
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
    const std::uint64_t max_text_length) {
  double tot_time = 0.0;
  // Print initial message.
  fprintf(stderr, "TEST, max_length = %lu\n",
      max_text_length);
  
  // Allocate text and SA.
  char_type * const text = new char_type[max_text_length];
  //Generate text
  for (std::uint64_t i = 0; i < max_text_length; ++i)
    text[i] = 'a' + utils::random_int<std::uint64_t>(0UL, 4);
  st_att<> * st_att_file = new st_att<>(2, text,  max_text_length);

  text_offset_type * const indexes = new text_offset_type[50000];
  for (std::uint64_t i = 0; i < 50000; ++i)
    indexes[i] = utils::random_int<std::uint64_t>(0UL, max_text_length-1);


  // Run tests.
  for (std::uint64_t testid = 0; testid < 20; ++testid) {  
    for(text_offset_type index=0; index< 50000; index++){
      double t1 = utils::wclock();
      char alpha  = st_att_file->query(indexes[index]);
      tot_time+=(utils::wclock()-t1);
      if(alpha != text[indexes[index]])
        fprintf(stderr,"Alphabet didnt match\n");
    }
  }
  delete(st_att_file);
  delete[] text;
  delete[] indexes;
  return tot_time;
}

int main() {

  // Init random number generator.
  srand(time(0) + getpid());

  static const std::uint64_t text_length_limit = (1 << 24);
  vector<double> time_to_construct(25,0);
  typedef std::uint8_t char_type;
  typedef long unsigned int text_offset_type;

  // Run tests.
  int ind = 0;
  for (std::uint64_t max_text_length = 1;
      max_text_length <= text_length_limit; max_text_length *= 2)
    time_to_construct[ind++] = test<char_type, text_offset_type>(max_text_length);

  // Print summary.
  for(long int i=0; i<=24;i++)
    fprintf(stderr,"Query time for string of size %ld is : %f\n",(long)1 << i , time_to_construct[i]);
}


