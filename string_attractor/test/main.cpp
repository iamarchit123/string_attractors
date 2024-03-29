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
void test(
    const std::uint64_t max_text_length,
    const std::uint64_t testcases) {

  // Print initial message.
  fprintf(stderr, "TEST, max_length = %lu, testcases = %lu\n",
      max_text_length, testcases);

  
  // Allocate text and SA.
  char_type * const text = new char_type[max_text_length];

  // Run tests.
  for (std::uint64_t testid = 0; testid < testcases; ++testid) {

    // Print progress message.
    if (testid % 10 == 0)
      fprintf(stderr, "%.2Lf%%\r", (100.L * testid) / testcases);

    // Generate the text.
    const std::uint64_t text_length =
      utils::random_int<std::uint64_t>(
          (std::uint64_t)1,
          (std::uint64_t)max_text_length);
    for (std::uint64_t i = 0; i < text_length; ++i)
      text[i] = 'a' + utils::random_int<std::uint64_t>(0UL, 4);
    
    /* Compute string attractor structure*/
    st_att<> * st_att_file = new st_att<>(2, text,  text_length);
    for(text_offset_type index=0; index< text_length; index++){
      char alpha = st_att_file->query(index);
      if(text[index]!=alpha){
        fprintf(stderr, "\nError:\n");
        fprintf(stderr, "  text_length = %lu\n", text_length);
        fprintf(stderr, "Wrong at index %lu\n",index);
        fprintf(stderr, "Actual alphabet %c Predicted Alphabet %c\n",text[index], alpha);
        fprintf(stderr, "  text = ");
        for (std::uint64_t i = 0; i < text_length; ++i)
          fprintf(stderr, "%c", text[index]);
        fprintf(stderr, "\n");
        
      std::exit(EXIT_FAILURE);
      }
    }
    delete(st_att_file);
  }
  delete[] text;
}

int main() {

  // Init random number generator.
  srand(time(0) + getpid());

  static const std::uint64_t text_length_limit = (1 << 20);
  static const std::uint64_t n_tests = 1000;

  typedef std::uint8_t char_type;
  typedef long unsigned int text_offset_type;

  // Run tests.
  for (std::uint64_t max_text_length = 1;
      max_text_length <= text_length_limit; max_text_length *= 2)
    test<char_type, text_offset_type>(max_text_length, n_tests);

  // Print summary.
  fprintf(stderr, "All tests passed.\n");
}


