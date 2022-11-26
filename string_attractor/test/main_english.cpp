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

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
void test_conversion(
    std::string parsing_filename) {

  // Read file
  double tot_time = 0.0;
  FILE *f = fopen(parsing_filename.c_str(), "rb");
  fseek(f, 0, SEEK_END);
  const std::uint64_t fsize = ftell(f);
  fseek(f, 0, SEEK_SET);
  char_type * const text = (char_type *)malloc(fsize);
  size_t result = fread(text, 1, fsize, f);
  if (result != fsize) {fputs ("Reading error",stderr); exit (3);}
  fclose(f);
  // Run tests.
  /* Compute string attractor structure*/
 /* for(int i=0; i<10;i++){
    fprintf(stderr,"1\n");
    double t1 = utils::wclock();
    st_att<> * st_att_file = new st_att<>(2, text,  fsize);
    tot_time+=(utils::wclock()-t1);
    fprintf(stderr,"10\n");
    delete(st_att_file);
  }*/
  
  st_att<> * st_att_file = new st_att<>(2, text,  fsize);
  fprintf(stderr,"Size of string %ld Time to construct : %fs\n",fsize , tot_time);
  text_offset_type * const indexes = new text_offset_type[50000];
  for (std::uint64_t i = 0; i < 50000; ++i)
    indexes[i] = utils::random_int<std::uint64_t>(0UL, fsize-1);
  for (std::uint64_t testid = 0; testid < 1; ++testid) {  
    int count = 0;
    for(text_offset_type index=0; index< 50000; index++){
      double t1 = utils::wclock();
      char_type alpha = st_att_file->query(indexes[index]);
      if(alpha!=text[indexes[index]]){
        fprintf(stderr,"actaul alphabet = %d predicted = %d\n",text[indexes[index]], alpha);
        count++;
      }
      tot_time+=(utils::wclock()-t1);
    }
    fprintf(stderr,"count of misclassified = %d\n",count);
  }
  delete(st_att_file);
  delete[] text;
  delete[] indexes;
  fprintf(stderr,"time taken for 1 million queries:%fs\n",tot_time);
}

int main(int argc, char **argv) {
  if (argc != 2)
    std::exit(EXIT_FAILURE);

  // Initialize runtime statistics.
  utils::initialize_stats();

  // Declare types.
  typedef std::uint8_t char_type;

  // Obtain filenames.
  std::string parsing_filename = argv[1];

  // Run the algorithm.
  test_conversion<char_type, long unsigned int>(
      parsing_filename);
}

