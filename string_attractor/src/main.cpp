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
#include <getopt.h>

#include "../include/utils.hpp"
#include "../include/compute_sa.hpp"
#include "../include/compute_lz77.hpp"
#include "../include/uint40.hpp"
#include "../include/compute_st_att.hpp"

std::uint8_t *generate_string(int alphabets, int length) {
        std::uint8_t * const text = new std::uint8_t[length];
        for(int i=0;i<length;i++)
                text[i] = rand()%alphabets + 'a';
        for(int i=0;i<length;i++)
                fprintf(stderr,"%c",text[i]);
        fprintf(stderr,"\n");
        return text;

}

//=============================================================================
// Computes the string atractor parsing of the file given
// as an argument and write to output file. Then converts it to string attractor
// and related sturctures
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::int64_t>
void test_st_att(void) 
  {
    //const char *z = "CDABCCDABCCA";
    //std::uint64_t text_length = 12;
    //char_type *text;
    //text = generate_string(z, text_length);
    int alphabet_size = 26, testcases  = 10,text_size = 20, text_length, alphabets,temp;
    int start ,end;
    char_type *text;

      text_length = 20000;
      alphabets = 2;
      text = generate_string(alphabets, text_length);
    while(testcases--) {
      start = rand()%text_length;
      end = rand()%text_length;
      if(start > end){
        temp = end;
        end = start;
        start = temp;
      }
      st_att<> * st_att_file = new st_att<>(2, text,  text_length);
      char *x = st_att_file->query(start,end - start + 1);
      for(int i=start; i<= end; i++)
      {
        if(x[i - start]!= text[i])
        {
          fprintf(stderr,"Testcase %d failed for query %d %d %d %c %c\n",testcases, start , end , i, x[i - start],text[i]);
          exit(-1);
        }
      }
      fprintf(stderr,"Testcase %d passed\n",testcases);
      
    }

    
  }

//=============================================================================
// Print usage instructions and exit.
//=============================================================================
void usage(
    const char * const program_name,
    const int status) {
  printf(

"Usage: %s [OPTION]... FILE\n"
"Construct the string attractor and data structures parsing of text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -h, --help              display this help and exit\n"
"  -o, --output=OUTFILE    specify output filename. Default: FILE.st_att\n",

    program_name);

  std::exit(status);
}

int main(int argc, char **argv) {

  // Initial setup.
  srand(time(0) + getpid());
  const char * const program_name = argv[0];

  // Declare flags.
  static struct option long_options[] = {
    {"help",     no_argument,       NULL, 'h'},
    {"output",   required_argument, NULL, 'o'},
    {NULL,       0,                 NULL, 0}
  };

  // Initialize output filename.
  std::string output_filename("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "ho:",
          long_options, NULL)) != -1) {
    switch(c) {
      case 'h':
        usage(program_name, EXIT_FAILURE);
        break;
      case 'o':
        output_filename = std::string(optarg);
        break;
      default:
        usage(program_name, EXIT_FAILURE);
        break;
    }
  }

  // Print error if there is not file.
  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(program_name, EXIT_FAILURE);
  }

  // Parse the text filename.
  const std::string text_filename = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  // Set default output filename (if not provided).
  if (output_filename.empty())
    output_filename = text_filename + ".st_att";

  // Check for the existence of text.
  if (!utils::file_exists(text_filename)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(program_name, EXIT_FAILURE);
  }

  // Check if output file exists.
  if (utils::file_exists(output_filename)) {

    // Output file exists, should we proceed?
    char *line = NULL;
    std::uint64_t buflen = 0;
    std::int64_t len = 0L;

    // Obtain the answer.
    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          output_filename.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        printf("\nError: failed to read answer\n\n");
        std::fflush(stdout);
        usage(program_name, EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    // If not, then exit.
    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }

    // Otherwise, we proceed.
    free(line);
  // Set types.
  }

  //typedef std::uint8_t char_type;
  //typedef uint40 text_offset_type;

  // Run the lz77 parsing algorithm.
  test_st_att<>();

}
