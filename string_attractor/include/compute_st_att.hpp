/**
 * @file    compute_sa.hpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
 *
 * Copyright (C) 2017-2021
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

#ifndef __COMPUTE_ST_ATT_HPP_INCLUDED
#define __COMPUTE_ST_ATT_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <math.h>
#include "uint40.hpp"
#include "uint48.hpp"
#include "sais.hxx"
#include "naive_compute_sa.hpp"
#include "compute_lz77.hpp"
#include "compute_sa.hpp"
#include "utils.hpp"
#include <cstring>

char *sliceString(unsigned char *str, int start, int end)
{

    int i;
    int size = (end - start) + 2;
    char *output = (char *)malloc(size * sizeof(char));

    for (i = 0 ; i < size - 1; i++)
    {
        output[i] = str[start + i];
    }

    output[size-1] = '\0';

    return output;
}

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
struct linked_indexes{
  typedef std::pair<text_offset_type, text_offset_type> pair_type;
  pair_type p;
  char* s;
  uint32_t level;
  
  linked_indexes(){
    s = NULL;
    p = std::make_pair<text_offset_type, text_offset_type>(UINT64_MAX,UINT64_MAX);
  }

  linked_indexes(uint32_t start, int attractor_loc_,int level_){
    level = level_;
    p = std::make_pair<text_offset_type, text_offset_type>(start, attractor_loc_);
    s = NULL;
  } 

  linked_indexes(uint32_t level_, char *s_){
    level = level_;
    s = strdup(s_);
    p = std::make_pair<text_offset_type, text_offset_type>(UINT64_MAX,UINT64_MAX);
  }

  void print(){
    fprintf(stderr,"offset_start = %lu, attractor start index = %lu, level = %u, string is %s\n"
    , p.first, p.second, level, s);
  }

};

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
int find_nearest(std::vector<text_offset_type>att_pos, std::uint64_t ind){
  uint64_t ind1 = 0;
  while(1){
    if(att_pos[ind1] >= ind)
      return ind1;
    ind1++;
  }
}

//To implement RMQ with binary search on top of SA later 
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
int find(char_type * text,int start,
 int len, std::uint32_t n){
   
    //fprintf(stderr,"Trying to find %s starting at index %d\n",text_to_find,start);
    for(std::uint32_t i = 0; i<n;i++){
      int j;
      for(j = 0; j< len;j++){
        if(text[i + j]!= text[start + j])
          break;
      }
      if(j==len)
        return i;
    }
   return -1;
}

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
class st_att{
  private:
    std::uint32_t tau;
    std::uint32_t gamma;
    std::uint32_t alpha;
    std::uint32_t n;
    std::vector<text_offset_type>att_pos;
    std::vector<std::vector<linked_indexes<>>>indexes;

  public:
    st_att(std::uint32_t m_tau, char_type * text, std::uint64_t text_length){
      std::uint32_t dim,temp, level;
      typedef std::pair<text_offset_type, text_offset_type> pair_type;
      tau = m_tau;
      n = text_length;
      // Allocate SA.
      text_offset_type * const sa = new text_offset_type[text_length];


      // Compute SA.
      {
        fprintf(stderr, "Compute SA... ");
        const long double start = utils::wclock();
        compute_sa(text, n, sa);
        const long double elapsed = utils::wclock() - start;
        fprintf(stderr, "%.2Lfs\n", elapsed);
      }
      
      // Compute parsing.
      std::vector<pair_type> parsing;
      {
        fprintf(stderr, "Compute LZ77... ");
        const long double start = utils::wclock();
        compute_lz77::kkp2n(text, text_length, sa, parsing);
        const long double elapsed = utils::wclock() - start;
        fprintf(stderr, "%.2Lfs\n", elapsed);
      }

      gamma = parsing.size();
      //TODO ::Discuss with prof or think about it how to set alpha
      alpha = 5;

      //Prepare att_pos
      att_pos.push_back(0);
      std::uint64_t ind = 0;
      for(uint32_t i=1; i < gamma; i++) {
        ind += parsing[i].second ? parsing[i].second : 1;
        att_pos.push_back(ind);
      }

      for(uint32_t i=0; i < gamma; i++)
        fprintf(stderr,"paramterer 1 = %lu, parameter 2 = %lu, att_pos = %lu\n", parsing[i].first, parsing[i].second, att_pos[i]);
      fprintf(stderr, "gamma = %u\n",gamma);
      fprintf(stderr, "tau = %u\n", tau);
      fprintf(stderr, "alpha = %u\n", alpha);
      fprintf(stderr, "Text length = %u\n", n);

      dim = n/gamma;
      temp = n/gamma;
      
      level = 0;
      
      //Make vector of vectors fill by level
      while(1){
        fprintf(stderr,"Calculating Level %d\n",level);
        if(!level){
          indexes.resize(level + 1);
          //Make first layer of n/gamma pointers ready
          for( std::uint32_t i = 0; i < n; i+=dim ){
            uint32_t len = std::min(dim,n-i);
            fprintf(stderr,"Before find\n");        
            ind = find(text, i, len,n);
            fprintf(stderr,"After find ind = %lu\n",ind);
            if(temp > alpha)
                indexes[0].push_back(linked_indexes<>(ind,find_nearest(att_pos,ind),0));
              else
                indexes[0].push_back(linked_indexes<>(0, sliceString(text,ind,ind + len -1)));
            fprintf(stderr,"After initializing vector for i = %u\n",i);
          }
        }
        else{
          indexes.resize(level + 1);
          indexes[level].resize((4*tau-1)*gamma);
          int skipper = 4*tau - 1;
          //make  (4*tau-1) * gamma pointers ready
          for(uint32_t i=0; i < gamma; i++) {
            int start = (int)att_pos[i] - (tau*temp) + 1, end;
            int tempstart = start;
            fprintf(stderr,"Filling index %d with start = %d\n",i,start);
            for(uint32_t j=1; j <= 2*tau; j++) { 
              if ((tempstart + (int)temp) <= 0 || tempstart >= (int)n){
                indexes[level][i*skipper + j-1] = linked_indexes<>();
                fprintf(stderr,"Hello from j = %d\n ",j);
                tempstart += temp;  
                continue;
              }  
              end = std::min(tempstart + temp - 1, n - 1);
              tempstart = std::max(tempstart, 0);
              ind = find(text, (uint32_t)tempstart, (end - tempstart + 1), n);
              
              if(temp > alpha)
                indexes[level][i*skipper + j-1] = linked_indexes<>(ind,find_nearest(att_pos,ind),level);
              else{
                fprintf(stderr,"Filling index %d with start =%d end =%d\n ",i*skipper+j-1,tempstart,end);
                indexes[level][i*skipper + j-1] = linked_indexes<>(level, sliceString(text,tempstart,end));
              }
              tempstart = end + 1;
          }
          start = (int)att_pos[i] - (tau*temp) + 1 + temp/2;
          tempstart = start;
          for(uint32_t j=1; j <= 2*tau - 1; j++){
              
              if ((tempstart + (int)temp) <= 0 || tempstart >= (int)n){
                indexes[level][i*skipper + 2* tau + j-1] = linked_indexes<>();
                fprintf(stderr,"Hello from j = %d\n ",j);
                tempstart+=temp;
                continue;
              }
              end = std::min(tempstart + temp - 1, n - 1);
              tempstart = std::max(tempstart, 0);
              ind = find(text, (uint32_t)tempstart, (end - tempstart + 1), n);
              if(temp > alpha)
                indexes[level][i*skipper + 2* tau + j-1] = linked_indexes<>(ind,find_nearest(att_pos,ind),level);
              else{
                fprintf(stderr,"Filling index %d with start =%d end =%d\n ",i*skipper + 2* tau + j-1,tempstart,end);
                indexes[level][i*skipper + 2* tau + j-1] = linked_indexes<>(level, sliceString(text,tempstart,end));
              }
              tempstart=end + 1;
            }
          }
        }

        fprintf(stderr,"Level %d calculated\n",level);
        if(temp <= alpha){
          alpha = temp;
          break;
        }
        
        temp/=tau;
        level++;
        
      }

      fprintf(stderr,"No of levels is %lu\n",indexes.size());
      for(std::uint32_t i=0; i < indexes.size(); i++) {
        fprintf(stderr,"Level %u size is %lu\n", i, indexes[i].size());
        for(std::uint32_t j = 0; j < indexes[i].size(); j++){
            fprintf(stderr,"index %u is ", j);
            indexes[i][j].print();
        }
      }
  }

  char * query_alpha(int start, int p_len){
    int level = 0;
    int start_cell = start/(n/gamma);
    int skipper = 4*tau - 1;
    int temp = n/gamma;
    int start_loc = start/temp * temp;
    while( indexes[level][start_cell].s==NULL ){
      temp/=tau;
      start = indexes[level][start_cell].p.first + (start - start_loc);
      int att_loc = indexes[level][start_cell].p.second;
      int start1 = (int)att_pos[att_loc] - ( tau * temp ) + 1;
      int start2 = start1 + temp/2;
      bool done = false;
      
      

      for(uint32_t i=0; i < (2*tau - 1); i++ ){
        fprintf(stderr,"start = %d start1 = %d start2 = %d\n",start, start1, start2);
        if(start >= start1 && (start + p_len) <= start1 + temp){
          start_cell = att_loc * skipper + i;
          start_loc  = start1;
          done = true;
          break;
         }
        else if(start >=start2 && (start + p_len) <= start2 + temp){
          start_cell = att_loc * skipper + 2* tau + i;
          start_loc = start2;
          done = true;
          break;
        }

        start1 += temp;
        start2 += temp;
      }

      if(!done){
        start_cell = att_loc * skipper + 2*tau - 1;
        start_loc = start1;
      }
      if(start_loc < 0)
        start_loc = 0;
      level++;
      fprintf(stderr, "level = %d start_cell = %d\n",level,start_cell);
    }
    return indexes[level][start_cell].s + start - start_loc;
  }

  //Store expilcit strings here
  char * query(int start, int len){
    int end = start + len -1;   
    int si = len + 1;
    char *output = (char *)malloc(si * sizeof(char));
    int tempend,templen;
    //fprintf(stderr,"alpha is %u\n", alpha);
    for(int i=start;i <= end; i+=alpha){
      tempend = i+ alpha - 1;
      tempend = std::min(tempend,end);
      //fprintf(stderr,"1 start = %d end =%d\n",i,tempend);
      if((i/(n/gamma)) != (tempend/(n/gamma))){
        //fprintf(stderr, "WTF\n");
        tempend = (i/(n/gamma) + 1) * (n/gamma) - 1;
      }
      tempend = std::min(tempend,end);
      fprintf(stderr,"start =%d end =%d\n",i,tempend);
      templen = tempend - i + 1;
      strncpy(output + (i - start), query_alpha(i,templen),templen);
      i += (templen - alpha);
      //fprintf(stderr,"i=%d\n",i);

    }
    output[si-1] = '\0';
    return output;
  }

};

#endif  // __COMPUTE_ST_ATT_HPP_INCLUDED