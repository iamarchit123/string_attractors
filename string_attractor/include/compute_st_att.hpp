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
#include <map>
using namespace std;

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::int64_t>
struct Block {
  text_offset_type start;
  text_offset_type len;
  text_offset_type end;

  Block(text_offset_type m_start, text_offset_type m_len)
  {
    start = m_start;
    len = m_len;
    end = m_start + m_len;
  }

};

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::int64_t>
struct linked_indexes{
  typedef std::pair<text_offset_type, text_offset_type> pair_type;
  pair_type p;
  
  linked_indexes(char_type* text, text_offset_type start,
                 text_offset_type end, std::vector<text_offset_type>&att_pos,
                 text_offset_type len)
  {
    text_offset_type attractor_loc_, offset;
    find(text, start, end, att_pos, &attractor_loc_, &offset, len);
    p = std::make_pair<text_offset_type, text_offset_type>( (long int)attractor_loc_, 
                                                            (long int)offset);
  }

  long unsigned int start()
  {
    return p.first - p.second;
  }

 static void find(char_type * text,text_offset_type start, text_offset_type end,
           std::vector<text_offset_type>&att_pos,text_offset_type *attractor_loc_, 
           text_offset_type *offset, text_offset_type len) {
   
    //This is very slow !! need to implement RMQ with binary search later on toop of it
    *attractor_loc_ = -1;
    *offset = -1;
    end = min(len-1,end);
    for (text_offset_type i = 0; i<len ; i++)
    {

      text_offset_type j;
      for(j = start; j <=end && j < len ;j++) {
        if (text[i + j - start] != text[j])
          break;
      }
      if(j==end+1 || j==len)
      {
        for(text_offset_type k = 0 ; k < att_pos.size(); k++)
        {
          if(att_pos[k]>=i && att_pos[k] <= (i + j - start))
          {
            *attractor_loc_ = att_pos[k];
            *offset = att_pos[k] - i;
            return;
          }
        }
        fprintf(stderr,"Couldn't find a match the program may not work ahead\n");
          exit(-1);
      }
    }
  }

  
};

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::int64_t>
std::vector<linked_indexes<>> make_linked_indexes(char_type *text,std::vector<Block<>>&b,
                                        std::vector<text_offset_type>&att_pos,text_offset_type len)
{
  std::vector<linked_indexes<>>v;
  for(text_offset_type i =0 ; i < b.size(); i++)
  {
    v.push_back(*(new linked_indexes<>(text,max((text_offset_type)0,b[i].start),b[i].end,att_pos,len)));
  }
  return v;
}

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::int64_t>
class st_att{
  private:
    text_offset_type tau;
    text_offset_type gamma;
    text_offset_type alpha;
    text_offset_type n;
    std::vector<text_offset_type>att_pos;
    std::vector<text_offset_type>b_si;
    std::vector<std::vector<linked_indexes<>>>indexes;
    std::map<text_offset_type,text_offset_type>att_map;
    std::vector<char *>v_s;

  public:
    st_att(text_offset_type m_tau, char_type * text, text_offset_type text_length){
      typedef std::pair<text_offset_type, text_offset_type> pair_type;
      n = text_length;
      text_offset_type block_len;
      tau = m_tau;
      std::vector<Block<>> v;
      // Allocate SA.
      text_offset_type * const sa = new text_offset_type[n];


      // Compute SA.
      {
        fprintf(stderr, "Compute SA... \n");
        const long double start = utils::wclock();
        compute_sa(text, n, sa);
        const long double elapsed = utils::wclock() - start;
        fprintf(stderr, "%.2Lfs\n", elapsed);
      }
      
      // Compute parsing.
      std::vector<pair_type> parsing;
      {
        fprintf(stderr, "Compute LZ77... \n");
        const long double start = utils::wclock();
        compute_lz77::kkp2n(text, text_length, sa, parsing);
        const long double elapsed = utils::wclock() - start;
        fprintf(stderr, "%.2Lfs\n", elapsed);
      }

      //TODO ::Discuss with prof or think about it how to set alpha

      

      //Prepare att_pos
      fprintf(stderr, "Compute LZ77... \n");
      std::uint64_t ind = -1;
      for(uint32_t i=0; i < parsing.size(); i++) {
        ind += parsing[i].second ? parsing[i].second : 1;
        att_map[ind] = i;
        att_pos.push_back(ind);
      }
      gamma = att_pos.size();
      //Make level 0 and assign alpha 
      block_len = n / gamma + (n % gamma !=0);
      for(text_offset_type i = 0; i < n; i+=block_len)
          v.push_back(Block<>(i, block_len));
      b_si.push_back(block_len);
      indexes.push_back(make_linked_indexes<>(text,v,att_pos,n));
      alpha = max((int)ceil(log(block_len)/log(tau)),1);
      //Now make all other levels
      while(block_len >= 2 * alpha)
      {
        block_len = block_len / tau + (block_len % tau !=0);
        b_si.push_back(block_len);
        v.clear();
        for(text_offset_type i=0; i < att_pos.size(); i++)
        {
          text_offset_type begin = att_pos[i] - tau * block_len;
          text_offset_type end = att_pos[i] + tau * block_len;
          for( text_offset_type j=begin; j < end; j+=block_len)
            v.push_back(Block<>(j,block_len));
          indexes.push_back(make_linked_indexes<>(text,v,att_pos,n));
        }
      }
      //Code for appending strings from v
      for(text_offset_type i =0 ; i < v.size(); i++)
      {
        Block<> t = v[i];
        char *s = new char[t.len];
        for(text_offset_type j = t.start; j < t.end; j++)
        {
          if(j>=0 && j<n)
            s[j - t.start] = text[j];
          else
          s[j - t.start] = '$';
        }
        v_s.push_back(s);
      }
      v.clear();
  }

  char query_alpha(text_offset_type off, uint32_t level, text_offset_type attractor)
  {
    text_offset_type block_position, offset, block_len = b_si[level];
    if (level == 0)
    {
      block_position = off/block_len;
      offset = off % block_len;
    }
    else
    { 
      block_position = att_map[attractor] * tau * 2 + tau + (off - attractor) / block_len;
      offset = (off - attractor) % block_len;
    }
    if (level == b_si.size() - 1)
      return v_s[block_position][offset];
    
    linked_indexes<> l = indexes[level][block_position];
      
    return query_alpha(
      l.start() + offset,
      level + 1,
      l.p.first
    );
  }

  //Store expilcit strings here
  char * query(text_offset_type start, text_offset_type len){
    char * s = new char[len+1];
    for(text_offset_type i = 0 ; i < len; i++)
    {
      s[i] = query_alpha(start + i,0,-1);
    }
    s[len]='\0';
    return s; 
  }

};

#endif  // __COMPUTE_ST_ATT_HPP_INCLUDED