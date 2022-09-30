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
#include "rmq.hpp"
#include <cstring>
#include <map>
using namespace std;

template <
    typename char_type = std::uint8_t,
    typename text_offset_type = std::int64_t>
struct Block
{
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

template <
    typename char_type = std::uint8_t,
    typename text_offset_type = std::int64_t>
struct linked_indexes
{
  typedef std::pair<text_offset_type, text_offset_type> pair_type;
  pair_type p;

  linked_indexes(char_type *text, text_offset_type start,
                 text_offset_type end, std::vector<text_offset_type> &att_pos,
                 text_offset_type len, rmq<> *sa_rmq,
                 text_offset_type * const sa)
  {
    text_offset_type attractor_loc_, offset;
    find(text, start, end, att_pos, &attractor_loc_, &offset, len,sa_rmq,sa);
    p = std::make_pair<text_offset_type, text_offset_type>((long int)attractor_loc_,
                                                           (long int)offset);
  }

  long unsigned int start()
  {
    return p.first - p.second;
  }

  static void find(char_type *text, text_offset_type start, text_offset_type end,
                   std::vector<text_offset_type> &att_pos, text_offset_type *attractor_loc_,
                   text_offset_type *offset, text_offset_type len, rmq<> *sa_rmq,
                   text_offset_type * const sa)
  {
    *attractor_loc_ = -1;
    *offset = -1;
    end = min(len-1,end);
    if (end  < 0 || start >=len)
      return;
    text_offset_type low=0, hi = len-1,r1,r2,mid,temp_st;
    while(low < hi)
    {
      mid = (low+hi)/2;
      temp_st = sa[mid];
      text_offset_type j = start;
      for(; j <=end;j++)
      {
        if((temp_st+j-start) == len || text[temp_st + j - start] < text[j])
        {
          low = mid + 1;
          break;
        }
        else if(text[temp_st + j - start] > text[j]){
          hi = mid - 1;
          break;
        }
      }
      if(j==(end+1))
        hi = mid;
    }

    r1 = low;
    low = 0;
    hi = len - 1;
    while(low < hi)
    {
      mid = (low+hi+1)/2;
      temp_st = sa[mid];
      text_offset_type j = start;
      for(; j <=end;j++)
      {
        if((temp_st+j-start) == len || text[temp_st + j - start] < text[j]){
          low = mid + 1;
          break;
        }
        else if(text[temp_st + j - start] > text[j]){
          hi = mid - 1;
          break;
        }
      }
      if(j==(end+1))
        low = mid;
    }
    r2 = hi;
    text_offset_type x = sa_rmq->rmq_min(r1,r2);
    for(text_offset_type k = 0 ; k < (long int)att_pos.size(); k++)
    {
      if(att_pos[k]>=x && att_pos[k] <= (x + end - start))
      {
            *attractor_loc_ = k;
            *offset = att_pos[k] - x;
            return;
      }
    }

    //This is very slow !! need to implement RMQ with binary search later on toop of it
  }
};

template <
    typename char_type = std::uint8_t,
    typename text_offset_type = std::int64_t>
std::vector<linked_indexes<> *> make_linked_indexes(char_type *text, std::vector<Block<>> &b,
                                                  std::vector<text_offset_type> &att_pos, text_offset_type len,
                                                  rmq<> *sa_rmq,
                                                  text_offset_type * const sa)
{
  std::vector<linked_indexes<>*> v;
  for (text_offset_type i = 0; i < (long int)b.size(); i++)
    v.push_back((new linked_indexes<>(text, max((text_offset_type)0, b[i].start), b[i].end, att_pos, len,sa_rmq,sa)));
  return v;
}

template <
    typename char_type = std::uint8_t,
    typename text_offset_type = std::int64_t>
class st_att
{
private:
  text_offset_type tau;
  text_offset_type gamma;
  text_offset_type alpha;
  text_offset_type n;
  std::vector<text_offset_type> att_pos;
  std::vector<text_offset_type> b_si;
  std::vector<std::vector<linked_indexes<> *>> indexes;
  std::vector<char *> v_s;
  char_type *t;
  rmq<> *sa_rmq;

public:
  st_att(text_offset_type m_tau, char_type *text, text_offset_type text_length)
  {
    typedef std::pair<text_offset_type, text_offset_type> pair_type;
    n = text_length;
    text_offset_type block_len;
    tau = m_tau;
    std::vector<Block<>> v;
    t = text;
    // Allocate SA.
    text_offset_type *const sa = new text_offset_type[n];

    // Compute SA.
    {
      compute_sa(text, n, sa);
      sa_rmq = new rmq<>(sa,n);
    }

    // Compute parsing.
    std::vector<pair_type> parsing;
      compute_lz77::kkp2n(text, text_length, sa, parsing);

    //TODO ::Discuss with prof or think about it how to set alpha

    //Prepare att_pos
    std::uint64_t ind = -1;
    for (uint32_t i = 0; i < parsing.size(); i++)
    {
      ind += parsing[i].second ? parsing[i].second : 1;
      att_pos.push_back(ind);
    }
    gamma = att_pos.size();
    //Make level 0 and assign alpha
    block_len = n / gamma + (n % gamma != 0);
    for (text_offset_type i = 0; i < n; i += block_len)
      v.push_back(Block<>(i, block_len));
    b_si.push_back(block_len);
    indexes.push_back(make_linked_indexes<>(text, v, att_pos, n,sa_rmq,sa));
    alpha = max((int)ceil(log(block_len) / log(tau)), 1);
    //Now make all other levels
    while (block_len >= 2 * alpha)
    {
      block_len = block_len / tau + (block_len % tau != 0);
      b_si.push_back(block_len);
      v.clear();
      for (text_offset_type i = 0; i < (long int)att_pos.size(); i++)
      {
        text_offset_type begin = att_pos[i] - tau * block_len;
        text_offset_type end = att_pos[i] + tau * block_len;
        for (text_offset_type j = begin; j < end; j += block_len)
          v.push_back(Block<>(j, block_len));
      }
      if (block_len >= 2 * alpha)
        indexes.push_back(make_linked_indexes<>(text, v, att_pos, n,sa_rmq,sa));
      else
        break;
    }
    
    for (text_offset_type i = 0; i < (long int)v.size(); i++)
    {
      Block<> bl = v[i];
      char *s = new char[bl.len];
      for (text_offset_type j = bl.start; j < bl.end; j++)
      {
        if (j >= 0 && j < n)
          s[j - bl.start] = text[j];
        else
          s[j - bl.start] = '$';
      }
      v_s.push_back(s);
    }
    v.clear();
    att_pos.clear();
    delete[] sa;
    delete sa_rmq;
  }

  char query(text_offset_type off, uint32_t level, text_offset_type attractor)
  {
    text_offset_type block_position, offset, block_len = b_si[level];
    if (level == 0)
    {
      block_position = off / block_len;
      offset = off % block_len;
    }
    else
    {
      block_position = attractor * tau * 2 + tau + (off - attractor) / block_len;
      offset = (off - attractor) % block_len;
      if (offset < 0)
      {
        block_position--;
        offset += block_len;
      }
    }
    if (level == b_si.size() - 1)
     return  v_s[block_position][offset];

    linked_indexes<> *l = indexes[level][block_position];

    return query(
        l->start() + offset,
        level + 1,
        l->p.first);
  }
  //Query alphabet at anindex
  char query(text_offset_type index){
    return query(index,0,-1);
  }
  ~st_att(){
    for(unsigned long i=0;i<indexes.size();i++){
      for(unsigned long j=0;j<indexes[i].size();j++)
        delete(indexes[i][j]);
    }
    indexes.clear();
    for(unsigned long i=0;i<v_s.size();i++)
      delete []v_s[i];
    v_s.clear();
  }

  
};

#endif // __COMPUTE_ST_ATT_HPP_INCLUDED