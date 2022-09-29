// C++ program for range minimum
// query using segment tree
#include <bits/stdc++.h>
using namespace std;

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::int64_t>
class rmq
{
	private:
		vector<text_offset_type> v;
		text_offset_type n;
		text_offset_type minVal(text_offset_type x, text_offset_type y) { return (x < y)? x: y; }
		text_offset_type getMid(text_offset_type s, text_offset_type e) { return s + (e -s)/2; }
		text_offset_type rmq_construct(text_offset_type ss, text_offset_type se, text_offset_type * const arr,text_offset_type si)
		{
			if (ss == se)
			{
				v[si] = arr[ss];
				return arr[ss];
			}

			text_offset_type mid = getMid(ss, se);
			v[si] = minVal(rmq_construct(ss, mid, arr, si*2+1),
					rmq_construct(mid+1, se, arr, si*2+2));
			return v[si];
		}

		text_offset_type RMQUtil(int ss, int se, int qs, int qe, int index)
		{
			if (qs <= ss && qe >= se)
				return v[index];

			if (se < qs || ss > qe)
				return INT_MAX;

			
			text_offset_type mid = getMid(ss, se);
			return minVal(RMQUtil(ss, mid, qs, qe, 2*index+1),
				RMQUtil(mid+1, se, qs, qe, 2*index+2));
		}
	public:
		rmq(text_offset_type * const arr,text_offset_type m_n)
		{
			n = m_n;
			text_offset_type x = (text_offset_type)(ceil(log2(n)));
			text_offset_type max_size = 2*(text_offset_type)pow(2, x) - 1;
			v.resize(max_size);
			rmq_construct(0, n-1, arr,0);			
		}
		
		text_offset_type rmq_min(text_offset_type qs, text_offset_type qe)
		{
			if (qs < 0 || qe > n-1 || qs > qe)
			{
				cout<<"Invalid Input ";
				return -1;
			}
			return RMQUtil(0, n-1, qs, qe, 0);
		}

};

