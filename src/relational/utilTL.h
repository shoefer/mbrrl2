/*  
    Copyright 2008-2012   Tobias Lang
    
    E-mail:    tobias.lang@fu-berlin.de
    
    This file is part of libPRADA.

    libPRADA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libPRADA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libPRADA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef UTIL_TL_h
#define UTIL_TL_h

#include <limits>
#include <float.h>
#ifndef MT_MSVC
	#include <sys/time.h>
	#include <sys/resource.h>
#endif
#ifndef MT_IMPLEMENT_TEMPLATES
#define MT_IMPLEMENT_TEMPLATES
#endif
#include <MT/array.h>
#include <MT/util.h>
#include <cmath>

#include <vector>

#define PRINT(x) cout << #x << "=" << x << endl;
#define PRINT2(x,y) y << #x << "=" << x << endl;
#define PRINT_(x) cout << #x << " " << x << endl;
#define PRINT1(x) cout << #x << "=" << endl << x << endl;
#define PRINTOUT(x, out) out << #x << "=" << x << endl;

#define TL_MAX(a,b) (a > b ? a : b)
#define TL_MIN(a,b) (a < b ? a : b)

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
// taken from Google, 19.03.09
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)


//----- check macros:
#ifndef MT_NOCHECK
#  define CHECK_(cond,code,msg) if(!(cond)) {cout<<endl<<endl<<endl<<endl;  code; HALT("CHECK failed: "<<msg);}
#else
#  define CHECK_(cond,code,msg)
#endif

#ifdef MT_MSVC
	#define pow(b, e) pow((double)b, (double)e)
#endif

namespace TL {
const uint UINT_NIL = std::numeric_limits<uint>::max();

const double TL_DOUBLE_NIL = -98765.43211234589;
const double TL_DOUBLE_MIN = -1. * std::numeric_limits<double>::max();

#ifdef MT_MSVC
    #define TL_INFINITY numeric_limits<double>::infinity();
#else
	const double TL_INFINITY = std::numeric_limits<double>::infinity();
#endif

inline int signOf(int a) { return (a == 0) ? 0 : (a<0 ? -1 : 1); }
inline double signOf(double a) { return a<0 ? -1 : 1; }
inline bool isZero(double a) {return fabs(a) < 10e-15;}
inline bool areEqual(double a, double b) {return isZero(a-b);}

// vary from left to right (left-most argument varies the fastest)
void allPermutations(MT::Array< uintA >& permutations, const uintA& arguments, uint length, bool withRepeat, bool returnEmptyIfNoneFound);
void allPermutations(MT::Array< uintA >& permutations, const MT::Array< uintA >& arguments_lists, bool returnEmptyIfNoneFound); // different arguments

void allSubsets(MT::Array< uintA >& subsets, const uintA& elements, uint length);
void allSubsets(MT::Array< uintA >& subsets, const uintA& elements, bool trueSubsets, bool withEmpty);

bool containsAllElements(const uintA& superlist, const uintA& list);

uint basic_sample(const arr& weights);

// some ALGOS

// descending
void sort_desc(uintA& values);
void sort_desc(uintA& sorted, const uintA& unsorted);
void sort_desc(uintA& sorted, uintA& sortedIndices, const uintA& unsorted);
void sort_desc(arr& sorted, uintA& sortedIndices, const arr& unsorted);
void sort_desc_keys(uintA& sortedIndices, const uintA& unsorted);
void sort_desc_keys(uintA& sortedIndices, const arr& unsorted);
// ascending
void sort_asc(uintA& values);
void sort_asc(uintA& sorted, const uintA& unsorted);
void sort_asc(uintA& sorted, uintA& sortedIndices, const uintA& unsorted);
void sort_asc(arr& sorted, uintA& sortedIndices, const arr& unsorted);
void sort_asc_keys(uintA& sortedIndices, const uintA& unsorted);
void sort_asc_keys(uintA& sortedIndices, const arr& unsorted);

// check for cycles
bool isAcyclic(boolA adjMatrix);

uint getIndex(const uintA& constants, const uintA& args);


double getcputime();


bool uint_compare(const uint& a, const uint& b);


// Reference Managing

template<class T> void del(T* p) {
	p->numRefs--;
  CHECK(p->numRefs>=0, "Negative #references");
	if (p->numRefs == 0) {
//     PRINT(p->numRefs);
		delete p;
  }
}

template<class T> T* getRef(T* p) {
	p->numRefs++;
	return p;
}



// for my ors simulator
double REPLACE_SIZE(double val);









//===========================================================================
//
// Rprop (taken from Marc Toussaint)
//

/*! Rprop, a fast gradient-based minimization */
class Rprop{
public:
  double incr;
  double decr;
  double dMax;
  double dMin;
  double rMax;
  double delta0;

  arr lastGrad; // last error gradient
  arr stepSize; // last update

  Rprop();

  void init(double _delta0);
  void step(arr& x,const arr& grad,uint *singleI=NULL);
  void step(double& x,const double& grad);
  bool done();
};

//===========================================================================
//
// More small helpers (added by S. Hoefer)
//

template <typename T>
int argmax(const std::vector<T> data) {
  if (data.empty()) return -1;
  int amx=0;
  for (int i=1; i < data.size(); i++){
    if (data[amx] < data[i]) {
      amx = i;
    }
  }
  return amx;
}

}


#endif
