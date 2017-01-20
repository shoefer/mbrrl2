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

#include "utilTL.h"

double TL::REPLACE_SIZE(double val) {
  // Blocks
  if (TL::areEqual(val,0.04)) return 1.0;
  if (TL::areEqual(val,0.06)) return 2.0;
  if (TL::areEqual(val,0.08)) return 3.0;
  
  // Balls
  if (TL::areEqual(val,0.03)) return 1.0;
  if (TL::areEqual(val,0.045)) return 2.0;
  
  // Table
  if (val > 0.9)
    return 4.0;
  
  // TODO
  return 5.0;
  if (TL::areEqual(val,0.1)) return 5.0; 
  
  HALT("unknown size "<<val);
  return -100000.;
}


typedef MT::Array< uintA > PermutationsList;



struct PermutationsMemory { 
  uint length;
  bool returnEmptyIfNoneFound;
  bool withRepeat;
  
  // only for array_t.cxx
  PermutationsMemory() {}
  
  void create_new_list(PermutationsList& permutationsList, const uintA& args) {
    uint DEBUG = 0;
    if (DEBUG>0) {cout<<"create_new_list [START]"<<endl;}
    uint id=0;
    uint a, i, _pow;
    // beim change ganz nach hinten rutschen und wieder alle zuruecksetzen
    while (id < pow(args.N, length)) {
      if (DEBUG>2) {PRINT(length);  PRINT(args.N);  PRINT(pow(args.N, length));}
      uintA nextPermutation(length);
      i = id;
      for (a=length; a>0; a--) {
        _pow = (uint) pow(args.N, a-1);
        nextPermutation(a-1) = args(i / _pow);
        // Auf gar keinen Fall hier die Reihenfolge aendern wollen, in der Liste bestueckt wird.
        // Andere Methoden nutzen explizit die Reihenfolge aus, in der's von links nach rechts variiert.
        i = i % _pow;
      }
      id++;
      if (withRepeat)
        permutationsList.append(nextPermutation);
      else {
        if (!nextPermutation.containsDoubles())  // teuer!
          permutationsList.append(nextPermutation);
      }
    }
    if (DEBUG>0) {
      FOR1D(permutationsList, i) {PRINT(permutationsList(i))}
    }
    if (DEBUG>0) {cout<<"create_new_list [END]"<<endl;}
  }
  
  
#if 1
// Version for arguments < 100
  MT::Array< short > indices;
  MT::Array<              // arg1-arg2
    MT::Array<            // various arguments
      PermutationsList
    > 
  > mem;
  MT::Array<              // arg1-arg2
    MT::Array<            // various arguments
      uintA
    >
  > mem_args;
  
  PermutationsMemory(uint _length, bool _withRepeat, bool _returnEmpty) : length(_length), returnEmptyIfNoneFound(_returnEmpty), withRepeat(_withRepeat) {
    indices.resize(100, 100);
    indices.setUni(100);
  }
  
  ~PermutationsMemory() {
  }
  
  void get(PermutationsList& permutationsList, const uintA& args) {
    uint DEBUG = 0;
//     if (length>=2) DEBUG = 2;
    if (DEBUG>0) {cout<<"PermutationsMemory::get [START]"<<endl;}
    if (DEBUG>0) {PRINT(length);  PRINT(returnEmptyIfNoneFound);  PRINT(withRepeat);  PRINT(args);}
    permutationsList.clear();
    if (length == 1) {
      uintA one(1);
      permutationsList.resize(args.N);
      permutationsList.setUni(one);
      uint i;
      FOR1D(args, i) {
        permutationsList(i)(0) = args(i);
      }
    }
    else if (length > 1  &&  args.N>0) {
      CHECK(withRepeat || args.N >= length, "args="<<args<<" length="<<length);
      short idx;
      if (args.N > 1)
        idx = indices(args(0), args(1));
      else
        idx = indices(args(0), args(0));
      if (DEBUG>0) {PRINT(idx);}
      if (idx == 100) {
        if (DEBUG>0) {cout<<"creating memory list for args(0)-args(1)"<<endl;}
        // extend memory
        MT::Array< PermutationsList > new_mem;
        mem.append(new_mem);
        MT::Array< uintA > new_mem_args;
        mem_args.append(new_mem_args);
        // update indices
        if (args.N > 1)
          indices(args(0), args(1)) = mem.N-1;
        else
          indices(args(0), args(0)) = mem.N-1;
        idx = mem.N-1;
      }
      if (DEBUG>1) {PRINT(mem(idx));  PRINT(mem_args(idx));}
      uint i;
      FOR1D(mem_args(idx), i) {
        if (mem_args(idx)(i) == args) {
          if (DEBUG>1) {cout<<"found existing entry"<<endl;}
          permutationsList = mem(idx)(i);
          break;
        }
      }
      if (i==mem_args(idx).N) {
        if (DEBUG>1) {cout<<"creating new entry"<<endl;}
        create_new_list(permutationsList, args);
        mem(idx).append(permutationsList);
        mem_args(idx).append(args);
      }
    }
    if (permutationsList.N == 0  &&  returnEmptyIfNoneFound) {
      uintA empty(0);
      permutationsList.append(empty);
    }
    if (DEBUG>0) {PRINT(permutationsList);}
    if (DEBUG>0) {cout<<"PermutationsMemory::get [END]"<<endl;}
  }
  
  
  
  std::map< uint, uint > map_1;
  MT::Array< std::map< uint, uint > > map_2;
  
  // 0 - argument1,  1 - argument2,  2 - lists
  MT::Array<  // arg 1
    MT::Array< // arg 2
      MT::Array< // various arguments (starting with arg1 and arg2)
        PermutationsList // lists
      >
    > 
  > mem__2plus;
  
  MT::Array<  // arg 1
    MT::Array< // arg 2
      MT::Array< uintA > // arguments lists
    > 
  > mem_args__2plus;
#else
  PermutationsMemory(uint _length, bool _withRepeat, bool _returnEmpty) : length(_length), returnEmptyIfNoneFound(_returnEmpty), withRepeat(_withRepeat) {
  }
  
  ~PermutationsMemory() {
  }
  
  void get(PermutationsList& permutationsList, const uintA& args) {
    uint DEBUG = 0;
    if (DEBUG>0) {PRINT(args);}
    permutationsList.clear();
    if (length == 1) {
      permutationsList.resize(args.N);
      uint i;
      FOR1D(args, i) {
        permutationsList.elem(i).resize(1);
        permutationsList.elem(i)(0) = args(i);
      }
    }
    else {
      if (map_1.count(args(0)) == 0) {
        // key
        std::map< uint, uint> new_map2;
        map_2.append(new_map2);
        map_1[args(0)] = map_2.N-1;
        // memory - permutations
        MT::Array< MT::Array< PermutationsList > > mem_list;
        mem__2plus.append(mem_list);
        CHECK(mem__2plus.N-1 == map_1[args(0)], "");
        // memory - arguments
        MT::Array< MT::Array< uintA > > mem_args_list;
        mem_args__2plus.append(mem_args_list);
      }
      int idx1 = map_1[args(0)];
      if (DEBUG>0) {PRINT(idx1);}
      if (map_2(idx1).count(args(1)) == 0) {
        // memory - permutations
        MT::Array< PermutationsList > mem_list;
        if (DEBUG>0) {PRINT(mem_list);}
        mem__2plus(idx1).append(mem_list);
        if (DEBUG>0) {PRINT(mem__2plus(idx1));}
        map_2(idx1)[args(1)] = mem__2plus(idx1).N-1;
        // memory - arguments
        MT::Array< uintA > mem_args_list;
        mem_args__2plus(idx1).append(mem_args_list);
      }
      int idx2 = map_2(idx1)[args(1)];
      if (DEBUG>0) {PRINT(idx2);}
      if (DEBUG>0) {PRINT(mem__2plus(idx1));}
      MT::Array< PermutationsList> & mem_prefix = mem__2plus(idx1)(idx2);
      if (DEBUG>0) {PRINT(mem_prefix);}
      MT::Array< uintA > & mem_args_prefix = mem_args__2plus(idx1)(idx2);
      if (DEBUG>0) {PRINT(mem_args_prefix);}
      uint i;
      FOR1D(mem_args_prefix, i) {
        if (mem_args_prefix(i) == args) {
          permutationsList = mem_prefix(i);
          break;
        }
      }
      // not found --> fill in
      if (mem_args_prefix.N == i) {
        mem_args_prefix.append(args);
        create_new_list(permutationsList, args);
        mem_prefix.append(permutationsList);
      }
    }
    if (permutationsList.N == 0  &&  returnEmptyIfNoneFound) {
      uintA empty(0);
      permutationsList.append(empty);
    }
    if (DEBUG>0) {PRINT(permutationsList);}
  }
#endif
};


struct PermutationsMemoryWrapper {
  MT::Array< PermutationsMemory* > all_mem;
  bool withRepeat;
  bool returnEmptyIfNoneFound;
  
  PermutationsMemoryWrapper(bool _withRepeat, bool _returnEmpty) : withRepeat(_withRepeat), returnEmptyIfNoneFound(_returnEmpty) {
  }
  
  ~PermutationsMemoryWrapper() {listDelete(all_mem);}
  
  void get(MT::Array< uintA >& output_lists, const uintA& args, uint length) {
    while (length >= all_mem.N) {
      all_mem.append(new PermutationsMemory(all_mem.N, withRepeat, returnEmptyIfNoneFound));
    } 
    all_mem(length)->get(output_lists, args);
  }
};


PermutationsMemoryWrapper mem__withRepeat_returnEmpty(true, true);
PermutationsMemoryWrapper mem__withRepeat_dontReturnEmpty(true, false);
PermutationsMemoryWrapper mem__withoutRepeat_returnEmpty(false, true);
PermutationsMemoryWrapper mem__withoutRepeat_dontReturnEmpty(false, false);


void TL::allPermutations(PermutationsList& permutationsList, const uintA& arguments, uint length, bool withRepeat, bool returnEmptyIfNoneFound) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"allPermutations [START]"<<endl;
  if (DEBUG>0) {PRINT(arguments); PRINT(length); PRINT(withRepeat);  PRINT(returnEmptyIfNoneFound);}
  
  permutationsList.clear();
  // Try to get from memory
  if (withRepeat && returnEmptyIfNoneFound) {
    mem__withRepeat_returnEmpty.get(permutationsList, arguments, length);
  }
  else if (withRepeat && !returnEmptyIfNoneFound) {
    mem__withRepeat_dontReturnEmpty.get(permutationsList, arguments, length);
  }
  else if (!withRepeat && returnEmptyIfNoneFound) {
    mem__withoutRepeat_returnEmpty.get(permutationsList, arguments, length);
  }
  else {
    mem__withoutRepeat_dontReturnEmpty.get(permutationsList, arguments, length);
  }
  if (DEBUG>0) {
    uint i;
    FOR1D(permutationsList, i) {PRINT(permutationsList(i));}
  }
  if (DEBUG>0) cout<<"allPermutations [END]"<<endl;
}


// different arguments
void TL::allPermutations(MT::Array< uintA >& lists, const MT::Array< uintA >& arguments_lists, bool returnEmptyIfNoneFound)  {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"allPermutations [START]"<<endl;}
  if (DEBUG>0) {PRINT(arguments_lists);  PRINT(returnEmptyIfNoneFound);}
  lists.clear();
  if (arguments_lists.N == 0) {
    if (returnEmptyIfNoneFound) {
      uintA empty(0);
      lists.append(empty);
//       if (DEBUG>0) cout<<"allPermutations [END]"<<endl;
      return;
    }
    else {
      MT_MSG("No lists returned!");
      return;
    }
  }
  uint id=0;
  uint a, i, k, _num_configs;
  // beim change ganz nach hinten rutschen und wieder alle zuruecksetzen
  uint num_configs = 1;
  FOR1D(arguments_lists, i) {num_configs *= arguments_lists(i).N;}
  if (DEBUG>0) {PRINT(arguments_lists.N);  PRINT(num_configs);}
  while (id < num_configs) {
    uintA nextList(arguments_lists.N);
    i = id;
    for (a=arguments_lists.N; a>0; a--) {
      _num_configs = 1;
      for (k=a; k>1; k--) {_num_configs *= arguments_lists(k-2).N;}
      nextList(a-1) = arguments_lists(a-1)(i / _num_configs);
	// Auf gar keinen Fall hier die Reihenfolge aendern wollen, in der Liste bestueckt wird.
	// Andere Methoden nutzen explizit die Reihenfolge aus, in der's von links nach rechts variiert.
      i = i % _num_configs;
    }
    lists.append(nextList);
    id++;
  }
  
  if (DEBUG>0) {PRINT(lists);}
  if (DEBUG>0) {cout<<"allPermutations [END]"<<endl;}
}


void TL::allSubsets(MT::Array< uintA >& lists, const uintA& elements, bool trueSubsets, bool withEmpty) {
  uint DEBUG = 0;
  uint length;
  uint max_length = elements.N;
  uint i, k;
  if (trueSubsets)
    max_length -= 1;
  // others
  for (length=1; length<=max_length; length++) {
    MT::Array< uintA > local_lists;
    allPermutations(local_lists, elements, length, false, false);
    FOR1D(local_lists, i) {
      for (k=0; k<local_lists(i).N-1; k++) {
        if (local_lists(i)(k) >= local_lists(i)(k+1))
          break;
      }
      if (k == local_lists(i).N-1)
        lists.append(local_lists(i));
    }
  }
  // empty list
  if (withEmpty) {
    uintA empty;
    lists.append(empty);
  }
  
  if (DEBUG > 0) {
    cout<<"allSubsets:"<<endl;
    PRINT(elements);
    PRINT(trueSubsets);
    PRINT(withEmpty);
    cout<<"Result:"<<endl;
    PRINT(max_length);
    PRINT(lists);
  }
}

void TL::allSubsets(MT::Array< uintA >& lists, const uintA& elements, uint length) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"allSubsets [START]"<<endl;}
  if (DEBUG>0) {
    PRINT(elements);
    PRINT(length);
  }
  CHECK(length<=elements.N, "too long subset size");
  if (length==0) {
    uintA empty;
    lists.append(empty);
    return;
  }
  uint i, k, l;
  uintA marker(length);
  FOR1D(marker, i) {marker(i)=i;}
  bool change;
  while(true) {
    uintA list;
    FOR1D(marker, i) {
      list.append(elements(marker(i)));
    }
    if (DEBUG>0) {PRINT(marker);cout<<"New list: "<<list<<endl;}
    lists.append(list);
    change = false;
    FOR1D_DOWN(marker, k) {
      if (marker(k) < elements.N-(length-k)) {
        marker(k)++;
        for(l=k+1;l<marker.N;l++)
          marker(l)=marker(k)+l-k;
        change = true;
        break;
      }
    }
    if (!change)
      break;
  }
  if (DEBUG>0) {cout<<"All subsets: "<<lists<<endl;}
  if (DEBUG>0) {cout<<"allSubsets [END]"<<endl;}
}

bool TL::containsAllElements(const uintA& superlist, const uintA& list) {
  uint i;
  FOR1D(list, i) {
    if (superlist.findValue(list(i)) < 0)
      return false;
  }
  return true;
}


uint TL::basic_sample(const arr& weights) {
    arr probs_normalized = weights;
    if (!TL::areEqual(sum(probs_normalized), 1.0)) {
        normalizeDist(probs_normalized);
    }
    double rndNum = rnd.uni();
    double massSoFar = 0.0;
    uint i;
    FOR1D(weights, i) {
      massSoFar += probs_normalized(i);
      if (massSoFar >= rndNum)
        break;
    }
    CHECK(weights.N > i, "sampling from dist failed");
    return i;
}


// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//    SORT


void TL::sort_desc(arr& sorted, uintA& sortedIndices, const arr& unsorted) {
  const double LOCAL_MIN = -333333.3;
  sortedIndices.clear();
  sorted.clear();
  arr un_copy = unsorted;
  uint i, k;
  double max;
  uint max_id;
  for (i=0; i<un_copy.N; i++) {
    max = LOCAL_MIN;
    max_id = TL::UINT_NIL;
    FOR1D(un_copy, k) {
      if (un_copy(k) > max) {
        max = un_copy(k);
        max_id = k;
      }
    }
    sorted.append(max);
    sortedIndices.append(max_id);
    un_copy(max_id) = LOCAL_MIN;
  }
  CHECK(sorted.N == unsorted.N, "");
}


void TL::sort_desc(uintA& sorted, uintA& sortedIndices, const uintA& unsorted) {
  sortedIndices.clear();
  sorted.clear();
  uint rank;
  uint i;
  uint maxID;
  boolA used(unsorted.N);
  used.setUni(false);
  for(rank=0; rank<unsorted.N; rank++) {
    maxID=0;
    while(used(maxID))  // initial set of maxID (= kleinstes verbleibendes Element)
      maxID++;
    for(i=maxID+1; i<unsorted.N; i++) {
      if (used(i))
        continue;
      if (unsorted(maxID) < unsorted(i))
        maxID = i;
    }
    sortedIndices.append(maxID);
    sorted.append(unsorted(maxID));
    used(maxID) = true;
  }
  CHECK(sorted.N == unsorted.N, "");
}

void TL::sort_desc(uintA& values) {
  uintA sorted, sortedIndices;
  sort_desc(sorted, sortedIndices, values);
  values = sorted;
}

void TL::sort_desc_keys(uintA& sortedIndices, const uintA& unsorted) {
  uintA sorted;
  sort_desc(sorted, sortedIndices, unsorted);
}

void TL::sort_desc_keys(uintA& sortedIndices, const arr& unsorted) {
  arr sorted;
  sort_desc(sorted, sortedIndices, unsorted);
}



void TL::sort_asc(uintA& values) {
  sort_desc(values);
  values.reverse();
}

void TL::sort_asc(uintA& sorted, const uintA& unsorted) {
  uintA sortedIndices;
  sort_asc(sorted, sortedIndices, unsorted);
}

void TL::sort_asc(uintA& sorted, uintA& sortedIndices, const uintA& unsorted) {
  sort_desc(sorted, sortedIndices, unsorted);
  sorted.reverse();
  sortedIndices.reverse();
}

void TL::sort_asc(arr& sorted, uintA& sortedIndices, const arr& unsorted) {
  sort_desc(sorted, sortedIndices, unsorted);
  sorted.reverse();
  sortedIndices.reverse();
}

void TL::sort_asc_keys(uintA& sortedIndices, const uintA& unsorted) {
  sort_desc_keys(sortedIndices, unsorted);
  sortedIndices.reverse();
}











// some testing wouldn't be all too bad
bool TL::isAcyclic(boolA adjMatrix) {
  // eliminate successively nodes without outgoing or without incoming 
  // edges since these cannot be part of a cycle
  CHECK(adjMatrix.nd==2, "adjMatrix must be square")
  CHECK(adjMatrix.d0==adjMatrix.d1, "adjMatrix must be square")
  boolA endangeredNodes(adjMatrix.d0);
  endangeredNodes.setUni(true);
  bool change;
  uint n_from, n_to;
  do {
    change=false;
//         PRINT(endangeredNodes);
      // without outgoing
      for(n_from=0; n_from<endangeredNodes.d0; n_from++) {
          if (!endangeredNodes(n_from)) continue;
          for (n_to=0; n_to<endangeredNodes.d0; n_to++) {
              if (!endangeredNodes(n_to)) continue;
              else if (adjMatrix(n_from, n_to)) break;
          }
          if (n_to == endangeredNodes.d0) {
              endangeredNodes(n_from) = false;
              change = true;
          }
      }
      // without incoming
      for(n_to=0; n_to<endangeredNodes.d0; n_to++) {
          if (!endangeredNodes(n_to)) continue;
          for (n_from=0; n_from<endangeredNodes.d0; n_from++) {
              if (!endangeredNodes(n_from)) continue;
              else if (adjMatrix(n_from, n_to)) break;
          }
          if (n_from == endangeredNodes.d0) {
              endangeredNodes(n_to) = false;
              change = true;
          }
      }
    } while (change);
    bool cyclic = sum(endangeredNodes);
    return !cyclic;
}



double TL::getcputime() {
  double t = 0;
#ifndef MT_MSVC
  struct timeval tim;
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  tim=ru.ru_utime;
  t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
  tim=ru.ru_stime;
  t+=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
#endif
  return t;
}

uint TL::getIndex(const uintA& constants, const uintA& args) {
  uint args_idx=0;
  uint i;
  FOR1D(args, i) {
    args_idx += ((uint) pow(constants.N, i)) * constants.findValue(args(i));
  }
//   cout<<"getIndex: constants="<<constants<<"  args="<<args<<"    args_idx="<<args_idx<<endl;
  return args_idx;
}


bool TL::uint_compare(const uint& a,const uint& b) { return a<b; }



//===========================================================================
//
// Rprop
//

int _sgn(double x){ if (x > 0) return 1; if (x < 0) return -1; return 0; }
double _mymin(double x,double y){ return x < y ? x : y; }
double _mymax(double x,double y){ return x > y ? x : y; }


TL::Rprop::Rprop(){
  incr   = 1.2;
  decr   = .33;
  dMax = 50;
  dMin = 1e-6;
  rMax = 0;
  delta0 = 1.;
}

void TL::Rprop::init(double _delta0){
  stepSize.resize(0);
  lastGrad.resize(0);
  delta0 = _delta0;
}

bool TL::Rprop::done(){
  double maxStep = stepSize(stepSize.maxIndex());
  return maxStep < incr*dMin;
}

void TL::Rprop::step(double& w,const double& grad){
  static arr W,GRAD;
  W.referTo(&w,1); GRAD.referTo(&grad,1);
  step(W,GRAD);
}

void TL::Rprop::step(arr& w,const arr& grad,uint *singleI){
  if(!stepSize.N){ //initialize
    stepSize.resize(w.N);
    lastGrad.resize(w.N);
    lastGrad.setZero();
    stepSize = delta0;
  }
  CHECK(grad.N==stepSize.N,"Rprop: gradient dimensionality changed!");
  CHECK(w.N==stepSize.N   ,"Rprop: parameter dimensionality changed!");

  uint i=0,I=w.N;
  if(singleI){ i=*(singleI); I=i+1; }
  for(;i<I;i++){
    if(grad.elem(i) * lastGrad(i) > 0){        //same direction as last time
      if(rMax) dMax=fabs(rMax*w.elem(i));
      stepSize(i) = _mymin(dMax, incr * stepSize(i)); //increase step size
      w.elem(i) += stepSize(i) * -_sgn(grad.elem(i)); //step in right direction
      lastGrad(i) = grad.elem(i);                    //memorize gradient
    }else if(grad.elem(i) * lastGrad(i) < 0){  //change of direction
      stepSize(i) = _mymax(dMin, decr * stepSize(i)); //decrease step size
      w.elem(i) += stepSize(i) * -_sgn(grad.elem(i)); //step in right direction
      lastGrad(i) = 0;                               //memorize to continue below next time
    }else{                                     //after change of direcion
      w.elem(i) += stepSize(i) * -_sgn(grad.elem(i)); //step in right direction 
      lastGrad(i) = grad.elem(i);                    //memorize gradient
    }
  }
}


