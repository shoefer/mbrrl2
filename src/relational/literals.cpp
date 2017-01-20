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
    GNU General Public License for more details->

    You should have received a copy of the GNU General Public License
    along with libPRADA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "symbols.h"
#include "reason.h" // for derive in state

#define NOISE_MARKER 'N'

namespace relational {

  
/************************************************
 * 
 *     LiteralStorage
 * 
 ************************************************/

//Maps non-binary symbols f(X)=[1..n] to their abstract counterpart f(X)=C
std::map<Symbol*, Symbol> varEqualitySymbols;

struct LiteralStorage {
  Symbol* s;
#if 1
  LitL mem_arity0;
  MT::Array< short > indices;
  MT::Array< LitL* > mem;
    
  LiteralStorage(Symbol* _s) : s(_s) {
    if (s->arity == 1) {
      indices.resize(100);
      indices.setUni(100);
    }
    else if (s->arity >= 2) {
      indices.resize(100, 100);
      indices.setUni(100);
    }
  }
  LiteralStorage() {} // only for array_t.cpp
  ~LiteralStorage() {
    uint i;
    if (s->arity == 0) {listDelete(mem_arity0);}
    else {
      FOR_ALL(mem, i) {
        listDelete((*mem(i)));
        delete mem(i);
      }
    }
  }
#else
  std::map<uint, uint> map_indices_arg1;
  MT::Array< std::map<uint, uint> > maps_indices_arg2;
  
  LitL mem_arity0;
  MT::Array< LitL >  mem_arity1;
  MT::Array< MT::Array< LitL> >  mem_arity2plus;
  
  
  LiteralStorage(Symbol* _s) : s(_s) {
  }
  LiteralStorage() {} // only for array_t.cpp
  ~LiteralStorage() {
    uint i, k;
    if (s->arity == 0) {listDelete(mem_arity0);}
    else if (s->arity == 1) {
      FOR_ALL(mem_arity1, i) {
        listDelete(mem_arity1(i));
      }
    }
    else {
      FOR_ALL(mem_arity2plus, i) {
        FOR_ALL(mem_arity2plus(i), k) {
          listDelete(mem_arity2plus(i)(k));
        }
      }
    }
  }
#endif
};


struct LiteralStorage_Container {
  MT::Array< LiteralStorage* > literal_storages;
  LiteralStorage_Container() {}
  ~LiteralStorage_Container() {listDelete(literal_storages);}
};

LiteralStorage_Container lsc;


Literal* __getLiteral(LitL& mem, const uintA& args, double value, Literal::ComparisonType comparison_type) {
  uint i;
  if (args.N <= 2) {
    FOR1D(mem, i) {
      if (TL::areEqual(mem(i)->value, value)  &&  mem(i)->comparison_type == comparison_type) {
        return mem(i);
      }
    }
  }
  else {
    FOR1D(mem, i) {
      if (mem(i)->args == args  &&  TL::areEqual(mem(i)->value, value)  &&  mem(i)->comparison_type == comparison_type) {
        return mem(i);
      }
    }
  }
  return NULL;
}


/************************************************
 * 
 *     Literal
 * 
 ************************************************/

Literal::Literal() {}

#if 1
// Version for arguments < 100
Literal* Literal::get(Symbol* s, const uintA& args, double value, ComparisonType comparison_type) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"Literal::get[START]"<<endl;}
  if (DEBUG>0) {cout<<"asking:  "<<s->name<<" args="<<args<<" value="<<value<<" comparison_type="<<comparison_type<<endl;}
  uint i;
  FOR1D(args, i) {if (args(i) > 99) HALT("arguments must be <=99");}
  if (s->arity != args.N) {HALT("s->arity!=args.N  s="<<*s<<"    args="<<args);}
  uint p_id = 0;
  FOR1D(lsc.literal_storages, i) {
    if (lsc.literal_storages.p[p_id++]->s == s)
      break;
  }
  if (i==lsc.literal_storages.N) {
    lsc.literal_storages.append(new LiteralStorage(s));
  }
  LiteralStorage* ls = lsc.literal_storages(i);
  Literal* l = NULL;
  if (s->arity == 0) {
    l = __getLiteral(ls->mem_arity0, args, value, comparison_type);
    if (l==NULL) {
      l = new Literal;
      l->s = s;
      l->args = args;
      l->value = value;
      l->comparison_type = comparison_type;
      ls->mem_arity0.append(l);
    }
  }
  else {
    uint idx1 = 10000;
    if (s->arity == 1) {
      idx1 = ls->indices(args(0));
      if (DEBUG>0) {PRINT(idx1);}
      if (idx1 == 100) {
        LitL* new_lits_list = new LitL;
        ls->mem.append(new_lits_list);
        ls->indices(args(0)) = ls->mem.N-1;
        idx1 = ls->mem.N-1;
      }
    }
    else {
      idx1 = ls->indices(args(0), args(1));
      if (DEBUG>0) {PRINT(idx1);}
      if (idx1 == 100) {
        LitL* new_lits_list = new LitL;
        ls->mem.append(new_lits_list);
        ls->indices(args(0), args(1)) = ls->mem.N-1;
        idx1 = ls->mem.N-1;
      }
    }
    l = __getLiteral(*ls->mem(idx1), args, value, comparison_type);
    if (DEBUG>0) {cout<<(l!=NULL?"found":"not found")<<endl;}
    if (l==NULL) {
      l = new Literal;
      l->s = s;
      l->args = args;
      l->value = value;
      l->comparison_type = comparison_type;
      ls->mem(idx1)->append(l);
    }
    if (DEBUG>0) {PRINT(ls->mem(idx1));}
  }
  CHECK(l!=NULL, "literal could not be made");
  if (DEBUG>0) {cout<<"getting:  "<<*l<<" "<<l<<endl;}
  if (l->s->range_type == Symbol::binary) {
    if (!TL::areEqual(l->value, 0)  &&  !TL::areEqual(l->value, 1))
      HALT("Literal "<<*l<<" is defined for a binary symbol "<<*l->s<<", but has value "<<l->value);
  }
  if (DEBUG>0) {CHECK(l->args.N==0 || l->args == args, "");}
  if (DEBUG>0) {cout<<"Literal::get[END]"<<endl;}
  return l;
}
#else
Literal* Literal::get(Symbol* s, const uintA& args, double value, ComparisonType comparison_type) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"Literal::get[START]"<<endl;}
  if (DEBUG>0) {cout<<"asking:  "<<s->name<<" args="<<args<<" value="<<value<<" comparison_type="<<comparison_type<<endl;}
  if (s->arity != args.N) {HALT("s->arity!=args.N  s="<<*s<<"    args="<<args);}
  uint i;
  uint p_id = 0;
  FOR1D(lsc.literal_storages, i) {
    if (lsc.literal_storages.p[p_id++]->s == s)
      break;
  }
  if (i==lsc.literal_storages.N) {
    lsc.literal_storages.append(new LiteralStorage(s));
  }
  LiteralStorage* ls = lsc.literal_storages(i);
  Literal* l = NULL;
  if (s->arity == 0) {
    l = __getLiteral(ls->mem_arity0, args, value, comparison_type);
    if (l==NULL) {
      l = new Literal;
      l->s = s;
      l->args = args;
      l->value = value;
      l->comparison_type = comparison_type;
      ls->mem_arity0.append(l);
    }
  }
  else if (s->arity == 1) {
    std::map<uint, uint>::iterator iter(ls->map_indices_arg1.lower_bound(args(0)));
    uint idx_1;
    if (iter == ls->map_indices_arg1.end()  ||  args(0) < iter->first) {
      LitL new_lits_list;
      ls->mem_arity1.append(new_lits_list);
      std::pair<uint,uint> pair = std::make_pair(args(0), ls->mem_arity1.N-1);
      ls->map_indices_arg1.insert(iter, pair);
      idx_1 = pair.second;
    }
    else
      idx_1 = iter->second;
    l = __getLiteral(ls->mem_arity1(idx_1), args, value, comparison_type);
    if (l==NULL) {
      l = new Literal;
      l->s = s;
      l->args = args;
      l->value = value;
      l->comparison_type = comparison_type;
      ls->mem_arity1(idx_1).append(l);
    }
  }
  // arity >= 2
  else {
    std::map<uint, uint>::iterator iter(ls->map_indices_arg1.lower_bound(args(0)));
    uint idx_1;
    if (iter == ls->map_indices_arg1.end()  ||  args(0) < iter->first) {
      // add memory
      MT::Array< LitL> new_outer_list;
      ls->mem_arity2plus.append(new_outer_list);
      // add to args_pos2
      std::map< uint, uint> new_arg2_map;
      ls->maps_indices_arg2.append(new_arg2_map);
      // insert in args_pos1
      std::pair<uint,uint> pair = std::make_pair(args(0), ls->mem_arity2plus.N-1);
      ls->map_indices_arg1.insert(iter, pair);
      idx_1 = pair.second;
    }
    else
      idx_1 = iter->second;
    std::map<uint, uint>::iterator iter2(ls->maps_indices_arg2(idx_1).lower_bound(args(1)));
    uint idx_2;
    if (iter2 == ls->maps_indices_arg2(idx_1).end()  ||  args(1) < iter2->first) {
      LitL new_lits_list;
      ls->mem_arity2plus(idx_1).append(new_lits_list);
      std::pair<uint,uint> pair = std::make_pair(args(1), ls->mem_arity2plus(idx_1).N-1);
      ls->maps_indices_arg2(idx_1).insert(iter2, pair);
      idx_2 = pair.second;
    }
    else
      idx_2 = iter2->second;
    if (DEBUG>1) {PRINT(idx_2);}
    l = __getLiteral(ls->mem_arity2plus(idx_1)(idx_2), args, value, comparison_type);
    if (l==NULL) {
      l = new Literal;
      l->s = s;
      l->args = args;
      l->value = value;
      l->comparison_type = comparison_type;
      ls->mem_arity2plus(idx_1)(idx_2).append(l);
    }
  }
  CHECK(l!=NULL, "literal could not be made");
  if (DEBUG>0) {cout<<"getting:  "<<*l<<endl;}
  if (DEBUG>0) {CHECK(l->args.N==0 || l->args == args, "");}
  if (DEBUG>0) {cout<<"Literal::get[END]"<<endl;}
  return l;
}
#endif

Literal* Literal::getVarComparison(Symbol* s, const uintA& args) {
  Symbol *varSymbol = Symbol::get(s->name, s->arity, s->symbol_type, s->range_type);
  return get(varSymbol, args, -1, Literal::comparison_variable);   //dummy value
}


uint variable__char2digit(const char c) {
  if (c == 'X') return 0;
  else if (c == 'Y') return 1;
  else if (c == 'Z') return 2;
  else if (c == 'V') return 3;
  else if (c == 'W') return 4;
  else if (c == 'U') return 5;
  else HALT("unknown variable "<<c);
}


char variable__digit2char(uint u) {
  CHECK(u<10, "only for digits;  u="<<u);
  if (u == 0) return 'X';
  else if (u == 1) return 'Y';
  else if (u == 2) return 'Z';
  else if (u == 3) return 'V';
  else if (u == 4) return 'W';
  else if (u == 5) return 'U';
  else return u + '0';
}


Literal* Literal::get(const char* text) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"Literal::get [START]"<<endl;}
  MT::String mt_string(text);
  if (DEBUG>0) {PRINT(mt_string);}
  
  double value = 1.;
  
  // Special negation with "-" for binary predicates
  if (MT::peerNextChar(mt_string) == '-') {
    value = 0.;
    skipOne(mt_string);
  }
  
  // Symbol
  MT::String name;
  name.read(mt_string, NULL, "(");
  if (DEBUG>0) PRINT(name);
  Symbol* s = Symbol::get(name);
  if (s == NULL) HALT("Symbol with name \""<<name<<"\" is unknown when reading literal \""<<text<<"\". Have you provided and read the correct symbols file?");
  if (DEBUG>0) {PRINT(*s);}
  
  // Arguments
  if (DEBUG>0) {cout<<"Reading arguments"<<endl;}
  uintA args(s->arity);
  uint i;
  FOR1D(args, i) {
    MT::String arg;
    arg.read(mt_string, NULL, "/, )");
    if (DEBUG>0) {PRINT(arg);}
    if (isdigit(MT::peerNextChar(arg))) {
      arg >> args(i);
    }
    else {  // for string variables
      if (arg.N != 1) HALT("non-digit argument in bad format: "<<arg);
      args(i) = variable__char2digit(arg(0));
    }
  }
  if (s->arity == 0) {
    MT::String arg;
    arg.read(mt_string, NULL, ")");
  }

  Literal* l = NULL;
  
  // Value
  Literal::ComparisonType comp_type;
  if (MT::peerNextChar(mt_string) == -1  &&  s->range_type == Symbol::binary) {
    comp_type = Literal::comparison_equal;
  }
  else {
    char op;
    mt_string >> op;
    if (DEBUG>0) {PRINT(op);}
    if (op == '>') {
      if (MT::peerNextChar(mt_string) == '=') {
        mt_string >> op;
        comp_type = Literal::comparison_greaterEqual;
      }
      else
        comp_type = Literal::comparison_greater;
    }
    else if (op == '<') {
      if (MT::peerNextChar(mt_string) == '=') {
        mt_string >> op;
        comp_type = Literal::comparison_lessEqual;
      }
      else
        comp_type = Literal::comparison_less;
    }
    else if (op == '+') {
      if (MT::peerNextChar(mt_string) == '=') {
        mt_string >> op;
        comp_type = Literal::comparison_offset;
      }
    }
    else
      comp_type = Literal::comparison_equal;

    if (MT::peerNextChar(mt_string) == '?')
      l = Literal::getVarComparison(s, args);
    else
      mt_string >> value;
    if (DEBUG>0) {PRINT(value);}
  }

  if (!l)  l = Literal::get(s, args, value, comp_type);
  if (DEBUG>0) {cout<<"==> "<<*l<<endl;}
  if (DEBUG>0) {cout<<"Literal::get [END]"<<endl;}
  return l;
}


void Literal::get(LitL& lits, const char* text) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"Literal::get(LitL) [START]"<<endl;}
  if (DEBUG>0) {PRINT(text);}
  lits.clear();
  MT::String mt_string(text);
  MT::String name;
  while (MT::skip(mt_string) != -1) {
    MT::String lit_text;
    lit_text.read(mt_string, NULL, ")", 1);
    if (MT::peerNextChar(mt_string) == ',') MT::skipOne(mt_string);
//     MT::skip(std::istream& is,const char *skipchars=" \n\r\t",bool skipCommentLines=true);
    char c = MT::peerNextChar(mt_string);
    if (c == '='  ||  c == '>'  ||  c == '<' || c == '+') {
      MT::String comparison;
      comparison.read(mt_string, NULL, " \n", 1);
      lit_text << ")" << comparison;
    }
    if (DEBUG>0) {PRINT(lit_text);}
    lits.append(Literal::get(lit_text));
  }
  if (DEBUG>0) {relational::write(lits); cout<<endl;}
  if (DEBUG>0) {cout<<"Literal::get(LitL) [END]"<<endl;}
}


bool Literal::operator==(Literal& l) const {
  if (s != l.s)
    return false;
  if (comparison_type != l.comparison_type)
    return false;
  if (value != l.value)
    return false;
  if (this->args.N==0 && l.args.N==0) // special case if no slot assignments
    return true;
  if (args == l.args)
    return true;
  else
    return false;
}


bool Literal::operator!=(Literal& lit) const {
  return !(*this == lit);
}


void Literal::write(ostream& os, bool withTypes) const {
  if (s->range_type == Symbol::binary  &&  TL::isZero(value))
    os << "-";
  os << s->name << "(";
  uint i;
  FOR1D(args, i) {
    if (args(i) < 10)
      os << variable__digit2char(args(i));
    else
      os << args(i);
    if (withTypes) {
      if (s->arg_types.N > 0) {
        os<<"/"<<s->arg_types(i)->name;
      }
    }
    if (i < args.N-1)
      os << " ";
  }
  os << ")";
  if (s->range_type != Symbol::binary) {
    switch(comparison_type) {
      case Literal::comparison_equal: os << "="; break;
      case Literal::comparison_less: os << "<"; break;
      case Literal::comparison_lessEqual: os << "<="; break;
      case Literal::comparison_greater: os << ">"; break;
      case Literal::comparison_greaterEqual: os << ">="; break;
      case Literal::comparison_offset: os << "+="; break;
      case Literal::comparison_variable: os << "=?"; break;
      default: HALT("Unknown comparison type")
    }
    if (comparison_type != Literal::comparison_variable)
      os << value;
  }
}


bool Literal::compareValueTo(Literal::ComparisonType compType, double a) {
  return compare(this->value, compType, a);
}


bool Literal::compare(double a, Literal::ComparisonType compType, double b) {
  switch(compType) {
    case Literal::comparison_equal: return TL::areEqual(a,b);
    case Literal::comparison_less: return a<b;
    case Literal::comparison_lessEqual: return a<=b;
    case Literal::comparison_greater: return a>b;
    case Literal::comparison_greaterEqual: return a>=b;
    default: HALT("Undefined comparison:  "<<compType)
  }
  return false;
}


bool Literal::isNegated() const {
  return s->range_type == Symbol::binary  &&  TL::isZero(value);
}


Literal* Literal::getNegated() {
  CHECK(s->range_type == Symbol::binary, "only defined for binary symbols")
  return Literal::get(s, args, (TL::isZero(value) ? 1. : 0.));
}


bool Literal::equivalent(const LitL& p1, const LitL& p2) {
  if (p1.N != p2.N)
    return false;
  uint i, j;
  // assumption: all predicate tuples in p1 are distinct
  FOR1D(p1, i) {
    FOR1D(p2, j) {
      if (*(p1(i))==*(p2(j)))
        break;
    }
    if (j== p2.N)
      return false;
  }
  return true;
}


bool Literal::nonContradicting(const LitL& l1, const LitL& l2) {
  uint i, k;
  FOR1D(l1, i) {
    FOR1D(l2, k) {
      if (l1(i)->s == l2(k)->s  &&  l1(i)->args == l2(k)->args) {
        if (!TL::areEqual(l1(i)->value, l2(k)->value))
          return false;
      }
    }
  }
  return true;
}

bool Literal::typeCheck() {
  if (s->arg_types.N == 0) return true;   //no types specified
  CHECK(args.N == s->arg_types.N, "Unexpected nmber of arguments!");
  uint i;
  FOR1D(args, i) {
    if (!s->arg_types(i)->subsumes(*reason::getArgumentTypeOfObject(args(i))))
      return false;
  }
  return true;
}


void Literal::getLiterals(LitL& lits, Symbol* s, const uintA& constants, double value, bool withRepeatingArguments) {
  lits.clear();
  MT::Array< uintA > args_lists;
  TL::allPermutations(args_lists, constants, s->arity, withRepeatingArguments, true);
  uint i;
  FOR1D(args_lists, i) {
    lits.append(Literal::get(s, args_lists(i), value));
  }
}


void Literal::getLiterals_state(LitL& lits, const uintA& constants, double value, bool binaryOnly) {
  lits.clear();
  SymL syms_state;
  Symbol::get_state(syms_state);
  uint i;
  FOR1D(syms_state, i) {
    if (syms_state(i)->range_type != Symbol::binary) continue;
    LitL lits_local;
    getLiterals(lits_local, syms_state(i), constants, value);
    lits.setAppend(lits_local);
  }
}


void Literal::getLiterals_state(LitL& lits, const uintA& constants, const uintA& constants_mustBeContained, double value, bool binaryOnly) {
  uint DEBUG=0;
  if (DEBUG>0) cout<<"getLiterals_state [START]"<<endl;
  lits.clear();
  LitL lits_total;
  Literal::getLiterals_state(lits_total, constants, value, binaryOnly);
  uint i;
  FOR1D(lits_total, i) {
    if (numberSharedElements(lits_total(i)->args, constants_mustBeContained) == constants_mustBeContained.N) {
      lits.append(lits_total(i));
    }
  }
  if (DEBUG>1) relational::write(lits);
  if (DEBUG>0) cout<<"getLiterals_state [END]"<<endl;
}


void Literal::getLiterals_actions(LitL& lits, const uintA& arguments) {
  lits.clear();
  SymL syms_action;
  Symbol::get_action(syms_action);
  uint i;
  FOR1D(syms_action, i) {
    LitL lits2;
    getLiterals(lits2, syms_action(i), arguments, 1.);
    lits.setAppend(lits2);
  }
}


double Literal::getValue(const LitL& literals, const Symbol* s, const uintA& args) {
  uint i;
  FOR1D(literals, i) {
    if (literals(i)->s == s  &&  literals(i)->args == args)
      return literals(i)->value;
  }
  HALT("value cannot be determined:  s="<<*s<<"  args="<<args<<"   in  literals="<<literals);
}


uint Literal::getArguments(uintA& args, const LitL& lits) {
  uint i;
  FOR1D(lits, i) {
    args.setAppend(lits(i)->args);
  }
  TL::sort_asc(args);
  return args.N;
}


bool Literal::negativeBinaryLiteralsLast(const LitL& lits) {
  bool negative_binary_started = false;
  uint i;
  FOR1D(lits, i) {
    bool is_positive = !lits(i)->isNegated();
    if (is_positive && negative_binary_started)
      return false;
    if (!is_positive && !negative_binary_started)
      negative_binary_started = true;
  }
  return true;
}


// positives first
void Literal::sort(LitL& lits) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"Literal::sort [START]"<<endl;}
  if (DEBUG>0) {PRINT(lits);}
  uint i;
  uintA keys(lits.N);
  SymL symbols;
  Symbol::get(symbols);
  // literals with small keys first
  // positions 1-4: arguments
  // positions 5-6: symbol
  // position 7: non-binary
  // position 8: negative
  FOR1D(lits, i) {
    uint key = 0;
    // positions 1-4: arguments
    if (lits(i)->s->arity > 0) {
      CHECK(lits(i)->args(0) < 100, "");
      key += 10e2 * lits(i)->args(0);
    }
    if (lits(i)->s->arity > 1) {
      CHECK(lits(i)->args(1) < 100, "");
      key += lits(i)->args(1);
    }
    // ATTENTION!!! When sorting, we only consider the first two arguments!
    
    // positions 5-6: symbols
    key += 10e4 * symbols.findValue(lits(i)->s);

    // position 7: non-binary?
    if (lits(i)->s->range_type != Symbol::binary) {
      key += 10e7;
    }

    // position 8: negated literals last
    if (lits(i)->isNegated()) {
      key += 10e8;
    }
    keys(i) = key;
    if (DEBUG>1) {cout<<"("<<i<<") "<<*lits(i)<<"\t"<<keys(i)<<endl;}
  }
  uintA sortedIndices;
//   PRINT(sortedIndices);
  TL::sort_asc_keys(sortedIndices, keys);
//   PRINT(sortedIndices);
  LitL lits_sorted;
  FOR1D(sortedIndices, i) {
    lits_sorted.append(lits(sortedIndices(i)));
  }
  if (DEBUG>0) {PRINT(lits_sorted);}
  lits = lits_sorted;
  if (DEBUG>0) {cout<<"Literal::sort [END]"<<endl;}
}


int Literal::findPattern(const LitL& actions, uint minRepeats) {
  uint length, repeat, pos;
  int max_repeat_length = 0;
  for (length=1; length<=actions.N / minRepeats; length++) {
    bool successful_repeat = true;
    for (repeat=1; repeat<minRepeats && successful_repeat; repeat++) {
      for (pos=0; pos<length && successful_repeat; pos++) {
        if (actions(repeat * length + pos) != actions(pos))
          successful_repeat = false;
      }
    }
    if (successful_repeat)
      max_repeat_length = length;
  }
  return max_repeat_length;
}


Literal* Literal::getLiteral_default_action() {
  uintA empty;
  return get(Symbol::getDefaultAction(), empty, 1.);
}


Literal* l_doNothing = NULL;
Literal* Literal::getLiteral_doNothing() {
  if (l_doNothing == NULL) {
    uintA empty;
    l_doNothing = get(Symbol::get(MT::String("doNothing"), 0, Symbol::action), empty, 1.);
  }
  return l_doNothing;
}


void write(const LitL& lits, ostream& os) {
  uint k;
  FOR1D(lits, k) {
    if (lits(k)!=NULL)
      lits(k)->write(os);
    else
      os<<"NULL";
    os << " ";
  }
}


void write(const MT::Array< LitL >& outcomes, ostream& os) {
  uint k;
  FOR1D(outcomes, k) {
    os << "(" << k << ") ";
    write(outcomes(k), os);
    os << endl;
  }
}



/************************************************
 * 
 *     SymbolicState
 * 
 ************************************************/

SymbolicState::SymbolicState() {
  derived_lits_are_calculated = false;
}


SymbolicState::SymbolicState(const MT::Array<Literal*>& _lits) {
  derived_lits_are_calculated = false;
  this->lits = _lits;
  reason::derive(this);
  Literal::getArguments(state_constants, lits);
}


void SymbolicState::write(ostream& os, bool primOnly) const {
  uint i;
  FOR1D(lits, i) {
    if (primOnly && lits(i)->s->symbol_type != Symbol::primitive)
      continue;
    os << *lits(i) << " ";
  }
}

std::string SymbolicState::toString(bool primOnly) const {
  std::ostringstream ss;
  this->write(ss, primOnly);
  return ss.str();
}


void SymbolicState::read(istream& in, bool read_constants) {
  if (read_constants) {
    in >> state_constants;
  }
  MT::String line;
  line.read(in, NULL, "\n");
  Literal::get(lits, line);
  if (!read_constants) {
    Literal::getArguments(state_constants, lits);
  }
  lits.memMove = true;
  uint i;
  FOR1D_DOWN(lits, i) {
    if (lits(i)->s->symbol_type != Symbol::primitive)
      lits.remove(i);
  }
  derived_lits_are_calculated = false;
  reason::derive(this);
//  write(cout, false);
//  cout << endl;
}


bool SymbolicState::operator==(const SymbolicState& s) const {
  LitL this__lits_prim, other__lits_prim;
  uint i;
  FOR1D(lits, i) {
    if (lits(i)->s->symbol_type == Symbol::primitive)
      this__lits_prim.append(lits(i));
  }
  FOR1D(s.lits, i) {
    if (s.lits(i)->s->symbol_type == Symbol::primitive)
      other__lits_prim.append(s.lits(i));
  }
  return Literal::equivalent(this__lits_prim, other__lits_prim);
}


bool SymbolicState::operator!=(const SymbolicState& s) const {
  return !(*this == s);
}


SymbolicState::~SymbolicState() {
  //  object deletion of concept-instances6 is managed ny the LogicEngine now!
//     uint i;
//     FOR1D(lits_prim, i) {
//         delete lits_prim(i);
//     }
//     FOR1D(lits_derived, i) {
//         delete lits_derived(i);
//     }
//     FOR1D(lits_comp, i) {
//         delete lits_comp(i);
//     }
//     FOR1D(fv_prim, i) {
//         delete fv_prim(i);
//     }
//     FOR1D(fv_derived, i) {
//         delete fv_derived(i);
//     }
}


void SymbolicState::calcDifferences(LitL& lits_diff_1to2, LitL& lits_diff_2to1, uintA& changedConstants, const SymbolicState& state1, const SymbolicState& state2) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcDifferences [START]"<<endl;}
  
  if (DEBUG>0) {
    cout<<"SymbolicState 1:"<<endl;  state1.write();  cout<<endl;
    cout<<"SymbolicState 2:"<<endl;  state2.write();  cout<<endl;
  }
  
  lits_diff_1to2.clear();
  lits_diff_2to1.clear();
  changedConstants.clear();
  
  uint i;
  FOR1D(state1.lits, i) {
    if (state2.lits.findValue(state1.lits(i)) < 0) {
      lits_diff_1to2.append(state1.lits(i));
      changedConstants.setAppend(state1.lits(i)->args);
    }
  }
  FOR1D(state2.lits, i) {
    if (state1.lits.findValue(state2.lits(i)) < 0) {
      lits_diff_2to1.append(state2.lits(i));
      changedConstants.setAppend(state2.lits(i)->args);
    }
  }
  
  if (DEBUG>0) {
    cout<<"lits_diff_1to2:    "<<lits_diff_1to2<<endl;
    cout<<"lits_diff_2to1:    "<<lits_diff_2to1<<endl;
    PRINT(changedConstants);
  }
  
  if (DEBUG>0) {cout<<"calcDifferences [END]"<<endl;}
}


void SymbolicState::filterState_full(SymbolicState& s_filtered, const SymbolicState& s_full, const uintA& filter_constants, bool primOnly) {
  uint i, k;
  // lits
  s_filtered.lits.clear();
  FOR1D(s_full.lits, i) {
    FOR1D(s_full.lits(i)->args, k) {
      if (filter_constants.findValue(s_full.lits(i)->args(k))<0)
        break;
    }
    if (k==s_full.lits(i)->args.N)
      s_filtered.lits.append(s_full.lits(i));
  }
}


void SymbolicState::filterState_atleastOne(SymbolicState& s_filtered, const SymbolicState& s_full, const uintA& filter_constants, bool primOnly) {
  uint i, k;
  // lits
  s_filtered.lits.clear();
  FOR1D(s_full.lits, i) {
    if (s_full.lits(i)->args.N == 0)
      s_filtered.lits.append(s_full.lits(i));
    FOR1D(s_full.lits(i)->args, k) {
      if (filter_constants.findValue(s_full.lits(i)->args(k))>=0) {
        s_filtered.lits.append(s_full.lits(i));
        break;
      }
    }
  }
}


uint SymbolicState::getArgument(const SymbolicState& state, const Symbol& s) {
  CHECK(s.arity == 1, "");
  uintA args;
  getArguments(args, state, s);
  if (args.N == 1)
    return args(0);
  else
    return TL::UINT_NIL;
}


void SymbolicState::getArguments(uintA& args, const SymbolicState& state, const Symbol& s) {
  args.clear();
  uint i;
  FOR1D(state.lits, i) {
    if (*state.lits(i)->s == s)
      args.setAppend(state.lits(i)->args);
  }
  TL::sort_asc(args);
}

void SymbolicState::getArguments(uintA& args, const SymbolicState& state) {
  args.clear();
  uint i;
  FOR1D(state.lits, i) {
    args.setAppend(state.lits(i)->args);
  }
  TL::sort_asc(args);
}

double SymbolicState::getValue(const Symbol* s, const SymbolicState& state) {
  if (s->arity != 0) HALT("only defined for 0-ary symbol");
  uintA empty;
  return Literal::getValue(state.lits, s, empty);
}


void SymbolicState::getValues(arr& values, const SymbolicState& state, const Symbol& s, const uintA& objs) {
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s != &s) continue;
    if (numberSharedElements(state.lits(i)->args, objs) == objs.N)
      values.append(state.lits(i)->value);
  }
}


void SymbolicState::getRelatedConstants(uintA& constants_related, uint id, bool id_covers_first, const Symbol& s, const SymbolicState& state) {
  constants_related.clear();
  uint i;
  FOR1D(state.lits, i) {
    if (*state.lits(i)->s == s) {
      if (id_covers_first) {
        if (state.lits(i)->args(0) == id)
          constants_related.append(state.lits(i)->args(1));
      }
      else {
        if (state.lits(i)->args(1) == id)
          constants_related.append(state.lits(i)->args(0));
      }
    }
  }
}






/************************************************
 * 
 *     StateTransition
 * 
 ************************************************/


StateTransition::StateTransition(const SymbolicState& pre, Literal* action, const SymbolicState& post, double reward) {
  this->pre = pre;
  this->action = action;
  this->post = post;
  this->reward = reward;
  calcChanges();
}


StateTransition::StateTransition() {}


StateTransition::~StateTransition() {}


void StateTransition::calcChanges() {
  // only look at primitives
  changedConstants.clear();
  changes.clear();
  uint i,j;
  //changed values from pre to post
  FOR1D(pre.lits, i) {
    if (pre.lits(i)->s->symbol_type != Symbol::primitive) continue;
    if (post.lits.findValue(pre.lits(i)) < 0) {     //Literal has changed
      changedConstants.setAppend(pre.lits(i)->args);
      if (pre.lits(i)->s->range_type == Symbol::binary)
        changes.append(pre.lits(i)->getNegated());      //predicate is negated in post state
      else {    //search for symbol with same arguments
        bool found = false;
        FOR1D(post.lits, j) {
          if (post.lits(j)->s == pre.lits(i)->s && post.lits(j)->args == pre.lits(i)->args) {
              changes.append(post.lits(j));    
              found = true;
              break;
          }
        }
        if (!found) { cout << "Warning: Undefined value for non-binary symbol: "; pre.lits(i)->s->write(); cout << endl; }
      }
    }
  }
  // pre-, post+
  //only consider binary symbols because we assume that other symbols are not added or removed but only change values
  FOR1D(post.lits, i) {
    if (post.lits(i)->s->symbol_type != Symbol::primitive || post.lits(i)->s->range_type != Symbol::binary) continue;
    if (pre.lits.findValue(post.lits(i)) < 0) {
      changes.append(post.lits(i));
      changedConstants.setAppend(post.lits(i)->args);
    }
  }
}


void StateTransition::write(ostream& os, uint detailLevel) const {
  os << "ACTION: " << *action << endl;
  if (detailLevel>0) {
    os << "PRE:    ";
    if (detailLevel == 1) this->pre.write(os, true);
    else this->pre.write(os, false);
    os << endl;
    os << "POST:   ";
    if (detailLevel == 1) this->post.write(os, true);
    else this->post.write(os, false);
    os<<endl;
    os<<"Changed constants: "<<changedConstants<<endl;
    os << "Diff: "<< changes.N << endl;
  }
  os<<reward<<endl;
  os<<"Changed literals ("<<changes.N<<"):    " << changes << endl;
}


bool StateTransition::noChange() {
  return changedConstants.N == 0  && changes.N == 0;
}

/*
void StateTransition::getConstants(uintA& constants) {
  constants.clear();

  if (pre.state_constants.N == 0)
    Literal::getArguments(pre.state_constants, pre.lits);
  constants.setAppend(pre.state_constants);

  if (post.state_constants.N == 0)
    Literal::getArguments(post.state_constants, post.lits);
  constants.setAppend(post.state_constants);

  constants.setAppend(action->args);
}
*/


void StateTransition::write(const StateTransitionL& exs, ostream& os, uint detailLevel) {
  uint i;
  FOR1D(exs, i) {
    os<<"# i="<<i<<" -----"<<endl;
    exs(i)->write(os, detailLevel);
  }
}


StateTransitionL& StateTransition::read_SASAS(const char* filename) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"StateTransition::read_SASAS [START]"<<endl;}
  if (DEBUG>0) {PRINT(filename);}
  static StateTransitionL experiences;
  ifstream in(filename);
  if (!in.is_open()) {HALT("File not found: "<<filename);}
  MT::skip(in);
  bool read_state = true;
  SymbolicState* state, *state_old = NULL;
  Literal* action = NULL;
  while (MT::skip(in) != -1) {
    if (DEBUG>2) {PRINT(MT::peerNextChar(in));}
    if (read_state) {
      state = new SymbolicState;
      state->read(in);
      if (state_old != NULL) {
        experiences.append(new StateTransition(*state_old, action, *state));
        if (DEBUG>0) {cout<<*experiences.last()<<endl;}
        if (!reason::isGround(*experiences.last())) HALT("Attention! State transition is not ground!"<<endl<<*experiences.last());
      }
      state_old = state;
    }
    else {
      MT::String line;
      line.read(in, NULL, "\n");
      if (DEBUG>1) {cout<<"READING ACTION:"<<endl; PRINT(line);}
      action = Literal::get(line);
    }
    read_state = !read_state;
  }
  if (DEBUG>0) {PRINT(experiences.N);}
  if (DEBUG>0) {cout<<"StateTransition::read_SASAS [END]"<<endl;}
  return experiences;
}


StateTransitionL StateTransition::read_SAS_SAS(const char* filename) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"StateTransition::read_SAS_SAS [START]"<<endl;}
  if (DEBUG>0) {PRINT(filename);}
  //static StateTransitionL experiences; // are you kidding me?
  StateTransitionL experiences;
  experiences.clear(); // bugfix by Sebastian
  ifstream in(filename);
  if (!in.is_open()) {HALT("File not found: "<<filename);}
  MT::skip(in);
  SymbolicState* state_pre, *state_post = NULL;
  Literal* action = NULL;
  while (MT::skip(in) != -1) {
    if (DEBUG>2) {PRINT(MT::peerNextChar(in));}
    state_pre = new SymbolicState;
    state_pre->read(in);
    MT::String line;
    line.read(in, NULL, "\n");
    // additionally remove trailing comments
    std::stringstream ss; ss << line;
    line.read(ss, NULL, "#");
    if (DEBUG>1) {cout<<"READING ACTION:"<<endl; PRINT(line);}
    action = Literal::get(line);
    state_post = new SymbolicState;
    state_post->read(in);
    double reward = 0.;
    char c = MT::skip(in);
    if (isdigit(c) || c == '-') {
      in >> reward;
    }
    experiences.append(new StateTransition(*state_pre, action, *state_post, reward));
    if (DEBUG>0) {cout<<"Read StateTransition #"<<experiences.N-1<<":"<<endl; experiences.last()->write(cout,0);}
    if (!reason::isGround(*experiences.last())) HALT("Attention! State transition is not ground!"<<endl<<*experiences.last());
  }
  if (DEBUG>0) {PRINT(experiences.N);}
  if (DEBUG>0) {cout<<"StateTransition::read_SAS_SAS [END]"<<endl;}
  return experiences;
}

bool StateTransition::operator==(const StateTransition& s) const {
  return this->reward == s.reward &&
      this->pre == s.pre && *(this->action) == *(s.action) &&
      this->post == s.post;
}

bool StateTransition::operator!=(const StateTransition& s) const {
  return ! (*this == s);
}

// void write(const StateTransitionL& exs, ostream& os) {
//   uint k;
//   FOR1D(exs, k) {
//     os << "[" << k << "] ("<<(k+1)<<"/"<<exs.N<<")"<<endl;
//     exs(k)->write(os);
//   }
// }


void getConstants(uintA& constants, const StateTransitionL& transitions) {
  constants.clear();
  uint i;
  FOR1D(transitions, i) {
    constants.setAppend(transitions(i)->pre.state_constants);
    constants.setAppend(transitions(i)->post.state_constants);
  }
}

void getStateSymbols(SymL& symbols, const StateTransitionL& transitions) {
  symbols.clear();
  uint i;
  FOR1D(transitions, i) {
    const StateTransition* s = transitions(i);
    uint j;
    FOR1D(s->pre.lits, j) {
      symbols.setAppend(s->pre.lits.elem(j)->s);
    }
    FOR1D(s->post.lits, j) {
      symbols.setAppend(s->post.lits.elem(j)->s);
    }
  }
}


}  // namespace PRADA



std::ostream& operator<<(std::ostream& os, const relational::SymbolicState& s) {
  s.write(os); return os;
}


std::ostream& operator<<(std::ostream& os, const relational::Literal& l) {
  l.write(os); return os;
}


std::ostream& operator<<(std::ostream& os, const relational::StateTransition& e) {
  e.write(os); return os;
}


std::ostream& operator<<(std::ostream& os, const relational::LitL& lits) {
  write(lits, os);
  return os;
}


std::ostream& operator<<(std::ostream& os, const relational::StateTransitionL& sl) {
  relational::StateTransition::write(sl, os, 0);
  return os;
}
