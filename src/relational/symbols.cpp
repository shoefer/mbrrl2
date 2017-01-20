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

#include <stdlib.h>
#include "symbols.h"
#include "relational/literals.h"


namespace relational {

  
/************************************************
 * 
 *     Symbol
 * 
 ************************************************/

MT::Array< Symbol* > Symbol::mem__all_symbols;
bool Symbol::symbols_sorted = false;

Symbol::Symbol() {
  mem__all_symbols.append(this);
  symbols_sorted = false;
}


Symbol* Symbol::get(const MT::String& name, uint arity, SymbolType symbol_type, RangeType range_type) {
  Symbol* potential_s = NULL;
  // Special treatment for "default" symbol
  if (name == "default") {
    uint i;
    FOR1D(mem__all_symbols, i) {
      if (mem__all_symbols(i)->name == name)
        return mem__all_symbols(i);
    }
    printf("Default action symbol does not exist\n");
  }
  else
    potential_s = get(name);
  if (potential_s != NULL) {
    if ( potential_s->arity != arity
      || potential_s->symbol_type != symbol_type
      || potential_s->range_type != range_type) {
      HALT("Symbol "<<name<<" already exists with different properties: "<<*potential_s);
    }
    else
      return potential_s;
  }
  Symbol* s = new Symbol;
  s->name = name;
  s->arity = arity;
  s->symbol_type = symbol_type;
  s->range_type = range_type;
  if (s->range_type == integer_set) {
    // default initialization of integer_set
    s->range.clear();
    uint i;
    for (i=0; i<=5; i++) s->range.append(i);
  }
  return s;
}


Symbol* Symbol::get(const MT::String& name) {
  if (name == "default") {
    return getDefaultAction();
  }
  uint i;
  FOR1D(mem__all_symbols, i) {
    if (mem__all_symbols(i)->name == name)
      return mem__all_symbols(i);
  }
  return NULL;
}


Symbol* Symbol::get(const char* name) {
  return get(MT::String(name));
}

Symbol* Symbol::getDefaultAction() {
  return Symbol::get(MT::String("default"), 0, Symbol::action);
}

void Symbol::getDefiningSymbols(SymL& symbols, bool only_direct_precessors) const {
  if (symbol_type == primitive)
    symbols.clear();
  else
    HALT("should be overwritten "<<*this);
}


void Symbol::write(ostream& os) const {
  os << name << " " << arity << " ";
  switch (symbol_type) {
    case action: os << "action"; break;
    case primitive: os << "primitive"; break;
    case conjunction: os << "conjunction"; break;
    case transclosure: os << "transclosure"; break;
    case count: os << "count"; break;
//     case function_count: os << "function_count"; break;
    case avg: os << "avg"; break;
    case max: os << "max"; break;
    case function_change: os << "function_change"; break;
    case sum: os << "sum"; break;
    case function_reward: os << "function_reward"; break;
    case function_difference: os << "function_difference"; break;
    default: NIY;
  }
  os << " ";
    switch (range_type) {
    case binary: os << "binary"; break;
    case integer_set: os << "integer_set"; break;
    case integers: os << "integers"; break;
    case reals: os << "reals"; break;
    default: HALT("wrong range_type="<<range_type<<" for "<<name); NIY;
  }
  os << " ";
  uint i;
  FOR1D(arg_types, i) {
    os << " " << arg_types(i)->type;
  }
}


bool Symbol::operator==(const Symbol& p) const {
  // not via id!!
  if (this->arity != p.arity  ||  this->symbol_type != p.symbol_type)
      return false;
  return (MT::String) this->name == (MT::String) p.name;
}


bool Symbol::operator!=(const Symbol& s) const {
  return !(*this == s);
}


void Symbol::get(SymL& all_symbols) {
  if (!symbols_sorted) {
    sort(mem__all_symbols);
    symbols_sorted = true;
  }
  all_symbols = mem__all_symbols;
}


void Symbol::get_state(SymL& symbols) {
  symbols.clear();
  SymL all_symbols;
  get(all_symbols);
  uint i;
  FOR1D(all_symbols, i) {
    if (all_symbols(i)->symbol_type != action)
      symbols.append(all_symbols(i));
  }
}


void Symbol::get_state_nonBinary(SymL& symbols) {
  SymL all_symbols;
  get_state(all_symbols);
  uint i;
  FOR1D(all_symbols, i) {
    if (all_symbols(i)->range_type != Symbol::binary)
      symbols.append(all_symbols(i));
  }
}


void Symbol::get_action(SymL& symbols) {
  symbols.clear();
  SymL all_symbols;
  get(all_symbols);
  uint i;
  FOR1D(all_symbols, i) {
    if (all_symbols(i)->symbol_type == action)
      symbols.append(all_symbols(i));
  }
}


void Symbol::sort(SymL& symbols) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"Symbol::sort [START]"<<endl;}
  if (DEBUG>0) {PRINT(symbols);}
  if (symbols.N==0) return;
  SymL sorted_symbols;
  uint i, k;
  FOR1D(symbols, i) {
    if (symbols(i)->symbol_type == Symbol::action)
      sorted_symbols.append(symbols(i));
  }
  FOR1D(symbols, i) {
    if (symbols(i)->symbol_type == Symbol::primitive  &&  symbols(i)->range_type == Symbol::binary)
      sorted_symbols.append(symbols(i));
  }
  FOR1D(symbols, i) {
    if (symbols(i)->symbol_type == Symbol::primitive  &&  symbols(i)->range_type != Symbol::binary)
      sorted_symbols.append(symbols(i));
  }
  SymL derived_symbols;
  FOR1D(symbols, i) {
    if (sorted_symbols.findValue(symbols(i)) < 0)
      derived_symbols.append(symbols(i));
  }
  while (sorted_symbols.N < symbols.N) {
    bool one_added = false;
    FOR1D(derived_symbols, i) {
      if (derived_symbols(i) != NULL) {
        if (DEBUG>0) {cout<<"trying to add "<<*derived_symbols(i)<<endl;}
        SymL defining_symbols;
        derived_symbols(i)->getDefiningSymbols(defining_symbols);
        if (DEBUG>0) {PRINT(defining_symbols);}
        bool add = true;
        FOR1D(defining_symbols, k) {
          if (!symbols.contains(defining_symbols(k))) {
            MT_MSG("WARNING: defining symbol "<<defining_symbols(k)->name<<" not among symbols to be sorted");
            continue;
          }
          else if (!sorted_symbols.contains(defining_symbols(k))) {
            add = false;
            break;
          }
        }
        if (DEBUG>0) {cout<<""<<(add ? "ADD" : "ignored")<<endl;}
        if (add) {
          sorted_symbols.append(derived_symbols(i));
          derived_symbols(i) = NULL;
          one_added = true;
        }
      }
    }
    if (!one_added) HALT("Symbol definitions are cyclic! "<<symbols);
  }
  symbols = sorted_symbols;
  if (DEBUG>0) {PRINT(symbols);}
  if (DEBUG>0) {cout<<"Symbol::sort [END]"<<endl;}
}



/************************************************
 * 
 *     ConjunctionSymbol
 * 
 ************************************************/


ConjunctionSymbol* ConjunctionSymbol::get(const MT::String& name, uint arity, MT::Array<Literal* > base_literals, bool free_vars_all_quantified) {
  Symbol* potential_s = Symbol::get(name);
  if (potential_s != NULL) {
    if ( potential_s->arity != arity
      || potential_s->symbol_type != conjunction
      || potential_s->range_type != binary
      || ((ConjunctionSymbol*) potential_s)->base_literals != base_literals
      || ((ConjunctionSymbol*) potential_s)->free_vars_all_quantified != free_vars_all_quantified
    ) {
      HALT("Symbol "<<name<<" already exists with different properties: "<<*potential_s);
    }
    else
      return (ConjunctionSymbol*) potential_s;
  }
  ConjunctionSymbol* s = new ConjunctionSymbol;
  s->name = name;
  s->arity = arity;
  s->base_literals = base_literals;
  s->free_vars_all_quantified = free_vars_all_quantified;
  
  uint i, k;
  s->redundant__base_literal_with_free_vars.resize(base_literals.N);
  s->redundant__base_literal_with_free_vars.setUni(false);
  FOR1D(s->base_literals, i) {
    FOR1D(s->base_literals(i)->args, k) {
      if (s->base_literals(i)->args(k) >= arity) {
        s->redundant__base_literal_with_free_vars(i) = true;
        break;
      }
    }
  }
  return s;
}


void ConjunctionSymbol::write(ostream& os) const {
  Symbol::write(os);
  uint i;
  os << " <-- ";
  uintA freeVars;
  this->getFreeVars(freeVars);
  if (freeVars.N>0) {
    if (!free_vars_all_quantified)
      os << "Ex ";
    else
      os << "All ";
    FOR1D(freeVars, i) {
      if (freeVars(i) == 0) os << "X";
      else if (freeVars(i) == 1) os << "Y";
      else if (freeVars(i) == 2) os << "Z";
      else if (freeVars(i) == 3) os << "V";
      else if (freeVars(i) == 4) os << "W";
      else if (freeVars(i) == 5) os << "U";
      if (freeVars.N-1 != i)
        os<<freeVars(i)<<" ";
    }
    os <<" ";
    os <<" ";
  }
  relational::write(base_literals, os);
}


void ConjunctionSymbol::getFreeVars(uintA& freeVars) const {
  uintA freeVars_pos, freeVars_neg;
  getFreeVars(freeVars_pos, freeVars_neg);
  freeVars.append(freeVars_pos);
  freeVars.append(freeVars_neg);
  TL::sort_asc(freeVars);
}


void ConjunctionSymbol::getFreeVars(uintA& freeVars_pos, uintA& freeVars_neg) const {
  uint i, k;
  FOR1D(base_literals, i) {
    if (!base_literals(i)->isNegated()) {
      FOR1D(base_literals(i)->args, k) {
        if (base_literals(i)->args(k) >= arity)
          freeVars_pos.setAppend(base_literals(i)->args(k));
      }
    }
  }
  FOR1D(base_literals, i) {
    if(base_literals(i)->isNegated()) {
      FOR1D(base_literals(i)->args, k) {
        if (base_literals(i)->args(k) >= arity  &&  freeVars_pos.findValue(base_literals(i)->args(k)) < 0)
          freeVars_neg.setAppend(base_literals(i)->args(k));
      }
    }
  }
}


void ConjunctionSymbol::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  uint i;
  FOR1D(base_literals, i) {
    if (!only_direct_predecessors) {
      SymL local_symbols;
      base_literals(i)->s->getDefiningSymbols(local_symbols, only_direct_predecessors);
      symbols.setAppend(local_symbols);
    }
    symbols.setAppend(base_literals(i)->s);
  }
}



/************************************************
 * 
 *     TransClosureSymbol
 * 
 ************************************************/


TransClosureSymbol::TransClosureSymbol() {symbol_type = transclosure;}


TransClosureSymbol* TransClosureSymbol::get(const MT::String& name, uint arity, Symbol* base_symbol) {
  Symbol* potential_s = Symbol::get(name);
  if (potential_s != NULL) {
    if ( potential_s->arity != arity
      || potential_s->symbol_type != transclosure
      || potential_s->range_type != binary
      || ((TransClosureSymbol*) potential_s)->base_symbol != base_symbol
    ) {
      HALT("Symbol "<<name<<" already exists with different properties: "<<*potential_s);
    }
    else
      return (TransClosureSymbol*) potential_s;
  }
  TransClosureSymbol* s = new TransClosureSymbol;
  s->name = name;
  s->arity = arity;
  s->base_symbol = base_symbol;
  return s;
}


void TransClosureSymbol::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  if (!only_direct_predecessors) {
    SymL local_symbols;
    base_symbol->getDefiningSymbols(local_symbols, only_direct_predecessors);
    symbols.setAppend(local_symbols);
  }
  symbols.setAppend(base_symbol);
}


void TransClosureSymbol::write(ostream& os) const {
  Symbol::write(os);
  os << " <-- + ";
  os << base_symbol->name;
}



/************************************************
 * 
 *     CountSymbol
 * 
 ************************************************/


CountSymbol::CountSymbol () {
  symbol_type = count;
  range_type = integers;  
}


CountSymbol* CountSymbol::get(const MT::String& name, uint arity, Literal* base_literal) {
  Symbol* potential_s = Symbol::get(name);
  if (potential_s != NULL) {
    if ( potential_s->arity != arity
      || potential_s->symbol_type != count
      || potential_s->range_type != integers
      || ((CountSymbol*) potential_s)->base_literal != base_literal
    ) {
      HALT("Symbol "<<name<<" already exists with different properties: "<<*potential_s);
    }
    else
      return (CountSymbol*) potential_s;
  }
  CountSymbol* s = new CountSymbol;
  s->name = name;
  s->arity = arity;
  s->base_literal = base_literal;
  return s;
}


void CountSymbol::write(ostream& os) const {
  Symbol::write(os);
  os << " <-- Num ";
  uintA freeVars;
  uint i;
  FOR1D(base_literal->args, i) {
    if (base_literal->args(i) >= arity) {
      if (base_literal->args(i) == 0) os << "X";
      else if (base_literal->args(i) == 1) os << "Y";
      else if (base_literal->args(i) == 2) os << "Z";
      else os<<base_literal->args(i)<<" ";
      os <<" ";
    }
  }
  os << " ";
  os<<*base_literal;
}


void CountSymbol::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  if (!only_direct_predecessors) {
    SymL local_symbols;
    base_literal->s->getDefiningSymbols(local_symbols, only_direct_predecessors);
    symbols.setAppend(local_symbols);
  }
  symbols.setAppend(base_literal->s);
}



/************************************************
 * 
 *     AverageFunction
 * 
 ************************************************/

AverageFunction::AverageFunction() {
  symbol_type = avg;
  range_type = reals;
}


void AverageFunction::write(ostream& os) const {
  NIY;
  Symbol::write(os);
//   os << " " << base_symbol->id;
}


void AverageFunction::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  if (!only_direct_predecessors) {
    SymL local_symbols;
    base_symbol->getDefiningSymbols(local_symbols, only_direct_predecessors);
    symbols.setAppend(local_symbols);
  }
  symbols.setAppend(base_symbol);
}



/************************************************
 * 
 *     SumFunction
 * 
 ************************************************/

SumFunction::SumFunction() {
  symbol_type = sum;
  range_type = integers;
}


SumFunction* SumFunction::get(const MT::String& name, uint arity, Symbol* base_symbol) {
  Symbol* potential_s = Symbol::get(name);
  if (potential_s != NULL) {
    if ( potential_s->arity != arity
      || potential_s->symbol_type != sum
      || potential_s->range_type != reals
      || ((SumFunction*) potential_s)->base_symbol != base_symbol
    ) {
      HALT("Symbol "<<name<<" already exists with different properties: "<<*potential_s);
    }
    else
      return (SumFunction*) potential_s;
  }
  SumFunction* s = new SumFunction;
  s->name = name;
  s->arity = arity;
  s->base_symbol = base_symbol;
  return s;
}


void SumFunction::write(ostream& os) const {
  Symbol::write(os);
  os << " <-- ";
  uint i;
  for(i=0;i<arity;i++) {
    if (i == arity-1)
      os << i;
    else
      os << i << " ";
  }
  os<<" Sum "<<base_symbol->name;
}


void SumFunction::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  if (!only_direct_predecessors) {
    SymL local_symbols;
    base_symbol->getDefiningSymbols(local_symbols, only_direct_predecessors);
    symbols.setAppend(local_symbols);
  }
  symbols.setAppend(base_symbol);
}



/************************************************
 * 
 *     MaxFunction
 * 
 ************************************************/


MaxFunction::MaxFunction() {
  symbol_type = max;
  range_type = integers;
}

void MaxFunction::write(ostream& os) const {
  NIY;
  Symbol::write(os);
//   os << " " << base_symbol->id;
}


void MaxFunction::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  if (!only_direct_predecessors) {
    SymL local_symbols;
    base_symbol->getDefiningSymbols(local_symbols, only_direct_predecessors);
    symbols.setAppend(local_symbols);
  }
  symbols.setAppend(base_symbol);
}



/************************************************
 * 
 *     ChangeFunction
 * 
 ************************************************/

ChangeFunction::ChangeFunction() {
  symbol_type = function_change;
}


void ChangeFunction::write(ostream& os) const {
  NIY;
  Symbol::write(os);
//   os << " " << base_symbol->id;
}


void ChangeFunction::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  if (!only_direct_predecessors) {
    SymL local_symbols;
    base_symbol->getDefiningSymbols(local_symbols, only_direct_predecessors);
    symbols.setAppend(local_symbols);
  }
  symbols.setAppend(base_symbol);
}



/************************************************
 * 
 *     RewardFunction
 * 
 ************************************************/

RewardFunction::RewardFunction() {
  symbol_type = function_reward;
}

void RewardFunction::write(ostream& os) const {
  NIY;
//   Symbol::write(os);
//   os << reward_value;
//   os<<" ( ";
//   uint i;
//   FOR1D(grounded_pis, i) {
//     os<<" ( ";
//     grounded_pis(i)->write_deprecated(os);
//     os<<" ) ";
//   }
//   os<<" )";
}

void RewardFunction::getDefiningSymbols(SymL & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  uint i;
  FOR1D(base_literals, i) {
    if (!only_direct_predecessors) {
      SymL local_symbols;
      base_literals(i)->s->getDefiningSymbols(local_symbols, only_direct_predecessors);
      symbols.setAppend(local_symbols);
    }
    symbols.setAppend(base_literals(i)->s);
  }
}



/************************************************
 * 
 *     DifferenceFunction
 * 
 ************************************************/

//f(X) - f(Y)
DifferenceFunction* DifferenceFunction::get(const MT::String& name, uint arity, Symbol* base_symbol, Literal *baseFunctionLit1, Literal *baseFunctionLit2, const MT::Array< Literal* > &restriction_lits) {
  Symbol* potential_s = Symbol::get(name);
  if (potential_s != NULL) {
    if ( potential_s->arity != arity
      || potential_s->symbol_type != function_difference
      || potential_s->range_type != integers
      || ((SumFunction*) potential_s)->base_symbol != base_symbol
    ) {
      HALT("Symbol "<<name<<" already exists with different properties: "<<*potential_s);
    }
    else
      return (DifferenceFunction*) potential_s;
  }
  DifferenceFunction* s = new DifferenceFunction;
  s->name = name;
  s->arity = arity;
  s->range_type = integers;
  s->baseSymbol = base_symbol;
  s->restrictionLits = restriction_lits;
  s->baseFunctionLit1 = baseFunctionLit1;
  s->baseFunctionLit2 = baseFunctionLit2;
  return s;
}

void DifferenceFunction::getFreeVars(uintA& freeVars) const {
  uint i, k;
  FOR1D(restrictionLits, i) {
    CHECK (!restrictionLits(i)->isNegated(), "Restriction Literals must be positive atm!")
      FOR1D(restrictionLits(i)->args, k) {
        if (restrictionLits(i)->args(k) >= arity)
          freeVars.setAppend(restrictionLits(i)->args(k));
    }
  }
}

void DifferenceFunction::write(ostream& os) const {
  Symbol::write(os);
  os << " <-- ";
  uint i;
  for(i=0;i<arity;i++) {
    if (i == arity-1)
      os << i;
    else
      os << i << " ";
  }
  os << baseSymbol->name << " " << restrictionLits;
}
void DifferenceFunction::getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors) const {
  symbols.clear();
  if (!only_direct_predecessors) {
    SymL local_symbols;
    baseSymbol->getDefiningSymbols(local_symbols, only_direct_predecessors);
    symbols.setAppend(local_symbols);

    for (int i = 0; i < restrictionLits.N; i++) {
      restrictionLits(i)->s->getDefiningSymbols(local_symbols, only_direct_predecessors);
      symbols.setAppend(local_symbols);
    }
  }
  symbols.setAppend(baseSymbol);
    for (int i = 0; i < restrictionLits.N; i++) {
      symbols.setAppend(restrictionLits(i)->s);
    }
}



/************************************************
 * 
 *     ArgumentType
 * 
 ************************************************/


ArgumentTypeL ArgumentType::mem_all_argument_types;


bool ArgumentType::subsumes(const ArgumentType& t) const {
  return this->name == t.name  ||  t.name == "any" || this->name == "any";
}


bool DisjunctionArgumentType::subsumes(const ArgumentType& t) const {
  bool subsumes = false;
  uint i;
  FOR1D(this->base_types, i) {
    if (this->base_types(i)->subsumes(t))
      subsumes = true;
  }
  return subsumes;
}


bool ArgumentType::operator==(const ArgumentType& t) const {
  return t.type == this->type;
}


bool ArgumentType::operator!=(const ArgumentType& t) const {
  return !(*this == t);
}


void ArgumentType::write(ostream& os) const {
  os << type << " " << name;
}

void DisjunctionArgumentType::write(ostream& os) const {
  ArgumentType::write(os);
  uint i;
  os << " = (";
  FOR1D(base_types, i) {
    if (i>0)
      os<<" OR";
    os << " "<<base_types(i)->name;
  }
  os << ")";
}


ArgumentType* ArgumentType::get(const MT::String& name) {
  uint i;
  FOR1D(mem_all_argument_types, i) {
    if (mem_all_argument_types(i)->name == name) {
      if (mem_all_argument_types(i)->type != ArgumentType::simple) {
        HALT("already exists with different type");
      }
      else
        return mem_all_argument_types(i);
    }
  }
  ArgumentType* t = new ArgumentType;
  t->name = name;
  mem_all_argument_types.append(t);
  return t;
}


ArgumentType* ArgumentType::read(ifstream& in) {
  uint DEBUG = 0;
  MT::String line;
  line.read(in, NULL, "\n");
  if (DEBUG>0) PRINT(line);
  if(line.N == 0) return NULL;
 
  uint type, typeI;
  MT::String name;
  line >> type;
  name.read(line, NULL, " ");
  line >> typeI;
  
  if (typeI == ArgumentType::simple) {
    ArgumentType* t = new ArgumentType;
    t->name = name;
    if (DEBUG>0) {t->write(); cout<<endl;}
    return t;
  }
  else if (typeI == ArgumentType::disjunction) {
    NIY;
    // TODO read list of strings and call "get"
#if 0
    ArgumentTypeL base_types;
    uint i, k;
    FOR1D(base_types__ids, i) {
      base_types.append((k));
      FOR1D(existing_types, k) {
        if (existing_types(k)->type == base_types__ids(i)) {
          
          break;
        }
      }
      CHECK(i < base_types__ids.N, "base type has not been found");
    }
    DisjunctionArgumentType* t = DisjunctionArgumentType::get(name, base_types);
    if (DEBUG>0) {cout<<*t<<endl;}
    return t;
#endif
  }
  else
    NIY;
  return NULL;
}


void ArgumentType::get(ArgumentTypeL& all_argumentTypes) {
  all_argumentTypes = mem_all_argument_types;
}


DisjunctionArgumentType* DisjunctionArgumentType::get(const MT::String& name, const ArgumentTypeL& base_types) {
  uint i;
  FOR1D(mem_all_argument_types, i) {
    if (mem_all_argument_types(i)->name == name) {
      if (mem_all_argument_types(i)->type != ArgumentType::disjunction) {
        HALT("already exists with different type");
      }
      else
        return (DisjunctionArgumentType*)mem_all_argument_types(i);
    }
  }
  DisjunctionArgumentType* t = new DisjunctionArgumentType();
  t->name = name;
  t->base_types = base_types;
  mem_all_argument_types.append(t);
  return t;
}

/************************************************
 * 
 *     Reading  &  Writing  of  Lists and Language
 * 
 ************************************************/


void writeSymbols(const SymL& symbols, ostream& os) {
  uint k;
  FOR1D(symbols, k) {
    if (symbols(k)!=NULL)
      symbols(k)->write(os);
    else
      os<<"NULL";
    os << " "<<endl;
  }
}


void writeSymbolsAndTypes(const SymL& syms, const ArgumentTypeL& types, const char* filename) {
  ofstream out;
  out.open(filename);
  writeSymbolsAndTypes(syms, types, out);
  out.close();
}


void writeSymbolsAndTypes(const SymL& syms, const ArgumentTypeL& types, ostream& out) {
  uint i;
  FOR1D(syms, i) {
    if (syms(i)->name == "default") continue;
    out << *syms(i) << endl;
  }
  out<<endl;  
  if (types.N > 0  &&  !(types.N==1 && types(0)->name == "any")) {
    out<<endl;
    out<<"["<<endl;
    FOR1D(types, i) {
      if (types(i)->name == "any") continue;
      out << *types(i) << endl;
    }
    out<<"]"<<endl;
  }
}


void writeSymbolsAndTypes(const char* filename) {
  SymL all_symbols;
  Symbol::get(all_symbols);
  ArgumentTypeL all_argumentTypes;
  ArgumentType::get(all_argumentTypes);
  writeSymbolsAndTypes(all_symbols, all_argumentTypes, filename);
}


void readSymbolsAndTypes(SymL& symbols, ArgumentTypeL& types, const char *filename) {
  ifstream in;
  in.open(filename);
  CHECK(in.is_open(), "Symbols file '"<<filename<<"' can't be opened!");
  readSymbolsAndTypes(symbols, types, in);
}


void readSymbolsAndTypes(SymL& symbols, ArgumentTypeL& types, istream& in) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"readSymbolsAndTypes [START]"<<endl;}
//   MT_MSG("ArgumentTypeL types IGNORED THUS FAR!!");
  symbols.clear();
  
  // READING
  // Ignore description beginning of file.
  // Symbols
  while (MT::skip(in) != -1  &&  MT::peerNextChar(in) != '[') {
    if (DEBUG>0) {cout<<"***** Reading next symbol..."<<endl;}
    MT::skip(in);
    MT::String line;
    line.read(in, NULL, "\n");
    if (DEBUG>1) {PRINT(line);}
    MT::String name;
    uint arity;
    name.read(line, NULL, " ");
    if (isupper(name(0))) HALT("Symbols need to start with a small letter. (Read symbol with name '"<<name<<"'.)");
    line >> arity;
    if (DEBUG>1) {PRINT(name);}
    if (DEBUG>1) {PRINT(arity);}
    
    // read type
    MT::String symbol_type__string;
    symbol_type__string.read(line, NULL, " ");
    if (DEBUG>1) {PRINT(symbol_type__string);}
    Symbol::SymbolType symbol_type;
    if (symbol_type__string == "action") symbol_type = Symbol::action;
    else if (symbol_type__string == "primitive") symbol_type = Symbol::primitive;
    else if (symbol_type__string == "conjunction") symbol_type = Symbol::conjunction;
    else if (symbol_type__string == "transclosure") symbol_type = Symbol::transclosure;
    else if (symbol_type__string == "count") symbol_type = Symbol::count;
    else if (symbol_type__string == "avg") symbol_type = Symbol::avg;
    else if (symbol_type__string == "max") symbol_type = Symbol::max;
    else if (symbol_type__string == "sum") symbol_type = Symbol::sum;
    else if (symbol_type__string == "function_change") symbol_type = Symbol::function_change;
    else if (symbol_type__string == "function_reward") symbol_type = Symbol::function_reward;
    else if (symbol_type__string == "function_difference") symbol_type = Symbol::function_difference;
    else HALT("unknown symbol_type:  "<<symbol_type__string);
    
    // read range
    MT::String range_type__string;
    range_type__string.read(line, NULL, " ");
    if (DEBUG>1) {PRINT(range_type__string);}
    Symbol::RangeType range_type;
    if (range_type__string == "binary") range_type = Symbol::binary;
    else if (range_type__string == "integer_set") range_type = Symbol::integer_set;
    else if (range_type__string == "integers") range_type = Symbol::integers;
    else if (range_type__string == "reals") range_type = Symbol::reals;
    else NIY;
    
    uintA range;
    if (range_type == Symbol::integer_set) {
      line >> range;
      if (DEBUG>1) {PRINT(range);}
      if (range.N == 0) {
          // default initialization of integer_set
          range.clear();
          uint i;
          for (i=0; i<=5; i++) range.append(i);
        }
      }
    
    Symbol* s;
    if (symbol_type == Symbol::primitive) {
      s = Symbol::get(name, arity, symbol_type, range_type);
    }
    else if (symbol_type == Symbol::action) {
      s = Symbol::get(name, arity, symbol_type, range_type);
    }
    else if (symbol_type == Symbol::conjunction) {
      MT::skip(line);
      MT::skipUntil(line, " ");  // skip "<--"
      while (MT::skip(line,"\n\r\t ,") != -1) {
        bool free_vars_all_quantified = false;
        while (isupper(MT::peerNextChar(line))) {
          if (MT::peerNextChar(line) == 'A') free_vars_all_quantified = true;
          else if (MT::peerNextChar(line) == 'E') free_vars_all_quantified = false;
          // else do nothing
          // --> leading variable can be ignored (only for visualization)
          MT::skipUntil(line, " ");
        }
        MT::String substring;
        line >> substring;
        LitL base_literals;
        if (DEBUG>0) {PRINT(substring);}
        Literal::get(base_literals, substring);
        if (DEBUG>0) {PRINT(base_literals);}
        s = ConjunctionSymbol::get(name, arity, base_literals, free_vars_all_quantified);
      }
    }
    else if (symbol_type == Symbol::transclosure) {
      MT::skip(line);
      MT::skipUntil(line, " ");  // skip "<--"
      if (MT::peerNextChar(line) != '+') HALT("incorrect specification of transclosure symbol (expecting '+')");
      MT::skipUntil(line, " ");  // skip "+"
      MT::String base_symbol_name;
      base_symbol_name.read(line, NULL, " ");
      if (DEBUG>2) {cout<<"base_symbol_name=\""<<base_symbol_name<<"\""<<endl;}
      Symbol* base_symbol = Symbol::get(base_symbol_name);
      if (base_symbol == NULL) {HALT("base symbol \""<<base_symbol_name<<"\" does not exist");}
      s = TransClosureSymbol::get(name, arity, base_symbol);
    }
    else if (symbol_type == Symbol::count) {
      MT::skip(line);
      MT::skipUntil(line, " ");  // skip "<--"
      if (MT::peerNextChar(line) != 'N') HALT("incorrect specification of transclosure symbol (expecting 'Num')");
      MT::skipUntil(line, " ");  // skip "num"
      // leading variables can be ignored (only for visualization)
      while (isupper(MT::peerNextChar(line))) {MT::skipUntil(line, " ");}
      MT::String remaining_string;
      line >> remaining_string;
      Literal* base_literal = Literal::get(remaining_string);
      if (base_literal == NULL) {HALT("base_literal "<<remaining_string<<" does not exist");}
      s = CountSymbol::get(name, arity, base_literal);
    }
    else if (symbol_type == Symbol::sum) {
      MT::skip(line);
      MT::skipUntil(line, " ");  // skip "<--"
      if (MT::peerNextChar(line) != 'S') HALT("incorrect specification of transclosure symbol (expecting 'Num')");
      MT::skipUntil(line, " ");  // skip "sum"
      // leading variables can be ignored (only for visualization)
      while (isupper(MT::peerNextChar(line))) {MT::skipUntil(line, " ");}
      MT::String remaining_string;
      line >> remaining_string;
      Symbol* base_symbol = Symbol::get(remaining_string);
      if (base_symbol == NULL) {HALT("base symbol "<<remaining_string<<" does not exist");}
      s = SumFunction::get(name, arity, base_symbol);
    }
    else if (symbol_type == Symbol::function_difference) {
      MT::skip(line);
      MT::skipUntil(line, " ");  // skip "<--"
      MT::String substring;
      line >> substring;
      LitL base_literals;
      base_literals.memMove = true;
      if (DEBUG>0) {PRINT(substring);}
      Literal::get(base_literals, substring);
      CHECK(base_literals.N >= 2, "Too few parameters!")
      CHECK(base_literals(0)->s == base_literals(1)->s, "Functions must be equal!")
      Symbol *baseSymbol = base_literals(0)->s;
      Literal *baseFLit1 = base_literals(0);
      Literal *baseFLit2 = base_literals(1);
      base_literals.remove(0);
      base_literals.remove(0);
      s = DifferenceFunction::get(name, arity, baseSymbol, baseFLit1, baseFLit2, base_literals);    //Really terrible code, must be improved!
    }
    else
      HALT("unknown symbol type "<<symbol_type);
    s->name = name;
    s->arity = arity;
    s->symbol_type = symbol_type;
    s->range_type = range_type;
    s->range = range;
    symbols.append(s);
    if (DEBUG>0) {cout<<"Read-in symbol:  "<<*symbols.last()<<endl;}
    MT::skip(in);
  }

  if (MT::peerNextChar(in) == '[') HALT("Reading types not implemented yet");
  // SORTING
  Symbol::sort(symbols);
  if (DEBUG>0) {cout<<"readSymbolsAndTypesSimpleFormat [END]"<<endl;}
}

  

}  // namespace PRADA



std::ostream& operator<<(std::ostream& os, const relational::Symbol& s) {
  s.write(os); return os;
}


std::ostream& operator<<(std::ostream& os, const relational::ArgumentType& t) {
  t.write(os); return os;
}


std::ostream& operator<<(std::ostream& os, const relational::SymL& symbols) {
  writeSymbols(symbols, os);
  return os;
}

