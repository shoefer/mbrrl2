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
#include <map>
#include "reason.h"
#include "plan.h"
#include "literals.h"

#include <set>
#include <algorithm>
#include <vector>

#define LE_fast 1



namespace relational {
  

/****************************************
  ABSTRACTING AND GROUNDING
 ***************************************/
  
 
uintA mem__constants;
ArgumentTypeL mem__constants_types;


bool reason::isConstant(uint obj) {
  if (getConstants().N > 0) {
    return getConstants().findValue(obj) >= 0;
  }
  else
    return obj >= 10;  // assumpation if no further information is given
}


bool reason::isGround(const Literal* lit) {
  uint i;
  FOR1D(lit->args, i) {
    if (!isConstant(lit->args(i)))
      return false;
  }
  return true;
}


bool reason::isGround(const SymbolicState& state) {
  uint i;
  FOR1D(state.lits, i) {
    if (!isGround(state.lits(i)))
      return false;
  }
  return true;
}


bool reason::isGround(const StateTransition& trans) {
  return isGround(trans.pre)  &&  isGround(trans.action)  &&  isGround(trans.post);
}


bool reason::isPurelyAbstract(const Literal* lit) {
  uint i;
  FOR1D(lit->args, i) {
    if (isConstant(lit->args(i)))
      return false;
  }
  return true;
}


void reason::setConstants(uintA& _constants) {
  ArgumentType* default_type = ArgumentType::get(MT::String("any"));
  ArgumentTypeL _constants_types;
  uint i;
  FOR1D(_constants, i) {_constants_types.append(default_type);}
  setConstants(_constants, _constants_types);
}


void reason::setConstants(uintA& _constants, const ArgumentTypeL& _constants_types) {
  uintA sortedIndices;
  TL::sort_asc(mem__constants, sortedIndices, _constants);
  if (_constants.N == _constants_types.N) {
    mem__constants_types.clear();
    uint i;
    FOR1D(sortedIndices, i) {
      mem__constants_types.append(_constants_types(sortedIndices(i)));
    }
  }
}


uintA& reason::getConstants() {
  return mem__constants;
}


ArgumentType* reason::getArgumentTypeOfObject(uint object_id) {
  int idx = mem__constants.findValue(object_id);
  if (idx >= 0) {
    return mem__constants_types(idx);
  }
  else {
    return ArgumentType::get(MT::String("any"));
  }
}





/****************************************
  CALCULATING DERIVED SYMBOLS
  ***************************************/
  

bool reason::derive_conjunction(LitL& lits_derived, ConjunctionSymbol& s, const LitL& lits_given, const uintA& constants) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_conjunction [START]"<<endl;}
//     double t_start = MT::cpuTime();
  if (DEBUG>1) {cout<<"ConjunctionSymbol: "<<s<<endl<<"lits_given: "<<lits_given<<endl;}
  uint i, j, k;
  lits_derived.clear();

  //filter lits_given
  LitL lits_filtered;
  FOR1D(lits_given, i) {
    FOR1D(s.base_literals, j) {
      // BUGFIX we must include base literals with opposing value too
      // TODO what with non-binary symbols?
      //if (lits_given(i)->s == s.base_literals(j)->s && TL::areEqual(lits_given(i)->value, s.base_literals(j)->value)) {
      if (lits_given(i)->s == s.base_literals(j)->s) {
        lits_filtered.append(lits_given(i));
        break;
      }
    }
  }

  // Existential
  if (!s.free_vars_all_quantified) {
    if (DEBUG>1) {cout<<"existential"<<endl;}
    // simply call cover
    SubstitutionSet subs;
    if (calcSubstitutions(subs, lits_filtered, s.base_literals, false)) {
      FOR1D_(subs, j) {
        if (s.arity == 2  &&  !subs.elem(j)->mapsToDistinct()) continue;
        uintA args;
        for (k=0; k<s.arity; k++) {args.append(subs.elem(j)->getSubs(k));}
        lits_derived.append(Literal::get(&s, args, 1.));
      }
    }
  }
  // all
  else {
    if (DEBUG>1) {cout<<"all"<<endl;}
    // We must treat positive vars with special care, since in cover(.) positives are always ex quantified.
    // calc free vars in base predicates
    uintA freeVars_pos, freeVars_neg;
    s.getFreeVars(freeVars_pos, freeVars_neg);
    if (DEBUG>1) {PRINT(freeVars_pos);  PRINT(freeVars_neg);}
    // only use constants which are provided in the state
    MT::Array< uintA > args_list;
    TL::allPermutations(args_list, constants, s.arity, false, true);
    if (DEBUG>0) {PRINT(args_list);}
    // investigate each possible arguments-tuple
    FOR1D(args_list, i) {
      if (DEBUG>2) {cout<<"++ Investigating "<<args_list(i)<<endl;}
      Substitution sub_args;
      FOR1D(args_list(i), j) {
        sub_args.addSubs(j, args_list(i)(j));
      }
      if (DEBUG>2) {PRINT(sub_args);}
      // For the remaining positive free variables, *all* instantiations must hold.
      MT::Array< uintA > args_freePos_lists;
      TL::allPermutations(args_freePos_lists, constants, freeVars_pos.d0, true, true);
      if (DEBUG>2) {PRINT(args_freePos_lists);}
      FOR1D(args_freePos_lists, j) {
        if (DEBUG>2) {cout<<"Checking "<<args_freePos_lists(j)<<endl;}
        Substitution sub_args__incl_freePos;
        sub_args__incl_freePos = sub_args;
        FOR1D(freeVars_pos, k) {
          sub_args__incl_freePos.addSubs(freeVars_pos(k), args_freePos_lists(j)(k));
        }
        if (DEBUG>2) {cout<<"sub_args__incl_freePos: "<<endl;  sub_args__incl_freePos.write();}
        if (!covers(lits_filtered, s.base_literals, true, &sub_args__incl_freePos)) {
          if (DEBUG>2) {cout<<"Fails."<<endl;}
          break;
        }
      }
      // if *all* instantiations held...
      if (args_freePos_lists.d0 == j) {
        uintA args;
        for (k=0; k<s.arity; k++) {args.append(sub_args.getSubs(k));}
        lits_derived.append(Literal::get(&s, args, 1.));
        if (DEBUG>2) {cout<<"Succeeds. New literal "<<*lits_derived.last()<<endl;}
      }
    }
  }
//     double t_finish = MT::cpuTime();
//     cout<<"derive_conjunction time = "<<(t_finish - t_start)<<endl;
  if (DEBUG>0) {cout<<"Derived literals ["<<lits_derived.N<<"]: "<<lits_derived<<endl;}
  if (DEBUG>0) {cout<<"derive_conjunction [END]"<<endl;}
  return lits_derived.N > 0;
}


// assumes acyclicity!
// seems standard graph problem i don't know standard solution of
bool reason::derive_transclosure(LitL& lits_derived, TransClosureSymbol& s, const LitL& lits_given) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_transclosure [START]"<<endl;}
  if (DEBUG>0) {cout<<s<<endl;}
  CHECK(s.arity==2, "transitive closure defined only for binary symbols with two arguments")
  lits_derived.clear();
  uint i;
  // (1) find edges
  std::map< uint, uintA > right;
  std::map< uint, uintA >::iterator iter;
  if (DEBUG>0) {cout<<"Given: ";}
  FOR1D(lits_given, i) {
    if (lits_given(i)->s == s.base_symbol) {
      if (DEBUG>0) {lits_given(i)->write(cout);cout<<" ";}
      right[lits_given(i)->args(0)].setAppend(lits_given(i)->args(1));
    }
  }
  if (DEBUG>0) {cout<<endl;}
  if (DEBUG>1) {
      cout<<"Direct right neighbors:"<<endl;
      for (iter = right.begin(); iter != right.end(); iter++) {
          cout<<iter->first<<": "<<iter->second<<endl;
      }
  }
  if (right.empty()) return false;
  // (2) build connections
  bool extended;
  uint num;
  do {
    extended = false;
    for (iter = right.begin(); iter != right.end(); iter++) {
      num = iter->second.d0;
      uintA newGuys_candidates;
      FOR1D(iter->second, i) {
        newGuys_candidates.setAppend(right[iter->second(i)]);
      }
      iter->second.setAppend(newGuys_candidates);
      if (iter->second.d0 > num)
        extended = true;
    }
  } while (extended);
  if (DEBUG>1) {
    cout<<"All right neighbors:"<<endl;
    for (iter = right.begin(); iter != right.end(); iter++) {
        cout<<iter->first<<": "<<iter->second<<endl;
    }
  }
  // (3) build TCP tuples
  if (DEBUG>0) {cout<<"logicObjectManager::p_derived: ";}
  for (iter = right.begin(); iter != right.end(); iter++) {
    FOR1D(iter->second, i) {
      uintA args(2);
      args(0)=iter->first;
      args(1)=iter->second(i);
      lits_derived.setAppend(Literal::get(&s, args, 1.));
      if (DEBUG>0) {cout<<*lits_derived.last()<<endl;}
    }
  }
  if (DEBUG>0) {cout<<endl;}
  if (DEBUG>0) {cout<<"derive_transclosure [END]"<<endl;}
  return true;
}



bool reason::derive_count(LitL& lits_derived, CountSymbol& s, const LitL& lits_given, const uintA& constants) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_count [START]"<<endl;}
  if (DEBUG>0) {cout<<s<<endl;  PRINT(lits_given);  PRINT(constants);}
  uint i, j, k;
  MT::Array< uintA > args_lists;
  TL::allPermutations(args_lists, constants, s.arity, true, true);

  uintA freeVarsPos(s.arity);
  FOR1D(s.base_literal->args, i) {
    if (s.base_literal->args(i) < s.arity)
      freeVarsPos(s.base_literal->args(i)) = i;
  }
  //filter lits_given
  LitL lits_filtered;
  FOR1D(lits_given, i) {
    if (lits_given(i)->s == s.base_literal->s && lits_given(i)->value == s.base_literal->value)
      lits_filtered.append(lits_given(i));
  }

  if (DEBUG>0) {cout<<s<<endl;  PRINT(lits_given); PRINT(lits_filtered); PRINT(constants); PRINT(freeVarsPos); }

  FOR1D(args_lists, i) {
    if (DEBUG>1) {PRINT(args_lists(i))}
    uint counter = 0;
    
    FOR1D(lits_filtered, j) {
      FOR1D(args_lists(i), k) {
        if (lits_filtered(j)->args(freeVarsPos(k)) == args_lists(i)(k))
          counter++;
      }
    }    

    /*Substitution sub_countsymbol;
    FOR1D(args_lists(i), k) {
      sub_countsymbol.addSubs(k, args_lists(i)(k));
    }
    if (DEBUG>1) {PRINT(sub_countsymbol)}
    uintA free_vars__base_literal;
    FOR1D(s.base_literal->args, k) {
      if (s.base_literal->args(k) >= s.arity) free_vars__base_literal.append(s.base_literal->args(k));
    }
    if (DEBUG>1) {PRINT(free_vars__base_literal);}
    SubstitutionSet subs_free_vars;
    SubstitutionSet::createAllPossibleSubstitutions(subs_free_vars, free_vars__base_literal, constants);
    FOR1D_(subs_free_vars, k) {
      Substitution* sub_total = Substitution::combine(sub_countsymbol, *subs_free_vars.elem(k));
      if (DEBUG>1) {PRINT(*sub_total);}
      Literal* lit_test = sub_total->apply(s.base_literal);
      if (holds(lits_filtered, lit_test)) counter++;
      delete sub_total;
    }*/
    lits_derived.append(Literal::get(&s, args_lists(i), counter));
    if (DEBUG>0) {cout<<" "<<*lits_derived.last()<<endl;}
  }
  if (DEBUG>0) {cout<<"derive_count [END]"<<endl;}
  return true; // always holds
}


bool reason::derive_functiondiff(LitL& lits_derived, DifferenceFunction& s, const LitL& lits_given, const uintA& constants) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_functiondiff [START]"<<endl;}

  SubstitutionSet subs;
  calcSubstitutions(subs, lits_given, s.restrictionLits, false);

  uint i,j;
  LitL baseSymbolLits;
  FOR1D(lits_given, i) {
    if (lits_given(i)->s == s.baseSymbol)
      baseSymbolLits.append(lits_given(i));
  }
  
  FOR1D_(subs, i) {
    double val1 = 0;
    FOR1D(baseSymbolLits, j) {
      if (baseSymbolLits(j)->args(0) == subs.elem(i)->getSubs(s.baseFunctionLit1->args(0)))
        val1 = baseSymbolLits(j)->value;
    }
    double val2 = 0;
    FOR1D(baseSymbolLits, j) {
      if (baseSymbolLits(j)->args(0) == subs.elem(i)->getSubs(s.baseFunctionLit2->args(0)))
        val2 = baseSymbolLits(j)->value;
    }
    double diff = val2-val1;
    //if (diff < 0) diff = 0;
    uintA args;
    args.append(subs.elem(i)->getSubs(0));
    lits_derived.append(Literal::get(&s, args, diff));
  }
  return true;
}


bool reason::derive_avg(LitL& lits_derived, AverageFunction& s, const LitL& lits_given, const uintA& constants) {
  NIY;
#if 0
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_avg [START]"<<endl;}
  CHECK(f.arity==0, "only implemented for zero-ary avg functions");
  uint i, k;
  MT::Array< uintA > combos;
  uintA local_constants;
  getConstants(s, local_constants);
  TL::allPermutations(combos, local_constants, f.f_base->arity, true, true);
  double avg = 0.0;
  FOR1D(combos, i) {
    if (DEBUG>1) {PRINT(combos(i))}
    if (f.f_base->category == category_primitive) {
      FOR1D(s.fv_prim, k) {
        if (s.fv_prim(k)->atom->f == f.f_base  
            &&  s.fv_prim(k)->args == combos(i)) {
          avg += s.fv_prim(k)->value;
          break;
        }
      }
    }
    else {
      FOR1D(s.fv_derived, k) {
        if (s.fv_derived(k)->atom->f == f.f_base  
            &&  s.fv_derived(k)->args == combos(i)) {
          avg += s.fv_derived(k)->value;
          break;
        }
      }
    }
  }
  avg /= combos.N;
  // function value bauen
  uintA empty;
  FunctionValue* fv = logicObjectManager::getFV(&f, empty, avg);
  s.fv_derived.setAppend(fv);
  if (DEBUG>0) {cout<<" ";fv->write(cout);}
  if (DEBUG>0) {cout<<"derive_avg [END]"<<endl;}
  return true; // always holds
#endif
}


bool reason::derive_sum(LitL& lits_derived, SumFunction& s, const LitL& lits_given, const uintA& constants) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_sum [START]"<<endl;}
  CHECK(s.arity==0, "only implemented for zero-ary avg functions");
  uint i, k;
  MT::Array< uintA > combos;
  TL::allPermutations(combos, constants, s.base_symbol->arity, true, true);
  double sum = 0.0;
  FOR1D(combos, i) {
    if (DEBUG>1) {PRINT(combos(i))}
    FOR1D(lits_given, k) {
      if (lits_given(k)->s == s.base_symbol
            &&  lits_given(k)->args == combos(i)) {
        sum += lits_given(k)->value;
        break;
      }
    }
  }
  uintA empty;
  Literal* lit_new = Literal::get(&s, empty, sum);
  if (DEBUG>0) {cout<<" "<<*lit_new<<endl;}
  lits_derived.setAppend(lit_new);
  if (DEBUG>0) {cout<<"derive_sum [END]"<<endl;}
  return true; // always holds
}


bool reason::derive_max(LitL& lits_derived, MaxFunction& s, const LitL& lits_given, const uintA& constants) {
  NIY;
#if 0
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_max [START]"<<endl;}
  CHECK(f.arity==0, "only implemented for zero-ary max functions");
  uint i, k;
  MT::Array< uintA > combos;
  uintA local_constants;
  getConstants(s, local_constants);
  TL::allPermutations(combos, local_constants, f.f_base->arity, true, true);
  arr values;
  FOR1D(combos, i) {
    if (DEBUG>1) {PRINT(combos(i))}
    if (f.f_base->category == category_primitive) {
      FOR1D(s.fv_prim, k) {
        if (s.fv_prim(k)->atom->f == f.f_base  
            &&  s.fv_prim(k)->args == combos(i)) {
          values.append(s.fv_prim(k)->value);
          break;
        }
      }
    }
    else {
      FOR1D(s.fv_derived, k) {
        if (s.fv_derived(k)->atom->f == f.f_base  
            &&  s.fv_derived(k)->args == combos(i)) {
          values.append(s.fv_derived(k)->value);
          break;
        }
      }
    }
  }
  double maxVal = values.max();
  // function value bauen
  uintA empty;
  FunctionValue* fv = logicObjectManager::getFV(&f, empty, maxVal);
  s.fv_derived.setAppend(fv);
  if (DEBUG>0) {cout<<" ";fv->write(cout);}
  if (DEBUG>0) {cout<<"derive_max [END]"<<endl;}
  return true; // always holds
#endif
}


bool reason::derive_reward(LitL& lits_derived, RewardFunction& s, const LitL& lits_given, const uintA& constants) {
  NIY;
#if 0
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"derive_reward [START]"<<endl;}
  CHECK(f.arity==0, "only zero-ary functions");
  double value = 0;
  if (reason::holds(s, f.ground_pis))
    value = f.reward_value;
  // function value bauen
  uintA empty;
  empty.resize(f.arity);
  FunctionValue* fv = logicObjectManager::getFV(&f, empty, value);
  s.fv_derived.setAppend(fv);
  if (DEBUG>0) {cout<<" ";fv->write(cout);}
  if (DEBUG>0) {cout<<"derive_reward [END]"<<endl;}
  return true; // always holds
#endif
}


// construct logicObjectManager::p_derived predicates and functions
// assumptions: no self recursion
// negation in base predicates is only allowed for primitive predicates
void reason::derive(LitL& lits_derived, const LitL& lits_given, const uintA& constants) {
  uint DEBUG = 0;
  if (DEBUG > 0) cout << "derive [START]" << endl;
  if (DEBUG > 1) {PRINT(lits_given);}  
  SymL symbols;
  Symbol::get(symbols);
  uint i;
  FOR1D(symbols, i) {
    LitL lits_all;  lits_all.append(lits_given);  lits_all.append(lits_derived);
    if (DEBUG>1) {cout<<"+++ Derive symbol: "<<*symbols(i)<<endl;}
    if (symbols(i)->symbol_type == Symbol::action) continue;
    else if (symbols(i)->symbol_type == Symbol::primitive) continue;
    else if (symbols(i)->symbol_type == Symbol::conjunction) {
      ConjunctionSymbol* s = dynamic_cast<ConjunctionSymbol*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_conjunction(lits_derived_local, *s, lits_all, constants);
      lits_derived.append(lits_derived_local);
    }
    else if (symbols(i)->symbol_type == Symbol::transclosure) {
      TransClosureSymbol* s = dynamic_cast<TransClosureSymbol*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_transclosure(lits_derived_local, *s, lits_all);
      lits_derived.append(lits_derived_local);
    }
    else if (symbols(i)->symbol_type == Symbol::count) {
      CountSymbol* s = dynamic_cast<CountSymbol*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_count(lits_derived_local, *s, lits_all, constants);
      lits_derived.append(lits_derived_local);
    }
    else if (symbols(i)->symbol_type == Symbol::avg) {
      AverageFunction* s = dynamic_cast<AverageFunction*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_avg(lits_derived_local, *s, lits_all, constants);
      lits_derived.append(lits_derived_local);
    }
    else if (symbols(i)->symbol_type == Symbol::max) {
      MaxFunction* s = dynamic_cast<MaxFunction*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_max(lits_derived_local, *s, lits_all, constants);
      lits_derived.append(lits_derived_local);
    }
    else if (symbols(i)->symbol_type == Symbol::sum) {
      SumFunction* s = dynamic_cast<SumFunction*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_sum(lits_derived_local, *s, lits_all, constants);
      lits_derived.append(lits_derived_local);
    }
    else if (symbols(i)->symbol_type == Symbol::function_reward) {
      RewardFunction* s = dynamic_cast<RewardFunction*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_reward(lits_derived_local, *s, lits_all, constants);
      lits_derived.append(lits_derived_local);
    }
    else if (symbols(i)->symbol_type == Symbol::function_difference) {
      DifferenceFunction *s = dynamic_cast<DifferenceFunction*>(symbols(i));
      CHECK(s!=NULL, "cast failed");
      LitL lits_derived_local;
      derive_functiondiff(lits_derived_local, *s, lits_all, constants);
      lits_derived.append(lits_derived_local);
    }
    else
      HALT("Unknown symbol type: "<<*symbols(i));    
  }
  Literal::sort(lits_derived);
  if (DEBUG > 0) {cout<<"Derived literals: "<<lits_derived<<endl;}
  if (DEBUG > 0) cout << "derive [END]" << endl;
}


void reason::derive(SymbolicState* state) {
  if (state->derived_lits_are_calculated)
    return;
  LitL lits_derived;
  if (state->state_constants.N == 0) {
    Literal::getArguments(state->state_constants, state->lits);
  }
  derive(state->lits, lits_derived, state->state_constants);
  state->lits.append(lits_derived);
  state->derived_lits_are_calculated = true;
}

void reason::dederive(SymbolicState* state) {
  uint i;
  FOR1D_DOWN(state->lits, i) {
    if (state->lits(i)->s->symbol_type != Symbol::primitive)
      state->lits.remove(i);
  }
  state->derived_lits_are_calculated = false;
}








/****************************************
    HOLDS (--> for ground)
***************************************/


bool reason::consistent(Literal* lit_given, Literal* lit_test) {
  uint DEBUG = 0;
  if (DEBUG>0) cout << "consistent [START]" << endl;
  if (DEBUG>0) {PRINT(*lit_given);  PRINT(*lit_test);}
  // special treatment for negated binary literals
  if (lit_test->isNegated()) {HALT("negated lit_test not allowed: "<<*lit_test);}
  bool does_hold = false;
  if (lit_given->s == lit_test->s
      &&  lit_given->args == lit_test->args) {
    if (lit_given->comparison_type == Literal::comparison_variable || lit_test->comparison_type == Literal::comparison_variable) {
      does_hold = true;     //Variable comparison always holds
    }
    else if (lit_given->comparison_type == Literal::comparison_equal) {
      does_hold = lit_given->compareValueTo(lit_test->comparison_type, lit_test->value);
    }
    else if (lit_test->comparison_type == Literal::comparison_equal) {
      does_hold = lit_test->compareValueTo(lit_given->comparison_type, lit_given->value);
    }
    else
      NIY;
  }
  if (DEBUG>0) cout<<" --> "<<does_hold<<endl;
  if (DEBUG>0) cout << "consistent [END]" << endl;
  return does_hold;
}


bool reason::holds(const LitL& lits_given, Literal* lit_test) {
  uint DEBUG = 0;
  if (DEBUG>0) cout << "holds [START]" << endl;
  if (DEBUG>0) {PRINT(lits_given);  PRINT(*lit_test);}
  bool does_hold = false;
  if (lit_test->isNegated()) {
    does_hold = !holds(lits_given, lit_test->getNegated());
  }
  else {
    uint i;
    FOR1D(lits_given, i) {
      if (consistent(lits_given(i), lit_test)) {
        does_hold = true;
        break;
      }
    }
  }
  if (DEBUG>0) cout<<" --> "<<does_hold<<endl;
  if (DEBUG>0) cout << "holds [END]" << endl;
  return does_hold;
}


bool reason::holds(const LitL& lits_given, const LitL& lits_test) {
  uint i;
  FOR1D(lits_test, i) {
    if (!holds(lits_given, lits_test(i))) return false;
  }
  return true;
}





/****************************************
    COVERAGE (--> for abstract)
  ***************************************/


bool reason::calcSubstitution(Substitution& sub, const uintA& ground_args, const uintA& other_args) {
  uint i;
  FOR1D(ground_args, i) {
    CHECK(isConstant(ground_args(i)), "ground_args contains not-ground "<<ground_args(i))
  }
  CHECK(ground_args.N==other_args.N, "wrong slot assignments:  ground_args="<<ground_args<<"  vs.  other_args="<<other_args);
  FOR1D(other_args, i) {
    if (isConstant(other_args(i))) {
      if (other_args(i) != ground_args(i))
        return false;
    }
    else {
      if (sub.hasSubs(other_args(i))) {
        if (sub.getSubs(other_args(i)) != ground_args(i))
          return false;
      }
      else {
        sub.addSubs(other_args(i), ground_args(i));
      }
    }
  }
  return true;
}


bool reason::calcSubstitution(Substitution& sub, const Literal* ground_lit, const Literal* other_lit) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcSubstitution [START]"<<endl;}
  if (DEBUG>0) {cout<<*ground_lit<<endl<<*other_lit<<endl;}
  CHECK(!ground_lit->isNegated(), "calcSubstitution only defined for non-negated literals, but ground_lit="<<*ground_lit);
  if (!TL::areEqual(ground_lit->value, other_lit->value)) {
    if (DEBUG>0) {cout<<"no coverage"<<endl<<"calcSubstitution [END]"<<endl;}
    return false;
  }
  if (!isGround(ground_lit)) {
    HALT("Literal is not ground:  "<<*ground_lit<<"  -- Have you set appropriate constants? constants="<<getConstants());
  }
  if (ground_lit->s != other_lit->s) {
    if (DEBUG>0) {cout<<"no coverage"<<endl<<"calcSubstitution [END]"<<endl;}
    return false;
  }
  if (ground_lit->args.N != other_lit->args.N) { 
      cout<<"ground_lit: "<<*ground_lit<<endl;
      cout<<"other_lit: "<<*other_lit<<endl;
      sub.write(cout);
      HALT("");
  }  
  bool covers = calcSubstitution(sub, ground_lit->args, other_lit->args);
  if (DEBUG>0) {PRINT(sub);  PRINT(covers);}
  if (DEBUG>0) {cout<<"calcSubstitution [END]"<<endl;}
  return covers;
}


bool reason::calcSubstitutions(SubstitutionSet& subs, const LitL& lits_ground, Literal* lit_input, bool free_neg_vars__all_quantified, Substitution* initSub) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"calcSubstitutions(..., const LitL& lits_ground, Literal* lit_input, ...) [START]"<<endl;
  CHECK(subs.num()==0, "There are already subs!")
  if (DEBUG>1) {
    cout << "lits_ground: "<<lits_ground<<endl;
    cout << "lit_input: "<<*lit_input<<endl;
    PRINT(free_neg_vars__all_quantified)
    cout << "initSub: "; if (initSub==NULL) cout<<" NULL"; else initSub->write(cout); cout<<endl;
  }
    
  // ------------------------------------------------------------
  // Grounded after initSub
  Literal* lit_input__ground_init;
  if (initSub != NULL)
    lit_input__ground_init = initSub->apply(lit_input);
  else
    lit_input__ground_init = lit_input;
  if (isGround(lit_input__ground_init)) {
    if (DEBUG>1) cout<<"lit is ground after initSub:  lit_input__ground_init="<<*lit_input__ground_init<<endl;
    if (lit_input__ground_init->isNegated()) {
      Literal* lit_input__ground_init__pos;
      lit_input__ground_init__pos = lit_input__ground_init->getNegated();
      if (!holds(lits_ground, lit_input__ground_init__pos)) {
        Substitution* s = new Substitution;
        if (initSub != NULL) *s = *initSub;
        subs.append(s);
      }
    }
    else {
      if (holds(lits_ground, lit_input__ground_init)) {
        Substitution* s = new Substitution;
        if (initSub != NULL) *s = *initSub;
        subs.append(s);
      }
    }
    if (DEBUG>1) cout<<"covered? " << (subs.num()>0) << endl;
    if (DEBUG>0) cout<<"cover(..., const LitL& lits_ground, Literal* lit_input, ...) [END]"<<endl;
    if (subs.num() > 0)
        return true;
    else
        return false;
  }
    
  // ------------------------------------------------------------
  // Non-ground after initSub
  if (DEBUG>0) {cout<<"non-ground after initSub"<<endl;}
  uint i;
  if (!lit_input__ground_init->isNegated()) {
    if (DEBUG>0) {cout<<"Positive"<<endl;}
    bool literalTrue;
    FOR1D(lits_ground, i) {
      Substitution sub1;
      literalTrue = calcSubstitution(sub1, lits_ground(i), lit_input__ground_init);
      if (literalTrue) {subs.append(Substitution::combine(sub1, *initSub));}
    }
  }
  else {
    if (DEBUG>0) {cout<<"Negative binary literal"<<endl;}
    // ------------------------------------------------------------
    // All-quantified
    // Trick: Algorithm checks that the positive version cannot be unified with state.
    if (free_neg_vars__all_quantified) {
      if (DEBUG>0) {cout<<"All quantification"<<endl;}
      bool literalTrue = false;
      uint i;
      Literal* lit_input__ground_init_pos = lit_input__ground_init->getNegated();
      FOR1D(lits_ground, i) {
        Substitution sub1;
        if (initSub != NULL) sub1 = *initSub;
        literalTrue = calcSubstitution(sub1, lits_ground(i), lit_input__ground_init_pos);
        if (literalTrue) break;
      }
      if (lits_ground.N == i) {  // lit_pos cannot be unified
        literalTrue = true;
        subs.append(initSub);
      }
    }
    // ------------------------------------------------------------
    // Existential-quantified
    // Trick: Algorithm tries all possible substitutions.
    else {
      if (DEBUG>0) {cout<<"Exists quantification"<<endl;}
      uintA freeVars;
      FOR1D(lit_input__ground_init->args, i) {
        if (!isConstant(lit_input__ground_init->args(i)))
          freeVars.setAppend(lit_input__ground_init->args(i));
      }
      // create all possible substitutions
      SubstitutionSet cand_subs;
      SubstitutionSet::createAllPossibleSubstitutions(cand_subs, freeVars, reason::getConstants());
      // check which lead to coverage by creating positive literals and unifying them with the given state
      uintA state_constants;
      Literal::getArguments(state_constants, lits_ground);
      FOR1D_(cand_subs, i) {
        // filter out constants which don't appear in the state
        uintA outs_distinct;  cand_subs.elem(i)->getOutsDistinct(outs_distinct);
        if (DEBUG>1) {cout<<"cand_subs.elem(i="<<i<<"):  ";  cand_subs.elem(i)->write();  cout<<endl;  PRINT(outs_distinct);}
        if (numberSharedElements(outs_distinct, state_constants) != outs_distinct.N) {
          if (DEBUG>1) {cout<<" -> invalid cand_subs"<<endl;}
          continue;
        }
        Literal* lit_input__ground_full = cand_subs.elem(i)->apply(lit_input__ground_init);
        Literal* lit_input__ground_full__pos = lit_input__ground_full->getNegated();
        // If positive version of ground literal does not hold,
        // then its negation _does_ hold.
        if (!holds(lits_ground, lit_input__ground_full__pos)) {
          subs.append(Substitution::combine(*initSub, *cand_subs.elem(i)));
          if (DEBUG>1) {cout<<" -> covers (positive literal does not hold!)"<<endl;}
        }
      }
    }
  }

  if (DEBUG>1) {
    if (subs.num()==0) cout << "Not covered." << endl;
    else {
      cout <<"Covered by the following "<<subs.num()<<" substitutions:"<<endl;
      subs.write();
    }
  }
  if (DEBUG>0) {cout<<"Covering? "<<(subs.num()>0)<<endl;}
  if (DEBUG>0) cout<<"calcSubstitutions(..., const LitL& lits_ground, Literal* lit_input, ...) [END]"<<endl;
  return subs.num() > 0;
}


bool reason::calcSubstitutions(SubstitutionSet& subs, const LitL& lits_ground, const LitL& lits_input, bool free_neg_vars__all_quantified, Substitution* initSub) {
  int DEBUG = 0;
  if (DEBUG>0) cout<<"calcSubstitutions (..., const LitL& lits_ground, const LitL& lits_input, ...) [START]"<<endl;
  if (DEBUG > 0) {
    cout<<"lits_input: "<<lits_input<<endl;
    cout<<"lits_ground: "<<lits_ground<<endl;
    cout << "initSub: ";
    if (initSub != NULL) {initSub->write(cout); cout<<endl;}
    else {cout<<"-"<<endl;}
  }
  CHECK(Literal::negativeBinaryLiteralsLast(lits_input), "Positive literals need to be first!");
  CHECK(subs.num()==0, "Subs has to be empty.");
  Substitution* initSub_copy = new Substitution;
  if (initSub != NULL)
    *initSub_copy = *initSub;
  subs.append(initSub_copy);
  uint i, l, k;
  FOR1D(lits_input, i) {
    SubstitutionSet all_next_subs;
    if (DEBUG > 0) {cout << "Inspecting lit #"<<i<<": "<<*lits_input(i)<<endl;}
    FOR1D_(subs, l) {
      if (DEBUG>2) {cout<<"Providing init sub: "<<*subs.elem(l)<<endl;}
      SubstitutionSet next_subs;
      if (calcSubstitutions(next_subs, lits_ground, lits_input(i), free_neg_vars__all_quantified, subs.elem(l))) {
        FOR1D_(next_subs, k) {
          all_next_subs.append(next_subs.elem(k));
        }
      }
    }
    subs = all_next_subs;
    // if no more subs found = if no more coverage
    if (all_next_subs.num() == 0) {
      break;
    }
    if (DEBUG > 2) {
      cout << "After lit #"<<i<<": "; PRINT(subs.num())
      subs.write();
    }
  }
  if (DEBUG>0) {cout << "Covered: " << (subs.num() > 0) << endl;}
  if (DEBUG>0) cout<<"calcSubstitutions (..., const LitL& lits_ground, const LitL& lits_input, ...) [END]"<<endl;
  return subs.num() > 0;
}


bool reason::covers(const LitL& lits_ground, const LitL& lits_input, bool free_neg_vars__all_quantified, Substitution* initSub) {
  SubstitutionSet subs_unneeded;
  return calcSubstitutions(subs_unneeded, lits_ground, lits_input, free_neg_vars__all_quantified, initSub);
}




/****************************************
        STATE UNIFICATION
  ***************************************/

bool reason::calcSubstitutions(SubstitutionSet& subs, const SymbolicState& state1, const SymbolicState& state2, Substitution* initSub) {
  uint minDiff = calcSubstitutions_asMuchAsPossible(subs, state1, state2, initSub);
  return minDiff == 0;
}

/**
 *
 */
bool reason::unifiableIgnoreCommonConstants(const SymbolicState& s1, const SymbolicState& s2) {
  SubstitutionSet subs;
  return unifiableIgnoreCommonConstants(subs, s1, s2);
}

bool reason::unifiableIgnoreCommonConstants(SubstitutionSet& subs, const SymbolicState& s1, const SymbolicState& s2_) {
  SymbolicState s2 = s2_;

  CHECK (s1.state_constants.N > 0, "Did you forget to compute the state constants for s1?");
  CHECK (s2.state_constants.N > 0, "Did you forget to compute the state constants for s2?");

  uintA common_constants;
  uint i;
  FOR1D(s1.state_constants, i){
      uint c = s1.state_constants(i);
      if (s2.state_constants.contains(c))
          common_constants.append(c);
  }
  if (common_constants.N > 0) {
      uintA constants;
      constants.setAppend(s1.state_constants);
      constants.setAppend(s2.state_constants);
      // arbitrarily resubstitute
      uint start = (uint) constants.absMax() + 1;
      Substitution sub;
      FOR1D(common_constants, i) {
          sub.addSubs(common_constants(i), start++);
      }
      SymbolicState* sub_s2 = sub.apply(s2);
      s2 = *sub_s2;
      Literal::getArguments(s2.state_constants, s2.lits);
      delete sub_s2;
  }


  using namespace std;

  // brute-force: generate all possible substitutions and check which of them works
  set<uint> vars1_set, vars2_set;

  // collect args
  FOR1D(s1.lits, i) {
      uint j;
      FOR1D(s1.lits(i)->args, j)
              vars1_set.insert(s1.lits(i)->args(j));
  }
  FOR1D(s2.lits, i) {
      uint j;
      FOR1D(s2.lits(i)->args, j)
              vars2_set.insert(s2.lits(i)->args(j));
  }

  // must have same number of variables; otherwise we have to
  // map on same variable which is not supported
  if (vars1_set.size() != vars2_set.size()) {
      return false;
  }

  vector<uint> vars1(vars1_set.begin(), vars1_set.end());
  vector<uint> vars2(vars2_set.begin(), vars2_set.end());

  sort(vars1.begin(), vars1.end());

  do {
      do {
          vector<uint>::iterator it1 = vars1.begin();
          vector<uint>::iterator it2 = vars2.begin();
          Substitution* s = new Substitution;
          for(;it1!=vars1.end() && it2!=vars2.end(); it1++, it2++) {
              s->addSubs(*it1, *it2);
          }
          //s->write(cout);
          SymbolicState* s1_sub = s->apply(s1);
          if (*s1_sub == s2) {
              delete s1_sub;
              subs.append(s);
              return true;
          }
          delete s1_sub;
          delete s;
      } while(std::next_permutation(vars2.begin(), vars2.end()));

  } while(std::next_permutation(vars1.begin(), vars1.end()));


  return false;
}



bool reason::unifiable(const SymbolicState& state1, const SymbolicState& state2) {
  SubstitutionSet subs;
  return calcSubstitutions(subs, state1, state2);
}


uint reason::calcSubstitutions_asMuchAsPossible(SubstitutionSet& subs, const SymbolicState& state1, const SymbolicState& state2, Substitution* initSub) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcSubstitutions_asMuchAsPossible [START]"<<endl;}
  
  if (DEBUG>0) {
    cout << "SymbolicState 1:   "; state1.write(cout, true);  cout<<endl;
    cout << "SymbolicState 2:   "; state2.write(cout, true);  cout<<endl;
    cout<<"initSub:  "<<flush; if (initSub == NULL) cout<<"NULL"<<endl;  else initSub->write();  cout<<endl;
  }
  
  // (1) Calc Differences (with initSub)
//   SymbolicState* state1_initSubbed = initSub->apply(state1);
//   makeOriginal(*state1_initSubbed);
  LitL lits_diff_1to2, lits_diff_2to1;
  uintA changedConstants;
  SymbolicState::calcDifferences(lits_diff_1to2, lits_diff_2to1, changedConstants, state1, state2);
  if (DEBUG>0) {
//     cout << "SymbolicState 1 (with init-sub):  "; state1_initSubbed->writeNice(cout, false, true);  cout<<endl;
    cout<<"Differences:"<<endl;
    cout<<"lits_diff_1to2:    "; write(lits_diff_1to2); cout<<endl;
    cout<<"lits_diff_2to1:    "; write(lits_diff_2to1); cout<<endl;
  }
//   delete state1_initSubbed;
  
  
  uint i, k, l;
  
  
  // (2) Determine substitutions leading from state1 to state2
  uintA in_ids, out_ids;
  SubstitutionSet potential_subs;
  if (initSub == NULL)
    potential_subs.append(new Substitution);
  else {
    potential_subs.append(initSub);
    initSub->getIns(in_ids);
    initSub->getIns(out_ids);
  }
  
  if (DEBUG>0) {PRINT(in_ids);}
  
  // (2.1) LOOK AT UNARIES
  // LitL
  FOR1D(lits_diff_1to2, i) {
    if (lits_diff_1to2(i)->s->arity == 1   &&  in_ids.findValue(lits_diff_1to2(i)->args(0)) < 0) {
      uintA possibleSubstitutes;
      FOR1D(lits_diff_2to1, k) {
        if (lits_diff_1to2(i)->s == lits_diff_2to1(k)->s    // same predicate
//                &&    substituted_ids.findValue(lits_diff_2to1(k)->args(0)) < 0      // substituting-id not already used
               ){
          possibleSubstitutes.append(lits_diff_2to1(k)->args(0));
        }
      }
      FOR1D(possibleSubstitutes, k) {
        FOR1D_(potential_subs, l) {
          potential_subs.elem(l)->addSubs(lits_diff_1to2(i)->args(0), possibleSubstitutes(k));
          in_ids.setAppend(lits_diff_1to2(i)->args(0));
        }
      }
    }
  }
  
  // (2.2) LOOK AT BINARIES
  FOR1D(lits_diff_1to2, i) {
    if (lits_diff_1to2(i)->s->arity == 2) {
      bool first_arg_covered = in_ids.findValue(lits_diff_1to2(i)->args(0)) >= 0;
      bool second_arg_covered = in_ids.findValue(lits_diff_1to2(i)->args(1)) >= 0;
      if (DEBUG>2) {lits_diff_1to2(i)->write();  cout<<"  "<<first_arg_covered<<"  "<<second_arg_covered<<endl;}
      if (!first_arg_covered) {
        uintA possibleSubstitutes;
        FOR1D(lits_diff_2to1, k) {
          if (lits_diff_1to2(i)->s == lits_diff_2to1(k)->s     // same predicate
//                 &&   substituted_ids.findValue(lits_diff_2to1(k)->args(0)) < 0         // substituting-id not already used
                &&   lits_diff_1to2(i)->args(0) != lits_diff_2to1(k)->args(0)) {  // different arguments --> true substitution
            possibleSubstitutes.append(lits_diff_2to1(k)->args(0));
          }
        }
        if (DEBUG>2) {cout<<"1:  "<<possibleSubstitutes<<endl;}
        FOR1D(possibleSubstitutes, k) {
          FOR1D_(potential_subs, l) {
            potential_subs.elem(l)->addSubs(lits_diff_1to2(i)->args(0), possibleSubstitutes(k));
            in_ids.setAppend(lits_diff_1to2(i)->args(0));
          }
        }
      }
      if (!second_arg_covered) {
        uintA possibleSubstitutes;
        FOR1D(lits_diff_2to1, k) {
          if (lits_diff_1to2(i)->s == lits_diff_2to1(k)->s     // same predicate
//                 &&   substituted_ids.findValue(lits_diff_2to1(k)->args(1)) < 0         // substituting-id not already used
                &&   lits_diff_1to2(i)->args(1) != lits_diff_2to1(k)->args(1)) {  // different arguments --> true substitution
            possibleSubstitutes.append(lits_diff_2to1(k)->args(1));
          }
        }
        if (DEBUG>2) {cout<<"2:  "<<possibleSubstitutes<<endl;}
        FOR1D(possibleSubstitutes, k) {
          FOR1D_(potential_subs, l) {
            potential_subs.elem(l)->addSubs(lits_diff_1to2(i)->args(1), possibleSubstitutes(k));
            in_ids.setAppend(lits_diff_1to2(i)->args(1));
          }
        }
      }
    }
  }
  
  if (DEBUG>0) {
    PRINT(in_ids);
  }
  
  
  // (3) RESUBSTITUTE ALL VARS which are not covered yet; e.g. falls nur 71-->70 bisher, fuehere dann auch 70-->71 ein
  FOR1D_(potential_subs, i) {
    if (DEBUG>0) {cout<<"Completing potential substitution: ";  potential_subs.elem(i)->write();  cout<<endl;}
//     CHECK(potential_subs.elem(i)->mapsToDistinct(), "should map to distinct");
    if (!potential_subs.elem(i)->mapsToDistinct())
      continue;
    
    uintA ins, outs;
    potential_subs.elem(i)->getIns(ins);
    potential_subs.elem(i)->getOuts(outs);
    uintA unsubstituted_outs;
    FOR1D(outs, k) {
      int idx = ins.findValue(outs(k));
      if (idx < 0) {
        unsubstituted_outs.append(outs(k));
      }
    }
    
    if (DEBUG>0) {PRINT(unsubstituted_outs);}
    while (unsubstituted_outs.N > 0) {
      potential_subs.elem(i)->getIns(ins);
      potential_subs.elem(i)->getOuts(outs);
      uintA free_ins;  // for backward substitution
      FOR1D(ins, k) {
        int idx = outs.findValue(ins(k));
        if (idx < 0) {
          free_ins.append(ins(k));
        }
      }
      CHECK(free_ins.N <= unsubstituted_outs.N, "");
      if (DEBUG>0) {PRINT(free_ins);}
      FOR1D_DOWN(free_ins, k) {
        potential_subs.elem(i)->addSubs(unsubstituted_outs(k), free_ins(k));
        unsubstituted_outs.remove(k);
      }
    }
    if (DEBUG>0) {
      cout<<"Final potential substitution: ";  potential_subs.elem(i)->write();  cout<<endl;
    }
  }
  
  
  // (4) TRY OUT POTENTIAL SUBS
  uint minDiff = 10000; // over all subsitutions
  FOR1D_(potential_subs, i) {
    // "potential_state2 should be state2
    SymbolicState* potential_state2 = potential_subs.elem(i)->apply(state1);
    if (DEBUG>0) {
      cout<<"Potential sub #" << i << ":"<<endl;
      potential_subs.elem(i)->write(cout);  cout<<endl;
      cout<<"Potential state 2 (built from state 1 using the substitution):"<<endl; potential_state2->write(cout, true); cout<<endl;
    }
    if (*potential_state2 == state2) {
      subs.append(potential_subs.elem(i));
      if (DEBUG>0) {cout<<"FITS!"<<endl;}
    }
    
    LitL local__pi_diff_1to2, local__pi_diff_2to1;
    uintA changedConstants;
    SymbolicState::calcDifferences(local__pi_diff_1to2, local__pi_diff_2to1, changedConstants, *potential_state2, state2);
    uint diff = local__pi_diff_1to2.N  +  local__pi_diff_2to1.N;
    
    if (DEBUG>0) {
      cout<<"Differences between potential state 2 and true state 2:"<<endl;
      cout<<"local__pi_diff_1to2:    "; write(local__pi_diff_1to2); cout<<endl;
      cout<<"local__pi_diff_2to1:    "; write(local__pi_diff_2to1); cout<<endl;
      PRINT(diff);
    }
    
    minDiff = TL_MIN(diff,  minDiff);
  }
  
  
  if (DEBUG>0) {
    cout<<"Final subs:"<<endl;
    FOR1D_(subs, i) {
      subs.elem(i)->write(cout);  cout<<endl;
    }
    PRINT(minDiff);
  }
  
  if (DEBUG>0) {cout<<"calcSubstitutions_asMuchAsPossible [END]"<<endl;}
  
  return minDiff;
}


uint reason::calcSubstitutions_asMuchAsPossible(SubstitutionSet& subs,
                                const SymbolicState& state1, const Literal& action1,
                                const SymbolicState& state2, const Literal& action2) {
  Substitution* action_sub = new Substitution;
  CHECK(action1.s == action2.s, "Try to unify different action symbols:  " << action2.s << " vs " << action1.s);
  uint i;
  FOR1D(action1.args, i) {
    action_sub->addSubs(action1.args(i), action2.args(i));
  }
  uint num_subs = calcSubstitutions_asMuchAsPossible(subs, state1, state2, action_sub);
  delete action_sub;
  return num_subs;
}


}
