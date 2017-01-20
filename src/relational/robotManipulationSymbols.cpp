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

#include "robotManipulationSymbols.h"
#include "reason.h"


namespace relational {


/************************************************
  * 
  *     Set up the symbols
  * 
  ************************************************/

// #define CLEARANCE
// #define TOWER

void RobotManipulationSymbols::initSymbols() {
  SymL symbols_state;
  // Primitives
  symbols_state.append(getSymbol_table());
  symbols_state.append(getSymbol_block());
  symbols_state.append(getSymbol_ball());
  symbols_state.append(getSymbol_on());
  symbols_state.append(getSymbol_inhand());
  symbols_state.append(getSymbol_upright());
  symbols_state.append(getSymbol_out());
  symbols_state.append(getFunction_size());
  // Derived
  symbols_state.append(getSymbol_clear());
  symbols_state.append(getSymbol_inhandNil());
  
  
  
  SymL symbols_action;
  symbols_action.append(getSymbol_action_grab());
  symbols_action.append(getSymbol_action_puton());
//   symbols_action.append(getSymbol_action_lift());
//   symbols_action.append(getSymbol_action_place());

  
#ifdef CLEARANCE
  symbols_state.append(getSymbol_above());    // CLEARANCE
  symbols_state.append(getSymbol_height());
#endif
  
#ifdef CLEARANCE
  symbols_state.append(getSymbol_homies());  // CLEARANCE
  symbols_state.append(getSymbol_above());    // CLEARANCE
  symbols_state.append(getSymbol_aboveNotable());    // CLEARANCE
  symbols_state.append(getSymbol_dirtyGuyBelow());    // CLEARANCE
  symbols_state.append(getSymbol_diffTower());   // CLEARANCE
  symbols_state.append(getSymbol_withoutHomies());    // CLEARANCE
  symbols_state.append(getSymbol_inorder());    // CLEARANCE
  
  symbols_state.append(getFunction_countInorder());    // CLEARANCE
#endif
  
#ifdef TOWER
  symbols_state.append(getSymbol_box());
  symbols_state.append(getSymbol_contains()); 
  symbols_state.append(getSymbol_closed()); 
  symbols_state.append(getSymbol_onBox());
  symbols_action.append(getSymbol_action_openBox());
  symbols_action.append(getSymbol_action_closeBox());
#endif
}



/************************************************
  * 
  *     Actions
  * 
  ************************************************/


Symbol* RobotManipulationSymbols::getSymbol_action_grab() {
  return Symbol::get(MT::String("grab"), 1, Symbol::action);
}

Symbol* RobotManipulationSymbols::getSymbol_action_puton() {
  return Symbol::get(MT::String("puton"), 1, Symbol::action);
}

Symbol* RobotManipulationSymbols::getSymbol_action_lift() {
  return Symbol::get(MT::String("lift"), 2, Symbol::action);
}

Symbol* RobotManipulationSymbols::getSymbol_action_place() {
  return Symbol::get(MT::String("place"), 2, Symbol::action);
}

Symbol* RobotManipulationSymbols::getSymbol_action_openBox() {
  return Symbol::get(MT::String("openBox"), 1, Symbol::action);
}

Symbol* RobotManipulationSymbols::getSymbol_action_closeBox() {
  return Symbol::get(MT::String("closeBox"), 1, Symbol::action);
}



/************************************************
  * 
  *     Primitives
  * 
  ************************************************/


Symbol* RobotManipulationSymbols::getSymbol_on() {
  return Symbol::get(MT::String("on"), 2);
}

Symbol* RobotManipulationSymbols::getSymbol_table() {
  return Symbol::get(MT::String("table"), 1);
}

Symbol* RobotManipulationSymbols::getSymbol_block() {
  return Symbol::get(MT::String("block"), 1);
}

Symbol* RobotManipulationSymbols::getSymbol_box() {
  return Symbol::get(MT::String("box"), 1);
}

Symbol* RobotManipulationSymbols::getSymbol_ball() {
  return Symbol::get(MT::String("ball"), 1);
}

Symbol* RobotManipulationSymbols::getSymbol_inhand() {
  return Symbol::get(MT::String("inhand"), 1);
}

Symbol* RobotManipulationSymbols::getSymbol_upright() {
  return Symbol::get(MT::String("upright"), 1);
}

Symbol* RobotManipulationSymbols::getSymbol_out() {
  return Symbol::get(MT::String("out"), 1);
}

Symbol* RobotManipulationSymbols::getSymbol_homies() {
  return Symbol::get(MT::String("homies"), 2);
}

Symbol* RobotManipulationSymbols::getSymbol_contains() {
  return Symbol::get(MT::String("contains"), 2);
}

Symbol* RobotManipulationSymbols::getSymbol_closed() {
  return Symbol::get(MT::String("closed"), 1);
}


Symbol* RobotManipulationSymbols::getFunction_size() {
  return Symbol::get(MT::String("size"), 1, Symbol::primitive, Symbol::integers);
}






/************************************************
  * 
  *     Derived
  * 
  ************************************************/

ConjunctionSymbol* RobotManipulationSymbols::getSymbol_clear() {
  uintA args;  args.append(1);  args.append(0);
  LitL base_literals;  base_literals.append(Literal::get(getSymbol_on(), args, 0.));
  return ConjunctionSymbol::get(MT::String("clear"), 1, base_literals, true);
}


TransClosureSymbol* RobotManipulationSymbols::getSymbol_above() {
  return TransClosureSymbol::get(MT::String("above"), 2, getSymbol_on());
}


ConjunctionSymbol* RobotManipulationSymbols::getSymbol_aboveNotable() {
  LitL base_literals;
  uintA args;
  args.append(0);  args.append(1);
  base_literals.append(Literal::get(getSymbol_above(), args, 1.));
  args.clear();  args.append(1);
  base_literals.append(Literal::get(getSymbol_table(), args, 0.));
  return ConjunctionSymbol::get(MT::String("aboveNotable"), 2, base_literals, true);
}


ConjunctionSymbol* RobotManipulationSymbols::getSymbol_dirtyGuyBelow() {
    NIY;
#if 0
    p_DIRTY_GUY_BELOW= new ConjunctionSymbol;
    p_DIRTY_GUY_BELOW->name = "dirtyGuyBelow";
//     p_DIRTY_GUY_BELOW->id = HAND_ID__PRED_DIRTY_GUY_BELOW;
    p_DIRTY_GUY_BELOW->arity = 1;
    p_DIRTY_GUY_BELOW->freeVarsAllQuantified = false;
    p_DIRTY_GUY_BELOW->base_symbols.resize(3);
    p_DIRTY_GUY_BELOW->base_symbols(0) = getSymbol_table();
    p_DIRTY_GUY_BELOW->base_symbols(1) = getSymbol_homies();
    p_DIRTY_GUY_BELOW->base_symbols(2) = getSymbol_on();
    p_DIRTY_GUY_BELOW->basePreds_mapVars2conjunction.resize(5);
    p_DIRTY_GUY_BELOW->basePreds_mapVars2conjunction(0) = 1;
    p_DIRTY_GUY_BELOW->basePreds_mapVars2conjunction(1) = 0;
    p_DIRTY_GUY_BELOW->basePreds_mapVars2conjunction(2) = 1;
    p_DIRTY_GUY_BELOW->basePreds_mapVars2conjunction(3) = 0;
    p_DIRTY_GUY_BELOW->basePreds_mapVars2conjunction(4) = 1;
    p_DIRTY_GUY_BELOW->base_symbols_positive.resize(3);
    p_DIRTY_GUY_BELOW->base_symbols_positive(0) = false;
    p_DIRTY_GUY_BELOW->base_symbols_positive(1) = false;
    p_DIRTY_GUY_BELOW->base_symbols_positive(2) = true;
#endif
}


ConjunctionSymbol* RobotManipulationSymbols::getSymbol_diffTower() {
    NIY;
#if 0
    p_DIFF_TOWER = new ConjunctionSymbol;
    p_DIFF_TOWER->name = "diffTower";
//     p_DIFF_TOWER->id = HAND_ID__PRED_DIFF_TOWER;
    p_DIFF_TOWER->arity = 2;
    p_DIFF_TOWER->base_symbols.resize(2);
    p_DIFF_TOWER->base_symbols(0) = getSymbol_above();
    p_DIFF_TOWER->base_symbols(1) = getSymbol_above();
    p_DIFF_TOWER->basePreds_mapVars2conjunction.resize(4);
    p_DIFF_TOWER->basePreds_mapVars2conjunction(0) = 0;
    p_DIFF_TOWER->basePreds_mapVars2conjunction(1) = 1;
    p_DIFF_TOWER->basePreds_mapVars2conjunction(2) = 1;
    p_DIFF_TOWER->basePreds_mapVars2conjunction(3) = 0;
    p_DIFF_TOWER->base_symbols_positive.resize(2);
    p_DIFF_TOWER->base_symbols_positive(0) = false;
    p_DIFF_TOWER->base_symbols_positive(1) = false;
    p_DIFF_TOWER->freeVarsAllQuantified = false;
#endif
}


ConjunctionSymbol* RobotManipulationSymbols::getSymbol_withoutHomies() {
    NIY;
#if 0
    p_WITHOUT_HOMIES = new ConjunctionSymbol;
    p_WITHOUT_HOMIES->name = "withoutHomies";
//     p_WITHOUT_HOMIES->id = HAND_ID__PRED_WITHOUT_HOMIES;
    p_WITHOUT_HOMIES->arity = 1;
    p_WITHOUT_HOMIES->base_symbols.resize(2);
    p_WITHOUT_HOMIES->base_symbols(0) = getSymbol_homies();
    p_WITHOUT_HOMIES->base_symbols(1) = getSymbol_diffTower();
    p_WITHOUT_HOMIES->basePreds_mapVars2conjunction.resize(4);
    p_WITHOUT_HOMIES->basePreds_mapVars2conjunction(0) = 0;
    p_WITHOUT_HOMIES->basePreds_mapVars2conjunction(1) = 1;
    p_WITHOUT_HOMIES->basePreds_mapVars2conjunction(2) = 0;
    p_WITHOUT_HOMIES->basePreds_mapVars2conjunction(3) = 1;
    p_WITHOUT_HOMIES->base_symbols_positive.resize(2);
    p_WITHOUT_HOMIES->base_symbols_positive(0) = true;
    p_WITHOUT_HOMIES->base_symbols_positive(1) = true;
    p_WITHOUT_HOMIES->freeVarsAllQuantified = false;
#endif
}


ConjunctionSymbol* RobotManipulationSymbols::getSymbol_inorder() {
    NIY;
#if 0
    p_INORDER = new ConjunctionSymbol;
    p_INORDER->name = "inorder";
//     p_INORDER->id = HAND_ID__PRED_INORDER;
    p_INORDER->arity = 1;
    NIY;
    p_INORDER->base_symbols.resize(3);
    p_INORDER->base_symbols(0) = getSymbol_table();
    p_INORDER->base_symbols(1) = getSymbol_withoutHomies();
    p_INORDER->base_symbols(2) = getSymbol_dirtyGuyBelow();
    p_INORDER->basePreds_mapVars2conjunction.resize(3);
    p_INORDER->basePreds_mapVars2conjunction(0) = 0;
    p_INORDER->basePreds_mapVars2conjunction(1) = 0;
    p_INORDER->basePreds_mapVars2conjunction(2) = 0;
    p_INORDER->base_symbols_positive.resize(3);
    p_INORDER->base_symbols_positive(0) = false;
    p_INORDER->base_symbols_positive(1) = false;
    p_INORDER->base_symbols_positive(2) = false;
#endif
}


ConjunctionSymbol* RobotManipulationSymbols::getSymbol_inhandNil() {
  uintA args;  args.append(0);
  LitL base_literals;  base_literals.append(Literal::get(getSymbol_inhand(), args, 0.));
  return ConjunctionSymbol::get(MT::String("inhandNil"), 0, base_literals, true);
}


ConjunctionSymbol* RobotManipulationSymbols::getSymbol_onBox() {
  NIY;
  return NULL;
#if 0
    p_ON_BOX = new ConjunctionSymbol;
    p_ON_BOX->arity = 1;
    p_ON_BOX->name = "onBox";
//     p_ON_BOX->id = HAND_ID__PRED_ON_BOX;
    p_ON_BOX->base_symbols.resize(2);
    p_ON_BOX->base_symbols(0) = getSymbol_on();
    p_ON_BOX->base_symbols(1) = getSymbol_box();
    p_ON_BOX->basePreds_mapVars2conjunction.resize(3);
    p_ON_BOX->basePreds_mapVars2conjunction(0) = 0;
    p_ON_BOX->basePreds_mapVars2conjunction(1) = 1;
    p_ON_BOX->basePreds_mapVars2conjunction(2) = 1;
    p_ON_BOX->base_symbols_positive.resize(2);
    p_ON_BOX->base_symbols_positive(0) = true;
    p_ON_BOX->base_symbols_positive(1) = true;
    p_ON_BOX->freeVarsAllQuantified = false;
#endif
}


CountSymbol* RobotManipulationSymbols::getSymbol_height() {
  uintA args;  args.append(0);  args.append(1);
  Literal* base_literal = Literal::get(getSymbol_aboveNotable(), args, 1);
  return CountSymbol::get(MT::String("height"), 1, base_literal);
}


AverageFunction* RobotManipulationSymbols::getFunction_avgheight() {
  NIY;
//   if (f_AVG_HEIGHT == NULL) {
//     f_AVG_HEIGHT = new AverageFunction;
//     f_AVG_HEIGHT->arity = 0;
//     f_AVG_HEIGHT->name = "avg_height";
//   //     f_AVG_HEIGHT->id = HAND_ID__FUNCTION_AVG_HEIGHT;
//     f_AVG_HEIGHT->base_symbol = getSymbol_height();
//   }
//   return f_AVG_HEIGHT;
}


SumFunction* RobotManipulationSymbols::getFunction_sumheight() {
  return SumFunction::get(MT::String("sum_height"), 0, getSymbol_height());
}


CountSymbol* RobotManipulationSymbols::getFunction_countInorder() {
  NIY;
    #if 0
    f_COUNT_INORDER = new CountSymbol;
//     f_COUNT_INORDER->id = HAND_ID__FUNCTION_COUNT_INORDER;
    f_COUNT_INORDER->name = "count_inorder";
    f_COUNT_INORDER->arity = 0;
    uint i;
    for (i=0;i<20;i++) {f_COUNT_INORDER->range.append(i);}
    f_COUNT_INORDER->counted_symbol = getSymbol_inorder();
    f_COUNT_INORDER->counted_symbol_mapVars2derived.resize(1);
    f_COUNT_INORDER->counted_symbol_mapVars2derived(0) = 0;
#endif
}


    
	









/************************************************
  * 
  *     RewardLibrary
  * 
  ************************************************/

Reward* RobotManipulationSymbols::RewardLibrary::on(uint o1, uint o2) {
  Symbol* p_ON = Symbol::get(MT::String("on"));
  uintA sa2(2);
  sa2(0)=o1;
  sa2(1)=o2;
  Literal* pt = Literal::get(p_ON, sa2, 1.);
  return new LiteralReward(pt);
}


Reward* RobotManipulationSymbols::RewardLibrary::inhand(uint o1) {
  Symbol* p_INHAND = Symbol::get(MT::String("inhand"));
  uintA sa1(1);
  sa1(0)=o1;
  Literal* pt = Literal::get(p_INHAND, sa1, 1.);
  return new LiteralReward(pt);
}


Reward* RobotManipulationSymbols::RewardLibrary::stack() {
//   ATTENTION:  Stack defined by sum (not max) over heights
  uintA empty;
  SymL symbols2add;
  if (!Symbol::get(MT::String("above"))) {
    TransClosureSymbol* p_ABOVE1 = getSymbol_above();
    p_ABOVE1->base_symbol = Symbol::get(MT::String("on")); // HACK da in regelfiles bis juni 2009 on andere id hat
    symbols2add.append(p_ABOVE1);
  }
  if (!Symbol::get(MT::String("aboveNotable"))) {
    symbols2add.append(getSymbol_aboveNotable());
  }
  if (!Symbol::get(MT::String("height"))) {
    symbols2add.append(getSymbol_height());
  }
  if (!Symbol::get(MT::String("sum_height"))) {
    symbols2add.append(getFunction_sumheight());
  }
  Literal* fa = Literal::get(Symbol::get(MT::String("sum_height")), empty, 1.);
  Reward* reward = new MaximizeReward(fa);
  return reward;
}


Reward* RobotManipulationSymbols::RewardLibrary::tower(uintA& objs) {
  Symbol* p_ON = Symbol::get(MT::String("on"));
  LitL pts;
  uint i;
  uintA sa2(2);
  FOR1D(objs, i) {
    sa2(0)=objs(i);
    if (i<objs.N-1)
      sa2(1)=objs(i+1);
    else
      sa2(1)=60;  // table id in my ors simulator
    pts.append(Literal::get(p_ON, sa2, 1.));
  }
  LiteralListReward* reward = new LiteralListReward(pts);
  return reward;
}


Reward* RobotManipulationSymbols::RewardLibrary::clearance() {
  SymL symbols2add;
  if (!Symbol::get(MT::String("above"))) {
    TransClosureSymbol* p_ABOVE1 = getSymbol_above();
    p_ABOVE1->base_symbol = Symbol::get(MT::String("on")); // HACK da in regelfiles bis juni 2009 on andere id hat
    symbols2add.append(p_ABOVE1);
  }
  if (!Symbol::get(MT::String("dirtyGuyBelow"))) {
    symbols2add.append(getSymbol_dirtyGuyBelow());
  }
  if (!Symbol::get(MT::String("differentTower"))) {
    symbols2add.append(getSymbol_diffTower());
  }
  if (!Symbol::get(MT::String("withoutHomies"))) {
    symbols2add.append(getSymbol_withoutHomies());
  }
  if (!Symbol::get(MT::String("inorder"))) {
    symbols2add.append(getSymbol_inorder());
  }
  if (!Symbol::get(MT::String("count_inorder"))) {
    symbols2add.append(getFunction_countInorder());
  }
  uintA empty;
  Literal* fi = Literal::get(getFunction_countInorder(), empty, 1.);
  Reward* reward = new MaximizeReward(fi);
  return reward;
}




/************************************************
  * 
  *     Helpers for SymbolicState
  * 
  ************************************************/

bool RobotManipulationSymbols::isBlock(uint id, const SymbolicState& state) {
  return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("block")), TUP(id), 1.0));
}


bool RobotManipulationSymbols::isOut(uint id, const SymbolicState& state) {
  return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("out")), TUP(id), 1.0));
}

bool RobotManipulationSymbols::isInhand(uint id, const SymbolicState& state) {
  return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("inhand")), TUP(id), 1.0));
}

bool RobotManipulationSymbols::isTable(uint id, const SymbolicState& state) {
  return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("table")), TUP(id), 1.0));
}

bool RobotManipulationSymbols::isBall(uint id, const SymbolicState& state) {
  return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("ball")), TUP(id), 1.0));
}

bool RobotManipulationSymbols::isBox(uint id, const SymbolicState& state) {
  if (Symbol::get(MT::String("box")) == NULL)
    return false;
  return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("box")), TUP(id), 1.0));
}

bool RobotManipulationSymbols::isClosed(uint id, const SymbolicState& state) {
  return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("closed")), TUP(id), 1.0));
}

bool RobotManipulationSymbols::isInorderGang(const uintA gang, const SymbolicState& state) {
  CHECK(gang.N > 0, "");
  if (Symbol::get(MT::String("inorder")) == NULL) {
    NIY;
  }
  else {
    return reason::holds(state.lits, Literal::get(Symbol::get(MT::String("inorder")), gang(0), 1.0));
  }
}


uint RobotManipulationSymbols::getBelow(uint id, const SymbolicState& state) {
  Symbol* p_ON = Symbol::get(MT::String("on"));
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s == p_ON) {
      if (state.lits(i)->args(0) == id)
        return state.lits(i)->args(1);
    }
  }
  return TL::UINT_NIL;
}

uint RobotManipulationSymbols::getAbove(uint id, const SymbolicState& state) {
  Symbol* p_ON = Symbol::get(MT::String("on"));
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s == p_ON) {
      if (state.lits(i)->args(1) == id)
        return state.lits(i)->args(0);
    }
  }
  return TL::UINT_NIL;
}

void RobotManipulationSymbols::getBelowObjects(uintA& ids, uint id_top, const SymbolicState& state) {
  SymbolicState::getRelatedConstants(ids, id_top, true, *Symbol::get(MT::String("above")), state);
}

void RobotManipulationSymbols::getAboveObjects(uintA& ids, uint id_top, const SymbolicState& state) {
  SymbolicState::getRelatedConstants(ids, id_top, false, *Symbol::get(MT::String("above")), state);
}

uint RobotManipulationSymbols::getInhand(const SymbolicState& state) {
  Symbol* p_INHAND = Symbol::get(MT::String("inhand"));
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s == p_INHAND) {
      return state.lits(i)->args(0);
    }
  }
  return TL::UINT_NIL;
}

void RobotManipulationSymbols::getBoxes(uintA& ids, const SymbolicState& state) {
  ids.clear();
  Symbol* p_BOX = Symbol::get(MT::String("box"));
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s == p_BOX) {
      ids.append(state.lits(i)->args(0));
    }
  }
}


uint RobotManipulationSymbols::getContainingBox(uint obj_id, const SymbolicState& state) {
  Symbol* p_CONTAINS = Symbol::get(MT::String("contains"));
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s == p_CONTAINS) {
      if (state.lits(i)->args(1) == obj_id)
        return state.lits(i)->args(0);
    }
  }
  return TL::UINT_NIL;
}


uint RobotManipulationSymbols::getContainedObject(uint box_id, const SymbolicState& state) {
  Symbol* p_CONTAINS = Symbol::get(MT::String("contains"));
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s == p_CONTAINS) {
      if (state.lits(i)->args(0) == box_id)
        return state.lits(i)->args(1);
    }
  }
  return TL::UINT_NIL;
}


void RobotManipulationSymbols::getHomieGangs(MT::Array< uintA >& homieGangs, const SymbolicState& state) {
  homieGangs.clear();
  Symbol* p_HOMIES = Symbol::get(MT::String("homies"));
  if (p_HOMIES == NULL)
    return;
  uint i, k;
  boolA constants_covered(reason::getConstants().N);
  constants_covered.setUni(false);
  FOR1D(reason::getConstants(), i) {
    if (reason::holds(state.lits, Literal::get(Symbol::get(MT::String("table")), reason::getConstants()(i), 1.0)))
      continue;
    if (constants_covered(i))
      continue;
    uintA homies;
    SymbolicState::getRelatedConstants(homies, reason::getConstants()(i), true, *p_HOMIES, state);
    homies.insert(0, reason::getConstants()(i));
    constants_covered(i) = true;
    FOR1D(homies, k) {
      constants_covered(reason::getConstants().findValue(homies(k))) = true;
    }
    homieGangs.append(homies);
  }
}





// Helpers for reward functions

double RobotManipulationSymbols::reward_buildTower(const SymbolicState& state) {
  uint DEBUG=0;
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::reward_buildTower [START]"<<endl;}
  uintA piles;
  uint id_table = TL::UINT_NIL;
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s->name == "table"  &&  state.lits(i)->value>0.) {
      id_table = state.lits(i)->args(0);
      break;
    }
  }
  CHECK(i<state.lits.N, "table id not found");
  calcPiles(state, piles, id_table);
  uint height;
  for(height=0; height<piles.d1; height++) {
    if (piles(0, height)==TL::UINT_NIL)
      break;
  }
  double reward = 1.0 * height;
  if (DEBUG>0) {
    PRINT(piles)
        PRINT(reward)
  }
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::reward_buildTower [END]"<<endl;}
  return reward;
}


void RobotManipulationSymbols::calcPiles(const SymbolicState& state, uintA& piles, uint sort_type) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<" RobotManipulationSymbols::calcPiles [START]"<<endl;}
  if (DEBUG>0) {cout<<"STATE:  ";  state.write(cout, true);  cout<<endl;}
  bool use_table = false;
  uint i;
  FOR1D(state.state_constants, i) {
    if (isTable(state.state_constants(i), state)) {
      use_table = true;
      break;
    }
  }
  // calc piles (unsorted)
  MT::Array< uintA > piles_unsorted;
  bool inserted;
  FOR1D(state.lits, i) {
    // on(A,B)
    if (state.lits(i)->s == Symbol::get(MT::String("on"))) {
      inserted = false;
      uint A_2_top;
      FOR1D(piles_unsorted, A_2_top) {
        if (piles_unsorted(A_2_top).N == 0)
          continue;
        // pile with [A, ..., top]  -->  put B below
        if (piles_unsorted(A_2_top)(0) == state.lits(i)->args(0)) {
          piles_unsorted(A_2_top).insert(0, state.lits(i)->args(1));
          inserted = true;
        }
      }
      uint table_2_B;
      FOR1D(piles_unsorted, table_2_B) {
        if (piles_unsorted(table_2_B).N == 0)
          continue;
        // pile with [table, lastBlock, ..., B]  -->  put A on top
        if (piles_unsorted(table_2_B).last() == state.lits(i)->args(1)) {
          if (inserted) {
            // when trying to insert a second time
            // find previous insertion point and combine both piles
            FOR1D(piles_unsorted, A_2_top) {
              if (piles_unsorted(A_2_top).N == 0)
                continue;
              // pile with [A, ..., top]  -->  put B below
              if (piles_unsorted(A_2_top)(0) == state.lits(i)->args(1)) {  // anderer check als oben, verdammt, da B ja jetzt schon eingefuegt!
                break;
              }
            }
            if (A_2_top == piles_unsorted.N) {
              cerr<<endl<<endl<<endl<<endl;
              MT_MSG("ERROR! A_2_top was not found. Sorting piles will fail.");
              PRINT2(A_2_top, cerr);
              PRINT2(piles_unsorted.N, cerr);
              PRINT2(piles_unsorted, cerr);
              PRINT2(use_table, cerr);
              PRINT2(*state.lits(i), cerr);
              cerr<<endl<<endl<<endl<<endl;
              break;
            }
            piles_unsorted(table_2_B).setAppend(piles_unsorted(A_2_top)); // schmeisst doppeleintrag raus
            piles_unsorted(A_2_top).clear();
          }
          else {
            piles_unsorted(table_2_B).append(state.lits(i)->args(0));
            inserted = true;
          }
        }
      }
      if (!inserted) {
        uintA newPile;
        newPile.append(state.lits(i)->args(1));
        newPile.append(state.lits(i)->args(0));
        piles_unsorted.append(newPile);
      }
    }
  }
  
  MT::Array< uintA > piles_unsorted2;
  FOR1D(piles_unsorted, i) {
    if (piles_unsorted(i).N > 0) {
      piles_unsorted2.append(piles_unsorted(i));
    }
  }
  piles_unsorted = piles_unsorted2;
  
  uint j;
    
  if (DEBUG>0) {
    cout <<"piles_unsorted (N="<<piles_unsorted.N<<"):"<<endl;
    PRINT(piles_unsorted);
    FOR1D(piles_unsorted, i) {
      FOR1D(piles_unsorted(i), j) {
        cout << piles_unsorted(i)(j) << " ";
      }
      cout << endl;
    }
  }
  if (piles_unsorted.N == 0) {
    if (DEBUG>0) {cout<<"No piles found."<<endl;}
    if (DEBUG>0) {cout<<" RobotManipulationSymbols::calcPiles [END]"<<endl;}
    return;
  }
  
  uintA heights;
  FOR1D(piles_unsorted, i) {
    heights.append(piles_unsorted(i).d0);
  }
  if (DEBUG>0) {
    FOR1D(heights, i) {
      PRINT(heights(i))
    }
  }
  uintA sortedIndices;
  TL::sort_desc_keys(sortedIndices, heights); // descending
  piles.resize(piles_unsorted.d0, heights(sortedIndices(0)));
  piles.setUni(TL::UINT_NIL);
  
  if (DEBUG>0) {PRINT(sort_type);}
  // sort by height
  if (sort_type == 1) {
    // reorder piles
    FOR1D(piles, i) {
      FOR1D(piles_unsorted(sortedIndices(i)), j) {
        piles(i,j) = piles_unsorted(sortedIndices(i))(j);
      }
    }
  }
  // sort by id
  else if (sort_type == 2) {
    uintA args_state;
    Literal::getArguments(args_state, state.lits);
    if (DEBUG>0) {PRINT(args_state);}
    uint next_pile_id = 0;
    FOR1D(args_state, i) {
      FOR1D(piles_unsorted, j) {
        uint first_inspected_id = 0;
        if (use_table)
          first_inspected_id = 1;
        if (piles_unsorted(j)(first_inspected_id) == args_state(i)) {   // (1) for table
          uint k;
          FOR1D(piles_unsorted(j), k) {
            piles(next_pile_id,k) = piles_unsorted(j)(k);
          }
          next_pile_id++;
          break;
        }
      }
    }
    if (next_pile_id != piles_unsorted.N) {
      cerr<<endl<<endl<<endl<<endl;
      PRINT(piles_unsorted);
      PRINT(piles);
      PRINT(args_state);
      PRINT(next_pile_id);
      MT_MSG("resorting failed");
      cerr<<endl<<endl<<endl<<endl;
      if (DEBUG>0) {cout<<" RobotManipulationSymbols::calcPiles [END]"<<endl;}
      return;
    }
  }
  else {
    FOR1D(piles_unsorted, i) {
      FOR1D(piles_unsorted(i), j) {
        piles(i,j) = piles_unsorted(i)(j);
      }
    }
  }
  
  if (DEBUG>0) {PRINT(piles);}
  if (DEBUG>0) {cout<<" RobotManipulationSymbols::calcPiles [END]"<<endl;}
}


// unfound objects get height = 0
void RobotManipulationSymbols::calcHeights(const uintA& objects, const uintA& piles, uintA& object_heights, uint id_table) {
  uint i, k, l;
  bool found;
  object_heights.resize(objects.N);
  FOR1D(objects, i) {
    found = false;
    FOR2D(piles, k, l) {
      if (piles(k,l)==objects(i)) {
        object_heights(i) = l+1;
        found = true;
        break;
      }
    }
    if (!found) {
      object_heights(i) = 0;
    }
  }
}


bool RobotManipulationSymbols::has_maximum_stack_value(const SymbolicState& state) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::has_maximum_stack_value [START]"<<endl;}
  if (DEBUG>0) {cout<<"STATE:   ";  state.write();  cout<<endl;}
  
  Symbol* p_BLOCK = Symbol::get(MT::String("block"));
  Symbol* p_BALL = Symbol::get(MT::String("ball"));
  Symbol* p_BOX = Symbol::get(MT::String("box"));
  
  uintA blocks, balls, boxes;
  SymbolicState::getArguments(blocks, state, *p_BLOCK);
  SymbolicState::getArguments(balls, state, *p_BALL);
  SymbolicState::getArguments(boxes, state, *p_BOX);
  
  if (DEBUG>0) {PRINT(blocks);  PRINT(balls);  PRINT(boxes);}
  
  uint i;
  double maximum_stack_value = 0.;
  if (boxes.N > 0) {
    FOR1D(blocks, i) {
      maximum_stack_value += i+1;
    }
  }
  else {
    FOR1D(blocks, i) {
      maximum_stack_value += i;
    }
  }
  
  FOR1D(balls, i) {
    // erster ball auf die bloecke
    if (i==0) {
      maximum_stack_value += blocks.N;
    }
    // weitere baelle auf die kisten
    else if (boxes.N > i-1) {
      maximum_stack_value += 1;
    }
    else break;
  }
  
  Symbol* f_SUM_HEIGHT = Symbol::get(MT::String("sum_height"));
  double real_stack_value = SymbolicState::getValue(f_SUM_HEIGHT, state);
  
  if (DEBUG>0) {PRINT(maximum_stack_value);  PRINT(real_stack_value);}
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::has_maximum_stack_value [END]"<<endl;}
  return TL::areEqual(maximum_stack_value, real_stack_value);
}


void RobotManipulationSymbols::calcSkyscraperWeights(const uintA& heights, double skyscraper_bias, arr& weights, bool highGood, uint id_table) {
  uint i;
  weights.resize(heights.N);
  weights.setUni(1.0);
  FOR1D(heights, i) {
    if (highGood)
      weights(i) += heights(i) * skyscraper_bias;
    else
      weights(i) -= heights(i) * skyscraper_bias;
  }
  if (!highGood) {
    double minValue = weights.min();
    FOR1D(weights, i) {
      weights(i) -= minValue - 1.0;
    }
  }
}


void RobotManipulationSymbols::calcSkyscraperWeights(const uintA& objects, const uintA& piles, double skyscraper_bias, arr& weights, bool highGood, uint id_table) {
  uintA object_heights;
//     PRINT(piles)
  calcHeights(objects, piles, object_heights, id_table);
  calcSkyscraperWeights(object_heights, skyscraper_bias, weights, highGood, id_table);
    // special care for table if putting on
  if (highGood) {
    uint i;
    FOR1D(objects, i) {
      if (objects(i)==id_table)
        break;
    }
    if (i<objects.N)
      weights(i) = skyscraper_bias;
  }
}


LiteralListReward* RobotManipulationSymbols::sampleGroundGoal__stack(const uintA& blocks, const uintA& balls, uint table_id, bool leave_existing_towers, SymbolicState* state) {
  uint DEBUG = 1;
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::sampleGroundGoal__stack [START]"<<endl;}
  if (DEBUG>1) {PRINT(blocks);  PRINT(balls);  PRINT(table_id);}
  uintA objs;
  objs.append(blocks);
  objs.append(balls);
  
  // -----------------------------------------------------------------------------
  // (1) Determine heights of towers
  
  arr weights_tower_heights;
  uint i, k;
  double basis = 1.2;
//   if (objs.N >= 10)
//     basis = 1.1;
  
  for (i=2; /*i<objs.N*/ i<5; i++) {
    weights_tower_heights.append(pow(basis, i-2));
  }
  weights_tower_heights /= sum(weights_tower_heights);
  uint target_height = TL::basic_sample(weights_tower_heights) + 2;
  if (DEBUG>1) {PRINT(weights_tower_heights);  PRINT(target_height);}
  
  
  // -----------------------------------------------------------------------------
  // (2) Fill towers with objects
  
  MT::Array< uintA > towers;
  if (!leave_existing_towers) {
    if (DEBUG>0) {cout<<"Dumb sampling"<<endl;}
    for (i=0; (i+1)*target_height <= objs.N; i++) {
      uintA tower;
      for (k=0; k<target_height; k++) {
        uint block;
        do {
          block = rnd.num(objs.N);
        }
        while (objs(block) == TL::UINT_NIL);
        tower.append(objs(block));
        objs(block) = TL::UINT_NIL;
      }
      tower.insert(0, table_id);
      towers.append(tower);
    }
    
    uintA last_tower;
    last_tower.append(table_id);
  
    FOR1D(objs, i) {
      if (objs(i) != TL::UINT_NIL)
        last_tower.append(objs(i));
    }
    if (last_tower.N > 1) {
      towers.append(last_tower);
    }
  }
  else {
    if (DEBUG>0) {cout<<"Preserving sampling"<<endl;}
    uintA piles;
    calcPiles(*state, piles, table_id);
    if (DEBUG>1) {PRINT(piles);}
    // leave high towers
    uintA used_objects;
    for (i=0; i<piles.d0; i++) {
      k=0;
      while (k < piles.d1 && piles(i,k) != TL::UINT_NIL) {
        k++;
      }
      if (k > 2) {
        uintA true_tower;
        k = 0;
        while (k < piles.d1 && piles(i,k) != TL::UINT_NIL) {
          true_tower.append(piles(i,k));
          used_objects.setAppend(piles(i,k));
          k++;
        }
        towers.append(true_tower);
      }
    }
    
    // randomly split other objects on top of these
    uintA remaining_objs;
    remaining_objs = objs;
    setMinus<uint>(remaining_objs, used_objects);
    if (DEBUG>1) {PRINT(used_objects);  PRINT(remaining_objs);}
    while (remaining_objs.N > 0) {
      uint obj = remaining_objs(rnd.num(remaining_objs.N));
      // lower chance for balls
      if (balls.findValue(obj) >= 0) {
        obj = remaining_objs(rnd.num(remaining_objs.N));
        if (balls.N >= 3 && balls.findValue(obj) >= 0) {
          obj = remaining_objs(rnd.num(remaining_objs.N));
            if (balls.N >= 4 && balls.findValue(obj) >= 0) {
              obj = remaining_objs(rnd.num(remaining_objs.N));
            }
        }
      }
      FOR1D(towers, i) {
        if (towers(i).N - 1 < target_height) {
          towers(i).append(obj);
          break;
        }
      }
      if (i == towers.N) {
        uintA new_tower;
        new_tower.append(table_id);
        new_tower.append(obj);
        towers.append(new_tower);
        if (DEBUG>1) {PRINT(new_tower);}
      }
      remaining_objs.removeValue(obj);
    }
  }

  if (DEBUG>0) {
    cout<<"****************"<<endl;
    cout<<"TOWERS:"<<endl;
    FOR1D(towers, i) {
      FOR1D(towers(i), k) {
        cout<<" ";
        cout<<towers(i)(k);
        if (balls.findValue(towers(i)(k)) >= 0) {
          cout<<"o";
        }
      }
      cout << endl;
    }
    cout<<"****************"<<endl;
    
    cerr<<"****************"<<endl;
    cerr<<"TOWERS:"<<endl;
    FOR1D(towers, i) {
      FOR1D(towers(i), k) {
        cerr<<" ";
        cerr<<towers(i)(k);
        if (balls.findValue(towers(i)(k)) >= 0) {
          cerr<<"o";
        }
      }
      cerr << endl;
    }
    cerr<<"****************"<<endl;
  }
  
  // -----------------------------------------------------------------------------
  // (3) Translate to logic
  
  LitL lits;
  FOR1D(towers, i) {
    FOR1D(towers(i), k) {
      if (k==0)
        continue;
      lits.append(Literal::get(Symbol::get(MT::String("on")), TUP(towers(i)(k), towers(i)(k-1)), 1.));
    }
  }
  
  if (DEBUG>0) {cout<<"LOGIC:"<<endl;  write(lits); cout<<endl;}
  
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::sampleGroundGoal__stack [END]"<<endl;}
  return new LiteralListReward(lits);
}


LiteralListReward* RobotManipulationSymbols::sampleGroundGoal__clearance(const SymbolicState& current_state, uint table_id) {
  uint DEBUG = 1;
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::sampleGroundGoal__clearance [START]"<<endl;}

  uint i, k;

  //-----------------------------------------------------------------------------
  // (1) Determine classes to be ordered
  
  MT::Array< uintA > gangs;
  getHomieGangs(gangs, current_state);

  // Determine gangs (= gangs of homies)
  if (DEBUG>1) {
    cout<<"GANGS:"<<endl;
    FOR1D(gangs, i) {
      cout << i << ":  " << gangs(i) << endl;
    }
  }
  
  uintA gang_ids;
  FOR1D(gangs, i) {
    if (!isInorderGang(gangs(i), current_state))
      gang_ids.append(i);
  }
  MT::Array< uintA > subsets;
  TL::allSubsets(subsets, gang_ids, false, false);
  
  arr sample_weights(subsets.N);
  FOR1D(subsets, i) {
    sample_weights(i) = subsets(i).N;
    if (sample_weights(i) > 2) // HACK 
      sample_weights(i) = 2;
  }
  sample_weights /= sum(sample_weights);
  
  uintA combo = subsets(TL::basic_sample(sample_weights));
  
  
  if (DEBUG>1) {
    PRINT(subsets);
    PRINT(sample_weights);
    PRINT(combo);
  }
  
  //-----------------------------------------------------------------------------
  // (2) Stack them together
  
  MT::Array< uintA > towers;
  FOR1D(combo, i) {
    uintA tower;
    uintA& local_gang = gangs(combo(i));
    FOR1D(local_gang, k) {
      uint candidate_position = rnd.num(tower.N+1);
      // lower chance for balls to be inside
//       if (isBall(local_gang(k), current_state)  &&  candidate_position != tower.N)
//         candidate_position = rnd.num(tower.N+1);
      tower.insert(candidate_position, local_gang(k));
    }
    tower.insert(0, table_id);
    towers.append(tower);
  }
  
  if (DEBUG>0) {
    cout<<"****************"<<endl;
    cout<<"TOWERS:"<<endl;
    FOR1D(towers, i) {
      FOR1D(towers(i), k) {
        cout<<" ";
        cout<<towers(i)(k);
        if (isBall(towers(i)(k), current_state)) {
          cout<<"o";
        }
      }
      cout << endl;
    }
    cout<<"****************"<<endl;
    
    cerr<<"****************"<<endl;
    cerr<<"TOWERS:"<<endl;
    FOR1D(towers, i) {
      FOR1D(towers(i), k) {
        cerr<<" ";
        cerr<<towers(i)(k);
        if (isBall(towers(i)(k), current_state)) {
          cerr<<"o";
        }
      }
      cerr << endl;
    }
    cerr<<"****************"<<endl;
  }
  

  //-----------------------------------------------------------------------------
  // (3) Translate to logic
  LitL lits;
  FOR1D(towers, i) {
    FOR1D(towers(i), k) {
      if (k==0)
        continue;
      lits.append(Literal::get(Symbol::get(MT::String("on")), TUP(towers(i)(k), towers(i)(k-1)), 1.));
    }
  }
  
  if (DEBUG>1) {cout<<"LOGIC:"<<endl;  write(lits); cout<<endl;}
  
  if (DEBUG>0) {cout<<"RobotManipulationSymbols::sampleGroundGoal__clearance [END]"<<endl;}
  return new LiteralListReward(lits);
}


void RobotManipulationSymbols::writeStateInfo(const SymbolicState& state, ostream& out) {
  out<<"--"<<endl;
  uint i, k;
  uint id_table = TL::UINT_NIL;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s->name == "table") {
      id_table = state.lits(i)->args(0);
      break;
    }
  }   
//   CHECK(id_table != TL::UINT_NIL, "");
  // Piles
  uintA piles;
  calcPiles(state, piles, 2);
//   PRINT(piles);
  
  for (i=0; i<piles.d0; i++) {
    k = 0;
    while (k < piles.d1 && piles(i,k) != TL::UINT_NIL) {
      if (k>0) out<<" ";
      if (isBox(piles(i,k), state))
        out << "b";
      out<<piles(i,k);
      if (isBall(piles(i,k), state))
        out<<"o";
      else if (isBox(piles(i,k), state)) {
        out<<"[ ";
        uint obj = getContainedObject(piles(i,k), state);
        if (obj != TL::UINT_NIL) {
          out<<obj;
          if (isBall(obj, state))
            out<<"o";
          out<<" ";
        }
        if (isClosed(piles(i,k), state))
          out<<"]";
      }
      k++;
    }
    out << endl;
  }
  
  out<<"--"<<endl;
  
  out<<"H ";
  uint id_inhand = getInhand(state);
  if (id_inhand != TL::UINT_NIL) {
    out<<id_inhand;
    if (isBall(id_inhand, state))
      out<<"o";
  }
  else
    out << "-";
  out<<endl;
  
  // Out objects
  uintA out_objects;
  SymbolicState::getArguments(out_objects, state, *Symbol::get(MT::String("out")));
  if (out_objects.N>0) {
    out<<"--"<<endl;
    out<<"OUT:  ";
    FOR1D(out_objects, i) {
      out<<out_objects(i);
      if (isBall(id_inhand, state))
        out<<"o";
      out<<" ";
    }
    out<<endl;
  }
    
  MT::Array< uintA > homieGangs;
  getHomieGangs(homieGangs, state);
  if (homieGangs.N > 0   &&   homieGangs(0).N > 1) {
    out<<"--"<<endl;
    out<<"Gangs:"<<endl;
    FOR1D(homieGangs, i) {
      out<<homieGangs(i)<<endl;
    }
  }
  
  out<<"--"<<endl;
  
  out<<"(existing objects: "<<reason::getConstants()<<")"<<endl;
}



}  // namespace PRADA
