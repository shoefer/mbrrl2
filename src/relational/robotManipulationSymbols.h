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

#ifndef RELATIONAL_robotManipulationSymbols_h
#define RELATIONAL_robotManipulationSymbols_h

#include "relational/reason.h"
#include "relational/plan.h"


namespace relational {
  
namespace RobotManipulationSymbols {
  

  /************************************************
  * 
  *     Set up the symbols
  * 
  ************************************************/
  
  void initSymbols();
  
  
  // actions
  Symbol* getSymbol_action_grab();
  Symbol* getSymbol_action_puton();
  Symbol* getSymbol_action_lift(); // 2 args
  Symbol* getSymbol_action_place(); // 2 args
  Symbol* getSymbol_action_openBox();
  Symbol* getSymbol_action_closeBox();
  
  // primitives
  Symbol* getSymbol_table();
  Symbol* getSymbol_block();
  Symbol* getSymbol_ball();
  Symbol* getSymbol_box();
  Symbol* getSymbol_on();
  Symbol* getSymbol_inhand();
  Symbol* getSymbol_upright();
  Symbol* getSymbol_out();
  Symbol* getSymbol_homies();
  Symbol* getSymbol_contains();
  Symbol* getSymbol_closed();
  Symbol* getFunction_size();
  
  // derived
  // --> look here to see how to define derived symbols
  ConjunctionSymbol* getSymbol_clear();
  TransClosureSymbol* getSymbol_above();
  ConjunctionSymbol* getSymbol_aboveNotable();
  ConjunctionSymbol* getSymbol_dirtyGuyBelow();
  ConjunctionSymbol* getSymbol_diffTower();
  ConjunctionSymbol* getSymbol_withoutHomies();
  ConjunctionSymbol* getSymbol_inorder();
  ConjunctionSymbol* getSymbol_inhandNil();
  ConjunctionSymbol* getSymbol_onBox();
  CountSymbol* getSymbol_height();
  AverageFunction* getFunction_avgheight();
  SumFunction* getFunction_sumheight();
  CountSymbol* getFunction_countInorder();
  

  
  
  
  /************************************************
  * 
  *     Reward
  * 
  ************************************************/
  
  namespace RewardLibrary {
    Reward* stack();
    Reward* on(uint o1, uint o2);
    Reward* inhand(uint o1);
    Reward* tower(uintA& tower_objects);
    Reward* clearance();
  };
  
  
  /************************************************
   * 
   *     Helpers for SymbolicState
   * 
   ************************************************/
  
  bool isBlock(uint id, const SymbolicState& s);
  bool isBall(uint id, const SymbolicState& s);
  bool isBox(uint id, const SymbolicState& s);
  bool isTable(uint id, const SymbolicState& s);
  bool isOut(uint id, const SymbolicState& s);
  bool isInhand(uint id, const SymbolicState& s);
  bool isClosed(uint box_id, const SymbolicState& s);
  bool isInorderGang(const uintA gang, const SymbolicState& s);
  void getBelowObjects(uintA& ids, uint id, const SymbolicState& s);
  void getAboveObjects(uintA& ids, uint id, const SymbolicState& s); // directly above!!!
  uint getBelow(uint id, const SymbolicState& s);
  uint getAbove(uint id, const SymbolicState& s); // directly above!!!
  uint getInhand(const SymbolicState& s);
  void getBoxes(uintA& ids, const SymbolicState& s);
  uint getContainingBox(uint obj_id, const SymbolicState& s);
  uint getContainedObject(uint box_id, const SymbolicState& s);
  void getHomieGangs(MT::Array< uintA >& homieGangs, const SymbolicState& s);
  void getOutObjects(uintA& outs, const SymbolicState& s);
  
  // Reward functions
  double reward_buildTower(const SymbolicState& s);
  // Piles;  sorted by heights! piles(0)=highest
  void calcPiles(const SymbolicState& s, uintA& piles, uint sort_type = 1);
  void calcHeights(const uintA& objects, const uintA& piles, uintA& object_heights, uint id_table);
  bool has_maximum_stack_value(const SymbolicState& s);
  void calcSkyscraperWeights(const uintA& objects, const uintA& piles, double skyscraper_bias, arr& weights, bool highGood, uint id_table);
  void calcSkyscraperWeights(const uintA& heights, double skyscraper_bias, arr& weights, bool highGood, uint id_table);

  // Goal state sampling
  LiteralListReward* sampleGroundGoal__stack(const uintA& blocks, const uintA& balls, uint table_id, bool leave_existing_towers = false, SymbolicState* = NULL);
  LiteralListReward* sampleGroundGoal__clearance(const SymbolicState& current_state, uint table_id);
  
  // Nice state information
  void writeStateInfo(const SymbolicState& s, ostream& out = std::cout);
}

}

#endif // RELATIONAL_robotManipulationSymbols_h
