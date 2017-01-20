/*  
    Copyright 2008-2012  Tobias Lang
    
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
    

#ifndef RELATIONAL_prada_h
#define RELATIONAL_prada_h

#include "plan.h"



namespace relational {


/************************************************
 * 
 *     PRADA_Planner
 * 
 ************************************************/

class PRADA_Reward;
class PRADA_DBN;

class PRADA_Planner : public NID_Planner {

  public:
    PRADA_Planner();
    ~PRADA_Planner();
    
    // *** central planning methods ***
    virtual Literal* plan_action(const SymbolicState& s, uint max_runs = 1);
    virtual void plan_full(LitL& plan, double& planValue, const SymbolicState& s, uint max_runs = 1);
    
    // setter methods
    virtual void setStartState(const SymbolicState& s0);
    void setNumberOfSamples(uint num_samples);
    virtual void setReward(Reward* reward);
    void setReward(Reward* reward, PRADA_Reward* prada_reward); // use only if you have domain specific knowledge!
    void setThresholdReward(double threshold_reward);
    static PRADA_Reward* create_PRADA_Reward(Reward* reward);
    void setNoiseSoftener(double noise_softener);
    void setRewardCalculation(bool reward_calculation__sum);
    void setActionChoice(bool action_choice__max);
    
    
    
    // -- All remaining methods not required for high-level planning from outside. -- 
    
    // samples "num_samples" plans and returns best
    void plan1_wrapper(LitL& best_plan, double& best_value, const SymbolicState& s, uint max_runs);
    bool plan1(LitL& best_plan, double& bestValue, uint num_samples);
    MT::Array< LitL > good_old_plans;
    
    // dynamic Bayesian network used for inference
    PRADA_DBN* dbn;
    SymL dbn_state_symbols;
    void calc_dbn_state_symbols();
    // DBN Construction
    void build_dbn(const uintA& constants, const SymL& preds, const SymL& actions);
    // Calculates from rules and rewards which predicates, functions, and actions to use
    void build_dbn(const uintA& constants);
    
    // Inference in DBN
    // calculate posterior over hidden state variables
    void infer(const LitL& plan); // one plan; infer until plan-length (no additional future action sampling)
    double inferStateRewards();
    double inferStateRewards(uint horizon);
    double inferStateRewards_limited_sum(uint horizon);
    double inferStateRewards_limited_max(uint horizon);
    double inferStateRewards_single(uint t);
    double calcRuleRewards(const LitL& actions);
    
    // writes variables in DBN
    void writeState(uint t, ostream& out = cout);
    void writeStateSparse(uint t, bool prim_only, ostream& out = cout);
    
    
  protected:
    // PRADA's parameters
    uint num_samples;
    double noise_softener;
    double threshold_reward;
    
    bool any_sensible_action;
    
    // true:   value = sum_t P(r_t)
    // false:  value = max_t P(r_t)
    bool reward_calculation__sum;
    
    // true:   first action of max action seq
    // false:  normalized sum over all action seqs starting with that action
    bool action_choice__max;
    
    // for statistics
    uintA action_choices;
    
    PRADA_Reward* prada_reward;
    static PRADA_Reward* convert_reward(LiteralReward* reward);
    static PRADA_Reward* convert_reward(MaximizeReward* reward);
    static PRADA_Reward* convert_reward(LiteralListReward* reward);
    static PRADA_Reward* convert_reward(DisjunctionReward* reward);
    static PRADA_Reward* convert_reward(NotTheseStatesReward* reward);
    static PRADA_Reward* convert_reward(RewardConjunction* reward);
    
    void setState(const SymbolicState& s, uint t);
    
    // Rather high-level methods
    void sampleActionsAndInfer(LitL& plan); // one plan
    // some actions are fixed
    void sampleActionsAndInfer(LitL& plan, const LitL& fixed_actions); // one plan
    // base function --> specify which net to sample on
    void sampleActionsAndInfer(LitL& plan, const LitL& fixed_actions, PRADA_DBN* net, uint local_horizon);
};




/************************************************
 * 
 *     A_PRADA
 * 
 ************************************************/

class A_PRADA : public PRADA_Planner {
  
  LitL last_seq;
  double last_value;
  
  // manipulations
  double shorten_plan(LitL& seq_best, const LitL& seq, double value_old);

  public:
    A_PRADA() {}
    void plan_full(LitL& plan, double& planValue, const SymbolicState& s, uint max_runs = 1);

    void reset();
};






/************************************************
 * 
 *     PRADA_DBN
 * 
 ************************************************/

class LiteralRV;
typedef MT::Array< LiteralRV* > RVL;

class PRADA_DBN {
  public:
    PRADA_DBN(const uintA& net_constants, const SymL& symbols, const SymL& actions, RuleSet& ground_rules, double noise_softener, uint horizon);
    ~PRADA_DBN();
    
    // Inference
    void inferRules(uint t);  // from time-slice at t
    void inferState(uint t, Literal* action);  // from action and time-slice at t-1
    void inferState(uint t, Literal* action, double given_action_weight);  // from action and time-slice at t-1
    void inferStates(const LitL& given_action_sequence); // for given_action_sequence
    
    // Set evidence
    void setAction(uint t, Literal* action);
    void setState(const LitL& lits, uint t);
    void setStateUniform(uint t);
    
    // Comparing probabilities
    double log_probability(uint t, const SymbolicState& state) const;
    double belief_difference(uint t, const arr& probs_p_prim, const arr& probs_f_prim) const;
    
    // Misc
    void calcDerived(uint t);
    void checkStateSoundness(uint t, bool omit_derived = false);
    void getBelief(uint t, arr& beliefs_p_prim, arr& beliefs_f_prim) const;
    
    // Writing
    void writeAllStates(bool prim_only = false, double threshold = 0.0, ostream& out = cout) const;
    void writeState(uint t, bool prim_only = false, double threshold = 0.0, ostream& out = cout) const;
    void writeStateSparse(uint t, bool prim_only, ostream& out = cout) const;
    void writeDAI(ostream& out = cout) const;
    void getActions(LitL& actions, uint horizon) const;
  
    
  public:
    
    void create_dbn_structure(const SymL& state_symbols, const SymL& actions);
    void create_dbn_params();
  
  RuleSet ground_rules;
  double noise_softener;
  uint horizon;
  uintA net_constants;  // ojects over which net is built; don't change them from outside!
  SymL net_symbols_state;
  SymL net_symbols_action;

  // Random variables
  uint start_id_rv__symbol_state;
  uint num_state_vars;
  RVL vars_state__prim;
  RVL vars_state__derived;
  RVL vars_action;
  RVL vars_context; // for faster code --> which vars are used as context
  
  arr vars_rules_simple;   // P(\phi_r | s)
  arr vars_rules; // P(\phi_r | -\phi_r', s)     2 dim: (1) timesteps,  (2) rules
  
  // Helping structures
  // Impacts only on primitives
  MT::Array< arr > impacts_V; // on vars as a whole;  2 dim: (1) var, (2) arr for rule x value
  MT::Array< arr > impacts_val; // on single values;  2 dim: (1) var, (2) arr for rule
  // dim_1 actions,  dim_2 arguments,  dim_3 rules
  uintA action2rules;
  uintA action2rules_num; // 2 dim; how many rules per constellation

  
  
  // --- Efficiency managing of random variables based on redundant data-structures ---
  // required for fast access by rewards and derive-methods
  void RVefficiency__init(const SymL& symbols);
  LiteralRV* RVefficiency__atom2var(Literal* lit) const;
  LiteralRV* RVefficiency__id2var(uint id_var) const;
  void RVefficiency__setAtom(Literal* lit, LiteralRV* var);
  SymL RVefficiency__symbols;
  RVL RVefficiency__LRV_A__structured;
  RVL RVefficiency__LRV_A__flat; // redundant container
  LitL RVefficiency__function_literals;
};



/************************************************
 * 
 *     PRADA_Reward
 * 
 ************************************************/

//   PRADA rewards  -->  working on beliefs
class PRADA_Reward {
  public :
    // Random variables --> Reals
    virtual double evaluate_prada_reward(const PRADA_DBN& net, uint t) = 0;
};
  


/************************************************
 * 
 *     LiteralRV
 * 
 ************************************************/

// -------------------------------------------------------------
// RANDOM VARIABLES


class LiteralRV {
  public:
    enum RVType {
      simple,
      expectation
    };
    
    Literal* lit;
    uint id;
    uint dim;
    uintA range;
    arr P; // over time;   dim-1: time,   dim-2: value
    bool changeable; // can change value over time
    RVType type;
    
    LiteralRV() {type = simple;}
    virtual ~LiteralRV() {}

    virtual void write(ostream& os = cout);
};

// Calculates only Expected Value of Random Variable.
// No true distribution is maintained!!!
class ExpectationRV : public LiteralRV {
  // in array P, expectation values are stored
  public:
    ExpectationRV() {type = expectation; dim=1;}
    void write(ostream& os = cout);
};



void write(const RVL& vars);



}



#endif // RELATIONAL_prada_h

