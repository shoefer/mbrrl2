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
    

#ifndef RELATIONAL_plan_h
#define RELATIONAL_plan_h

#include "reason.h"


namespace relational {
  

/************************************************
 * 
 *     NID-Planner
 * 
 ************************************************/

// WorldAbstraction provides basic knowledge about the action space and sampling state transitions.
class WorldAbstraction {
  public:
    LitL ground_actions;
    // returns TL::TL_DOUBLE_NIL if action not applicable
    virtual double sampleSuccessorState(SymbolicState& s_suc, uint& flag, const SymbolicState& s_prev, Literal* action) const  = 0;
    virtual double postprocessValue(double value, uint flag) const {return value;}
};


class Reward;

class NID_Planner : public WorldAbstraction {
  public:
    NID_Planner();
    virtual ~NID_Planner();
    
    // ***** central planning method: state --> action *****
    virtual Literal* plan_action(const SymbolicState& current_state, uint max_runs = 1) = 0;
    
    // setter methods
    void setGroundRules(RuleSet& ground_rules);
    void setDiscount(double discount);
    virtual void setHorizon(uint horizon); // planning horizon
    virtual void setReward(Reward* reward);
    void setNoiseScalingFactor(double noise_scaling_factor);
    Reward* getReward() {return reward;}
    
    // returns TL::TL_DOUBLE_NIL if action not applicable
    double sampleSuccessorState(SymbolicState& s_suc, uint& flag, const SymbolicState& s_prev, Literal* action) const;
    double postprocessValue(double value, uint flag) const;
    
    
  #define DEFAULT__NID_PLANNER__HORIZON 1
  #define DEFAULT__NID_PLANNER__DISCOUNT 0.95
  #define DEFAULT__NID_PLANNER__USE_RULE_OUTCOME_REWARDS false
  #define DEFAULT__NID_PLANNER__NOISE_SCALING_FACTOR 0.1
  protected:
    double discount;
    uint horizon;
    Reward* reward;
    double noise_scaling_factor;
    RuleSet ground_rules;
    boolA is_manipulating_rule;  // predicts non-noise changes?
    bool use_ruleOutcome_rewards;
    arr expected_rule_rewards;
    
    arr discount_pow; // for faster code
};




/************************************************
 * 
 *     NID-SST
 * 
 ************************************************/

namespace SST {
  Literal* plan_action(double& value, const SymbolicState& current_state, const Reward& reward, uint branch, uint T, double discount, const WorldAbstraction& wa);
}


class NID_SST : public NID_Planner {
  uint branch;
  public :
    NID_SST(uint branch);
    Literal* plan_action(const SymbolicState& current_state, uint max_runs = 1);
};




/************************************************
 * 
 *     NID-UCT
 * 
 ************************************************/

struct StateActionValues {
  public:
    const SymbolicState& s;
    arr values;
    uintA visits;

    StateActionValues(const SymbolicState& s, uint num_actions);
    ~StateActionValues();
    
    // counts: c(s,a)
    uint getVisits();
    uint getVisits(uint action_id);
    void increaseVisits(uint action_id);
    
    // Q-values: Q(s,a)
    double getQvalue(uint action_id);
    void setQvalue(uint action_id, double value);
};


class NID_UCT : public NID_Planner {
public:
  NID_UCT();
  ~NID_UCT();
  
  Literal* plan_action(const SymbolicState& current_state, uint max_runs = 1);
  
  void setC(double c);
  void setNumEpisodes(uint numEpisodes);
  
private:
  #define DEFAULT__NID_UCT__C 1.0
  #define DEFAULT__NID_UCT__NUM_EPISODES 500
  double c;
  uint numEpisodes;
  MT::Array< StateActionValues* > s_a_values;
  StateActionValues* getStateActionValues(const SymbolicState& s);
  void killStateActionValues();
  void runEpisode(double& reward, const SymbolicState& s, uint t);
};


  

/************************************************
 * 
 *     Reward
 * 
 ************************************************/

class Reward {  
  public:
    enum RewardType { 
                      reward_literal,  //  literal such as on(a,b)=true
                      reward_literalList,  // list ofliterals
                      reward_maximize_function,  // try to make function as big as possible
                                                 // --> Advantage: can work on expectations in case of beliefs over states
                      reward_not_these_states, 
                      reward_one_of_literal_list,
                      conjunction_of_rewards
    };
    
    Reward();
    Reward(RewardType reward_type);
    virtual ~Reward() {}
    
    virtual double evaluate(const SymbolicState& s) const = 0;
    virtual bool satisfied(const SymbolicState& s) const = 0;
    virtual bool possible(const SymbolicState& s) const = 0;
    
    virtual void getRewardConstants(uintA& constants, const SymbolicState* s) const = 0;
    
    virtual void write(ostream& out = cout) const {};
    virtual void write(const char* filename) const = 0;
    
    RewardType reward_type;
    
    static Reward* read(const char* filename);
};


class LiteralReward : public Reward {
  public:
    Literal* lit;
    
    LiteralReward(Literal* lit);
    
    double evaluate(const SymbolicState& s) const ;
    bool satisfied(const SymbolicState& s) const;
    bool possible(const SymbolicState& s) const;
    
    void getRewardConstants(uintA& constants, const SymbolicState* s = NULL) const;
    
    void write(ostream& out = cout) const;
    void write(const char* filename) const;
};


class LiteralListReward : public Reward {
  public:
    LitL lits;
    
    LiteralListReward(LitL& lits);
    
    double evaluate(const SymbolicState& s) const ;
    bool satisfied(const SymbolicState& s) const;
    bool possible(const SymbolicState& s) const;
    
    void getRewardConstants(uintA& constants, const SymbolicState* s) const;
    
    void write(ostream& out = cout) const;
    void write(const char* filename) const;
};


class MaximizeReward : public Reward {
  public:
    Literal* literal_to_be_maximized;
    double offset;
    
    MaximizeReward(); // --> maximize function
    MaximizeReward(Literal* function_literal); // --> maximize function
    
    double evaluate(const SymbolicState& s) const ;
    bool satisfied(const SymbolicState& s) const;
    bool possible(const SymbolicState& s) const;
    
    void getRewardConstants(uintA& constants, const SymbolicState* s) const;
    
    void write(ostream& out = cout) const;
    void write(const char* filename) const;
};


class DisjunctionReward : public Reward {
  public:
    LitL lits;
    arr weights;

    DisjunctionReward(LitL& lits);
    DisjunctionReward(LitL& lits, arr& weights);
    
    double evaluate(const SymbolicState& s) const;
    bool satisfied(const SymbolicState& s) const;
    bool possible(const SymbolicState& s) const;
    
    void getRewardConstants(uintA& constants, const SymbolicState* s) const;
    
    void write(ostream& out = cout) const;
    void write(const char* filename) const;
};


class NotTheseStatesReward : public Reward {
  public:
    SymbolicStateL undesired_states;
    
    NotTheseStatesReward(const SymbolicStateL& undesired_states); // --> maximize function
    ~NotTheseStatesReward() {}
    
    double evaluate(const SymbolicState& s) const ;
    bool satisfied(const SymbolicState& s) const;
    bool possible(const SymbolicState& s) const;
    
    void getRewardConstants(uintA& constants, const SymbolicState* s) const;
    
    void write(ostream& out = cout) const;
    void write(const char* filename) const;
};


class RewardConjunction : public Reward {
  public:
    MT::Array<Reward*> rewards;

    RewardConjunction() : Reward(conjunction_of_rewards) {}
    virtual ~RewardConjunction() {}

    void addReward(Reward *reward);

    double evaluate(const SymbolicState& s) const ;
    bool satisfied(const SymbolicState& s) const;
    bool possible(const SymbolicState& s) const;

    void getRewardConstants(uintA& constants, const SymbolicState* s = NULL) const;

    void write(ostream& out = cout) const;
    void write(const char* filename) const;
};


}

std::ostream& operator<<(std::ostream& os, const relational::Reward& p);

#endif // RELATIONAL_plan_h
