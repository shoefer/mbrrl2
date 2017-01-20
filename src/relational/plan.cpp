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


#include "plan.h"
#include "reason.h"


namespace relational {

  
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
//    NID Planner
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


NID_Planner::NID_Planner() {
  this->horizon = DEFAULT__NID_PLANNER__HORIZON;
  this->discount = DEFAULT__NID_PLANNER__DISCOUNT;
  this->use_ruleOutcome_rewards = DEFAULT__NID_PLANNER__USE_RULE_OUTCOME_REWARDS;
  this->noise_scaling_factor = DEFAULT__NID_PLANNER__NOISE_SCALING_FACTOR;
}


NID_Planner::~NID_Planner() {
}


void NID_Planner::setGroundRules(RuleSet& ground_rules) {
  this->ground_rules = ground_rules;
  // fill actions list
  ground_actions.clear();
  is_manipulating_rule.clear();
  Literal* last_action = NULL;
  uint i;
  FOR1D_(ground_rules, i) {
    if (last_action != ground_rules.elem(i)->action  &&
        ground_rules.elem(i)->action != Literal::getLiteral_default_action()) {
      last_action = ground_rules.elem(i)->action;
      ground_actions.setAppend(last_action);
    }
    bool manipulates = false;
    uint o;
    for (o=0; o<ground_rules.elem(i)->outcomes.N-1; o++) {
      if (ground_rules.elem(i)->outcomes(o).N > 0) {
        manipulates = true;
        break;
      }
    }
    is_manipulating_rule.append(manipulates);
  }
  CHECK(is_manipulating_rule.N == this->ground_rules.num(), "");
  
  FOR1D_(ground_rules, i) {
    if (ground_rules.elem(i)->outcome_rewards.N > 0) {
      this->use_ruleOutcome_rewards = true;
      break;
    }
  }
  
  if (this->use_ruleOutcome_rewards) {
    expected_rule_rewards.resize(ground_rules.num());
    uint k;
    FOR1D_(ground_rules, i) {
      Rule* rule = ground_rules.elem(i);
      if (rule->outcome_rewards.N > 0) {
        double expected_reward = 0.;
//         rule->writeNice();
        FOR1D(rule->outcomes, k) {
          expected_reward += rule->probs(k) * rule->outcome_rewards(k);
        }
        expected_rule_rewards(i) = expected_reward;
      }
    }
  }
}


void NID_Planner::setDiscount(double discount) {
  this->discount = discount;
  discount_pow.resize(horizon+1);
  uint i=0;
  for (i=0;i<=horizon; i++) {
    discount_pow(i) = pow(discount, i);
  }
}


void NID_Planner::setHorizon(uint horizon) {
  this->horizon = horizon;
  if (discount_pow.N < horizon+1) {
    discount_pow.resize(horizon+1);
    uint i=0;
    for (i=0;i<=horizon; i++) {
      discount_pow(i) = pow(discount, i);
    }
  }
}


void NID_Planner::setReward(Reward* reward) {
  this->reward = reward;
}


void NID_Planner::setNoiseScalingFactor(double noise_scaling_factor) {
  this->noise_scaling_factor = noise_scaling_factor;
}


double NID_Planner::sampleSuccessorState(SymbolicState& s_suc, uint& flag, const SymbolicState& s_prev, Literal* action) const {
  return reason::sampleSuccessorState_groundRules(s_suc, s_prev, ground_rules, action, flag);
}


double NID_Planner::postprocessValue(double value, uint flag) const {
  double value_processed;
  value_processed = value;
  if (flag == STATE_TRANSITION__NOISE_OUTCOME) {
    value_processed *= noise_scaling_factor;
  }
  return value_processed;
}

  


/************************************************
 * 
 *     NID-SST
 * 
 ************************************************/

LitL DEBUG_actions;

Literal* SST::plan_action(double& value, const SymbolicState& s0, const Reward& reward, uint branch, uint tau, double discount, const WorldAbstraction& wa) {
  uint DEBUG = 0;
  if (DEBUG_actions.N == 0)
    DEBUG_actions.resize(tau);
  if (DEBUG>0) {
    cout<<"+ SST - start  tau="<<tau<<endl;
    s0.write(); cout<<endl;
  }
  if (tau == 5)
    cout<<"."<<std::flush;
  double reward_s0 = reward.evaluate(s0);
  double reward_tree = 0.;
  uint i, b;
  Literal* action = NULL;
  bool action_is_applicable;
  if (tau>0) {
    arr action_values(wa.ground_actions.N);
    FOR1D(wa.ground_actions, i) {
      double action_value = 0.;
      action_is_applicable = true;
      if (DEBUG>0) {DEBUG_actions(DEBUG_actions.N-tau) = wa.ground_actions(i);}
      for(b=0; b<branch; b++) {
        uint flag;
        SymbolicState s_suc;
        double ruleOutcome_reward = wa.sampleSuccessorState(s_suc, flag, s0, wa.ground_actions(i));
        if (DEBUG>1) {PRINT(ruleOutcome_reward);}
        action_is_applicable = (ruleOutcome_reward != TL::TL_DOUBLE_NIL);
        if (!action_is_applicable) {
          if (DEBUG>1) {
            cout<<"++ Omitting at "<<tau<<" "; wa.ground_actions(i)->write(); cout<<" b="<<b<<endl;
          }
          action_value = -10000.;
          break;
        }
        else {
          double tree_value;
          if (DEBUG>1) {
            cout<<"++ going down tau="<<tau<<" "; wa.ground_actions(i)->write(); cout<<" b="<<b<<endl;
          }
          plan_action(tree_value, s_suc, reward, branch, tau-1, discount, wa);
          tree_value = wa.postprocessValue(tree_value, flag);
          tree_value += ruleOutcome_reward;
          action_value += tree_value;
        }
      }
      action_value /= branch;
      action_values(i) = action_value;
    }
    uint max_id = action_values.maxIndex();
    action = wa.ground_actions(max_id);
    reward_tree = action_values(max_id);
    if (DEBUG>1) {
      cout<<"Best Action for "; uint t; for (t=0; t<DEBUG_actions.N-tau; t++){DEBUG_actions(t)->write(); cout<<" ";}
      cout<<"  -->  ";action->write(); cout<<endl;
    }
  }
  value = reward_s0 + discount * reward_tree;
  if (DEBUG>0) {
    PRINT(reward_s0);
    PRINT(reward_tree);
    PRINT(value);
    cout<<"+ SST - end  tau="<<tau<<endl;
  }
  return action;
}


NID_SST::NID_SST(uint branch) : NID_Planner() {
  this->branch = branch;
}


Literal* NID_SST::plan_action(const SymbolicState& current_state, uint max_runs) {
  double value;
  uint i;
  for (i=0; i<max_runs; i++) {
     Literal* action = SST::plan_action(value, current_state, *reward, branch, horizon, discount, *this);
#define SST_THRESHOLD 0.005
     if (value > SST_THRESHOLD) {
      return action;
    }
  }
  MT_MSG("NID_SST: No good action found.");
  return NULL;
}






/************************************************
 * 
 *     NID-UCT
 * 
 ************************************************/

StateActionValues::StateActionValues(const SymbolicState& _s, uint num_actions) : s(_s) {
  values.resize(num_actions);
  values.setUni(0.);
  visits.resize(num_actions);
  visits.setUni(0.);
}


StateActionValues::~StateActionValues() {
}


uint StateActionValues::getVisits() {
  return sum(visits);
}


uint StateActionValues::getVisits(uint action_id) {
  return visits(action_id);
}


void StateActionValues::increaseVisits(uint action_id) {
  visits(action_id)++;
}


double StateActionValues::getQvalue(uint action_id) {
  return values(action_id);
}


void StateActionValues::setQvalue(uint action_id, double value) {
  values(action_id) = value;
}


NID_UCT::NID_UCT() : NID_Planner() {
  c = DEFAULT__NID_UCT__C;
  numEpisodes = DEFAULT__NID_UCT__NUM_EPISODES;
}


NID_UCT::~NID_UCT() {
  killStateActionValues();
}


LitL DEBUG__UCT_ACTIONS;
Literal* NID_UCT::plan_action(const SymbolicState& s, uint max_runs) {
  uint DEBUG = 0;
  killStateActionValues(); // full replanning... comment if not desired
  uint i, k;
  DEBUG__UCT_ACTIONS.resize(horizon);
  for (k=0; k<max_runs; k++) {
    for (i=0; i<numEpisodes; i++) {
      DEBUG__UCT_ACTIONS.setUni(NULL);
      if (i%10 == 0) cerr<<"."<<std::flush;
      double dummy_reward;
      runEpisode(dummy_reward, s, 0);
      if (DEBUG>1) {cout<<i<<":  "; write(DEBUG__UCT_ACTIONS); cout<<endl;}
    }
    // get maximum q value for s
    StateActionValues* s_a_info = getStateActionValues(s);
    if (DEBUG>1) {
      cout<<"Q values for starting state for states tried more than 0 times:"<<endl;
      FOR1D(ground_actions, i) {
        if (s_a_info->getVisits(i) > 0) {
          ground_actions(i)->write(); cout<<": "<<s_a_info->getQvalue(i)<<"   ("<<s_a_info->getVisits(i)<<" visits)"<<endl;
        }
      }
    }
    uint max_id = s_a_info->values.maxIndex();
#define THRESHOLD_UCT 0.0005
    if (s_a_info->values(max_id) > THRESHOLD_UCT)
      return ground_actions(max_id);
    MT_MSG("NID_UCT: No good action found --> Retry!");
  }
  MT_MSG("NID_UCT: Still no good action found. I'll give up :-(.");
  return NULL;
}


void NID_UCT::setC(double c) {
  this->c = c;
}


void NID_UCT::setNumEpisodes(uint numEpisodes) {
  this->numEpisodes = numEpisodes;
}


StateActionValues* NID_UCT::getStateActionValues(const SymbolicState& state) {
  uint i;
  FOR1D(s_a_values, i) {
    if (s_a_values(i)->s == state)
      return s_a_values(i);
  }
  // create new one
  StateActionValues* s_a_info = new StateActionValues(state, ground_actions.N);
  s_a_values.append(s_a_info);
  return s_a_info;
}


void NID_UCT::killStateActionValues() {
  listDelete(s_a_values);
}


void NID_UCT::runEpisode(double& value, const SymbolicState& s, uint t) {
  uint DEBUG = 0;
  if (t==0)
    DEBUG = 0;
  else
    DEBUG = 0;
  if (DEBUG>0) {
    cout<<"NID_UCT::runEpisode() [START]"<<endl;
    PRINT(t);
    cout<<"SymbolicState: "; s.write(cout, true); cout<<endl;
  }
  double reward_s = reward->evaluate(s);
//   if (reward_s > 0.) {
//     cout<<"DERBE GESCHICHTE!!!!!"<<endl;
//     writeNice(DEBUG__UCT_ACTIONS); cout<<endl;
//     PRINT(t);
//   }
  if (DEBUG>0) {PRINT(reward_s);}
  if (t==horizon) {
    value = reward_s;
  }
  else if (reward->satisfied(s)) {
    value = reward_s;
//     cout<<"finished:"<<endl;
//     cout<<"SymbolicState: "; s.writeNice(); cout<<endl;
//     uint t2;
//     cout<<"Actions:  ";
//     for (t2=0; t2<t; t2++) {
//       DEBUG__UCT_ACTIONS(t2)->writeNice(cout); cout<<" ";
//     }
//     cout<<endl;
  }
  else {
    StateActionValues* s_a_info = getStateActionValues(s);
    uint i;
    arr UCB(ground_actions.N);
    RuleSet rules;
    uintA untried_action_ids;
    if (DEBUG>0) {cout<<"visits(s)="<<s_a_info->getVisits()<<endl;}
    FOR1D(ground_actions, i) {
      SymbolicState dummy_state = s;
      Rule* r = ground_rules.elem(reason::calc_uniqueCoveringRule_groundRules_groundAction(ground_rules, s, ground_actions(i)));
      if (!Rule::isDefaultRule(r)) {
        rules.append(r);
        if (s_a_info->getVisits(i) == 0) {
          untried_action_ids.append(i);
          UCB(i) = -11.;
        }
        else
          UCB(i) = s_a_info->getQvalue(i)  +  c * sqrt(log((double)s_a_info->getVisits()) / (1.0 * s_a_info->getVisits(i)));
      }
      else {
        rules.append(Rule::getDoNothingRule()); // just as a hack
        UCB(i) = -22.;
      }
      if (DEBUG>1) {
        ground_actions(i)->write(); cout<<" UCB="<<UCB(i)<<"   (q="<<s_a_info->getQvalue(i)<<",  visits="<<s_a_info->getVisits(i)<<")"<<endl;
      }
    }
    uint id_opt;
    if (untried_action_ids.N > 0)
      id_opt = untried_action_ids(rnd.num(untried_action_ids.N));
    else
      id_opt = UCB.maxIndex();
    if (DEBUG>0) {
      cout<<" --> Chosen action: "; ground_actions(id_opt)->write(); cout<<endl;
    }
    DEBUG__UCT_ACTIONS(t) = ground_actions(id_opt);
    uint flag;
    SymbolicState s_suc;
    double ruleOutcome_value = reason::sampleSuccessorState_groundRule(s_suc, s, rules.elem(id_opt), flag);
  
    double reward_suc = 0.;
    runEpisode(reward_suc, s_suc, t+1);  // recursive call
    reward_suc = postprocessValue(reward_suc, flag);
    value = reward_s +  discount * reward_suc;
    if (use_ruleOutcome_rewards) {
      value += ruleOutcome_value;
    }
    // update Q-value
    s_a_info->increaseVisits(id_opt);
    double old_q = s_a_info->getQvalue(id_opt);
    double new_q = old_q + (1 / (1.0 * s_a_info->getVisits(id_opt))) * (value - old_q);
    s_a_info->setQvalue(id_opt, new_q);
    if (DEBUG>0) {
      PRINT(reward_suc);
      PRINT(ruleOutcome_value);
      PRINT(value);
      PRINT(old_q);
      PRINT(new_q);
    }
  }
  
  if (DEBUG>0) {
    cout<<"NID_UCT::runEpisode() [END]"<<endl;
  }
}










/************************************************
 * 
 *     Reward
 * 
 ************************************************/

Reward::Reward() {
}


Reward::Reward(RewardType _reward_type) {
  reward_type = _reward_type;
}


Reward* Reward::read(const char* filename) {
  ifstream in(filename);
  if (!in.is_open()) HALT("File cannot be opened.");
  MT::skip(in);
  uint type;
  in >> type;
  MT::skip(in);
  if (type == Reward::reward_literal) {
    MT::String line;
    line.read(in, NULL, "\n");
    Literal* lit = Literal::get(line);
    uint i;
    FOR1D(lit->args, i) {
      if (reason::getConstants().findValue(lit->args(i)) < 0)
        HALT("Reward uses unknown argument "<<lit->args(i));
    }
    return new LiteralReward(lit);
  }
  else if (type == Reward::reward_literalList) {
    LitL lits;
    while (MT::skip(in) != -1) {
      MT::String line;
      line.read(in, NULL, "\n");
      Literal* lit = Literal::get(line);
      uint i;
      FOR1D(lit->args, i) {
        if (reason::getConstants().findValue(lit->args(i)) < 0)
          HALT("Reward uses unknown argument "<<lit->args(i));
      }
      lits.append(lit);
    }
    return new LiteralListReward(lits);
  }
  else if (type == Reward::reward_maximize_function) {
    MT::String line;
    line.read(in, NULL, "\n");
    NIY;
//     FunctionAtom* fa = logicObjectManager::getFA(line);
//     uint i;
//     FOR1D(fa->args, i) {
//       if (logicObjectManager::constants.findValue(fa->args(i)) < 0)
//         HALT("Reward uses unknown argument "<<fa->args(i));
//     }
//     fa = logicObjectManager::getFAorig(fa);
//     return new MaximizeReward(fa);
  }
  else
    HALT("Unknown reward type " << type << " in file "<< filename);
  return NULL;
}


/************************************************
 * 
 *     LiteralReward
 * 
 ************************************************/

LiteralReward::LiteralReward(Literal* _lit) : Reward(reward_literal) {
  //CHECK(_lit->value, "");
  lit = _lit;
}


double LiteralReward::evaluate(const SymbolicState& state) const {
  if (reason::holds(state.lits, lit))
    return 1.0;
  else
    return 0.0;
}


bool LiteralReward::satisfied(const SymbolicState& state) const {
  double result = evaluate(state);
  return TL::areEqual(1., result);
}


bool LiteralReward::possible(const SymbolicState& state) const {
  // special implementation for robot manipulation domain
  Symbol* p_OUT = Symbol::get(MT::String("out"));
  if (p_OUT == NULL) return true;
  uint i;
  FOR1D(lit->args, i) {
    uintA args(1);  args(0)=lit->args(i);
    Literal* pi_out = Literal::get(p_OUT, args, 1.);
    if (reason::holds(state.lits, pi_out))
      return false;
  }
  return true;
}


void LiteralReward::getRewardConstants(uintA& constants, const SymbolicState* s) const {
  constants.clear();
  constants.setAppend(lit->args);
}


void LiteralReward::write(ostream& out) const {
  lit->write(out);
}


void LiteralReward::write(const char* filename) const {
  ofstream out(filename);
  out<<reward_type<<endl;
  out<<"# LiteralReward"<<endl;
  out<<"# "; lit->write(out); out<<endl;
  lit->write(out); out<<endl;
  out.close();
}



/************************************************
 * 
 *     LiteralListReward
 * 
 ************************************************/

LiteralListReward::LiteralListReward(LitL& _lits) : Reward(reward_literalList) {
  lits = _lits;
}


double LiteralListReward::evaluate(const SymbolicState& state) const {
  if (reason::holds(state.lits, lits))
    return 1.0;
  else
    return 0.0;
}


bool LiteralListReward::satisfied(const SymbolicState& state) const {
  double result = evaluate(state);
  if (TL::areEqual(1., result))
    return true;
  else
    return false;
}


bool LiteralListReward::possible(const SymbolicState& state) const {
  // BRING IN DOMAIN KNOWLEDGE
  
  // Desktop world domain
  Symbol* p_OUT = Symbol::get(MT::String("out"));
  uint i, k;
  if (p_OUT != NULL) {
    FOR1D(lits, k) {
      FOR1D(lits(k)->args, i) {
        uintA args(1);  args(0)=lits(k)->args(i);
        Literal* pt_out = Literal::get(p_OUT, args, 1.);
        if (reason::holds(state.lits, pt_out))
          return false;
      }
    }
  }
  
  // Ex-Blocksworld domain
  Symbol* p_NO_DESTROYED_TABLE = Symbol::get(MT::String("no-destroyed-table"));
  Symbol* p_ON_TABLE = Symbol::get(MT::String("on-table"));
  Symbol* p_NO_DESTROYED = Symbol::get(MT::String("no-destroyed"));
  if (p_NO_DESTROYED_TABLE != NULL  &&  p_ON_TABLE != NULL) {
    uintA empty;
    Literal* pi_no_destroyed_table = Literal::get(p_NO_DESTROYED_TABLE, empty, 1.);
    FOR1D(lits, k) {
      // (1)  Impossible if one object X still has to be put "on-table(X)" (but is not yet!)
      // and it already does NOT hold "no-destroyed-table()".
      if (lits(k)->s == p_ON_TABLE  &&  lits(k)->value) {
//         PRINT(*lits(k));
//         PRINT(!reason::holds(s, lits(k)));
//         PRINT(!reason::holds(s, pi_no_destroyed_table));
        if (!reason::holds(state.lits, lits(k)) &&  !reason::holds(state.lits, pi_no_destroyed_table)) {
          cout<<"Impossible as "<<*pi_no_destroyed_table<<" does not hold and we still require "<<*lits(k)<<endl;
          cerr<<"Impossible as "<<*pi_no_destroyed_table<<" does not hold and we still require "<<*lits(k)<<endl;
          return false;
        }
      }
      // (2)  Impossible if object to be moved is destroyed.
      if (lits(k)->args.N == 2) {
        if (!reason::holds(state.lits, lits(k))) {
          uintA helper;  helper.append(lits(k)->args(0));
          Literal* pi_helper = Literal::get(p_NO_DESTROYED, helper, 1.);
          if (!reason::holds(state.lits, pi_helper)) {
            cout<<"Impossible as "<<lits(k)->args(0)<<" is already destroyed; i.e., it does not hold that "<<*pi_helper<<endl;
            cerr<<"Impossible as "<<lits(k)->args(0)<<" is already destroyed; i.e., it does not hold that "<<*pi_helper<<endl;
            return false;
          }
        }
      }
    }
  }
  
  // Triangle-tireworld domain
  Symbol* p_VEHICLE_AT = Symbol::get(MT::String("vehicle-at"));
  Symbol* p_SPARE_IN = Symbol::get(MT::String("spare-in"));
  Symbol* p_NOT_FLATTIRE = Symbol::get(MT::String("not-flattire"));
  Symbol* p_HASSPARE = Symbol::get(MT::String("hasspare"));
  if (p_VEHICLE_AT != NULL  &&  p_NOT_FLATTIRE != NULL) {
    uintA empty;
    Literal* pi_not_flattire = Literal::get(p_NOT_FLATTIRE, empty, 1.);
    Literal* pi_hasspare = Literal::get(p_HASSPARE, empty, 1.);
    FOR1D(lits, k) {
      // (1)  Impossible if   (i) not in goal-location,  (ii) not not-flattire,  (iii) not has-spare   and (iv) not spare-in current location
      if (lits(k)->s == p_VEHICLE_AT  &&  lits(k)->value) {
        // (i)
        if (!reason::holds(state.lits, lits(k))) {
          // (ii), (iii)
          if (!reason::holds(state.lits, pi_not_flattire)  &&  !reason::holds(state.lits, pi_hasspare) ) {
            // (iv)
            uint current_location = SymbolicState::getArgument(state, *p_VEHICLE_AT);
            uintA wrapper;  wrapper.append(current_location);
            Literal* pi_spare_in_current_location = Literal::get(p_SPARE_IN, wrapper, 1.);
            if (!reason::holds(state.lits, pi_spare_in_current_location)) {
              cout<<"Impossible as -" << *lits(k) << " (REWARD), but -"<<*pi_not_flattire<<",  current_location="
                        <<current_location<<", -"<<*pi_hasspare<<" and -"
                        <<*pi_spare_in_current_location<< " (CURRENT STATE)."<<endl;
              return false;
            }
          }
        }
      }
    }
  }
  return true;
}


void LiteralListReward::getRewardConstants(uintA& constants, const SymbolicState* s) const {
  constants.clear();
  uint i;
  FOR1D(lits, i) {
    constants.setAppend(lits(i)->args);
  }
}


void LiteralListReward::write(ostream& out) const {
  uint i;
  out<<"[" << lits.N  << "]  ";
  FOR1D(lits, i) {
    lits(i)->write(out);out<<" ";
  }
}


void LiteralListReward::write(const char* filename) const {
  ofstream out(filename);
  out<<reward_type<<endl;
  out<<"# LiteralListReward"<<endl;
  uint i;
  FOR1D(lits, i) {
    out<<"# "; lits(i)->write(out); out<<endl;
  }
  out<<lits.N<<endl;
  FOR1D(lits, i) {
    lits(i)->write(out); out<<endl;
  }
  out.close();
}



/************************************************
 * 
 *     MaximizeReward
 * 
 ************************************************/

MaximizeReward::MaximizeReward() : Reward(reward_maximize_function) {
  literal_to_be_maximized = NULL;
  offset = 0;
}

MaximizeReward::MaximizeReward(Literal* _function_literal) : Reward(reward_maximize_function) {
  offset = 0;
  literal_to_be_maximized = _function_literal;
}

double MaximizeReward::evaluate(const SymbolicState& state) const {
  uint i;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s == literal_to_be_maximized->s
      &&  state.lits(i)->args == literal_to_be_maximized->args) {
      return state.lits(i)->value+offset;
    }
  }
  HALT("failed");
  return -10000.;
}


bool MaximizeReward::satisfied(const SymbolicState& state) const {
  // Satisfied = maximum value
  // For count function with 1-arity has maximium value if all constants (except table) fullfill it
  if (literal_to_be_maximized->s->symbol_type == Symbol::count) {
    CountSymbol* cs = (CountSymbol*) literal_to_be_maximized->s;
    if (cs->arity == 0  &&  cs->base_literal->s->arity == 1) {
      double value = SymbolicState::getValue(literal_to_be_maximized->s, state);
      if (state.state_constants.N > 0) {
        return TL::areEqual(value, state.state_constants.N);
      }
      else
        return TL::areEqual(value, reason::getConstants().N);
    }
  }
  else if (literal_to_be_maximized->s->symbol_type == Symbol::function_reward) {
    RewardFunction* rf = dynamic_cast<RewardFunction*>(literal_to_be_maximized->s);
    return reason::holds(state.lits, rf->base_literals);
  }
  return false; // is never satisfied...
}


bool MaximizeReward::possible(const SymbolicState& state) const {
  return true; // is always true
}


void MaximizeReward::getRewardConstants(uintA& constants, const SymbolicState* s) const {
  constants.clear();
   if (literal_to_be_maximized->s->symbol_type == Symbol::count) {
     // Reward constants are those for which predicate instances don't hold yet
     LitL lits;
     Literal::getLiterals(lits, ((CountSymbol*)literal_to_be_maximized->s)->base_literal->s, reason::getConstants(), 1.0);
     lits.memMove = true;
     uint i;
     FOR1D_DOWN(lits, i) {
       if (reason::holds(s->lits, lits(i)))
         lits.remove(i);
     }
     FOR1D(lits, i) {
       constants.setAppend(lits(i)->args);
     }
   }
   else NIY;
//   else if (fa->f->type == Function::function_reward) {
//     RewardFunction* rf = (RewardFunction*) fa->f;
//     reason::getConstants(rf->grounded_pis, constants);
//   }
//   else if (fa->f->type == Function::function_sum) {
//     // Reward constants are all constants in state.
//     reason::getConstants(*s, constants);
//   }
//   else 
//     HALT("Warning: get reward constants has not been implemented in a korrekt way yet");
}


void MaximizeReward::write(ostream& out) const {
  out<<"maximize "<<*literal_to_be_maximized<<endl;
}


void MaximizeReward::write(const char* filename) const {
  ofstream out(filename);
  out<<reward_type<<endl;
  out<<"# MaximizeReward"<<endl;
  out<<*literal_to_be_maximized<<endl;
  out<<reward_maximize_function<<endl;
  out.close();
}




/************************************************
 * 
 *     DisjunctionReward
 * 
 ************************************************/

DisjunctionReward::DisjunctionReward(LitL& _lits) : Reward(reward_one_of_literal_list) {
  this->lits = _lits;
  this->weights.resize(this->lits.N);
  this->weights.setUni(1.0);
}


DisjunctionReward::DisjunctionReward(LitL& _lits, arr& _weights) : Reward(reward_one_of_literal_list) {
  this->lits = _lits;
  this->weights = _weights;
}


double DisjunctionReward::evaluate(const SymbolicState& state) const {
  uint i;
  double max = 0.0;
  FOR1D(lits, i) {
    if (reason::holds(state.lits, lits(i))) {
      max = TL_MAX(max, weights(i));
    }
  }
  return max;
}


bool DisjunctionReward::satisfied(const SymbolicState& state) const {
  uint i;
  FOR1D(lits, i) {
    if (reason::holds(state.lits, lits(i)))
      return true;
  }
  return false;
}


bool DisjunctionReward::possible(const SymbolicState& state) const {
  // TODO This is domain knowledge!
  Symbol* p_OUT = Symbol::get(MT::String("out"));
  uint i, k;
  if (p_OUT != NULL) {
    FOR1D(lits, k) {
      FOR1D(lits(k)->args, i) {
        uintA args(1);  args(0)=lits(k)->args(i);
        Literal* pt_out = Literal::get(p_OUT, args, 1.);
        if (reason::holds(state.lits, pt_out))
          return false;
      }
    }
  }
  return true;
}


void DisjunctionReward::getRewardConstants(uintA& constants, const SymbolicState* s) const {
  constants.clear();
  uint i;
  FOR1D(lits, i) {
    constants.setAppend(lits(i)->args);
  }
}


void DisjunctionReward::write(const char* filename) const {
  ofstream out(filename);
  out<<reward_type<<endl;
  out<<"# DisjunctionReward"<<endl;
  uint i;
  FOR1D(lits, i) {
    out<<"# "; lits(i)->write(out); out<<endl;
  }
  out<<lits.N<<endl;
  FOR1D(lits, i) {
    lits(i)->write(out); out<<endl;
  }
  out.close();
}

void DisjunctionReward::write(ostream& out) const {
  uint i;
  out<<"[" << lits.N  << "]  OR  ";
  FOR1D(lits, i) {
    out<<weights(i)<<":";lits(i)->write(out); out << "  ";
  }
}





/************************************************
 * 
 *     NotTheseStatesReward
 * 
 ************************************************/

NotTheseStatesReward::NotTheseStatesReward(const SymbolicStateL& _undesired_states) : Reward(reward_not_these_states) {
  undesired_states = _undesired_states;
}


double NotTheseStatesReward::evaluate(const SymbolicState& state) const {
  uint i;
  FOR1D(undesired_states, i) {
    if (state == *undesired_states(i))
      return 0.0;
  }
  return 1.0;
}


bool NotTheseStatesReward::satisfied(const SymbolicState& state) const {
  uint i;
  FOR1D(undesired_states, i) {
    if (state == *undesired_states(i))
      return false;
  }
  return true;
}


bool NotTheseStatesReward::possible(const SymbolicState& state) const {
  return true; // is always true
}


void NotTheseStatesReward::getRewardConstants(uintA& constants, const SymbolicState* s) const {
  constants.clear();
  NIY;
}


void NotTheseStatesReward::write(ostream& out) const {
  out<<"reward_not_these_states  "<<endl;
  uint i;
  out<<undesired_states.N<<" undesired states:"<<endl;
  FOR1D(undesired_states, i) {
    out<<i<<": ";
    undesired_states(i)->write(out, true);
    out<<endl;
  }
}

void NotTheseStatesReward::write(const char* filename) const {
  ofstream out(filename);
  out<<reward_type<<endl;
  out<<"# NotTheseStatesReward"<<endl;
  out<<reward_not_these_states<<endl;
  out<<"# Number of undesired states = "<<undesired_states.N << endl;
  uint i;
  FOR1D(undesired_states, i) {
    undesired_states(i)->write(out);
    out<<endl;
  }
  out.close();
}




/************************************************
 * 
 *     RewardConjunction
 * 
 ************************************************/

void RewardConjunction::addReward(Reward *reward) {
  rewards.append(reward);
}

double RewardConjunction::evaluate(const SymbolicState& s) const {
  double product = 1;
  for (uint i = 0; i < rewards.N; i++) {
    product *= rewards(i)->evaluate(s);
  }
  return product;
}

bool RewardConjunction::satisfied(const SymbolicState& s) const {
  for (uint i = 0; i < rewards.N; i++) {
    if (!rewards(i)->satisfied(s)) return false;
  }
  return true;
}
bool RewardConjunction::possible(const SymbolicState& s) const {
  for (uint i = 0; i < rewards.N; i++) {
    if (!rewards(i)->possible(s)) return false;
  }
  return true;
}

void RewardConjunction::getRewardConstants(uintA& constants, const SymbolicState* s) const {
  for (uint i = 0; i < rewards.N; i++) {
    uintA rewardConstants;
    rewards(i)->getRewardConstants(rewardConstants, s);
    constants.setAppend(rewardConstants);
  }
}

void RewardConjunction::write(ostream& out) const {
  cout << "Conjunction of rewards:" << endl;
  for (uint i = 0; i < rewards.N; i++)
    rewards(i)->write(out);
}
void RewardConjunction::write(const char* filename) const {
}






}  // namespace PRADA

std::ostream& operator<<(std::ostream& os, const relational::Reward& p) {
   p.write(os);
   return os;
}
