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


#include "prada.h"
#include "symbols.h"

#define FAST 1


#define MAX_FUNCTION_VALUE 10
#define RULE_MIN_PROB 0.03



namespace relational {

inline uint getIndex__PRADA_Planner(const uintA& constants, const uintA& args) {
  uint args_idx=0;
  uint i;
  FOR1D(args, i) {
    args_idx += ((uint) pow(constants.N, i)) * constants.findValue(args(i));
  }
//   cout<<"getIndex: constants="<<constants<<"  args="<<args<<"    args_idx="<<args_idx<<endl;
  return args_idx;
}

  
  
/************************************************
 * 
 *     PRADA_Planner
 * 
 ************************************************/
  
  
PRADA_Planner::PRADA_Planner() : NID_Planner(), num_samples(100), noise_softener(1.0)  {
  this->dbn = NULL;
  this->prada_reward = NULL;
  this->threshold_reward = 0.05;
  this->reward_calculation__sum = true;
  this->action_choice__max = true;
}


PRADA_Planner::~PRADA_Planner() {
  if (dbn!=NULL) delete dbn;
  if (prada_reward!=NULL) delete prada_reward;
}


void PRADA_Planner::plan_full(LitL& best_plan, double& best_value, const SymbolicState& s, uint max_runs) {
  plan1_wrapper(best_plan, best_value, s, max_runs);
}


Literal* PRADA_Planner::plan_action(const SymbolicState& s, uint max_runs) {
  LitL best_plan;
  double best_value;
  plan_full(best_plan, best_value, s, max_runs);
  if (best_plan.N == 0) {
    return NULL;
  }
  else
    return best_plan(0);
}


void PRADA_Planner::setStartState(const SymbolicState& s1) {
  setState(s1, 0);
}


void PRADA_Planner::setNumberOfSamples(uint num_samples) {
  this->num_samples = num_samples;
}


void PRADA_Planner::setReward(Reward* reward) {
//   MT_MSG("Automatic conversion of reward!");
  this->reward = reward;
  LiteralReward* pg = dynamic_cast<LiteralReward*>(reward);
  if (pg!= NULL) {
    this->prada_reward = convert_reward((LiteralReward*) reward);
  }
  else {
    LiteralListReward* plg = dynamic_cast<LiteralListReward*>(reward);
    if (plg!= NULL) {
      this->prada_reward = convert_reward((LiteralListReward*) reward);
    }
    else {
      MaximizeReward* fg = dynamic_cast<MaximizeReward*>(reward);
      if (fg!= NULL) {
        this->prada_reward = convert_reward((MaximizeReward*) reward);
      }
      else {
        NotTheseStatesReward* ntsg = dynamic_cast<NotTheseStatesReward*>(reward);
        if (ntsg!= NULL)
          this->prada_reward = convert_reward((NotTheseStatesReward*) reward);
        else {
          DisjunctionReward* dg = dynamic_cast<DisjunctionReward*>(reward);
          if (dg!= NULL)
            this->prada_reward = convert_reward((DisjunctionReward*) reward);
          else {
            RewardConjunction *rc = dynamic_cast<RewardConjunction*>(reward);
            if (rc != NULL)
              this->prada_reward = convert_reward(rc);
            else 
              NIY;
          }
        }
      }
    }
  }
}


void PRADA_Planner::setReward(Reward* reward, PRADA_Reward* prada_reward) {
  this->reward = reward;
  this->prada_reward = prada_reward;
}


void PRADA_Planner::setThresholdReward(double threshold_reward) {
  this->threshold_reward = threshold_reward;
}


PRADA_Reward* PRADA_Planner::create_PRADA_Reward(Reward* reward) {
  switch(reward->reward_type) {
    case Reward::reward_literal: return convert_reward((LiteralReward*) reward);
    case Reward::reward_literalList: return convert_reward((LiteralListReward*) reward);
    case Reward::reward_maximize_function: return convert_reward((MaximizeReward*) reward);
    case Reward::reward_not_these_states: return convert_reward((NotTheseStatesReward*) reward);
    case Reward::reward_one_of_literal_list: return convert_reward((DisjunctionReward*) reward);    
    case Reward::conjunction_of_rewards: return convert_reward((RewardConjunction*) reward);    
    default: NIY;
  }
}


void PRADA_Planner::setNoiseSoftener(double noise_softener) {
  if (noise_softener > 1.0 || noise_softener < 0.) {
    MT_MSG("noise_softener has to be in [0,1]");
  }
  else
    this->noise_softener = noise_softener;
}


void PRADA_Planner::setRewardCalculation(bool reward_calculation__sum) {
  this->reward_calculation__sum = reward_calculation__sum;
}


void PRADA_Planner::setActionChoice(bool action_choice__max) {
  this->action_choice__max = action_choice__max;
}



/************************************************
 * 
 *     central run method 
 * 
 ************************************************/

void my_handmade_plan(LitL& plan) {
//   NIY;
  Literal::get(plan, "grab(69) puton(67) grab(66) puton(69)");
#if 0
  uint i;
  Predicate* p_MOVE = logicObjectManager::getPredicate(MT::String("move"));
  uintA args(3);
  FOR1D(reason::getConstants(), i) {
    if (i < 2) continue;
    args(0) = reason::getConstants()(i);
    args(1) = 60;
    args(2) = reason::getConstants()(i-1);
    plan.append(logicObjectManager::getLiteral(p_MOVE, true, args));
  }
#endif
  
  // search-and-rescue
//   char* geiler_plan = "takeoff(11) takeoff(14) goto(12) explore(12) land(12) takeoff(12) goto(11) explore(11) land(11) end-mission()";
  // Triangle-tireworld
//   char* geiler_plan = "move-car(11 14) changetire() move-car(14 17) loadtire(17) changetire() move-car(17 15) loadtire(15) changetire() move-car(15 13)";
  // Ex-Blocksworld 1
//   char* geiler_plan = "pick-up(13 12) put-down(13) pick-up(11 14) put-on-block(11 13) pick-up(14 15) put-down(14) pick-up-from-table(12) put-on-block(12 14)";
  // Ex-Blocksworld 5
//   char* geiler_plan = "pick-up(15 13) put-on-block(15 17) pick-up(13 14) put-on-block(13 16) pick-up(15 17) put-on-block(15 13)";
//   logicObjectManager::getLiterals(plan, geiler_plan);
//   cout<<endl<<endl<<endl<<"USING HAND-MADE PLAN!!"<<endl<<endl<<endl;
//   MT_MSG("USING HAND-MADE PLAN!");
}


void PRADA_Planner::plan1_wrapper(LitL& best_plan, double& best_value, const SymbolicState& s, uint max_runs) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"PRADA_Planner::plan1_wrapper [START]"<<endl;}
  if (DEBUG>0) {
    cout<<"Reward:  ";  if (reward != NULL) reward->write();  else cout<<"NULL";  cout<<endl;
    PRINT(num_samples);
    int size_search_space = pow(ground_actions.N, horizon);
    cout<<endl<<endl<<"########"<<endl
      <<"size_search_space="<<size_search_space<<" (number of possible plans)"<<endl
      <<"(number of actions="<<ground_actions.N<<", horizon="<<horizon<<", number of possible plans="<<ground_actions.N<<"^"<<horizon<<")"<<endl
      <<"With "<<num_samples<<" samples you'll cover potentially "<<(num_samples * 100.0 / size_search_space)<<"%."<<endl
      <<"########"<<endl<<endl;
  }
  
  // Check 1 whether planning makes sense at all:  are there changing concepts?
  SymL changing_symbols;
  ground_rules.changingSymbols(changing_symbols);
  if (changing_symbols.N == 0) {
    best_value = -100000.;
    best_plan.clear();
    if (DEBUG>0) {cout<<"Planning does not make sense here: no changing concepts"<<endl;}
    if (DEBUG>0) {cout<<"PRADA_Planner::plan1_wrapper [END]"<<endl;}
    return;
  }
  
  if (DEBUG>0) {cout<<"SymbolicState:  ";  s.write();  cout<<endl;  PRINT(max_runs);  cout<<"Trying to set state..."<<endl;}
  setStartState(s);
  if (DEBUG>0) {cout<<"SymbolicState has been set."<<endl;}
  
  // Check 2 whether planning makes sense at all:  are reward concepts among random variables?
  if (reward != NULL) {
    LitL lits_reward;
    LiteralReward* pg = dynamic_cast<LiteralReward*>(reward);
    if (pg!= NULL) {
      lits_reward.append(pg->lit);
    }
    else {
      LiteralListReward* plg = dynamic_cast<LiteralListReward*>(reward);
      if (plg!= NULL) {
        lits_reward.append(plg->lits);
      }
    }
    if (DEBUG>0) {cout<<"lits_reward:  ";  write(lits_reward);  cout<<endl;}
    uint i;
    FOR1D(lits_reward, i) {
      LiteralRV* var = dbn->RVefficiency__atom2var(lits_reward(i));
      if (var == NULL) {
        best_value = -100000.;
        best_plan.clear();
        if (DEBUG>0) {cout<<"Planning does not make sense here: reward concept " << *lits_reward(i) << " not among changing concepts"<<endl;}
        if (DEBUG>0) {cout<<"PRADA_Planner::plan1_wrapper [END]"<<endl;}
        return;
      }
    }
  }
  
  
  uint idx_run = 0;
  LitL run_plan;
  double run_value = 0.0;
  do {
    bool improve = plan1(run_plan, run_value, num_samples);
    if (improve) {
      best_plan = run_plan;
      best_value = run_value;
      break;
    }
    if (DEBUG > 0) {
      cerr<<endl<<"PRADA run "<< (idx_run+1) << " of "<<max_runs<<": no good sample found."<<endl;
    }
  } while (++idx_run < max_runs);
  
  if (idx_run == max_runs) {
    if (DEBUG > 0) {
      cerr<<"PRADA: Can't find a good plan, sire!"<<endl;
      cout<<"PRADA: Can't find a good plan, sire!"<<endl;
    }
    best_value = -100000.;
    best_plan.clear();
  }
  if (DEBUG>0) {cout<<"PRADA_Planner::plan1_wrapper [END]"<<endl;}
}


bool PRADA_Planner::plan1(LitL& best_plan, double& bestValue, uint num_samples) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"PRADA_Planner::plan1 [START]"<<endl;}
  
  uint t;

  if (DEBUG>1) {
    PRINT(num_samples);
    PRINT(threshold_reward);
    cout<<"STARTING STATE:  "<<endl;
    dbn->writeState(0, true, 0.9);
    cout<<endl;
    cout<<"REWARD:  ";  reward->write();  cout<<endl;
  }
  
//   FOR1D_(ground_rules, t) {
//     if (ground_rules.elem(t)->action == logicObjectManager::getLiteral("grab(69)"))
//       ground_rules.elem(t)->write();
//   }
 
  // doNothing-value [START]
  uint v, val;
  for (t=1; t<=horizon; t++) {
    FOR1D(dbn->vars_state__prim, v) {
      for (val=0; val<dbn->vars_state__prim(v)->dim; val++) {
        dbn->vars_state__prim(v)->P(t,val) = dbn->vars_state__prim(v)->P(0,val);
      }
    }
    FOR1D(dbn->vars_state__derived, v) {
      for (val=0; val<dbn->vars_state__derived(v)->dim; val++) {
        dbn->vars_state__derived(v)->P(t,val) = dbn->vars_state__derived(v)->P(0,val);
      }
    }
    // Rules -- need to be inferred first!!!  (depending on reward!)
    dbn->inferRules(0);
    if (t<dbn->vars_rules_simple.d0) {
      for (v=0; v<dbn->vars_rules_simple.d1; v++) {
        dbn->vars_rules_simple(t,v) = dbn->vars_rules_simple(0,v);
        dbn->vars_rules(t,v) = dbn->vars_rules(0,v);
      }
    }
  }
  double do_nothing_reward = inferStateRewards();
  if (DEBUG>0) {cout<<"Calculated do_nothing_reward="<<do_nothing_reward<<endl;}
  // doNothing-value [END]
 
 
 
  action_choices.resize(horizon, ground_actions.N);
  action_choices.setUni(0);

  uint s;
  MT::Array< LitL > plans;
  arr reward_values;
  arr action_values;
  arr values;
//   MT::Array< arr > reward_probs(num_samples);
  
  any_sensible_action = false;
  bool good_old_plan__is_good = false;
  for(s=0; s<num_samples; s++) {
    if (DEBUG>2) cout <<"--- SAMPLING ROUND "<<s<<"  (/" << num_samples << ") ---"<<endl;
    if (DEBUG>0 && s%10==0) cerr <<"."<<std::flush;
    LitL sampled_actions;
    
    if (good_old_plans.N > 0 && s==0) {
      LitL old_inspiration_plan = good_old_plans(0);
      if (DEBUG>0) {cout<<"Candidate good old (complete) - s==0:  "<<old_inspiration_plan<<endl;}
      old_inspiration_plan.memMove = true;
      FOR1D_DOWN(old_inspiration_plan, t) {
        if (ground_actions.findValue(old_inspiration_plan(t)) < 0) {
          old_inspiration_plan.remove(t);
        }
      }
      if (DEBUG>0) {cout<<"Retrying filtered good old (complete) - s==0:  "<<old_inspiration_plan<<endl;}
      sampleActionsAndInfer(sampled_actions, old_inspiration_plan);
    }
    else if (good_old_plans.N > 0 && s==1) {
      uint q;
      LitL old_inspiration_plan;
      FOR1D(good_old_plans(0), q) {
        if (q==0) continue;
        if (ground_actions.findValue(good_old_plans(0)(q)) < 0) continue;
        old_inspiration_plan.append(good_old_plans(0)(q));
      }
      if (DEBUG>0) {cout<<"Retrying good old (kill action 0) - s==1:  "<<old_inspiration_plan<<endl;}
      sampleActionsAndInfer(sampled_actions, old_inspiration_plan);
    }
    else if (good_old_plans.N > 0 && s==2) {
      uint q;
      LitL old_inspiration_plan;
      FOR1D(good_old_plans(0), q) {
        if (q==1) continue;
        if (ground_actions.findValue(good_old_plans(0)(q)) < 0) continue;
        old_inspiration_plan.append(good_old_plans(0)(q));
      }
      if (DEBUG>0) {cout<<"Retrying good old (kill action 1) - s==2:  "<<old_inspiration_plan<<endl;}
      sampleActionsAndInfer(sampled_actions, old_inspiration_plan);
    }
    else if (good_old_plans.N > 0 && s==3) {
      uint q;
      LitL old_inspiration_plan;
      FOR1D(good_old_plans(0), q) {
        if (q==0 || q==1) continue;
        if (ground_actions.findValue(good_old_plans(0)(q)) < 0) continue;
        old_inspiration_plan.append(good_old_plans(0)(q));
      }
      if (DEBUG>0) {cout<<"Retrying good old (kill actions 0 and 1) - s==3:  "<<old_inspiration_plan<<endl;}
      sampleActionsAndInfer(sampled_actions, old_inspiration_plan);
    }
    else sampleActionsAndInfer(sampled_actions);
    
#if 0
    if (s==0) {
      LitL hand_plan;
      my_handmade_plan(hand_plan);
      // ensure that all actions are possible
      uint q;
//       FOR1D(hand_plan, q) {
//         if (dbn->RVefficiency__atom2var(hand_plan(q)) == NULL)
//           HALT("handmade plan undoable / not recognized");
//       }
      while (hand_plan.N < horizon) {
        hand_plan.append(NULL);
      }
      cout<<"Hand-made plan:  ";  write(hand_plan);  cout<<endl;
      sampleActionsAndInfer(sampled_actions, hand_plan, dbn, horizon);
      cout<<"---> Reward = " << inferStateRewards(hand_plan.N) << endl;
      exit(0);
    }
#endif
    
    plans.append(sampled_actions);
    if (DEBUG>2) {cout<<"Sampled actions: ";write(sampled_actions); cout<<endl;}
    // bis zur ersten NULL evaluieren...
    FOR1D(sampled_actions, t) {
      if (sampled_actions(t) == NULL) {
        break;
      }
    }
    reward_values.append(inferStateRewards(t));  // t ist hier horizont.
    action_values.append(0.);
    if (use_ruleOutcome_rewards) {
      action_values(s) = calcRuleRewards(sampled_actions);
    }
    values.append(reward_values(s) + action_values(s));
    // statistics
    FOR1D(sampled_actions, t) {
      uint action_id = ground_actions.findValue(sampled_actions(t));
      action_choices(t, action_id)++;
    }
    
    CHECK(plans.N == reward_values.N   &&  plans.N == action_values.N  &&  plans.N == values.N, "");
    
    if (DEBUG>2) {cout<<"--> values(s)="<<values(s)<<"  (reward_values(s)="<<reward_values(s)<<",  action_values(s)="<<action_values(s)<<")"<<endl;}
    
    
    if (!use_ruleOutcome_rewards  &&  reward->reward_type != Reward::reward_maximize_function) {
//       if (any_sensible_action  &&  s>300) {
//         if (DEBUG>0) {cout<<"Great plan already found (s="<<s<<")!"<<endl;}
//         break;
//       }
    }
    else if (!use_ruleOutcome_rewards  &&  reward->reward_type == Reward::reward_maximize_function) {
      if (values(s) > do_nothing_reward + threshold_reward  &&  s<4) {
        good_old_plan__is_good = true;
      }
      else if (good_old_plans.N > 0  &&  values(s) > do_nothing_reward + threshold_reward  &&  s>200  &&  good_old_plan__is_good) {
//       else if (good_old_plans.N > 0  &&  values(s) > do_nothing_reward + threshold_reward  &&  s>500  &&  good_old_plan__is_good) {
// MT_MSG("STACK version!");
        if (DEBUG>0) {cout<<"Great good old plan (s="<<s<<")!"<<endl;}
        break;
      }
      else if (values(s) > do_nothing_reward + threshold_reward  &&  s>800) {
        if (DEBUG>0) {cout<<"Great plan already found (s="<<s<<")!"<<endl;}
        break;
      }
    }
  }
  if (DEBUG>0) {cerr<<s<<" samples taken."<<endl;}

  int max_id = values.maxIndex();
  double best_value = values(max_id);
  Literal* max_action = NULL;
  if (plans(max_id).N > 0)
   max_action = plans(max_id)(0);
  if (DEBUG>0) {
    if (DEBUG>2) {
      cout<<"Sample values:  "<<endl;
      FOR1D(plans, s) {
        cout<<"("<<s<<") "<<values(s)<<":  ";write(plans(s));cout<<endl;
      }
    }
    arr sorted_values;
    uintA sortedIndices;
    TL::sort_desc(sorted_values, sortedIndices, values);
    cout<<"All plans sorted (horizon="<< horizon << "):"<<endl;
    FOR1D(sortedIndices, s) {
      if (DEBUG<3 && s>10) break;
      printf("#%4u",s);
      cout<<":   ";
      printf("(%4u)", sortedIndices(s));
      cout<< "  ";
      printf("%5.3f", values(sortedIndices(s)));
      printf(" (%5.3f + %5.3f)", reward_values(sortedIndices(s)), action_values(sortedIndices(s)));
      cout<< "  ";
      write(plans(sortedIndices(s)));
      cout<<endl;
    }
    
    PRINT(max_id);
    if (!any_sensible_action) {cout<<"NO SENSIBLE ACTION WAS FOUND!"<<endl;}
    cout<<"Max action:  ";
    if (max_action != NULL) {
      max_action->write();cout<<"  v="<<values(max_id)<<""<<endl;
    }
    else
      cout<<" -"<<endl;
    cout<<"Max sequence: ("<<max_id<<") ";write(plans(max_id));cout<<endl;
    reward->write();cout<<endl;
  }
  
   // statistics
  if ((DEBUG>0 && ground_actions.N < 20)
       || DEBUG>1 ) {
//   if (true) {
    cout<<"Action application statistics:"<<endl;
    printf("%-14s","Time");
    for(t=0;t<horizon;t++) printf("%5u",t);
    cout<<endl;
    uint i;
    FOR1D(ground_actions, i) {
      uint number = 0;
      for(t=0;t<action_choices.d0;t++) {
        number += action_choices(t,i);
      }
      if (number == 0)
        continue;
      MT::String name;  name << *ground_actions(i);
      printf("%-14s",(char*)name);
      for(t=0;t<action_choices.d0;t++) printf("%5u", action_choices(t,i));
      cout<<endl;
    }
  }
  
  // ACTION CHOICE
  // (1)  Choose maximizing action sequence
  if (action_choice__max) {
    best_plan = plans(max_id);
    bestValue = best_value;
  }
  // (2)  Choose action with maximal normalized sum
  else {
    MT::Array< arr > action_values__all(ground_actions.N);
    FOR1D(plans, s) {
      if (plans(s)(0) != NULL)
        action_values__all(ground_actions.findValue(plans(s)(0))).append(values(s));
    }
    uint a;
    arr action_values(ground_actions.N);
    FOR1D(ground_actions, a) {
      if (action_values__all(a).N > 0) {
        #if 0
        action_values(a) = sum(action_values__all(a));
        action_values(a) /= 1.0 * action_values__all(a).N;
        #else
        // taking only X best into account
        uint X = 10;
        uint q;
        arr sorted;
        uintA sortedIndices;
        TL::sort_desc(sorted, sortedIndices, action_values__all(a));
        for (q=0; q<TL_MIN(sorted.N, X); q++) {
          action_values(a) += sorted(q);
        }
        action_values(a) /= 1.0 * TL_MIN(sorted.N, X);
        #endif
        
      }
      else
        action_values(a) = -10000.;
    }
    
    if (DEBUG>0) {
      FOR1D(ground_actions, a) {
        MT::String name;  name << *ground_actions(a);
        printf("%-14s",(char*)name);
        cout << "  " << action_values(a) <<  "  "<<action_values__all(a).N;
        cout<<action_values__all(a)<<endl;
      }
    }
    
    uint max_id2 = action_values.maxIndex();
    bestValue = action_values(max_id2);
    // Random choice of plan: just use first plan with the best action
    FOR1D(plans, s) {
      if (plans(s)(0) == ground_actions(max_id2)) {
        best_plan = plans(s);
        break;
      }
    }
  }


  // Return value: depends on type of reward
  if (use_ruleOutcome_rewards) {
    if (DEBUG>0) {cout<<"PRADA_Planner::plan1 [END]"<<endl;}
    return true;
  }
  else if (reward->reward_type == Reward::reward_maximize_function) {
    if (DEBUG>0) {
      PRINT(bestValue);  PRINT(do_nothing_reward);  PRINT(threshold_reward);  PRINT(do_nothing_reward+threshold_reward);
      PRINT((bestValue > do_nothing_reward + threshold_reward));
      cout<<"Returning-1 "<<(bestValue > do_nothing_reward + threshold_reward)<<endl;
    }
    if (DEBUG>0) {cout<<"PRADA_Planner::plan1 [END]"<<endl;}
    return (bestValue > do_nothing_reward + threshold_reward);
  }
  else {
    if (DEBUG>0) cout<<"Returning02 "<<any_sensible_action<<endl;
    if (DEBUG>0) {cout<<"PRADA_Planner::plan1 [END]"<<endl;}
    return any_sensible_action;
  }
}



/************************************************
 * 
 *     DBN construction
 * 
 ************************************************/


void calc_dbn_state_symbols_for_rewards(SymL& reward_symbols, Reward *r) {
  reward_symbols.clear();
  uint k;
  if (r->reward_type == Reward::reward_literal) {
    reward_symbols.setAppend(((LiteralReward*) r)->lit->s);
    SymL defining_symbols;
    reward_symbols.last()->getDefiningSymbols(defining_symbols, false);
    reward_symbols.setAppend(defining_symbols);
  }
  else if (r->reward_type == Reward::reward_literalList) {
    FOR1D(((LiteralListReward*) r)->lits, k) {
      reward_symbols.setAppend(((LiteralListReward*) r)->lits(k)->s);
      SymL defining_symbols;
      reward_symbols.last()->getDefiningSymbols(defining_symbols, false);
      reward_symbols.setAppend(defining_symbols);
    }
  }
  else if (r->reward_type == Reward::reward_maximize_function) {
    MaximizeReward* mfr = (MaximizeReward*) r;
    if (mfr->literal_to_be_maximized != NULL) {
      reward_symbols.setAppend(mfr->literal_to_be_maximized->s);
      SymL defining_symbols;
      mfr->literal_to_be_maximized->s->getDefiningSymbols(defining_symbols, false);
      reward_symbols.setAppend(defining_symbols);
    }
  }
  else if (r->reward_type == Reward::conjunction_of_rewards) {
    for (uint i = 0; i < ((RewardConjunction*)r)->rewards.N; i++) {
      SymL reward_symbols__inner;
      calc_dbn_state_symbols_for_rewards(reward_symbols__inner, ((RewardConjunction*)r)->rewards(i));
      reward_symbols.setAppend(reward_symbols__inner);
    }
  }
  else
    NIY;
}



// The DBN used by (A-)PRADA consists of grounded predicate and 
// function instances which are used 
//       (i)  in the (ground) rules and/or 
//       (ii) in the reward description.
void PRADA_Planner::calc_dbn_state_symbols() {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"PRADA_Planner::calc_dbn_state_symbols [START]"<<endl;}
  CHECK(reward != NULL, "Reward has to be set first");
  CHECK(ground_rules.num() > 0, "Ground rules have to be set first");
  dbn_state_symbols.clear();
  uint i, k, l;
  // (i-a) rule concepts - primitive
  FOR1D_(ground_rules, i) {
    FOR1D(ground_rules.elem(i)->context, k) {
      if (ground_rules.elem(i)->context(k)->s->symbol_type == Symbol::primitive) {
        dbn_state_symbols.setAppend(ground_rules.elem(i)->context(k)->s);
      }
    }
    FOR1D(ground_rules.elem(i)->outcomes, k) {
      FOR1D(ground_rules.elem(i)->outcomes(k), l) {
        if (ground_rules.elem(i)->outcomes(k)(l)->s->symbol_type == Symbol::primitive) {
          dbn_state_symbols.setAppend(ground_rules.elem(i)->outcomes(k)(l)->s);
        }
      }
    }
  }
  // (i-b) rule concepts - derived
  FOR1D_(ground_rules, i) {
    FOR1D(ground_rules.elem(i)->context, k) {
      if (ground_rules.elem(i)->context(k)->s->symbol_type != Symbol::primitive) {
        dbn_state_symbols.setAppend(ground_rules.elem(i)->context(k)->s);
        SymL defining_symbols;
        dbn_state_symbols.last()->getDefiningSymbols(defining_symbols, false);
        dbn_state_symbols.setAppend(defining_symbols);
      }
    }
  }
  // (ii) reward concepts
  SymL symbols_reward;
  calc_dbn_state_symbols_for_rewards(symbols_reward, reward);
  dbn_state_symbols.setAppend(symbols_reward);

  Symbol::sort(dbn_state_symbols);
  if (DEBUG>0) {
    cout<<"(A-) PRADA's DBN will be built with random variables for the following state symbols:"<<endl;
    writeSymbols(dbn_state_symbols);
  }
  if (DEBUG>0) {cout<<"PRADA_Planner::calc_dbn_state_symbols [END]"<<endl;}
}


void PRADA_Planner::build_dbn(const uintA& constants, const SymL& preds, const SymL& actions) {
  dbn = new PRADA_DBN(constants, preds, actions, ground_rules, noise_softener, horizon);
}


void PRADA_Planner::build_dbn(const uintA& constants) {
  calc_dbn_state_symbols();
  SymL symbols_action;
  Symbol::get_action(symbols_action);
  build_dbn(constants, dbn_state_symbols, symbols_action);
}




/************************************************
 * 
 *     Inference in DBN
 * 
 ************************************************/

void PRADA_Planner::infer(const LitL& plan) {
  LitL sampled_actions;
  sampleActionsAndInfer(sampled_actions, plan, dbn, plan.N);
}


double PRADA_Planner::inferStateRewards() {
  return inferStateRewards(horizon);
}


double PRADA_Planner::inferStateRewards(uint horizon) {
  if (reward_calculation__sum) {
    return inferStateRewards_limited_sum(horizon);
  }
  else {
    return inferStateRewards_limited_max(horizon);
  }
}


double PRADA_Planner::inferStateRewards_limited_max(uint horizon) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"inferStateRewards_limited_max [START]"<<endl;}
  uint t;
  if (DEBUG>0) {
    reward->write();cout<<endl;
  }
  double value_t, value_0, value_max=-50.;
  
  for (t=0; t<=horizon; t++) {
//     if (t>=T) break;
    value_t = discount_pow(t) * inferStateRewards_single(t);
    if (t==0) value_0 = value_t;
    else {
      if (value_t - value_0 > threshold_reward) {
        any_sensible_action = true;
      }
    }
    if (value_t > value_max)
      value_max = value_t;
    if (DEBUG>1) {cout<<"t="<<t<<":  "<<value_t<<endl;}
  }
  
  if (DEBUG>0) {
    PRINT(value_0);
    PRINT(any_sensible_action);
  }
  if (DEBUG>0) {cout<<"--> value="<<value_max<<endl;}
  if (DEBUG>0) {cout<<"inferStateRewards_limited_max [END]"<<endl;}
  return value_max;
}


double PRADA_Planner::inferStateRewards_limited_sum(uint horizon) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"inferStateRewards_limited_sum [START]"<<endl;}
  uint t;
  if (DEBUG>1) {
    reward->write();cout<<endl;
    PRINT(threshold_reward);
    PRINT(any_sensible_action);
  }
  double value_t, value_0=0.0, value_total=0.;
  arr debug_values;
  for (t=0; t<=horizon; t++) {
//     if (t>=T) break;
    value_t = inferStateRewards_single(t);
    if (t==0) value_0 = value_t;
    else {
      if (value_t - value_0 > threshold_reward)
        any_sensible_action = true;
    }
    value_total += discount_pow(t) * value_t;
    debug_values.append(value_t);
    if (DEBUG>1) {cout<<"t="<<t<<":  "<<value_t<<endl;}
  }
  
  if (DEBUG>0) {
    LitL actions;
    dbn->getActions(actions, horizon);
    write(actions);  cout<<endl;
    PRINT(debug_values);
    PRINT(value_0);
    PRINT(any_sensible_action);
  }
  if (DEBUG>0) {cout<<"--> value="<<value_total<<endl;}
  if (DEBUG>0) {cout<<"inferStateRewards_limited_sum [END]"<<endl;}
  return value_total;
}


double PRADA_Planner::inferStateRewards_single(uint t) {
  return prada_reward->evaluate_prada_reward(*dbn, t);
}


double PRADA_Planner::calcRuleRewards(const LitL& actions) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcRuleRewards [START]"<<endl;}
  double total_reward = 0.;
  uint t, r;
  FOR1D(actions, t) {
    if (DEBUG>0) {cout<<"+++++  t="<<t<<"  ";  actions(t)->write();  cout<<endl;}
    NIY;
    uint idx_symbol = 11111111;
//     uint idx_symbol = logicObjectManager::p_actions.findValue(actions(t)->s);
    uint idx_args = getIndex__PRADA_Planner(reason::getConstants(), actions(t)->args);
    double t_reward = 0.;
    for (r=0; r<dbn->action2rules_num(idx_symbol, idx_args); r++) {
      uint rule_number = dbn->action2rules(idx_symbol, idx_args, r);
      Rule* rule = ground_rules.elem(rule_number);
      double expected_reward = this->expected_rule_rewards(rule_number);
      double probability = dbn->vars_rules(t, rule_number);
      t_reward += discount_pow(t) * probability * expected_reward;
      if (DEBUG>1) {
        if (dbn->vars_rules(t, rule_number) > 0.01) {
          cout<<"Rule #"<<r<<endl;
          rule->write(cout);
          cout<<" expected_reward="<<expected_reward << "   *   probability="<<probability<<"    *   discount="<<discount_pow(t)<<endl;
          cout<<" -> Brings "<< (discount_pow(t) * probability * expected_reward) << endl;
        }
      }
    }
    if (DEBUG>0) {cout<<"==> t_reward="<<t_reward<<endl;}
    total_reward += t_reward;
  }
  if (DEBUG>0) {PRINT(total_reward);}
  if (DEBUG>0) {cout<<"calcRuleRewards [END]"<<endl;}
  return total_reward;
}



/************************************************
 * 
 *     Write DBN
 * 
 ************************************************/

void PRADA_Planner::writeState(uint t, ostream& out) {
  dbn->writeState(t, out);
}


void PRADA_Planner::writeStateSparse(uint t, bool prim_only, ostream& out) {
  dbn->writeStateSparse(t, prim_only, out);
}



/************************************************
 * 
 *     PRADA_Reward
 * 
 ************************************************/

class PRADA_Reward_Literal : public PRADA_Reward {
  Literal* lit;
  LiteralRV* rv;
  uint value;
  
  public:
    PRADA_Reward_Literal(Literal* lit) {
      this->lit = lit;
      rv = NULL;
      if (lit->s->range_type == Symbol::binary) {
        if (lit->value > 0.)
          value = 1;
        else
          value = 0;
      }
      else {
        uint rangeIndex = lit->s->range.findValue((uint)lit->value);
        CHECK(rangeIndex >= 0, "Value out of range!")
        value = rangeIndex;
      }
    }
    
    double evaluate_prada_reward(const PRADA_DBN& dbn, uint t) {
      if (rv == NULL) {
        rv = dbn.RVefficiency__atom2var(lit);
      }
      return rv->P(t,value);
    }
};


class PRADA_Reward_LiteralList : public PRADA_Reward {
  LitL lits;
  RVL rvs;
  
  public:
    PRADA_Reward_LiteralList(LitL& lits) {
      this->lits = lits;
      rvs.resize(this->lits.N);
      rvs.setUni(NULL);
    }
    
    double evaluate_prada_reward(const PRADA_DBN& dbn, uint t) {
      uint i;
      if (rvs(0) == NULL) {
        FOR1D(lits, i) {
          rvs(i) = dbn.RVefficiency__atom2var(lits(i));
        }
      }
      double prob = 1.0;
      FOR1D(rvs, i) {
        if (lits(i)->s->range_type == Symbol::binary) {
          if (lits(i)->value > 0.)
            prob *= rvs(i)->P(t,1);
          else
            prob *= rvs(i)->P(t,0);
        }
        else {
          uint rangeIndex = lits(i)->s->range.findValue((uint)lits(i)->value);
          CHECK(rangeIndex >= 0, "Value out of range!")
          prob *= rvs(i)->P(t, rangeIndex);
        }
      }
      return prob;
    }
};


class PRADA_Reward_Disjunction : public PRADA_Reward {
  LitL lits;
  RVL rvs;
  arr weights;
  
  public:
    PRADA_Reward_Disjunction(LitL& lits, arr& weights) {
      this->lits = lits;
      rvs.resize(this->lits.N);
      rvs.setUni(NULL);
      this->weights = weights;
    }
    
    double evaluate_prada_reward(const PRADA_DBN& dbn, uint t) {
      uint i;
      if (rvs(0) == NULL) {
        FOR1D(lits, i) {
          rvs(i) = dbn.RVefficiency__atom2var(lits(i));
        }
      }
      double max_prob = 0.0;
      FOR1D(rvs, i) {
        if (lits(i)->value > 0.)
          max_prob = weights(i) * TL_MAX(max_prob, rvs(i)->P(t,1));
        else
          max_prob = weights(i) * TL_MAX(max_prob, rvs(i)->P(t,0));
      }
      return max_prob;
    }
};

class PRADA_Reward_ConjunctionOfRewards : public PRADA_Reward {
  MT::Array<PRADA_Reward*> pradaRewards;
  
  public:
    PRADA_Reward_ConjunctionOfRewards(const MT::Array<PRADA_Reward*> &rewards) {
      pradaRewards = rewards;
    }
    
    double evaluate_prada_reward(const PRADA_DBN& dbn, uint t) {
      double product = 1;
      for (uint i = 0; i < pradaRewards.N; i++) {
        product *= pradaRewards(i)->evaluate_prada_reward(dbn, t);
      }
      return product;
    }
};


class PRADA_Reward_Maximize : public PRADA_Reward {
  Literal* fl;
  LiteralRV* rv;
  double offset;
  
  public:
    PRADA_Reward_Maximize(Literal* fl, double offset) {
      this->fl = fl;
      this->offset = offset;
      rv = NULL;
    }
    
    double evaluate_prada_reward(const PRADA_DBN& dbn, uint t) {
      if (rv == NULL) {
        rv = dbn.RVefficiency__atom2var(fl);
      }
      return rv->P(t,0)+offset;	//quick and dirty, need this for reward conjunction
    }
};


class PRADA_Reward_S : public PRADA_Reward {
  SymbolicStateL undesired_states;
  double penalty;
  
  public:
    PRADA_Reward_S(SymbolicStateL& _undesired_states) {
      undesired_states = _undesired_states;
      penalty = -1.;
    }
    
    double evaluate_prada_reward(const PRADA_DBN& dbn, uint t) {
      uint DEBUG = 0;
      double value = 10.;
      uint i;
      if (DEBUG>0) {cout<<"=== t="<<t<<" ===="<<endl;   dbn.writeState(t, true, 0.1);}
      FOR1D(undesired_states, i) {
        if (DEBUG>0) {cout<<"-----"<<endl; undesired_states(i)->write(cout); cout<<endl;}
        double state_log_probability = dbn.log_probability(t, *undesired_states(i));
        if (DEBUG>0) {PRINT(state_log_probability);  PRINT(exp(state_log_probability));  PRINT(penalty * exp(state_log_probability));}
        value += penalty * exp(state_log_probability);
      }
      if (DEBUG>0) {cout<<"-->  ";  PRINT(value);}
      return value;
    }
};


PRADA_Reward* PRADA_Planner::convert_reward(LiteralReward* reward) {
  return new PRADA_Reward_Literal(reward->lit);
}


PRADA_Reward* PRADA_Planner::convert_reward(LiteralListReward* reward) {
  return new PRADA_Reward_LiteralList(reward->lits);
}


PRADA_Reward* PRADA_Planner::convert_reward(MaximizeReward* reward) {
  return new PRADA_Reward_Maximize(reward->literal_to_be_maximized, reward->offset);
}


PRADA_Reward* PRADA_Planner::convert_reward(NotTheseStatesReward* reward) {
  return new PRADA_Reward_S(reward->undesired_states);
}


PRADA_Reward* PRADA_Planner::convert_reward(DisjunctionReward* reward) {
  return new PRADA_Reward_Disjunction(reward->lits, reward->weights);
}


PRADA_Reward* PRADA_Planner::convert_reward(RewardConjunction* reward) {
  MT::Array<PRADA_Reward*> pradaRewards;
  for (uint i = 0; i < reward->rewards.N; i++) {
    pradaRewards.append(create_PRADA_Reward(reward->rewards(i)));
  }
  return new PRADA_Reward_ConjunctionOfRewards(pradaRewards);
}



/************************************************
 * 
 *     setState
 * 
 ************************************************/

void PRADA_Planner::setState(const SymbolicState& s, uint t) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"PRADA_Planner::setState [START]"<<endl;}
  // Construct DBN if necessary
  if (dbn == NULL) {
    if (DEBUG>0) {cout<<"Constructing dbn."<<endl;}
    build_dbn(reason::getConstants());
  }
  if (DEBUG>0) {cout<<"Creating lits_prim_filtered and fv_prim_filtered."<<endl;}
    if (DEBUG>0) {cout<<"Setting state: "<<s<<endl;}
  
  // Set
  LitL lits_prim_filtered;
  uint i;
  FOR1D(s.lits, i) {
    if (dbn_state_symbols.findValue(s.lits(i)->s) >= 0)
      lits_prim_filtered.append(s.lits(i));
  }
  
  if (DEBUG>0) {cout<<"Setting in dbn."<<endl;}
  dbn->setState(lits_prim_filtered, t);
  if (DEBUG>0) {cout<<"PRADA_Planner::setState [END]"<<endl;}
}


/************************************************
 * 
 *     sampleActionsAndInfer
 * 
 ************************************************/

void PRADA_Planner::sampleActionsAndInfer(LitL& sampled_actions) {
  LitL fixed_actions(horizon);
  fixed_actions.setUni(NULL);
  sampleActionsAndInfer(sampled_actions, fixed_actions);
}


void PRADA_Planner::sampleActionsAndInfer(LitL& sampled_actions, const LitL& fixed_actions) {
  sampleActionsAndInfer(sampled_actions, fixed_actions, dbn, horizon);
}


void PRADA_Planner::sampleActionsAndInfer(LitL& sampled_actions, const LitL& fixed_actions1, PRADA_DBN* local_net, uint local_horizon) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"sampleActionsAndInfer [START]"<<endl;}
  LitL fixed_actions;  fixed_actions = fixed_actions1;
  if (DEBUG>0) {cout<<"fixed_actions: " << fixed_actions <<endl;}
  if (fixed_actions.N < local_horizon) {
    if (DEBUG>>0) {cout<<"fixed_actions have wrong number " << fixed_actions.N << " instead of " << local_horizon << " and will be refilled"<<endl;}
    do {
      fixed_actions.append(NULL);
    } while (fixed_actions.N < local_horizon);
  }
  bool TIMING = false;
  clock_t start,finish;
  double time;
  sampled_actions.clear();
  uint A = ground_actions.N;
  arr action_weights(local_horizon, A);
  uint t, r, a;
  for (t=0; t<local_horizon; t++) {
    if (DEBUG>1) {cout<<"+++++ TIME-STEP "<<t<<" +++++"<<endl;}
    if (DEBUG>1) {local_net->writeState(t, false, 0.05);cout<<endl;}
    
    // (1) INFER rule random variables at t
    if (TIMING) {start = clock();}
    local_net->inferRules(t);
    if (TIMING) {
      finish = clock();
      time = (double(finish)-double(start))/CLOCKS_PER_SEC;
      cerr<<"inferPhis"<<time<<endl;
    }
    
    // (2) INFER ACTION-WEIGHT at t,  cf. paper Eq. (24)
    if (TIMING) {start = clock();}
    arr action_weights__sampling(ground_actions.N);  // --> need this for sampling distribution which shall ignore non-manipulating rules
    FOR1D(ground_actions, a) {
      action_weights(t,a) = 0.0;
      action_weights__sampling(a) = 0.0;
      uint idx_symbol = dbn->net_symbols_action.findValue(ground_actions(a)->s);
      uint idx_args = getIndex__PRADA_Planner(reason::getConstants(), ground_actions(a)->args);
      for (r=0; r<local_net->action2rules_num(idx_symbol, idx_args); r++) {
        action_weights(t,a) += local_net->vars_rules(t, local_net->action2rules(idx_symbol, idx_args, r));
        int idx = local_net->action2rules(idx_symbol, idx_args, r); // this is a fix because sometimes idx is not in is_manipulating_rule
        if (Literal::getLiteral_doNothing() == ground_actions(a)  || (is_manipulating_rule.N > idx && is_manipulating_rule(idx) ))  // only consider manipulating rules
          action_weights__sampling(a) += local_net->vars_rules(t, local_net->action2rules(idx_symbol, idx_args, r));
      }
      // Manually forbid doNothing-action at first time-step if free choice 
      if (t==0  &&  Literal::getLiteral_doNothing() == ground_actions(a)  &&  fixed_actions(t) == NULL) {
        action_weights(t,a) = 0.0;
        action_weights__sampling(a) = 0.0;
      }
    }
    if (DEBUG>2) {
      cout<<"action_weights__sampling (action applicability):";
      for(a=0; a<A; a++) {
        cout<<"  "<<action_weights__sampling(a);
      }
      cout<<endl;
      cout<<"APPLICABLE ACTIONS:   ";
      for(a=0; a<A; a++) {
        if (action_weights(t,a)>0) {
          cout<<"   "<<action_weights__sampling(a)<<" ("<<action_weights(t,a)<<") ";
          ground_actions(a)->write();
        }
      }
      cout<<endl;
    }
    if (TL::isZero(sum(action_weights__sampling))) {
      if (DEBUG>0) {cout<<"No more action applicable. --> Stop!"<<endl;};
      break;
    }
    arr P_action = action_weights__sampling / sum(action_weights__sampling);
        
    // (3) SELECT ACTION
    int sampled_action_id;
    Literal* sampled_action;
    // (3a) take fixed action
    if (fixed_actions(t) != NULL) {
      sampled_action_id = ground_actions.findValue(fixed_actions(t));
      sampled_action = fixed_actions(t);
      if (sampled_action_id < 0) {
        HALT("fixed_actions(t)="<<*fixed_actions(t)<<"  has not been found in ground_actions "<<ground_actions);
      }
      if (DEBUG>1) {cout<<"*** SIMULATED ACTION (fixed):   ";  sampled_action->write(cout);cout<<endl;}
    }
    // (3b) or sample action according to sampling distribution
    else {
      sampled_action_id = TL::basic_sample(P_action);
      sampled_action = ground_actions(sampled_action_id);
      if (DEBUG>1) {cout<<"*** SIMULATED ACTION (sampled):   ";  sampled_action->write(cout);cout<<endl;}
    }
    sampled_actions.append(sampled_action);
    FOR1D(local_net->vars_action, a) {
      if (local_net->vars_action(a)->lit == sampled_action) {
        local_net->vars_action(a)->P(t,0) = 0.;
        local_net->vars_action(a)->P(t,1) = 1.;
      }
      else {
        local_net->vars_action(a)->P(t,0) = 1.;
        local_net->vars_action(a)->P(t,1) = 0.;
      }
    }
    if (DEBUG>2) {
      cout<<"P_action:  ";
      uint i;
      FOR1D(P_action, i) {
        if (P_action(i) > 0.) {
          cout<<P_action(i)<< " ";  ground_actions(i)->write();  cout<<"   ";
        }
      }
      cout<<endl;
      PRINT(sampled_action_id);
      if (fixed_actions(t) != NULL) {cout<<"Fixed ";}
      else {cout<<"Sampled ";}
      cout<<"action: ";sampled_action->write();cout<<endl;
//       cerr<<"Action: ";sampled_action->write(cerr);cerr<<endl;
    }
    if (TIMING) {
      finish = clock();
      time = (double(finish)-double(start))/CLOCKS_PER_SEC;
      cerr<<"TIME calcaction_weights and sample action "<<time<<endl;
    }
    
    
    // (4) INFER STATE t+1   (--> propagate action effects)
    if (TIMING) {start = clock();}
    local_net->inferState(t+1, sampled_action, action_weights(t,sampled_action_id));
    if (TIMING) {
      finish = clock();
      time = (double(finish)-double(start))/CLOCKS_PER_SEC;
      cerr<<"TIME infer state t+1 "<<time<<endl;
    }
  }

  if (DEBUG>1) {cout<<"+++++ TIME-STEP "<<t<<" +++++ (sans action)"<<endl;}
  if (DEBUG>1) {local_net->writeState(t, false, 0.05);}
  if (DEBUG>0) {cout<<"sampleActionsAndInfer [END]"<<endl;}
}


















/************************************************
 * 
 *     A_PRADA
 * 
 ************************************************/

double A_PRADA::shorten_plan(LitL& seq_best, const LitL& seq_old, double value_old) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"A_PRADA::shorten_plan [START]"<<endl;}
  if (DEBUG>0) {
    cout<<"seq_old: ";write(seq_old);cout<<endl;
    PRINT(value_old);
  }
  if (seq_old.N == 0) {
    seq_best.clear();
    if (DEBUG>0) {cout<<"No sensible input sequence"<<endl;  cout<<"A_PRADA::shortenSeq [END]"<<endl;}
    return -5999999999.;
  }
  double value_best = value_old;
  double value = 0.;
  LitL seq_new;
  seq_best = seq_old;
  uint DEBUG_TOTAL_R=0;
  // always leave out t
  uint t=0, t2;
  while (t<horizon-1 && t<seq_best.N && seq_best.N > 1) {
//   int t=T-2, t2;
//   while (t>=0) {
//     PRINT(DEBUG_TOTAL_R);
    if (DEBUG_TOTAL_R++>1000) {
      cerr<<"too many round"<<endl;
      break;
    }
    if (DEBUG>1) {cout<<"*** t="<<t<<" ***"<<endl;}
    if (DEBUG>1) {cout<<"seq_best:   ";write(seq_best);cout<<endl;}
    #if 0
    // A-PRADA samplet auch noch die letzten Sequenzen
    LitL seq_fixed(horizon);
    seq_fixed.setUni(NULL);
    for (t2=0; t2<seq_best.N; t2++) {
      if (t2 < t)
        seq_fixed(t2) = seq_best(t2);
      else if (t2 > t)
        seq_fixed(t2-1) = seq_best(t2);
    }
    if (DEBUG>2) {cout<<"t="<<t<<"  --  seq_fixed:  ";write(seq_fixed);cout<<endl;}
    //     sampleActionsAndInfer(seq_new, seq_fixed, dbn, horizon)
    if (DEBUG>2) {cout<<" seq_fixed + sampled:  ";write(seq_new);cout<<endl;}
    #else
    LitL seq_new;
    for (t2=0; t2<seq_best.N; t2++) {
      if (t2 < t)
        seq_new.append(seq_best(t2));
      else if (t2 > t)
        seq_new.append(seq_best(t2));
    }
    infer(seq_new);
    if (DEBUG>1) {cout<<"kill:    "; seq_best(t)->write(cout); cout<<endl<<"seq_new:    ";write(seq_new);cout<<endl;}
    #endif
    // bis zur ersten NULL evaluieren...
    FOR1D(seq_new, t2) {
      if (seq_new(t2) == NULL)
        break;
    }
    value = inferStateRewards(t2-1);
    double ruleOutcome_rewards_sum = 0.;
    if (use_ruleOutcome_rewards) {
      ruleOutcome_rewards_sum = calcRuleRewards(seq_new);
      value += ruleOutcome_rewards_sum;
    }
    if (DEBUG>1) {cout<<"Shortened value = "<<value<<"  (use_ruleOutcome_rewards="<<use_ruleOutcome_rewards<<", ruleOutcome_rewards_sum="<<ruleOutcome_rewards_sum<<")"<<endl;}
    cout<<std::flush;
    if (value > value_best) {
      if (DEBUG>1) {cout<<"ACCEPTED."<<endl;}
      value_best = value;
      seq_best = seq_new;
    }
    else {
      if (DEBUG>1) {cout<<"Rejected."<<endl;}
      t++;
    }
  }
  if (DEBUG>0) {
    cout<<"Result:"<<endl;
    cout<<"seq_best: ";write(seq_best);cout<<endl;
    if (seq_best.N > 0   &&  seq_best(0) != seq_old(0)) cout<<"  --> Shortened!"<<endl;
    PRINT(value_best);
  }
  if (DEBUG>0) {cout<<"A_PRADA::shorten_plan [END]"<<endl;}
  return value_best;
}


void A_PRADA::plan_full(LitL& best_plan, double& best_value, const SymbolicState& s, uint max_runs) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"A_PRADA::plan_full [START]"<<endl;}
  setStartState(s);
  
  // (1) take into account last best sequence: NOT IMPLEMENTED YET
//   if (last_seq.N > 0) {
//   uint i;
//     LitL last_seq_without_first;
//     FOR1D(last_seq, i) {
//       if (i>0)
//         last_seq_without_first.append(last_seq(i));
//     }
//     last_seq_without_first.append(NULL);
//     LitL last_seq_without_first__full;
//     DEBUG_SINGLE_SAMPLE = 0; // 100
//     infer(last_seq_without_first__full); //  fehlt ne aktion
//     DEBUG_SINGLE_SAMPLE = 0;
//     double v_last_new = inferStateRewards();
//     if (DEBUG>1) {
//       cout<<"Checking last sequence: ";
//       cout<<"Last:  ";write(last_seq);cout<<endl;
//       PRINT(last_value);
//       cout<<"Last minus first:  ";write(last_seq_without_first__full);cout<<endl;
//       PRINT(v_last_new);
//     }
//     double ESTIMATED_NOISE_GAIN = 0.03;
//     // take only if last action has sensible effects
//     // hacked: +1000
//     if (v_last_new >=  last_value / discount + 1000  &&  v_last_new >= horizon * ESTIMATED_NOISE_GAIN) {
//       if (DEBUG>1) {cout<<" --> taken!"<<endl;}
//       last_seq = last_seq_without_first__full;
//       last_value = v_last_new;
//       return last_seq_without_first__full(0);
//     }
//   }
  
  // (2) do usual PRADA planning and thereafter try to shorten returned plan
  uint i;
  LitL plan_short;
  for (i=0; i<3; i++) {
    LitL provisory__best_plan;
    double provisory__best_value;
    plan1_wrapper(provisory__best_plan, provisory__best_value, s, max_runs);
    if (DEBUG>1) {cout<<"best_plan found: "<<provisory__best_plan<<endl;}
    LitL plan_short__local;
    last_value = shorten_plan(plan_short__local, provisory__best_plan, provisory__best_value);
    if (DEBUG>1) {cout<<"plan_short__local (value="<<last_value<<"): "<<plan_short__local<<endl;}
    if (plan_short__local.N > 0 
      && reason::calc_uniqueCoveringRule_groundRules_groundAction(this->ground_rules, s, plan_short__local(0)) == 0) {
      MT_MSG("stupid action which we can't use");
      if (DEBUG>1) {cout<<"stupid action which we can't use"<<endl;}
    }
    else if (plan_short__local.N > 0 
      && reason::calc_uniqueCoveringRule_groundRules_groundAction(this->ground_rules, s, plan_short__local(0)) != 0) {
      plan_short = plan_short__local;
      last_seq = plan_short;
      break;
    }
    if (i==2) {// nur fuer den Abbruch
      plan_short = plan_short__local;
      last_seq = plan_short;
      break;
    }
  }

  best_plan = plan_short;
  best_value = last_value;
  
  if (DEBUG>1) {
    PRINT(best_plan);
    PRINT(best_value);
  }
  
  if (DEBUG>0) {cout<<"A_PRADA::plan_full [END]"<<endl;}
}


void A_PRADA::reset() {
  last_seq.clear();
  last_value = -5000.;
}







/************************************************
 * 
 *     PRADA_DBN
 * 
 ************************************************/

uintA __prada_dbn__constantsIndices(100);
uint __prada_dbn__num_constants = 10000;

inline uint getIndex__PRADA_DBN(const uintA& args) {
  uint args_idx=0;
  uint i;
  FOR1D(args, i) {
    args_idx += ((uint) pow(__prada_dbn__num_constants, i)) * __prada_dbn__constantsIndices(args(i));
  }
//   cout<<"getIndex: constants="<<constants<<"  args="<<args<<"    args_idx="<<args_idx<<endl;
  return args_idx;
}


PRADA_DBN::PRADA_DBN(const uintA& _net_constants, const SymL& _net_symbols_state, const SymL& _net_symbols_action, RuleSet& _ground_rules, double noise_softener, uint horizon) {
  this->net_constants = _net_constants;
  uint i;
  FOR1D(this->net_constants, i) {
    __prada_dbn__constantsIndices(this->net_constants(i)) = i;
  }
  __prada_dbn__num_constants = this->net_constants.N;
  this->net_symbols_state = _net_symbols_state;
  Symbol::sort(this->net_symbols_state);
  this->net_symbols_action = _net_symbols_action;
  this->ground_rules = _ground_rules;
  this->noise_softener = noise_softener;
  this->horizon = horizon;
  create_dbn_structure(this->net_symbols_state, this->net_symbols_action);
  create_dbn_params();
}


PRADA_DBN::~PRADA_DBN() {
  uint i;
  FOR1D(vars_state__prim, i) {
    delete vars_state__prim(i);
  }
  FOR1D(vars_state__derived, i) {
    delete vars_state__derived(i);
  }
  FOR1D(vars_action, i) {
    delete vars_action(i);
  }
}





/************************************************
 * 
 *     PRADA_DBN  -  infer rules
 * 
 ************************************************/


#ifdef FAST
double inferPhi(const Rule& grounded_rule, uint rule_id, uint t, const RVL& vars_context) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"inferPhi [START]"<<endl;}
  if (DEBUG>0) {
    PRINT(t);
    grounded_rule.write();
  }
  double phi = 1.0;
  uint i;
  double prob;
  uint v_current=0;
  FOR1D(grounded_rule.context, i) {
    prob = 0.0;
    if (grounded_rule.context(i)->s->range_type != Symbol::binary) {
      //NIY;

      uint val;
      uint id_2D__val1, id_2D__val2;
      Literal* clit = grounded_rule.context(i);
      // f(x) comp constant
      {
        LiteralRV* var = vars_context(rule_id, v_current++);
        if (var->type == LiteralRV::expectation) {
          if (Literal::compare(var->P(t,0), clit->comparison_type, clit->value))
            prob += 1.;
        }
        else {
          FOR1D(var->range, val) {
            if (Literal::compare(var->range(val), clit->comparison_type, clit->value))
              prob += var->P(t, val);
          }
        }
      }
#if 0
      // f(x) < f(y)
      else {
        LiteralRV* variable1 = vars_context(rule_id, v_current++);
        LiteralRV* variable2 = vars_context(rule_id, v_current++);
        if (DEBUG>2) {
          variable1->write(cout);
          variable2->write(cout);
        }
        if (variable1->type == RV_TYPE__FUNC_EXPECT) {
          if (compare(variable1->P(t,0), variable2->P(t,0), ((ComparisonAtom*)clit->atom)->comparisonType))
            prob +=  1.0;
        }
        else {
          id_2D__val1 = variable1->P.d1 * t; // Wird unten am Schleifenende je 1 hochgezaehlt.
          FOR1D(variable1->range, val) {
            if (((ComparisonAtom*)clit->atom)->comparisonType == comparison_equal) {
              prob += variable1->P.p[id_2D__val1] * variable2->P.p[id_2D__val1];
            }
            else if (((ComparisonAtom*)clit->atom)->comparisonType == comparison_less) {
              id_2D__val2 = variable2->P.d1 * t + val + 1;
              for (val2 = val+1; val2 < variable2->range.d0; val2++) {
                prob += variable1->P.p[id_2D__val1] * variable2->P.p[id_2D__val2];
                id_2D__val2++;
              }
            }
            else if (((ComparisonAtom*)clit->atom)->comparisonType == comparison_lessEqual) {
              id_2D__val2 = variable2->P.d1 * t + val;
              for (val2 = val; val2 < variable2->range.d0; val2++) {
                prob += variable1->P.p[id_2D__val1] * variable2->P.p[id_2D__val2];
                id_2D__val2++;
              }
            }
            else if (((ComparisonAtom*)clit->atom)->comparisonType == comparison_greater) {
              id_2D__val2 = variable2->P.d1 * t + val - 1;
              if (val>0) {
                for (val2 = val-1; val2--; ) {
                  prob += variable1->P.p[id_2D__val1] * variable2->P.p[id_2D__val2];
                  id_2D__val2--;
                }
              }
            }
            else if (((ComparisonAtom*)clit->atom)->comparisonType == comparison_greaterEqual) {
              id_2D__val2 = variable2->P.d1 * t + val;
              for (val2 = val; val2--; ) {
                prob += variable1->P.p[id_2D__val1] * variable2->P.p[id_2D__val2];
                id_2D__val2--;
              }
            }
            
            id_2D__val1++;
          }
//           FOR1D(var2->range, val2) {
// //               if (((ComparisonAtom*)clit->atom)->compare(var1->range(val), var2->range(val2))) prob += var1->P(t,val) * var2->P(t,val2);
// //               if (((ComparisonAtom*)clit->atom)->compare(var1->range(val), var2->range(val2))) prob += var1->P.p[id_2D__val1] * var2->P.p[id_2D__val2];
//             if (((ComparisonAtom*)clit->atom)->compare(var1->range.p[id_1D__val1], var2->range.p[id_1D__val2])) prob += var1->P.p[id_2D__val1] * var2->P.p[id_2D__val2];
//             id_2D__val2++;
//             id_1D__val2++;
//           }
          
        }
      }
#endif
    }
    else {
      // Can be safely assumed to be binary (since function values only appear in p_comp).
      LiteralRV* var = vars_context(rule_id, v_current++);
      if (grounded_rule.context(i)->value)
        prob = var->P(t,1);
      else
        prob = var->P(t,0);
    }
    if (DEBUG>0) {cout<<i<<":"<<prob<<"  ";}
    phi *= prob;
    if (phi < RULE_MIN_PROB) {phi=0.0; break;}
  }
  if (DEBUG == 1 && phi > 0) {
    grounded_rule.write();
    cout<<"  -->  phi="<<phi<<endl;
    if (phi>0.) {grounded_rule.action->write(); cout<<" *****"<<endl;}
    cout<<"inferPhi [END]"<<endl;
  }
  else if (DEBUG > 1) {
    cout<<"  -->  phi="<<phi<<endl;
    if (phi>0.) {grounded_rule.action->write(); cout<<" *****"<<endl;}
    cout<<"inferPhi [END]"<<endl;
  }
  return phi;
}
#endif
#ifndef FAST
double inferPhi(const Rule& grounded_rule, uint t, RV_Manager* rvm) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"inferPhi [START]"<<endl;}
  if (DEBUG>0) {
    PRINT(t);
    grounded_rule.write();
  }
  double phi = 1.0;
  uint i, val, val2;
  double prob;
  FOR1D(grounded_rule.context, i) {
    prob = 0.0;
    if (grounded_rule.context(i)->s->type == Predicate::predicate_comparison) {
      ComparisonLiteral* clit = (ComparisonLiteral*) grounded_rule.context(i);
      // f(x) < c
      if (((ComparisonPredicate*) clit->s)->constantBound) {
        FunctionAtom* fi = rvm->getFVW(((ComparisonAtom*)clit->atom)->f, ((ComparisonAtom*)clit->atom)->args);
        LiteralRV* var = rvm->fi2v(fi);
        if (var->type == RV_TYPE__FUNC_EXPECT) {
          if (((ComparisonAtom*)clit->atom)->compare(var->P(t,0))) prob += 1.;
        }
        else {
          FOR1D(var->range, val) {
            if (((ComparisonAtom*)clit->atom)->compare(var->range(val))) prob += var->P(t,val);
          }
        }
      }
      // f(x) < f(y)
      else {
        uintA sa_1, sa_2;
        FOR1D(((ComparisonAtom*)clit->atom)->args, val) {
          if (val<((ComparisonAtom*)clit->atom)->args.N/2)
            sa_1.append(((ComparisonAtom*)clit->atom)->args(val));
          else
            sa_2.append(((ComparisonAtom*)clit->atom)->args(val));
        }
        FunctionAtom* fvw1 = rvm->getFVW(((ComparisonAtom*)clit->atom)->f, sa_1);
        FunctionAtom* fvw2 = rvm->getFVW(((ComparisonAtom*)clit->atom)->f, sa_2);
        LiteralRV* var1 = rvm->fi2v(fvw1);
        LiteralRV* var2 = rvm->fi2v(fvw2);
        if (var1->type == RV_TYPE__FUNC_EXPECT) {
          if (((ComparisonAtom*)clit->atom)->compare(var1->P(t,0), var2->P(t,0))) prob +=  1.0;
        }
        else {
          FOR1D(var1->range, val) {
            FOR1D(var2->range, val2) {
              if (((ComparisonAtom*)clit->atom)->compare(var1->range(val), var2->range(val2))) prob += var1->P(t,val) * var2->P(t,val2);
            }
          }
        }
      }
    }
    else {
      // Can be safely assumed to be binary (since function values only appear in p_comp).
      if (grounded_rule.context(i)->positive)
        prob = RVefficiency__atom2var(grounded_rule.context(i))->P(t,1);
      else
        prob = RVefficiency__atom2var(grounded_rule.context(i))->P(t,0);
    }
    if (DEBUG>0) {cout<<i<<":"<<prob<<"  ";}
    alpha *= prob;
    if (alpha < RULE_MIN_PROB) {
      alpha=0.0; 
      if (DEBUG>0) {cout<<"Pruning rule to 0."<<endl;}
      break;
    }
  }
  if (DEBUG>0) {cout<<"  -->  phi="<<phi<<endl<<"inferPhi [END]"<<endl;}
  return phi;
}
#endif



// vars_rules_simple:  P(\phi_r | s)
// vars_rules: P(\phi_r | -\phi_r', s)
void PRADA_DBN::inferRules(uint t) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"inferRules [START]"<<endl;}
  if (DEBUG>0) {PRINT(t);}
  uint r, r2;
  
  // vars_rules_simple
  uint p_id = t * vars_rules_simple.d1;
  FOR1D_(ground_rules, r) {
#ifdef FAST
    vars_rules_simple.p[p_id++] = inferPhi(*ground_rules.elem(r), r, t, vars_context);
#endif
#ifndef FAST
    vars_rules_simple(t,r) = inferPhi(*ground_rules.elem(r), t);
#endif
  }
  if (DEBUG>1) {
    cout<<"phi  P(phi_r | s): (for >0)";
    for(r=0; r<ground_rules.num(); r++) {
      if (vars_rules_simple(t,r)>0.0) {
        ground_rules.elem(r)->write();
        cout << " -->  P(\\phi_r | s) = " << vars_rules_simple(t,r) << endl;
      }
    }
    cout<<endl;
  }
  
  // vars_rules
  p_id = t * vars_rules_simple.d1;
  FOR1D_(ground_rules, r) {
#ifdef FAST
    if (r>0) p_id++;
//     PRINT(p_id);
//     PRINT(vars_rules.d1 * t + r);
    vars_rules.p[p_id] = vars_rules_simple.p[p_id];
    if (vars_rules.p[p_id] < RULE_MIN_PROB) {vars_rules.p[p_id]=0.0; continue;}
    uint idx_symbol = net_symbols_action.findValue(ground_rules.elem(r)->action->s);
    uint idx_args = getIndex__PRADA_DBN(ground_rules.elem(r)->action->args);
    for (r2=0; r2<action2rules_num(idx_symbol, idx_args); r2++) {
      if (action2rules(idx_symbol, idx_args, r2) == r) continue;
      else {
        vars_rules.p[p_id] *= (1. - vars_rules_simple(t,action2rules(idx_symbol, idx_args, r2)));
        if (vars_rules.p[p_id] < RULE_MIN_PROB) {vars_rules.p[p_id]=0.0; break;}
      }
    }
#endif
#ifndef FAST
    vars_rules(t,r) = alpha(t,r);
    if (vars_rules(t,r) < RULE_MIN_PROB) {vars_rules(t,r)=0.0; continue;}
    uint idx_symbol = logicObjectManager::p_actions.findValue(ground_rules.elem(r)->action->s);
    uint idx_args = getIndex(constants, ground_rules.elem(r)->action->args);
    for (r2=0; r2<action2rules_num(idx_symbol, idx_args); r2++) {
      if (action2rules(idx_symbol, idx_args,r2) == r) continue;
      else {
        vars_rules(t,r) *= (1. - alpha(t,action2rules(idx_symbol, idx_args,r2)));
        if (vars_rules(t,r) < RULE_MIN_PROB) {vars_rules(t,r)=0.0; break;}
      }
    }
#endif
  }
  if (DEBUG>0) {
    cout<<"vars_rules  P(phi_r | s, phi_r'):";
    for(r=0; r<ground_rules.num(); r++) {
      cout<<"     "<<*ground_rules.elem(r)->action<<" "<<vars_rules(t,r)<<" (" << vars_rules_simple(t,r) << ")";
    }
    cout<<endl;
    if (DEBUG>1) {
      cout<<"Rules with prob>0.0:"<<endl;
      for(r=0; r<ground_rules.num(); r++) {
        if (vars_rules(t,r)>0.0) {
          ground_rules.elem(r)->write();
          cout << " --> P(phi_r | s, phi_r') = " << vars_rules(t,r) << endl;
        }
      }
      cout<<endl;
    }
  }
  
  if (DEBUG>0) {cout<<"inferRules [END]"<<endl;}
}



/************************************************
 * 
 *     PRADA_DBN  -  infer states
 * 
 ************************************************/

void PRADA_DBN::inferState(uint t, Literal* given_action) {
  NIY;
//   uint idx_symbol = logicObjectManager::p_actions.findValue(given_action->s);
//   uint idx_args = getIndex(reason::getConstants(), given_action->args);
//   double action_weight = 0.;
//   uint r;
//   for (r=0; r<action2rules_num(idx_symbol, idx_args); r++) {
//     action_weight += vars_rules(t-1, action2rules(idx_symbol, idx_args, r));
//   }
//   inferState(t, given_action, action_weight);
}


// from state and action at t-1
void PRADA_DBN::inferState(uint t, Literal* given_action, double given_action_weight) {
  uint DEBUG = 0;
  
  if (DEBUG>0) {cout<<"inferState [START]"<<endl;}
  CHECK(t>0, "only works for state > 0");
  if (DEBUG>0) {
    PRINT(t);
    given_action->write(cout); cout<<endl;
    PRINT(given_action_weight);
  }
   
  uint r, r2, v, val;
  
  // Determine which action has been taken
  uint action_idx_symbol = net_symbols_action.findValue(given_action->s);
  uint action_idx_args = getIndex__PRADA_DBN(given_action->args);
//   PRINT(given_action->args);
//   PRINT(action_idx_symbol);
//   PRINT(action_idx_args);
  
  uintA rules_with_sampled_action;
  uint p_id_start = (action_idx_symbol * action2rules.d1  +  action_idx_args) * action2rules.d2;
  for(r=0; r<action2rules_num(action_idx_symbol, action_idx_args); r++) {
    rules_with_sampled_action.append(action2rules.p[p_id_start+r]);
  }
  if (DEBUG>1) {PRINT(rules_with_sampled_action);}
  if (DEBUG>1) {
    cout<<"Important rules: "<<endl;
    FOR1D(rules_with_sampled_action, r) {
      ground_rules.elem(rules_with_sampled_action(r))->write(cout);
      cout<<"  --> applicability "<<vars_rules(t-1, rules_with_sampled_action(r));
      if (vars_rules(t-1, rules_with_sampled_action(r)) > 0.05) {
        cout<<" ***";
      }
      cout<<endl;
    }
  }
  uint p_id, p_id2;  
  // (1) calc P_t+1(v_i|r, a_t, s)
  MT::Array< arr* > p_P_v__r_val;  // dim 1: variable;  dim 2: arr for rule x value
  FOR1D(vars_state__prim, v) {
    LiteralRV* var = vars_state__prim(v);
    arr* p_local_P_r_val = new arr(ground_rules.num(), var->dim);
    arr& reference__impacts_val_v = impacts_val(v);
    if (DEBUG>1) {PRINT(rules_with_sampled_action);}
    FOR1D(rules_with_sampled_action, r) {
      r2 = rules_with_sampled_action(r);
      double one_minus_impacts_V__r2_v = 1. - impacts_V(v)(r2);
      p_id = r2 * p_local_P_r_val->d1;
      p_id2 = var->P.d1 * (t-1);
      FOR1D(var->range, val) {
        p_local_P_r_val->p[p_id] = reference__impacts_val_v.p[p_id] + one_minus_impacts_V__r2_v * var->P.p[p_id2];
        p_id++;
        p_id2++;
      }
    }
    p_P_v__r_val.append(p_local_P_r_val);
  }
  if (DEBUG>2) {
    cout<<"P_v_r =  P(v_i | r, a_t, s):"<<endl;
    for(r=0; r<ground_rules.num(); r++) {
      cout<<"rule "<<r;
      if (rules_with_sampled_action.findValue(r) < 0) {
        cout << "  not important --> ignored"<<endl;  continue;
      }
      cout<<endl;
      for(v=0; v<vars_state__prim.N; v++) {
        cout<<*vars_state__prim(v)->lit<<":  ";
        FOR1D(vars_state__prim(v)->range, val) {
          cout<<vars_state__prim(v)->range(val)<<":"<<(*p_P_v__r_val(v))(r,val);
          cout<<"  ";
        }
        cout<<endl;
      }
    }
  }
  // (2) calc P_t+1(v_i|a,s) = next time-step distribution
#ifdef FAST
  FOR1D(vars_state__prim, v) {
    LiteralRV* var = vars_state__prim(v);
    FOR1D(var->range, val) {
      var->P(t,val) = (1.-given_action_weight) * var->P(t-1,val);
    }
  }
  for(r=0; r<action2rules_num(action_idx_symbol, action_idx_args); r++) {
    r2 = action2rules(action_idx_symbol, action_idx_args, r);
    double r2_prob = vars_rules(t-1, r2);
    FOR1D(vars_state__prim, v) {
      arr* p_P_r_val = p_P_v__r_val(v);
      LiteralRV* var = vars_state__prim(v);
      p_id = r2 * p_P_r_val->d1;
      p_id2 = t * var->P.d1;
      FOR1D(var->range, val) {
//         var->P(t,val) += p_P_r_val->p[p_id] * r2_prob;
        var->P.p[p_id2] += p_P_r_val->p[p_id] * r2_prob;
        p_id++;
        p_id2++;
      }
    }
  }
#else
  FOR1D(vars_state__prim, v) {
    arr& reference_P_v_r__v = P_v_r(v);
    LiteralRV* var = vars_state__prim(v);
    FOR1D(var->range, val) {
      var->P(t,val) = 0.0;
      for(r=0; r<action2rules_num(action_idx_symbol, action_idx_args); r++) {
        r2 = action2rules(action_idx_symbol, action_idx_args, r);
        var->P(t,val) += reference_P_v_r__v(r2,val) * vars_rules(t-1, r2);
      }
      var->P(t,val) += (1.-given_action_weight) * var->P(t-1,val);
    }
  }
#endif
  listDelete(p_P_v__r_val);
    
  // INFER DERIVED SYMBOLS
  if (DEBUG>1) {cout<<"After inference over primitive symbols:" << endl; writeState(t);}
  calcDerived(t);
  checkStateSoundness(t);


  if (DEBUG>1) {writeState(t);}
  if (DEBUG>0) {cout<<"inferState [END]"<<endl;}
}


void PRADA_DBN::inferStates(const LitL& given_action_sequence) {
//   CHECK(given_action_sequence.N == horizon, "");
  uint t;
  FOR1D(given_action_sequence, t) {
    inferRules(t);
    inferState(t+1, given_action_sequence(t));
  }
}

  
  
/************************************************
 * 
 *     PRADA_DBN  -  Set evidence
 * 
 ************************************************/
  

void PRADA_DBN::setAction(uint t, Literal* action) {
  LiteralRV* action_rv = RVefficiency__atom2var(Literal::get(action->s, action->args, 1.));
  uint i;
  FOR1D(vars_action, i) {
    if (action_rv == vars_action(i)) {
      vars_action(i)->P(t, 0) = 0.0;
      vars_action(i)->P(t, 1) = 1.0;
    }
    else {
      vars_action(i)->P(t, 0) = 1.0;
      vars_action(i)->P(t, 1) = 0.0;
    }
  }
}


void PRADA_DBN::setState(const LitL& lits, uint t) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"PRADA_DBN::setState [START]"<<endl;}
  uint v, val;
  if (DEBUG>0) {cout<<"SymbolicState: "<<lits<<endl;  PRINT(t);}
  // Everything not specified is set to false!
  FOR1D(vars_state__prim, v) {
    vars_state__prim(v)->P(t,0) = 1.0;
    for (val=1; val<vars_state__prim(v)->dim; val++) {
      vars_state__prim(v)->P(t,val) = 0.0;
    }
  }
  // Incorporating evidence
  FOR1D(lits, v) {
    if (DEBUG>0) {cout<<"setState:  "<<*lits(v)<<endl;}
    LiteralRV* rv = RVefficiency__atom2var(lits(v));
    if (rv == NULL) {
      if (DEBUG>0) {cout<<"no random variable --> omitted"<<endl;}
//       MT_MSG("Omitting setting state lit "<<*lits(v)<<" in DBN^(t=0).");
      continue;
    }
    if (rv->type == LiteralRV::simple) {
      for (val=0; val<rv->dim; val++) {
        if (TL::areEqual(lits(v)->value, rv->range(val))) {
          rv->P(t,val) = 1.0;
        }
        else
          rv->P(t,val) = 0.0;
      }
    }
    else {
      rv->P(t,0) = lits(v)->value;
    }
    if (DEBUG>0) {cout<<" --> "; rv->write(); cout<<endl;}
  }
  
  calcDerived(t);
  if (DEBUG>1){writeState(0);}
  checkStateSoundness(t);
  if (DEBUG>0) {cout<<"PRADA_DBN::setState [END]"<<endl;}
}


void PRADA_DBN::setStateUniform(uint t) {
  uint v, val;
  FOR1D(vars_state__prim, v) {
    FOR1D(vars_state__prim(v)->range, val) {
      vars_state__prim(v)->P(t,val) = 1. / (1. * vars_state__prim(v)->range.N);
    }
  }
  
  checkStateSoundness(t, true);
}

  
  
  
/************************************************
 * 
 *     PRADA_DBN  -  Comparing probabilities
 * 
 ************************************************/

double PRADA_DBN::log_probability(uint t, const SymbolicState& state) const {
  uint DEBUG = 0;
  if (DEBUG>0) {
    cout<<"log_probability [START]"<<endl;
    PRINT(t);
    state.write(cout); cout<<endl;
  }
  uint i;
  double total_log_prob = 0.;
  double prob_variable = 0.;
  FOR1D(state.lits, i) {
    if (state.lits(i)->s->symbol_type != Symbol::primitive)
      continue;
    LiteralRV* var = RVefficiency__atom2var(state.lits(i));
    int val_idx = var->range.findValue(state.lits(i)->value);
    CHECK(val_idx>=0, "value not found");
    prob_variable = var->P(t, val_idx);
    if (DEBUG>1) {
      MT::String name;  name << var->lit;
      printf("%-14s",(char*)name);
      cout << "  has prob=" << prob_variable << "   log-prob=" << log(prob_variable);
    }
    total_log_prob += log(prob_variable);
    if (DEBUG>1) {
      cout << "  -->  total_log_prob=" << total_log_prob << endl;
    }
  }
  
  if (DEBUG>0) {
    PRINT(total_log_prob);
    cout<<"log_probability [END]"<<endl;
  }
  
  return total_log_prob; 
}


double PRADA_DBN::belief_difference(uint t, const arr& probs_p_prim, const arr& probs_f_prim) const {
  uint DEBUG = 0;
  uint i, val;
  double total_diff = 0.;
  double local_diff;
  if (DEBUG>0) {cout<<"+++ belief_difference at "<<t<<endl;}
  FOR1D(vars_state__prim, i) {
    local_diff = 0.;
    if (DEBUG>0) {vars_state__prim(i)->lit->write(cout);cout<<": ";}
    FOR1D(vars_state__prim(i)->range, val) {
      local_diff += fabs(probs_p_prim(i,val) - vars_state__prim(i)->P(t,val));
      if (DEBUG>0) {cout << " (" << probs_p_prim(i,val) << " - " << vars_state__prim(i)->P(t,val)<<")";}
    }
    if (DEBUG>0) {cout<<"   --> "<<local_diff<<endl;}
    local_diff /= 1.0 * vars_state__prim(i)->range.N;
    total_diff += local_diff;
  }
//   FOR1D(vars_state__prim, i) {
//     local_diff = 0.;
//     if (DEBUG>0) {vars_state__prim(i)->fi->write(cout);cout<<": ";}
//     FOR1D(vars_state__prim(i)->range, val) {
//       local_diff += fabs(probs_f_prim(i,val) - vars_state__prim(i)->P(t,val));
//       if (DEBUG>0) {cout << " (" << probs_f_prim(i,val) << " - " << vars_state__prim(i)->P(t,val)<<")";}
//     }
//     if (DEBUG>0) {cout<<"   --> "<<local_diff<<endl;}
//     local_diff /= 1.0 * vars_state__prim(i)->range.N;
//     total_diff += local_diff;
//   }
  if (DEBUG>0) {PRINT(total_diff);}
  return total_diff;
}




/************************************************
 * 
 *     PRADA_DBN  -  calcDerived
 * 
 ************************************************/

void calcDerived1(ConjunctionSymbol* s, uint t, const uintA& constants, PRADA_DBN* dbn);
void calcDerived1(TransClosureSymbol* s, uint t, const uintA& constants, PRADA_DBN* dbn);
void calcDerived1(CountSymbol* s, uint t, const uintA& constants, PRADA_DBN* dbn);
void calcDerived1(SumFunction* f, uint t, const uintA& constants, PRADA_DBN* dbn);
void calcDerived1(RewardFunction* f, uint t, const uintA& constants, PRADA_DBN* dbn);
void calcDerived1(DifferenceFunction* s, uint t, const uintA& constants, PRADA_DBN* dbn);


void PRADA_DBN::calcDerived(uint t) {
  uint i;
  FOR1D(net_symbols_state, i) {
    if (net_symbols_state(i)->symbol_type == Symbol::primitive) continue;
    else if (net_symbols_state(i)->symbol_type == Symbol::conjunction)
      calcDerived1((ConjunctionSymbol*) net_symbols_state(i), t, net_constants, this);
    else if (net_symbols_state(i)->symbol_type == Symbol::transclosure)
      calcDerived1((TransClosureSymbol*) net_symbols_state(i), t, net_constants, this);
    else if (net_symbols_state(i)->symbol_type == Symbol::count)
      calcDerived1((CountSymbol*) net_symbols_state(i), t, net_constants, this);
    else if (net_symbols_state(i)->symbol_type == Symbol::sum)
      calcDerived1((SumFunction*) net_symbols_state(i), t, net_constants, this);
    else if (net_symbols_state(i)->symbol_type == Symbol::function_reward)
      calcDerived1((RewardFunction*) net_symbols_state(i), t, net_constants, this);
    else if (net_symbols_state(i)->symbol_type == Symbol::function_difference)
      calcDerived1((DifferenceFunction*) net_symbols_state(i), t, net_constants, this);
    else
      HALT("unknown symbol type: "<<*net_symbols_state(i))
  }
}


// With free vars!
void calcDerived1(ConjunctionSymbol* s, uint t, const uintA& constants, PRADA_DBN* dbn) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcDerived [START]"<<endl;}
  if (DEBUG>0) {cout<<*s<<endl;}
  uintA freeVars;
  s->getFreeVars(freeVars);
  uint i;
  MT::Array< uintA > combos_args__conjunction_lit;
  TL::allPermutations(combos_args__conjunction_lit, constants, s->arity, true, true);
  uint c1, c2;
  double prob;  
  uintA sa(s->arity);
  FOR1D(combos_args__conjunction_lit, c1) {
    Literal* target = Literal::get(s, combos_args__conjunction_lit(c1), 1.);
    if (DEBUG>1) {PRINT(*target);}
    // Free Vars EXISTENTIAL
    // P(p) = 1 - PRODUCT[over free-var combos c](1 - P(basePTs[sub=c]))
    // Intuition: Predicate true if not all base-pred combinations are false.
    if (!s->free_vars_all_quantified) {
      MT::Array< uintA > combos_freevars___base_lits;
      uintA constants_freevars = constants;
      setMinus(constants_freevars, combos_args__conjunction_lit(c1));
      TL::allPermutations(combos_freevars___base_lits, constants_freevars, freeVars.N, false, true);
      arr probs__combos_freevars___base_lits(combos_freevars___base_lits.N);
      FOR1D(combos_freevars___base_lits, c2) {
        Substitution sub;
        for(i=0;i<s->arity;i++) {
          sub.addSubs(i, combos_args__conjunction_lit(c1)(i));
        }
        FOR1D(freeVars, i) {
          sub.addSubs(freeVars(i), combos_freevars___base_lits(c2)(i));
        }
        prob = 1.0;
        FOR1D(s->base_literals, i) {
          Literal* base_lit_ground = sub.apply(s->base_literals(i));
          LiteralRV* var = dbn->RVefficiency__atom2var(base_lit_ground);
          int val_idx = var->range.findValue(base_lit_ground->value);
          prob *= var->P(t,val_idx);
          if (prob < 0.01) {
            if (DEBUG>0) {cout<<"Stop calculation as probablity is already very small - 1: "<<prob<<endl;}
            break;
          }
        }
        probs__combos_freevars___base_lits(c2) = prob;
      }
      prob = 1.0;
      FOR1D(probs__combos_freevars___base_lits, c2) {
        prob *= 1-probs__combos_freevars___base_lits(c2);
      }
      prob = 1.-prob;
    }
    // Free Vars ALL
    // Intuition: Predicate true if all base-pred combinations are true.
    else {
      Substitution sub;
      for(i=0;i<s->arity;i++) {
        sub.addSubs(i, combos_args__conjunction_lit(c1)(i));
      }
      // base literals WITHOUT additional variables
      prob = 1.0;
      FOR1D(s->base_literals, i) {
        if (!s->redundant__base_literal_with_free_vars(i)) {
          Literal* base_lit_ground = sub.apply(s->base_literals(i));
          LiteralRV* var = dbn->RVefficiency__atom2var(base_lit_ground);
          CHECK(var!=NULL, "variable is missing for base_lit_ground "<<*base_lit_ground);
          int val_idx = var->range.findValue(base_lit_ground->value);
          prob *= var->P(t,val_idx);
          if (prob < 0.01)
            break;
        }
      }
      if (prob > 0.99) {
        // base literals WITH additional variables
        uintA constants_freevars = constants;
        setMinus(constants_freevars, combos_args__conjunction_lit(c1));
        MT::Array< uintA > combos_freevars___base_lits;
        TL::allPermutations(combos_freevars___base_lits, constants_freevars, freeVars.N, false, true);
        Substitution sub2;
        sub2 = sub;
        FOR1D(combos_freevars___base_lits, c2) {
          FOR1D(freeVars, i) {
            sub2.addSubs(freeVars(i), combos_freevars___base_lits(c2)(i));
          }
          FOR1D(s->base_literals, i) {
            if (s->redundant__base_literal_with_free_vars(i)) {
              Literal* base_lit_ground = sub2.apply(s->base_literals(i));
              LiteralRV* var = dbn->RVefficiency__atom2var(base_lit_ground);
              int val_idx = var->range.findValue(base_lit_ground->value);
              prob *= var->P(t,val_idx);
              if (prob < 0.01)
                break;
            }
          }
          if (prob < 0.01) break;
        }
      }
      else {
        if (DEBUG>0) {cout<<"Stop calculation as probablity is already very small - 2: "<<prob<<endl;}
      }
    }
    dbn->RVefficiency__atom2var(target)->P(t,0) = 1.-prob;
    dbn->RVefficiency__atom2var(target)->P(t,1) = prob;
    if (DEBUG>0) {cout<<" ---> "<<prob<<endl;}
  }
  if (DEBUG>0) {cout<<"calcDerived [END]"<<endl;}
}



uintA __auxiliary_transclosure__constants;
uintA __auxiliary_transclosure__constantsIndices(100);

// FAST version
void calcDerived_tcp_dfs(arr& probs, const uintA& remaining_constants, const arr& probs_last, double prev_prob, uint prev_constant, const uintA& constants, const arr& PROBS_BASE_TABLE) {
  uint DEBUG = 0;
  uint i;
  double prob_final, prob_single;
  int idx1, idx2;
  idx1 = __auxiliary_transclosure__constantsIndices(prev_constant);
  FOR1D(remaining_constants, i) {
    idx2 = __auxiliary_transclosure__constantsIndices(remaining_constants(i));
    prob_single = PROBS_BASE_TABLE(idx1, idx2);
    if (DEBUG>0) {PRINT(prev_constant);PRINT(remaining_constants(i));PRINT(prob_single);PRINT(prev_prob);PRINT(probs_last(remaining_constants(i)));}
    prob_final = prob_single * prev_prob * probs_last(remaining_constants(i));
    probs.append(prob_final);
    double prev_prob_new = prob_single * prev_prob;
    #define TRANS_CLOSURE_STOP_PROB 0.01
    if (prev_prob_new > TRANS_CLOSURE_STOP_PROB) {
      uintA remaining_constants2 = remaining_constants;
      remaining_constants2.removeValueInSorted(remaining_constants(i), TL::uint_compare);
      calcDerived_tcp_dfs(probs, remaining_constants2, probs_last, prev_prob_new, remaining_constants(i), constants, PROBS_BASE_TABLE);
    }
    else {
//       if (DEBUG>0) {cerr<<"  give-up "<<(constants.N - remaining_constants.N)/*<<endl*/;}
    }
  }
}


MT::Array< TransClosureSymbol* > __auxiliary_transclosure__symbols;
MT::Array< LitL > __auxiliary_transclosure__baseLits;
MT::Array< LitL > __auxiliary_transclosure__transclosureLits;
MT::Array< RVL > __auxiliary_transclosure__transclosureVars;

void calcDerived1(TransClosureSymbol* s, uint t, const uintA& constants, PRADA_DBN* dbn) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcDerived - TransClosureSymbol [START]"<<endl;}
  if (DEBUG>0) {
    s->write(cout);cout<<endl;
    PRINT(constants);
  }
  CHECK(constants.isSorted(TL::uint_compare), "");
  uint i, k;
  CHECK(s->arity==2, "TransClosureSymbol has to be 2dim");
  if (__auxiliary_transclosure__constants != constants) {
    __auxiliary_transclosure__symbols.clear();
    __auxiliary_transclosure__baseLits.clear();
    __auxiliary_transclosure__transclosureVars.clear();
    __auxiliary_transclosure__constants = constants;
    FOR1D(constants, i) {
      __auxiliary_transclosure__constantsIndices(constants(i)) = i;
    }
  }
  LitL baseLits;
  LitL transclosureLits;
  RVL transclosureVars;
  FOR1D(__auxiliary_transclosure__symbols, i) {
    if (__auxiliary_transclosure__symbols(i) == s) {
      baseLits = __auxiliary_transclosure__baseLits(i);
      transclosureLits = __auxiliary_transclosure__transclosureLits(i);
      transclosureVars = __auxiliary_transclosure__transclosureVars(i);
      break;
    }
  }
  if (i==__auxiliary_transclosure__symbols.N) {
    __auxiliary_transclosure__symbols.append(s);
    Literal::getLiterals(baseLits, s->base_symbol, constants, 1.);
    __auxiliary_transclosure__baseLits.append(baseLits);
    Literal::getLiterals(transclosureLits, s, constants, 1.);
    __auxiliary_transclosure__transclosureLits.append(transclosureLits);
    FOR1D(transclosureLits, i) {
      transclosureVars.append(dbn->RVefficiency__atom2var(transclosureLits(i)));
    }
    __auxiliary_transclosure__transclosureVars.append(transclosureVars);
  }
  // precompute all probs
  arr PROBS_BASE_TABLE(constants.N, constants.N);
  PROBS_BASE_TABLE.setUni(0.);
  int idx1, idx2;
  FOR1D(baseLits, i) {
    idx1 = __auxiliary_transclosure__constantsIndices(baseLits(i)->args(0));
    idx2 = __auxiliary_transclosure__constantsIndices(baseLits(i)->args(1));
    PROBS_BASE_TABLE(idx1, idx2) = dbn->RVefficiency__atom2var(baseLits(i))->P(t,1);
  }
  if (DEBUG>0) {PRINT(PROBS_BASE_TABLE);}
  double prob;
  FOR1D(transclosureLits, i) {
    if (DEBUG>0) {cout<<"-- "<<*transclosureLits(i)<<endl;}
    // A-B
    idx1 = __auxiliary_transclosure__constantsIndices(baseLits(i)->args(0));
    idx2 = __auxiliary_transclosure__constantsIndices(baseLits(i)->args(1));
    arr probs(1);
    probs(0) = PROBS_BASE_TABLE(idx1, idx2);
    if (DEBUG>0) {cout<<*baseLits(i)<<" has start prob = "<<probs.last()<<endl;}
    if (idx1!=idx2) {
      // A-...X...-B
      uintA constants_remaining(constants.N-2);
      int __l = 0;
      FOR1D(constants, k) {
        uint idx_constant = __auxiliary_transclosure__constantsIndices(constants(k));
        if (idx_constant != idx1  &&  idx_constant != idx2)
          constants_remaining(__l++) = constants(k);
      }
      if (DEBUG>0) {PRINT(constants_remaining);}
      arr probs_last(constants_remaining.max()+1);
      uint num_nonzero_probs__arg1 = 0, num_nonzero_probs__arg2 = 0;
      FOR1D(constants_remaining, k) {
        uint idx_remaining = __auxiliary_transclosure__constantsIndices(constants_remaining(k));
        double prob = PROBS_BASE_TABLE(idx_remaining, idx2);
        probs_last(constants_remaining(k)) = prob;
        if (DEBUG>1) {cout<<"prob for ("<<constants(idx_remaining)<<","<<constants(idx2)<<") = "<<prob<<endl;}
        if (prob>0.01) num_nonzero_probs__arg2++;
        prob = PROBS_BASE_TABLE(idx1, idx_remaining);
        if (DEBUG>1) {cout<<"prob for ("<<constants(idx1)<<","<<constants(idx_remaining)<<") = "<<prob<<endl;}
        if (prob > 0.01) num_nonzero_probs__arg1++;
      }
      if (DEBUG>0) {PRINT(num_nonzero_probs__arg1);  PRINT(num_nonzero_probs__arg2);}
      if (num_nonzero_probs__arg1>0  &&  num_nonzero_probs__arg2>0) {
        if (DEBUG>0) {cout<<"starting depth search"<<endl;}
        calcDerived_tcp_dfs(probs, constants_remaining, probs_last, 1.0, baseLits(i)->args(0), constants, PROBS_BASE_TABLE);
      }
    }
    // true if not all possible chains are false
    prob = 1.;
    FOR1D(probs, k) {
      prob *= (1-probs(k));
    }
    prob = 1 - prob;
    if (DEBUG>0) {
      PRINT(probs);
      PRINT(prob);
    }
    LiteralRV* var = transclosureVars(i);
    var->P(t,0) = 1.-prob;
    var->P(t,1) = prob;
  }
  if (DEBUG>0) {cout<<"calcDerived - TransClosureSymbol [END]"<<endl;}
}


void calcDerived1(CountSymbol* s, uint t, const uintA& constants, PRADA_DBN* dbn) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcDerived - CountSymbol [START]"<<endl;}
  if (DEBUG>0) {
    s->write(cout);cout<<endl;
    PRINT(t);
  }
  uint i, c, c2;
  MT::Array< uintA > combos;
  TL::allPermutations(combos, constants, s->arity, true, true);
  uintA args_base(s->base_literal->s->arity);
  double expect;
  uintA freeVars;
  FOR1D(s->base_literal->args, i) {
    if (s->base_literal->args(i) >= s->arity)
      freeVars.setAppend(s->base_literal->args(i));
  }
  if (DEBUG>0) {PRINT(freeVars);}
  FOR1D(combos, c) {
    if (DEBUG>0) {cout<<"-- "<<combos(c)<<endl;}
    Literal* lit = Literal::get(s, combos(c), 1.);
    LiteralRV* var = dbn->RVefficiency__atom2var(lit);
    FOR1D(s->base_literal->args, i) {
      if (s->base_literal->args(i) < s->arity) {
        args_base(i) = combos(c)(s->base_literal->args(i));
      }
    }
    uintA local_constants = constants;
    setMinus(local_constants, combos(c));
    MT::Array< uintA > inner_combos;
    TL::allPermutations(inner_combos, local_constants, freeVars.N, true, true);
    // (1) Expectation
    if (var->type == LiteralRV::expectation) {
      expect = 0.0;
      FOR1D(inner_combos, c2) {
        FOR1D(s->base_literal->args, i) {
          if (s->base_literal->args(i) >= s->arity) {
            args_base(i) = inner_combos(c2)(s->base_literal->args(i) - s->arity);
          }
        }
        Literal* lit_base_grounded = Literal::get(s->base_literal->s, args_base, 1.);
        LiteralRV* var__lit_base = dbn->RVefficiency__atom2var(lit_base_grounded);
        CHECK(var__lit_base != NULL, "var__lit_base is missing!");
        FOR1D(var__lit_base->range, i) {
          expect += var__lit_base->range(i) * var__lit_base->P(t,i);
        }
        if (DEBUG>2) {
          cout<<"inner combo: "<<inner_combos(c2)<<"  "<<*lit_base_grounded
              <<" "<<var__lit_base->P(t,1)<<endl;
        }
      }
      if (DEBUG>1) {PRINT(expect);}
      var->P(t,0) = expect;
    }
    // or (2) True Distribution
    else {
      NIY;
    }
    if (DEBUG>0) {cout<<"Updated variable: ";var->write();}
  }
  if (DEBUG>0) {cout<<"calcDerived - CountSymbol [END]"<<endl;}
}


void calcDerived1(SumFunction* s, uint t, const uintA& constants, PRADA_DBN* dbn) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcDerived - SumFunction [START]"<<endl;}
  if (DEBUG>0) {
    s->write(cout);cout<<endl;
    PRINT(t);
  }
  CHECK(s->arity == 0, "so far implemented only for zero-ary");
  uintA empty;
  Literal* lit = Literal::get(s, empty, 1);
  if (DEBUG>0) {PRINT(*lit);}
  LiteralRV* var = dbn->RVefficiency__atom2var(lit);
  if (DEBUG>0) {var->write();}
  CHECK(var->type == LiteralRV::expectation, "only defined for expectation random vars");
  uint c;
  MT::Array< uintA > combos;
  TL::allPermutations(combos, constants, s->base_symbol->arity, true, true);
  double sum=0.;
  FOR1D(combos, c) {
    if (DEBUG>1) {cout<<"-- "<<combos(c)<<endl;}
    Literal* lit_base = Literal::get(s->base_symbol, combos(c), 1.);
    if (DEBUG>1) {cout<<*lit_base;}
    LiteralRV* var_base = dbn->RVefficiency__atom2var(lit_base);
    if (DEBUG>1) {cout<<"  "<<var_base->P(t,0)<<endl;}
    CHECK(var_base->type == LiteralRV::expectation, "only defined for expectation random vars");
    sum += var_base->P(t,0);
  }
  var->P(t,0) = sum;
  if (DEBUG>0) {var->write();}
  if (DEBUG>0) {cout<<"calcDerived - SumFunction [END]"<<endl;}
}


void calcDerived1(RewardFunction* f, uint t, const uintA& constants, PRADA_DBN* dbn) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcDerived - RewardFunction [START]"<<endl;}
  if (DEBUG>0) {
    f->write(cout);cout<<endl;
    PRINT(t);
  }
  CHECK(f->arity == 0, "so far implemented only for zero-ary");
  uintA empty;
  NIY;
#if 0
  FunctionAtom* fi = rvm->getFVW(f, empty);
  LiteralRV* var = rvm->fi2v(fi);
  double prob = 1.0;
  uint i;
  FOR1D(f->grounded_pis, i) {
    PredicateRV* sub_var = RVefficiency__atom2var(f->grounded_pis(i));
    if (f->grounded_pis(i)->positive)
      prob *= sub_var->P(t,1);
    else
      prob *= sub_var->P(t,0);
  }
  var->P(t,0) = prob * f->reward_value;
  if (DEBUG>0) {PRINT(prob);  var->write();}
  if (DEBUG>0) {cout<<"calcDerived - RewardFunction [END]"<<endl;}
#endif
}


void calcDerived1(DifferenceFunction* s, uint t, const uintA& constants, PRADA_DBN* dbn) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calcDerived - DifferenceFunction [START]"<<endl;}
  if (DEBUG>0) {
    s->write(cout);cout<<endl;
    PRINT(t);
  }

  uintA freeVars;
  s->getFreeVars(freeVars);
  uint i;
  MT::Array< uintA > combos_args__conjunction_lit;
  TL::allPermutations(combos_args__conjunction_lit, constants, s->arity, true, true);
  uint c1, c2;
  double prob;  
  uintA sa(s->arity);
  FOR1D(combos_args__conjunction_lit, c1) {
    Literal* target = Literal::get(s, combos_args__conjunction_lit(c1), 1.);
    if (DEBUG>1) {PRINT(*target);}
    // Free Vars EXISTENTIAL
    // P(p) = 1 - PRODUCT[over free-var combos c](1 - P(basePTs[sub=c]))
    // Intuition: Predicate true if not all base-pred combinations are false.
    MT::Array< uintA > combos_freevars___base_lits;
    uintA constants_freevars = constants;
    setMinus(constants_freevars, combos_args__conjunction_lit(c1));
    TL::allPermutations(combos_freevars___base_lits, constants_freevars, freeVars.N, false, true);
    arr probs__combos_freevars___base_lits(combos_freevars___base_lits.N);
    FOR1D(combos_freevars___base_lits, c2) {
      Substitution sub;
      for(i=0;i<s->arity;i++) {
        sub.addSubs(i, combos_args__conjunction_lit(c1)(i));
      }
      FOR1D(freeVars, i) {
        sub.addSubs(freeVars(i), combos_freevars___base_lits(c2)(i));
      }
      prob = 1.0;
      FOR1D(s->restrictionLits, i) {
        Literal* base_lit_ground = sub.apply(s->restrictionLits(i));
        LiteralRV* var = dbn->RVefficiency__atom2var(base_lit_ground);
        int val_idx = var->range.findValue(base_lit_ground->value);
        prob *= var->P(t,val_idx);
        if (prob < 0.01) {
          if (DEBUG>1) {cout<<"Stop calculation as probablity is already very small - 1: "<<prob<<endl;}
          break;
        }
      }
      if (prob > 0.99) {
        LiteralRV *rv1 = dbn->RVefficiency__atom2var(sub.apply(s->baseFunctionLit1));
        LiteralRV *rv2 = dbn->RVefficiency__atom2var(sub.apply(s->baseFunctionLit2));
        CHECK(rv1->type == LiteralRV::expectation && rv2->type == LiteralRV::expectation, "Both functions of a function difference symbol must be of type expectation!")
        if (DEBUG > 0) cout << *target << " expectation: " << rv1->P(t,0) - rv2->P(t,0) << endl;
        dbn->RVefficiency__atom2var(target)->P(t,0) = rv1->P(t,0) - rv2->P(t,0);
        return;
      }
    }
  }
  CHECK(false, "Restriction literals do not cover!")
}






/************************************************
 * 
 *     PRADA_DBN  -  Misc
 * 
 ************************************************/

void PRADA_DBN::checkStateSoundness(uint t, bool omit_derived) {
  uint i, val;
  double sum;
  FOR1D(vars_state__prim, i) {
    sum = 0.;
    FOR1D(vars_state__prim(i)->range, val) {
      if (vars_state__prim(i)->P(t,val) < -0.01) {
        cerr<<"invalid distribution for rv with id="<<vars_state__prim(i)->id<<"  for value="<<val<<": "<< vars_state__prim(i)->P(t,val) <<endl;
        HALT("");
      }
      sum += vars_state__prim(i)->P(t,val);
    }
    if (!TL::areEqual(1.0, sum))
      HALT("invalid distribution for rv with id="<<vars_state__prim(i)->id << " and literal " << *vars_state__prim(i)->lit);
  }
  if (!omit_derived) {
    FOR1D(vars_state__derived, i) {
      if (vars_state__derived(i)->type == LiteralRV::expectation)
        continue;
      sum = 0.;
      FOR1D(vars_state__derived(i)->range, val) {
        if (vars_state__derived(i)->P(t,val) < -0.01) {
          cerr<<"invalid distribution for rv with id="<<vars_state__derived(i)->id<<"  for value="<<val<<": "<<vars_state__derived(i)->P(t,val) <<endl;
          HALT("");
        }
        sum += vars_state__derived(i)->P(t,val);
      }
      if (!TL::areEqual(1.0, sum))
        HALT("invalid distribution for rv with id="<<vars_state__derived(i)->id
            <<" for literal "<<*vars_state__derived(i)->lit<<"   sum="<<sum);
    }
  }
}


void PRADA_DBN::getBelief(uint t, arr& beliefs_p_prim, arr& beliefs_f_prim) const {
  uint i, val;
  beliefs_p_prim.resize(vars_state__prim.N, 2);
  beliefs_p_prim.setZero();
  FOR1D(vars_state__prim, i) {
    FOR1D(vars_state__prim(i)->range, val) {
      beliefs_p_prim(i, val) += vars_state__prim(i)->P(t,val);
    }
  }
  NIY;
//   beliefs_f_prim.resize(vars_state__prim.N, MAX_FUNCTION_VALUE);
//   beliefs_f_prim.setZero();
//   FOR1D(vars_state__prim, i) {
//     FOR1D(vars_state__prim(i)->range, val) {
//       beliefs_f_prim(i, val) += vars_state__prim(i)->P(t,val);
//     }
//   }
}




/************************************************
 * 
 *     PRADA_DBN  -  write
 * 
 ************************************************/

void PRADA_DBN::writeAllStates(bool prim_only, double threshold, ostream& out) const {
  uint t;
  for(t=0; t<=horizon; t++) {
    out<<"+++++  t=" << t << "  +++++"<<endl;
    writeState(t, prim_only, threshold, out);
  }
}


void PRADA_DBN::writeState(uint t, bool prim_only, double threshold, ostream& out) const {
  uint v, val;
  FOR1D(vars_state__prim, v) {
    if (vars_state__prim(v)->dim == 2) {
      if (vars_state__prim(v)->P(t,1) >= threshold) {
        out<<*vars_state__prim(v)->lit<<" "<<vars_state__prim(v)->P(t,1)<<"  ";
      }
    }
    else {
      out<<*vars_state__prim(v)->lit<<"  ";
      FOR1D(vars_state__prim(v)->range, val) {
        out<<vars_state__prim(v)->range(val)<<":"<<vars_state__prim(v)->P(t,val)<<"  ";
      }
    }
  }
  if (!prim_only) {
    FOR1D(vars_state__derived, v) {
      if (vars_state__derived(v)->type == LiteralRV::expectation) {
        out<<*vars_state__derived(v)->lit<<" "<<vars_state__derived(v)->P(t,0)<<" [E]   ";
      }
      else {
        if (vars_state__derived(v)->dim == 2) {
          if (vars_state__derived(v)->P(t,1) >= threshold) {
            out<<*vars_state__derived(v)->lit<<" "<<vars_state__derived(v)->P(t,1)<<"  ";
          }
        }
        else {
          out<<*vars_state__derived(v)->lit<<"  ";
          FOR1D(vars_state__derived(v)->range, val) {
            out<<vars_state__derived(v)->range(val)<<":"<<vars_state__derived(v)->P(t,val)<<"  ";
          }
        }
      }
    }
  }
}


void PRADA_DBN::writeStateSparse(uint t, bool prim_only, ostream& out) const {
  uint v;
  FOR1D(vars_state__prim, v) {
    if (TL::areEqual(vars_state__prim(v)->P(t,1), 1.0)) {
      vars_state__prim(v)->lit->write(out); out<<"  ";
    }
  }
  if (!prim_only) {
    FOR1D(vars_state__derived, v) {
      if (TL::areEqual(vars_state__derived(v)->P(t,1), 1.0)) {
        vars_state__derived(v)->lit->write(out); out<<"  ";
      }
    }
  }
  NIY;
//   uint val;
//   FOR1D(vars_state__prim, v) {
//     if (vars_state__prim(v)->type == RV_TYPE__FUNC_EXPECT) {
//       vars_state__prim(v)->fi->write(out);
//       out<<"="<<vars_state__prim(v)->P(t,0)<<"  ";
//     }
//     else {
//       FOR1D(vars_state__prim(v)->range, val) {
//         if (TL::areEqual(vars_state__prim(v)->P(t,val), 1.0)) {
//           vars_state__prim(v)->fi->write(out);
//           out<<"="<<vars_state__prim(v)->range(val)<<"  ";
//         }
//       }
//     }
//   }
//   if (!prim_only) {
//     FOR1D(rvs_state__f_derived, v) {
//       if (rvs_state__f_derived(v)->type == RV_TYPE__FUNC_EXPECT) {
//         rvs_state__f_derived(v)->fi->write(out);
//         out<<"="<<rvs_state__f_derived(v)->P(t,0)<<"  ";
//       }
//       else {
//         FOR1D(rvs_state__f_derived(v)->range, val) {
//           if (TL::areEqual(rvs_state__f_derived(v)->P(t,val), 1.0)) {
//             rvs_state__f_derived(v)->fi->write(out);
//             out<<"="<<rvs_state__f_derived(v)->range(val)<<"  ";
//           }
//         }
//       }
//     }
//   }
  out<<endl;
}


void PRADA_DBN::getActions(LitL& actions, uint horizon) const {
  actions.clear();
  uint a, t;
  CHECK(horizon <= vars_action(0)->P.d0, "");
  for (t=0; t<horizon; t++) {
    FOR1D(vars_action, a) {
      if (TL::isZero(vars_action(a)->P(t,1) - 1.0)) {
        actions.append(vars_action(a)->lit);
        break;
      }
      if (vars_action.N == a) {
        HALT("action has not been found");
      }
    }
  }
}






/************************************************
 * 
 *     PRADA_DBN  -  create DBN
 * 
 ************************************************/


void PRADA_DBN::create_dbn_structure(const SymL& symbols_state, const SymL& actions) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"PRADA_DBN::create_dbn_structure [START]"<<endl;}
  if (DEBUG>0) {PRINT(symbols_state);  PRINT(actions);}
  CHECK(!net_constants.containsDoubles(), "constants ain't distinct, brother! look here: "<<net_constants);
  
  uint i;
  for (i=0; symbols_state.N>0 && i<symbols_state.N-1; i++) {
    if (symbols_state(i+1)->symbol_type == Symbol::primitive  &&  symbols_state(i)->symbol_type != Symbol::primitive)
      HALT("Symbols need to be ordered:  "<<symbols_state);
  }
  
  vars_state__prim.clear();
  vars_state__derived.clear();
  vars_action.clear();
  // Init efficiency structure
  SymL symbols_all;  symbols_all.append(symbols_state);  symbols_all.append(actions);
  CHECK(RVefficiency__symbols.N == 0, "DBN has already been built before!");
  RVefficiency__init(symbols_all);
  
  // Build mapping: action --> rules
  // Calc max #rules per grounded action
  uint max_rules_per_action=0, rules_per_action = 0;
  Literal* last_action = NULL;
  FOR1D_(ground_rules, i) {
    if (ground_rules.elem(i)->action != last_action) {
      rules_per_action = 1;
      last_action = ground_rules.elem(i)->action;
    }
    else
      rules_per_action++;
//     ground_rules.elem(i)->action->write(cout); cout<<"   rules_per_action="<<rules_per_action<<endl;
    if (rules_per_action > max_rules_per_action)
      max_rules_per_action = rules_per_action;
  }
  // max action arity
  uint max_arity_actions = 0;
  FOR1D(actions,i) {
    if (actions(i)->arity > max_arity_actions)
      max_arity_actions = actions(i)->arity;
  }
  if (DEBUG>0) {PRINT(max_rules_per_action);  PRINT(max_arity_actions); PRINT(actions.N);  PRINT(net_constants.N);}
  
  action2rules.resize(actions.N, pow(net_constants.N, max_arity_actions), max_rules_per_action);
  action2rules.setUni(9999); // dummy value
  action2rules_num.resize(actions.N, pow(net_constants.N, max_arity_actions));
  action2rules_num.setZero();
  FOR1D_(ground_rules, i) {
    uint idx_symbol = actions.findValue(ground_rules.elem(i)->action->s);
    uint idx_args = getIndex__PRADA_DBN(ground_rules.elem(i)->action->args);
    action2rules_num(idx_symbol, idx_args)++;
//     ground_rules.elem(i)->action->write(cout);   cout<<"  idx_symbol="<<idx_symbol<<"  idx_args="<<idx_args<<"  ";   cout<<"   action2rules_num(idx_symbol, idx_args)="<<action2rules_num(idx_symbol, idx_args)<<endl;
    action2rules(idx_symbol, idx_args, action2rules_num(idx_symbol, idx_args)-1) = i;
  }
  
  // Determine changeable symbols
  SymL changing_symbols;
  ground_rules.changingSymbols(changing_symbols);
  if (changing_symbols.N == 0)
    HALT("no changing concepts! would be stupid to plan");
  
  // Determine literals
  LitL lits_state;
  FOR1D(symbols_state, i) {
    LitL lits;
    Literal::getLiterals(lits, symbols_state(i), net_constants, 1.0);
    lits_state.setAppend(lits);
    if (DEBUG>3) {cout<<"Lits for "<<symbols_state(i)->name << ":  "<<lits<<endl;}
  }
  // Filter PIs with ground-rules
  //   Guideline: If no noise outcome, then only PIs appearing in ground rules are interesting and we can filter.
  bool do_filter = true; // filter if no noise outcomes
  FOR1D_(ground_rules, i) {
    if (!TL::isZero(ground_rules.elem(i)->probs.last())) {
      do_filter = false;
      break;
    }
  }
  
//   cout<<endl<<endl<<endl<<endl;
//   cout<<"GROUND RULES"<<endl;
//   ground_rules.write(cout);
//   cout<<endl<<endl<<endl<<endl;

  if (do_filter) {
    uint k;
    // (1) Use all derived predicates and their precessors.
    //     Only a hack: to account for rewards (which are often specified by means of derived predicates)
    SymL necessary_symbols;
    FOR1D(symbols_state, i) {
      if (symbols_state(i)->symbol_type != Symbol::primitive) {
        SymL defining_symbols;
        symbols_state(i)->getDefiningSymbols(defining_symbols);
        FOR1D(defining_symbols, k) {
          if (changing_symbols.findValue(defining_symbols(k)) < 0) // only take non-changeable into account for sure; OB DAS SO PASST??
            necessary_symbols.setAppend(defining_symbols(k));
        }
      }
    }
    FOR1D(symbols_state, i) {
      if (symbols_state(i)->symbol_type != Symbol::primitive)
        necessary_symbols.setAppend(symbols_state(i));
    }
    
    // (2) Use all zero-ary preds
    //     Another HACK: since these might be used for the rewards
    FOR1D(symbols_state, i) {
      if (symbols_state(i)->arity == 0)
        necessary_symbols.append(symbols_state(i));
    }
    
    // (3) Filter
    LitL lits_filtered;
    bool lit_used;
    FOR1D(lits_state, i) {
      // (i)  keep all literals from "necessary_symbols"
      if (necessary_symbols.findValue(lits_state(i)->s) >= 0) { // only primitives taken into account
        lit_used = true;
      }
      // (ii) keep all literals appearing in ground rules
      else {
        lit_used = false;
        FOR1D_(ground_rules, k) {
          if (ground_rules.elem(k)->usesAtom(lits_state(i))) {
            lit_used = true;
            break;
          }
        }
      }
      if (lit_used)
        lits_filtered.append(lits_state(i));
    }
    lits_state = lits_filtered;
  }
  
  LitL lits_action;
  FOR1D_(ground_rules, i) {
    lits_action.setAppend(ground_rules.elem(i)->action);
    CHECK(lits_action.last()->s->arity < 1  
          ||  reason::getConstants().findValue(lits_action.last()->args(0)) >= 0,
          "reason::getConstants() "<<reason::getConstants()<<" don't contain action arg for action "<<*lits_action.last());
  }
  
  
  if (DEBUG>0) {
    cout<<"+++++   FINAL LITERALS (after filtering)   +++++"<<endl;
    cout<<"SymbolicState lits["<<lits_state.N<<"]:  "<<lits_state<<endl;
    cout<<"Action lits["<<lits_action.N<<"]:  "<<lits_action<<endl;
  }
  uint a;
  uint rv_ids = 1;
  
  
  // Build random variables for state literals
  start_id_rv__symbol_state = rv_ids;
  FOR1D(lits_state, i) { // assumed to be in order p_prim to p_derived
    LiteralRV* var = NULL;
    if (lits_state(i)->s->range_type == Symbol::binary  ||  lits_state(i)->s->range_type == Symbol::integer_set)
      var = new LiteralRV;
    else
      var = new ExpectationRV;
    var->lit = lits_state(i);
    if (lits_state(i)->s->range_type == Symbol::binary) {
      var->dim = 2;
      var->range.resize(2); var->range(0)=0; var->range(1)=1;
    }
    else if (lits_state(i)->s->range_type == Symbol::integer_set) {
      if (lits_state(i)->s->range.N == 0) {
        for (a=0; a<MAX_FUNCTION_VALUE; a++) {
          lits_state(i)->s->range.append(a);
        }
      }
      var->dim = lits_state(i)->s->range.N;
      var->range = lits_state(i)->s->range;
    }
    var->P.resize(horizon+1, var->dim);
    var->P.setZero();
    var->id = rv_ids++;
    if (changing_symbols.findValue(lits_state(i)->s) >= 0)
      var->changeable = true;
    else
      var->changeable = false;
    if (lits_state(i)->s->symbol_type == Symbol::primitive)
      vars_state__prim.append(var);
    else
      vars_state__derived.append(var);
    RVefficiency__setAtom(lits_state(i), var);
  }
  // Build random variables for action literals
  FOR1D(lits_action, a) {
    if (lits_action(a) == Literal::getLiteral_default_action()) continue; // Noisy default rule
    LiteralRV* var = new LiteralRV;
    var->lit = Literal::get(lits_action(a)->s, lits_action(a)->args, 1.);
    var->dim = 2; // binary variables
    var->range.resize(2); var->range(0)=0; var->range(1)=1;
    var->P.resize(horizon+1, var->dim);
    var->P.setZero();
    var->id = rv_ids++;
    vars_action.append(var);
    RVefficiency__setAtom(Literal::get(lits_action(a)->s, lits_action(a)->args, 1.), var);
  }
  if (DEBUG>1) {
    cout<<"SymbolicState rvs prim:  "<<endl;write(vars_state__prim);cout<<endl;
    cout<<"SymbolicState rvs p_derived:  "<<endl;write(vars_state__derived);cout<<endl;
    cout<<"Action rvs:  "<<endl;write(vars_action);cout<<endl;
  }
  
  // (4) Pseudo-RVs for rules
  vars_rules_simple.resize(horizon, ground_rules.num());
  vars_rules.resize(horizon, ground_rules.num());
  
  if (DEBUG) {cout<<"PRADA_DBN::create_dbn_structure [END]"<<endl;}
}


void PRADA_DBN::create_dbn_params() {
  uint DEBUG = 0;
  if (DEBUG) {cout<<"create_dbn_params [START]"<<endl;}
  uint r, v, o, k;
  int val;
  // impacts_val
  impacts_val.clear();
  uint no_changeable_rvs = 0;
  FOR1D(vars_state__prim, v) {
    if (vars_state__prim(v)->changeable) no_changeable_rvs++;
  }
  FOR1D(vars_state__prim, v) {
    arr local_impacts_val(ground_rules.num(), vars_state__prim(v)->dim);
    local_impacts_val.setZero();
    LiteralRV* var = vars_state__prim(v);
    FOR1D_(ground_rules, r) {
      Rule* r_grounded = ground_rules.elem(r);
      double prob_noise_change = r_grounded->noise_changes / no_changeable_rvs * noise_softener;
      FOR1D(r_grounded->outcomes, o) {
        // noise outcome (= last outcome)
        if (o==r_grounded->outcomes.N-1) {
          // touches every changeable random variable in the state!
          if (var->changeable) {
            FOR1D(var->range, val) {
              local_impacts_val(r,val) += prob_noise_change * r_grounded->probs(o);
            }
          }
        }
        // non-noise
        else {
          Literal* lit_predicted = NULL;
          // check that literal appears somehow
          FOR1D(r_grounded->outcomes(o), k) {
            if (r_grounded->outcomes(o)(k)->s == var->lit->s  &&  r_grounded->outcomes(o)(k)->args == var->lit->args) {
              lit_predicted = r_grounded->outcomes(o)(k);
              break;
            }
          }
          if (lit_predicted == NULL)
            continue;
          int val = var->range.findValue(lit_predicted->value);
          CHECK(val>= 0, "value has not been found");
          local_impacts_val(r, val) += r_grounded->probs(o);
        }
      }
    }
    impacts_val.append(local_impacts_val);
  }
  // impacts_V
  impacts_V.clear();
  FOR1D(vars_state__prim, v) {
    arr local_impacts_V(ground_rules.num());
    local_impacts_V.setZero();
    FOR1D_(ground_rules, r) {
      FOR1D(vars_state__prim(v)->range, val) {
        local_impacts_V(r) += impacts_val(v)(r, val);
      }
    }
    impacts_V.append(local_impacts_V);
  }
  
  if (DEBUG>1) {
    cout<<"+++++   IMPACTS   +++++"<<endl;
    for (r=0; r<ground_rules.num(); r++) {
      cout<<"["<<r<<"]"<<endl;
      ground_rules.elem(r)->write(cout);
      FOR1D(vars_state__prim, v) {
        vars_state__prim(v)->lit->write(cout);cout<<":  ";
        for (val=0; val<vars_state__prim(v)->dim; val++) {
          cout<<vars_state__prim(v)->range(val)<<":"<<impacts_val(v)(r,val)<<"  "<<std::flush;
        }
        cout<<" -->  "<< impacts_V(v)(r) << endl;
      }
    }
  }
  
  // build helper array vars 
  
  vars_context.resize(ground_rules.num(), 20);
  FOR1D_(ground_rules, r) {
    uint v_current = 0;
    FOR1D(ground_rules.elem(r)->context,v) {
      if (v_current == vars_context.d1) {
        HALT("Too large context for ground rule r="<<r<<" (PRADA only copes with maximum 20 literals in rule contexts)");
      }
      vars_context(r, v_current++) = RVefficiency__atom2var(ground_rules.elem(r)->context(v));
    }
  } 
  
  if (DEBUG) {cout<<"create_dbn_params [END]"<<endl;}
}






/************************************************
 * 
 *     PRADA_DBN  -  RVefficiency
 * 
 ************************************************/

void PRADA_DBN::RVefficiency__init(const SymL& symbols) {
  RVefficiency__symbols = symbols;
  uint max_arity = 0;
  uint i;
  FOR1D(symbols, i) {
    if (symbols(i)->arity > max_arity)
      max_arity = symbols(i)->arity;
  }
  RVefficiency__LRV_A__structured.resize(symbols.N, (uint) pow(net_constants.N, max_arity));
  RVefficiency__LRV_A__structured.setUni(NULL);
  RVefficiency__LRV_A__flat.resize(RVefficiency__LRV_A__structured.N+10);
  RVefficiency__LRV_A__flat.setUni(NULL);
}


LiteralRV* PRADA_DBN::RVefficiency__atom2var(Literal* lit) const {
  int s_idx = RVefficiency__symbols.findValue(lit->s);
  if (s_idx<0) {
    MT_MSG("Symbol '"<<lit->s->name<<"' for lit="<<*lit<<" not in RVefficiency.");
    return NULL;
  }
  LiteralRV* rv = RVefficiency__LRV_A__structured(s_idx, getIndex__PRADA_DBN(lit->args));
  if (rv == NULL) {
//     MT_MSG("RV_Manager::l2v -- No LiteralRV for "<<*lit<<" with calculated s_idx="<<s_idx);
    // Some methods expect that NULL is returned in that case.
  }
  return rv;
}


LiteralRV* PRADA_DBN::RVefficiency__id2var(uint id_var) const {
  return RVefficiency__LRV_A__flat(id_var);
}


void PRADA_DBN::RVefficiency__setAtom(Literal* lit, LiteralRV* var) {
  int idx_symbol = RVefficiency__symbols.findValue(lit->s);
  uint idx_args = getIndex__PRADA_DBN(lit->args);
  RVefficiency__LRV_A__structured(idx_symbol, idx_args) = var;
  RVefficiency__LRV_A__flat(var->id) = var;
}













/************************************************
 * 
 *     LiteralRV
 * 
 ************************************************/

void LiteralRV::write(ostream& os) {
  os<<"Var over "<<*lit<<endl;
  os<<"id="<<id<<endl;
  os<<"dim="<<dim<<endl;
  os<<"range="<<range<<endl;
  os<<"changeable="<<changeable<<endl;
  uint d,t;
  for(d=0; d<dim; d++) {
    FOR1D(P, t) {
      os<<"  "<<P(t, d);
    }
    os<<endl;
  }
}


void ExpectationRV::write(ostream& os) {
  os<<"Var over "<<*lit<<endl;
  os<<"id="<<id<<endl;
  os<<"changeable="<<changeable<<endl;
  os<<"Expectation"<<endl;
  uint t;
  FOR1D(P, t) {
    os<<"  "<<P(t,0);
  }
  os<<endl;
}


void write(const RVL& vars) {
  uint i;
  FOR1D(vars, i) {
    vars(i)->write(cout);
  }
}


}
