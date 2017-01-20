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

#include "reason.h"
#include "reason.h"
#include <stdlib.h>

namespace relational {

/****************************************
    Grounding
  ***************************************/


bool reason::isGround(const Rule& r) {
  CHECK(!Rule::isDefaultRule(&r), "don't use this method for default rule");
  CHECK(r.action != Literal::getLiteral_doNothing(), "don't use this method for doNothing rule");
  uint i, j;
  FOR1D(r.context, i) {
    if (!reason::isGround(r.context(i)))
      return false;
  }
  if (!reason::isGround(r.action))
    return false;
  FOR1D(r.outcomes, i) {
    FOR1D(r.outcomes(i), j) {
      if (!reason::isGround(r.outcomes(i)(j)))
        return false;
    }
  }
  return true;
}


bool reason::isPurelyAbstract(const Rule& r) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"reason::isPurelyAbstract [START]"<<endl;}
  if (DEBUG>1) {cout<<r<<endl;  PRINT(getConstants());}
  if (Rule::isDefaultRule(&r)) return true;
  uint i, j;
  FOR1D(r.context, i) {
    if (!reason::isPurelyAbstract(r.context(i))) {
      if (DEBUG>0) {cout<<"No abstract context #i="<<i<<": "<<*r.context(i)<<endl;   cout<<"reason::isPurelyAbstract [END]"<<endl;}
      return false;
    }
  }
  if (!reason::isPurelyAbstract(r.action)) {
    if (DEBUG>0) {cout<<"No abstract action: "<<*r.action<<endl;   cout<<"reason::isPurelyAbstract [END]"<<endl;}
    return false;
  }
  FOR1D(r.outcomes, i) {
    FOR1D(r.outcomes(i), j) {
      if (!reason::isPurelyAbstract(r.outcomes(i)(j))) {
        if (DEBUG>0) {cout<<"No abstract outcome #i="<<i<<", j="<<j<<": "<<*r.outcomes(i)(j)<<endl;   cout<<"reason::isPurelyAbstract [END]"<<endl;}
        return false;
      }
    }
  }
  if (DEBUG>0) {cout<<"Rule is abstract."<<endl;   cout<<"reason::isPurelyAbstract [END]"<<endl;}
  return true;
}


bool reason::isGround_positives(Rule* r) {
  uint i, j;
  FOR1D(r->context, i) {
    if (r->context(i)->value>0  &&  !reason::isGround(r->context(i)))
      return false;
  }
  if (!reason::isGround(r->action))
    return false;
  FOR1D(r->outcomes, i) {
    FOR1D(r->outcomes(i), j) {
      if (!reason::isGround(r->outcomes(i)(j)))
        return false;
    }
  }
  return true;
}


void reason::calcGroundDeicticReferences(uintA& ground_drefs, const SymbolicState& state, const Literal* groundAction, const Rule* rule) {
  ground_drefs.clear();
  uintA drefs;
  rule->getDeicticRefs(drefs);  
  SubstitutionSet subs;
  bool covers = reason::calcSubstitutions_rule_groundAction(subs, state, groundAction, rule);
  if ((!covers) || ((rule->action->args.N > 0 || drefs.N > 0)  &&  subs.num() != 1)) {
    cerr<<endl<<endl<<endl;
    cerr<<"Bad rule:"<<endl;  rule->write(cerr);
    cerr<<"groundAction="<<*groundAction<<endl;
    if (!covers) HALT("Rule does not cover!");
    HALT("Bad ground action: subs.num()="<<subs.num());
  }

  uint i;
  FOR1D(drefs, i) {
    ground_drefs.append(subs.elem(0)->getSubs(drefs(i)));
  }
}




/****************************************
  TRANSITION KNOWLEDGE
***************************************/
  
void reason::calcSuccessorState(SymbolicState& suc, const SymbolicState& pre, const LitL& outcome) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"calcSuccessorState [START]"<<endl;
  if (DEBUG>0) {cout<<"pre (N="<<pre.lits.N<<"): "<<pre<<endl<<"Outcome: "<<outcome<<endl;}
  uint i, k;
	// keep literals from predecessor that are not specified in outcome
  FOR1D(pre.lits, i) {
    if (pre.lits(i)->s->symbol_type != Symbol::primitive) continue;
    FOR1D(outcome, k) {
      if (outcome(k)->s == pre.lits(i)->s
      && outcome(k)->args == pre.lits(i)->args) {
        if (outcome(k)->s->range_type != Symbol::binary && outcome(k)->comparison_type == Literal::comparison_offset) {    //apply offset
          double sumValue = pre.lits(i)->value + outcome(k)->value;
          if (outcome(k)->s->range_type == Symbol::integer_set) {   //clamp   
            uint rangeMin = outcome(k)->s->range.min();
            uint rangeMax = outcome(k)->s->range.max();  
            if (sumValue < rangeMin) sumValue = rangeMin;
            if (sumValue > rangeMax) sumValue = rangeMax;
          }
          suc.lits.append(Literal::get(pre.lits(i)->s, pre.lits(i)->args, sumValue));
        }
        break;
      }
    }
    if (k==outcome.N) { // add if not found
      suc.lits.append(pre.lits(i));
      if (DEBUG>2) {cout<<*pre.lits(i)<<" --> keep"<<endl;}
    }
    else {
      if (DEBUG>1) {cout<<*pre.lits(i)<<" --> DEL"<<endl;}
    }
  }

	// add outcome literals
  FOR1D(outcome, i) {
    if (outcome(i)->s->range_type != Symbol::binary && outcome(i)->comparison_type == Literal::comparison_offset) continue;
    if (outcome(i)->s->range_type == Symbol::binary && outcome(i)->value == 0) continue;  //do not add negated binary symbols, they are implicitly negative
    suc.lits.setAppend(outcome(i));
  }

  Literal::sort(suc.lits);
  reason::derive(&suc);
  if (DEBUG>0) {cout<<"suc (N="<<suc.lits.N<<"):  "<<suc<<endl;}
  if (DEBUG>0) cout<<"calcSuccessorState [END]"<<endl;
}


double reason::sampleSuccessorState_groundRule(SymbolicState& successor, const SymbolicState& predecessor, Rule* ground_rule, uint& flag) {
  uint outcome_id = TL::basic_sample(ground_rule->probs);
  if (outcome_id == ground_rule->outcomes.N - 1) { // noise outcome
    flag = STATE_TRANSITION__NOISE_OUTCOME;
    successor = predecessor;
  }
  else {
    flag = 0;
    reason::calcSuccessorState(successor, predecessor, ground_rule->outcomes(outcome_id));
  }
  if (ground_rule->outcome_rewards.N > 0)
    return ground_rule->outcome_rewards(outcome_id);
  else
    return 0.0;
}


double reason::sampleSuccessorState_groundRules(SymbolicState& successor, const SymbolicState& predecessor, const RuleSet& ground_rules, Literal* action, uint& flag) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"sampleSuccessorState_groundRules [START]"<<endl;}
  if (DEBUG>0) {cout<<"Predecessor state:" <<predecessor<<endl<<"Action:  "<<*action<<endl;}
  RuleSet ground_rules_cov;
  calc_coveringGroundRules_groundAction(ground_rules_cov, ground_rules, predecessor, action);
  if (DEBUG>0) {cout<<"Covering rules (#="<<ground_rules_cov.num()<<"):"<<endl;  ground_rules_cov.write();}
  CHECK(ground_rules_cov.num()>0, "No covering rules (at least default should cover)!");
  double value = 0.;
  if (ground_rules_cov.num() > 2  ||  ground_rules_cov.num() == 0) {
    if (DEBUG>0) {cout<<"At least two non-default, or zero, covering rules."<<endl;}
    value = TL::TL_DOUBLE_NIL;
  }
  else if (ground_rules_cov.num() == 2
    &&   !(Rule::isDefaultRule(ground_rules_cov.elem(0))  ||  Rule::isDefaultRule(ground_rules_cov.elem(1)))) { 
    if (DEBUG>0) {cout<<"Two non-default covering rules."<<endl;}
    value = TL::TL_DOUBLE_NIL;
  }
  // -> unique non-default covering rule
  else if (ground_rules_cov.num() == 2) {
    CHECK(!Rule::isDefaultRule(ground_rules_cov.elem(1)), "default rule should be first");
    // ground_rules_cov.elem(0)  =  noisy default rule
    // ground_rules_cov.elem(1)  =  unique non-default covering rule
    if (DEBUG>0) {cout<<"Unique non-default covering rule."<<endl;}
    value = sampleSuccessorState_groundRule(successor, predecessor, ground_rules_cov.elem(1), flag);
  }
  else {
    if (DEBUG>0) {cout<<"Unique covering rule."<<endl;}
    value = sampleSuccessorState_groundRule(successor, predecessor, ground_rules_cov.elem(0), flag);
  }
  if (DEBUG>0) {PRINT(value);  PRINT(successor);}
  if (DEBUG>0) {cout<<"sampleSuccessorState_groundRules [END]"<<endl;}
  return value;
}


double reason::sampleSuccessorState_groundRules(SymbolicState& successor, const SymbolicState& predecessor, const RuleSet& ground_rules, Literal* action) {
  uint flag;
  return sampleSuccessorState_groundRules(successor, predecessor, ground_rules, action, flag);
}


double reason::sampleSuccessorState_abstractRules(SymbolicState& successor, const SymbolicState& predecessor, const RuleSet& abstract_rules, Literal* action) {
  RuleSet coveringGroundRules;
  calc_coveringRules_groundAction(coveringGroundRules, abstract_rules, predecessor, action);
  return reason::sampleSuccessorState_groundRules(successor, predecessor, coveringGroundRules, action);
}


double reason::probability_groundRule(Rule* r_ground, const SymbolicState& pre, const SymbolicState& post, double noiseStateProbability) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"probability [START]"<<endl;
  if (DEBUG>1) {
    r_ground->write(cout);
    cout<<"PRE: ";pre.write(cout);cout<<endl;
    cout<<"POST: ";post.write(cout);cout<<endl;
  }
  CHECK(isGround_positives(r_ground), "Rule is not ground.");
	
  uintA covering_outcomes;
  calc_coveringOutcomes(covering_outcomes, r_ground, pre, post);
  uint o;
  double succ_prob;
  double totalProb = 0.0;
  FOR1D(covering_outcomes, o) {
		// noisy outcome
    if (covering_outcomes(o)==r_ground->outcomes.N-1)
      succ_prob = noiseStateProbability;
    else
      succ_prob = 1.0;
    totalProb += succ_prob * r_ground->probs(covering_outcomes(o));
    if (DEBUG>1) {
      cout<<"Outcome: ";write(r_ground->outcomes(covering_outcomes(o)));
      cout<<" --> P(s'|o,s,a,r)="<<succ_prob<<" * P(o|s,a,r)="<<r_ground->probs(covering_outcomes(o))<<endl;
    }
  }
	
  if (DEBUG>0) PRINT(totalProb);
  if (DEBUG>0) cout<<"probability [END]"<<endl;
	
  return totalProb;
}


double reason::probability_abstractRule(Rule* abstractRule, const SymbolicState& pre, Literal* groundAction, const SymbolicState& post, double noiseStateProbability, Substitution* sub) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"probability [START]"<<endl;
  CHECK(reason::isGround(groundAction), "Action is not ground");
  CHECK(reason::isPurelyAbstract(abstractRule->action), "Rule action has to be abstract!");
	
  SubstitutionSet subs;
  if (sub==NULL) {
    if (calcSubstitutions_rule_groundAction(subs, pre, groundAction, abstractRule))
      sub = subs.elem(0);
    else
      return 0.0;
  }
	
  if (DEBUG>1) {cout<<"Sub: ";sub->write(cout);cout<<endl;}
	
  Rule* r_ground = sub->apply(*abstractRule);
  double prob = probability_groundRule(r_ground, pre, post, noiseStateProbability);
  delete r_ground;
		
  if (DEBUG>0) PRINT(prob);
  if (DEBUG>0) cout<<"probability [END]"<<endl;
  return prob;
}


double reason::probability_defaultRule(Rule* defaultRule, const SymbolicState& pre, const SymbolicState& post, double noiseStateProbability) {
  if (pre == post)
    return defaultRule->probs(0) * (100.0 * noiseStateProbability);
  else
    return defaultRule->probs(1) * noiseStateProbability;
}





  
    
/****************************************
  COVERAGE
 ***************************************/

bool reason::calcSubstitutions_context(SubstitutionSet& subs, const SymbolicState& state, const Rule* rule, Substitution* actionSub, bool consider_non_unique_covering) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"calcSubstitutions_context [START]"<<endl;
  CHECK(actionSub!=NULL, "Action substitution must be provided!");
  CHECK(actionSub->num() <= rule->action->args.N, "");
  if (DEBUG>0) {
    cout<<*rule<<state<<endl;
    cout<<"Action sub: ";actionSub->write(cout);cout<<endl;
  }
  CHECK(subs.num()==0, "Already subs given!");
  CHECK(rule->action->args.N==actionSub->num(), "Incomplete actionSub.");
  if (!state.derived_lits_are_calculated) {MT_MSG("Derived symbols have not been derived in state!");}
  uintA actionSub_outs;
  actionSub->getOuts(actionSub_outs);
  //  The deep question of QUANTIFICATION with NEGATED FREE DEICTIC VARS
  //  Example rule:
  //    ACTION:  act(X)   CONTEXT:  p2(X), -p1(Y)
  //  Two possibilities:
  //   (i)  Representing all quantification:  \forall Y: -p1(Y)
  //        --> BAD BAAAAAAAAAAD  Don't do that!
  //            This is like introducing a formula "for all constants, the respective predicate does not hold".
  //   (ii) Representing exists quantification:  \exists Y: -p1(Y)
  //        --> does not make intuitive sense
  //            Example:   ACTION:  puton(X)    CONTEXT:  -on(Y,X)
  //                       Rule applicable if exactly one arbitrary object not on X
  //   Baseline: Should be forbidden in NID rules, in general.
  //             However, if there comes such a case nonetheless, we use exists quantification.
  bool mightbecovering = reason::calcSubstitutions(subs, state.lits, rule->context, false, actionSub);
  uint i, j;
  if (DEBUG>0) {
    PRINT(mightbecovering);
    cout<<"Substitutions (#="<<subs.num()<<"):"<<endl;
    FOR1D_(subs, i) {
      subs.elem(i)->write(cout);cout<<endl;
    }
    PRINT(actionSub_outs);
  }
  bool covering;
// 	// NON-DEICTIC [start]
// 	covering = mightbecovering;
// 	// NON-DEICTIC [end]
	// DEICTIC VARIANT [start]
  if (mightbecovering) {
    uintA drefs;
    rule->getDeicticRefs(drefs);
    if (DEBUG>0) {PRINT(drefs);}
    uint dref__subs_id;
    covering = true;
    // Filter substitutions
    FOR1D_DOWN_(subs, j) {
      FOR1D(drefs, i) {
        dref__subs_id = subs.elem(j)->getSubs(drefs(i));
        CHECK_(reason::isConstant(dref__subs_id),
               {rule->write(cout);  state.write(cout); cout<<endl;  cout<<"Action sub: "; actionSub->write(cout);cout<<endl;
                PRINT(dref__subs_id);  PRINT(drefs);
                cout<<"subs(j="<<j<<"):  ";  subs.elem(j)->write();  cout<<endl;},
                "No substitution for deictic reference.");
        // Deictic references must be different from action arguments
        if (actionSub_outs.findValue(dref__subs_id) >= 0) {
          if (DEBUG>0) {cout<<"Removing subs(j="<<j<<")  as  DR=" << dref__subs_id << " is action arg"<<endl;}
          subs.remove(j);
          break;
        }
      }
    }
    if (DEBUG>0) {PRINT(subs.num());}
    if (consider_non_unique_covering) { covering = subs.num() >= 1; }
    else { covering = subs.num() == 1; }
  }
  else
    covering = false;
	// DEICTIC VARIANT [end]
  if (DEBUG>0) PRINT(covering);
  if (DEBUG>0) cout<<"calcSubstitutions_context [END]"<<endl;
  return covering;
}


bool reason::calcSubstitutions_rule_groundAction(SubstitutionSet& subs, const SymbolicState& s, const Literal* groundAction, const Rule* rule, bool consider_non_unique_covering) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"calcSubstitutions_rule_groundAction [START]"<<endl;
  if (DEBUG>0) {
    cout<<"STATE:  "<<s<<endl;
    cout<<"GROUNDED ACTION:  "<<*groundAction<<endl;
    cout<<"RULE:"<<endl<<*rule;
    cout<<"SUBS (N="<<subs.num()<<"):"<<endl;  uint i; FOR1D_(subs, i) {cout<<"("<<i<<")  ";  subs.elem(i)->write();}
  }
  CHECK(subs.num()==0, "There are already substitutions!");
  CHECK(groundAction->s->symbol_type == Symbol::action, "Not an action symbol in groundAction-literal!");
	// applied world knowledge: default rule always covers
  if (Rule::isDefaultRule(rule)) {
    if (subs.num() == 0)
      subs.append(new Substitution);
    if (DEBUG>0) cout<<"isDefaultRule --> return true"<<endl<<"calcSubstitutions_rule_groundAction [END]"<<endl;
    return true;
  }
  // doNothing rule
  if (Literal::getLiteral_doNothing() == rule->action) {
    if (rule->action == groundAction)
      return true;
    else
      return false;
  }
  CHECK(isPurelyAbstract(*rule), "Rule has to be abstract:  "<<endl<<*rule);
  Substitution* actionSub = new Substitution;
  bool covers;
  if (reason::calcSubstitution(*actionSub, groundAction, rule->action))
    // ###### ESSENTIAL CALL
    covers = calcSubstitutions_context(subs, s, rule, actionSub, consider_non_unique_covering);
  else
    covers = false;
  delete actionSub;
	
  if (DEBUG>0) {
    cout<<"covers? "<<covers;
    if(covers){cout<<"  Sub: ";subs.elem(0)->write(cout);}
    cout<<endl;
  }
  if (DEBUG>0) cout<<"calcSubstitutions_rule_groundAction [END]"<<endl;
  return covers;
}


bool reason::calcSubstitutions_rule(SubstitutionSet& subs, const SymbolicState& state, const Rule* rule) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"calcSubstitutions_rule [START]"<<endl;
  if (DEBUG>0) {
    /*cout<<"Rule:"<<endl;*/rule->write(cout);
    cout<<"SymbolicState: "<<state<<endl;
  }
  MT::Array< uintA > list_actionArguments;
  TL::allPermutations(list_actionArguments, reason::getConstants(), rule->action->s->arity, true, true);
  uint i;
  FOR1D(list_actionArguments, i) {
    Literal* groundAction = Literal::get(rule->action->s, list_actionArguments(i), 1.);
    SubstitutionSet subs_inner;
    if (calcSubstitutions_rule_groundAction(subs_inner, state, groundAction, rule)) {
      subs.append(subs_inner);
    }
  }
  if (DEBUG>0) {cout<<"Covers: "<<(subs.num()>0)<<endl;}
  if (DEBUG>1) {
    FOR1D_(subs,i) {
      subs.elem(i)->write(cout);cout<<endl;
    }
  }
  if (DEBUG>0) cout<<"calcSubstitutions_rule [END]"<<endl;
  return subs.num()>0;
}


void reason::calc_coveringRules(RuleSet& r_grounds, const RuleSet& all_abstract_rules, const SymbolicState& s) {
  r_grounds.clear();
  uint i, k;
  for (i=0; i<all_abstract_rules.num(); i++) {
    SubstitutionSet subs;
    if (calcSubstitutions_rule(subs, s, all_abstract_rules.elem(i))) {
      FOR1D_(subs, k) {
        r_grounds.append(subs.elem(k)->apply(*all_abstract_rules.elem(i)));
      }
    }
  }
}


void reason::calc_coveringRules_groundAction(RuleSet& r_grounds, const RuleSet& all_abstract_rules, const SymbolicState& s, Literal* groundAction) {
  r_grounds.clear();
  uint i, k;
  for (i=0; i<all_abstract_rules.num(); i++) {
    SubstitutionSet subs;
    if (calcSubstitutions_rule_groundAction(subs, s, groundAction, all_abstract_rules.elem(i))) {
      FOR1D_(subs, k) {
        r_grounds.append(subs.elem(k)->apply(*all_abstract_rules.elem(i)));
      }
    }
  }
}


void reason::calc_coveringRules_groundAction(uintA& coveringRuleIDs, const RuleSet& all_abstract_rules, const SymbolicState& s, Literal* groundAction) {
  coveringRuleIDs.clear();
  uint i;	
  for (i=0; i<all_abstract_rules.num(); i++) {
    SubstitutionSet subs;
    if (calcSubstitutions_rule_groundAction(subs, s, groundAction, all_abstract_rules.elem(i)))
      coveringRuleIDs.append(i);
  }
}


void reason::calc_coveringGroundRules_groundAction(uintA& coveringRuleIDs, const RuleSet& all_ground_rules, const SymbolicState& s, Literal* groundAction) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"calc_coveringGroundRules_groundAction [START]"<<endl;}
  coveringRuleIDs.clear();
  uint i; 
  for (i=0; i<all_ground_rules.num(); i++) {
    if (DEBUG>0) {cout<<"Candidate rule:"<<endl;  all_ground_rules.elem(i)->write();}
    if (Rule::isDefaultRule(all_ground_rules.elem(i))) {
      coveringRuleIDs.append(i);
      if (DEBUG>0) {cout<<"--> Covers"<<endl;}
    }
    else if (all_ground_rules.elem(i)->action == groundAction  &&  reason::holds(s.lits, all_ground_rules.elem(i)->context)) {
      coveringRuleIDs.append(i);
      if (DEBUG>0) {cout<<"--> Covers"<<endl;}
    }
  }
  if (DEBUG>0) {cout<<"calc_coveringGroundRules_groundAction [END]"<<endl;}
}


void reason::calc_coveringGroundRules_groundAction(RuleSet& coveringGroundRules, const RuleSet& all_ground_rules, const SymbolicState& s, Literal* groundAction) {
  uintA coveringRuleIDs;
  reason::calc_coveringGroundRules_groundAction(coveringRuleIDs, all_ground_rules, s, groundAction);
  uint i;
  coveringGroundRules.clear();
  FOR1D(coveringRuleIDs, i) {
    coveringGroundRules.append(all_ground_rules.elem(coveringRuleIDs(i)));
  }
}


void reason::calc_coveringRules(uintA& coveringRulesIDs, const RuleSet& abstract_rules, const LitL& ground_actions, const SymbolicState& s) {
  coveringRulesIDs.resize(ground_actions.N);
  uint i;
  FOR1D(ground_actions, i) {
    uintA ids_covering_rules__action;
    reason::calc_coveringRules_groundAction(ids_covering_rules__action, abstract_rules, s, ground_actions(i));
    if (ids_covering_rules__action.N == 2)
      coveringRulesIDs(i) = ids_covering_rules__action(1);
    else
      coveringRulesIDs(i) = 0;   // default rule
  }
}


void reason::calc_coveringGroundRules(uintA& coveringRulesIDs, const RuleSet& ground_rules, const LitL& ground_actions, const SymbolicState& s) {
  coveringRulesIDs.resize(ground_actions.N);
  uint i;
  FOR1D(ground_actions, i) {
    uintA ids_covering_rules__action;
    reason::calc_coveringGroundRules_groundAction(ids_covering_rules__action, ground_rules, s, ground_actions(i));
    if (ids_covering_rules__action.N == 2)
      coveringRulesIDs(i) = ids_covering_rules__action(1);
    else
      coveringRulesIDs(i) = 0;   // default rule
  }
}


uint reason::calc_uniqueCoveringRule_groundRules_groundAction(const RuleSet& all_ground_rules, const SymbolicState& s, Literal* groundAction) {
  uintA coveringGroundRules_ids;
  calc_coveringGroundRules_groundAction(coveringGroundRules_ids, all_ground_rules, s, groundAction);
  if (coveringGroundRules_ids.N == 2) {
    CHECK(coveringGroundRules_ids(0) == 0, "first rule should be noisy default rule");
    return coveringGroundRules_ids(1);
  }
  else if (coveringGroundRules_ids.N == 1  &&  !Rule::isDefaultRule(all_ground_rules.elem(coveringGroundRules_ids(0)))) {
    return coveringGroundRules_ids(0);
  }
  else
    return 0;
}


uint reason::calc_uniqueAbstractCoveringRule_groundAction(const RuleSet& all_abstract_rules, const SymbolicState& s, Literal* groundAction) {
  uintA coveringRuleIDs;
  calc_coveringRules_groundAction(coveringRuleIDs, all_abstract_rules, s, groundAction);
  if (coveringRuleIDs.N == 2) {
    CHECK(coveringRuleIDs(0) == 0, "first rule should be noisy default rule");
    return coveringRuleIDs(1);
  }
  else
    return 0;
}


void reason::calc_coveringOutcomes(uintA& covering_outcomes, Rule* r_ground, const SymbolicState& pre, const SymbolicState& post) {
  covering_outcomes.clear();
  uint o;
  FOR1D(r_ground->outcomes, o) {
		// noisy outcome
    if (o==r_ground->outcomes.N-1) {
      covering_outcomes.append(o);
    }	
    else {
      SymbolicState suc;
      calcSuccessorState(suc, pre, r_ground->outcomes(o));
      if (suc == post)
        covering_outcomes.append(o);
    }
  }
}


void reason::calc_coveringOutcomes(uintA& covering_outcomes, Rule* abstractRule, const SymbolicState& pre, Literal* groundAction, const SymbolicState& post) {
  CHECK(reason::isPurelyAbstract(abstractRule->action), "Rule action has to be abstract!");
  SubstitutionSet subs;
  calcSubstitutions_rule_groundAction(subs, pre, groundAction, abstractRule);
  CHECK(subs.num()==1, "rule coverage only in case of exactly one sub");
  Rule* r_ground = subs.elem(0)->apply(*abstractRule);
  calc_coveringOutcomes(covering_outcomes, r_ground, pre, post);
  delete r_ground;
}




/****************************************
   SCORING
 ***************************************/
  
double reason::log_likelihood(const RuleSet& rules, const StateTransitionL& experiences, double p_min) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"log_likelihood [START]"<<endl;}
  if (DEBUG>0)  {PRINT(p_min);}
  double loglik = 0.;
  uint i;
  uint noise_predictions = 0;
  FOR1D(experiences, i) {
    if (DEBUG>1) {cout<<"+++++ Ex "<<i<<" +++++"<<endl;}
    if (DEBUG>2) {experiences(i)->write(cout, 1);}
    if (DEBUG>1) {cout<<"Ex-Changed: "<<experiences(i)->changes<<endl;}
    RuleSet coveringGroundRules;
    calc_coveringRules_groundAction(coveringGroundRules, rules, experiences(i)->pre, experiences(i)->action);
    CHECK(coveringGroundRules.num() > 0, "at least default rule should cover");
    if (DEBUG>1) {PRINT(coveringGroundRules.num());}
    Rule* explaining_rule = NULL;
    double lik = 0.;
    CHECK(Rule::isDefaultRule(coveringGroundRules.elem(0)), "first covering rule needs to be default rule");
    // Explain as non-noise
    if (coveringGroundRules.num() == 2) {explaining_rule = coveringGroundRules.elem(1);}
    // Explain as noise
    else {explaining_rule = coveringGroundRules.elem(0);  noise_predictions++;}
    if (DEBUG>1) {cout<<"Explaining rule:"<<endl<<*explaining_rule;}
    lik = probability_groundRule(explaining_rule, experiences(i)->pre, experiences(i)->post, p_min);
    if (explaining_rule == coveringGroundRules.elem(0)  &&  experiences(i)->changes.N == 0) {
      lik *= (100.0 * p_min);  // arbitrary choice: punish modeling persistence by means of default rule
      if (DEBUG>1) {cout<<"Modeling via default rule always get punished,"
                    <<" even when the no-change outcome applies instead of the noise outcome "
                    <<" (likelihood multiplied with p_min="<<p_min<<")."<<endl;}
    }
    if (DEBUG>1) {PRINT(lik);  PRINT(log(lik));}
    loglik += log(lik);
  }
  if (DEBUG>0) {
    cout<<"Non-noise predictions:  "<<((experiences.N - noise_predictions) * 1.0)/experiences.N<<"  ("<<experiences.N - noise_predictions<<")"<<endl;
    cout<<"Noise predictions:      "<<(noise_predictions * 1.0)/experiences.N<<"  ("<<noise_predictions<<")"<<endl;
    PRINT(loglik);
  }
  if (DEBUG>0) {cout<<"log_likelihood [END]"<<endl;}
  return loglik;
}



double reason::score(const RuleSet& rules, const StateTransitionL& experiences, double p_min, double alpha_pen) {
  uint DEBUG = 3;
  if (DEBUG>0) {cout<<"score [START]"<<endl;}
  // Loglik
  double loglik = log_likelihood(rules, experiences, p_min);
  if (DEBUG>0) {PRINT(loglik);}
  // Penalty
  double penalty = 0.0;
  uint i;
  FOR1D_(rules, i) {
    penalty += rules.elem(i)->numberLiterals();
  }
  penalty *= alpha_pen * (-1.0);
  if (DEBUG>0) {PRINT(penalty);}
  // Score
  double score = loglik + penalty;
  if (DEBUG>0) {PRINT(score);}
  if (DEBUG>0) {cout<<"score [END]"<<endl;}
  return score;
}





} // namespace PRADA



