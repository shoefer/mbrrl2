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


#include <math.h>
#include "reason.h"
#include "learn.h"

#include <set>

namespace relational {

/************************************************
 * 
 *     SearchOperator
 * 
 ************************************************/

uint SearchOperator::experience_idx_min = 0;

SearchOperator::SearchOperator() {
  this->approximative = false;
}


void SearchOperator::integrateNewRules(const RuleSetContainer& rulesC_old, const RuleSetContainer& rulesC_2add, 
                                       const StateTransitionL& experiences, RuleSetContainer& rulesC_new) {
  uint DEBUG = 0;
  if (DEBUG > 0) {cout << "integrateNewRules [START]" << endl;}
  if (DEBUG > 1) {
    cout << "***** Really old rules *****" << endl;
    rulesC_old.write(cout, true);
    cout << "***** rulesC_2add *****" << endl;
    rulesC_2add.write(cout, true);
    cout << "------" << endl;
  }
  rulesC_new.clear();
  uint i, j;
  // Create a copy of the input rule-set
  rulesC_new = rulesC_old;
  
  // For each new rule r'
  FOR1D_(rulesC_2add.rules, i) {
    // Remove rules in R' that cover any experiences r' covers
    uintA& covered_new = rulesC_2add.experiences_per_rule(i);
    // Clean-up rules
    rulesC_2add.rules.elem(i)->cleanup();
    // SPECIAL: default rule (at pos 0) is never removed
    for (j=rulesC_new.rules.num()-1; j>0; j--) {
      uintA& covered_old = rulesC_new.experiences_per_rule(j);
      if (numberSharedElements(covered_new, covered_old) > 0) {
        rulesC_new.remove(j);
      }
    }
    // Add r' to R'
    rulesC_new.append(rulesC_2add.rules.elem(i), rulesC_2add.experiences_per_rule(i), rulesC_2add.experiences_per_ruleOutcome(i));
  }
  // Recompute the set of experiences that the default rule in R' covers and the parameters of this default ruleLearner
  rulesC_new.recomputeDefaultRule();
  
  if (DEBUG > 1) {
    if (DEBUG>2) {
      cout << "***** Old rules *****" << endl;
      rulesC_old.write(cout, true);
    }
    cout << "***** New rules *****" << endl;
    rulesC_new.write(cout, true);
  }
  
  if (DEBUG > 0) {cout << "integrateNewRules [END]" << endl;}
}


const char* SearchOperator::getName() {
  return name;
}


// Algorithm of Figure 4 in Pasula et al. (2007)
void SearchOperator::createRuleSets(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights,
        MT::Array< RuleSetContainer >& set_of_new_rulesC) {
  uint DEBUG = 0;
  set_of_new_rulesC.clear();
  reset();
  uint i;
  while (true) {
    RuleSetContainer rulesC_2add(&experiences);
    // this is where the local knowledge of the individual search operators comes in
    findRules(rulesC_old, experiences, experience_weights, rulesC_2add);
    if (DEBUG>1) {
      FOR1D_(rulesC_2add.rules, i) {
        rulesC_2add.rules.elem(i)->write(cout);
      }
    }
    if (rulesC_2add.rules.num() == 0)
      break;
    RuleSetContainer rulesC_new(&experiences);
    integrateNewRules(rulesC_old, rulesC_2add, experiences, rulesC_new);
    set_of_new_rulesC.append(rulesC_new);
    if (approximative) {
      if (set_of_new_rulesC.d0 >= APPROXIMATOR__RULES_PER_ROUND)
        break;
    }
  }
}















// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    E X P L A I N   E X A M P L E S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

uint ExplainExperiences::maxRelationalActionDistance = 0; // <sebastian>

ExplainExperiences::ExplainExperiences(bool slimContext, bool comparingValues) : SearchOperator() {
  if (!slimContext) {
    if (!comparingValues)
      name = "ExplainExperiences";
    else
      name = "ExplainExperiences_comparingValues";
  }
  else {
    if (!comparingValues)
      name = "ExplainExperiences_slim";
    else
      name = "ExplainExperiences_slim_comparingValues";
  }
  //nextPotentialExperience = 0;
  nextPotentialExperience = experience_idx_min;
  this->slimContext = slimContext;
  this->comparingValues = comparingValues;
}


// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void ExplainExperiences::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"ExplainExperiences::findRules [START]"<<endl;
  uint i;
//  PRINT(Rule::globalRuleCounter)
  for (i=nextPotentialExperience; i<experiences.N; i++) {
    uintA& covering_rules = rulesC_old.nonDefaultRules_per_experience(i);
    if (DEBUG>1) {
      write(rulesC_old.rules);
      cout<<"#Covering non-default rules for example (" << i << ") = " << covering_rules.N << endl;
    }
    if (covering_rules.N == 0) {
      // Create new rule by explaining current example (sets context and action)
      if (DEBUG>0) cout << "findRules: let's explain #" << i << endl;
      if (DEBUG>2) {experiences(i)->write(cout, 2);}
      Rule* newRule = explainExperience(experiences(i));
      // calc experience coverage
      StateTransitionL covered_experiences;
      uintA covered_experiences_ids;
      arr covered_experience_weights;
      learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
      if (DEBUG>2) {cout<<"Learned new rule without outcomes:"<<endl<<*newRule;}
      // Estimate new outcomes for r'
      CHECK(covered_experiences.N>0  &&  covered_experiences_ids.findValue(i)>=0, "At least the explained example should be covered."<<endl<<*newRule)
      MT::Array< uintA > experiences_per_outcome;
      learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
      if (DEBUG>0) {cout<<"New Rule:"<<endl; newRule->write(cout);}
      if (DEBUG>0) {cout<<"#Covering experiences: "<<covered_experiences.N<<"   "<<covered_experiences_ids<<endl;}
      rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
//       rulesC_2add.sanityCheck();
      nextPotentialExperience = i+1;
      break;
    }
  }
//  PRINT(rules2_add.num());
//  PRINT(Rule::globalRuleCounter)
	if (DEBUG>0) cout<<"ExplainExperiences::findRules [END]"<<endl;
}


Rule* ExplainExperiences::explainExperience(StateTransition* ex) {
//   return explainExperience_deictic(ex);
  return explainExperience_deictic_ALL_DRs(ex);
}


// TODO Should be replaced by an explicit vocabulary that shall be used for the rules...
void filter_language(LitL& lits) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout << "filter_language before: [N="<<lits.N<<"] "<<lits<<endl;}
  // determine predicates to be filtered
  SymL filtered_symbols;
  uint i;
  FOR1D(lits, i) {
    if (lits(i)->s->name == "above") {
      filtered_symbols.append(lits(i)->s);
    }
  }
  if (DEBUG>0) {PRINT(filtered_symbols);}
  
  lits.memMove = true;
  FOR1D_DOWN(lits, i) {
    if (filtered_symbols.findValue(lits(i)->s) >= 0) {
      if (DEBUG>0) {cout<<"Removing "<<*lits(i)<<endl;}
      lits.remove(i);
    }
  }
  
  if (DEBUG>0) {cout << "filter_language after: [N="<<lits.N<<"] "<<lits<<endl;}
}


// TODO YOU MAY WANT TO FILTER SOME LITERALS WHICH ARE NEVER TRUE.
//  Exploit prior domain knowledge to make learning faster.
//  Some day in the future, should be replaced by more principled language constraints...
void trim_hack(LitL& lits) {
  uint DEBUG = 0;
  uint k;
  FOR1D(lits, k) {
    if (lits(k)->s->name == "table")
      break;
  }
  if (lits.N == k) // --> only for desktop world
    return;
  if (DEBUG>0) {cout << "trim_hack before: [N="<<lits.N<<"] "; write(lits); cout<<endl;}
  
  uint __TABLE_OBJECT_ID = 60; // ATTENTION: this has to be ensured!
  uint i;
  lits.memMove = true;
  
  // no self-referencing binary predicates
  // on(X,X) etc.
  FOR1D_DOWN(lits, i) {
    if (lits(i)->args.N == 2) {
      if (lits(i)->args(0) == lits(i)->args(1)) {
        if (DEBUG>0) {cout<<"Removing ";lits(i)->write(cout);cout<<endl;}
        lits.remove(i);
      }
    }
  }
  if (DEBUG>1) {cout << "trim_hack now: "<<lits<<endl;}
  
  // for X>__TABLE_OBJECT_ID: not table(X) [--> cares only about constants]
  Symbol* s_TABLE = Symbol::get(MT::String("table"));
  if (s_TABLE != NULL) {
    FOR1D_DOWN(lits, i) {
      if (lits(i)->s == s_TABLE) { // table
        if (lits(i)->args(0) > __TABLE_OBJECT_ID) {
          if (DEBUG>0) {cout<<"Removing ";lits(i)->write(cout);cout<<endl;}
            lits.remove(i);
        }
      }
    }
    if (DEBUG>1) {cout << "trim_hack now: "<<lits<<endl;}
  }
  
  // if on(X,Y) then -on(Y,X) is redundant
  Symbol* s_ON = Symbol::get(MT::String("on"));
  if (s_ON != NULL) {
    FOR1D_DOWN(lits, i) {
      if (!lits(i)->value > 0.) {
        if (lits(i)->s == s_ON) {
          uint k;
          bool remove = false;
          FOR1D_DOWN(lits, k) {
            if (lits(k)->value > 0.) {
              if (lits(k)->s) { // primitive
                if (lits(i)->args(0) == lits(k)->args(1)
                    &&     lits(i)->args(1) == lits(k)->args(0)) {
                  remove = true;
                  break;
                }
              }
            }
          }
          if (remove) {
            if (DEBUG>0) {cout<<"Removing ";  lits(i)->write(cout);cout<<endl;}
            lits.remove(i);
          }
        }
      }
    }
    if (DEBUG>1) {cout << "trim_hack now: "<<lits<<endl;}
  }
  
  // sort out specific concepts
  filter_language(lits);
  if (DEBUG>1) {cout << "trim_hack now: "<<lits<<endl;}
  
  if (DEBUG>0) {cout << "trim_hack after: [N="<<lits.N<<"] "; write(lits); cout<<endl;}
}


// Algorithm Pasula et al. (2007) p. 330
#if 0
Rule* ExplainExperiences::explainExperience_deictic(StateTransition* ex) {
  NIY;
  // TODO Do we still use this one at all??
  
  uint DEBUG = 0;
  if (DEBUG>0) cout << "explainExperience_deictic [START]" << endl;
  if (DEBUG>1) ex->write(cout);
  uint i, k;
  Rule* newRule = new Rule;
  
  // ensure that all complex are derived
  reason::derive(&ex->pre);
  reason::derive(&ex->post);
  uintA arguments, arguments_mustBeContained;
  
  // Step 1.1: Create an action and context
  // create action
  Substitution invSub;
  FOR1D(ex->action->args, i) {
    if (reason::isConstant(ex->action->args(i)) && !invSub.hasSubs(ex->action->args(i))) {
      invSub.addSubs2Variable(ex->action->args(i));
    }
  }
  newRule->action = invSub.apply(ex->action);
  if (DEBUG>1) {cout<<"New action: ";newRule->action->write(cout);cout<<endl;}
  // create context
  // create normal literals; first only with action arguments
  invSub.getIns(arguments);
  // (also accounts for negations)
  LitL context_candidates, context_candidates_pos, context_candidates_neg;
  Literal::getLiterals_state(context_candidates_pos, arguments, 1.0);
  Literal::getLiterals_state(context_candidates_neg, arguments, 0.0);
  context_candidates.append(context_candidates_pos);
  context_candidates.append(context_candidates_neg);
  
  // hack -- don't use complex reward-concepts [START]
  context_candidates.memMove = true;
  FOR1D_DOWN(context_candidates, k) {
    if ( context_candidates(k)->s->symbol_type == Symbol::function_reward
      || context_candidates(k)->s->name == "homies"
      )
      context_candidates.remove(k);
  }
  // hack -- don't use complex reward-concepts [END]
  
  if (DEBUG>2) {cout<<"Context literal candidates (based on action arguments, w./o. comparisons): ";write(context_candidates);cout<<endl;}
  FOR1D(context_candidates, i) {
    if (reason::holds(ex->pre.lits, context_candidates(i))) {
      if (DEBUG>3) {cout<<"Accepting context_candidates(i="<<i<<")="<<*context_candidates(i)<<endl;}
      newRule->context.append(invSub.apply(context_candidates(i)));
    }
    else {
      if (DEBUG>3) {cout<<"Rejecting context_candidates(i="<<i<<")="<<*context_candidates(i)<<endl;}
    }
  }
  
  // create comparison literals
  LitL equalityLiterals;
  FOR1D(ex->pre.lits, i) {
    if (ex->pre.lits(i)->s->range_type != Symbol::binary)
      equalityLiterals.append(ex->pre.lits(i));
  }
  
  // hack -- don't use complex reward-concepts [START]
  FOR1D_DOWN(equalityLiterals, k) {
    if ( equalityLiterals(k)->s->symbol_type == Symbol::count
      || equalityLiterals(k)->s->symbol_type == Symbol::function_count
      || equalityLiterals(k)->s->symbol_type == Symbol::avg
      || equalityLiterals(k)->s->symbol_type == Symbol::max
      || equalityLiterals(k)->s->symbol_type == Symbol::function_change
      || equalityLiterals(k)->s->symbol_type == Symbol::sum
      || equalityLiterals(k)->s->symbol_type == Symbol::function_reward
    ) {
      equalityLiterals.remove(k);
    }
  }
  // hack -- don't use complex reward-concepts [END]
  
  FOR1D(equalityLiterals, i) {
      newRule->context.append(invSub.apply(equalityLiterals(i)));
  }

  // order by positives first
  Literal::sort(newRule->context);
  if (DEBUG>1) {cout<<"Context (preliminary): ";write(newRule->context);cout<<endl;}
  
  // Step 1.2: Create deictic references and their literals
  if (DEBUG > 2) {
    cout << "ex->del: "; write(ex->del); cout << endl;
    cout << "ex->add: "; write(ex->add); cout << endl;
    PRINT(ex->changedConstants)
  }
  LitL newContext;
  FOR1D(ex->changedConstants, i) {
    if (!invSub.hasSubs(ex->changedConstants(i))) {
      if (DEBUG>1) {cout<<"Deictic candidate: "<<ex->changedConstants(i)<<endl;}
      Substitution newInvSub = invSub;
      Substitution sub;
      invSub.getInverse(sub);
      newInvSub.addSubs2Variable(ex->changedConstants(i));
      newInvSub.getIns(arguments);
      arguments_mustBeContained.clear();
      arguments_mustBeContained.append(ex->changedConstants(i));
      // create normal predicates (also accounts for negations)
      LitL context_candidates__pos, context_candidates__neg;
      Literal::getLiterals_state(context_candidates__pos, arguments, arguments_mustBeContained, 1.0);
      Literal::getLiterals_state(context_candidates__neg, arguments, arguments_mustBeContained, 0.0);
      context_candidates.append(context_candidates__pos);  context_candidates.append(context_candidates__neg);
      
      // hack -- don't use complex reward-concepts [START]
      context_candidates.memMove = true;
      FOR1D_DOWN(context_candidates, k) {
        if ( context_candidates(k)->s->range_type != Symbol::binary
          || context_candidates(k)->s->symbol_type == Symbol::count
          || context_candidates(k)->s->symbol_type == Symbol::function_count
          || context_candidates(k)->s->symbol_type == Symbol::avg
          || context_candidates(k)->s->symbol_type == Symbol::max
          || context_candidates(k)->s->symbol_type == Symbol::function_change
          || context_candidates(k)->s->symbol_type == Symbol::sum
          || context_candidates(k)->s->symbol_type == Symbol::function_reward
        ) {
          context_candidates.remove(k);
        }
      }
      // hack -- don't use complex reward-concepts [END]
            
      // create constant-bound comparison literals
      LitL equalityLiterals;
      uintA changedConstantWrapper;
      changedConstantWrapper.append(ex->changedConstants(i));
      NIY;
//       logicObjectManager::getCompLiterals_constantBound(equalityLiterals, changedConstantWrapper, ex->pre, 0);
      
      // hack -- don't use complex reward-concepts [START]
      equalityLiterals.memMove = true;
      FOR1D_DOWN(equalityLiterals, k) {
        NIY;
//         ComparisonLiteral* clit = (ComparisonLiteral*) equalityLiterals(k);
//         if (((ComparisonLiteral*)clit->atom)->fa1->f->category >= category_derived)
//           equalityLiterals.remove(k);
      }
      // hack -- don't use complex reward-concepts [END]
      
      context_candidates.append(equalityLiterals);
      if (comparingValues) {
        // create dynamic-bound comparison literals
        uintA vars_sofar;
        newInvSub.getIns(vars_sofar);
        LitL dynamic_comparisons;
        NIY;
//         logicObjectManager::getCompLiterals_dynamicBound(dynamic_comparisons, vars_sofar, ex->pre, 0);
//         if (DEBUG>3) {cout<<"Dynamic-bound comparison literals:  ";  write(dynamic_comparisons); cout<<endl;}
//         FOR1D(dynamic_comparisons, k) {
//           // hack -- don't use complex reward-concepts [START]
//           if (((ComparisonLiteral*) dynamic_comparisons(k)->atom)->fa1->f->category >= category_derived)
//             continue;
//           if (dynamic_comparisons(k)->args.findValue(ex->changedConstants(i)) >= 0) {
//             context_candidates.append(dynamic_comparisons(k));
//           }
//         }
      }
      if (DEBUG>2) {cout<<"Context literal candidates (based on deictic candidate): ";write(context_candidates);cout<<endl;}
      // create possible newContext
      newContext = newRule->context;
      FOR1D(context_candidates, k) {
        if (DEBUG>3) {cout<<"context_candidates(k)="<<*context_candidates(k)<<" accepted? ";}
        if (reason::holds(ex->pre.lits, context_candidates(k))) {
          newContext.append(newInvSub.apply(context_candidates(k)));
          if (DEBUG>3) {cout<<" yes"<<endl;}
        }
        else {
          if (DEBUG>3) {cout<<" no"<<endl;}
        }
      }
      Literal::sort(newContext);
      if (DEBUG>1) {cout<<"Context with deic ref (potential): ";write(newContext);cout<<endl;}
      // trim literals
      trim_hack(newContext);
      if (slimContext) {NIY;}
      if (DEBUG>0) {
        cout << "newContext after trimming: ";
        write(newContext);
        cout << endl;
      }
      // check whether new variables refers uniquely to s
      SubstitutionSet subs;
      // check whether truly deictic ref (only one sub)
      bool covers = reason::cover(subs, ex->pre.lits, newContext, true, &sub);
      if (covers && subs.num()==1) {
        // check for neg free DRs
        Rule helper_rule;
        helper_rule.action = newRule->action;
        helper_rule.context = newContext;
        uintA drefs_pos, drefs_neg;
        helper_rule.getDeicticRefs(drefs_pos, drefs_neg);
        if (drefs_neg.N == 0) {
          invSub = newInvSub;
          newRule->context = newContext;
          if (DEBUG>1) {cout<<"Accepted ("<<ex->changedConstants(i)<<")"<<endl;}
        }
        else {
          if (DEBUG>1) {cout << "Not accepted ("<<ex->changedConstants(i)<<"):  drefs_neg=" << drefs_neg << endl;}
        }
      }
      else {
        if (DEBUG>1) {cout << "Not accepted ("<<ex->changedConstants(i)<<"):  covers="<<covers<<"  subs.num()=" << (subs.num()) << "=/=1"<< endl;}
      }
    }
  }

  if (DEBUG>0) {
    newRule->write(cout);
  }
  if (DEBUG>0) cout << "explainExperience_deictic [END]" << endl;
  
  return newRule;
}
#endif



Rule* ExplainExperiences::explainExperience_deictic_ALL_DRs(StateTransition* ex) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout << "explainExperience_deictic_ALL_DRs [START]" << endl;}
  if (DEBUG>1) {cout<<"Experience:"<<endl; ex->write(cout, 2);}
  uint i;
  Rule* newRule = new Rule;
  
  // ensure that all complex are derived
  reason::derive(&ex->pre);
  reason::derive(&ex->post);
  
  // Step 1.1: Create an action and context
  // create action
  Substitution invSub;
  FOR1D(ex->action->args, i) {
    if (reason::isConstant(ex->action->args(i)) && !invSub.hasSubs(ex->action->args(i))) {
      invSub.addSubs2Variable(ex->action->args(i));
    }
  }
  newRule->action = invSub.apply(ex->action);
  if (DEBUG>1) {cout<<"New action: ";newRule->action->write(cout);cout<<endl;}

  /////////////////////////////
  // <sebastian>
  // add and substitute constants from action up to a certain relational "distance"
  // in the pre state (TODO post state might also be interesting?)
  // TODO add relational distance > 1
  std::set<uint> addArgCandidates;

  if (maxRelationalActionDistance > 0) {
    uint l;
    FOR1D(ex->pre.lits, l) {
        FOR1D(ex->action->args, i) {
            Literal* lit = ex->pre.lits(l);
            // does literal contain this action arg? if not ignore
            if (!lit->args.contains(ex->action->args(i))) continue;

            // literal contains action args. add all of the literals args as candidates
            uint j;
            FOR1D(lit->args, j) {
                addArgCandidates.insert(lit->args(j));
            }
        }
    }

    if (DEBUG>0) {
        cout << "Distance 1";
        for(std::set<uint>::iterator it=addArgCandidates.begin(); it!=addArgCandidates.end(); ++it)
            cout << *it << " ";
        cout << endl;
    }

    if (maxRelationalActionDistance > 1) {
      // FIXME relational distance = 2
      FOR1D(ex->pre.lits, l) {
        for(std::set<uint>::iterator it=addArgCandidates.begin(); it!=addArgCandidates.end(); ++it) {
          uint arg = *it;
          Literal* lit = ex->pre.lits(l);
          // does literal contain this action arg? if not ignore
          if (!lit->args.contains(arg)) continue;

          // literal contains action args. add all of the literals args as candidates
          uint j;
          FOR1D(lit->args, j) {
              addArgCandidates.insert(lit->args(j));
          }
        }
      }
      if (DEBUG>0) {
        cout << "Distance 1";
        for(std::set<uint>::iterator it=addArgCandidates.begin(); it!=addArgCandidates.end(); ++it)
            cout << *it << " ";
        cout << endl;
      }
    }

    // now add all addArgCandidates which are not substituted yet to substitution
    for(std::set<uint>::iterator it=addArgCandidates.begin(); it!=addArgCandidates.end(); ++it) {
      uint arg = *it;
      if (reason::isConstant(arg) && !invSub.hasSubs(arg)) {
        invSub.addSubs2Variable(arg);
      }
    }
  }
  // </sebastian>
  /////////////////////////////

  // create context
 
  // Step 1.2: Create deictic references and their literals
  Substitution newInvSub = invSub;
  Substitution sub_action;
  invSub.getInverse(sub_action);
  if (DEBUG>0) {
    PRINT(invSub);
    PRINT(newInvSub);
    PRINT(sub_action);
  }

  // now add the changed constants
  if (ex->changedConstants.N < 6) {
    FOR1D(ex->changedConstants, i) {
      if (newInvSub.hasSubs(ex->changedConstants(i))) continue;
      newInvSub.addSubs2Variable(ex->changedConstants(i));
    }
  }


  LitL newContext;
//   LitL newContext_ground;

  FOR1D(ex->pre.lits, i) {
    Literal *substituted = newInvSub.apply(ex->pre.lits(i));
    if (reason::isPurelyAbstract(substituted)
      &&  ex->pre.lits(i)->s->symbol_type != Symbol::count
      &&  ex->pre.lits(i)->s->symbol_type != Symbol::avg
      &&  ex->pre.lits(i)->s->symbol_type != Symbol::max
      &&  ex->pre.lits(i)->s->symbol_type != Symbol::sum
      &&  ex->pre.lits(i)->s->symbol_type != Symbol::function_change
      &&  ex->pre.lits(i)->s->symbol_type != Symbol::function_reward
    ) {
      newContext.append(substituted);
//       newContext_ground.append(ex->pre.lits(i));
    }
  }
  
  //negated literals are not specified in the context, they must be added seperately
  uintA context_vars;  context_vars.setAppend(ex->action->args);
  if (ex->changedConstants.N < 6) context_vars.setAppend(ex->changedConstants);
  TL::sort_asc(context_vars);
  if (DEBUG>1) {PRINT(context_vars);}
  LitL ground_context_candidates_neg;
  Literal::getLiterals_state(ground_context_candidates_neg, context_vars, 0.0);

  FOR1D(ground_context_candidates_neg, i) {
    if (DEBUG>3) {cout<<"ground_context_candidates_neg(k)="<<*ground_context_candidates_neg(i)<<" accepted? ";}
    if (reason::holds(ex->pre.lits, ground_context_candidates_neg(i))) {
      newContext.append(newInvSub.apply(ground_context_candidates_neg(i)));
//       newContext_ground.append(ground_context_candidates_neg(i));
      if (DEBUG>3) {cout<<" yes"<<endl;}
    }
    else {
      if (DEBUG>3) {cout<<" no"<<endl;}
    }
  }
//   Literal::sort(newContext_ground);
  Literal::sort(newContext);
      
  if (DEBUG>0) {
//     PRINT(newContext_ground);
    PRINT(newContext);
  }

#ifdef NO_DEICTICREFS_BY_NONBINARY
  SubstitutionSet subsNB;
  //build rule that has the context of r with all non-binary symbols removed
  relational::LitL contextWithoutNonBinaries;
  bool containsNonBinaries = false;
  FOR1D(newContext, i) {
    if (newContext(i)->s->range_type == Symbol::binary)
      contextWithoutNonBinaries.append(newContext(i));
    else containsNonBinaries = true;
  }
  if (containsNonBinaries) {
    bool covers = reason::calcSubstitutions(subsNB, ex->pre.lits, contextWithoutNonBinaries, true, &sub_action);
    if (covers && subsNB.num() > 1) {
      if (DEBUG>1) { cout << "Rule not accepted because NO_DEICTICREFS_BY_NONBINARY is specified!" << endl; }
      return newRule;
    }
  }
#endif

  // check whether new variables refers uniquely to s
  SubstitutionSet subs;
    // check whether truly deictic ref (only one sub)
  bool covers = reason::calcSubstitutions(subs, ex->pre.lits, newContext, true, &sub_action);
  if (covers && subs.num()==1) {
    // check for neg free DRs
    Rule helper_rule;
    helper_rule.action = newRule->action;
    helper_rule.context = newContext;
    uintA drefs_neg, drefs_pos, drefs_nonBinary;
    helper_rule.getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
    if (drefs_neg.N == 0 && drefs_nonBinary.N == 0) {
      invSub = newInvSub;
      newRule->context = newContext;
      if (DEBUG>1) {cout<<"Rule accepted"<<endl;}
    }
    else {
      if (DEBUG>1) {cout << "Rule not accepted:  drefs_neg=" << drefs_neg << endl;}
    }
  }
  else {
    if (DEBUG>1) {
      cout << "Rule not accepted:  covers="<<covers<<"  subs.num()=" << (subs.num()) << "=/=1"<< endl;
      subs.write();
    }
  }
    
  if (DEBUG>0) {newRule->write(cout);}
  if (DEBUG>0) cout << "explainExperience_deictic_ALL_DRs [END]" << endl;
  return newRule;
}



void ExplainExperiences::reset() {
  nextPotentialExperience = experience_idx_min;
}







// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    D R O P   C O N D I T I O N S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

DropContextLiterals::DropContextLiterals() : SearchOperator() {
  name = "DropContextLiterals";
  nextRule = 0;
  nextContextLiteral = 0;
}


// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void DropContextLiterals::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"DropContextLiterals::findRules [START]"<<endl;
  uint r, p, i;
  Rule* newRule;
  bool stop = false;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    for (p=nextContextLiteral; p<rulesC_old.rules.elem(r)->context.N; p++) {
      newRule = new Rule;
      FOR1D(rulesC_old.rules.elem(r)->context, i) {
        if (i!=p)
          newRule->context.append(rulesC_old.rules.elem(r)->context(i));
      }
      newRule->action = rulesC_old.rules.elem(r)->action;
      // check for neg free DRs
      uintA drefs_neg, drefs_pos, drefs_nonBinary;
      newRule->getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
      if (drefs_neg.N > 0 || drefs_nonBinary.N > 0) {
        delete newRule;
        continue;
      }
      // Ensure that new rule covers experiences
      StateTransitionL covered_experiences;
      uintA covered_experiences_ids;
      arr covered_experience_weights;
      // FIXME could remove
      learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
      if (covered_experiences.N > 0) {
          if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
        if (DEBUG>3) {
          cout<<"Covered experiences:"<<endl;
          uint k;
          FOR1D(covered_experiences, k) {
            covered_experiences(k)->write(cout, 2);
          }
        }
        MT::Array< uintA > experiences_per_outcome;
        learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
        rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
        stop = true;
        nextContextLiteral += 1;
        break;
      }
      else {
        delete newRule;
      }
    }
    if (stop)
      break;
    
    nextRule=r+1;
    nextContextLiteral=0;
  }
  if (DEBUG>0) cout<<"DropContextLiterals::findRules [END]"<<endl;
}


void DropContextLiterals::reset() {
	nextRule = 0;
	nextContextLiteral = 0;
}






// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    D R O P   C O N D I T I O N S approximative Version
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

DropContextLiterals_approximativeVersion::DropContextLiterals_approximativeVersion() : SearchOperator() {
  name = "DropContextLiterals_approximativeVersion";
  approximative = true;
  prepareTotalNewSearch = false;
}


// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void DropContextLiterals_approximativeVersion::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"DropPre_approx::findRules [START]"<<endl;
  rulesC_2add.clear();
  uint i, k;
  if (prepareTotalNewSearch) {
    usableContextLiterals.clear();
    uint numContextLiterals=0;
    FOR1D_(rulesC_old.rules, i) {
      numContextLiterals += rulesC_old.rules.elem(i)->context.N;
    }
    usableContextLiterals.resize(numContextLiterals);
    usableContextLiterals.setUni(1);
    prepareTotalNewSearch=false;
    if (DEBUG>0) cout<<"Prepared new usableContextLiterals-Array with N="<<usableContextLiterals.N<<endl;
  }
  uint random_contextLiteral, id_contextLiteral, id_rule, starting_id_for_current_rule;
  uint num_open_context;
  Rule* newRule;
  uint resamples=0;
  while (rulesC_2add.rules.num()==0) {
    num_open_context = sum(usableContextLiterals);
    if (num_open_context==0)
        break;
    // randomly choose which context literal to delete
    random_contextLiteral = rnd.num(num_open_context);
    if (DEBUG>1) {cout<<"New rule-finding try:"<<endl; PRINT(num_open_context); PRINT(random_contextLiteral); PRINT(usableContextLiterals);}
    id_contextLiteral=0;
    id_rule=0;
    starting_id_for_current_rule=0;
    // find the correct rule and context for the random_contextLiteral
    FOR1D(usableContextLiterals, k) {
      // If target surely not in current rule "rulesC_old.rules.elem(id_rule)"...
      // Counter for id_rule
      if (k >= rulesC_old.rules.elem(id_rule)->context.N+starting_id_for_current_rule) {
        id_rule++;
        while (rulesC_old.rules.elem(id_rule)->context.N == 0)  // account for empty rules thereafter
          id_rule++;
        starting_id_for_current_rule=k;
      }
      // Counter for id_contextLiteral
      //    found
      if (usableContextLiterals(k) && id_contextLiteral==random_contextLiteral) {
        usableContextLiterals(k)=0;
        id_contextLiteral=k-starting_id_for_current_rule;
        break;
      }
      //   searching
      if (usableContextLiterals(k))
        id_contextLiteral++;
    }
    if (DEBUG>1) {
      cout<<"Thinking about deleting in rule "<<id_rule<<" context literal "<<id_contextLiteral<<endl;
      rulesC_old.rules.elem(id_rule)->write(cout);
      cout<<"Context literal of interest:  "; rulesC_old.rules.elem(id_rule)->context(id_contextLiteral)->write(cout); cout<<endl;
    }
    
    // bias for positive
    if (rulesC_old.rules.elem(id_rule)->context(id_contextLiteral)->value > 0. && resamples < DROP_NEGATIVE_BIAS) {
      if (DEBUG>1) {cout<<" (Don't wanna delete positive lit: "; rulesC_old.rules.elem(id_rule)->context(id_contextLiteral)->write(cout); cout<<")"<<endl;}
      usableContextLiterals(k)=1; // set back
      resamples++;
      continue;
    }
    else
      resamples=0;
    newRule = new Rule;
    FOR1D(rulesC_old.rules.elem(id_rule)->context, k) {
      if (k!=id_contextLiteral)
        newRule->context.append(rulesC_old.rules.elem(id_rule)->context(k));
    }
    if (DEBUG>1) {cout<<"Deletion of "<<*rulesC_old.rules.elem(id_rule)->context(id_contextLiteral)<<" executed."<<endl;}
    newRule->action = rulesC_old.rules.elem(id_rule)->action;
    // check for neg free DRs
    uintA drefs_neg, drefs_pos, drefs_nonBinary;
    newRule->getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
    if (drefs_neg.N > 0 || drefs_nonBinary.N > 0) {
      delete newRule;
      continue;
    }
    // calc experience coverage
    StateTransitionL covered_experiences;
    uintA covered_experiences_ids;
    arr covered_experience_weights;
//    learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);

    // <FIX> added by sebastian
    // we might actually have dropped a rule that potentially covers other experiences non-uniquely
    static const bool consider_non_unique_covering = true;
    learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights, consider_non_unique_covering);
    // no check that all experiences that were covered are uniquely covered
    FOR1D(covered_experiences, i) {
      SubstitutionSet subs;
      bool covered = reason::calcSubstitutions_rule_groundAction(subs, covered_experiences(i)->pre, covered_experiences(i)->action, newRule);
      if (!covered) {
        if (DEBUG>0) { cout << "WARN: experience " << covered_experiences_ids.elem(i) << " was non-uniquely covered. Discarding!" << endl; }
        covered_experiences.clear();
        break;
      }
    }
    // </FIX>

    if (DEBUG>1) {cout<<"Potential new rule:"<<endl;newRule->write(cout);}
    if (covered_experiences.N > 0) {
      if (DEBUG>1) cout<<"+++++ Covers "<<covered_experiences.N<<" experiences " << covered_experiences_ids << " and will be kept."<<endl;
      if (DEBUG>3) {
        cout<<"Covered experiences:"<<endl;
        uint k;
        FOR1D(covered_experiences, k) {
          covered_experiences(k)->write(cout);
        }
      }
      MT::Array< uintA > experiences_per_outcome;
      learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
      rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
    }
    else {
      delete newRule;
    }
  }
  
  if (DEBUG>1) cout<<"# rules found = "<<rulesC_2add.rules.num()<<endl;
  if (DEBUG>0) cout<<"DropPre_approx::findRules [END]"<<endl;
}


void DropContextLiterals_approximativeVersion::reset() {
}


void DropContextLiterals_approximativeVersion::reset_total_approximator() {
  prepareTotalNewSearch = true;
}






// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    D R O P   R E F E R E N C E S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------


DropReferences::DropReferences() : SearchOperator() {
  name = "DropReferences";
  nextReference = 0;
  nextRule=0;
}

// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void DropReferences::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"DropReferences::findRules [START]"<<endl;
  uint r, p, i;
  bool stop = false;
  Rule* newRule = NULL;
  uintA drefs;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    rulesC_old.rules.elem(r)->getDeicticRefs(drefs);
    for (p=nextReference; p<drefs.N; p++) {
      if (DEBUG>0) {
        cout<<"Removing deictic reference "<<drefs(nextReference)<<" in rule:"<<endl;
        rulesC_old.rules.elem(r)->write(cout);
      }
      newRule = new Rule;
      newRule->action = rulesC_old.rules.elem(r)->action;
      FOR1D(rulesC_old.rules.elem(r)->context, i) {
        Literal* lit = rulesC_old.rules.elem(r)->context(i);
        if (DEBUG>4) {PRINT(lit->args);}
        if (lit->s->range_type != Symbol::binary) {
          //(ANDREAS): Funktioniert hier nicht der selbe Code wie für Binaries?
          if (lit->args.findValue(drefs(nextReference))<0)
            newRule->context.append(rulesC_old.rules.elem(r)->context(i));
        }
        else if (lit->args.findValue(drefs(nextReference))<0)
          newRule->context.append(rulesC_old.rules.elem(r)->context(i));
      }
      if (DEBUG>0) {cout<<"Yielding the new context: "; write(newRule->context); cout<<endl;}
      // check for neg free DRs
      uintA drefs_neg, drefs_pos, drefs_nonBinary;
      newRule->getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
      if (drefs_neg.N > 0 || drefs_nonBinary.N > 0) {
        delete newRule;
        continue;
      }
      // calc experience coverage
      StateTransitionL covered_experiences;
      uintA covered_experiences_ids;
      arr covered_experience_weights;
//      learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);

      // <FIX> added by sebastian
      // we might actually have dropped a rule that potentially covers other experiences non-uniquely
      static const bool consider_non_unique_covering = true;
      learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights, consider_non_unique_covering);

      // no check that all experiences that were covered are uniquely covered
      FOR1D(covered_experiences, i) {
          SubstitutionSet subs;
          bool covered = reason::calcSubstitutions_rule_groundAction(subs, covered_experiences(i)->pre, covered_experiences(i)->action, newRule);
          if (!covered) {
              if (DEBUG>0) { cout << "WARN: experience " << covered_experiences_ids.elem(i) << " was non-uniquely covered. Discarding!" << endl; }
              covered_experiences.clear();
              break;
          }
      }
      // </FIX>

      if (covered_experiences.N > 0) {        
        MT::Array< uintA > experiences_per_outcome;
        learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);

        bool keep = true;
        // <FIX> not working
        // if we use experiences for crossvalidation also see if they are covered
        /*
        if (rulesC_old.p_experiences_crossvalidation) {
          StateTransitionL cv_covered_experiences;
          uintA cv_covered_experiences_ids;
          arr cv_covered_experience_weights;

          arr p_experiences_crossvalidation_weights(rulesC_old.p_experiences_crossvalidation->N);
          p_experiences_crossvalidation_weights.setUni(1.);
          learn::calcCoverage(cv_covered_experiences, cv_covered_experiences_ids, cv_covered_experience_weights,
                              newRule,
                              * (rulesC_old.p_experiences_crossvalidation), p_experiences_crossvalidation_weights, consider_non_unique_covering);

          // now check whether outcome also coverss

          // TODO if outcome does NOT cover (i.e. different
          if (cv_covered_experiences.N > 0) {
            //if (DEBUG>0)
            { cout << "WARN: rule covers cross-validation experiences. Discarding!" << endl; }
            covered_experiences.clear();
          }
        } */
        // </FIX>

        if (keep) {
          if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
          if (DEBUG>3) {
            cout<<"Covered experiences:"<<endl;
            uint k;
            FOR1D(covered_experiences, k) {
              covered_experiences(k)->write(cout);
            }
          }

          rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
          stop = true;
          nextReference += 1;
          break;
        }
      }
      else {
        delete newRule;
      }
    }
    if (stop)
      break;

    nextRule=r+1;
    nextReference=0;
  }
  if (DEBUG>0) {if (newRule!=NULL) {cout<<"Found rule: ";newRule->write(cout);cout<<endl;}}
  
  if (DEBUG>0) cout<<"DropReferences::findRules [END]"<<endl;
}


void DropReferences::reset() {
	nextRule=0;
	nextReference = 0;
}









// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    D R O P   R U L E S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------


DropRules::DropRules() : SearchOperator() {
  name = "DropRules";
}
    

void DropRules::createRuleSets(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights,
          MT::Array< RuleSetContainer >& one_rulesC_new) {
  uint DEBUG = 0;
  one_rulesC_new.clear();
  uint i, j;
  if (DEBUG>0) {cout<<"Old rule-set:"<<endl;  rulesC_old.write(cout, true); cout<<endl<<endl;}
  for (i=1; i<rulesC_old.rules.num(); i++) { // default rule must always be in
    if (DEBUG>0) {cout<<"Dropping rule #" << i << endl;}
    RuleSetContainer rulesC_new(&experiences);
//     rulesC_new = rulesC_old;
    FOR1D_(rulesC_old.rules, j) {
      if (i!=j) {
        if (DEBUG>0) {cout<<"j="<<j<<":  "; PRINT(rulesC_old.experiences_per_rule(j));}
        rulesC_new.append(rulesC_old.rules.elem(j), rulesC_old.experiences_per_rule(j), rulesC_old.experiences_per_ruleOutcome(j));
      }
    }
    rulesC_new.recomputeDefaultRule();
//     rulesC_new.sanityCheck();
    one_rulesC_new.append(rulesC_new);
  }
}

void DropRules::reset() {}
void DropRules::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {}











// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    S P L I T   O N   L I T E R A L S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

SplitOnLiterals::SplitOnLiterals() : SearchOperator() {
  name = "SplitOnLiterals";
  nextRule=1; // ignore default rule
  nextLiteral=0;
}


void SplitOnLiterals::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"SplitOnLiterals::findRules [START]"<<endl;
  rulesC_2add.clear();
  uint r, k;
  Rule* newRule_pos = NULL;
  Rule* newRule_neg = NULL;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    if (absentLiterals.N == 0) { // first round
      rulesC_old.rules.elem(r)->getAbsentLiterals(absentLiterals, true);
      // use refs [START]
      uintA rule_args;
      rulesC_old.rules.elem(r)->getArguments(rule_args);
      newVar = 0;
      FOR1D(rule_args, k) {
        if (rule_args(k) == newVar)
          newVar++;
      }
      uintA lit_args;
      lit_args.append(rule_args);
      lit_args.append(newVar);
      uintA wrapper;
      wrapper.append(newVar);
      LitL lits_dr;
      Literal::getLiterals_state(lits_dr, lit_args, wrapper, 1.0, true);
      // hack -- don't use complex reward-concepts [START]
      FOR1D_DOWN(lits_dr, k) {
        if ( lits_dr(k)->s->range_type != Symbol::binary
          || lits_dr(k)->s->symbol_type == Symbol::count
          || lits_dr(k)->s->symbol_type == Symbol::avg
          || lits_dr(k)->s->symbol_type == Symbol::max
          || lits_dr(k)->s->symbol_type == Symbol::function_change
          || lits_dr(k)->s->symbol_type == Symbol::sum
          || lits_dr(k)->s->symbol_type == Symbol::function_reward
        ) {
          lits_dr.remove(k);
        }
      }
      // hack -- don't use complex reward-concepts [END]
      if (DEBUG>2) {
        cout << "Calculated restriction literals for rule:"<<endl<<*rulesC_old.rules.elem(r);
        cout<<"Restriction literals: "<<lits_dr<<endl;
      }
      
      absentLiterals.append(lits_dr);
      // use refs [END]
      
      if (DEBUG>2) {
        cout << "Calculated absent literals for rule:"<<endl;
        rulesC_old.rules.elem(r)->write(cout);
        cout<<"Absent literals: ";write(absentLiterals);cout<<endl;
      }
    }
    while (nextLiteral<absentLiterals.N) {
      // positive version
      LitL wrapper;
      wrapper.append(absentLiterals(nextLiteral));
      if (DEBUG>1) {
        cout<<"Trying to insert ";absentLiterals(nextLiteral)->write(cout);cout<<"   into   ";
        write(rulesC_old.rules.elem(r)->context);cout<<endl;
      }
      if (Literal::nonContradicting(wrapper, rulesC_old.rules.elem(r)->context)
             && rulesC_old.rules.elem(r)->context.findValue(absentLiterals(nextLiteral)) < 0) {
        if (DEBUG>1) cout<<" --> Feasible and will be done."<<endl;
        newRule_pos = new Rule;
        newRule_pos->action = rulesC_old.rules.elem(r)->action;
        newRule_pos->context = rulesC_old.rules.elem(r)->context;
        newRule_pos->insertContext(absentLiterals(nextLiteral));
        StateTransitionL covered_experiences;
        uintA covered_experiences_ids;
        arr covered_experience_weights;
        learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule_pos, experiences, experience_weights);
        if (covered_experiences.N > 0) {
          if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
          MT::Array< uintA > experiences_per_outcome;
          learn::learn_outcomes(newRule_pos, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
          rulesC_2add.append(newRule_pos, covered_experiences_ids, experiences_per_outcome);
        }
        else {
          if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
          delete newRule_pos;
        }
      }
      else {
        if (DEBUG>1) cout<<" --> Impossible."<<endl;
      }
      // negative version
      Literal* nextLiteral_neg = absentLiterals(nextLiteral)->getNegated();
      wrapper.clear();
      wrapper.append(nextLiteral_neg);
      if (DEBUG>1) {
        cout<<"Trying to insert ";nextLiteral_neg->write(cout);cout<<"   into   ";
        write(rulesC_old.rules.elem(r)->context);cout<<endl;
      }
      // check that negated literal does not introduce a negative deictic reference
      uintA args;
      rulesC_old.rules.elem(r)->getArguments(args);
      bool no_negated_dref = true;
      FOR1D(nextLiteral_neg->args, k) {
        if (args.findValue(nextLiteral_neg->args(k)) < 0) {
          no_negated_dref = false;
          break;
        }
      }
      // insert
      if (Literal::nonContradicting(wrapper, rulesC_old.rules.elem(r)->context)
             && rulesC_old.rules.elem(r)->context.findValue(nextLiteral_neg) < 0
             && no_negated_dref
             && nextLiteral_neg->args.findValue(newVar) < 0) {    //don't add negated deictic references
        if (DEBUG>1) cout<<" --> Feasible and will be done."<<endl;
        newRule_neg = new Rule;
        newRule_neg->action = rulesC_old.rules.elem(r)->action;
        newRule_neg->context = rulesC_old.rules.elem(r)->context;
        newRule_neg->insertContext(nextLiteral_neg);
        StateTransitionL covered_experiences;
        uintA covered_experiences_ids;
        arr covered_experience_weights;
        learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule_neg, experiences, experience_weights);
        if (covered_experiences.N > 0) {
          if (DEBUG>1) {cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;}
          if (DEBUG>3) {
            cout<<"Covered experiences:"<<endl;
            uint k;
            FOR1D(covered_experiences, k) {
              covered_experiences(k)->write(cout);
            }
          }
          MT::Array< uintA > experiences_per_outcome;
          learn::learn_outcomes(newRule_neg, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
          rulesC_2add.append(newRule_neg, covered_experiences_ids, experiences_per_outcome);
        }
        else {
          if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
          delete newRule_neg;
        }
      }
      else {
        if (DEBUG>1) cout<<" --> Impossible."<<endl;
      }
      
      nextLiteral++;
      if (rulesC_2add.rules.num() > 0)
        break;
    }
    if (rulesC_2add.rules.num() > 0)
      break;
    else {
      nextRule = r+1;
      nextLiteral = 0;
      absentLiterals.clear();
    }
  }
  if (DEBUG>0) {if (rulesC_2add.rules.num() > 0) {write(rulesC_2add.rules);}}
  if (DEBUG>0) cout<<"SplitOnLiterals::findRules [END]"<<endl;
}

void SplitOnLiterals::reset() {
  nextRule = 1;
  nextLiteral = 0;
}


// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    S P L I T  O N  L I T E R A L  C O N J U N C T I O N S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

SplitOnLiteralConjunctions::SplitOnLiteralConjunctions() : SearchOperator() {
  name = "SplitOnLiteralConjunctions";
  nextRule=1; // ignore default rule
  nextLiteral1=nextLiteral2=0;
}


void SplitOnLiteralConjunctions::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"SplitOnLiteralConjunctions::findRules [START]"<<endl;
  rulesC_2add.clear();
  uint r, k;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    if (absentLiterals.N == 0) { // first round
      rulesC_old.rules.elem(r)->getAbsentLiterals(absentLiterals, true);
      // use refs [START]
      uintA rule_args;
      rulesC_old.rules.elem(r)->getArguments(rule_args);
      newVar = 0;
      FOR1D(rule_args, k) {
        if (rule_args(k) == newVar)
          newVar++;
      }
      uintA lit_args;
      lit_args.append(rule_args);
      lit_args.append(newVar);
      uintA wrapper;
      wrapper.append(newVar);
      LitL lits_dr;
      Literal::getLiterals_state(lits_dr, lit_args, wrapper, 1.0, true);
      // hack -- don't use complex reward-concepts [START]
      FOR1D_DOWN(lits_dr, k) {
        if ( lits_dr(k)->s->range_type != Symbol::binary
          || lits_dr(k)->s->symbol_type == Symbol::count
          || lits_dr(k)->s->symbol_type == Symbol::avg
          || lits_dr(k)->s->symbol_type == Symbol::max
          || lits_dr(k)->s->symbol_type == Symbol::function_change
          || lits_dr(k)->s->symbol_type == Symbol::sum
          || lits_dr(k)->s->symbol_type == Symbol::function_reward
        ) {
          lits_dr.remove(k);
        }
      }
      // hack -- don't use complex reward-concepts [END]

      if (DEBUG>2) {
        cout << "Calculated restriction literals for rule:"<<endl<<*rulesC_old.rules.elem(r);
        cout<<"Restriction literals: "<<lits_dr<<endl;
      }

      absentLiterals.append(lits_dr);
      // use refs [END]

      if (DEBUG>0) {
        cout << "Calculated absent literals for rule:"<<endl;
        rulesC_old.rules.elem(r)->write(cout);
        cout<<"Absent literals: ";write(absentLiterals);cout<<endl;
      }
    }

    /// FIXME
//    absentLiterals.clear();
//    absentLiterals.append( Literal::get("rightMost(Z)")  );
//    absentLiterals.append( Literal::get("isRev(Z)")  );

    uintA already_covered_experiences;
      if (DEBUG>0) cout << "===" << endl;
    for (; nextLiteral1<absentLiterals.N;nextLiteral1++) {
      if (absentLiterals(nextLiteral1)->s->arity > 1) continue;

      bool inner_loop_break=false;
      if (DEBUG>0) cout << "--" << endl;
      for (; nextLiteral2<nextLiteral1;nextLiteral2++) {
        if (absentLiterals(nextLiteral2)->s->arity > 1) continue;
        if (nextLiteral1 == nextLiteral2) continue;

        // forbid identical predicates
        if (absentLiterals(nextLiteral1)->s == absentLiterals(nextLiteral2)->s)
          continue;

        // forbid predicates with different args
        if (absentLiterals(nextLiteral1)->args.elem(0) != absentLiterals(nextLiteral2)->args.elem(0))
          continue;

        if (DEBUG>0) cout << ".." << endl;
        // 4 combinations
        for(uint l1=0; l1<=1; l1++) {

          for(uint l2=0; l2<=1; l2++) {
            if (l1 == 0 && l2 == 0) continue;

            bool one_negated = l1 == 0 || l2 == 0;

            LitL wrapper;
//            wrapper.clear();

            Literal *lit1, *lit2;
            lit1=Literal::get(absentLiterals(nextLiteral1)->s, absentLiterals(nextLiteral1)->args, absentLiterals(nextLiteral1)->value);
            lit2=Literal::get(absentLiterals(nextLiteral2)->s, absentLiterals(nextLiteral2)->args, absentLiterals(nextLiteral2)->value);
            if (l1 == 0) lit1=lit1->getNegated();
//            else lit1 = absentLiterals(nextLiteral1);

            if (l2 == 0) lit2=lit2->getNegated();
//            else lit2 = absentLiterals(nextLiteral2);

            wrapper.append(lit1);
            wrapper.append(lit2);

            if (DEBUG>0) {
                cout<<"Trying to insert ";lit1->write(cout);cout<<", ";lit2->write(cout);cout<<"   into   ";
                write(rulesC_old.rules.elem(r)->context);cout<<endl;
            }

            uintA args;
            rulesC_old.rules.elem(r)->getArguments(args);
            bool no_negated_dref = true;
            bool no_in_context = true;
            uint j;
            FOR1D(wrapper, j) {
              // avoid neg deictive references
              if (one_negated && wrapper.elem(j)->isNegated() ){
                // collect other wrapper args in positive literals
                uint jj,kk;
                uintA other_wrapper_args;
                FOR1D(wrapper, jj) {
                  if (jj == j) continue;
                  if (wrapper.elem(jj)->isNegated()) continue;
                  FOR1D(wrapper.elem(jj)->args, kk) {
                    other_wrapper_args.append(wrapper.elem(jj)->args.elem(kk));
                  }
                }
                // for each arg of current literal
                FOR1D(wrapper.elem(j)->args, k) {
                  // is in context or in another wrapper?
                  if (args.findValue(wrapper.elem(j)->args(k)) < 0
                      && other_wrapper_args.findValue(wrapper.elem(j)->args(k))) {
                    no_negated_dref = false;
                    break;
                  }
                }

                if (!no_negated_dref)
                  break;
              }
              // make sure the literals are not already in the context
              if (rulesC_old.rules.elem(r)->context.findValue(wrapper.elem(j)) >= 0) {
                no_in_context = false;
                break;
              }
            }
            if (DEBUG>1) cout<<" no_negated_dref --> " << no_negated_dref <<endl;
            if (DEBUG>1) cout<<" no_in_context --> " << no_in_context <<endl;

            if (Literal::nonContradicting(wrapper, rulesC_old.rules.elem(r)->context)
                  && no_in_context && no_negated_dref) {
              if (DEBUG>1) cout<<" --> Feasible and will be done."<<endl;
              Rule* newRule = new Rule;
              newRule->action = rulesC_old.rules.elem(r)->action;
              newRule->context = rulesC_old.rules.elem(r)->context;
              newRule->insertContext(lit1);
              newRule->insertContext(lit2);

              if (DEBUG>1) {
                cout<<" ** NEW RULE: "<<endl;
                newRule->write(cout);
              }

              StateTransitionL covered_experiences;
              uintA covered_experiences_ids;
              arr covered_experience_weights;
              learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
              if (covered_experiences.N > 0) {
                bool covers_already_covered_experiences=false;
                uint l;
                FOR1D(already_covered_experiences,l) {
                  if (covered_experiences_ids.contains(already_covered_experiences.elem(l))) {
                    covers_already_covered_experiences=true;
                    break;
                  }
                }
                if (covers_already_covered_experiences) {
                  if (DEBUG>0) cout<<"  Covers already covered experience "<< l <<" and will be dropped."<<endl;
                  delete newRule;
                } else {
                  if (DEBUG>0) cout<<"  Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
                  MT::Array< uintA > experiences_per_outcome;
                  learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
                  rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);

                  already_covered_experiences.append(covered_experiences_ids);
                }
              }
              else {
                  if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
                  delete newRule;
              }
            }
            else {
                if (DEBUG>1) cout<<" --> Impossible."<<endl;
            }
          }
        }
        // now check that the four combinations actually generated at least two rules which
        // cover different things
        if (rulesC_2add.rules.num() > 1) {
          if (DEBUG>0) cout<< rulesC_2add.rules.num() << " rules found. breaking! (inner loop)"<<endl;
          nextLiteral2++;
          inner_loop_break =true;
          break;
        } else {
          already_covered_experiences.clear();
          rulesC_2add.clear();
        }
      }
      if (inner_loop_break)
        break;
      else
        nextLiteral2 = 0;

      if (rulesC_2add.rules.num() > 1) {
        if (DEBUG>0) cout<< rulesC_2add.rules.num() << "rules found. breaking!"<<endl;
        nextLiteral1++;
        break;
      }
    }

    if (rulesC_2add.rules.num() > 0) {
      break;
    }
    else {
      nextRule = r+1;
      nextLiteral1 = nextLiteral2 = 0;
      absentLiterals.clear();
    }
  }

  if (DEBUG>0) {if (rulesC_2add.rules.num() > 0) {write(rulesC_2add.rules);}}
  if (DEBUG>0) cout<<"SplitOnLiteralConjunctions::findRules [END]"<<endl;
}

void SplitOnLiteralConjunctions::reset() {
  nextRule = 1;
  nextLiteral1 = nextLiteral2 = 0;
}






// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    A D D   L I T E R A L S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

AddLiterals::AddLiterals() : SearchOperator() {
  name = "AddLiterals";
  nextRule=1; // ignore default rule
  nextLiteral=0;
}


void AddLiterals::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"AddLits::findRules [START]"<<endl;
  uint r;
  Rule* newRule = NULL;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    if (absentLiterals.N == 0) { // first round
      rulesC_old.rules.elem(r)->getAbsentLiterals(absentLiterals);
      // hack -- don't use complex reward-concepts [START]
      absentLiterals.memMove = true;
      uint k;
      // hack -- don't use complex reward-concepts [START]
      FOR1D_DOWN(absentLiterals, k) {
        if ( absentLiterals(k)->s->range_type != Symbol::binary
          || absentLiterals(k)->s->symbol_type == Symbol::count
          || absentLiterals(k)->s->symbol_type == Symbol::avg
          || absentLiterals(k)->s->symbol_type == Symbol::max
          || absentLiterals(k)->s->symbol_type == Symbol::function_change
          || absentLiterals(k)->s->symbol_type == Symbol::sum
          || absentLiterals(k)->s->symbol_type == Symbol::function_reward
        ) {
          absentLiterals.remove(k);
        }
      }
      // hack -- don't use complex reward-concepts [END]
      if (DEBUG>2) {
        cout << "Calculated absent literals for rule:"<<endl;
        rulesC_old.rules.elem(r)->write(cout);
        cout<<"Absent literals: ";write(absentLiterals);cout<<endl;
      }
    }
    while (nextLiteral<absentLiterals.N) {
      LitL wrapper;
      wrapper.append(absentLiterals(nextLiteral));
      if (DEBUG>1) {
        cout<<"Trying to insert ";absentLiterals(nextLiteral)->write(cout);cout<<"   into   ";
        write(rulesC_old.rules.elem(r)->context);cout<<endl;
      }
      if (Literal::nonContradicting(wrapper, rulesC_old.rules.elem(r)->context)) {
        if (DEBUG>1) cout<<" --> Feasible and will be done."<<endl;
        newRule = new Rule;
        newRule->action = rulesC_old.rules.elem(r)->action;
        newRule->context = rulesC_old.rules.elem(r)->context;
        newRule->insertContext(absentLiterals(nextLiteral));
        StateTransitionL covered_experiences;
        uintA covered_experiences_ids;
        arr covered_experience_weights;
        learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
        if (covered_experiences.N > 0) {
          if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
          if (DEBUG>3) {
            cout<<"Covered experiences:"<<endl;
            uint k;
            FOR1D(covered_experiences, k) {
              covered_experiences(k)->write(cout);
            }
          }
          MT::Array< uintA > experiences_per_outcome;
          learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
          rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
        }
        else {
          if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
          delete newRule;
        }
      }
      else {
        if (DEBUG>1) cout<<" --> Impossible."<<endl;
      } 
      nextLiteral++;
      if (rulesC_2add.rules.num() > 0)
        break;  
    }
    if (rulesC_2add.rules.num() > 0)
      break;
    else {
      nextRule = r+1;
      nextLiteral = 0;
      absentLiterals.clear();
    }
  }
  if (DEBUG>0) {if (rulesC_2add.rules.num() > 0) {newRule->write(cout);}}
  if (DEBUG>0) cout<<"AddLits::findRules [END]"<<endl;
}


void AddLiterals::reset() {
	nextRule = 1;
	nextLiteral = 0;
}









// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    A D D   R E F E R E N C E S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

AddReferences::AddReferences() : SearchOperator() {
  name = "AddReferences";
  nextRule=1; // ignore default rule
  nextLiteral=0;
}


void AddReferences::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"AddReferences::findRules [START]"<<endl;
  uint r, i;
  Rule* newRule = NULL;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    if (restrictionLiterals.N == 0) { // first round
      uintA rule_args;
      rulesC_old.rules.elem(r)->getArguments(rule_args);
      uint newVar = 0;
      FOR1D(rule_args, i) {
        if (rule_args(i) == newVar)
          newVar++;
      }
      uintA arguments;
      arguments.append(rule_args);
      arguments.append(newVar);
      uintA wrapper;
      wrapper.append(newVar);
      Literal::getLiterals_state(restrictionLiterals, arguments, wrapper, 1.);
      // hack -- don't use complex reward-concepts [START]
      restrictionLiterals.memMove = true;
      FOR1D_DOWN(restrictionLiterals, i) {
        if ( restrictionLiterals(i)->s->range_type != Symbol::binary
          || restrictionLiterals(i)->s->symbol_type == Symbol::count
          || restrictionLiterals(i)->s->symbol_type == Symbol::avg
          || restrictionLiterals(i)->s->symbol_type == Symbol::max
          || restrictionLiterals(i)->s->symbol_type == Symbol::function_change
          || restrictionLiterals(i)->s->symbol_type == Symbol::sum
          || restrictionLiterals(i)->s->symbol_type == Symbol::function_reward
        ) {
          restrictionLiterals.remove(i);
        }
      }
      // hack -- don't use complex reward-concepts [END]
      if (DEBUG>2) {
        cout << "Calculated restriction literals for rule "<<rulesC_old.experiences_per_rule(r)<<":"<<endl;
        rulesC_old.rules.elem(r)->write(cout);
        cout<<"Restriction literals: ";write(restrictionLiterals);cout<<endl;
      }
    }
    while (nextLiteral<restrictionLiterals.N) {
      if (DEBUG>1) {
        cout<<"Inserting ";restrictionLiterals(nextLiteral)->write(cout);cout<<"   into   ";
        write(rulesC_old.rules.elem(r)->context);cout<<endl;
      }
      newRule = new Rule;
      newRule->action = rulesC_old.rules.elem(r)->action;
      newRule->context = rulesC_old.rules.elem(r)->context;
      newRule->insertContext(restrictionLiterals(nextLiteral));
      StateTransitionL covered_experiences;
      uintA covered_experiences_ids;
      arr covered_experience_weights;
      learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
      if (covered_experiences.N > 0) {
        if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences "<<covered_experiences_ids<<" and will be kept."<<endl;
        if (DEBUG>3) {
          cout<<"Covered experiences:"<<endl;
          uint k;
          FOR1D(covered_experiences, k) {
            covered_experiences(k)->write(cout);
          }
        }
        MT::Array< uintA > experiences_per_outcome;
        learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
        rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
      }
      else {
        if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
        delete newRule;
      }
      nextLiteral++;
      if (rulesC_2add.rules.num() > 0)
        break;
    }
    if (rulesC_2add.rules.num() > 0)
      break;
    else {
      nextRule = r+1;
      nextLiteral = 0;
      restrictionLiterals.clear();
    }
  }
  if (DEBUG>0) {if (rulesC_2add.rules.num() > 0) {newRule->write(cout);}}
  if (DEBUG>0) cout<<"AddReferences::findRules [END]"<<endl;
}

void AddReferences::reset() {
  nextRule = 1;
  nextLiteral = 0;
}















// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    G E N E R A L I Z E   E Q U A L I T Y
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

GeneralizeEquality::GeneralizeEquality() : SearchOperator() {
  name = "GeneralizeEquality";
  nextRule=1; // ignore default rule
  nextLiteral=0;
  doneLess=false;
}

// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void GeneralizeEquality::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"GeneralizeEquality::findRules [START]"<<endl;
  uint r, p, i;
  Rule* newRule;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    for (p=nextLiteral; p<rulesC_old.rules.elem(r)->context.N; p++) {
      if (rulesC_old.rules.elem(r)->context(p)->s->range_type != Symbol::binary
                &&  rulesC_old.rules.elem(r)->context(p)->comparison_type == Literal::comparison_equal) {
        if (DEBUG>1) {cout<<"Generalizing literal #"<<p<<": ";rulesC_old.rules.elem(r)->context(p)->write(cout);cout<<endl;}
        newRule = new Rule;
        newRule->action = rulesC_old.rules.elem(r)->action;
        FOR1D(rulesC_old.rules.elem(r)->context, i) {
          if (i!=p) {
            newRule->context.append(rulesC_old.rules.elem(r)->context(i));
          }
          else {
            Literal* lit_comp__old = rulesC_old.rules.elem(r)->context(p);
            Literal* lit_comp__new = NULL;
//             create new context literal
            if (!doneLess) {
              doneLess = true; // try greater-equal thereafter
              nextLiteral = p;
              lit_comp__new = Literal::get(lit_comp__old->s, lit_comp__old->args, lit_comp__old->value, Literal::comparison_lessEqual);
            }
            else {
              doneLess=false;
              nextLiteral = p+1; // try next literal thereafter
              lit_comp__new = Literal::get(lit_comp__old->s, lit_comp__old->args, lit_comp__old->value, Literal::comparison_greaterEqual);
            }
            if (DEBUG>1) {cout<<"Generalized literal: "<<*lit_comp__new<<endl;}
            newRule->context.append(lit_comp__new);
          }
        }
        Literal::sort(newRule->context);
        if (DEBUG>0) {cout<<"Potential new rule:"<<endl;newRule->write(cout);}
        StateTransitionL covered_experiences;
        uintA covered_experiences_ids;
        arr covered_experience_weights;
        learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
        if (covered_experiences.N > 0) {
          if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
          if (DEBUG>3) {
            cout<<"Covered experiences:"<<endl;
            uint k;
            FOR1D(covered_experiences, k) {
              covered_experiences(k)->write(cout);
            }
          }
          MT::Array< uintA > experiences_per_outcome;
          learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
          rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
          break;
        }
        else {
          delete newRule;
        }
      }
    }
    if (rulesC_2add.rules.num()>0) {
      break;
    }
    else {
      nextRule=r+1;
      nextLiteral=0;
      doneLess=false;
    }
  }
  if (DEBUG>0) cout<<"GeneralizeEquality::findRules [END]"<<endl;
}


void GeneralizeEquality::reset() {
    nextRule = 1;
    nextLiteral = 0;
    doneLess= false;
}












// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    S P L I T   O N   E Q U A L I T I E S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------


SplitOnEqualities::SplitOnEqualities() : SearchOperator() {
  name = "SplitOnEqualities";
  nextRule=1; // ignore default rule
  nextVar=0;
  nextFunc=0;
  SymL function_symbols;
  Symbol::get_state_nonBinary(function_symbols) ;
  usedFunctions.append(function_symbols);
  // Achtung DON'T used derived: usedFunctions.append(logicObjectManager::f_derived);!!!!
}

void SplitOnEqualities::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"SplitOnEqualities::findRules [START]"<<endl;
  rulesC_2add.clear();
  uint r, v, f, i;
  
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    if (vars.N == 0) { // first round
      // determine vars
      rulesC_old.rules.elem(r)->getArguments(vars);
      if (DEBUG>2) {
        cout << "Inspecting rule:"<<endl;
        rulesC_old.rules.elem(r)->write(cout);
        cout<<"findValue the following vars: "<<vars<<endl;
      }
    }
    for (v=nextVar; v<vars.d0; v++) {
      for (f=nextFunc; f<usedFunctions.d0; f++) {
        if (usedFunctions(f)->arity != 1)  // nur fuer unary functions!!
          continue;
        if (DEBUG>1) {cout<<"Checking var "<<v<<" with function "<<usedFunctions(f)->name<<endl;}
        // check whether function has not been used with this variable
        bool alreadyUsed = false;
        FOR1D(rulesC_old.rules.elem(r)->context, i) {
          if (rulesC_old.rules.elem(r)->context(i)->s == usedFunctions(f)) {
            if (rulesC_old.rules.elem(r)->context(i)->args.findValue(vars(v)) >= 0)
              alreadyUsed = true;
          }
        }
        if (DEBUG>1) { PRINT(alreadyUsed) }
        if (!alreadyUsed) {
          uintA tosplit_args;
          tosplit_args.append(vars(v));

          arr usedValues = (*usedFVs)[usedFunctions(f)];
          if (DEBUG>2) {cout<<"Using function values: "<<usedValues<<endl;}
          FOR1D(usedValues, i) {
            Rule* newRule = new Rule;
            newRule->action = rulesC_old.rules.elem(r)->action;
            newRule->context = rulesC_old.rules.elem(r)->context;
            Literal* lit_tosplit = Literal::get(usedFunctions(f), tosplit_args, usedValues(i), Literal::comparison_equal);
            newRule->insertContext(lit_tosplit);
            if (DEBUG>0) {cout<<"Potential new rule:"<<endl<<*newRule;}
            StateTransitionL covered_experiences;
            uintA covered_experiences_ids;
            arr covered_experience_weights;
            learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
            if (covered_experiences.N > 0) {
              if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
              MT::Array< uintA > experiences_per_outcome;
              learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
              rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
            }
            else {
              if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
              delete newRule;
            }
          }
          if (rulesC_2add.rules.num() > 0) {
            nextFunc = f+1;
            break;
          }

        }
      }
      if (rulesC_2add.rules.num() > 0) {
        break;
      }
      else {
        nextFunc = 0;
        nextVar = v+1;
      }
    }
    if (rulesC_2add.rules.num() > 0)
      break;
    else {
      nextRule = r+1;
      nextVar = 0;
      nextFunc = 0;
      vars.clear();
    }
  }
  if (DEBUG>0) {if (rulesC_2add.rules.num() > 0) {write(rulesC_2add.rules);}}
  if (DEBUG>0) cout<<"SplitOnEqualities::findRules [END]"<<endl;
}

void SplitOnEqualities::reset() {
  nextRule = 1;
  nextVar = 0;
  nextFunc = 0;
}


// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    S P L I T   O N   I N E Q U A L I T I E S
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------


SplitOnInequalities::SplitOnInequalities() : SearchOperator() {
  name = "SplitOnInequalities";
  nextRule=1; // ignore default rule
  nextVar=0;
  nextFunc=0;
  nextValue = 0;
  SymL function_symbols;
  Symbol::get_state_nonBinary(function_symbols) ;
  usedFunctions.append(function_symbols);
  // Achtung DON'T used derived: usedFunctions.append(logicObjectManager::f_derived);!!!!
}

void SplitOnInequalities::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"SplitOnInequalities::findRules [START]"<<endl;
  rulesC_2add.clear();
  uint i;
  
  for (; nextRule<rulesC_old.rules.num(); nextRule++) {
    if (vars.N == 0) { // first round
      // determine vars
      rulesC_old.rules.elem(nextRule)->getArguments(vars);
      if (DEBUG>2) {
        cout << "Inspecting rule:"<<endl;
        rulesC_old.rules.elem(nextRule)->write(cout);
        cout<<"findValue the following vars: "<<vars<<endl;
      }
    }
    for (; nextVar<vars.d0; nextVar++) {
      for (; nextFunc<usedFunctions.d0; nextFunc++) {
        if (usedFunctions(nextFunc)->arity != 1)  // nur fuer unary functions!!
          continue;
        if (DEBUG>1) {cout<<"Checking var "<<nextVar<<" with function "<<usedFunctions(nextFunc)->name<<endl;}
        // check whether function has not been used with this variable
        bool alreadyUsed = false;
        FOR1D(rulesC_old.rules.elem(nextRule)->context, i) {
          if (rulesC_old.rules.elem(nextRule)->context(i)->s == usedFunctions(nextFunc)) {
            if (rulesC_old.rules.elem(nextRule)->context(i)->args.findValue(vars(nextVar)) >= 0)
              alreadyUsed = true;
          }
        }
        if (DEBUG>1) { PRINT(alreadyUsed) }
        if (!alreadyUsed) {
          uintA tosplit_args;
          tosplit_args.append(vars(nextVar));

          arr usedValues = (*usedFVs)[usedFunctions(nextFunc)];
          if (DEBUG>2) {cout<<"Using function values: "<<usedValues<<endl;}
          for (; nextValue < usedValues.d0; nextValue++) {
            //Rule with <=
            Rule* newRule = new Rule;
            newRule->action = rulesC_old.rules.elem(nextRule)->action;
            newRule->context = rulesC_old.rules.elem(nextRule)->context;
            Literal* lit_tosplit = Literal::get(usedFunctions(nextFunc), tosplit_args, usedValues(nextValue), Literal::comparison_lessEqual);
            newRule->insertContext(lit_tosplit);
            if (DEBUG>0) {cout<<"Potential new rule:"<<endl<<*newRule;}
            StateTransitionL covered_experiences;
            uintA covered_experiences_ids;
            arr covered_experience_weights;
            learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
            if (covered_experiences.N > 0) {
              if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
              MT::Array< uintA > experiences_per_outcome;
              learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
              rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
            }
            else {
              if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
              delete newRule;
            }

            //Rule with >
            newRule = new Rule;
            newRule->action = rulesC_old.rules.elem(nextRule)->action;
            newRule->context = rulesC_old.rules.elem(nextRule)->context;
            lit_tosplit = Literal::get(usedFunctions(nextFunc), tosplit_args, usedValues(nextValue), Literal::comparison_greater);
            newRule->insertContext(lit_tosplit);
            if (DEBUG>0) {cout<<"Potential new rule:"<<endl<<*newRule;}
            covered_experiences.clear();
            covered_experiences_ids.clear();
            covered_experience_weights.clear();
            learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
            if (covered_experiences.N > 0) {
              if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
              MT::Array< uintA > experiences_per_outcome;
              learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
              rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
            }
            else {
              if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
              delete newRule;
            }

            if (rulesC_2add.rules.num() > 0) {
              nextValue++;
              if (nextValue == usedValues.d0) {
                nextValue = 0;
                nextFunc++;
                if (nextFunc == usedFunctions.d0) {
                  nextFunc = 0;
                  nextVar++;
                  if (nextVar == vars.d0) {
                    nextVar = 0;
                    vars.clear();
                    nextRule++;
                  }
                }
              }
              if (DEBUG>0) { write(rulesC_2add.rules); }
              return;
            }
          }
          nextValue = 0;
        }
      }
      nextFunc = 0;
    }
    vars.clear();
    nextVar = 0;
  }
  if (DEBUG>0) cout<<"SplitOnInequalities::findRules [END]"<<endl;
}

void SplitOnInequalities::reset() {
  nextRule = 1;
  nextVar = 0;
  nextFunc = 0;
  nextValue = 0;
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    ADD ABSTRACT EQUALITY
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------


AddAbstractEquality::AddAbstractEquality() : SearchOperator() {
  name = "AddAbstractEquality";
  nextRule=1; // ignore default rule
  nextVar=0;
  nextFunc=0;
  SymL function_symbols;
  Symbol::get_state_nonBinary(function_symbols) ;
  usedFunctions.append(function_symbols);
}

void AddAbstractEquality::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"AddAbstractEquality::findRules [START]"<<endl;
  rulesC_2add.clear();
  uint r, v, f, i;
  
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    if (vars.N == 0) { // first round
      // determine vars
      rulesC_old.rules.elem(r)->getArguments(vars);
      if (DEBUG>2) {
        cout << "Inspecting rule:"<<endl;
        rulesC_old.rules.elem(r)->write(cout);
        cout<<"findValue the following vars: "<<vars<<endl;
      }
    }
    for (v=nextVar; v<vars.d0; v++) {
      for (f=nextFunc; f<usedFunctions.d0; f++) {
        if (usedFunctions(f)->arity != 1)  // nur fuer unary functions!!
          continue;
        if (DEBUG>1) {cout<<"Checking var "<<v<<" with function "<<usedFunctions(f)->name<<endl;}
        // check whether function has not been used with this variable
        bool alreadyUsed = false;
        FOR1D(rulesC_old.rules.elem(r)->context, i) {
          if (rulesC_old.rules.elem(r)->context(i)->s == usedFunctions(f)) {
            if (rulesC_old.rules.elem(r)->context(i)->args.findValue(vars(v)) >= 0)
              alreadyUsed = true;
          }
        }
        if (DEBUG>1) { PRINT(alreadyUsed) }
        if (!alreadyUsed) {
          uintA newLitArgs;
          newLitArgs.append(vars(v));

          Rule* newRule = new Rule;
          newRule->action = rulesC_old.rules.elem(r)->action;
          newRule->context = rulesC_old.rules.elem(r)->context;

          Literal *newLit = Literal::getVarComparison(usedFunctions(f), newLitArgs);
          newRule->insertContext(newLit);

          if (DEBUG>0) {cout<<"Potential new rule:"<<endl<<*newRule;}
          StateTransitionL covered_experiences;
          uintA covered_experiences_ids;
          arr covered_experience_weights;
          learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
          if (covered_experiences.N > 0) {      
            MT::Array< uintA > experiences_per_outcome;
            learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
            if (newRule->existsInOutcome(newLit->s, newLitArgs)) {    //check whether newLit was used to predict something in the outcome
              if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;                 
              rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);   
              nextFunc = f+1;
              break;
            }
            else if (DEBUG>1) cout<<"Nothing useful was predicted, rule will be dropped."<<endl;
          }
          if (DEBUG>1) cout<<"Covers 0 experiences and will be dropped."<<endl;
          delete newRule;

        }
      }
      if (rulesC_2add.rules.num() > 0)
        break;
      else {
        nextFunc = 0;
        nextVar = v+1;
      }
    }
    if (rulesC_2add.rules.num() > 0)
      break;
    else {
      nextRule = r+1;
      nextVar = 0;
      nextFunc = 0;
      vars.clear();
    }
  }
  if (DEBUG>0) {if (rulesC_2add.rules.num() > 0) {write(rulesC_2add.rules);}}
  if (DEBUG>0) cout<<"AddAbstractEquality::findRules [END]"<<endl;
}

void AddAbstractEquality::reset() {
  nextRule = 1;
  nextVar = 0;
  nextFunc = 0;
}



// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    C H A N G E   R A N G E
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

ChangeRange::ChangeRange() : SearchOperator() {
    name = "ChangeRange";
    nextRule=1; // ignore default rule
    nextLiteral=0;
    nextPossibleValue=0;
}

// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
// only for CONSTANT-bound comparison predicates
void ChangeRange::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"ChangeRange::findRules [START]"<<endl;
  uint r, p, v;
  Rule* newRule;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    for (p=nextLiteral; p<rulesC_old.rules.elem(r)->context.N; p++) {
      Literal* clit = rulesC_old.rules.elem(r)->context(p);
      if (clit->s->range_type == Symbol::integers || clit->s->range_type == Symbol::integer_set) {
        if (DEBUG>1) {cout<<"Changing range of literal #"<<p<<": ";clit->write(cout);cout<<endl;}
        if (nextPossibleValue == 0) {
          // collect possibleValues 
          possibleValues = (*usedFVs)[clit->s];
          possibleValues.removeValueSafe(clit->value);
        }
        if (DEBUG>2) {cout<<"Using values: "<<possibleValues<<endl;}
        
        for (v=nextPossibleValue; v<possibleValues.d0; v++) {
          newRule = new Rule;
          newRule->action = rulesC_old.rules.elem(r)->action;
          uint n;
          FOR1D(rulesC_old.rules.elem(r)->context, n) {
            if (n!=p) {
              newRule->context.append(rulesC_old.rules.elem(r)->context(n));
            }
            else {
              // create new context literal
              Literal* eq = Literal::get(clit->s, clit->args, possibleValues(v), clit->comparison_type);
              newRule->insertContext(eq);
            }
          }
          if (DEBUG>0) {cout<<"Potential new rule:"<<endl;newRule->write(cout);}
          StateTransitionL covered_experiences;
          uintA covered_experiences_ids;
          arr covered_experience_weights;
          learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
          if (covered_experiences.N > 0) {
            if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
            if (DEBUG>3) {
              cout<<"Covered experiences:"<<endl;
              uint k;
              FOR1D(covered_experiences, k) {
                covered_experiences(k)->write(cout);
              }
            }
            MT::Array< uintA > experiences_per_outcome;
            learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
            rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
            nextPossibleValue = v+1;
            break;
          }
          else {
            if (DEBUG>1) cout<<"Covers 0 experiences and will be discarded."<<endl;
            delete newRule;
          }
        } // for values
      }
      if (rulesC_2add.rules.num()>0) {
        break;
      }
      else {
        nextLiteral=p+1;
        nextPossibleValue = 0;
      }
    } // for literals
    if (rulesC_2add.rules.num()>0) {
      nextLiteral = p;
      break;
    }
    else {
      nextRule=r+1;
      nextLiteral=0;
    }
  } // for rules
  if (DEBUG>0) cout<<"ChangeRange::findRules [END]"<<endl;
}


void ChangeRange::reset() {
  nextRule=1; // ignore default rule
  nextLiteral=0;
  nextPossibleValue=0;
}





// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    M A K E   I N T E R V A L
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

MakeInterval::MakeInterval() : SearchOperator() {
  name = "MakeInterval";
  nextRule=1; // ignore default rule
  nextLiteral=0;
  nextPossibleValue=0;
}

// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void MakeInterval::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"MakeInterval::findRules [START]"<<endl;
  uint r, p, v;
  Rule* newRule;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    for (p=nextLiteral; p<rulesC_old.rules.elem(r)->context.N; p++) {
      Literal *clit = rulesC_old.rules.elem(r)->context(p);
      if (clit->s->range_type == Symbol::integer_set || clit->s->range_type == Symbol::integers) {
        // TODO CHECK FOR DOUBLE INTERVALS!!!!!!!!!!!!
        // omit equality literals
        if (!(clit->comparison_type == Literal::comparison_less || clit->comparison_type == Literal::comparison_lessEqual
          || clit->comparison_type == Literal::comparison_greater || clit->comparison_type == Literal::comparison_greaterEqual))
          continue;
        if (DEBUG>1) {cout<<"Making interval for literal #"<<p<<": ";clit->write(cout);cout<<endl;}
        if (nextPossibleValue == 0) {
          // collect possibleValues 
          possibleValues = (*usedFVs)[clit->s];
          if (clit->comparison_type == Literal::comparison_less 
              || clit->comparison_type == Literal::comparison_lessEqual) {
            uint pv;
            FOR1D_DOWN(possibleValues, pv) {
              if (possibleValues(pv) >= clit->value)
                possibleValues.remove(pv);
            }
          }
          else if (clit->comparison_type == Literal::comparison_greater 
                    || clit->comparison_type == Literal::comparison_greaterEqual) {
            uint pv;
            FOR1D_DOWN(possibleValues, pv) {
              if (possibleValues(pv) <= clit->value)
                possibleValues.remove(pv);
            }
          }
          else {HALT("Don't know this comparison type, digger.");}
        }
        if (DEBUG>2) {cout<<"Using values: "<<possibleValues<<endl;}
        
        for (v=nextPossibleValue; v<possibleValues.d0; v++) {
          newRule = new Rule;
          newRule->action = rulesC_old.rules.elem(r)->action;
          newRule->context = rulesC_old.rules.elem(r)->context;
          
          // create new context literal
          Literal::ComparisonType newCompType;
          if (clit->comparison_type == Literal::comparison_less 
              || clit->comparison_type == Literal::comparison_lessEqual) {
            newCompType = Literal::comparison_greaterEqual;
          }
          else if (clit->comparison_type == Literal::comparison_greater 
                      || clit->comparison_type == Literal::comparison_greaterEqual) {
            newCompType = Literal::comparison_lessEqual;
          }
          else
              HALT("Don't know this comparison type, digger.")
          Literal* new_clit = Literal::get(clit->s, clit->args, possibleValues(v), newCompType);
          newRule->insertContext(new_clit);

          if (DEBUG>0) {cout<<"Potential new rule:"<<endl;newRule->write(cout);}
          StateTransitionL covered_experiences;
          uintA covered_experiences_ids;
          arr covered_experience_weights;
          learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
          if (covered_experiences.N > 0) {
            if (DEBUG>1) cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
            if (DEBUG>3) {
                cout<<"Covered experiences:"<<endl;
                uint k;
                FOR1D(covered_experiences, k) {
                    covered_experiences(k)->write(cout);
                }
            }
            MT::Array< uintA > experiences_per_outcome;
            learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
            rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
            nextPossibleValue = v+1;
            break;
          }
          else {
            if (DEBUG>1) cout<<"Covers 0 experiences and will be discarded."<<endl;
            delete newRule;
          }
        } // for values
      }
      if (rulesC_2add.rules.num()>0) {
          break;
      }
      else {
          nextLiteral=p+1;
          nextPossibleValue = 0;
      }
    } // for literals
    if (rulesC_2add.rules.num()>0) {
      nextLiteral = p;
      break;
    }
    else {
      nextRule=r+1;
      nextLiteral=0;
    }
  } // for rules
  if (DEBUG>0) cout<<"MakeInterval::findRules [END]"<<endl;
}


void MakeInterval::reset() {
  nextRule=1; // ignore default rule
  nextLiteral=0;
  nextPossibleValue=0;
}



// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    CompareFunctionValues
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

CompareFunctionValues::CompareFunctionValues() : SearchOperator() {
  name = "CompareFunctionValues";
  nextRule=1; // ignore default rule
  nextFunction=0;
  nextComparisonType=0;
  nextTermCombination=0;
  SymL function_symbols;
  Symbol::get_state_nonBinary(function_symbols) ;
  usedFunctions.append(function_symbols);
  // collect comparison types SAVE THE ORDER!
  comparisonTypes.append(Literal::comparison_equal);comparisonTypes.append(Literal::comparison_less);comparisonTypes.append(Literal::comparison_lessEqual);
  comparisonTypes.append(Literal::comparison_greater);comparisonTypes.append(Literal::comparison_greaterEqual);
}


// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void CompareFunctionValues::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"CompareFunctionValues::findRules [START]"<<endl;
#if 0
  uint r, f, c, t;
  Rule* newRule;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    for (f=nextFunction; f<usedFunctions.d0; f++) {
      // calc new term combos if needed
      if (nextTermCombination==0) {
        termCombos.clear();
        uintA terms;
        rulesC_old.rules.elem(r)->getArguments(terms);
        MT::Array< uintA > termCombos_unfiltered;
        TL::allPermutations(termCombos_unfiltered, terms, 2 * usedFunctions(f)->arity, false, true);
        // filter such that only terms are only combined once (i.e. [X,Y] and not also [Y,X]
        uint i,j;
        FOR1D(termCombos_unfiltered,i) {
          FOR1D(termCombos_unfiltered(i), j) {
            if (j==0) continue;
            if (termCombos_unfiltered(i)(j-1) > termCombos_unfiltered(i)(j))
              break;
          }
          if (termCombos_unfiltered(i).N == j)
            termCombos.append(termCombos_unfiltered(i));
        }
        if (DEBUG>2) cout<<"termCombos for function "<<usedFunctions(f)->name<<": "<<termCombos<<"   (unfiltered: "<<termCombos_unfiltered<<")"<<endl;
      }
      for (t=nextTermCombination; t<termCombos.d0; t++) {
        for (c=nextComparisonType; c<comparisonTypes.d0; c++) {
          if (c==0) {
              coveredExIDsPerComparisonType.clear();
              if (DEBUG>1) {cout<<"Extending OLD RULE:"<<endl;rulesC_old.rules.elem(r)->write(cout);}
          }
          // FILTER FOR SENSIBLE COMPARISONS THAT DO NOT CONTRADICT OTHER STUFF?
          // add comparison to rule
          newRule = new Rule;
          newRule->action = rulesC_old.rules.elem(r)->action;
          newRule->context = rulesC_old.rules.elem(r)->context;
          uintA args1, args2;
          uint u;
          for (u=0; u<termCombos(t).N/2; u++) args1.append(termCombos(t)(u));
          for (u=termCombos(t).N/2; u<termCombos(t).N; u++) args2.append(termCombos(t)(u));
          ComparisonLiteral* new_clit = logicObjectManager::getCompLiteral_dynamic(usedFunctions(f), comparisonTypes(c), args1, args2);
          reason::insert(*newRule, *new_clit);
          if (DEBUG>1) {cout<<"Potential new rule:"<<endl;newRule->write(cout);}
          StateTransitionL covered_experiences;
          uintA covered_experiences_ids;
          calcCoverage(covered_experiences, covered_experiences_ids, newRule, experiences);
          if (covered_experiences.N > 0) {
            if (DEBUG>1) {
//               cout<<"Covers "<<coveredExperiences.N<<" experiences and will be kept."<<endl;
//               uintA coveredExperiencesIDs;
//               calcCoverage(newRule, experiences, coveredExperiencesIDs);
//               coveredExIDsPerComparisonType.append(coveredExperiencesIDs);
//               cout<<"Example IDs: "<<coveredExperiencesIDs<<endl;
//               if (DEBUG>3) {
//                 cout<<"Covered experiences:"<<endl;
//                 uint k;
//                 FOR1D(coveredExperiences, k) {
//                   cout<<coveredExperiencesIDs(k)<<":"<<endl;
//                   coveredExperiences(k)->write(cout);
//                 }
//               }
            }
            MT::Array< uintA > experiences_per_outcome;
            learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, alpha_PEN, p_min, pen_sum, pen_pos);
            rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
            break;
          }
          else {
            delete newRule;
            if (DEBUG>1) cout<<"Covers 0 experiences and will be discarded."<<endl;
          }
        } // for comparison types
        if (rulesC_2add.rules.num()>0) {
          nextComparisonType = c+1;
          break;
        }
        else {
          nextComparisonType=0;
          if (DEBUG>2) {
//             cout<<"Recap on coverages:"<<endl;
//             StateTransitionL coveredExperiences;
//             calcCoverage(rulesC_old.rules.elem(r), experiences, coveredExperiences);
//             cout<<"Old rule: "<<endl;rulesC_old.rules.elem(r)->write(cout);cout<<" covers "<<coveredExperiences.N<<endl;
//             cout <<"while the new rules with comparisons cover: "<<endl;
//             uint xx;
//             FOR1D(coveredExIDsPerComparisonType, xx) {
//               cout<<"comptype #"<<xx<<"="<<comparisonTypes(xx)<<": "<<coveredExIDsPerComparisonType(xx).d0<<" "<<coveredExIDsPerComparisonType(xx)<<endl;
//             }
            // several tests for correct calculation
            // These tests are intuitive, but wrong!!
//                         CHECK(coveredExIDsPerComparisonType.d0==comparisonTypes.N, "Too many new rules!!!")
//                         CHECK(coveredExIDsPerComparisonType(0).d0+coveredExIDsPerComparisonType(1).d0==coveredExIDsPerComparisonType(2).d0, "equal + less != less-equal!!")
//                         CHECK(coveredExIDsPerComparisonType(0).d0+coveredExIDsPerComparisonType(3).d0==coveredExIDsPerComparisonType(4).d0, "equal + less != less-equal!!")
//                         CHECK(numberSharedElements(coveredExIDsPerComparisonType(0), coveredExIDsPerComparisonType(2)) == coveredExIDsPerComparisonType(0).N, "<= does not subsume ==!")
//                         CHECK(numberSharedElements(coveredExIDsPerComparisonType(0), coveredExIDsPerComparisonType(4)) == coveredExIDsPerComparisonType(0).N, ">= does not subsume ==!")
//                         CHECK(numberSharedElements(coveredExIDsPerComparisonType(1), coveredExIDsPerComparisonType(3)) == 0, "< and > share elements!")
          }
        }
      } // for term combos
      if (rulesC_2add.rules.num()>0) {
        nextTermCombination = t;
        break;
      }
      else {
        nextTermCombination=0;
      }
    } // for functions
    if (rulesC_2add.rules.num()>0) {
        nextFunction = f;
        break;
    }
    else {
        nextRule=r+1;
        nextFunction=0;
    }
  }
#endif
  if (DEBUG>0) cout<<"CompareFunctionValues::findRules [END]"<<endl;
}


void CompareFunctionValues::reset() {
  nextRule=1; // ignore default rule
  nextFunction=0;
  nextComparisonType=0;
  nextTermCombination=0;
}




// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    SplitOnCompareFunctionValues
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

SplitOnCompareFunctionValues::SplitOnCompareFunctionValues() : SearchOperator() {
  name = "SplitOnCompareFunctionValues";
  nextRule=1; // ignore default rule
  nextFunction=0;
  nextTermCombination=0;
  SymL function_symbols;
  Symbol::get_state_nonBinary(function_symbols);
  // collect comparison types SAVE THE ORDER!
  comparisonTypes.append(Literal::comparison_equal);comparisonTypes.append(Literal::comparison_less);comparisonTypes.append(Literal::comparison_greater);
}

// creates possible new rules for the given rule-set
// newRules are potential additional rules which are all intended to be added to the SAME rule-set!
void SplitOnCompareFunctionValues::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rulesC_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"SplitOnCompareFunctionValues::findRules [START]"<<endl;
#if 0
  uint r, f, c, t;
  Rule* newRule;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    for (f=nextFunction; f<usedFunctions.d0; f++) {
      // calc new term combos if needed
      if (nextTermCombination==0) {
        termCombos.clear();
        uintA terms;
        rulesC_old.rules.elem(r)->getArguments(terms);
        MT::Array< uintA > termCombos_unfiltered;
        TL::allPermutations(termCombos_unfiltered, terms, 2 * usedFunctions(f)->arity, false, true);
        // filter such that only terms are only combined once (i.e. [X,Y] and not also [Y,X]
        uint i,j;
        FOR1D(termCombos_unfiltered,i) {
          FOR1D(termCombos_unfiltered(i), j) {
            if (j==0) continue;
            if (termCombos_unfiltered(i)(j-1) > termCombos_unfiltered(i)(j))
              break;
          }
          if (termCombos_unfiltered(i).N == j)
            termCombos.append(termCombos_unfiltered(i));
        }
        if (DEBUG>2) cout<<"termCombos for function "<<usedFunctions(f)->name<<": "<<termCombos<<"   (unfiltered: "<<termCombos_unfiltered<<")"<<endl;
      }
      for (t=nextTermCombination; t<termCombos.d0; t++) {
        if (DEBUG>1) {cout<<"Extending OLD RULE:"<<endl;rulesC_old.rules.elem(r)->write(cout);coveredExIDsPerComparisonType.clear();}
        FOR1D(comparisonTypes, c) {
          // add comparison to rule
          newRule = new Rule;
          newRule->action = rulesC_old.rules.elem(r)->action;
          newRule->context = rulesC_old.rules.elem(r)->context;
          uintA args1, args2;
          uint u;
          for (u=0; u<termCombos(t).N/2; u++) args1.append(termCombos(t)(u));
          for (u=termCombos(t).N/2; u<termCombos(t).N; u++) args2.append(termCombos(t)(u));
          ComparisonLiteral* new_clit = logicObjectManager::getCompLiteral_dynamic(usedFunctions(f), comparisonTypes(c), args1, args2);
          reason::insert(*newRule, *new_clit);
          if (DEBUG>1) {cout<<"Potential new rule:"<<endl;newRule->write(cout);}
          StateTransitionL covered_experiences;
          uintA covered_experiences_ids;
          calcCoverage(covered_experiences, covered_experiences_ids, newRule, experiences);
          if (covered_experiences.N > 0) {
            if (DEBUG>1) {
              cout<<"Covers "<<covered_experiences.N<<" experiences and will be kept."<<endl;
              coveredExIDsPerComparisonType.append(covered_experiences_ids);
              cout<<"Example IDs: "<<covered_experiences_ids<<endl;
              if (DEBUG>3) {
                cout<<"Covered experiences:"<<endl;
                uint k;
                FOR1D(covered_experiences, k) {
                  cout<<covered_experiences_ids(k)<<":"<<endl;
                  covered_experiences(k)->write(cout);
                }
              }
            }
            MT::Array< uintA > experiences_per_outcome;
            learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, alpha_PEN, p_min, pen_sum, pen_pos);
            rulesC_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
          }
          else {
            delete newRule;
            if (DEBUG>1) cout<<"Covers 0 experiences and will be discarded."<<endl;
          }
        } // for comparison types
        if (rulesC_2add.rules.num()>0) {
          break;
        }
      } // for term combos
      if (rulesC_2add.rules.num()>0) {
          nextTermCombination = t+1;
          break;
      }
      else {
          nextTermCombination=0;
      }
    } // for functions
    if (rulesC_2add.rules.num()>0) {
      nextFunction = f;
      break;
    }
    else {
      nextRule=r+1;
      nextFunction=0;
    }
  }
#endif
  
  if (DEBUG>0) cout<<"SplitOnCompareFunctionValues::findRules [END]"<<endl;
}


void SplitOnCompareFunctionValues::reset() {
  nextRule=1; // ignore default rule
  nextFunction=0;
  nextTermCombination=0;
}



// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//    AbstractEquality
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

AbstractEquality::AbstractEquality() : SearchOperator() {
  name = "AbstractEquality";
  nextRule = 1; // ignore default rule
  nextLiteral = 0;
}

void AbstractEquality::findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"AbstractEquality::findRules [START]"<<endl;
  uint r, p, i;
  Rule* newRule;
  for (r=nextRule; r<rulesC_old.rules.num(); r++) {
    for (p=nextLiteral; p<rulesC_old.rules.elem(r)->context.N; p++) {
      Literal *contextLit = rulesC_old.rules.elem(r)->context(p);
      if (contextLit->comparison_type == Literal::comparison_equal 
        && (contextLit->s->range_type == Symbol::integers || contextLit->s->range_type == Symbol::integer_set)
        && (rulesC_old.rules.elem(r)->existsInOutcome(contextLit->s, contextLit->args))) {
          if (DEBUG>1) {cout<<"Generalizing literal #" << p << ": " << *rulesC_old.rules.elem(r)->context(p) << " in rule:" << endl << *rulesC_old.rules.elem(r) << endl;}
          newRule = new Rule();
          newRule->action = rulesC_old.rules.elem(r)->action;
          Literal *varComp = NULL;
          FOR1D(rulesC_old.rules.elem(r)->context, i) {
            if (i!=p)
              newRule->insertContext(rulesC_old.rules.elem(r)->context(i));
            else { // create new comparison literal of the form f(X)=C
              varComp = Literal::getVarComparison(rulesC_old.rules.elem(r)->context(i)->s, rulesC_old.rules.elem(r)->context(i)->args);
              newRule->insertContext(varComp);
              if (DEBUG > 1) cout << "Inserting " << *varComp << endl;
            }
          }
          if (DEBUG > 1) cout << "Generalized context: " << newRule->context << endl;
          StateTransitionL covered_experiences;
          uintA covered_experiences_ids;
          arr covered_experience_weights;
          learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, newRule, experiences, experience_weights);
          if (covered_experiences.N > 0) {
            MT::Array< uintA > experiences_per_outcome;
            learn::learn_outcomes(newRule, experiences_per_outcome, covered_experiences, covered_experiences_ids, covered_experience_weights);
            if (newRule->existsInOutcome(varComp->s, varComp->args)) {
              if (DEBUG > 0) { cout << "Keep rule" << endl << *newRule << endl; }
              rules_2add.append(newRule, covered_experiences_ids, experiences_per_outcome);
              break;
            }
            else cout << "No prediction learned, rule will be dropped" << endl;
          }
          delete newRule;
      }
    }
    if (rules_2add.rules.num()>0) {
      nextLiteral = p+1;
      break;
    }
    else {
      nextRule=r+1;
      nextLiteral=0;
    }
  }
  if (DEBUG>0) cout<<"AbstractEquality::findRules [END]"<<endl;
}

void AbstractEquality::reset() {
    nextRule = 1;
    nextLiteral = 0;
}










}
