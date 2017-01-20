/*  
    Copyright 2008-2012   Tobias Lang
    
    E-mail:    tobias.lang@fu-berlin.de
    
    This file is part of libPRADA.

    libPRADA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your olition) any later version.

    libPRADA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libPRADA.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef RELATIONAL_reason_h
#define RELATIONAL_reason_h

#include <relational/rules.h>


namespace relational {


namespace reason {
  
  // *********************************************************************************************
  // *********************************************************************************************
  //   Basic Reasoning
  // *********************************************************************************************
  // *********************************************************************************************
	
	
  /****************************************
    ABSTRACTING AND GROUNDING
  ***************************************/
  
  bool isConstant(uint id);
  bool isGround(const Literal* lit);
  bool isGround(const SymbolicState& state);
  bool isGround(const StateTransition& trans);
  bool isPurelyAbstract(const Literal* lit);
  
  void setConstants(uintA& constants);
  void setConstants(uintA& constants, const ArgumentTypeL& constants_types);
  uintA& getConstants();
  ArgumentType* getArgumentTypeOfObject(uint object_id);

  
  
  
  /****************************************
    CALCULATING DERIVED SYMBOLS
   ***************************************/
  
  bool derive_conjunction(LitL& lits_derived, ConjunctionSymbol& s, const LitL& lits_given, const uintA& constants);
  bool derive_transclosure(LitL& lits_derived, TransClosureSymbol& s, const LitL& lits_given);
  bool derive_count(LitL& lits_derived, CountSymbol& s, const LitL& lits_given, const uintA& constants);
  bool derive_avg(LitL& lits_derived, AverageFunction& s, const LitL& lits_given, const uintA& constants);
  bool derive_max(LitL& lits_derived, MaxFunction& s, const LitL& lits_given, const uintA& constants);
  bool derive_sum(LitL& lits_derived, SumFunction& s, const LitL& lits_given, const uintA& constants);
  bool derive_reward(LitL& lits_derived, RewardFunction& s, const LitL& lits_given, const uintA& constants);
  bool derive_functiondiff(LitL& lits_derived, DifferenceFunction& s, const LitL& lits_given, const uintA& constants);
    
  // If constants are empty, then constants are calculated from literals.
  void derive(LitL& lits_derived, const LitL& lits_given, const uintA& constants);
  void derive(SymbolicState* s);
  
  void dederive(SymbolicState* s);
  
    
    
  /****************************************
      HOLDS (--> for ground)
  ***************************************/

  // GOAL: Coverage for ground literals (-> simple existence check)
  // In contrast to "cover(.)" which covers abstract / unground literals (see below)
  bool consistent(Literal* lit_given, Literal* lit_test);
  bool holds(const LitL& lits_given, Literal* lit_test);
  bool holds(const LitL& lits_given, const LitL& lits_test);
	
  
 
  /****************************************
      COVERAGE (--> for abstract)
   ***************************************/
  
  // "calcSubstitution" unifies literals, that is, calculates substitutions mapping from one to the other.
  bool calcSubstitution(Substitution& sub, const uintA& ground_args, const uintA& other_args);
  bool calcSubstitution(Substitution& sub, const Literal* ground_lit, const Literal* other_lit);
  
  // GOAL: Coverage for abstract literal arrays
  
  // SEMANTICS
  // Free variables:
  // In "cover", free variables are dealt with depending upon the literal they are involved in:
  // (i)  Positive binary literals  and  non-binary literals:
  //      Existential quantification, i.e., we need to find a single substitution for the literal to hold
  // (ii) Negative binary literals:
  //      There are two possibilities specified by the "free_neg_vars__all_quantified"-flag.
  //      In either case, free variables in negative literals lead to a substantial computational burden
  //      and should be avoided. Therefore, we usually order Literal arrays such that
  //      the positive literals appear first (which bind many variables).
  //      (A) All-quantification:   "forall  free-variables X_i:  -p(X1, ..., Xn)"
  //          A negative literal -p(X1,...,Xn) is only true if p(X1, ..., Xn) does not hold for all
  //          possible substitutions of the free variables.
  //      (B) Existential-quantification:  "exists  free-variables X_i:  -p(X1, ..., Xn)"
  //          A negative literal -p(X1,...,Xn) is true if there is at least one
  //          substitution S={X1->c1,...,Xn-->cn} such that p(c1,...,cn) does not hold.
  
  // USAGE
  // Object instances of data-structure "Substitution":
  // - All substitutions in "subs" arrays found by "cover" are newly created Substition instances. It is
  //   the responsibility of the calling methods to delete them!
  // - If a cover method returns false, the "subs" array is supposed to be empty. (So, the caller does
  //   not have to care about deleting object instances.)
    
  // BASIC FUNCTION 1
	// Calculates substitutions "subs" (extensions of "initSub") which
	// substitute the arguments in "lit" such that the resulting
	// ground literal is satisfied in "state".
  bool calcSubstitutions(SubstitutionSet& subs, const LitL& lits_ground, Literal* lit_input, bool free_neg_vars__all_quantified, Substitution* initSub = NULL);
  // BASIC FUNCTION 2
	// Calculates substitutions "subs" (extensions of "initSub") which
	// substitute the arguments in the literal list "lits" such that the resulting
	// ground literals are satisfied in "state".
  bool calcSubstitutions(SubstitutionSet& subs, const LitL& lits_ground, const LitL& lits_input, bool free_neg_vars__all_quantified, Substitution* initSub = NULL);
  bool covers(const LitL& lits_ground, const LitL& lits_input, bool free_neg_vars__all_quantified, Substitution* initSub = NULL);
    
    
    
	
  /****************************************
            STATE UNIFICATION 
   ***************************************/
  
  bool calcSubstitutions(SubstitutionSet& subs, const SymbolicState& state1, const SymbolicState& state2, Substitution* initSub = NULL);
  bool unifiable(const SymbolicState& state1, const SymbolicState& state2);

  // this method makes sure that state1 and state2 have disjoint constants before
  // trying unification. this is helpful if you just care whether states are structurally
  // equal, not about the actual constants
  bool unifiableIgnoreCommonConstants(SubstitutionSet& subs, const SymbolicState& state1, const SymbolicState& state2);
  bool unifiableIgnoreCommonConstants(const SymbolicState& state1, const SymbolicState& state2);

  // returns min-difference in number of literals
  uint calcSubstitutions_asMuchAsPossible(SubstitutionSet& subs, const SymbolicState& state1, const Literal& action1, const SymbolicState& state2, const Literal& action2);
  uint calcSubstitutions_asMuchAsPossible(SubstitutionSet& subs, const SymbolicState& state1, const SymbolicState& state2, Substitution* initSub = NULL);

	
  
  
  // *********************************************************************************************
  // *********************************************************************************************
  //   Rule Reasoning
  // *********************************************************************************************
  // *********************************************************************************************
  
  
   /****************************************
     Grounding
   ***************************************/
  
  bool isGround(const Rule& rule);
  bool isPurelyAbstract(const Rule& rule);
  // we allow for unground free vars in negative context literals, e.g. X in -on(a,X)
  bool isGround_positives(Rule* rule);  
  void calcGroundDeicticReferences(uintA& ground_drefs, const SymbolicState& state, const Literal* groundAction, const Rule* rule);
  
  
  /****************************************
    TRANSITION KNOWLEDGE
   ***************************************/
  
  #define STATE_TRANSITION__NOISE_OUTCOME 11
  
  // Assumption: outcome contains primitives only
  void calcSuccessorState(SymbolicState& successor, const SymbolicState& precessor, const LitL& outcome);
  double sampleSuccessorState_groundRule(SymbolicState& successor, const SymbolicState& precessor, Rule* ground_rules, uint& flag);
  // returns true iff unique covering rule; otherwise returns false (--> undefined successor state)
  double sampleSuccessorState_groundRules(SymbolicState& successor, const SymbolicState& precessor, const RuleSet& ground_rules, Literal* action, uint& flag);
  double sampleSuccessorState_groundRules(SymbolicState& successor, const SymbolicState& precessor, const RuleSet& ground_rules, Literal* action);
  double sampleSuccessorState_abstractRules(SymbolicState& successor, const SymbolicState& precessor, const RuleSet& abstract_rules, Literal* action);
  
  // assumes that rule is indeed applicable and rule is ground
  double probability_groundRule(Rule* groundRule, const SymbolicState& pre, const SymbolicState& post, double noiseStateProbability = 0.0);
  // grounds rule first (checks whether context and action hold)
  double probability_abstractRule(Rule* abstractRule, const SymbolicState& pre, Literal* groundAction, const SymbolicState& post, double noiseStateProbability = 0.0, Substitution* sub = NULL);
  double probability_defaultRule(Rule* defaultRule, const SymbolicState& pre, const SymbolicState& post, double noiseStateProbability = 0.0);
  
  
  
  /****************************************
    STATE-ACTION COVERAGE
   ***************************************/
  
  // Calculates the substitutions "subs" of context variables
  // such that the context of "rule" are satisified in state "s".
  // Ensures that deictic references are different from action arguments!
  bool calcSubstitutions_context(SubstitutionSet& subs, const SymbolicState& s, const Rule* rule, Substitution* actionSub, bool consider_non_unique_covering=false);
   
  // RULE COVERAGE
  // action prescribed to be "groundAction"
  // (Calculates the substitutions "subs" which ground the _action_ and the context
  // of "rule" such that they are satisified in "action" and state "s".)
  bool calcSubstitutions_rule_groundAction(SubstitutionSet& subs, const SymbolicState& s, const Literal* groundAction, const Rule* rule, bool consider_non_unique_covering=false);
  // General version.
  // Returns all possible instantiations for this rule.
  // (Remember: action arguments are allowed
  // to have different substitutions, while deictic reference can only have one sub once the
  // subs for the action arguments are given.)
  bool calcSubstitutions_rule(SubstitutionSet& subs, const SymbolicState& s, const Rule* rule);
  
  
  // WRAPPERS FOR RULE-SETS
  
  // action prescribed to be "groundAction"
  void calc_coveringRules_groundAction(RuleSet& coveringGroundRules, const RuleSet& allAbstractRules, const SymbolicState& s, Literal* groundAction);
  void calc_coveringRules_groundAction(uintA& coveringRuleIDs, const RuleSet& allAbstractRules, const SymbolicState& s, Literal* groundAction);
  void calc_coveringGroundRules_groundAction(RuleSet& coveringGroundRules, const RuleSet& allGroundRules, const SymbolicState& s, Literal* groundAction);
  void calc_coveringGroundRules_groundAction(uintA& coveringRuleIDs, const RuleSet& allGroundRules, const SymbolicState& s, Literal* groundAction);
  uint calc_uniqueCoveringRule_groundRules_groundAction(const RuleSet& allGroundRules, const SymbolicState& s, Literal* groundAction);
  uint calc_uniqueAbstractCoveringRule_groundAction(const RuleSet& allAbstractRules, const SymbolicState& s, Literal* groundAction);
  
  // General version for all actions
  void calc_coveringRules(RuleSet& coveringGroundRules, const RuleSet& allAbstractRules, const SymbolicState& s);
  void calc_coveringGroundRules(RuleSet& coveringGroundRules, const RuleSet& allGroundRules, const SymbolicState& s);
  
  // For several actions
  void calc_coveringRules(uintA& coveringRulesIDs, const RuleSet& abstract_rules, const LitL& ground_actions, const SymbolicState& s);
  void calc_coveringGroundRules(uintA& coveringRulesIDs, const RuleSet& ground_rules, const LitL& ground_actions, const SymbolicState& s);
  
  // explaining state transitions
  void calc_coveringOutcomes(uintA& covering_outcomes, Rule* groundRule, const SymbolicState& pre, const SymbolicState& post);
  void calc_coveringOutcomes(uintA& covering_outcomes, Rule* abstractRule, const SymbolicState& pre, Literal* groundAction, const SymbolicState& post);
  
  
  /****************************************
    SCORING
   ***************************************/
  
  double log_likelihood(const RuleSet& rules, const StateTransitionL& experiences, double p_min);
  double score(const RuleSet& rules, const StateTransitionL& experiences, double p_min, double alpha_pen);


  
    //      MORE THOUGHTS W.R.T. DEICTIC REFERENCES IN RULES:
  //         The deep question of QUANTIFICATION with NEGATED FREE DEICTIC VARS
  //         Example rule:
  //            ACTION:  act(X)   CONTEXT:  p2(X), -p1(Y)
  //         Two possibilities:
  //           (i)  Representing all quantification:  \forall Y: -p1(Y)
  //                --> BAD BAAAAAAAAAAD  Don't do that!
  //                    This is like introducing a formula "for all constants, the respective predicate does not hold".
  //           (ii) Representing exists quantification:  \exists Y: -p1(Y)
  //                --> does not make intuitive sense
  //                    Example:   ACTION:  puton(X)    CONTEXT:  -on(Y,X)
  //                               Rule applicable if exactly one arbitrary object not on X
  //           Baseline: Should be forbidden in NID rules, in general.
  //                     However, if there comes such a case nonetheless, we use exists quantification.
  
};


}



#endif // RELATIONAL_reason_h
