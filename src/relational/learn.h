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

#ifndef RELATIONAL_learn_h
#define RELATIONAL_learn_h

#include <relational/reason.h>




namespace relational {
  
/************************************************
 * 
 *     learn -- Major interface
 * 
 ************************************************/

class RuleSetContainer;

namespace learn {
  void learn_rules(RuleSetContainer& rulesC, StateTransitionL& experiences, bool incremental=false, const char* logfile = "learn.log", const char* verbose_logfile="verbose_learn.log", uint debug=2);
  void learn_rules(RuleSetContainer& rulesC, StateTransitionL& experiences, arr& experience_weights, bool incremental=false, const char* logfile = "learn.log", const char* verbose_logfile="verbose_learn.log", uint debug=2);

  /*
   *  Regularization penalty
   *     Heuristics for setting alpha:
   *     --  For non-noise outcomes: for alpha<0.1 it pays off to use 2 rules with p=1.0 each instead of
   *         1 rule with two outcomes with p=0.5 each (i.e., not to merge several experiences).
   *         (assumption: the complexity of rules is not more than about 10 literals)
   */
  void set_penalty(double alpha_PEN);
  /* 
   *  Noise outcome probabilities
   *    Heuristics:
   *    -- For strange experiences, which would require a separate rule with very, very many literals,
   *       it does *not* pay off to include this separate rule. Rather, this should be modeled by the noise
   *       outcome of a different rule (potentially the default rule).
   *    -- Given alpha=1.0 and p_min < 10e-5, it pays off to include a separate rule for 10 rule literals.
   *       --> Thus, it almost always pays off to include a separate rule for "normal" 
   *       experiences (i.e., which require a normal-sized rule)
   *
   */
  void set_p_min(double p_min);
  // separate p_min choice for the noisy default rule
  void set_p_min(double p_min, double p_min_noisyDefaultRule);
  
  
  double score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, std::ostream& out=std::cout, uint DEBUG=0);
  double score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, arr& experience_weights, std::ostream& out=std::cout, uint DEBUG=0);
  
  
  // Choice of search operators
  #define RULE_LEARNER__OP_CHOICE__LINEAR 1
  #define RULE_LEARNER__OP_CHOICE__RANDOM 2
  void set_ChoiceTypeSearchOperators(uint choice_type);
  // Consideration horizon of previous successful search operator applications:
  #define SEARCH_OP_CHOICE__PAST_HORIZON 20
  #define SEARCH_OP_CHOICE__PAST_WEIGHT 0.5

  //If NO_DEICTICREFS_BY_NONBINARY is specified deictic references must have a unique covering substitution without considering non-binary symbols.
  #define NO_DEICTICREFS_BY_NONBINARY
  
  // statistics
  extern uintA num_so_improvements;
  extern uintA num_so_applied;
  extern arr so_improvements;
  extern uintA so_successfulUsageHistory;
  extern uintA so_UsageHistory;
  extern arr scores;
}




/************************************************
 * 
 *     RuleSetContainer
 * 
 *     Efficiency wrapper for rule-sets
 *     --> stores the coverage of experiences
 * 
 ************************************************/

struct RuleSetContainer {
  // First rule = default rule
  // Don't set rules directly, but use append(.) below.
  RuleSet rules;
  const StateTransitionL* p_experiences;
  bool p_experiences_copied;
  //StateTransitionL* p_experiences;
  arr experience_weights;

  // use additional experiences to crossvalidate the rule set
  // -> provide sync across actions
  const StateTransitionL* p_experiences_crossvalidation;

  // redundant memories
  MT::Array< uintA > nonDefaultRules_per_experience;  // only non-default rules!
  MT::Array< uintA > experiences_per_rule;
  
  MT::Array< MT::Array < uintA > > experiences_per_ruleOutcome;

  bool rules_copied;

  RuleSetContainer(const StateTransitionL* _p_experiences, const arr&  experience_weights = arr());
  RuleSetContainer(); // use init() later when using this constructor
  RuleSetContainer& operator= (const RuleSetContainer& other);
  virtual void copyWithRules(const RuleSetContainer&);

  virtual ~RuleSetContainer();

  virtual void init(const StateTransitionL* _p_experiences, const arr& experience_weights = arr());
  virtual void append(Rule* rule, uintA& experiences_of_this_rule, MT::Array< uintA >& experiences_per_outcome_of_this_rule);
  virtual void remove(uint id);
  virtual void clear();
  virtual void recomputeDefaultRule();
  virtual void sort();

  virtual void addExperiences(const StateTransitionL* _p_experiences, const arr&  experience_weights=arr());

  virtual void getResponsibilities(arr& responsibilities, MT::Array< uintA >& covered_experiences, uintA& covered_experiences_num) const;
  virtual void getPartitionsForAction(MT::Array< uintA >& partititions, Literal* action) const;
  
  virtual void write(ostream& out = std::cout, bool only_action = false, bool additional_experience_info = true) const;
  virtual void write(const char* filename, bool only_action = false, bool additional_experience_info = true) const;
  virtual void write_experiencesWithRules(ostream& os = std::cout) const;
  virtual void write_rulesWithExperiences(ostream& os = std::cout) const;
  
  virtual void sanityCheck(bool ignore_default_rule = false) const;  // EXPENSIVE!!!
};
  


/************************************************
 * 
 *     learn -- details
 * 
 *     (only of interest for programmers)
 * 
 ************************************************/
  
namespace learn {
  void calcCoverage(StateTransitionL& covered_experiences, uintA& covered_experiences_ids, const Rule* r, const StateTransitionL& experiences, bool consider_non_unique_covering=false);
  void calcCoverage(StateTransitionL& covered_experiences, uintA& covered_experiences_ids, arr& covered_experience_weights, const Rule* r, const StateTransitionL& experiences, const arr& experience_weights, bool consider_non_unique_covering=false);
  
  // *** Outcome learning ***
  // Main method which is safe to call from outside (namely from the SearchOperators).
  void learn_outcomes(Rule* rule, MT::Array< uintA >& coveredExperiences_per_outcome, const StateTransitionL& covered_experiences, const uintA& covered_experiences_ids, const arr& covered_experience_weights, uint DEBUG=0);
  
  // *** Parameter learning ***
  // (outcome probabilities)  (For efficiency, it's interwoven with "learn_outcomes" in the implementation details.)
  double learn_parameters(const MT::Array< LitL >& outcomes, doubleA& probs, uint DEBUG = 0);
  // Penalties for adapting gradients to learn sound probability distribution
  // - pen_sum: probs sum to 1
  void setProbabilitiesLearningPenalty_sum(double pen_sum);
  // - pen_pos: probs non-negative
  void setProbabilitiesLearningPenalty_pos(double pen_pos);
  
  // *** Cost function ***
  // (for rule-set evaluation)
  namespace CostFunction {
    void setRuleCoveredExperiences(const StateTransitionL& coveredEx);
    void setOutcomesCoverage(const boolA& coverage);
    double loglikelihood(const arr& probs);
    
    // C = -LOG_LIK + PENALTIES
    // --> to minimize
    // need knowledge of covered experiences and corresponding covered outcomes
    double calc(const arr& in); 
    void calc_grad(arr& out, const arr& in);
  };
};






/************************************************
 * 
 *     SearchOperator
 * 
 ************************************************/



// Search operator weights determine how often each operator is tried.
// (The larger the weight, the more often.)
#define SO_WEIGHT__EXPLAIN_EXPERIENCES 4.0
#define SO_WEIGHT__EXPLAIN_EXPERIENCES_SLIM 4.0
// incl. dynamic bound cpts
#define SO_WEIGHT__EXPLAIN_EXPERIENCES_SLIM_AND_COMPARING 0.0  // 4.0

#define SO_WEIGHT__DROP_CONTEXT_LITERALS 0.0
#define SO_WEIGHT__DROP_CONTEXT_LITERALS_APPROX 100.0  // 48.0 vorher
#define APPROXIMATOR__RULES_PER_ROUND 5
#define SO_WEIGHT__DROP_REFS 2.0

#define SO_WEIGHT__DROP_RULES 3.0

#define SO_WEIGHT__SPLIT_ON_LITS 3.0
#define SO_WEIGHT__ADD_LITS 2.0
#define SO_WEIGHT__ADD_REFS 3.0

// constant bounds
#define SO_WEIGHT__SPLIT_ON_EQS 0.0 // 1.0
#define SO_WEIGHT__SPLIT_ON_INEQUALITIES 0.0 // 0.5
#define SO_WEIGHT__CHANGE_RANGE 0.0 // 1.0
#define SO_WEIGHT__MAKE_INTVL 0.0 // 1.0
// dynamic bounds
#define SO_WEIGHT__COMPARE_FUNCTIONVALUES 0.0
#define SO_WEIGHT__SPLIT_ON_COMPARE_FUNCTIONVALUES 0.0
// both bounds
#define SO_WEIGHT__GENERALIZE_EQS 0//3.0

#define SO_WEIGHT__ABSTRACT_EQS 2.0
#define SO_WEIGHT__ADD_ABSTRACT_EQS 3.0
#define SO_WEIGHT__ADD_REFS_AND_ADD_LITS 0//0.5
#define SO_WEIGHT__ADD_REFS_INDIRECT 0//0.1

#define SO_WEIGHT__SPLIT_ON_LIT_CONJUNCTIONS 0.5 //3.0

class SearchOperator {
  
protected:

  MT::String name;
  bool approximative;
        
  // if we do incremental learning new experiences are
  // appended at the end and search operators can ignore them
  static uint experience_idx_min;

  // takes the rulelist rules2add and integrates it into existing ruleset
  void integrateNewRules(const RuleSetContainer& rulesC_old, const RuleSetContainer& rules_2add,
                          const StateTransitionL& experiences, RuleSetContainer& rules_new);
          
public:
  SearchOperator();
  virtual ~SearchOperator() {}

  virtual void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, RuleSetContainer& rules_2add) {
    arr experience_weights(experiences.N);
    experience_weights.setUni(1);
    this->findRules(rulesC_old, experiences, experience_weights, rules_2add);
  }

  // Creates possible new rules for the given rule-set.
  // "rules_2add" are potential additional rules which are all supposed to become part of the SAME rule-set!
  // I.e., they will be integrated into the old rule-set.
  virtual void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add) = 0;

  // resets search operator for next application round
  virtual void reset() = 0;
  // central method which is called by the RuleLearner
  virtual void createRuleSets(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights,
                  MT::Array< RuleSetContainer >& sets_of_new_rules);
  const char* getName();
  bool isApproximator() {return approximative;}
  virtual void reset_total_approximator() {}

  static void setExperienceIndexMin(uint idx) { experience_idx_min = idx; }
  static void resetExperienceIndexMin() { experience_idx_min = 0; }
};





class ExplainExperiences : public SearchOperator {
  uint nextPotentialExperience;
  bool slimContext;
  bool comparingValues;
  public:
    ExplainExperiences(bool slimContext, bool comparingValues);
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    Rule* explainExperience(StateTransition* ex);
    Rule* explainExperience_straightforward(StateTransition* ex);
    Rule* explainExperience_deictic(StateTransition* ex);
    Rule* explainExperience_deictic_ALL_DRs(StateTransition* ex);
    // resets field "nextPotentialExperience
    void reset();

    static uint maxRelationalActionDistance; // <sebastian>
};


class DropContextLiterals : public SearchOperator {
  uint nextRule;
  uint nextContextLiteral;
  public:
    DropContextLiterals();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};


class DropContextLiterals_approximativeVersion : public SearchOperator {
  uintA usableContextLiterals;
  const static uint DROP_NEGATIVE_BIAS = 3;
  bool prepareTotalNewSearch;
  public:
    DropContextLiterals_approximativeVersion();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
    void reset_total_approximator();
};


class DropReferences : public SearchOperator {
  uint nextRule;
  uint nextReference;
  public:
    DropReferences();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};


class DropRules : public SearchOperator {
  public:
    DropRules();
    void createRuleSets(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights,
              MT::Array< RuleSetContainer >& sets_of_new_rules);
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};


// considers only symbols with binary ranges (= predicates)
class SplitOnLiterals : public SearchOperator {
  uint nextRule;
  uint nextLiteral;
  LitL absentLiterals;
  uint newVar;
  public:
    SplitOnLiterals();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();

    virtual void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, RuleSetContainer& rules_2add) {
      arr experience_weights(experiences.N);
      experience_weights.setUni(1);
      this->findRules(rulesC_old, experiences, experience_weights, rules_2add);
    }
};

// considers only symbols with binary ranges (= predicates)
class SplitOnLiteralConjunctions : public SearchOperator {
  uint nextRule;
  uint nextLiteral1;
  uint nextLiteral2;
  LitL absentLiterals;
  uint newVar;
  public:
    SplitOnLiteralConjunctions();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};

class AddLiterals : public SearchOperator {
  uint nextRule;
  uint nextLiteral;
  LitL absentLiterals;
  public:
    AddLiterals();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};


class AddReferences : public SearchOperator {
  uint nextRule;
  uint nextLiteral;
  LitL restrictionLiterals;
  public:
    AddReferences();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};


// ---------------------------------------------------------
// ---------------------------------------------------------
//  SearchOperators on functions

class GeneralizeEquality : public SearchOperator {
  uint nextRule;
  uint nextLiteral;
  bool doneLess;

  public:
    GeneralizeEquality();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};

class AbstractEquality : public SearchOperator {
  uint nextRule;
  uint nextLiteral;

  public:
    AbstractEquality();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};

class AddAbstractEquality : public SearchOperator {
  uint nextRule;
  uint nextVar;
  uint nextFunc;
  uintA vars;
  SymL usedFunctions;
  
  public:
      AddAbstractEquality();
      void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
      void reset();
};

class AddIndirectReferences: public SearchOperator {
  uint nextRule;
  uint nextRefLiteral;
  LitL restrictionRefLiterals;

  uint newVar;
  uint nextNewVarLit;
  LitL restrictionNewVarLiterals;

  public:
    AddIndirectReferences();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};

class AddReferencesAndAddLits: public SearchOperator {
  uint nextRule;
  uint nextRefLiteral;
  LitL restrictionRefLiterals;

  uint newVar;
  uint nextNewVarLit;
  LitL restrictionNewVarLiterals;

  public:
    AddReferencesAndAddLits();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};

// for each variable v for each function f for which a comparison predicate with v is not
// used yet, we introduce equality literals for all existing values for v.
class SplitOnEqualities : public SearchOperator {
  uint nextRule;
  uint nextVar;
  uint nextFunc;
  uintA vars;
  SymL usedFunctions;

  //cached used function values accross all experiences
  std::map<relational::Symbol*, MT::Array<double> > *usedFVs;
  
  public:
      SplitOnEqualities();
      void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
      void reset();

      void setUsedFunctionValues(std::map<relational::Symbol*, MT::Array<double> > &usedFunctionValues) { usedFVs = &usedFunctionValues; }
};

class SplitOnInequalities : public SearchOperator {
  uint nextRule;
  uint nextVar;
  uint nextFunc;
  uint nextValue;
  uintA vars;
  SymL usedFunctions;

  //cached used function values accross all experiences
  std::map<relational::Symbol*, MT::Array<double> > *usedFVs;
  
  public:
      SplitOnInequalities();
      void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
      void reset();

      void setUsedFunctionValues(std::map<relational::Symbol*, MT::Array<double> > &usedFunctionValues) { usedFVs = &usedFunctionValues; }
};

// changes ranges of all comp preds
class ChangeRange : public SearchOperator {
  uint nextRule;
  uint nextLiteral;
  uint nextPossibleValue;
  arr possibleValues;

  //cached used function values accross all experiences
  std::map<relational::Symbol*, MT::Array<double> > *usedFVs;

  public:
    ChangeRange();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();

    void setUsedFunctionValues(std::map<relational::Symbol*, MT::Array<double> > &usedFunctionValues) { usedFVs = &usedFunctionValues; }
};


// for all less/greater and less/greater-equal compPreds, we introduce
// less/greater-equal compPreds for all possible values
// only defined for constant comparisons!
class MakeInterval : public SearchOperator {
  uint nextRule;
  uint nextLiteral;
  uint nextPossibleValue;
  arr possibleValues;

  //cached used function values accross all experiences
  std::map<relational::Symbol*, MT::Array<double> > *usedFVs;

  public:
    MakeInterval();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();

    void setUsedFunctionValues(std::map<relational::Symbol*, MT::Array<double> > &usedFunctionValues) { usedFVs = &usedFunctionValues; }
};


// Compares values of the same function, applied to all possible term combinations.
// (1) for all functions (2) for all term combinations (3) for all comparison types
// example: (1) size(.),  (2) size(x1)Xsize(x2), size(x1)Xsize(x3)...,
//   (3) size(x1)<size(x2), size(x1)<=size(x2), size(x1)==size(x2)...
// [CPT]: ComparisonPredicates dynamic
// 
// Be aware of the following subtle behavior: cover(CPT with ==) and cover(CPT with <)
// does not necessarily imply cover(CPT with <=). Explanation: Rule with (CPT==) may cover with sub1
// and rule with (CPT<) may cover with sub2 (=/= sub1), so that rule with (CPT<=) could cover with
// sub1 and sub2. But since we are using deictic referencing, the rule with (CPT<=) 
// does not apply in case in that case since we have more than one possible substitution for the rule.
class CompareFunctionValues : public SearchOperator {
  uint nextRule;
  uint nextFunction;
  uint nextComparisonType;
  uint nextTermCombination;
  
  MT::Array<uintA> termCombos;
  SymL usedFunctions;
  
  MT::Array<uintA> coveredExIDsPerComparisonType; // only for debugging
  
  MT::Array< Literal::ComparisonType > comparisonTypes;

  public:
    CompareFunctionValues();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};


// Introduces a comparison between two functions values and defines 3 corresponding rules,
// namely with <,==,>, that are introduced to a single new rule-set.
//
// (1) for all functions (2) for all term combinations
// example: (1) size(.),  (2) size(x1)Xsize(x2), size(x1)Xsize(x3)...,
// [CPT]: ComparisonPredicates dynamic
// 
// Be aware of the following subtle behavior: cover(CPT with ==) and cover(CPT with <)
// does not necessarily imply cover(CPT with <=). Explanation: Rule with (CPT==) may cover with sub1
// and rule with (CPT<) may cover with sub2 (=/= sub1), so that rule with (CPT<=) could cover with
// sub1 and sub2. But since we are using deictic referencing, the rule with (CPT<=) 
// does not apply in case in that case since we have more than one possible substitution for the rule.
class SplitOnCompareFunctionValues : public SearchOperator {
  uint nextRule;
  uint nextFunction;
  uint nextTermCombination;
  MT::Array< Literal::ComparisonType > comparisonTypes;
  MT::Array<uintA> termCombos;
  SymL usedFunctions;
  
  MT::Array<uintA> coveredExIDsPerComparisonType; // only for debugging

  public:
    SplitOnCompareFunctionValues();
    void findRules(const RuleSetContainer& rulesC_old, const StateTransitionL& experiences, const arr& experience_weights, RuleSetContainer& rules_2add);
    void reset();
};



}  // namespace relational

#endif // RELATIONAL_learn_h
