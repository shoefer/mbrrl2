#ifndef RULELEARNER2_H_
#define RULELEARNER2_H_

#include "learn.h"
#include "score_functions.h"
#include "symbols.h"
#include "rules.h"
#include "literals.h"

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <relational/learn.h>

namespace relational {

class RuleLearner2 {
public:

  typedef std::map<int, std::map<int, std::vector<int> > > SingleRuleSetCoveringResponsibilitiesMap;
  typedef std::map<std::string, SingleRuleSetCoveringResponsibilitiesMap > RuleCoveringResponsibilitiesMap;

  RuleLearner2(const StateTransitionL* _p_experiences);
  RuleLearner2(const StateTransitionL* _p_experiences, const arr& experience_weights);
  RuleLearner2(); // use init() later when using this constructor

  /**
   * @brief Initialize the rule set container with a set of experiences
   */
  virtual void init(const StateTransitionL* _p_experiences, const arr& weights, const SymL ignore_symbols=SymL());

  /**
   * @brief Reset and clear all members of the rule set container
   */
  virtual void clear(SymL ignore_symbols=SymL());

  /**
   * @brief Set a rule set (e.g. loaded from file)
   *
   * If you also have the experiences that the rule set was learnt from
   * (or other experiences) first call "init" for this rule set container
   * with the experiences and then set the rule set
   */
  virtual void setRuleSet(const RuleSet& rules, const RuleCoveringResponsibilitiesMap& responsibilities);

  /**
   * @brief Get the stored experiences for an action symbol
   */
  virtual bool getExperiences(Symbol* s, StateTransitionL& experiences, arr& weights);

  virtual bool getExperiences(Symbol* s, StateTransitionL& experiences);

  /**
   * @brief Get the covering responsibilities as a map of maps
   */
  virtual void getRuleCoveringResponsibilities(const RuleSetContainer& rsc, SingleRuleSetCoveringResponsibilitiesMap& responsibilities) const;

  static bool loadRuleCoveringResponsibilities(std::string filename, RuleCoveringResponsibilitiesMap& map);
  static void saveRuleCoveringResponsibilities(std::string log_name, const RuleCoveringResponsibilitiesMap& responsibilities_full);

  /**
   * @brief Write the symbolic transitions in this container to a stream
   */
  virtual void write(ostream& out = std::cout, bool only_action = false, bool additional_experience_info = true) const;

  /**
   * @brief Write the symbolic transitions in this container to a file called "filename"
   */
  virtual void write(const char* filename, bool only_action = false, bool additional_experience_info = true) const;

  /**
   * @brief Write the symbolic transitions in this container to a file to
   * folder log_root
   */
  virtual void writeSymbolicTransitions(std::string log_root, std::map<int, int> experience_index_map) const;

  virtual void writeRuleCoveringResponsibilities(std::string log_root) const;

  virtual void writeExperienceWeights(std::string log_name) const;
  static void writeExperienceWeights(std::string log_name, const std::vector<double>& weights);

  static bool loadExperienceWeights(std::string log_name, arr& experience_weights);

  /**
   * @brief Learn the rule set from scratch with a given set of experiences
   */
  virtual bool learn(StateTransitionL& experiences, arr& weights, std::string log_root, uint repeat=1, SymL ignore_symbols=SymL(), uint debug=2);

  virtual bool learn(StateTransitionL& experiences, std::string log_root, uint repeat=1, SymL ignore_symbols=SymL(), uint debug=2);

  //  virtual bool learn(StateTransitionL& experiences, arr& experience_weights, std::string log_root, uint repeat=1, SymL ignore_symbols=SymL(), uint debug=2);

  /**
   * @brief Learn the rule set incrementally, i.e. by adding just a few new experiences. Only the symbols for which
   *  experiences have been added will be learnt.
   *
   * If you set from_scratch to true, the previous rule set will be discarded and learning will be non-incremental
   *
   * FIXME repeat is not implemented
   */
  virtual bool learnIncremental(StateTransitionL& new_experiences, std::string log_root, uint repeat=1, uint debug=2, bool from_scratch=false);

  /**
   * @brief Get the score of the container for the currently stored experiences
   */
  virtual double score(double cutting_threshold=TL::TL_DOUBLE_MIN, std::ostream& out=std::cout, uint DEBUG=0);
  
  /**
   * @brief Get the score of the container for the currently stored experiences using a particular score function
   */
  virtual double score(relational::learn::ScoreFunction* sf, double cutting_threshold=TL::TL_DOUBLE_MIN, std::ostream& out=std::cout, uint DEBUG=0);

  /**
   * @brief Get the score of the container for the currently stored experiences
   *  (ideally the ones that the container was learnt from), but as a map of
   * individual scores per action symbol
   */
  virtual double score(std::map<Symbol*,double>& scores, double cutting_threshold=TL::TL_DOUBLE_MIN, std::ostream& out=std::cout, uint DEBUG=0,
	relational::learn::ScoreFunction* sf=NULL);

  /**
   * @brief Get the score of the container for the currently stored experiences
   *  (ideally the ones that the container was learnt from), but only for one action symbol
   */
  virtual double score(Symbol *s, double cutting_threshold=TL::TL_DOUBLE_MIN, std::ostream &out=std::cout, uint DEBUG=0);

  /**
   * @brief Get the score of the container for some arbitrary experiences
   *  but only for one action symbol
   */
  virtual double score(Symbol *s, StateTransitionL &exp, arr& weights, double cutting_threshold=TL::TL_DOUBLE_MIN, std::ostream &out=std::cout, uint DEBUG=0);

  /**
   * @brief Check if container is trivial
   */
  virtual bool isTrivial(const RuleSetContainer& rsc, bool ignore_outcomes=true) const;

  /**
   * @brief Check if container is trivial
   */
  virtual bool isTrivial(Symbol* s, bool ignore_outcomes=true) const;

  /**
   * @brief Check if container has experiences in the noise outcome of the default rule
   */
  virtual bool usesNoisyDefaultOutcome(const RuleSetContainer& rsc) const;

  /**
   * @brief Check if container has experiences in the noise outcome of the default rule
   */
  virtual bool usesNoisyDefaultOutcome(Symbol* s) const;

  /**
   * @brief Get the individual rule set container which only contains the
   *  the rules for one action symbol + associated experiences
   *
   *  Attention: there is no copying - you get a reference to the internal RuleSetContainer!
   *
   *  Returns false if no container for the action symbol is found.
   */
  virtual bool getRuleSetContainer(Symbol* symbol, RuleSetContainer& rsc, StateTransitionL& rsc_experiences) const;
  virtual bool getRuleSetContainer(Symbol* symbol, RuleSetContainer& rsc, StateTransitionL& rsc_experiences, uintA& old_ids) const;

  /**
   * @brief Add a new individual rule set container which only contains the
   *  the rules for one action symbol
   */
  virtual void addRuleSetContainer(Symbol* symbol, RuleSetContainer rsc, StateTransitionL rsc_experiences);

  /**
   * @brief Get the rule set container as a class of type RuleSet, where
   *  all actions are merged into it
   */
  virtual void getMergedRuleSet(RuleSet& rules, bool copy_rules=false);

  virtual double coveringProbability(Symbol* symbol, const StateTransition& st, Literal * action);
  static double coveringProbability(RuleSetContainer& rsc, const StateTransition& st, Literal * action);
  static void calc_coveringRules_groundAction(RuleSet& r_grounds, std::vector<uint>& abstract_rules_idx, std::vector<uintA>& abstract2GroundMap,
                                              const RuleSet& all_abstract_rules, const SymbolicState& s, Literal* groundAction);

  static void calcCoverage(relational::StateTransitionL& covered_experiences, uintA& covered_experiences_ids, MT::Array<relational::SubstitutionSet>& subsets, const relational::Rule* r, const relational::StateTransitionL& experiences);

  virtual bool coverage(Symbol* symbol, const StateTransition& st,
                Literal* action,
                std::vector<uint>& covering_rules,
                std::vector<uintA>& covering_outcomes_per_rule,
                bool ignore_default_rule=true);
  virtual void coverage(RuleSetContainer& rsc, const StateTransition& st,
                Literal * action,
                std::vector<uint>& covering_rules,
                std::vector<uintA>& covering_outcomes,
                bool ignore_default_rule=true);


  /**
   * @brief Print out some information for the container, e.g. which symbols
   * and how many rules and experiences it contains
   */
  virtual void info(std::ostream& out=std::cout);

  /**
   * @brief Perform a sanity check of the container
   */
  virtual void sanityCheck(bool ignore_default_rule = false) const;  // EXPENSIVE!!!


  bool approximateScoreAfterAdding(const StateTransition &experience, const double weight, const std::vector<uint> &covering_rules, const std::vector<uintA> &covering_outcomes,
                                     double& probability, double& delta);

  bool approximateScoreChangeAfterRemoving(Symbol* s, uint experienceId, double& delta_minus);

  void setCompressBeforeLearning(bool v) {
    compress_before_learning = v;
  }

  static void disableActionSymbols(SymL& symbols_to_keep);
  static void reenableActionSymbols();

  virtual bool computeRelationalChangeDistribution();

  virtual bool getRelationalChangeProbability(Symbol* symbol, const StateTransition& st,
                                              double & probability, std::vector<int>& effect_experience_ids) const;

//  virtual bool getRuleSetContainers(std::map <Symbol*, RuleSetContainer>& c) const {
//    c = containers;
//    return true;
//  }

  virtual bool getSymbols(std::vector<Symbol*>& s) const {
    for (std::map <Symbol*, arr>::const_iterator it = experience_weights.begin();
         it != experience_weights.end(); it++) {
      s.push_back(it->first);
    }
    return true;
  }

protected:

  /**
   * @brief Helper function to save the symbolic transitions of a single action
   */
  virtual void saveSymbolicTransitions(
    const std::string& file_path,
    const relational::StateTransitionL& res_st,
    const std::map<int, int> experience_index_map) const;


  double approximateNewRulePenalty(const RuleSetContainer &rsc, const StateTransition &experience, double additionalLiterals =0., double alpha =-1);

  void getExperiencesCoveredByRule(StateTransitionL &exp_rule, arr& weights_rule, const RuleSetContainer &rsc_old, const uint rule);

  void createRuleSetContainerFromSingleRule(RuleSetContainer &rsc, StateTransitionL &experiences_rule, arr &weights_rule, const RuleSetContainer &rsc_old, const uint rule);
  void createRuleSetContainerFromSingleRule(RuleSetContainer &rsc, StateTransitionL &experiences_rule, arr &weights_rule, Rule *rule, const MT::Array< uintA > &experiences_per_ruleOutcome);
  void createRuleSetContainerFromSingleRule(RuleSetContainer &rsc, const RuleSetContainer &rsc_old, const uint rule);

  void calculateUpdatedRuleOutcomeScoreDifference(RuleSetContainer &rsc, const StateTransition &experience, double weight,
                                                  const uint &covering_rule,
                                                  double& probability, double& delta, uint DEBUG=0);

  void approximateScoreAfterAddingOneRuleCovers(RuleSetContainer &rsc, const StateTransition &experience, double weight, const uint &covering_rule, const uintA &covering_outcomes,
                                                  double& probability, double& delta);
  void approximateScoreAfterAddingSeveralRulesCover(RuleSetContainer &rsc, const StateTransition &experience, double weight, const std::vector<uint> &covering_rules, const std::vector<uintA> &covering_outcomes,
                                                      double& probability, double& delta);
  void approximateScoreAfterAddingNoRuleCovers(RuleSetContainer &rsc, StateTransitionL& experiences, arr& weights, const StateTransition &experience, double weight, double& probability, double& delta);

  void approximateScoreAfterAdding(RuleSetContainer &rsc, StateTransitionL& experiences, arr& experiences_weights, const StateTransition &experience, double weight,
                                   const std::vector<uint> &covering_rules, const std::vector<uintA> &covering_outcomes,
                                   double& probability, double& delta);

  static Rule *findDistinctiveRuleAddLiterals(const Rule *const rule, const StateTransitionL &experiences, const StateTransition &distinctive_experience);

  /**
   * @brief This version of learn_outcomes deals with the "double covering" bug
   *  This bug happens if experiences with changes A and experiences without changes B
   *  have equivalent post conditions. The method removes experiences in B from the
   *  covering list of the outcomes which covers A
   */
  static void learnUniquelyCoveringOutcomes(Rule* r, MT::Array< uintA >& coveredExperiences_per_outcome,
                                            const StateTransitionL& coveredExperiences, const arr& coveredExperience_weights,
                                            const uintA& covered_experiences_ids, uint DEBUG=0);

  /**
   * @brief Sets the rule set container along with the symbolic transitions
  // @deprecated
   */
  virtual void setRuleSet(const Symbol* s, RuleSetContainer& rsc, StateTransitionL symbol_transitions);

  virtual void setRuleSet(const Symbol* s, RuleSetContainer& rsc,
                          const SingleRuleSetCoveringResponsibilitiesMap& rule_covering_responsibilities);

  virtual void compressExperiences(
      const StateTransitionL& experiences, const arr& original_weights,
      StateTransitionL &experiences_filtered, arr& experience_weights_filtered,
      uintA& expid_to_equivalence_class,
      MT::Array<uintA>& equivalence_class_to_expid
      );


  virtual void decompressRuleSetContainer(const uintA& expid_to_equivalence_class,
      const MT::Array<uintA>& equivalence_class_to_expid,
      const StateTransitionL& experiences, const arr& original_weights,
      RuleSetContainer& rsc);

  virtual bool computeRelationalChangeDistribution(Symbol* symbol);

  std::map <Symbol*, RuleSetContainer> containers;
  // action symbol -> transition set with this action
  std::map <Symbol*, arr> experience_weights;
  // action symbol -> transition set with this action
  std::map <Symbol*, StateTransitionL > experience_map;
  // original experience id -> action symbol -> new experience id
  std::map <int, std::pair<Symbol*, int > > experience_index_map;
  // action symbol -> new experience id -> original experience id
  std::map <Symbol*, std::map<int, int> > reverse_experience_index_map;

  // action symbol -> [effect, effect probability]
  std::map <Symbol*, std::vector< std::pair<SymbolicState, double> > > relational_change_probabilities;
  // action symbol -> [effect, [(new) experience ids] ]
  std::map <Symbol*, std::vector< std::pair<SymbolicState, std::vector<int> > > > relational_change_experience_map;

  bool compress_before_learning;

  bool crossvalidate_accross_actions;

  static MT::Array< relational::Symbol* > mem__all_symbols_reenable;

};

}

#endif // RULELEARNER2_H_
