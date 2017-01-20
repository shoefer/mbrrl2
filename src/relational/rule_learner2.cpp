// RuleLearner2 requires boost!
#ifdef RULELEARNER2

#include "rule_learner2.h"
#include "logging_macros.h"
#include "utilTL.h"

#include <boost/filesystem.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/unordered/unordered_set.hpp>

#include <iostream>
#include <string>
#include <sstream>

#ifdef ROS_LOGGING
#include <log4cxx/logger.h>
#endif

namespace relational {

std::string MTString2String(MT::String str) {
  std::stringstream s;
  s << str;
  return s.str();
}

using namespace std;

bool compareUint(const uint& a, const uint& b) {
  return a > b ? 1 : (a < b ? - 1 : 0);
}

MT::Array< relational::Symbol* > RuleLearner2::mem__all_symbols_reenable;

RuleLearner2::RuleLearner2() : compress_before_learning(false)
{
  compress_before_learning = false;
  crossvalidate_accross_actions = false;
  arr weights;
  init(NULL, weights);
}

RuleLearner2::RuleLearner2(const StateTransitionL *_p_experiences)
{
  compress_before_learning = false;
  crossvalidate_accross_actions = false;
  arr weights(_p_experiences->N);
  weights.setUni(1.);
  init (_p_experiences, weights);
}

RuleLearner2::RuleLearner2(const StateTransitionL *_p_experiences, const arr& weights)
{
  compress_before_learning = false;
  crossvalidate_accross_actions = false;
  init (_p_experiences, weights);
}

void RuleLearner2::init(const StateTransitionL *_p_experiences, const arr& weights, const SymL ignore_symbols)
{
  //containers.clear();
  //experience_map.clear();
  //experience_index_map.clear();
  //reverse_experience_index_map.clear();
  clear(ignore_symbols);

  if (_p_experiences == NULL || _p_experiences->N == 0) {
    RLOG_WARN_STR_NAMED("RuleLearner2.init", "Experiences are NULL or empty");
    return;
  }

  if(weights.N != _p_experiences->N) {
    RLOG_ERROR_STR_NAMED("RuleLearner2.init", "#weights does not match #experiences");
    throw 0;
  }

  uint i;
  // store experiences
  FOR1D( (*_p_experiences), i) {
    StateTransition* t = (*_p_experiences).elem(i);
    Symbol* action_sym = t->action->s;

    if (ignore_symbols.contains(action_sym)) {
      RLOG_DEBUG_STR_NAMED("RuleLearner2.init", "Ignoring experience " << i << " with symbol " << action_sym->name);
//      assert (containers.find(action_sym) != containers.end());
      continue;
    }

    if (containers.find(action_sym) == containers.end()) {
      //containers[action_sym] = StateTransitionL();
      containers[action_sym];
      experience_weights[action_sym];
      experience_map[action_sym];
      relational_change_probabilities[action_sym];
      experience_index_map[i];
      RLOG_DEBUG_STR_NAMED("RuleLearner2.init",
                             "Inserting symbol " << *action_sym);
    }
    StateTransitionL &cur = experience_map[action_sym];
    experience_index_map[i] = make_pair(action_sym, cur.N);
    reverse_experience_index_map[action_sym][cur.N] = i;
    cur.append(t);
    experience_weights[action_sym].append(weights(i));

    RLOG_DEBUG_STR_NAMED("RuleLearner2.init",
                           "Experience " << i << " [w=" << weights(i) << "] becomes " << (cur.N-1) << " of " << *action_sym);
  }

  // init containers
  std::map <Symbol*, StateTransitionL >::iterator it;
  std::map <Symbol*, arr >::iterator itw=experience_weights.begin();
  for (it = experience_map.begin(); it != experience_map.end(); it++, itw++) {
    Symbol* action_sym = it->first;

    if (ignore_symbols.contains(action_sym)) {
      RLOG_INFO_STR_NAMED("RuleLearner2.init", "Ignoring initialization of container " << action_sym->name);
      continue;
    }

    RuleSetContainer& rsc = containers[action_sym];
    rsc.init(&it->second, itw->second);

    RLOG_INFO_STR_NAMED("RuleLearner2.init", action_sym->name
                          << " contains " << it->second.N << " experiences");
  }



  crossvalidate_accross_actions = false; // FIXME

  if (crossvalidate_accross_actions) {
    std::map <Symbol*, RuleSetContainer >::iterator it;
    for (it = containers.begin(); it != containers.end(); it++) {
      StateTransitionL* cv_transitions = new StateTransitionL;

      std::map <Symbol*, StateTransitionL >::iterator ite;
      for(ite=experience_map.begin(); ite !=experience_map.end(); ite++) {
        if (ite->first == it->first) continue;
        uint i;
        FOR1D(ite->second, i) {
          StateTransition* e = ite->second(i);
          StateTransition* e_new = new StateTransition;
          *e_new = *e;

          // THIS IS NOT GOOD! rather create a new literal!
          e_new->action->s = it->first;
          cv_transitions->append(e_new);
        }
      }
      RuleSetContainer& rsc = it->second;
      delete rsc.p_experiences_crossvalidation;
      rsc.p_experiences_crossvalidation = cv_transitions;
    }

  }
}

void RuleLearner2::setRuleSet(const Symbol *s, RuleSetContainer& rsc,
                                                 const SingleRuleSetCoveringResponsibilitiesMap& rule_covering_responsibilities) {
  if (rsc.rules.num() == 0) {
    RLOG_ERROR_STR_NAMED("RuleLearner2.setRuleSet", *s << " RuleSetContainer does not contain any rules!");
    return;
  }

  rsc.experiences_per_rule.resize(rsc.rules.num());
  rsc.experiences_per_ruleOutcome.resize(rsc.rules.num());
  //  rsc.nonDefaultRules_per_experience is resized in RuleSetContainer::init

  // Load all experiences from

  // special care for default rule
  rsc.experiences_per_ruleOutcome.elem(0).resize(2);

  SingleRuleSetCoveringResponsibilitiesMap::const_iterator it;
  for (it = rule_covering_responsibilities.begin(); it != rule_covering_responsibilities.end(); it++) {
    uint r = it->first;
    const std::map<int, std::vector<int> > & outcomes = it->second;
    if (rsc.experiences_per_ruleOutcome.elem(r).N == 0)
      rsc.experiences_per_ruleOutcome.elem(r).resize(outcomes.size());

    for (std::map<int, std::vector<int> >::const_iterator ito = outcomes.begin(); ito != outcomes.end(); ito++) {
      uint o = ito->first;
      for (std::vector<int>::const_iterator ite = ito->second.begin(); ite !=  ito->second.end(); ite++) {
         rsc.experiences_per_rule.elem(r).append(*ite);
         rsc.experiences_per_ruleOutcome.elem(r).elem(o).append(*ite);
      }
    }
  }

  //rsc.nonDefaultRules_per_experience.resize(rsc.p_experiences->N);

  std::stringstream log;
  rsc.write(log);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.setRuleSet", log.str());
  log.flush();

}

void RuleLearner2::setRuleSet(const Symbol *s, RuleSetContainer& rsc, StateTransitionL symbol_experiences) {

  if (rsc.rules.num() == 0) {
    RLOG_ERROR_STR_NAMED("RuleLearner2.setRuleSet", *s << " RuleSetContainer does not contain any rules!");
    return;
  }

  rsc.experiences_per_rule.resize(rsc.rules.num());
  rsc.experiences_per_ruleOutcome.resize(rsc.rules.num());
  //  rsc.nonDefaultRules_per_experience is resized in RuleSetContainer::init

  // the default rule experiences
  uintA default_rule_experiences_ids;

  // all experiences cached; will be reduced step by step
  uintA all_experiences_ids;
  for (uint i = 0; i < symbol_experiences.N; i++) {
    all_experiences_ids.append(i);
  }

  uint i;
  FOR1D_(rsc.rules, i) {
    Rule* r = rsc.rules.elem(i);

//    cout << "----------" << endl;
//    PRINT (*r);

    // resize
    rsc.experiences_per_ruleOutcome(i).resize(r->outcomes.N);

    // default rule
    if (i == 0) {
      continue;
    }

    // which experiences are covered by rule
    StateTransitionL covered_experiences;
    uintA covered_experiences_ids;
    MT::Array<SubstitutionSet> subsets;
    calcCoverage(covered_experiences, covered_experiences_ids, subsets, r, symbol_experiences);

//    PRINT(covered_experiences_ids);

    // remove these from all_experiences / add to default rule if non uniquely covered
    uint k;
    FOR1D(covered_experiences_ids, k) {
      // experience -> rules
      //rsc.nonDefaultRules_per_experience.insertInSorted(i);
      rsc.nonDefaultRules_per_experience(covered_experiences_ids(k)).append(i);

      RLOG_DEBUG_STR_NAMED("RuleLearner2.setRuleSet",
                             "rule " << *r
                             << " covers "
                             << covered_experiences_ids(k) << endl
                             );

      // manage all_experience_ids / default_rule_experiences_ids
      if (all_experiences_ids.findValue(covered_experiences_ids(k)) >= 0) {
        all_experiences_ids.removeAllValues(covered_experiences_ids(k));
      } else {
        // not in all_experiences -> was removed before -> non uniquely covered
        default_rule_experiences_ids.append(covered_experiences_ids(k));
      }

      // rule -> experiences
      rsc.experiences_per_rule(i).append(covered_experiences_ids(k));

      // rule outcomes -> experiences
      // TODO does this really work?
      bool experience_covered=false;

      uint j;
      FOR1D(r->outcomes, j) {

        if (j == r->outcomes.N-1) // noise outcome
          break;

        SubstitutionSet& sub = subsets(k);
        uint l;
        FOR1D_(sub, l){
          //PRINT(*sub.elem(l));
          if (reason::covers(covered_experiences(k)->post.lits, r->outcomes(j), true, sub.elem(l)) ) {
            //rsc.experiences_per_ruleOutcome(i)(j).insertInSorted(covered_experiences_ids(k));
//            PRINT(covered_experiences(k)->post.lits);
            rsc.experiences_per_ruleOutcome(i)(j).append(covered_experiences_ids(k));
            experience_covered = true;
            break;
          }
        }
        // we need this hack because reason::covers always covers empty outcome
        // since the empty outcome is always the last before the noise outcome
        // this should not be harmful
        if (experience_covered)
          break;
      }
      if (!experience_covered) { // is noise outcome?
        rsc.experiences_per_ruleOutcome(i)(r->outcomes.N-1).append(covered_experiences_ids(k));
      }
    }

  }

  // default rules
  FOR1D(all_experiences_ids, i) {
    uint bin;
    if (symbol_experiences(all_experiences_ids(i))->changes.N == 0) {
      bin = 0;
    } else {
      bin = 1;
      RLOG_WARN_STR_NAMED("RuleLearner2.setRuleSet",
                             "Symbol " << *s << ", experience " << all_experiences_ids(i) << " explained as noise; not yet implemented");
    }
    rsc.experiences_per_ruleOutcome(0)(bin).insertInSorted(all_experiences_ids(i), compareUint);
//    rsc.experiences_per_rule(0).insertInSorted(all_experiences_ids(i));
//    rsc.experiences_per_ruleOutcome(0)(bin).append(all_experiences_ids(i));
    rsc.experiences_per_rule(0).append(all_experiences_ids(i));
  }

  // TODO finally remove all non-uniquely covered experiences from rule outcome covering

  //PRINT(rsc.experiences_per_rule);
//  PRINT(rsc.experiences_per_ruleOutcome);
}

void RuleLearner2::setRuleSet(const RuleSet& rules, const RuleCoveringResponsibilitiesMap& responsibilities)
{
  std::map <Symbol*, RuleSetContainer>::iterator itr = containers.begin();

  for (; itr != containers.end(); itr++) {
    Symbol* s = itr->first;

    // it is safer to put this here in case less containers than clusters have been initialized
    std::map <Symbol*, StateTransitionL >::iterator ite = experience_map.find(s);
    assert (ite != experience_map.end());

    RuleSetContainer& rsc = itr->second;
    StateTransitionL st = ite->second;

    // remove all rules
    rsc.rules.clear();
    // add default rule
    rsc.rules.append(Rule::generateDefaultRule());

    // filter rules for this symbol
    uint i;
    FOR1D_(rules, i) {
      Rule* r = rules.elem(i);

      if (Rule::isDefaultRule(r)) {
        // skip default rule
        continue;
      }

      if (*(r->action->s) == *s) {
        Rule* r_copy = new Rule;
        r_copy->copyBody(*r);

        rsc.rules.append(r_copy);
      }
    }
    if (rsc.rules.num() == 0)
      // ignore; no rule set for this action type
      continue;

    if (responsibilities.empty()) {
      setRuleSet(s, rsc, st);
    } else {
      std::string sname = MTString2String(s->name);
      RuleCoveringResponsibilitiesMap::const_iterator rit = responsibilities.find(sname);
      if (rit == responsibilities.end()) {
        if (sname[0] != 'c') {
          RLOG_WARN_STR_NAMED("RuleLearner2.setRuleSet", "No responsibilities for symbol " << sname << " - ignoring");
          continue;
        } else {
          RLOG_ERROR_STR_NAMED("RuleLearner2.setRuleSet", "No responsibilities for symbol " << sname);
          throw 0;
        }
      }
      setRuleSet(s, rsc, rit->second);
    }
  }
}

void RuleLearner2::clear(SymL ignore_symbols)
{
  SymL cleanup_symbols;

  std::map <Symbol*, RuleSetContainer >::iterator it;
  for (it = containers.begin(); it != containers.end(); it++) {
      Symbol* s = it->first;
    // it is safer to put this here in case less containers than clusters have been initialized
    std::map <Symbol*, arr >::iterator itw = experience_weights.find(s);
    assert (itw != experience_weights.end());

    if (ignore_symbols.contains(it->first)) {
      RLOG_INFO_STR_NAMED("RuleLearner2.clear", "Ignoring symbol " << it->first->name);
      continue;
    }

    cleanup_symbols.append(it->first);
    (it->second).clear();
    (itw->second).clear();
  }

  if (ignore_symbols.N == 0) {
    // clean up everything
    containers.clear();
    experience_map.clear();
    experience_weights.clear();
    experience_index_map.clear();
    reverse_experience_index_map.clear();
    relational_change_probabilities.clear();
    return;
  }

  uint s;
  FOR1D(cleanup_symbols, s) {
    RLOG_INFO_STR_NAMED("RuleLearner2.clear", "Cleaning symbol " << cleanup_symbols(s)->name);

    std::map <Symbol*, RuleSetContainer >::iterator itc = containers.find(cleanup_symbols(s));
    assert(itc != containers.end());
    containers.erase(itc);

    std::map <Symbol*, StateTransitionL >::iterator ite = experience_map.find(cleanup_symbols(s));
    assert(ite != experience_map.end());
    experience_map.erase(ite);

    std::map <Symbol*, arr >::iterator itw = experience_weights.find(cleanup_symbols(s));
    assert(itw != experience_weights.end());
    experience_weights.erase(itw);

    // index map should be overwritten
    //std::map <Symbol*, RuleSetContainer >::iterator itim = experience_index_map.find(cleanup_symbols(s));

    std::map <Symbol*, std::map<int, int> >::iterator iteim = reverse_experience_index_map.find(cleanup_symbols(s));
    assert(iteim != reverse_experience_index_map.end());
    reverse_experience_index_map.erase(iteim);

    std::map <Symbol*, std::vector< std::pair<SymbolicState, double> > >::iterator itrc = relational_change_probabilities.find(cleanup_symbols(s));
    assert(itrc != relational_change_probabilities.end());
    relational_change_probabilities.erase(itrc);
  }

}

bool RuleLearner2::getExperiences(Symbol* s, StateTransitionL& experiences) {
  arr weights;
  return getExperiences(s, experiences, weights);
}

bool RuleLearner2::getExperiences(Symbol* s, StateTransitionL& experiences, arr& weights) {
  std::map <Symbol*, StateTransitionL >::iterator its = experience_map.find(s);
  std::map <Symbol*, arr >::iterator itw = experience_weights.find(s);
  if (its == experience_map.end()) {
    return false;
  }
  assert (itw != experience_weights.end());
  experiences = its->second;
  weights = itw->second;
  return true;
}

void RuleLearner2::write(ostream &out, bool only_action, bool additional_experience_info) const
{
  std::map <Symbol*, RuleSetContainer >::const_iterator it;
  for (it = containers.begin(); it != containers.end(); it++) {
    /*
    std::map <Symbol*, StateTransitionL >::const_iterator emit = experience_map.find(it->first);
    if (emit->second.N == 0) {
      RLOG_DEBUG_STR_NAMED("RuleLearner2.write",
                          "RuleSet " << it->first->name << " is empty!");
      continue;
    }*/
    out << "# ============================ " << endl;
    out << it->first->name << endl;
    (it->second).write(out, only_action, additional_experience_info);
  }
}

void RuleLearner2::write(const char *filename, bool only_action, bool additional_experience_info) const
{
  ofstream out(filename);
  write(out, only_action, additional_experience_info);
}

void RuleLearner2::writeSymbolicTransitions(std::string log_name, std::map<int, int> old_index_map) const
{
//  std::map <Symbol*, RuleSetContainer> containers;
//  std::map <Symbol*, StateTransitionL > experience_map;
//  std::map <int, std::pair<Symbol*, int > > experience_index_map;
//  std::map <Symbol*, std::map<int, int> > reverse_experience_index_map;

  map <Symbol*, std::map<int, int> >::const_iterator it = reverse_experience_index_map.begin();
  for (; it != reverse_experience_index_map.end(); it++) {

    map<int, int> experience_index_remap;
    map<int, int>::const_iterator mit = it->second.begin();
    for (; mit != it->second.end(); mit++) {
      // local index in transitions -> index in old_index_map
      experience_index_remap[mit->first] = old_index_map[mit->second];
    }

    map <Symbol*, StateTransitionL >::const_iterator tit = experience_map.find(it->first);
//    const StateTransitionL transitions = experience_map.at(s);

    stringstream l_name;
    const Symbol* s = it->first;
    l_name << log_name << s->name << ".log";
    
    saveSymbolicTransitions(l_name.str(),
                            tit->second,
                            experience_index_remap);
  }

}

void RuleLearner2::saveSymbolicTransitions(
    const std::string& file_path,
    const relational::StateTransitionL& res_st,
    const std::map<int, int> experience_index_map) const {
  ofstream st_out;
  st_out.open( file_path.c_str() );
  uint i;
  FOR1D(res_st, i) {
    st_out << "#" << i << endl;
    st_out << res_st(i)->pre << endl;
    st_out << * (res_st(i)->action);

    map<int, int>::const_iterator new_index = experience_index_map.find(i);
    if (new_index != experience_index_map.end()) {
      st_out << " # was " << new_index->second;
    }
    st_out << endl;
    st_out << res_st(i)->post << endl;
    st_out << res_st(i)->reward << endl << endl;
  }
  st_out.close();
}


void RuleLearner2::getRuleCoveringResponsibilities(const RuleSetContainer& rsc, SingleRuleSetCoveringResponsibilitiesMap& responsibilities) const {
    uint ruleId, outcomeId, experienceId;
    FOR1D(rsc.experiences_per_ruleOutcome, ruleId) {
      const MT::Array < uintA >& eprO = rsc.experiences_per_ruleOutcome.elem(ruleId);
      responsibilities[ruleId]; // create map for rule (even if empty)
      FOR1D(eprO, outcomeId) {
        const uintA & epO = eprO.elem(outcomeId);
        responsibilities[ruleId][outcomeId].resize(epO.N);  // create map for outcome (even if empty)
        FOR1D(epO, experienceId) {
          responsibilities[ruleId][outcomeId][experienceId] = epO.elem(experienceId);
        }
      }
    }
}

bool RuleLearner2::loadRuleCoveringResponsibilities(std::string filename, RuleCoveringResponsibilitiesMap& responsibilities) {
  responsibilities.clear();
  using namespace std;

  if (!boost::filesystem::exists(filename)) {
    RLOG_WARN_STR_NAMED("RuleLearner2.loadRuleCoveringResponsibilities", " Responsibility file does not exist: " << filename);
    return false;
  }

  ifstream in;
  in.open(filename.c_str(), ios::in);

  boost::archive::text_iarchive iarch(in);
  iarch >> responsibilities;

  std::stringstream sslog;
  RuleCoveringResponsibilitiesMap::iterator it = responsibilities.begin();
  for (; it != responsibilities.end(); it++) {
    sslog << " " << it->first << ":" << endl;
    SingleRuleSetCoveringResponsibilitiesMap::iterator itk = it->second.begin();
    for (; itk != it->second.end(); itk++) {
      sslog << "\tR" << itk->first << ": " << endl;
      std::map<int, std::vector<int> >::iterator its = itk->second.begin();
      for (; its != itk->second.end(); its++) {
        sslog << "\t\to" << its->first << ":  [";
        std::vector<int>::iterator itw = its->second.begin();
        for (; itw != its->second.end(); itw++) {
          sslog << " " << *itw;
        }
        sslog << "]" << endl;
      }
      sslog << endl;
    }
  }
  RLOG_DEBUG_STR_NAMED("RuleLearner2.loadRuleCoveringResponsibilities", " Responsibilities: " << endl << sslog.str());


  in.close();

  return true;
}

void RuleLearner2::saveRuleCoveringResponsibilities(std::string log_name, const RuleCoveringResponsibilitiesMap& responsibilities_full) {
  using namespace std;

  ofstream out;
  out.open(log_name.c_str(), ofstream::out);
  boost::archive::text_oarchive oarch(out);
  oarch << responsibilities_full;
  out.close();
}

void RuleLearner2::writeRuleCoveringResponsibilities(std::string log_name) const {
  map <Symbol*, RuleSetContainer>::const_iterator it;

  RuleCoveringResponsibilitiesMap responsibilities_full;

  for (it = containers.begin(); it != containers.end(); it++) {
    Symbol* s = it->first;
    const RuleSetContainer& rsc = it->second;

    SingleRuleSetCoveringResponsibilitiesMap& responsibilities = responsibilities_full[MTString2String(s->name)];
    getRuleCoveringResponsibilities(rsc, responsibilities);
  }

  saveRuleCoveringResponsibilities(log_name, responsibilities_full);

}

void RuleLearner2::writeExperienceWeights(std::string log_name) const {
  using namespace std;

  vector<double> weights;

  int N = experience_index_map.size();
  for (int i=0; i < N; i++) {
    const pair<Symbol*, int >& p = experience_index_map.at(i);
    double w = experience_weights.at(p.first).elem(p.second);
    weights.push_back(w);;
  }

  writeExperienceWeights(log_name, weights);
}

void RuleLearner2::writeExperienceWeights(std::string log_name, const vector<double>& weights)  {
  using namespace std;

  ofstream out;
  out.open(log_name.c_str(), ofstream::out);
  boost::archive::text_oarchive oarch(out);
  oarch << weights;
  out.close();
}


bool RuleLearner2::loadExperienceWeights(std::string filename, arr& experience_weights) {
  experience_weights.clear();
  using namespace std;

  if (!boost::filesystem::exists(filename)) {
    RLOG_WARN_STR_NAMED("RuleLearner2.loadExperienceWeights", " Weights file does not exist: " << filename);
    return false;
  }

  vector<double> weights;

  ifstream in;
  in.open(filename.c_str(), ios::in);

  boost::archive::text_iarchive iarch(in);
  iarch >> weights;

  in.close();

  for (vector<double>::iterator it = weights.begin(); it != weights.end(); it++) {
    experience_weights.append(*it);
  }

  return true;
}

/*
void RuleLearner2::write_experiencesWithRules(ostream &os) const
{
  std::map <Symbol*, RuleSetContainer >::const_iterator it;
  for (it = containers.begin(); it != containers.end(); it++) {
    (it->second).write_experiencesWithRules(os);
  }
}

void RuleLearner2::write_rulesWithExperiences(ostream &os) const
{
  std::map <Symbol*, RuleSetContainer >::const_iterator it;
  for (it = containers.begin(); it != containers.end(); it++) {
    (it->second).write_rulesWithExperiences(os);
  }
}*/

void RuleLearner2::sanityCheck(bool ignore_default_rule) const
{
  std::map <Symbol*, RuleSetContainer >::const_iterator it;
  for (it = containers.begin(); it != containers.end(); it++) {
    (it->second).sanityCheck(ignore_default_rule);
  }
}


bool RuleLearner2::getRuleSetContainer(Symbol* s, RuleSetContainer& rsc, StateTransitionL& rsc_experiences) const {
  uintA old_ids;
  return getRuleSetContainer(s, rsc, rsc_experiences, old_ids);
}
bool RuleLearner2::getRuleSetContainer(Symbol* s, RuleSetContainer& rsc, StateTransitionL& rsc_experiences, uintA& old_ids) const {
  assert (s != NULL);
  std::map <Symbol*, RuleSetContainer >::const_iterator itc = containers.find(s);
  if (itc == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.getRuleSetContainer",
                       "RuleSetContainer for symbol " << s->name << " not found");
    return false;
  }
  std::map <Symbol*, StateTransitionL>::const_iterator its = experience_map.find(s);
  assert (its != experience_map.end());

  rsc = itc->second;
  rsc_experiences = its->second;

  // append experience indexes
  stringstream old_ids_log;
  old_ids_log << " Relative to rule set idx -> global idx" << endl;
  old_ids.clear();
  std::map  <Symbol*, std::map<int, int> > ::const_iterator itm = reverse_experience_index_map.find(s);
  assert (itm != reverse_experience_index_map.end());
  int i = 0;
  for (map<int,int>::const_iterator it = itm->second.begin(); it != itm->second.end(); it++, i++) {
    old_ids_log << " " << i << ") " << it->first << " -> " << it->second << endl;
    old_ids.insertInSorted(it->second,compareUint);
//    old_ids.append(it->second);
  }
  old_ids.reverse();

  PRINTOUT(old_ids, old_ids_log);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.getRuleSetContainer",
                     "RuleSetContainer for symbol " << s->name << " oldIds: " << endl << old_ids_log.str());


  if (rsc_experiences.N != old_ids.N) {
    RLOG_ERROR_STR_NAMED("RuleLearner2.getRuleSetContainer",
                        "Error collecting global IDs of rule set " << s->name
                        << " missed " << (rsc_experiences.N - old_ids.N) << " indices");
    throw "";
  }

  return true;
}

void RuleLearner2::addRuleSetContainer(Symbol* s, RuleSetContainer rsc, StateTransitionL rsc_experiences) {
  assert (s != NULL);
  std::map <Symbol*, RuleSetContainer >::iterator itc = containers.find(s);
  if (itc != containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.getRuleSetContainer",
                       "RuleSetContainer for symbol " << s->name << " not empty");
  }
  containers[s] = rsc;
  experience_map[s] = rsc_experiences;
}

bool RuleLearner2::learnIncremental(StateTransitionL &new_experiences, std::string log_name, uint, uint debug, bool from_scratch)
{
  // sort new experiences according to symbols
  set<Symbol*> updated_symbols;
  uint i;

  std::map <Symbol*, StateTransitionL > experience_map_new;
  FOR1D(new_experiences, i){
    StateTransition* st = new_experiences(i);
    Symbol* s = st->action->s;
    experience_map_new[s].append(st);
    updated_symbols.insert(s);
  }

  // retrain updated containers
  for (set<Symbol*>::iterator it=updated_symbols.begin(); it != updated_symbols.end(); it++) {
    Symbol* s = *it;
    RLOG_INFO_STR_NAMED("RuleLearner2.learnIncremental", "Training " << s->name);
    stringstream l_name, vl_name;
    l_name << log_name << s->name << ".log";
    vl_name << log_name << s->name << "_FULL.log";
    if (!from_scratch) {
      // learn really incrementally
      learn::learn_rules(containers[s], experience_map_new[s], true, l_name.str().c_str(), vl_name.str().c_str(), debug);
    }
    // append to experience cache
    FOR1D(experience_map_new[s], i){
      experience_map[s].append(experience_map_new[s].elem(i));
    }
    if (from_scratch) {
      // re-learn rules from scratch including the new experience
      learn::learn_rules(containers[s], experience_map[s], false, l_name.str().c_str(), vl_name.str().c_str(), debug);
    }
  }

  return true;
}

//bool RuleLearner2::learn(StateTransitionL &experiences, std::string log_name, uint repeat, SymL ignore_symbols, uint debug) {
//  arr experience_weights(experiences.N);
//  experience_weights.setUni(1);
//  return learn(experiences, experience_weights, log_name, repeat, ignore_symbols, debug);
//}

void RuleLearner2::compressExperiences(
    const StateTransitionL& experiences,
    const arr& original_weights,
    StateTransitionL &experiences_filtered, arr& experience_weights_filtered,
    uintA& expid_to_equivalence_class,
    MT::Array<uintA>& equivalence_class_to_expid
    ) {

  experiences_filtered.clear();
  experience_weights_filtered.clear();

  expid_to_equivalence_class.clear();
  expid_to_equivalence_class.resize(experiences.N);
  equivalence_class_to_expid.clear();

  uint i,j;
  FOR1D(experiences, i) {
    StateTransition *st = experiences(i);
    double w = original_weights(i);
    bool contains=false;
    FOR1D(experiences_filtered, j) {
      StateTransition *stf = experiences_filtered(j);
      if ( *st == *stf ) {
        contains = true;
        //experience_weights_filtered(j)++;
        experience_weights_filtered(j)+=w;
        expid_to_equivalence_class(i) = j;
        equivalence_class_to_expid(j).append(i);
        break;
      }
    }

    if (!contains) {
//      StateTransition* st_new = new StateTransition;
//      *st_new = *st; // copy
//      experiences_filtered.append(st_new);
      experiences_filtered.append(st);
      //experience_weights_filtered.append(1);
      experience_weights_filtered.append(w);
      expid_to_equivalence_class(i) = experiences_filtered.N-1;
      equivalence_class_to_expid.append(uintA());
      equivalence_class_to_expid(experiences_filtered.N-1).append(i);
    }
  }

  // now normalize weights
//  double min_w = experience_weights_filtered.min();
//  cout << (min_w) <<endl;
//  FOR1D(experience_weights_filtered,i){
//    experience_weights_filtered(i) /= min_w;
//  }

  #ifdef ROS_LOGGING
  log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger("ros.RuleLearner2.compressExperiences");
  if (logger && logger->getLevel() == log4cxx::Level::getDebug()) {

    std::stringstream ss;
    PRINTOUT(equivalence_class_to_expid, ss);

    ss << endl;
    FOR1D(experiences_filtered, i) {
      StateTransition* s = experiences_filtered(i);
      ss << "#" << i << " - subsumes [";
      FOR1D(equivalence_class_to_expid(i),j) {
        uint local_id = equivalence_class_to_expid(i)(j);
        ss << " " << local_id;
        ss << " (g=" << this->reverse_experience_index_map[s->action->s][local_id] << ")";
      }
      ss<< "]" << endl;
      s->pre.write(ss);
      ss << endl;
      ss << *(s->action) << endl;
      s->post.write(ss);
      ss << s->reward;
      ss << endl << endl;
    }

    RLOG_DEBUG_STR_NAMED("RuleLearner2.compressExperiences",
                         std::endl << ss.str());
  }
  #endif
}

void RuleLearner2::decompressRuleSetContainer(
    const uintA& expid_to_equivalence_class,
    const MT::Array<uintA>& equivalence_class_to_expid,
    const StateTransitionL& experiences,
    const arr& original_weights,
    RuleSetContainer& rsc
    ) {

  //
  MT::Array< uintA > nDRpe = rsc.nonDefaultRules_per_experience;
  MT::Array< uintA > epr = rsc.experiences_per_rule;
  MT::Array< MT::Array < uintA > > eprO = rsc.experiences_per_ruleOutcome;

  rsc.nonDefaultRules_per_experience.clear();
  rsc.nonDefaultRules_per_experience.resize(experiences.N);
  uint i,j;

  FOR1D(experiences, i) {
    uint eqc = expid_to_equivalence_class(i);
    rsc.nonDefaultRules_per_experience(i).append(nDRpe(eqc));
  }

  // fix experiences_per_rule
  FOR1D(epr, i) {
    rsc.experiences_per_rule(i).clear();
    FOR1D(epr(i), j) {
      uintA & ids = equivalence_class_to_expid (epr(i)(j));
      uint id;
      FOR1D(ids, id) {
        rsc.experiences_per_rule(i).append(ids(id));
      }
    }

    // fix experiences_per_ruleOutcome
    FOR1D(eprO(i), j) {
      rsc.experiences_per_ruleOutcome(i)(j).clear();
      uint k;
      FOR1D(eprO(i)(j), k) {
        uintA& ids = equivalence_class_to_expid (eprO(i)(j)(k));
        uint id;
        FOR1D(ids, id) {
          rsc.experiences_per_ruleOutcome(i)(j).append(ids(id));
        }
      }
    }
  }

//  PRINT(eprO);
//  PRINT(rsc.experiences_per_ruleOutcome);

  // overwrite internal experiences
  StateTransitionL* new_p_experiences = new StateTransitionL;
  *new_p_experiences = experiences;
  //delete rscc1.p_experiences;
  rsc.p_experiences = new_p_experiences;

  rsc.experience_weights.clear();
//  rsc.experience_weights.setUni(1., experiences.N);
  rsc.experience_weights = original_weights;

}

void RuleLearner2::disableActionSymbols(relational::SymL& symbols_to_keep) {

  // begin remove symbols
  // remove unnecessary symbols to speed up learning
  Symbol* defaultAction = Symbol::getDefaultAction(); // get it before erasing mem__all_symbols
  mem__all_symbols_reenable.clear();
  mem__all_symbols_reenable = Symbol::mem__all_symbols;

//      PRINTOUT(mem__all_symbols_old, cout);

  Symbol::mem__all_symbols.clear();
  uint i;
  FOR1D(mem__all_symbols_reenable, i) {
    Symbol* sym = mem__all_symbols_reenable(i);
    // remove all but current action symbols
    if (sym->symbol_type == Symbol::action && sym != defaultAction && symbols_to_keep.findValue(sym) < 0) {
      RLOG_DEBUG_STR_NAMED("RuleLearner2.disableActionSymbols", " Temporally removing " << sym->name);
      continue;
    }
    // keep symbol
    Symbol::mem__all_symbols.append(sym);
  }

  // <logging>
  std::stringstream ss;
  ss << endl;
  PRINTOUT(mem__all_symbols_reenable, ss);
  ss << endl;
  ss << " Removed " << (mem__all_symbols_reenable.N - Symbol::mem__all_symbols.N) << " symbols:" << endl;
  PRINTOUT(Symbol::mem__all_symbols, ss);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.disableActionSymbols", ss.str() << endl);
  // </logging>
  // end remove symbols
}

void RuleLearner2::reenableActionSymbols() {
  if (mem__all_symbols_reenable.N == 0) {
    RLOG_WARN_STR_NAMED("RuleLearner2.reenableActionSymbols", "Did you disable symbols before?");
    return;
  }
  Symbol::mem__all_symbols.clear();
  Symbol::mem__all_symbols = mem__all_symbols_reenable;

}

bool RuleLearner2::learn(StateTransitionL &experiences, std::string log_name, uint repeat, SymL ignore_symbols, uint debug) {
  arr weights(experiences.N);
  weights.setUni(1.);
  return learn(experiences, weights, log_name, repeat, ignore_symbols, debug);
}

bool RuleLearner2::learn(StateTransitionL &experiences, arr& weights, std::string log_name, uint repeat, SymL ignore_symbols, uint debug)
{
  this->clear(ignore_symbols);
  this->init(&experiences, weights, ignore_symbols);

  // gather all constants from experiences only only use the present ones
  // to speed up learning
  uintA old_constants, old_constants_;
  old_constants_ = relational::reason::getConstants(); // returns a reference...
  uint i_;
  FOR1D(old_constants_, i_) {
	old_constants.append(old_constants_(i_));
  }
  
  // cache old constants to not interfere with surrounding code
  uintA constants;
  relational::getConstants(constants, experiences);
  constants.sort(compareUint);
  constants.reverse();
  reason::setConstants(constants);

//  PRINT(constants);

  if (repeat == 0) {
    RLOG_ERROR_STR_NAMED("RuleLearner2.learn", "repeat must not be 0! Setting to 1");
    repeat = 1;
  }

  std::map <Symbol*, RuleSetContainer >::iterator itc = containers.begin();
  for (; itc != containers.end(); itc++) {
    Symbol* s = itc->first;

    // it is safer to put this here in case less containers than clusters have been initialized
    std::map <Symbol*, StateTransitionL >::iterator its = experience_map.find(s);
    assert (its != experience_map.end());
    std::map <Symbol*, arr >::iterator itw = experience_weights.find(s);
    assert (itw != experience_weights.end());

    if (ignore_symbols.contains(s)) {
      RLOG_INFO_STR_NAMED("RuleLearner2.learn", "NOT training " << s->name << " because it is on ignore list");
    } else {
      RLOG_INFO_STR_NAMED("RuleLearner2.learn", "Training " << s->name);
      StateTransitionL* exp = & its->second;
      arr* exp_weights = & itw->second;

      // begin remove symbols
      // remove unnecessary symbols to speed up learning
      Symbol* defaultAction = Symbol::getDefaultAction(); // get it before erasing mem__all_symbols
      MT::Array< Symbol* > mem__all_symbols_old = Symbol::mem__all_symbols;

//      PRINTOUT(mem__all_symbols_old, cout);

      SymL state_symbols;
      relational::getStateSymbols(state_symbols, *exp);

      Symbol::mem__all_symbols.clear();
      uint i;
      FOR1D(mem__all_symbols_old, i) {
        Symbol* sym = mem__all_symbols_old(i);
        // remove all but current action symbols
        if (sym->symbol_type == Symbol::action && sym != defaultAction && sym != s) {
          RLOG_DEBUG_STR_NAMED("RuleLearner2.learn", " Temporally removing " << sym->name);
          continue;
        }
        // state symbols
        if (sym->symbol_type != Symbol::action && state_symbols.findValue(sym) < 0) {
          RLOG_DEBUG_STR_NAMED("RuleLearner2.learn", " Temporally removing " << sym->name);
          continue;
        }
        // keep symbol
        Symbol::mem__all_symbols.append(sym);
      }

      // <logging>
      std::stringstream ss;
      ss << endl;
      PRINTOUT(mem__all_symbols_old, ss);
      ss << endl;
      ss << " Removed " << (mem__all_symbols_old.N - Symbol::mem__all_symbols.N) << " symbols:" << endl;
      PRINTOUT(Symbol::mem__all_symbols, ss);
      RLOG_DEBUG_STR_NAMED("RuleLearner2.learn", ss.str() << endl);
      // </logging>
      // end remove symbols

      // compress
      uintA expid_to_equivalence_class;
      MT::Array<uintA> equivalence_class_to_expid;

      if (compress_before_learning) {
        RLOG_INFO_STR_NAMED("RuleLearner2.learn", "  Compressing " << its->second.N << " experiences");

        StateTransitionL* experience_filtered = new StateTransitionL;
        arr* experience_weights_filtered = new arr;
        compressExperiences(*exp, *exp_weights, *experience_filtered, *experience_weights_filtered, expid_to_equivalence_class, equivalence_class_to_expid);

        exp = experience_filtered;
        exp_weights = experience_weights_filtered;

        RLOG_INFO_STR_NAMED("RuleLearner2.learn", "  Compressed " << its->second.N << " down to " << exp->N);
      }

      vector<double> scores;
      vector<RuleSetContainer> rscs(repeat);
      for (uint i = 0; i < repeat; i++) {
        rscs[i].init(exp, *exp_weights);
        if (crossvalidate_accross_actions) {
          rscs[i].p_experiences_crossvalidation = itc->second.p_experiences_crossvalidation;
        }
        stringstream l_name, vl_name;
        l_name << log_name << s->name << "_" << i << " .log";
        vl_name << log_name << s->name << "_FULL_" << i << ".log";
        learn::learn_rules(rscs[i], *exp, *exp_weights, false, l_name.str().c_str(), vl_name.str().c_str(), debug);
        double score = learn::score(rscs[i], *exp, TL::TL_DOUBLE_MIN, *exp_weights, cout, 0);
        RLOG_INFO_STR_NAMED("RuleLearner2.learn", "  Repetition " << i << " gives score: " << score);
        scores.push_back(score);
      }
      int winner = TL::argmax(scores);
      double winner_score = scores[winner];
      assert(winner >= 0 && winner < static_cast<int>(scores.size()));

      RLOG_INFO_STR_NAMED("RuleLearner2.learn", "  Winner: " << winner << " gives score: " << scores[winner]);


      // sanity check
      if (isTrivial(rscs[winner]) && usesNoisyDefaultOutcome(rscs[winner])) {
          RLOG_WARN_STR_NAMED("RuleLearner2.learn",  "Symbol " << s->name << " resulted in a pure NID ruleset");

          /*
          // FIXME
          PRINT(*exp_weights);

          string output_root("/tmp/test_");

          stringstream learn_log_name;
          learn_log_name << output_root << "/DEBUG" << "_" << "_ruleset_learn.log";
          stringstream verbose_learn_log_name;
          verbose_learn_log_name << output_root << "/DEBUG" << "_" << "_ruleset_learn_FULL.log";

          RuleSetContainer rsc_debug;

          // HACK 1
          StateTransitionL exp2;
          for (int i = 0; i < exp->N; i++) {
              StateTransition* newE = new StateTransition;
              *newE = * (exp->elem(i));
              //newE->calcChanges();
              exp2.append(newE);
          }
          // hACK 2
          // LOAD symbolic experience
          // HACK make sure we do not have not enough constants
          uintA constants;
          for (uint i = 20; i < 40; i++) // FIXME too many
              constants.append(i);
          reason::setConstants(constants);

          // symbols
          relational::SymL symbols;
          relational::ArgumentTypeL types;
          relational::readSymbolsAndTypes(symbols, types, (output_root+"/symbols.log").c_str());
          //StateTransitionL exp2;
          //exp2 = StateTransition::read_SAS_SAS( (output_root+"/experiences.log").c_str() );

          //rsc_debug.init(exp, *exp_weights);
          rsc_debug.init(&exp2, *exp_weights);
//          if (crossvalidate_accross_actions) {
//              rsc_debug.p_experiences_crossvalidation = itc->second.p_experiences_crossvalidation;
//          }

          srand(12345);

          rlts_sim::World::saveSymbolicTransitions(output_root+"/experiences.log", exp2);
          vector<double> weights;
          uint j;
          FOR1D( (*exp_weights), j){
              weights.push_back(exp_weights->elem(j));
          }
          writeExperienceWeights(output_root+"/experience_weights.log", weights);

          //arr dummy_weights(exp->N);
//          dummy_weights.setUni(1.);
          learn::learn_rules(rsc_debug, exp2, *exp_weights, false,
          //learn::learn_rules(rsc_debug, *exp, *exp_weights, false,
          //learn::learn_rules(rsc_debug, *exp, dummy_weights, false,
                             learn_log_name.str().c_str(),
                             verbose_learn_log_name.str().c_str(),
                             3);

          learn::ScoreFunction* sf = relational::learn::ScoreFunction::getCurrentScoreFunction();
          cout << sf->name() << endl;
          learn::PasulaScoreFunction *psf  = dynamic_cast<learn::PasulaScoreFunction*>(sf);
          printf(" alpha=%.20f", psf->getAlpha());
          printf(" pmin=%.20f", psf->getPmin());
          printf(" pminNID=%.20f\n", psf->getPminNID());


          double score_debug = learn::score(rsc_debug, *exp, TL::TL_DOUBLE_MIN, *exp_weights, cout, 0);
          std::cout << "  DEBUG Repetition gives score: " << score_debug << std::endl;
          std::cout << " " << std::endl;

          cout << "abort?" << endl;
          getchar();
          throw 0;
          */
      }


      if (compress_before_learning) {
        RLOG_INFO_STR_NAMED("RuleLearner2.learn", "  Decompressing winner");
        decompressRuleSetContainer(expid_to_equivalence_class, equivalence_class_to_expid, its->second, itw->second, rscs[winner]);
        double score = learn::score(rscs[winner], its->second, TL::TL_DOUBLE_MIN, itw->second, cout, 0);
        if (fabs(winner_score - score) > 1e-5) {
          RLOG_ERROR_STR_NAMED("RuleLearner2.learn", "  Winner decompressed score " << score
                              << " does not match compressed score " << winner_score);
          throw 0;
        }
        delete exp;
        delete exp_weights;
      }

      itc->second.clear();
      itc->second.copyWithRules(rscs[winner]);
  //    itc->second = rscs[winner];

      // copy winner logs
      if (repeat > 1) {
        stringstream l_name, vl_name;
        l_name << log_name << s->name << "_" << winner << " .log";
        vl_name << log_name << s->name << "_FULL_" << winner << ".log";

        stringstream l_name_winner, vl_name_winner;
        l_name_winner << log_name << s->name << "_WIN_" << winner << ".log";
        vl_name_winner << log_name << s->name << "_FULL_WIN_" << winner << ".log";

        //assert(boost::filesystem::is_regular_file(from) && "Assume file is present");
        boost::filesystem::copy_file(l_name.str(),l_name_winner.str());
        boost::filesystem::copy_file(vl_name.str(),vl_name_winner.str());
      }

      // reset symbols list
      Symbol::mem__all_symbols.clear();
      Symbol::mem__all_symbols = mem__all_symbols_old;

    }
  }

  // set back to old constants
  reason::setConstants(old_constants);

  return true;
}

double RuleLearner2::score(Symbol* s, double cutting_threshold, std::ostream& out, uint DEBUG) {
  return score(s, experience_map[s], experience_weights[s], cutting_threshold, out, DEBUG);
}



bool RuleLearner2::usesNoisyDefaultOutcome(const RuleSetContainer& rsc) const {
  // get default rule
  int nid_exps = rsc.experiences_per_ruleOutcome.elem(0).elem(1).N;

  /*
  std::stringstream log;
  Rule* rule = rsc.rules.elem(0);
  log << "default rule:" << std::endl;
  rule->write(log);
  log << endl << " nid outcomes " << nid_exps << " (non-nid " << rsc.experiences_per_ruleOutcome.elem(0).elem(0).N << ")";
  RLOG_DEBUG_STR_NAMED("RuleLearner2.usesNoisyDefaultOutcome", log.str());
  */

  return nid_exps > 0;
}


bool RuleLearner2::isTrivial(Symbol* s, bool ignore_outcomes) const {
  std::map <Symbol*, StateTransitionL >::const_iterator ite = experience_map.find(s);
  std::map <Symbol*, RuleSetContainer >::const_iterator itc = containers.find(s);

  if (itc == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.isTrivial",
                       "RuleSetContainer for symbol " << s->name << " not found");
    return false;
  }

  const RuleSetContainer& rsc = itc->second;

  if (ite->second.N == 0) {
    RLOG_WARN_STR_NAMED("RuleLearner2.isTrivial",
                       "RuleSetContainer for symbol " << s->name << " has no experiences");
    return false;
  }

  return isTrivial(rsc, ignore_outcomes);

}

double RuleLearner2::score(Symbol* s,
                                              StateTransitionL& exp,
                                              arr& weights,
                                              double cutting_threshold,
                                              std::ostream& out, uint DEBUG) {
  if (exp.N == 0) {
    RLOG_WARN_STR_NAMED("RuleLearner2.score",
                       "No experiences given for " << s->name << "");
    return 0.0;
  }

  std::map <Symbol*, RuleSetContainer >::iterator itc = containers.find(s);
  //std::map <Symbol*, StateTransitionL >::iterator its = experience_map.find(s);
  //if (exp.N == 0) {  }

  if (itc == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.score",
                       "RuleSetContainer for symbol " << s->name << " not found");
    return -std::numeric_limits<double>::max();
  }

  // copy container (does not WORK!)
  // RuleSetContainer rsc = itc->second;
  RuleSetContainer& rsc = itc->second;

  double score = learn::score(rsc, exp, cutting_threshold, weights, out, DEBUG);

  return score;
}


bool RuleLearner2::usesNoisyDefaultOutcome(Symbol* s) const {

  std::map <Symbol*, RuleSetContainer >::const_iterator itc = containers.find(s);

  if (itc == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.isTrivial",
                       "RuleSetContainer for symbol " << s->name << " not found");
    return false;
  }

  const RuleSetContainer& rsc = itc->second;

  return usesNoisyDefaultOutcome(rsc);
}

bool RuleLearner2::isTrivial(const RuleSetContainer& rsc, bool ignore_outcomes) const {
  // more than default + 1 rule
  if (rsc.rules.num() > 2)
    return false;

  // only default rule
  if (rsc.rules.num() == 1)
    return true;

  // only one rule and context is empty?
  Rule* rule = rsc.rules.elem(1);
  if (rule->context.N == 0) {
    if (ignore_outcomes)
      return true;
    else {
      // if only one non-noise outcome which has no literals -> trivial
      return (rule->outcomes.N <= 2 && rule->outcomes(0).N == 0);
    }
  }

  return false;
}


double RuleLearner2::score(double cutting_threshold, std::ostream& out, uint DEBUG) {
  std::map<Symbol*,double> scores;
  return score(scores, cutting_threshold, out, DEBUG);
}

double RuleLearner2::score(relational::learn::ScoreFunction* sf, 
	double cutting_threshold, std::ostream& out, uint DEBUG) {
  std::map<Symbol*,double> scores;
  return score(scores, cutting_threshold, out, DEBUG, sf);
}

double RuleLearner2::score(map<Symbol*,double>& scores,
                           double cutting_threshold,
                           std::ostream& out, uint DEBUG,
                           relational::learn::ScoreFunction* sf) {

  double score=0.;

  std::map <Symbol*, RuleSetContainer >::iterator itc = containers.begin();
  for (; itc != containers.end(); itc++) {
    Symbol* s = itc->first;

    // it is safer to put this here in case less containers than clusters have been initialized
    std::map <Symbol*, StateTransitionL >::iterator its = experience_map.find(s);
    assert (its != experience_map.end());
    std::map <Symbol*, arr >::iterator itw = experience_weights.find(s);
    assert (itw != experience_weights.end());

    RLOG_INFO_STR_NAMED("RuleLearner2.score", "Scoring " << s->name);

    // no rules - not trained?
    if (itc->second.rules.num() == 0) {
      RLOG_WARN_STR_NAMED("RuleLearner2.score", s->name << " has no rules - has it been trained?");
      continue;
    }

    if (DEBUG >= 1)
      out << "Scoring " << s->name << endl;
    StateTransitionL& exp = its->second;
    if (sf != NULL) {
		scores[s] = sf->score(itc->second, exp, cutting_threshold, itw->second, out, DEBUG);
	} else {
		scores[s] = learn::score(itc->second, exp, cutting_threshold, itw->second, out, DEBUG);
	}
    score += scores[s];
  }
  return score;
}


void RuleLearner2::getMergedRuleSet(RuleSet& rs, bool copy_rules){
  // merge single rulesets into one
//  RuleSet rs;
  rs.clear();

  // add default rule
  rs.append(Rule::generateDefaultRule());

  RLOG_DEBUG_STR_NAMED("RuleLearner2.getMergedRuleSet",
                      "experience_map.size() " << experience_map.size());
  RLOG_DEBUG_STR_NAMED("RuleLearner2.getMergedRuleSet",
                      "containers.size() " << containers.size());
  RLOG_DEBUG_STR_NAMED("RuleLearner2.getMergedRuleSet",
                      "reverse_experience_index_map.size() " << reverse_experience_index_map.size());

  std::map <Symbol*, RuleSetContainer >::iterator itc = containers.begin();

  // only necessary if we want to build a combined rule set container
  for (; itc != containers.end(); itc++) {
    Symbol* s = itc->first;

    // it is safer to put this here in case less containers than clusters have been initialized
    std::map <Symbol*, StateTransitionL >::iterator its = experience_map.find(s);
    assert (its != experience_map.end());
    std::map <Symbol*, std::map<int, int> >::iterator itr = reverse_experience_index_map.find(s);
    assert (itr != reverse_experience_index_map.end());

//    Symbol* s = its->first;
    RuleSetContainer& rsc = itc->second;
//    StateTransitionL& exp = its->second;
    uint r;
    FOR1D_(rsc.rules,r) {
      if (r == 0) continue;
      Rule* newr;
      if (copy_rules) {
        newr = new Rule;
        newr->copyBody(* rsc.rules.elem(r));
      } else {
        newr = rsc.rules.elem(r);
      }
      rs.append(newr);
    }
  }
}

void RuleLearner2::calc_coveringRules_groundAction(RuleSet& r_grounds, std::vector<uint>& abstract_rules_idx,
                                                                      std::vector<uintA>& abstract2GroundMap,
                                                                      const RuleSet& all_abstract_rules, const SymbolicState& s, Literal* groundAction) {
  r_grounds.clear();
  abstract_rules_idx.clear();
  abstract2GroundMap.resize(all_abstract_rules.num());

  uint i, k;
  for (i=0; i<all_abstract_rules.num(); i++) {
    SubstitutionSet subs;
    if (reason::calcSubstitutions_rule_groundAction(subs, s, groundAction, all_abstract_rules.elem(i))) {
      abstract_rules_idx.push_back(i);
      FOR1D_(subs, k) {
        r_grounds.append(subs.elem(k)->apply(*all_abstract_rules.elem(i)));
        abstract2GroundMap[i].append(r_grounds.num()-1);
      }
    }
  }
}

void RuleLearner2::calcCoverage(relational::StateTransitionL& covered_experiences, uintA& covered_experiences_ids, MT::Array<relational::SubstitutionSet>& subsets, const relational::Rule* r, const relational::StateTransitionL& experiences) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"learn::calcCoverage [START]"<<endl;
  if (DEBUG>0) r->write(cout);
  covered_experiences.clear();
  covered_experiences_ids.clear();
  uint i;
#ifdef NO_DEICTICREFS_BY_NONBINARY
  //build rule that has the context of r with all non-binary symbols removed
  Rule ruleWithoutNonBinaries;
  ruleWithoutNonBinaries.action = r->action;
  bool containsNonBinaries = false;
  FOR1D(r->context, i) {
    if (r->context(i)->s->range_type == Symbol::binary)
      ruleWithoutNonBinaries.context.append(r->context(i));
    else containsNonBinaries = true;
  }
  if (DEBUG>0) {cout<<"ruleWithoutNonBinaries:"<<endl<<ruleWithoutNonBinaries;}
#endif

  FOR1D(experiences, i) {
    if (DEBUG>0) cout<<"ex "<<i<< " " << endl;
    if (DEBUG>1) experiences(i)->write(cout);

#ifdef NO_DEICTICREFS_BY_NONBINARY
    //Deictic refs may be ambigous due to the missing non-binary symbols.
    //If NO_DEICTICREFS_BY_NONBINARY is set these references are not permitted.
    if (containsNonBinaries) {
      SubstitutionSet subsNB;
      if (!reason::calcSubstitutions_rule_groundAction(subsNB, experiences(i)->pre, experiences(i)->action, &ruleWithoutNonBinaries)) {
        if (DEBUG>0) cout<<" --> 0  (NO_DEICTICREFS_BY_NONBINARY)"<<endl;
        continue;
      }
    }
#endif
    SubstitutionSet subs;
    if (reason::calcSubstitutions_rule_groundAction(subs, experiences(i)->pre, experiences(i)->action, r)) {
      covered_experiences.append(experiences(i));
      covered_experiences_ids.append(i);
      subsets.append(subs);
      if (DEBUG>0) cout<<" --> 1"<<endl;
    }
    else {
      if (DEBUG>0) cout<<" --> 0"<<endl;
    }
  }
  if (DEBUG>0) cout<<"learn::calcCoverage [END]"<<endl;
}

void RuleLearner2::info(ostream &out)
{
  std::map <Symbol*, RuleSetContainer >::iterator itc = containers.begin();
  for (; itc != containers.end(); itc++) {

    Symbol* s = itc->first;
    std::map <Symbol*, StateTransitionL >::iterator its = experience_map.find(s);
    assert(its != experience_map.end());

    out << "RuleSetContainer "  //<< (& (itc->second))
        << ", symbol: " << itc->first->name
        << ", rules:" << itc->second.rules.num()
        << "  experiences: " << its->second.N
        << endl;
  }
}

bool RuleLearner2::computeRelationalChangeDistribution()
{
  std::map <Symbol*, RuleSetContainer>::iterator it;
  for (it = containers.begin(); it != containers.end(); it++) {
    if (!computeRelationalChangeDistribution(it->first))
      return false;
  }

  return true;
}

bool RuleLearner2::computeRelationalChangeDistribution(Symbol* symbol)
{
  if (symbol == NULL) {
    RLOG_WARN_STR_NAMED("RuleLearner2.computeRelationalChangeDistribution",
                       "Requesting probability for NULL symbol");
    return false;
  }
  if (containers.find(symbol) == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.computeRelationalChangeDistribution",
                       "Requesting probability for unknown symbol " << symbol
                       << " (name=" << symbol->name << ")");
    return false;
  }
  assert (experience_map.find(symbol) != experience_map.end());

  // collect all outcomes
  Symbol* fakeUnaryActionSymbol = Symbol::get(MT::String("fakeUnaryAction"), 1, Symbol::primitive, Symbol::binary);

  /*
  RuleSetContainer& rsc = (containers.find(symbol))->second;
  RuleSet& rules = rsc.rules;

  vector< pair<SymbolicState, int> > relCh2numExp;

  for (uint j=1; j < rules.num(); j++) { // ignore NID rule
    Rule* r = rules.elem(j);
    uint o;
    FOR1D(r->outcomes, o) {
      if (r->outcomes(o).N == 0) continue;

      // collect the changes
      SymbolicState ss(r->outcomes(o));

      // create a fake action to make unification more effective
      uintA action_arg;
      action_arg.append(r->action->args(0)); // using only the first arg because the second is "invented" by rel change
      Literal* fakeAction = Literal::get(fakeUnaryActionSymbol, action_arg, 1 );
      ss.lits.append(fakeAction);
      // necessary call for data sanity
      Literal::getArguments(ss.state_constants, ss.lits);

      // count experience covered by this rule set
      int ne = rsc.experiences_per_ruleOutcome(j)(o).N;
      relCh2numExp.push_back( make_pair(ss, ne) );

      stringstream log;
      ss.write(log);
      log << " (#" << ne << " experiences)";
      RLOG_DEBUG_STR_NAMED("RuleLearner2.computeRelationalChangeDistribution",
                          "[" << symbol->name << "] Adding " << log.str());

    }
  }*/

  StateTransitionL& experiences = (experience_map.find(symbol))->second;

  vector< pair<SymbolicState, int> > relCh2numExp;
  vector< int > relCh2numExp_ids;

  for (uint i=0; i < experiences.N; i++) {
    StateTransition* st = experiences(i);

    if (st->changes.N == 0) continue;

    // collect the changes
    SymbolicState ss(st->changes);

    // create a fake action to make unification more effective
    uintA action_arg;
    action_arg.append(st->action->args(0)); // using only the first arg because the second is "invented" by rel change
    Literal* fakeAction = Literal::get(fakeUnaryActionSymbol, action_arg, 1 );
    ss.lits.append(fakeAction);
    // necessary call for data sanity
    Literal::getArguments(ss.state_constants, ss.lits);

    // count experience covered by this rule set
    //int ne = rsc.experiences_per_ruleOutcome(j)(o).N;
    int ne = 1;
    relCh2numExp.push_back( make_pair(ss, ne) );
    relCh2numExp_ids.push_back(i);

    stringstream log;
    ss.write(log);
    log << " (#" << ne << " experiences)";
    RLOG_DEBUG_STR_NAMED("RuleLearner2.computeRelationalChangeDistribution",
                        "[" << symbol->name << "] Adding " << log.str());

  }

  // compute equivalence classes
  int n = relCh2numExp.size();

  typedef vector<int> VecInt;
  typedef boost::unordered_set<int> SetInt;

  SetInt elements;
  int set_cnt;

  VecInt rank (n);
  VecInt parent (n);
  boost::disjoint_sets<int*,int*> ds(&rank[0], &parent[0]);

  for (int i=0; i<n; ++i) {
    ds.make_set(i);
    elements.insert(i);
  }

  for (int i=0; i < n; i++) {
    SymbolicState& s1 = relCh2numExp[i].first;
    for (int j=0; j < i; j++) {
      SymbolicState& s2 = relCh2numExp[j].first;

      SubstitutionSet subs;
      bool equal = reason::unifiableIgnoreCommonConstants(subs, s1, s2);

      if (equal) {
        ds.union_set(i, j);
      }
    }
  }
  set_cnt = (int)ds.count_sets(elements.begin(), elements.end());
  RLOG_DEBUG_STR_NAMED("RuleLearner2.computeRelationalChangeDistribution",
                  "[" << symbol->name << "] DONE! Number of sets: " << set_cnt);
  ds.normalize_sets(elements.begin(), elements.end()); // normalize set so that parent is always the smallest number

  map<int, int> relChEq2numExp; //  representative (= effect id) -> no experiences covered by this effect
  map<int, SymbolicState> relChEq2ss;  // representative (= effect id) -> effect
  map<int, vector<int> > relChEq2expIds;  // representative (= effect id) -> list of experience ids
  int sumExp=0;
  // iterate over all elements
  for (SetInt::const_iterator i = elements.begin(); i != elements.end(); ++i) {
    // find representative
    int repr_id = ds.find_set(*i);
    if (relChEq2numExp.find(repr_id) == relChEq2numExp.end()) {
      // representative has not yet occurred -> init relCh... maps
      relChEq2numExp[repr_id] = 0;
      assert (relChEq2ss.find(repr_id) == relChEq2ss.end());
      relChEq2ss[repr_id] = relCh2numExp[*i].first;
      relChEq2expIds[repr_id] = vector<int>();
    }
    // increment experience count
    relChEq2numExp[repr_id] += relCh2numExp[*i].second;
    sumExp += relCh2numExp[*i].second;
    // add *global* experience id
    RLOG_DEBUG_STR_NAMED("RuleLearner2.computeRelationalChangeDistribution", relCh2numExp[*i].first
                        << " ID local " << relCh2numExp_ids[*i] << " -> ID global " << reverse_experience_index_map[symbol][relCh2numExp_ids[*i]]);
    relChEq2expIds[repr_id].push_back( reverse_experience_index_map[symbol][relCh2numExp_ids[*i]] );
  }


  std::vector< std::pair<SymbolicState, double> > &rcp = relational_change_probabilities[symbol];
  rcp.clear();

  std::vector< std::pair<SymbolicState, vector<int> > > &rcexp = relational_change_experience_map[symbol];
  rcexp.clear();

//  rlts_sim::ActionComparator* ac = rlts_sim::ActionComparator::getInstance();
//  ac->getPrototype();

  stringstream log;

  map<int, int>::iterator rit = relChEq2numExp.begin();
  map<int, SymbolicState>::iterator sit = relChEq2ss.begin();
  map<int, vector<int> >::iterator eit = relChEq2expIds.begin();

  assert (relChEq2numExp.size() == relChEq2ss.size());
  assert (relChEq2numExp.size() == relChEq2expIds.size());

  for (; rit != relChEq2numExp.end() && sit != relChEq2ss.end() && eit != relChEq2expIds.end(); rit++, sit++, eit++) {
    assert(rit->second == static_cast<int>(eit->second.size()));
    // normalize count to obtain probability
    double p = rit->second/(double)sumExp;
    rcp.push_back(make_pair(sit->second, p));
    log << " " << sit->second.toString() << " (abs: " << rit->second << ") " << " -> " << p << endl;
    // copy experience idds
    rcexp.push_back(make_pair(sit->second, eit->second));
  }
  RLOG_DEBUG_STR_NAMED("RuleLearner2.computeRelationalChangeDistribution",
                      "RESULTING Distribution: " << log.str());

  return true;
}

bool RuleLearner2::getRelationalChangeProbability(Symbol* symbol, const StateTransition& st,
            double & probability, vector<int>& effect_experience_ids) const {

  static const double no_change_probability = 0.5;

  if (symbol == NULL) {
    RLOG_WARN_STR_NAMED("RuleLearner2.getRelationalChangeProbability",
                       "Requesting probability for NULL symbol");
    return false;
  }
  if (relational_change_probabilities.find(symbol) == relational_change_probabilities.end()) {
    RLOG_ERROR_STR_NAMED("RuleLearner2.getRelationalChangeProbability",
                       "Requesting probability for unknown symbol " << symbol
                       << " (name=" << symbol->name << ") - have you computed the probabilities already?");
    return false;
  }

  if (st.changes.N == 0) {
    probability = no_change_probability;
    return true;
  }

  if (isTrivial(symbol) && !usesNoisyDefaultOutcome(symbol)) {
    // make trivial clusters attractive for badly represented experiences
    probability = no_change_probability;
    return true;
  }

  // default: action not covered, i.e. 0 probability
  probability = 0.;
  effect_experience_ids.clear();

  Symbol* fakeUnaryActionSymbol = Symbol::get(MT::String("fakeUnaryAction"), 1, Symbol::primitive, Symbol::binary);
  Literal* fakeAction = Literal::get(fakeUnaryActionSymbol, TUP(st.action->args(0)), 1 );

  SymbolicState ss;
  ss.lits.append(st.changes);
  ss.lits.append(fakeAction);
  Literal::getArguments(ss.state_constants, ss.lits);

  std::vector< std::pair<SymbolicState, double> > rcp = relational_change_probabilities.find(symbol)->second;
  std::vector< std::pair<SymbolicState, double> >::iterator it = rcp.begin();

  std::vector< std::pair<SymbolicState, vector<int> > > rcexp = relational_change_experience_map.find(symbol)->second;
  std::vector< std::pair<SymbolicState, vector<int> > >::iterator eit = rcexp.begin();

  for (; it != rcp.end() && eit != rcexp.end(); it++, eit++) {
    SubstitutionSet subs;
    bool equal = reason::unifiableIgnoreCommonConstants(subs, it->first, ss);
    if (equal) {
      probability = it->second;
      effect_experience_ids = eit->second;

      // debug logging
      stringstream sseei;
      for (vector<int>::iterator eeit=effect_experience_ids.begin(); eeit != effect_experience_ids.end(); eeit++) {
        //effect_experience_ids.push_back(*eeit);
        sseei << *eeit << " ";
      }
      RLOG_DEBUG_STR_NAMED("RuleLearner2.getRelationalChangeProbability",
                          it->first.toString() << " matches (experiences = " << sseei.str() << ")");
      break;
    }
  }

  RLOG_DEBUG_STR_NAMED("RuleLearner2.getRelationalChangeProbability",
                      "Probability: " << probability);

  return true;

}

bool RuleLearner2::coverage(Symbol* symbol, const StateTransition& st,
                                                 Literal* action,
                                                 std::vector<uint>& covering_rules,
                                                 std::vector<uintA>& covering_outcomes_per_rule,
                                                 bool ignore_default_rule){
  if (symbol == NULL) {
    RLOG_WARN_STR_NAMED("RuleLearner2.coverage",
                       "Requesting probability for NULL symbol");
    return false;
  }
  if (containers.find(symbol) == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.coverage",
                       "Requesting probability for unknown symbol " << symbol
                       << " (name=" << symbol->name << ")");
    return false;
  }

  RuleSetContainer& rsc = (containers.find(symbol))->second;

  std::stringstream ss;
  rsc.write(ss);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.coverage",
                      "RuleSetContainer " << &rsc << "; "
                      << endl << ss.str() << endl
                      << "Symbol " << symbol->name << ", rules: " << rsc.rules.num());

  coverage(rsc, st, action, covering_rules, covering_outcomes_per_rule, ignore_default_rule);
  return true;
}

void RuleLearner2::coverage(RuleSetContainer& rsc, const StateTransition& st,
                                                 Literal * action,
                                                 std::vector<uint>& covering_rules,
                                                 std::vector<uintA>& covering_outcomes_per_rule,
                                                 bool ignore_default_rule) {

  RLOG_DEBUG_STR_NAMED("RuleLearner2.coverage",
                      "Action: " << *action);

  RuleSet coveringGroundRules;
  std::vector<uintA> abstract2GroundMap;
  calc_coveringRules_groundAction(coveringGroundRules, covering_rules,
                                  abstract2GroundMap,
                                  rsc.rules,
                                  st.pre,
                                  action);

  // remove default rule (always the first one)
  if (covering_rules.empty()) {

    calc_coveringRules_groundAction(coveringGroundRules, covering_rules,
                                    abstract2GroundMap,
                                    rsc.rules,
                                    st.pre,
                                    action);
    throw 0;
  }
  if (ignore_default_rule && covering_rules[0] == 0) {
    covering_rules.erase(covering_rules.begin());
  }

//  uint unique = relational::reason::calc_uniqueCoveringRule_groundRules_groundAction
//                (coveringGroundRules, st.pre, action);

//  if (unique == 0) {
//    RLOG_DEBUG_STR_NAMED("RuleLearner2.coverage",
//                        "Not unique (" << non_def_covering_rules.size() << " covering)");

//    // noisy default rule
//    return 0.0;
//  }

  covering_outcomes_per_rule.resize(covering_rules.size());
  int i=0;
  for (vector<uint>::iterator it = covering_rules.begin();
       it != covering_rules.end(); it++, i++) {
    if (*it == 0) continue; // ignore default rule
    uint sub;
    FOR1D(abstract2GroundMap[*it], sub) {
      uintA covering_outcomes;
      Rule* ground_rule = coveringGroundRules.elem(abstract2GroundMap[*it](sub));
      relational::reason::calc_coveringOutcomes(covering_outcomes,
                                                ground_rule,
                                                st.pre,
//                                                action,
                                                st.post);

      // last outcome is noise outcome -> remove
      assert (covering_outcomes.last() == ground_rule->outcomes.N-1);
      covering_outcomes.popLast();

      if (covering_outcomes.N == 0)
        RLOG_DEBUG_STR_NAMED("RuleLearner2.coverage", "No outcomes covers the experience");

//      SymbolicState suc;
//      calcSuccessorState(suc, st.pre, ground_rule->outcomes(o));

      // what if there are several substitutions yielding different outcomes?
      // -> this is not possible because if there are several substitutions
      // the rule is assumed to be not covering
      covering_outcomes_per_rule[i] = covering_outcomes;
      // TODO store substitution
    }
  }
}


double RuleLearner2::coveringProbability(RuleSetContainer& rsc, const StateTransition& st,
                                                            Literal * action)
{
  RLOG_DEBUG_STR_NAMED("RuleLearner2.coveringProbability",
                      "Action: " << *action);

  RuleSet coveringGroundRules;
  std::vector<uint> non_def_covering_rules;
  std::vector<uintA> abstract2GroundMap;
  calc_coveringRules_groundAction(coveringGroundRules, non_def_covering_rules,
                                  abstract2GroundMap,
                                  rsc.rules,
                                  st.pre,
                                  action);

  /*
  uint i;
  stringstream cgr_log;
  FOR1D_(coveringGroundRules, i) {
    coveringGroundRules.elem(i)->write(cgr_log);
  }
  RLOG_DEBUG_STR_NAMED("RuleLearner2.coveringProbability",
                      "coveringGroundRules: " << cgr_log.str());
  */

  uint unique = relational::reason::calc_uniqueCoveringRule_groundRules_groundAction
                (coveringGroundRules, st.pre, action);

  if (unique == 0) {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.coveringProbability",
                        "Not unique (" << non_def_covering_rules.size() << " covering)");

    // noisy default rule
    return 0.0;
  }

  stringstream rule_log;
  coveringGroundRules.elem(unique)->write(rule_log, false);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.coveringProbability",
                      "Uniquely covering rules: "
                      << rule_log.str());

  uintA covering_outcomes;
  relational::reason::calc_coveringOutcomes(covering_outcomes,
                                            coveringGroundRules.elem(unique),
                                            st.pre, st.post);

  if (covering_outcomes.N == 2) {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.coveringProbability",
                        "Unique outcome! ("
                        << coveringGroundRules.elem(unique)->probs(covering_outcomes(0))
                        << " "
                        << coveringGroundRules.elem(unique)->probs(covering_outcomes(1))
                        << ")"
                        );
    double p=coveringGroundRules.elem(unique)->probs(covering_outcomes(0));
    return p;
  } else {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.coveringProbability",
                        "Non-unique outcome! ("
                        << covering_outcomes.N
                        << " covering)");
    // no unique rule
    return 0.0;
  }

}

double RuleLearner2::coveringProbability(Symbol* symbol, const StateTransition& st,
                                                            Literal* action){
  if (symbol == NULL) {
    RLOG_WARN_STR_NAMED("RuleLearner2.coveringProbability",
                       "Requesting probability for NULL symbol");
    return 0.0;
  }
  if (containers.find(symbol) == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.coveringProbability",
                       "Requesting probability for unknown symbol " << symbol
                       << " (name=" << symbol->name << ")");
    return 0.0;
  }

  RuleSetContainer& rsc = (containers.find(symbol))->second;

  RLOG_DEBUG_STR_NAMED("RuleLearner2.coveringProbability",
                      "RuleSetContainer " << &rsc << "; "
                      "Symbol " << symbol->name << ", rules: " << rsc.rules.num());

  // "invent" action symbol
  //Literal* action = Literal::get(symbol, st.action->args, 1.);

  return coveringProbability(rsc, st, action);
}


/* ============================================================================================== */


double RuleLearner2::approximateNewRulePenalty(const RuleSetContainer& rsc, const StateTransition& experience, double additionalLiterals, double alpha) {
  // approximate Delta by multiplying the average context length+the outcome length with alpha_pen
  if (alpha == -1) {
    learn::PasulaScoreFunction* sf =
        dynamic_cast<learn::PasulaScoreFunction*>(learn::ScoreFunction::getCurrentScoreFunction());
    assert(sf);
    alpha = sf->getAlpha();
  }

  // how many changed lits do we have (will appear in outcome)
  int changes = experience.changes.N;

  // what is the average context length of our rules
  double avg_context_len=0;
  uint r;
  FOR1D_(rsc.rules, r) {
    Rule* rule = rsc.rules.elem(r);
    avg_context_len += rule->context.N;
  }
  avg_context_len /= rsc.rules.num()-1; // ignore default rule

  double pen = - alpha * (changes + avg_context_len + additionalLiterals) * 2.; // yes, times 2! be conservative man!
  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateNewRulePenalty",
                      "alpha=" << alpha << ", changes=" << changes << ", avg(context)=" << avg_context_len
                      << ", add="<< additionalLiterals);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateNewRulePenalty",
                      "RESULT: " << pen);

  return pen;
}

void RuleLearner2::getExperiencesCoveredByRule(StateTransitionL& exp_rule, arr& weights_rule,
                                                                  const RuleSetContainer& rsc_old, const uint rule) {
  uint e;
  FOR1D(rsc_old.experiences_per_rule(rule), e) {
    StateTransition* st = new StateTransition;
    *st = * rsc_old.p_experiences->elem(rsc_old.experiences_per_rule(rule)(e));
    double w = rsc_old.experience_weights.elem(rsc_old.experiences_per_rule(rule)(e));
    exp_rule.append(st);
    weights_rule.append(w);
  }
}

void RuleLearner2::createRuleSetContainerFromSingleRule(RuleSetContainer& rsc,
    StateTransitionL& experiences_rule, arr& weights_rule, const RuleSetContainer& rsc_old, const uint rule) {

  // default rule
  rsc.rules.append(Rule::generateDefaultRule());
  rsc.experiences_per_rule.append(uintA());
  rsc.experiences_per_ruleOutcome.append(MT::Array < uintA >());
  rsc.experiences_per_ruleOutcome(0).append(uintA()); // default empty
  rsc.experiences_per_ruleOutcome(0).append(uintA()); // default noise

  uintA experiences_per_rule(experiences_rule.N);
  uintA experience_id_map(rsc_old.p_experiences->N);
  for (uint i=0; i < experiences_rule.N; i++) {
    experiences_per_rule(i) = i;
    experience_id_map(rsc_old.experiences_per_rule(rule)(i)) = i;
  }
  // remap experiences to their new IDs; assuming relative ordering has not changed!
  MT::Array < uintA > experiences_per_ruleOutcome( rsc_old.experiences_per_ruleOutcome(rule).N );
  uint o, i;
  FOR1D(rsc_old.experiences_per_ruleOutcome(rule), o) {
    FOR1D(rsc_old.experiences_per_ruleOutcome(rule)(o), i) {
      experiences_per_ruleOutcome(o).append( experience_id_map(rsc_old.experiences_per_ruleOutcome(rule)(o)(i)) );
    }
  }

  rsc.experiences_per_rule.append(experiences_per_rule);
  rsc.experiences_per_ruleOutcome.append(experiences_per_ruleOutcome);

  rsc.nonDefaultRules_per_experience.resize(experiences_rule.N);
  for (uint i=0; i < experiences_rule.N; i++) {
    uintA lst; lst.append(1);
    rsc.nonDefaultRules_per_experience(i) = lst; // only one non-def rule
  }

  rsc.init(&experiences_rule, weights_rule);
  Rule* r = new Rule;
  r->copyBody(* (rsc_old.rules.elem(rule)));
  rsc.rules.append(r);
}

void RuleLearner2::createRuleSetContainerFromSingleRule(RuleSetContainer& rsc,
                                          StateTransitionL& experiences_rule, arr& weights_rule, Rule* rule,
                                          const MT::Array < uintA >& experiences_per_ruleOutcome) {
  // default rule
  rsc.rules.append(Rule::generateDefaultRule());
  rsc.experiences_per_rule.append(uintA());
  rsc.experiences_per_ruleOutcome.append(MT::Array < uintA >());
  rsc.experiences_per_ruleOutcome(0).append(uintA()); // default empty
  rsc.experiences_per_ruleOutcome(0).append(uintA()); // default noise

  uintA experiences_per_rule(experiences_rule.N);
  for (uint i=0; i < experiences_rule.N; i++) {
    experiences_per_rule(i) = i;
  }

  rsc.experiences_per_rule.append(experiences_per_rule);
  rsc.experiences_per_ruleOutcome.append(experiences_per_ruleOutcome);

  rsc.nonDefaultRules_per_experience.resize(experiences_rule.N);
  for (uint i=0; i < experiences_rule.N; i++) {
    uintA lst; lst.append(1);
    rsc.nonDefaultRules_per_experience(i) = lst; // only one non-def rule
  }

  rsc.init(&experiences_rule, weights_rule);
//  Rule* r = new Rule;
//  r->copyBody(* (rsc_old.rules.elem(rule)));
  rsc.rules.append(rule);
}

void RuleLearner2::createRuleSetContainerFromSingleRule(RuleSetContainer& rsc, const RuleSetContainer& rsc_old, const uint rule) {
  StateTransitionL* experiences_rule = new StateTransitionL;
  arr weights_rule;
  getExperiencesCoveredByRule(*experiences_rule, weights_rule, rsc_old, rule);
  createRuleSetContainerFromSingleRule(rsc, *experiences_rule, weights_rule, rsc_old, rule);
  delete experiences_rule;
}


void RuleLearner2::calculateUpdatedRuleOutcomeScoreDifference(RuleSetContainer& rsc,
                                                                                 const StateTransition& experience, double weight,
                                                                                 const uint& covering_rule,
                                                                                 double& probability, double& delta,
                                                                                 uint DEBUG) {

  stringstream ss;

  // add a new outcome
  StateTransitionL rule_experiences;
  arr rule_weights;
  getExperiencesCoveredByRule(rule_experiences, rule_weights, rsc, covering_rule);
  // FIXME
  ss << "rule_weights: "; rule_weights.write(ss); ss << endl;

  // old score
  RuleSetContainer rsc_old_rule;
  createRuleSetContainerFromSingleRule(rsc_old_rule, rule_experiences, rule_weights, rsc, covering_rule);
  double old_score = learn::score(rsc_old_rule, rule_experiences, TL::TL_DOUBLE_MIN, rule_weights);
  ss << "Old rule score: " << old_score << endl;

  // copy new experience and append to rule experiences
  StateTransition* new_exp = new StateTransition;
  * new_exp = experience;
  rule_experiences.append(new_exp);
  rule_weights.append(weight);

  // copy rule
  Rule* rule = new Rule;
  rule->copyBody(* (rsc.rules.elem(covering_rule)));
  uintA covered_experiences_ids(rule_experiences.N);
  for (uint i=0; i < rule_experiences.N; i++) covered_experiences_ids(i) = i;

  // = rsc.experiences_per_rule(covering_rule);
  //covered_experiences_ids.append(rule_experiences.N-1);

  // induce outcomes
  MT::Array< uintA > cEpo;
  //learn::learn_outcomes(rule, cEpo, rule_experiences, covered_experiences_ids, DEBUG);
  learnUniquelyCoveringOutcomes(rule, cEpo, rule_experiences, rule_weights, covered_experiences_ids, DEBUG);
  PRINTOUT(cEpo, ss);

  // new score
  RuleSetContainer rsc_new_rule;
  createRuleSetContainerFromSingleRule(rsc_new_rule, rule_experiences, rule_weights, rule, cEpo);
  stringstream debug_ruleset_single;
  rsc_new_rule.write(debug_ruleset_single);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.calculateUpdatedRuleOutcomeScoreDifference",
                      "SingleRule container: " << endl <<debug_ruleset_single.str());

  ss << " Induced Outcomes: " << endl;
  rsc_new_rule.rules.elem(1)->write(ss, false);

  double new_score = learn::score(rsc_new_rule, rule_experiences, TL::TL_DOUBLE_MIN, rule_weights);

  ss << " Induced Outcomes: " << endl;
  rsc_new_rule.rules.elem(1)->write(ss, false);

  ss << "New rule score: " << new_score << endl;

  // outcome double covering bug. assign correctly and recompute score
  /*
  uint o;
  // find empty outcome
  int emptyOutcome = -1;

//  rsc_new_rule.write(cout);
  FOR1D(rsc_new_rule.experiences_per_ruleOutcome(1), o) {
    if (o == rule->outcomes.N-1) // noise outcome
      break;
    if (rule->outcomes(o).N == 0) {
      emptyOutcome = 0;
      break;
    }
  }
  if (emptyOutcome != -1) {
    // no check for all non-empty outcomes whether experiences of empty outcome appear
    // and remove them
    FOR1D(rsc_new_rule.experiences_per_ruleOutcome(1), o) {
      if (emptyOutcome == o)
        continue;

      uintA& emptyOutcomeExps = rsc_new_rule.experiences_per_ruleOutcome(1)(emptyOutcome);
      uintA& exps = rsc_new_rule.experiences_per_ruleOutcome(1)(o);

      uint i;
      FOR1D(emptyOutcomeExps, i) {
        if (exps.findValue(emptyOutcomeExps(i)) != -1) {
          exps.removeValue(emptyOutcomeExps(i));
        }
      }
    }
//    rsc_new_rule.write(cout);

    ss << "OUTCOME induction is wrong. Correcting" << endl;

    ss << " Corrected Outcomes: " << endl;
    rule->write(ss, false);

    new_score = learn::score(rsc_new_rule, rule_experiences, TL::TL_DOUBLE_MIN);
    ss << "CORRECTED New rule score: " << new_score << endl;
  }

  */

  RLOG_DEBUG_STR_NAMED("RuleLearner2.calculateUpdatedRuleOutcomeScoreDifference",
                      ss.str());

  // check out the probability of the outcome that has covered the experience
  probability = -1;

  uint o;
  FOR1D(cEpo, o) {
    uintA& ce = cEpo(o);
    if (ce.findValue(rule_experiences.N-1) != -1) {
      probability = rule->probs(o);
      // special case: outcome is noise --> we set it to 0 (we do the same in learning based DeltaComputation)
      if (o == cEpo.N-1) {
        RLOG_DEBUG_STR_NAMED("RuleLearner2.calculateUpdatedRuleOutcomeScoreDifference",
                            "New outcome is noise outcome --> ignore real probability and set to 0");
        probability = 0.;
      }
      break;
    }
  }
  assert (probability != -1);

  // compute the delta between the rule scores
  delta = new_score - old_score;

  // clean up
  listDelete(rule_experiences);
}



/* ============================================================================================== */


void RuleLearner2::approximateScoreAfterAddingOneRuleCovers(RuleSetContainer& rsc,
                                                const StateTransition& experience, double weight,
                                                const uint& covering_rule, const uintA& covering_outcomes,
                                                double& probability, double& delta) {
  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                      "Case 1: Exactly one rule (" << covering_rule << ") context covers (s,a)");

  if (covering_outcomes.N >= 1) {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                        "Possibility a: Exactly one outcome (" << covering_outcomes(0) << ") covers s'");

    if (covering_outcomes.N > 1) {
      stringstream ss;
      uint i;
      FOR1D(covering_outcomes, i) {
        ss << covering_outcomes(i) << " ";
      }
      ss << endl;
      ss << "RULE: " << endl;
      rsc.rules.elem(covering_rule)->write(ss);
      ss << "EXPERIENCE: " << endl;
      experience.write(ss, 2);
      ss << " WEIGHT: " << weight << endl;

      if (experience.changes.N == 0) {
        // delete the outcome which has
      }

      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "More than one outcome is covering: " << ss.str());

    }

    stringstream ss;
    rsc.rules.elem(covering_rule)->write(ss);
    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                        "Rule: " << endl << ss.str());

    // recompute outcome after adding new experience and compute delta
    calculateUpdatedRuleOutcomeScoreDifference(rsc, experience, weight, covering_rule, probability, delta);

    // approximate Delta by log-likelihood of old rule
    //    return log(rsc.rules.elem(covering_rule)->probs(covering_outcomes(0)));

    return;
  }

  if (covering_outcomes.N == 0) {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                        "Possibility b: No outcome covers s'");

    // try changing the rule by adding literals
    StateTransitionL exp_for_rule;
    arr weights_for_rule;
    getExperiencesCoveredByRule(exp_for_rule, weights_for_rule, rsc, covering_rule);

    // copy new experience because experience is const
    StateTransition* experience_new = new StateTransition;
    *experience_new = experience;

    // show deltas
    delta = -std::numeric_limits<double>::max();
    probability = 0.;

    // score the old rule set container
    RuleSetContainer rsc_single_rule;
    createRuleSetContainerFromSingleRule(rsc_single_rule, exp_for_rule, weights_for_rule, rsc, covering_rule);
//    rsc_single_rule.write(cout); // FIXME

    double old_score = learn::score(rsc_single_rule, exp_for_rule, TL::TL_DOUBLE_MIN, weights_for_rule);

    // ---------------------------------------------
    // first attempt: try to find a distinctive rule
    Rule* new_rule = findDistinctiveRuleAddLiterals(rsc.rules.elem(covering_rule), exp_for_rule, experience);

    if (new_rule != NULL) {
      delta = approximateNewRulePenalty(rsc, experience, 1);
      probability = 1.0;
      stringstream log;
      new_rule->write(log, false);
      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "[DistinctiveRule] Rule differentiates between old experience and new one: "
                          << endl << log.str());
      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "[DistinctiveRule] best score so far." << endl
                          << "   delta_ee = " << delta << ", p=1.0");

    }


    // -------------------------------------------
    // second attempt: try to explain experience
    {
//      ExplainExperiences::setExperienceIndexMin(new_rule_id);
//      ExplainExperiences ee(true, false);
//      uint new_rule_id = exp_for_rule.N-1;
//      RuleSetContainer rules_2add(&exp_for_rule);
//      ee.findRules(rsc_single_rule, exp_for_rule, rules_2add);

      ExplainExperiences ee(true, false);
      Rule* new_rule = ee.explainExperience_deictic_ALL_DRs(experience_new);
      if (new_rule != NULL) {
        // check that other experiences are not covered by this new rule
        StateTransitionL covered_experiences;
        uintA covered_experiences_ids;
        arr covered_experiences_weights;
        //arr weights_for_rule(exp_for_rule.N); // assume all equal=1
        //exp_weights_for_rule.setUni(1);
        learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experiences_weights, new_rule, exp_for_rule, weights_for_rule);

        stringstream debug_rule;
        new_rule->write(debug_rule);
        RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                            "[ExplainExperiences] suggested rule " << endl
                            << debug_rule.str());

        if (covered_experiences.N == 0) {
          // rule is accepted

          // compare scores
          double delta_ee = approximateNewRulePenalty(rsc, experience, 1.);

          RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                              "[ExplainExperiences] delta_ee = " << delta_ee << ", p=1.0");

          if (delta_ee > delta) {
            RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                                "[ExplainExperiences] best score so far." << endl
                                << "   delta_ee = " << delta_ee << ", p=1.0");

            delta = delta_ee;
            probability = 1.0;
          }
        } else {
          RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                              "[ExplainExperiences] rule covers some old experiences. discarded");

        }
      }
    }

    // ---------------------------------------------
    // second attempt: try SplitOnLiterals operator
    {
    SplitOnLiterals splitOnLiterals;


    // append to experiences that need to be explained
    exp_for_rule.append(experience_new);
    weights_for_rule.append(weight);
    RuleSetContainer rules_2add(&exp_for_rule, weights_for_rule);
    splitOnLiterals.findRules(rsc_single_rule, exp_for_rule, weights_for_rule, rules_2add);

    double new_score=std::numeric_limits<double>::max();
    double delta_split=std::numeric_limits<double>::max();
    if (rules_2add.rules.num() == 0) {
      // FIXME
      // I think this can happen if the rule was not pruned
      // and contains all possible literals

      RLOG_ERROR_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "[SplitOnLiterals] Error applying SplitOnLiterals operator for computing delta.");

      stringstream lg;
      lg << "[SplitOnLiterals] " << endl;
      lg << "Old Rule Set: " << endl;
      rsc_single_rule.write(lg);
      lg << "rules2add: " << endl;
      rules_2add.write(lg);
      lg << "Experience: " << endl;
      experience_new->write(lg, 10);
      lg << "Weight: " << weight << endl;

      RLOG_ERROR_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          lg.str());

      RuleSetContainer rules_2add2(&exp_for_rule, weights_for_rule);
      SplitOnLiterals splitOnLiterals2;
      splitOnLiterals2.findRules(rsc_single_rule, exp_for_rule, weights_for_rule, rules_2add2);


    } else {
      // score the rule set container
      new_score = learn::score(rules_2add, exp_for_rule, TL::TL_DOUBLE_MIN, weights_for_rule);

      // FIXME
      //    uint learn_outcomes_DEBUG=0;
      //    if (new_score < -60000) {
      //      RLOG_WARN_STR_NAMED("RuleLearner2.approximateScoreAfterAdding", "kacke");
      //      learn::score(rules_2add, exp_for_rule, TL::TL_DOUBLE_MIN);
      //      learn_outcomes_DEBUG=10;
      //    }

      // compare scores
      delta_split = new_score - old_score;

      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "[SplitOnLiterals] new_score (" << new_score <<")"
                          << " - old_score (" << old_score << ") = " << delta_split);
    }

    if (delta_split > delta && rules_2add.rules.num() > 0) {

      stringstream log;
      rules_2add.write(log);
      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "[SplitOnLiterals] best score so far. Rules:"
                          << endl << log.str());

      delta = delta_split;
      probability = -1;
      // TODO compute the probability of the new rule
      uint r;
      uint new_rule_id = exp_for_rule.N-1;
      FOR1D_(rules_2add.rules, r) {
        if (rules_2add.experiences_per_rule(r).findValue(new_rule_id) != -1) {
          uint o;
          FOR1D(rules_2add.experiences_per_ruleOutcome(r), o) {
            if (rules_2add.experiences_per_ruleOutcome(r)(o).findValue(new_rule_id) != -1) {
              probability = rules_2add.rules.elem(r)->probs(o);
            }
          }
        }
      }
      assert (probability!=-1);
    }
    }

    // -------------------------------------------
    // last possibility: create a new (probabilistic) outcome and use its score
//    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
//                        "Alter outcomes only and check score");

    double delta_outcome;
    double probability_outcome;

    calculateUpdatedRuleOutcomeScoreDifference(rsc, experience, weight, covering_rule,
                                               probability_outcome, delta_outcome);

    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                        "[NewOutcome] Calculating new outcomes gives delta_outcome=" << delta_outcome <<","
                        << " prob_outcome=" << probability_outcome);

    if (delta_outcome > delta) {
      stringstream ss;


      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "[NewOutcome] best score so far. updated rule:" << ss.str());

      delta = delta_outcome;
      probability = probability_outcome;
    }

    // clean up
    listDelete(exp_for_rule);
    delete new_rule;

    return;
  }

  // these are all cases; (covering_outcomes.N > 0) cannot happen
  throw "(covering_outcomes.N > 0)?";

}

void RuleLearner2::approximateScoreAfterAddingSeveralRulesCover(RuleSetContainer& rsc, const StateTransition& experience,
    double weight, const std::vector<uint>& covering_rules, const std::vector<uintA>& covering_outcomes, double& probability, double& delta) {
  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                      "Case 2: Several rule contexts cover (s,a)");

  //vector<Rule*> specialized_rules(covering_rules.size(), NULL);
  vector<uint> specializable_rules;
  int non_specializable_rules=covering_rules.size();
  // for all rules check which of them can be specialized in order not to cover e
  for (uint r = 0; r < covering_rules.size(); r++) {
    StateTransitionL exp_for_rule;
    arr weights_for_rule;
    getExperiencesCoveredByRule(exp_for_rule, weights_for_rule, rsc, covering_rules[r]);
    if (exp_for_rule.N == 0)
      RLOG_ERROR_STR_NAMED("RuleLearner2.approximateScoreAfterAddingSeveralRulesCover",
                          "Rule " << covering_rules[r] << " does not cover any experiences - cannot be!");

    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                        "Trying to alter rule " << covering_rules[r]);
    Rule* rule = findDistinctiveRuleAddLiterals(rsc.rules.elem(covering_rules[r]), exp_for_rule, experience);

    if (rule != NULL) {
      specializable_rules.push_back(covering_rules[r]);
      non_specializable_rules--;
    }

    listDelete(exp_for_rule);
    delete rule;
  }

  // get alpha
  learn::PasulaScoreFunction* sf =
      dynamic_cast<learn::PasulaScoreFunction*>(learn::ScoreFunction::getCurrentScoreFunction());
  assert(sf);
  double alpha_pen = sf->getAlpha();


  if (non_specializable_rules <= 1) {

    int approxNewRuleAddLits=0;
    if (non_specializable_rules == 0) {
      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "Possibility a: all rules can be specialized to not cover e -> Alter one of them and create a new rule for e");
      approxNewRuleAddLits = 1;
    } else if (non_specializable_rules == 1) {
      RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                          "Possibility b: one rule r can be specialized to not cover e -> Alter all the others and assign e to r");
      approxNewRuleAddLits = covering_rules.size()-1;
    }

    // what would a new rule cost
    double delta1 = approximateNewRulePenalty(rsc, experience, approxNewRuleAddLits);

    // what would the outcome rearrangment cose
    double delta2 = -numeric_limits<double>::max();
    double delta2_probability = 0.;
    for(vector<uint>::iterator it = specializable_rules.begin(); it != specializable_rules.end(); it++) {
      double d, p;
      calculateUpdatedRuleOutcomeScoreDifference(rsc, experience, weight, *it, d, p);
      if (d > delta2) {
        delta2 = d;
        delta2_probability = p;
      }
    }
    // we have to add the penalty for the literals which are added to the specialized rules
    delta2 -= (specializable_rules.size()-1) * alpha_pen;

    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                        "delta(new rule)=" << delta1 << ", delta(new outcomes)=" << delta2);

    if (delta1 >= delta2) {
      // new rule gives better score
      delta = delta1;
      probability = 1.0;
    } else {
      // outcome rearrangement gives better score
      delta = delta2;
      probability = delta2_probability;
    }
  }
  else {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                      "Possibility c: at least 2 rules cover e. Subsume e under default rule");

    calculateUpdatedRuleOutcomeScoreDifference(rsc, experience, weight, 0, probability, delta);
  }
}

void RuleLearner2::approximateScoreAfterAddingNoRuleCovers(RuleSetContainer& rsc, StateTransitionL& old_experiences,
  arr& old_weights, const StateTransition& experience, double weight, double& probability, double& delta) {

  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                      "Case 3: No rule context covers (s,a)");

  // make sure that ExplainExperiences can build a rule which does not cover any
  // old experiences

  StateTransition new_experience = experience;
  double new_weight = weight;

  ExplainExperiences ee(true, false);
  Rule* new_rule = ee.explainExperience_deictic_ALL_DRs(&new_experience);
  if (new_rule == NULL) {
    RLOG_ERROR_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                        "ExplainExperiences cannot build a rule for 1 experience - that must not happen!");
    delta = -std::numeric_limits<double>::max();
    probability = 0.0;
    return;
  }

  stringstream debug_rule;
  new_rule->write(debug_rule);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                      "ExplainExperiences suggests rule " << endl
                      << debug_rule.str());

  // check that other experiences are not covered by this new rule
//  uint e;
//  FOR1D(experiences, e) {
  StateTransitionL covered_experiences;
  arr covered_weights;
  uintA covered_experiences_ids;
  //learn::calcCoverage(covered_experiences, covered_experiences_ids, new_rule, old_experiences);
  learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_weights, new_rule, old_experiences, old_weights);

  if (covered_experiences.N == 0) {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                        "ExplainExperiences' new rule covers only new experience. Accepted!");

    ///////////////////////////////////
    ///////////////////////////////////
    // check if there are more experiences in old default rule
    /*
    if (rsc.experiences_per_rule(0).N > 0) {
      // create a set of default experiences
      StateTransitionL defExperiences;
      uint i;
      FOR1D(rsc.experiences_per_rule(0), i) {
        defExperiences.append( old_experiences ( rsc.experiences_per_rule(0)(i)) );
      }

      // create a dummy set which contains the new experience
      StateTransitionL newExperienceSet;
      newExperienceSet.append(&new_experience);

      // create a single rule rule set container
      RuleSetContainer rsc_tmp;
      rsc_tmp.copyWithRules(rsc);
      rsc_tmp.addExperiences(&newExperienceSet);
      // learn outcomes for new rule
      MT::Array< uintA > cEpo;
      uintA covered_experiences_ids;
      covered_experiences_ids.append(rsc_tmp.p_experiences->N-1);
      learnUniquelyCoveringOutcomes(new_rule, cEpo, newExperienceSet, covered_experiences_ids);
      rsc_tmp.append(new_rule, covered_experiences_ids, cEpo);
      rsc_tmp.write(cout);

      // compute score of this rule set container
      StateTransitionL score_exps = *(rsc_tmp.p_experiences);
      double old_score = learn::score(rsc_tmp, score_exps, TL::TL_DOUBLE_MIN);

      MT::Array< RuleSetContainer > candidate_rulesC;

      double best_score = old_score;

      RuleSetContainer current;
      current.copyWithRules(rsc_tmp);
      current.p_experiences_copied = false; // memory leak
      while (true) {
        candidate_rulesC.clear();

        // try to drop literals until we can remove experiences from default rule
        // now invoke dropcontextliterals
        DropContextLiterals dropLits;
        dropLits.createRuleSets(current, score_exps, candidate_rulesC);

        //best_score = -numeric_limits<double>::max();
        int best_score_idx = -1;
        for (uint i = 0; i < candidate_rulesC.N; i++) {
          double score = learn::score(candidate_rulesC.elem(i), score_exps, TL::TL_DOUBLE_MIN);
          if (score > best_score) {
            best_score = score;
            best_score_idx = i;
          }
        }

        if (candidate_rulesC.N == 0 || best_score_idx == -1) {
          break;
        }

        current = candidate_rulesC(best_score_idx);

//        stringstream test;
//        current.write(test);
//        RLOG_WARN_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCover", test.str());
      }

      if (best_score > old_score) {
        stringstream test;
        current.write(test);
        RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCover", test.str());
      }

      delta = best_score - old_score;

      // TODO probability
    }
    ///////////////////////////////////
    ///////////////////////////////////
    else {
      delta = approximateNewRulePenalty(rsc, experience);
    }
    */

    delta = approximateNewRulePenalty(rsc, experience);
    probability = 1.0;

    return;
  }

  // --------------------------------------
  // Now it gets dirty: covered_experiences.N > 0

  stringstream debug_ce;
  PRINTOUT(covered_experiences_ids, debug_ce);

  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                      "ExplainExperiences coveres old experiences: " << debug_ce.str());

  // score old rule set container
  double old_score = learn::score(rsc, old_experiences, TL::TL_DOUBLE_MIN, old_weights);

  // get copy of old rsc
  RuleSetContainer rsc_tmp;
  rsc_tmp.copyWithRules(rsc);
  // add experience to default rule of old rsc
  StateTransitionL new_experience_set;
  new_experience_set.append(&new_experience);
  arr new_weights;
  new_weights.append(new_weight);
  rsc_tmp.addExperiences(&new_experience_set, new_weights);

  // are covered old experiences all in default rule?
  bool allCoveredOldExpAreInDefaultRule = true;
  uint i;
  FOR1D(covered_experiences, i) {
    if (!rsc.experiences_per_rule(0).contains(covered_experiences_ids(i))) {
      allCoveredOldExpAreInDefaultRule = false;
    }
  }
  // --> add these old experiences to new rule and recompute outcomes
  if (allCoveredOldExpAreInDefaultRule) {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                        "ExplainExperiences created a rule which covers NEW experience AND old default covered experiences. Adding rule to set and scoring.");

    int new_id = old_experiences.N;


    MT::Array< uintA > cEpo;
    covered_experiences.append(&new_experience);
    covered_experiences_ids.append( new_id );
    covered_weights.append(new_weight);
    learnUniquelyCoveringOutcomes(new_rule, cEpo, covered_experiences, covered_weights, covered_experiences_ids);

    rsc_tmp.append(new_rule, covered_experiences_ids, cEpo);
    rsc_tmp.recomputeDefaultRule();

    probability = -1;
    uint o;
    FOR1D(cEpo, o) {
      if (cEpo(o).contains(new_id)) {
        probability = new_rule->probs(o);
        break;
      }
    }
    assert (probability != -1);
  }
  // otherwise: we assume worst case, namely that new experience goes to default rule
  // --> rule is NOT accepted
  else {
    RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                        "ExplainExperiences creates a rule which covers old uniquely covered experiences. Rejecting rule.");

    if (experience.changes.N == 0)
      probability = rsc_tmp.rules.elem(0)->probs(0);
    else
      probability = rsc_tmp.rules.elem(0)->probs(1);

  }

  stringstream debug_rsc;
  rsc_tmp.write(debug_rsc);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                      debug_rsc.str());

  StateTransitionL all_experiences = * (rsc_tmp.p_experiences);
  arr all_weights = rsc_tmp.experience_weights;
  double new_score = learn::score(rsc_tmp, all_experiences, TL::TL_DOUBLE_MIN, all_weights);

  delta = new_score - old_score;

  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterAddingNoRuleCovers",
                      "old_score=" << old_score << ", new_score=" << new_score
                      << endl << "delta=" << delta << ", p=" << probability);

  // clean up
//  covered_experiences.popLast(); // new_experience will be deleted when destroying new_experience_set
//  listDelete(covered_experiences);

}


bool RuleLearner2::approximateScoreAfterAdding(const StateTransition& experience,
    double weight,
    const std::vector<uint>& covering_rules, const std::vector<uintA>& covering_outcomes,
    double& probability, double& delta) {

  // Case 0: weight is zero
  if (weight == 0) {
    RLOG_WARN_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                       "weight is 0 - will not contribute score");
    return false;
  }

  std::map <Symbol*, RuleSetContainer >::iterator itc = containers.find(experience.action->s);

  if (itc == containers.end()) {
    RLOG_WARN_STR_NAMED("RuleLearner2.approximateScoreAfterAdding",
                       "RuleSetContainer for symbol " << experience.action->s->name << " not found");
    return false;
  }

  std::map <Symbol*, StateTransitionL >::iterator ite = experience_map.find(experience.action->s);
  assert (ite != experience_map.end());

  std::map <Symbol*, arr >::iterator itw = experience_weights.find(experience.action->s);
  assert (itw != experience_weights.end());

  RuleSetContainer& rsc = itc->second;
  approximateScoreAfterAdding(rsc, ite->second, itw->second, experience, weight, covering_rules, covering_outcomes, probability, delta);

  return true;
}

void RuleLearner2::approximateScoreAfterAdding(RuleSetContainer& rsc, StateTransitionL& experiences,
    arr& weights,
    const StateTransition& experience, double weight,
    const std::vector<uint>& covering_rules, const std::vector<uintA>& covering_outcomes, double& probability, double& delta) {

  // Case 1: exactly one rule context covers the (s,a)
  if (covering_rules.size() == 1)
    approximateScoreAfterAddingOneRuleCovers(rsc, experience, weight, covering_rules[0], covering_outcomes[0], probability, delta );

  // Case 2: Several rule contexts cover (s,a)
  else if (covering_rules.size() > 1)
    approximateScoreAfterAddingSeveralRulesCover(rsc, experience, weight, covering_rules, covering_outcomes, probability, delta);

  // Case 3: no rule context covers (s,a)
  else if (covering_rules.size() == 0)
    approximateScoreAfterAddingNoRuleCovers(rsc, experiences, weights, experience, weight, probability, delta);

}


/* ============================================================================================== */
bool RuleLearner2::approximateScoreChangeAfterRemoving(Symbol* s, uint experienceId, double& delta_minus) {
  std::map <int, std::pair<Symbol*, int > >::iterator iteim = experience_index_map.find(experienceId);
  assert(iteim != experience_index_map.end());
  std::pair<Symbol*, int> p = iteim->second;
  assert (p.first == s); // make sure that you are checking removal of an experience which is actually covered by the model

  // get the "local" experience ID
  int e = p.second;

  RuleSetContainer& rsc = this->containers[s];

  RuleSetContainer rscWoExp;
  rscWoExp.copyWithRules(rsc);

  // manually copy experiences
  StateTransitionL* reduced_experiences = new StateTransitionL;
  arr reduced_weights;

  std::map <Symbol*, StateTransitionL >::iterator it = experience_map.find(s);
  assert (it != experience_map.end());
  const StateTransitionL& all_experiences = it->second;
  std::map <Symbol*, arr >::iterator itw = experience_weights.find(s);
  assert (itw != experience_weights.end());
  const arr& all_weights = itw->second;

  rscWoExp.nonDefaultRules_per_experience.clear();

  // copy experiences
  uint i;
  FOR1D(all_experiences, i) {
    if (static_cast<int>(i) != e) {
//      StateTransition* st = new StateTransition;
//      *st = * (all_experiences(i));
      StateTransition* st = all_experiences(i);
      reduced_experiences->append(st);
      reduced_weights.append(all_weights(i));
      rscWoExp.nonDefaultRules_per_experience.append(rsc.nonDefaultRules_per_experience(i));
    }
  }

  rscWoExp.p_experiences = reduced_experiences;
  rscWoExp.experience_weights = reduced_weights;
  rscWoExp.p_experiences_copied = false;

  // remove e from redundant lists
//  rscWoExp.nonDefaultRules_per_experience.remove(e);

  FOR1D(rscWoExp.experiences_per_rule, i) {
    // remove e from experiences_per_rule
    int idx = rscWoExp.experiences_per_rule(i).findValue(e);
    if (idx >= 0) {
      rscWoExp.experiences_per_rule.elem(i).remove(idx);
    }

    // now adapt indices of all experiences with index>e
    uint j;
    FOR1D(rscWoExp.experiences_per_rule(i), j) {
      if (rscWoExp.experiences_per_rule(i)(j) > e) {
        rscWoExp.experiences_per_rule(i)(j) = rscWoExp.experiences_per_rule(i)(j)-1;
      }
    }

    FOR1D(rscWoExp.experiences_per_ruleOutcome(i), j) {
      // remove e from experiences_per_ruleOutcome
      int idx = rscWoExp.experiences_per_ruleOutcome(i)(j).findValue(e);
      if (idx >= 0) {
        rscWoExp.experiences_per_ruleOutcome(i)(j).remove(idx);
      }

      uint k;
      FOR1D(rscWoExp.experiences_per_ruleOutcome(i)(j), k) {
        if (rscWoExp.experiences_per_ruleOutcome(i)(j)(k) > e) {
          rscWoExp.experiences_per_ruleOutcome(i)(j)(k) = rscWoExp.experiences_per_ruleOutcome(i)(j)(k)-1;
        }
      }
    }
  }

//  PRINT(rsc.experiences_per_rule);
//  PRINT(rsc.experiences_per_ruleOutcome);

//  PRINT(rscWoExp.experiences_per_rule);
//  PRINT(rscWoExp.experiences_per_ruleOutcome);

  double score_all = this->score(s, it->second, itw->second, TL::TL_DOUBLE_MIN, std::cout, 0);
  double score_reduced = learn::score(rscWoExp, *reduced_experiences, TL::TL_DOUBLE_MIN, reduced_weights, std::cout, 0);

//  double delta_minus_alt = score(s, removed_experience);
//  double score_all = score(s, all_experiences);
//  double score_reduced = score(s, reduced_experiences);
  delta_minus = score_reduced - score_all;
//  delta_minus = -score(s, removed_experience);


  RLOG_DEBUG_STR_NAMED("RuleLearner2.approximateScoreAfterRemoving",
                      s->name << std::endl
                      << "score with all experiences " << score_all << std::endl
                      << "score without experience " << experienceId << " " << score_reduced << std::endl
                      << "deltaMk_minus(" << experienceId << "->" << e << ") " << delta_minus << std::endl
//                      << "delta_minus_alt " << delta_minus_alt
                      );

  delete reduced_experiences;

  return true;

}

/* ============================================================================================== */

void RuleLearner2::learnUniquelyCoveringOutcomes(Rule* rule, MT::Array< uintA >& coveredExperiences_per_outcome, const StateTransitionL& coveredExperiences, const arr& coveredExperience_weights, const uintA& covered_experiences_ids, uint DEBUG) {

  //arr coveredExperience_weights(coveredExperiences.N);
  //coveredExperience_weights.setUni(1);
  learn::learn_outcomes(rule, coveredExperiences_per_outcome, coveredExperiences, covered_experiences_ids, coveredExperience_weights, DEBUG);
  stringstream ss;
  ss << "learn::learn_outcomes: ";
  PRINTOUT(coveredExperiences_per_outcome, ss);
  RLOG_DEBUG_STR_NAMED("RuleLearner2.learnUniquelyCoveringOutcomes", ss.str());

  // outcome double covering bug. assign correctly and recompute score
  uint o;
  // find empty outcome
  int emptyOutcome = -1;

//  rsc_new_rule.write(cout);
  FOR1D(coveredExperiences_per_outcome, o) {
    if (o == rule->outcomes.N-1) // noise outcome
      break;
    if (rule->outcomes(o).N == 0) {
      emptyOutcome = o;
      break;
    }
  }
  if (emptyOutcome != -1) {
    bool correction_needed = false;

    // no check for all non-empty outcomes whether experiences of empty outcome appear
    // and remove them
    FOR1D(coveredExperiences_per_outcome, o) {
      if (emptyOutcome == o)
        continue;

      uintA& emptyOutcomeExps = coveredExperiences_per_outcome(emptyOutcome);
      uintA& exps = coveredExperiences_per_outcome(o);

      uint i;
      FOR1D(emptyOutcomeExps, i) {
        if (exps.findValue(emptyOutcomeExps(i)) != -1) {
          exps.removeValue(emptyOutcomeExps(i));
          correction_needed = true;
        }
      }
    }

    // FIXME recompute outcome probabilities

//
//    boolA coverage;
//    arr probs;

//    __get_trimmed_outcomes(MT::Array< LitL >& outcomes, arr& probs, coverage, coveredExperiences, *rule);
//    rule->probs = probs;
//    learn::learn_parameters(, rule->probs);

//    stringstream ss;
//    ss << "OUTCOME induction is wrong. Corrected Outcomes: " << endl;
//    rule->write(ss, false);

//    RLOG_DEBUG_STR_NAMED("RuleLearner2.learnUniquelyCoveringOutcomes",
//                        ss.str());

    stringstream ss;
    if (correction_needed)
      ss << "OUTCOME induction is wrong. Corrected Outcomes: " << endl;

    PRINTOUT(coveredExperiences_per_outcome, ss);
    RLOG_DEBUG_STR_NAMED("RuleLearner2.learnUniquelyCoveringOutcomes",
                        ss.str());

  }

}

/**
 * Find an extension of the rule such that it covers
 * experiences but not distinctive_experience
 */
Rule* RuleLearner2::findDistinctiveRuleAddLiterals(const Rule* const rule, const StateTransitionL& experiences, const StateTransition& distinctive_experience) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"AddLits::findRules [START]"<<endl;

  if (DEBUG>0) {
    rule->write();
  }
  if (DEBUG>1) {
    uint i;
    FOR1D(experiences, i) {
      experiences(i)->write(cout, 2);
    }
  }

  LitL absentLiterals;
  Rule* newRule = NULL;

  // copy experience (because destructor of StateTransitionL=MT::Array deletes content)
  StateTransitionL d_experiences;
  StateTransition* st = new StateTransition;
  *st = distinctive_experience;
  d_experiences.append(st);

  // compute absent literals
  rule->getAbsentLiterals(absentLiterals);
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
//  absentLiterals.append(Literal::get("handle(X)"));

  // hack -- don't use complex reward-concepts [END]
  if (DEBUG>2) {
    cout << "Calculated absent literals for rule:"<<endl;
    rule->write(cout);
    cout<<"Absent literals: "; relational::write(absentLiterals);cout<<endl;
  }
  for (uint nextLiteral=0; nextLiteral<absentLiterals.N; nextLiteral++) {
    LitL wrapper;
    wrapper.append(absentLiterals(nextLiteral));
    if (DEBUG>1) {
      cout<<"Trying to insert ";absentLiterals(nextLiteral)->write(cout);cout<<"   into   ";
      relational::write(rule->context);cout<<endl;
    }
    if (Literal::nonContradicting(wrapper, rule->context)) {
      if (DEBUG>2) cout<<" --> Feasible and will be done."<<endl;
      newRule = new Rule;
      newRule->action = rule->action;
      newRule->context = rule->context;
      newRule->insertContext(absentLiterals(nextLiteral));

      // verify that old experiences are still covered
      StateTransitionL covered_experiences;
      uintA covered_experiences_ids;
      learn::calcCoverage(covered_experiences, covered_experiences_ids, newRule, experiences);
      if (covered_experiences.N < experiences.N) {
        if (DEBUG>1) cout<<"Does not cover all "<<covered_experiences.N<<" experiences and will be discarded."<<endl;
        delete newRule;
        newRule=NULL;
        continue;
      }

      // verify that new experience is NOT covered
      covered_experiences.clear();
      covered_experiences_ids.clear();
      learn::calcCoverage(covered_experiences, covered_experiences_ids, newRule, d_experiences);
      if (covered_experiences.N != 1) {
        // newRule is accepted!
        if (DEBUG>1) cout<<"  --> Will be kept! Does not cover new distinctive experience"<<endl;
        if (DEBUG>0) {newRule->write(cout);}
        break;
      }
      if (DEBUG>1) cout<<"Also covera distinctive experience and will be discarded."<<endl;
      delete newRule;
      newRule=NULL;
    }
  }

  // clean up
  listDelete(d_experiences);

  return newRule;
  if (DEBUG>0) cout<<"AddLits::findRules [END]"<<endl;
}



/* ============================================================================================== */


}

#endif //RULELEARNER2
