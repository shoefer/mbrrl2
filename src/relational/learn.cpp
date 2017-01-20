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
#include "learn.h"
#include "reason.h"

#include "score_functions.h"

namespace relational {

/************************************************
 *
 *     learn -- Major interface
 *
 ************************************************/


namespace learn {

// major learning parameters
double __alpha_PEN;
double __p_min;
double __p_min_noisyDefaultRule;  

// statistics
uintA num_so_improvements;
uintA num_so_applied;
arr so_improvements;
uintA so_successfulUsageHistory;
uintA so_UsageHistory;
arr scores;

MT::Array< SearchOperator* > __searchOperators;
// choosing the next operator
uint __so_choice_type = RULE_LEARNER__OP_CHOICE__RANDOM;
// int __so_current;
arr __so_priorWeights;

//caches values for a symbol that occur accross all experiences
std::map<relational::Symbol*, MT::Array<double> > __usedFunctionValues;

}


// **************************************************************** //

namespace learn {

double PasulaScoreFunction::score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, arr& experience_weights, std::ostream& out, uint DEBUG) {
//  DEBUG=3;//  uint DEBUG = 0;
  if (DEBUG > 0) {out << "SCORE [start]" << endl;}
  if (DEBUG > 1) {PRINTOUT(cutting_threshold, out);  rulesC.write(); }
  uint i;

  // (1) Penalty
  double penalty = 0.0;
  FOR1D_(rulesC.rules, i) {
    penalty += rulesC.rules.elem(i)->numberLiterals();
  }
  penalty *= __alpha_PEN;
  if (DEBUG > 0) {out << "alpha=" << __alpha_PEN << ", PEN=" << (penalty) << std::endl;}

  // (2) Log-Likelihood
  double loglikelihood = 0.0, exLik;
  if (DEBUG> 0) {out<<"Calculating likelihoods..."<<endl;}
  FOR1D(experiences, i) {
    if (DEBUG > 1) {out << "+++ ex " << i << ": ";  experiences(i)->action->write(out); out << " (weight=" << experience_weights(i) << ")"; out<<endl;}

    if (experience_weights(i) == 0.) {
      if (DEBUG > 1) {out << "Ignoring ex with w=0" <<endl;}
      continue;
    }

    const uintA& covering_rules = rulesC.nonDefaultRules_per_experience(i);
    if (DEBUG>1) {PRINTOUT(covering_rules, out);}
    // only one non-default rule covers
    if (covering_rules.N == 1) {
      Rule* rule = rulesC.rules.elem(covering_rules(0));
      if (DEBUG > 2) {
        out << "Use rule #" << covering_rules(0) << "  " << rulesC.experiences_per_rule(covering_rules(0)) << endl;
        rule->write(out);
      }
      uint o;
      MT::Array< uintA >& exs_per_out = rulesC.experiences_per_ruleOutcome(covering_rules(0));
      if (DEBUG>2) {PRINTOUT(exs_per_out, out);}
      exLik = 0.0;
      FOR1D(exs_per_out, o) {
        // Non-noise outcome
        if (o<exs_per_out.N-1  &&  exs_per_out(o).findValue(i) >= 0) {
          exLik += rule->probs(o);
//          exLik += pow(rule->probs(o), experience_weights(i));
          if (DEBUG>2) {out << "after o=" << o << ":  " << exLik << std::endl;}
        }
        // Noise-outcome:  always covers
        if (o == exs_per_out.N-1) {
          exLik += __p_min * rule->probs(o);
//          exLik += pow(__p_min * rule->probs(o), experience_weights(i));
          if (DEBUG>2) {out << "after o=" << o << ":  " << exLik << std::endl;}
        }
      }
      if (__p_min > 0)
        CHECK(exLik>0., "bad referencing  exLik="<<exLik);
    }
    // only default covers or more than one non-default covers
    // --> apply default rule
    else {
      if (DEBUG > 2) {
        out << " Using default rule." << endl;
      }
      Rule* default_rule = rulesC.rules.elem(0);
      if (rulesC.experiences_per_ruleOutcome(0)(0).findValue(i) >= 0)
        exLik = default_rule->probs(0) * (100.0 * __p_min_noisyDefaultRule);   // Avoid modeling persistence by means of noise default rule.
      else
        exLik = default_rule->probs(1) * __p_min_noisyDefaultRule;
      //exLik = pow(exLik, experience_weights(i));

      CHECK(exLik>0., "bad referencing  exLik="<<exLik);
    }
    loglikelihood += experience_weights(i) * log(exLik); // weight the experiences differently
    //loglikelihood += log(exLik); // weight the experiences differently
    //loglikelihood += log(experience_weights(i)) + log(exLik); // weight the experiences differently
    if (DEBUG>1) out<<" --> lik=" << exLik<<endl;
    if (DEBUG > 1) {PRINTOUT(loglikelihood, out);  }

    CHECK(!isnan(loglikelihood), "Pasula loglikelihood is NaN - that should not happen");
    CHECK(!isnan(penalty), "Pasula penalty is NaN - that should not happen");

    if (loglikelihood - penalty  <  cutting_threshold) {
      if (DEBUG>0) {
        PRINTOUT(loglikelihood - penalty, out);
        PRINTOUT(cutting_threshold, out);
        out<<"Score becomes too small... giving up."<<endl;
        out << "SCORE [end]" << endl;
      }
      return TL::TL_DOUBLE_MIN;
    }
  }

  if (DEBUG > 0) {out << "loglikelihood=" << loglikelihood << std::endl;  }

  double score = loglikelihood - penalty;
  if (DEBUG > 0) {
    out << (score) << std::endl;
    out << "SCORE [end]" << endl;
  }

  return score;
}

double PasulaScoreFunction::getAlpha() const {
    return __alpha_PEN;
}

double PasulaScoreFunction::getPmin() const {
    return __p_min;
}

double PasulaScoreFunction::getPminNID() const {
    return __p_min_noisyDefaultRule;
}

/*********/

double WeightVarianceScoreFunction::score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, arr& experience_weights, std::ostream& out, uint DEBUG) {
//  DEBUG = 3;
  if (DEBUG > 0) {out << "WeightedVarianceScoreFunction [start]" << endl;}

  double pasula_score = PasulaScoreFunction::score(rulesC, experiences, cutting_threshold, experience_weights, out, DEBUG);

  uint i;

  // (3) Average weight (reward/weight) variance of rules
  double penalty = 0.0;

  double weight_min=std::numeric_limits<double>::max();
  double weight_max=std::numeric_limits<double>::min();

  CHECK(rulesC.nonDefaultRules_per_experience.N == rulesC.p_experiences->N,
        "WeightVarianceScoreFunction can only be evaluated if rule set contains the correct experiences");

  // -------------------
  // [START] reward variance per rule

  if (this->__beta > 0) {
      FOR1D_(rulesC.rules, i) {
          //      Rule* rule = rulesC.rules.elem(i);
          uintA& exs_per_rule = rulesC.experiences_per_rule (i);

          // default rule might not cover anything
          if (exs_per_rule.N == 0)
              continue;

          if (DEBUG > 2) {
              out<<"-------------"<<endl;
              //			PRINT(rule->context);
              PRINT(exs_per_rule);
          }
          // collect statistics
          double mean=0., sum_sq=0.;

          uint j;

          double no_exs_per_rule=0;
          FOR1D(exs_per_rule, j) {
              //        double weight = experience_weights(j);

              uint e = exs_per_rule.elem(j);
              const StateTransition* st = rulesC.p_experiences->elem(e);
              double weight = st->reward;
              if (DEBUG > 2) {
                  out << "experience " << e << " --> " << weight << endl;
              }
              mean+=weight * experience_weights(e);
              sum_sq+=weight*weight * experience_weights(e);

              // store weight range
              if (weight < weight_min) weight_min = weight;
              if (weight > weight_max) weight_max = weight;

              no_exs_per_rule += experience_weights(e);
          }

          double var= sum_sq/no_exs_per_rule - pow(mean/no_exs_per_rule,2);
          if (DEBUG > 2)
              out << "rule " << i << " ---> " << (var) << endl;

          // default rule --> penalize worse
          if (i == 0) {
              // weigh NID rule variance by no of experiences?
              double exp_weight = 1;
              if (__weigh_nid_by_no_exp)
                  exp_weight = no_exs_per_rule;
              penalty+= (__nid_penalty*var * exp_weight) ;
          }
          else
              penalty+=var;
      }
      penalty *= this->__beta;
  }
  // [END]

  // -------------------
  // [START] reward variance per outcome
  double penalty_outcome = 0.0;
  if (__beta_outcome > 0) {
    FOR1D_(rulesC.rules, i) {
      Rule* rule = rulesC.rules.elem(i);

      uint o;
      FOR1D(rule->outcomes, o) {
        // ignore noise outcome
        //if (o == 0) continue;
        if (o == rule->outcomes.N-1) continue;

        uintA& exs_per_outcome = rulesC.experiences_per_ruleOutcome(i)(o);

        // default rule outcome might not cover anything
        if (exs_per_outcome.N == 0)
          continue;

        if (DEBUG > 2) {
          out<<"-------------"<<endl;
          //			PRINT(rule->context);
          PRINT(exs_per_outcome);
        }
        // collect statistics
        double mean=0., sum_sq=0.;

        double no_exs_per_outcome = 0;
        uint j;
        FOR1D(exs_per_outcome, j) {
          uint e = exs_per_outcome.elem(j);

          if (experience_weights(e) == 0) continue;

          const StateTransition* st = rulesC.p_experiences->elem(e);
          double weight = st->reward;
          if (DEBUG > 2) {
            out << "experience " << e << " --> " << weight << endl;
          }
          mean+=weight * experience_weights(e);
          sum_sq+=weight*weight * experience_weights(e);

          // store weight range
          if (weight < weight_min) weight_min = weight;
          if (weight > weight_max) weight_max = weight;

          no_exs_per_outcome += experience_weights(e);
        }

        if (no_exs_per_outcome == 0) continue;

        double var= sum_sq/no_exs_per_outcome - pow(mean/no_exs_per_outcome,2);
        if (DEBUG > 2)
          out << "rule " << i << " ---> " << (var) << endl;

        CHECK(!isnan(var), "var is nan! that should not happen");

        // default rule --> penalize worse
        if (i == 0) {
          // weigh NID rule variance by no of experiences?
          double exp_weight = 1;
          if (__weigh_nid_by_no_exp)
            exp_weight = no_exs_per_outcome;
          penalty_outcome+= (__nid_penalty*var * exp_weight) ;
        }
        else
          penalty_outcome+=var;
      }
    }
  }
  penalty_outcome *= this->__beta_outcome;
  // [END]

  // -------------------
  // [START] low reward likelihood
  double penalty2 = 0;
  if (__beta_negrew > 0) {
    if (DEBUG > 1) {
      out << " weight range: [" << weight_min << ", " << weight_max
          << "]" << endl;
    }

    FOR1D_(rulesC.rules, i) {
      uint o,j;

      //uintA& exs_per_rule = rulesC.experiences_per_rule (i);
      MT::Array < uintA > & exps_per_outcome = rulesC.experiences_per_ruleOutcome(i);

      FOR1D(exps_per_outcome, o) {
        double o_prob = rulesC.rules.elem(i)->probs.elem(o);
        double out_mean=0.;

        // compute mean
        if (exps_per_outcome.elem(o).N == 0)
          continue;
        double no_exps_per_outcome=0;
        FOR1D(exps_per_outcome.elem(o), j) {
          uint e = exps_per_outcome.elem(o).elem(j);
          const StateTransition* st = rulesC.p_experiences->elem(e);
          double weight = st->reward * experience_weights(e);
          if (DEBUG > 2) {
            out << "experience " << e << " --> " << weight << endl;
          }
          // normalize reward to 0..1
          out_mean+= fabs ((weight - weight_min)/(weight_max - weight_min));
          //            out_mean+= weight - weight_min;
          //cout << out_mean << endl;
          no_exps_per_outcome += experience_weights(e);
        }
        out_mean /= no_exps_per_outcome;

        penalty2 += (1-o_prob) * out_mean;
      }
    }
    penalty2 *= __beta_negrew;
  }
  // [END]

  //double score = pasula_score - penalty;
  double score = pasula_score - penalty - penalty2 - penalty_outcome;

  if (isnan(score)) {
    HALT("Weighted variance score is NaN - that should not happen");
  }

  if (DEBUG > 0) {
    out << "pasula score: " << pasula_score << endl;
    out << "weight variance (rule) penalty: " << penalty << " (weight=" << __beta << ")" << endl;
    out << "weight variance (outcome) penalty: " << penalty_outcome << " (weight=" << __beta_outcome << ")" << endl;
    out << "reward likelihood penalty: " << penalty2 << " (weight=" << __beta_negrew << ")" << endl;
    out << "SCORE: " << score << endl;
    out << "SCORE [end]" << endl;
  }

  return score;
}

/*********/
// use a default scoring function
ScoreFunction* __scoring_function = new PasulaScoreFunction;

void ScoreFunction::activateScoreFunction(ScoreFunction* sf, bool delete_old) {
  std::cout << "Activating scoring function: " << sf->name() << std::endl;
  if (__scoring_function == sf)
    return;
  if (delete_old)
    delete __scoring_function;
  __scoring_function = sf;
}

ScoreFunction* ScoreFunction::getCurrentScoreFunction() {
  return __scoring_function;
}

} // learn

// **************************************************************** //



void __initUsedFunctionValues(StateTransitionL& experiences) {
  learn::__usedFunctionValues.clear();

  uint i, j;
  //cache used function values
  FOR1D(experiences, i) {
    reason::derive(&experiences(i)->pre);
    reason::derive(&experiences(i)->post);
    FOR1D(experiences(i)->pre.lits, j) {
      if (experiences(i)->pre.lits(j)->s->range_type == Symbol::integers || experiences(i)->pre.lits(j)->s->range_type == Symbol::integer_set)
        learn::__usedFunctionValues[experiences(i)->pre.lits(j)->s].setAppend(experiences(i)->pre.lits(j)->value);
    }
    FOR1D(experiences(i)->post.lits, j) {
      if (experiences(i)->post.lits(j)->s->range_type == Symbol::integers || experiences(i)->post.lits(j)->s->range_type == Symbol::integer_set)
        learn::__usedFunctionValues[experiences(i)->post.lits(j)->s].setAppend(experiences(i)->post.lits(j)->value);
    }
  }
}

void __init_search_operators() {
  listDelete(learn::__searchOperators);
  learn::__so_priorWeights.clear();

  if (SO_WEIGHT__EXPLAIN_EXPERIENCES>0.0) {
    ExplainExperiences* explainEx = new ExplainExperiences(false, false);
    learn::__searchOperators.append(explainEx);
    learn::__so_priorWeights.append(SO_WEIGHT__EXPLAIN_EXPERIENCES);
  }
  
  if (SO_WEIGHT__EXPLAIN_EXPERIENCES_SLIM>0.0) {
    ExplainExperiences* explainEx_slim = new ExplainExperiences(true, false);
    learn::__searchOperators.append(explainEx_slim);
    learn::__so_priorWeights.append(SO_WEIGHT__EXPLAIN_EXPERIENCES_SLIM);
  }
  
  if (SO_WEIGHT__EXPLAIN_EXPERIENCES_SLIM>0.0) {
    ExplainExperiences* explainEx_slim_compare = new ExplainExperiences(true, true);
    learn::__searchOperators.append(explainEx_slim_compare);
    learn::__so_priorWeights.append(SO_WEIGHT__EXPLAIN_EXPERIENCES_SLIM_AND_COMPARING);
  }

  if (SO_WEIGHT__DROP_CONTEXT_LITERALS>0.0) {
    DropContextLiterals* dropPre = new DropContextLiterals();
    learn::__searchOperators.append(dropPre);
    learn::__so_priorWeights.append(SO_WEIGHT__DROP_CONTEXT_LITERALS);
  }
  
  if (SO_WEIGHT__DROP_CONTEXT_LITERALS_APPROX>0.0) {
    DropContextLiterals_approximativeVersion* dropPre_approx = new DropContextLiterals_approximativeVersion();
    learn::__searchOperators.append(dropPre_approx);
    learn::__so_priorWeights.append(SO_WEIGHT__DROP_CONTEXT_LITERALS_APPROX);
  }

  if (SO_WEIGHT__DROP_REFS>0.0) {
    DropReferences* dropRef = new DropReferences();
    learn::__searchOperators.append(dropRef);
    learn::__so_priorWeights.append(SO_WEIGHT__DROP_REFS);
  }
  
  if (SO_WEIGHT__DROP_RULES>0.0) {
    DropRules* dropRules = new DropRules();
    learn::__searchOperators.append(dropRules);
    learn::__so_priorWeights.append(SO_WEIGHT__DROP_RULES);
  }

  if (SO_WEIGHT__SPLIT_ON_LITS>0.0) {
    SplitOnLiterals* splitOnLits = new SplitOnLiterals();
    learn::__searchOperators.append(splitOnLits);
    learn::__so_priorWeights.append(SO_WEIGHT__SPLIT_ON_LITS);
  }

  if (SO_WEIGHT__ADD_LITS>0.0) {
    AddLiterals* addLit = new AddLiterals();
    learn::__searchOperators.append(addLit);
    learn::__so_priorWeights.append(SO_WEIGHT__ADD_LITS);
  }

  if (SO_WEIGHT__ADD_REFS>0.0) {
    AddReferences* addRef = new AddReferences();
    learn::__searchOperators.append(addRef);
    learn::__so_priorWeights.append(SO_WEIGHT__ADD_REFS);
  }
  
  if (SO_WEIGHT__SPLIT_ON_EQS>0.0) {
    SplitOnEqualities* splitOnEqs = new SplitOnEqualities();
    splitOnEqs->setUsedFunctionValues(learn::__usedFunctionValues);
    learn::__searchOperators.append(splitOnEqs);
    learn::__so_priorWeights.append(SO_WEIGHT__SPLIT_ON_EQS);
  }

  if (SO_WEIGHT__SPLIT_ON_INEQUALITIES > 0.0) {
    SplitOnInequalities* splitOnIeqs = new SplitOnInequalities();
    splitOnIeqs->setUsedFunctionValues(learn::__usedFunctionValues);
    learn::__searchOperators.append(splitOnIeqs);
    learn::__so_priorWeights.append(SO_WEIGHT__SPLIT_ON_INEQUALITIES);
  }
  
  if (SO_WEIGHT__GENERALIZE_EQS>0.0) {
    GeneralizeEquality* generalizeEq = new GeneralizeEquality();
    learn::__searchOperators.append(generalizeEq);
    learn::__so_priorWeights.append(SO_WEIGHT__GENERALIZE_EQS);
  }
  
  if (SO_WEIGHT__CHANGE_RANGE>0.0) {
    ChangeRange* changeRange = new ChangeRange();
    changeRange->setUsedFunctionValues(learn::__usedFunctionValues);
    learn::__searchOperators.append(changeRange);
    learn::__so_priorWeights.append(SO_WEIGHT__CHANGE_RANGE);
  }
  
  if (SO_WEIGHT__MAKE_INTVL>0.0) {
    MakeInterval* makeIntvl = new MakeInterval();
    makeIntvl->setUsedFunctionValues(learn::__usedFunctionValues);
    learn::__searchOperators.append(makeIntvl);
    learn::__so_priorWeights.append(SO_WEIGHT__MAKE_INTVL);
  }
  
  if (SO_WEIGHT__COMPARE_FUNCTIONVALUES>0.0) {
    CompareFunctionValues* compFVs = new CompareFunctionValues();
    learn::__searchOperators.append(compFVs);
    learn::__so_priorWeights.append(SO_WEIGHT__COMPARE_FUNCTIONVALUES);
  }
  
  if (SO_WEIGHT__SPLIT_ON_COMPARE_FUNCTIONVALUES>0.0) {
    SplitOnCompareFunctionValues* split_compFVs = new SplitOnCompareFunctionValues();
    learn::__searchOperators.append(split_compFVs);
    learn::__so_priorWeights.append(SO_WEIGHT__SPLIT_ON_COMPARE_FUNCTIONVALUES);
  }

  if (SO_WEIGHT__ADD_ABSTRACT_EQS > 0.0) {
    AddAbstractEquality* addAbstractEq = new AddAbstractEquality();
    learn::__searchOperators.append(addAbstractEq);
    learn::__so_priorWeights.append(SO_WEIGHT__ADD_ABSTRACT_EQS);
  }

  if (SO_WEIGHT__ABSTRACT_EQS > 0.0) {
    AbstractEquality* abstractEq = new AbstractEquality();
    learn::__searchOperators.append(abstractEq);
    learn::__so_priorWeights.append(SO_WEIGHT__ABSTRACT_EQS);
  }

  if (SO_WEIGHT__ADD_REFS_AND_ADD_LITS > 0.0) {
    AddReferencesAndAddLits* addReferencesAndAddLits = new AddReferencesAndAddLits();
    learn::__searchOperators.append(addReferencesAndAddLits);
    learn::__so_priorWeights.append(SO_WEIGHT__ADD_REFS_AND_ADD_LITS);
  }
  if (SO_WEIGHT__ADD_REFS_INDIRECT > 0.0) {
    AddIndirectReferences* addIndirectReferencesAndSplit = new AddIndirectReferences();
    learn::__searchOperators.append(addIndirectReferencesAndSplit);
    learn::__so_priorWeights.append(SO_WEIGHT__ADD_REFS_INDIRECT);
  }
  if (SO_WEIGHT__SPLIT_ON_LIT_CONJUNCTIONS > 0.0) {
    SplitOnLiteralConjunctions* splitOnLiteralConjunctions = new SplitOnLiteralConjunctions();
    learn::__searchOperators.append(splitOnLiteralConjunctions);
    learn::__so_priorWeights.append(SO_WEIGHT__SPLIT_ON_LIT_CONJUNCTIONS);
  }

  learn::num_so_applied.resize(learn::__searchOperators.N);
  learn::num_so_applied.setUni(0);
  learn::num_so_improvements.resize(learn::__searchOperators.N);
  learn::num_so_improvements.setUni(0);
  learn::so_improvements.resize(learn::__searchOperators.N);
  learn::so_improvements.setUni(0.0);
}


int chooseNextOperator(boolA& op_applicable) {
  if (sum(op_applicable) == 0)
    return -1;
  uint op = -1;
  // LINEAR VARIANT
  if (learn::__so_choice_type == RULE_LEARNER__OP_CHOICE__LINEAR) {
    static int op_linear = -1;
    op_linear++;
    op_linear %= learn::__searchOperators.N;
    op = op_linear;
  }
  // RANDOM VARIANT
  else if (learn::__so_choice_type == RULE_LEARNER__OP_CHOICE__RANDOM) {
    arr so_weights;
    so_weights = learn::__so_priorWeights;
    uint i;
    FOR1D_DOWN(learn::so_successfulUsageHistory, i) {
      if (learn::so_successfulUsageHistory.d0-SEARCH_OP_CHOICE__PAST_HORIZON>i)
        break;
      so_weights(learn::so_successfulUsageHistory(i)) += SEARCH_OP_CHOICE__PAST_WEIGHT;
    }
    FOR1D(op_applicable, i) {
      if (!op_applicable(i))
        so_weights(i) = 0.0;
    }
    op = TL::basic_sample(so_weights);
  }
  else {
    HALT("Undefined search operator choice procedure: " << learn::__so_choice_type)
  }
  // FINAL
  return op;
}


void learn::learn_rules(RuleSetContainer& rulesC, StateTransitionL& experiences, bool incremental, const char* logfile, const char* verbose_logfile, uint debug) {
  arr experience_weights(experiences.N);
  experience_weights.setUni(1.);
  learn_rules(rulesC, experiences, experience_weights, incremental, logfile, verbose_logfile, debug);
}


// Algorithm in Pasula, Zettlemoyer, Kaelbling, JAIR (2007), Figure 2
void learn::learn_rules(RuleSetContainer& rulesC, StateTransitionL& exp, arr& ew, bool incremental, const char* logfile, const char* verbose_logfile, uint debug) {
  uint DEBUG = debug; //  2 is a good choice
  uint i, k;
  FOR1D(exp, i) {
    Literal::sort(exp(i)->pre.lits);
  }
  Literal::sort(exp.last()->post.lits);
  
  StateTransitionL& experiences = exp;
  arr experience_weights = ew;
  if (incremental) {
    rulesC.addExperiences(&experiences, experience_weights);
    experience_weights = rulesC.experience_weights;
    experiences = *(rulesC.p_experiences);

    SearchOperator::setExperienceIndexMin(experiences.N - exp.N);

  } else {
    rulesC.clear();
    rulesC.init(&experiences, experience_weights);

    SearchOperator::resetExperienceIndexMin();
  }
  
  __initUsedFunctionValues(experiences);
  __init_search_operators();

  // verbose log writing
  std::ofstream ss;
  MT::open(ss, verbose_logfile);

  // Init default rule
  if (!rulesC.rules.hasDefaultRule()) {
    CHECK(rulesC.rules.num() == 0, "No default rule must mean that there is no rule at all");
    rulesC.rules.append(Rule::generateDefaultRule());
  }
  rulesC.recomputeDefaultRule();
  bool betterRulesFound = true;
  uint round = 0;
  double bestscore = score(rulesC, experiences, TL::TL_DOUBLE_MIN, experience_weights);
  scores.append(bestscore);
  if (DEBUG > 0) {
    ss << "RuleLearner::learn_rules [START]" << endl;
    //    PRINT (__alpha_PEN);
    //    PRINT(__p_min);
    //    PRINT(__p_min_noisyDefaultRule);
    ss << (__alpha_PEN);
    ss << (__p_min);
    ss << (__p_min_noisyDefaultRule);
    ss <<
          ss << "Number of training experiences " << experiences.N << endl;
    ss << "Default rule:" << endl;
    //     rulesC.rules.elem(0)->write(ss);
    rulesC.write(ss);
    ss << "SCORE = " << bestscore << endl;
  }
  if (DEBUG > 2) {ss << "Experiences:" << endl << experiences << endl;}
  // log writing
  std::ofstream log;
  MT::open(log, logfile);

  log<<"# Search Operators:"<<endl;
  FOR1D(__searchOperators, i) {
    log<<"# "<<i<<" "<<__searchOperators(i)->getName()<<endl;
  }
  log<<"#"<<endl;
  log << "# round  op  bestscore  #newRulesets  improvement"<<endl;
  
  uint MAX_STEPS = 100000;
  //   MT_MSG("RuleLearner: MAX_STEPS = "<<MAX_STEPS);
  boolA op_applicable;
  op_applicable.resize(__searchOperators.d0);
  op_applicable.setUni(true);
  bool so_useAgain=false;
  int op = 0;
  while (round++ < MAX_STEPS) {
    if (!so_useAgain) {
      op = chooseNextOperator(op_applicable);
      if (op < 0)
        break;
      if(__searchOperators(op)->isApproximator())
        __searchOperators(op)->reset_total_approximator();
    }
    
    if (op < 0)
      break;
    if (DEBUG > 0) {
      ss << "========== LEARN RULE-SET ROUND " << round << " ==========" << endl;
    }
    if (DEBUG > 2) {
      cout << "========== LEARN RULE-SET ROUND " << round << " ==========" << endl;
    }
    betterRulesFound = false;
    if (DEBUG > 1) {ss<<">>> Search operator ***"<<__searchOperators(op)->getName()<<"*** gives it a try. <<<"<<endl;}
    if (DEBUG > 1) {if(so_useAgain) ss<<"Using op again."<<endl; else ss<<"Using fresh operator."<<endl;}
    so_UsageHistory.append(op);
    num_so_applied(op)++;
    MT::Array< RuleSetContainer > set_of__rulesC_new;
    try {
      __searchOperators(op)->createRuleSets(rulesC, experiences, experience_weights, set_of__rulesC_new);
    } catch (...) {
        std::cerr << "Error creating rule set of search operator #" << op << std::endl;
        set_of__rulesC_new.clear();
    }

    if (set_of__rulesC_new.N == 0) {
      op_applicable(op) = false;
      if (so_useAgain) {
        so_useAgain = false;
        if (DEBUG>1) ss<<"Turning off search operator."<<endl;
      }
      //write log
      log << round << " " << op << " " << bestscore << " 0" << endl;
      if (DEBUG>1) {ss << "No new rules found." << endl;}
      continue;
    }
    //     if (DEBUG>1) {ss<<"Sanity checks of "<<set_of__rulesC_new.N<<" new rule-sets"<<endl;}
    //     FOR1D(set_of__rulesC_new, j) {
    //       if (DEBUG>1) {ss<<"Sanity check for new candidate rule-set #"<<j<<endl;}
    //       set_of__rulesC_new(k).sanityCheck();
    //     }
    arr new_scores;
    FOR1D(set_of__rulesC_new, k) {
      new_scores.append(
            score(set_of__rulesC_new(k), experiences, bestscore, experience_weights, ss, (DEBUG <= 1 ? 0 : DEBUG-2)));
    }
    if (DEBUG > 2) {
      ss << "Search operator found the following " << set_of__rulesC_new.N << " new rule-sets:" << endl;
      FOR1D(set_of__rulesC_new, k) {
        ss << "+++ New rule-set " << k << ": +++" << endl;
        set_of__rulesC_new(k).write(ss);
        ss << " --> score=" << new_scores(k);
        if (TL::areEqual(new_scores(k), TL::TL_DOUBLE_MIN))
          ss <<" (calculation had been aborted as intermediate score is already really bad...)";
        ss << endl;
      }
      ss << "Old best score: " << bestscore << endl;
      ss << "Scores of the new rule-sets: " << new_scores << endl;
    }
    uint maxIdx = new_scores.maxIndex();
    bool betterRulesFound_thisSearchOperator = false;
    // write log
    log << round << " " << op << " " << (new_scores(maxIdx) > bestscore ? new_scores(maxIdx) : bestscore) << " " << new_scores.N << " " << (new_scores(maxIdx) > bestscore ? (new_scores(maxIdx)-bestscore) : 0) << endl;
    // We did improve!
    if (new_scores(maxIdx) > bestscore) {
      // some statistic
      num_so_improvements(op)++;
      so_improvements(op) += new_scores(maxIdx)-bestscore;
      // algorithmic part
      rulesC = set_of__rulesC_new(maxIdx);
      bestscore = new_scores(maxIdx);
      betterRulesFound = true;
      betterRulesFound_thisSearchOperator = true;
      op_applicable.setUni(true); // now, all sos are applicable again since we have a rule-set change
      learn::so_successfulUsageHistory.append(op);
      so_useAgain = false; // here we go, efficiency! we don't need to use this SO again.
    }
    // Found rule-sets but not better ones
    else {
      if (!__searchOperators(op)->isApproximator())
        op_applicable(op) = false;
      else
        so_useAgain = true; // might produce nice results again
    }
    if (DEBUG > 0) {
      if (betterRulesFound_thisSearchOperator) {
        ss << "A new best rule-set was found:" << endl;
        rulesC.write(ss, false, false);
        ss << "SCORE = " << bestscore << endl;
        //rulesC.sanityCheck();  // TODO remove
        ss<<"Sanity check successful."<<endl;
      }
      else {ss << "No better rule-set was found (although I did my very best)." << endl;}
      if (round%10==0) {
        uintA covered_experiences_num;
        MT::Array< uintA > covered_experiences;
        arr responsibilities;
        rulesC.getResponsibilities(responsibilities, covered_experiences, covered_experiences_num);
        ss<<"Responsibilities after round "<<round<<":"<<endl;
        FOR1D(responsibilities, k) {
          ss << "[" << k << "] " << covered_experiences_num(k) << " " << responsibilities(k) << " " << covered_experiences(k) << endl;
        }
        ss << "-> Explanation of non-default rules: " << (1.0 - responsibilities(0)) << endl;
        uint total_noise = 0;
        FOR1D(rulesC.experiences_per_ruleOutcome, k) {
          total_noise += rulesC.experiences_per_ruleOutcome(k).last().N;
        }
        total_noise += responsibilities(0);
        ss <<"-> Non-noise explanations: "<<(1. - total_noise * 1.0 / experiences.N)<<endl;
      }
    }
    scores.append(bestscore);
    
    if (round % 50 == 0) rulesC.write("current_learned_rules.dat", false, false);
  }
  //   rulesC.sanityCheck();
  rulesC.recomputeDefaultRule();
  if (DEBUG>0) {ss<<"Trying to sort rules..."<<endl;}
  rulesC.sort();

  if (DEBUG>0) ss<<"Puh, that was it, now I can't find any better rules."<<endl;
  
  // LOGFILE WRITING [start]
  std::ofstream log_info;
  MT::String logfile_info;
  logfile_info << logfile << ".info";
  open(log_info, logfile_info);
  
  log_info<<"--- SYMBOLS ---"<<endl;
  log_info<<"*State symbols*"<<endl;
  SymL symbols_state;
  Symbol::get_state(symbols_state);
  writeSymbols(symbols_state, log_info);
  log_info<<"--- RULES ---"<<endl;
  rulesC.rules.write(log_info);
  log_info<<endl;
  // responsibilities
  uintA covered_experiences_num;
  MT::Array< uintA > covered_experiences;
  arr responsibilities;
  rulesC.getResponsibilities(responsibilities, covered_experiences, covered_experiences_num);
  uint num_experiences_explained_as_noise = 0;
  FOR1D(rulesC.experiences_per_ruleOutcome, k) {
    num_experiences_explained_as_noise += rulesC.experiences_per_ruleOutcome(k).last().N;
  }
  num_experiences_explained_as_noise += responsibilities(0);
  
  log_info<<"Responsibilities:"<<endl;
  FOR1D(responsibilities, k) {
    log_info << "[" << k << "] " << covered_experiences_num(k) << " " << responsibilities(k) << endl;
  }
  log_info << "-> Coverage of non-default rules: " << (100. * (1.0 - responsibilities(0))) << "%" << endl;
  log_info << "-> Non-noise explanations: "<<(100. * (1. - num_experiences_explained_as_noise * 1.0 / experiences.N))<<"% (1. - "<<num_experiences_explained_as_noise<<"/"<<experiences.N<<")"<<endl;
  log_info<<endl;
  log_info <<"--- STATISTICS ---"<<endl;
  log_info << "#rounds = " << (round-1) << endl;
  //     log_info << "Scores: "<<scores<<endl;
  log_info <<"SEARCH OPERATORS:  (#applied  #improve  ratio     improveTotal  improveStep)"<<endl;
  FOR1D(__searchOperators, i) {
    log_info<<"["<<i<<"] "<<"  "<<num_so_applied(i)<<"  "<<num_so_improvements(i)
           << "  " <<  ((1.0 * num_so_improvements(i))/num_so_applied(i))
           <<"  "<<so_improvements(i);
    if (so_improvements(i)>0.0)log_info<<"  "<<((1.0*so_improvements(i))/num_so_improvements(i));
    log_info<<"  "<<__searchOperators(i)->getName();
    log_info<<endl;
  }
  log_info << "History of successful SO applications: "<<endl<<learn::so_successfulUsageHistory<<endl;
  log_info << "History of SO applications: "<<endl<<so_UsageHistory<<endl;
  // LOGFILE WRITING [end]

  if (DEBUG>0) {
    ss<<"==================================================="<<endl;
    ss<<"BEST RULE-SET:"<<endl;
    rulesC.write(ss);
    ss<<"Responsibilities:"<<endl;
    FOR1D(responsibilities, k) {
      ss << "[" << k << "] " << covered_experiences_num(k) << " " << responsibilities(k) << " " << covered_experiences(k) << endl;
    }
    ss << "-> Coverage of non-default rules: " << (100. * (1.0 - responsibilities(0))) << "%" << endl;
    ss << "-> Non-noise explanations: "<<(100. * (1. - num_experiences_explained_as_noise * 1.0 / experiences.N))<<"% (1. - "<<num_experiences_explained_as_noise<<"/"<<experiences.N<<")"<<endl;
    ss<<"STATISTICS:"<<endl;
    ss << "#rounds = " << (round-1) << endl;
    ss << "Scores: "<<scores<<endl;
    ss<<"SEARCH OPERATORS:  (name  #applied  #improve  improveTotal  improveStep)"<<endl;
    FOR1D(__searchOperators, i) {
      ss<<"["<<i<<"] "<<"  "<<num_so_applied(i)<<"  "<<num_so_improvements(i)<<"  "<<so_improvements(i);
      if (so_improvements(i)>0.0) ss<<"  "<<(so_improvements(i)/num_so_improvements(i));
      ss<<"  "<<__searchOperators(i)->getName()<<endl;
    }
    ss << "History of successful SO applications: "<<learn::so_successfulUsageHistory<<endl;
    ss << "History of SO applications: "<<so_UsageHistory<<endl;
  }

  if (DEBUG > 0) ss << "RuleLearner::learn_rules [END]" << endl;
}


void learn::set_penalty(double alpha_PEN) {
  __alpha_PEN = alpha_PEN;
}


void learn::set_p_min(double p_min) {
  set_p_min(p_min, p_min);
}


void learn::set_p_min(double p_min, double p_min_noisyDefaultRule) {
  __p_min = p_min;
  __p_min_noisyDefaultRule = p_min_noisyDefaultRule;
}


void learn::set_ChoiceTypeSearchOperators(uint choice_type) {
  __so_choice_type = choice_type;
}


// changed by shoefer
double learn::score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, std::ostream& out, uint DEBUG) {
  return __scoring_function->score(rulesC, experiences, cutting_threshold, out, DEBUG);
}

double learn::score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, arr& experience_weights, std::ostream& out, uint DEBUG) {
  return __scoring_function->score(rulesC, experiences, cutting_threshold, experience_weights, out, DEBUG);
}

/*
double learn::score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, std::ostream& out, uint DEBUG) {
  arr experience_weights(experiences.N);
  experience_weights.setUni(1.0);
  return score(rulesC, experiences, cutting_threshold, experience_weights, out, DEBUG);
}


double learn::score(RuleSetContainer& rulesC, StateTransitionL& experiences, double cutting_threshold, arr& experience_weights, std::ostream& out, uint DEBUG) {
//  uint DEBUG = 0;
  if (DEBUG > 0) {out << "SCORE [start]" << endl;}
  if (DEBUG > 1) {PRINTOUT(cutting_threshold, out);  rulesC.write(); }
  uint i;
  
  // (1) Penalty
  double penalty = 0.0;
  FOR1D_(rulesC.rules, i) {
    penalty += rulesC.rules.elem(i)->numberLiterals();
  }
  penalty *= __alpha_PEN;
  if (DEBUG > 0) {out << "alpha=" << __alpha_PEN << ", PEN=" << (penalty) << std::endl;}
  
  // (2) Log-Likelihood
  double loglikelihood = 0.0, exLik;
  if (DEBUG> 0) {out<<"Calculating likelihoods..."<<endl;}
  FOR1D(experiences, i) {
    if (DEBUG > 1) {out << "+++ ex " << i << ": ";  experiences(i)->action->write(out);  out<<endl;}
    const uintA& covering_rules = rulesC.nonDefaultRules_per_experience(i);
    if (DEBUG>1) {PRINTOUT(covering_rules, out);}
    // only one non-default rule covers
    if (covering_rules.N == 1) {
      Rule* rule = rulesC.rules.elem(covering_rules(0));
      if (DEBUG > 2) {
        out << "Use rule #" << covering_rules(0) << "  " << rulesC.experiences_per_rule(covering_rules(0)) << endl;
        rule->write(out);
      }
      uint o;
      MT::Array< uintA >& exs_per_out = rulesC.experiences_per_ruleOutcome(covering_rules(0));
      if (DEBUG>2) {PRINTOUT(exs_per_out, out);}
      exLik = 0.0;
      FOR1D(exs_per_out, o) {
        // Non-noise outcome
        if (o<exs_per_out.N-1  &&  exs_per_out(o).findValue(i) >= 0) {
          exLik += rule->probs(o);
          if (DEBUG>2) {out << "after o=" << o << ":  " << exLik << std::endl;}
        }
        // Noise-outcome:  always covers
        if (o == exs_per_out.N-1) {
          exLik += __p_min * rule->probs(o);
          if (DEBUG>2) {out << "after o=" << o << ":  " << exLik << std::endl;}
        }
      }
      CHECK(exLik>0., "bad referencing  exLik="<<exLik);
    }
    // only default covers or more than one non-default covers
    // --> apply default rule
    else {
      if (DEBUG > 2) {
        out << " Using default rule." << endl;
      }
      Rule* default_rule = rulesC.rules.elem(0);
      if (rulesC.experiences_per_ruleOutcome(0)(0).findValue(i) >= 0)
        exLik = default_rule->probs(0) * (100.0 * __p_min_noisyDefaultRule);   // Avoid modeling persistence by means of noise default rule.
      else
        exLik = default_rule->probs(1) * __p_min_noisyDefaultRule;
    }
    loglikelihood += experience_weights(i) * log(exLik); // weight the experiences differently
    if (DEBUG>1) out<<" --> lik=" << exLik<<endl;
    if (DEBUG > 1) {PRINTOUT(loglikelihood, out);  }

    
    if (loglikelihood - penalty  <  cutting_threshold) {
      if (DEBUG>0) {
        PRINTOUT(loglikelihood - penalty, out);
        PRINTOUT(cutting_threshold, out);
        out<<"Score becomes too small... giving up."<<endl;
        out << "SCORE [end]" << endl;
      }
      return TL::TL_DOUBLE_MIN;
    }
  }

  if (DEBUG > 0) {out << "loglikelihood=" << loglikelihood << std::endl;  }

  double score = loglikelihood - penalty;
  if (DEBUG > 0) {
    out << (score) << std::endl;
    out << "SCORE [end]" << endl;
  }
  
  return score;
}
*/






/************************************************
 *
 *     RuleSetContainer
 *
 *     Efficiency wrapper for rule-sets
 *     --> stores the coverage of experiences
 *
 ************************************************/


RuleSetContainer::RuleSetContainer() : p_experiences(NULL), p_experiences_copied(false), rules_copied(false),p_experiences_crossvalidation(NULL) {
  init(NULL);
}


RuleSetContainer::RuleSetContainer(const StateTransitionL* _p_experiences, const arr& experience_weights)
    : p_experiences(NULL), p_experiences_copied(false), rules_copied(false),p_experiences_crossvalidation(NULL) {
  init(_p_experiences, experience_weights);
}

RuleSetContainer& RuleSetContainer::operator= (const RuleSetContainer& other){
  this->rules = other.rules;
  this->nonDefaultRules_per_experience = other.nonDefaultRules_per_experience;
  this->experiences_per_rule = other.experiences_per_rule;
  this->experiences_per_ruleOutcome = other.experiences_per_ruleOutcome;
  this->experience_weights = other.experience_weights;

  // deep copy experiences
  //this->p_experiences = new StateTransitionL;
  //*(this->p_experiences) = * (other.p_experiences);
  // shallow copy experiences
  if (this->p_experiences != other.p_experiences) {
    if (this->p_experiences_copied)
      delete this->p_experiences;
    this->p_experiences = other.p_experiences;
    this->p_experiences_copied = false;
  }
  this->p_experiences_crossvalidation = other.p_experiences_crossvalidation;

  return *this;
}

void RuleSetContainer::copyWithRules(const RuleSetContainer& other) {
  this->nonDefaultRules_per_experience = other.nonDefaultRules_per_experience;
  this->experiences_per_rule = other.experiences_per_rule;
  this->experiences_per_ruleOutcome = other.experiences_per_ruleOutcome;
  this->experience_weights = other.experience_weights;

  // deep copy experiences
  //this->p_experiences = new StateTransitionL;
  //*(this->p_experiences) = * (other.p_experiences);
  // shallow copy experiences
  if (this->p_experiences != other.p_experiences) {
    if (this->p_experiences_copied)
      delete this->p_experiences;
    this->p_experiences = other.p_experiences;
    this->p_experiences_copied = false;
  }

  this->p_experiences_crossvalidation = other.p_experiences_crossvalidation;

  // copy manually
  uint r;
  FOR1D_(other.rules, r) {
    Rule* newRule = new Rule;
    newRule->copyBody(* other.rules.elem(r));
    rules.append(newRule);
  }
}


RuleSetContainer::~RuleSetContainer() {
  if (this->p_experiences_copied)
    delete this->p_experiences;

  if (rules_copied) {
      uint r;
      FOR1D_(rules, r) {
          delete rules.elem(r);
      }
  }
}

void RuleSetContainer::init(const StateTransitionL* _p_experiences, const arr&  _experience_weights) {
  //this->p_experiences = _p_experiences;
  if (_p_experiences != NULL) {
    experience_weights = _experience_weights;
    if (experience_weights.N == 0 && _p_experiences->N > 0) {
      experience_weights.resize(_p_experiences->N);
      experience_weights.setUni(1.);
    }
    if (_p_experiences->N != experience_weights.N) {
      HALT("RuleSetContainer.init expects same number of experiences and weights");
    }
    //this->p_experiences = new StateTransitionL;
    //*(this->p_experiences) = *_p_experiences; // copy
    this->p_experiences = _p_experiences; // do not copy
    //this->experience_weights = experience_weights;
    this->nonDefaultRules_per_experience.resize(this->p_experiences->N);
    this->p_experiences_copied = false;
  }
}

void RuleSetContainer::addExperiences(const StateTransitionL* _new_experiences, const arr& _experience_weights) {
  if (this->p_experiences == NULL || this->p_experiences->N == 0) {
    HALT("RuleSetContainer::addExperiences only applicable to already learnt rule sets");
  }
  if (_new_experiences == NULL || _new_experiences->N == 0)
    HALT("RuleSetContainer::addExperiences did not get new experiences");

  arr ew =_experience_weights;
  if (ew.N == 0 && _new_experiences->N > 0) {
    ew.resize(_new_experiences->N);
    ew.setUni(1.);
  }

  uint no_old_experiences = this->p_experiences->N;

  nonDefaultRules_per_experience.resizeCopy(no_old_experiences + _new_experiences->N);

  // copy experiences
  StateTransitionL* new_p_experiences = new StateTransitionL;
  *new_p_experiences = *p_experiences;
  new_p_experiences->append(*_new_experiences);
  p_experiences_copied = true;
  p_experiences = static_cast<const StateTransitionL*>(new_p_experiences);
  CHECK(p_experiences != NULL, "error in addExperiences");

  //this->p_experiences->append(*_new_experiences);
  this->experience_weights.append(ew);

  // now check whether new experiences can be uniquely explained by existing rules
  uint r;
  FOR1D_(rules, r) {
    if (r == 0) continue; // ignore default rule
    StateTransitionL covered_new_experiences; // ignored
    uintA covered_new_experiences_ids;
    arr covered_experience_weights;
    learn::calcCoverage(covered_new_experiences, covered_new_experiences_ids, covered_experience_weights, rules.elem(r), *_new_experiences, this->experience_weights);

    //      calc_coveringOutcomes(uintA& covering_outcomes, Rule* abstractRule, const SymbolicState& pre, Literal* groundAction, const SymbolicState& post)

    uint i;
    FOR1D(covered_new_experiences_ids,i) {
      // adjust ID
      uint id = i+no_old_experiences;

      SubstitutionSet subs;

      // if there are more than one substitution, reason::calc_coveringOutcomes will crash -> check subs
      reason::calcSubstitutions_rule_groundAction(subs, _new_experiences->elem(i)->pre, _new_experiences->elem(i)->action, rules.elem(r));

      if (subs.num() != 1) {
          continue;
      }

      experiences_per_rule(r).append(id);
      nonDefaultRules_per_experience(id).append(r);

      // find the covering outcome
      uintA covering_outcomes;
      reason::calc_coveringOutcomes(covering_outcomes, rules.elem(r),
                                    _new_experiences->elem(i)->pre, _new_experiences->elem(i)->action, _new_experiences->elem(i)->post);

      if (covering_outcomes.N == 2) { // one regular + one noise -> ignore noise
          experiences_per_ruleOutcome(r)(covering_outcomes(0)).append(id);
      } else {
          uint o;
          FOR1D(covering_outcomes, o) {
              experiences_per_ruleOutcome(r)(covering_outcomes(o)).append(id);
          }
      }

    }

    // recompute the outcome probabilities of the rules
    if (covered_new_experiences_ids.N > 0) {
      Rule* rule = rules.elem(r);
      // rule->outcomes does not have to be deleted because LitL contains
      // Literal "singletons"
      StateTransitionL covered_experiences;
      arr covered_experience_weights;
      uint i;
      FOR1D(experiences_per_rule(r), i) {
          covered_experiences.append(p_experiences->elem(experiences_per_rule(r)(i)));
          covered_experience_weights.append(this->experience_weights(experiences_per_rule(r)(i)));
      }

      MT::Array< uintA > experiences_per_outcome;
      learn::learn_outcomes(rule, experiences_per_outcome, covered_experiences, experiences_per_rule(r), covered_experience_weights);
      experiences_per_ruleOutcome(r) = experiences_per_outcome;
    }

  }

  // invoke recomputeDefaultRule to collected uncovered / double covered experiences
  recomputeDefaultRule();
  // done

}


void RuleSetContainer::append(Rule* rule, uintA& experiences_of_this_rule, MT::Array< uintA >& experiences_per_outcome_of_this_rule) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"RuleSetContainer::append [START]"<<endl;}
  if (DEBUG>0) {rule->write(cout);  PRINT(experiences_of_this_rule);  PRINT(experiences_per_outcome_of_this_rule);}
  if (Rule::isDefaultRule(rule)  &&   rules.num() > 0) {
    HALT("don't append default rule");
  }
  // (1) update rules
  rules.append(rule);

  if (!Rule::isDefaultRule(rule)) {
    // (2 - A) update experiences_per_rule
    experiences_per_rule.append(experiences_of_this_rule);
    // (3) update nonDefaultRules_per_experience
    uint i;
    FOR1D(experiences_of_this_rule, i) {
      nonDefaultRules_per_experience(experiences_of_this_rule(i)).append(rules.num()-1);
    }
    // (4) update experiences_per_ruleOutcome
    experiences_per_ruleOutcome.append(experiences_per_outcome_of_this_rule);
  }
  else {
    // (2 - B) update experiences_per_rule
    uintA empty;
    experiences_per_rule.append(empty);
    // (4 - B) update experiences_per_ruleOutcome
    MT::Array< uintA > empty_outcome;
    experiences_per_ruleOutcome.append(empty_outcome);
  }
  if (DEBUG>0) {PRINT(nonDefaultRules_per_experience);  PRINT(experiences_per_rule);  PRINT(experiences_per_ruleOutcome);}
  if (DEBUG>0) {cout<<"RuleSetContainer::append [END]"<<endl;}
}


void RuleSetContainer::remove(uint id) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"RuleSetContainer::remove [START]"<<endl;}
  if (DEBUG>0) {PRINT(id);  cout<<"RULES BEFORE REMOVE:"<<endl;  write(cout, true);}
  
  // (1) rules
  rules.remove(id);
  // (2) nonDefaultRules_per_experience: account for new ids of rules behind rule #id
  uint i, k;
  FOR1D(experiences_per_rule(id), i) {
    nonDefaultRules_per_experience(experiences_per_rule(id)(i)).removeValue(id);
  }
  FOR1D(nonDefaultRules_per_experience, i) {
    FOR1D(nonDefaultRules_per_experience(i), k) {
      if (nonDefaultRules_per_experience(i)(k) > id)
        nonDefaultRules_per_experience(i)(k)--;
    }
  }
  // (3) experiences_per_rule  and  experiences_per_ruleOutcome
  MT::Array< uintA > experiences_per_rule__new(rules.num()); // rules have already new size
  MT::Array< MT::Array < uintA > > experiences_per_ruleOutcome__new(rules.num());
  FOR1D_(rules, i) {
    if (i<id) {
      experiences_per_rule__new(i) = experiences_per_rule(i);
      experiences_per_ruleOutcome__new(i) = experiences_per_ruleOutcome(i);
    }
    else {
      experiences_per_rule__new(i) = experiences_per_rule(i+1);
      experiences_per_ruleOutcome__new(i) = experiences_per_ruleOutcome(i+1);
    }
  }
  experiences_per_rule = experiences_per_rule__new;
  experiences_per_ruleOutcome = experiences_per_ruleOutcome__new;
  
  if (DEBUG>0) {cout<<"RULES AFTER REMOVE:"<<endl;  write(cout, true);}
  //   sanityCheck();
  if (DEBUG>0) {cout<<"RuleSetContainer::remove [END]"<<endl;}
}


void RuleSetContainer::clear() {
  rules.clear();
  uint i, k;
  FOR1D(nonDefaultRules_per_experience, i) {
    nonDefaultRules_per_experience(i).clear();
  }
  FOR1D(experiences_per_rule, i) {
    experiences_per_rule(i).clear();
  }
  FOR2D(experiences_per_ruleOutcome, i, k) {
    experiences_per_ruleOutcome(i)(k).clear();
  }
  FOR1D(experiences_per_ruleOutcome, i) {
    experiences_per_ruleOutcome(i).clear();
  }
  experiences_per_rule.clear();
  experiences_per_ruleOutcome.clear();
}


void RuleSetContainer::recomputeDefaultRule() {
  uint DEBUG = 0;

  if (DEBUG>0) cout<<"recomputeDefault [START]"<<endl;
  CHECK(Rule::isDefaultRule(rules.elem(0)), "First rule should be default rule, digger");
  // prepare data-fields
  if (experiences_per_rule.N == 0) {
    uintA empty;
    experiences_per_rule.append(empty);
    MT::Array< uintA > empty_outcome(2);
    experiences_per_ruleOutcome.append(empty_outcome);
  }
  experiences_per_rule(0).clear();
  experiences_per_ruleOutcome(0).clear();
  MT::Array< uintA > empty_outcome(2);
  experiences_per_ruleOutcome(0) = empty_outcome;
  // experiences_per_rule
  uint i;
  FOR1D((*p_experiences), i) {
    if (nonDefaultRules_per_experience(i).N != 1) {
      experiences_per_rule(0).append(i);
      if (((*p_experiences)(i))->noChange())
        experiences_per_ruleOutcome(0)(0).append(i);
      else
        experiences_per_ruleOutcome(0)(1).append(i);
    }
  }
  if (DEBUG>1) {PRINT(experiences_per_rule(0));}
  // finalize rule
  double DEFAULT_CHANGE_PROB = 0.5;
  if ((*p_experiences).N == 0   ||   experiences_per_rule(0).N == 0) {
    rules.overwrite(0, Rule::generateDefaultRule(DEFAULT_CHANGE_PROB));
  }
  else {
    uint changes = 0;
    double no_experiences_per_rule_0 = 0;
    FOR1D(experiences_per_rule(0), i) {
      no_experiences_per_rule_0 += experience_weights(experiences_per_rule(0)(i));
      if ((*p_experiences)(experiences_per_rule(0)(i))->pre != (*p_experiences)(experiences_per_rule(0)(i))->post)
        //changes++;
        changes += experience_weights(experiences_per_rule(0)(i));
    }
    if (no_experiences_per_rule_0 == 0.) {
      rules.overwrite(0, Rule::generateDefaultRule(DEFAULT_CHANGE_PROB));
    } else {
      //double prob_change = ((double) changes) / experiences_per_rule(0).N;
      double prob_change = ((double) changes) / no_experiences_per_rule_0;
      if (DEBUG>1) {PRINT(prob_change);}
      Rule* new_default_rule = Rule::generateDefaultRule(prob_change, 0.05);
      rules.overwrite(0, new_default_rule);
    }
  }

  if (DEBUG>1) {rules.elem(0)->write(cout);}
  if (DEBUG>0) {cout<<"recomputeDefault [END]"<<endl;}
}


void RuleSetContainer::sort() {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"sort [START]"<<endl;}
  RuleSet new__rules;
  MT::Array< uintA > new__nonDefaultRules_per_experience;
  MT::Array< uintA > new__experiences_per_rule;
  MT::Array< MT::Array < uintA > > new__experiences_per_ruleOutcome;
  
  SymL sym_actions;
  uint i, k;
  FOR1D_(rules, i) {
    sym_actions.setAppend(rules.elem(i)->action->s);
  }
  Symbol::sort(sym_actions);
  
  Substitution sub;
  FOR1D(sym_actions, i) {
    if (DEBUG>0) {cout<<"Dealing with action (i="<<i<<")="<<*sym_actions(i)<<endl;}
    uintA rules_ids_for_this_action;
    uintA num_experiences;
    FOR1D_(rules, k) {
      if (rules.elem(k)->action->s == sym_actions(i)) {
        rules_ids_for_this_action.append(k);
        num_experiences.append(experiences_per_rule(k).N);
        // hack for better sorting [start]
        num_experiences.last() = num_experiences.last() * 1000;
        if (experiences_per_rule(k).N > 0)  num_experiences.last() += experiences_per_rule(k)(0);
        // hack for better sorting [end]
      }
    }
    uintA sorted_indices;
    TL::sort_desc_keys(sorted_indices, num_experiences);
    //     PRINT(sorted_indices);
    FOR1D(sorted_indices, k) {
      new__rules.append(rules.elem(rules_ids_for_this_action(sorted_indices(k))));
      sub.addSubs(rules_ids_for_this_action(sorted_indices(k)), new__rules.num()-1);
      new__experiences_per_rule.append(experiences_per_rule(rules_ids_for_this_action(sorted_indices(k))));
      new__experiences_per_ruleOutcome.append(experiences_per_ruleOutcome(rules_ids_for_this_action(sorted_indices(k))));
    }
  }
  new__nonDefaultRules_per_experience = nonDefaultRules_per_experience;
  FOR1D(new__nonDefaultRules_per_experience, i) {
    FOR1D(new__nonDefaultRules_per_experience(i), k) {
      new__nonDefaultRules_per_experience(i)(k) = sub.getSubs(new__nonDefaultRules_per_experience(i)(k));
    }
  }
  
  rules = new__rules;
  experiences_per_rule = new__experiences_per_rule;
  experiences_per_ruleOutcome = new__experiences_per_ruleOutcome;
  nonDefaultRules_per_experience = new__nonDefaultRules_per_experience;
  
  //   sanityCheck();
  if (DEBUG>0) {cout<<"sort [END]"<<endl;}
}


// takes into account that 2 or more non-default rules may cover an experience
void RuleSetContainer::getResponsibilities(arr& responsibilities, MT::Array< uintA >& covered_experiences, uintA& covered_experiences_num) const {
  responsibilities.clear();
  covered_experiences.clear();
  covered_experiences_num.clear();
  responsibilities.resize(rules.num());
  covered_experiences.resize(rules.num());
  covered_experiences_num.resize(rules.num());
  covered_experiences_num.setUni(0);
  uint i;
  FOR1D((*p_experiences), i) {
    if (nonDefaultRules_per_experience(i).N == 1) {
      covered_experiences(nonDefaultRules_per_experience(i)(0)).append(i);
      covered_experiences_num(nonDefaultRules_per_experience(i)(0))++;
    }
    else {
      covered_experiences(0).append(i);
      covered_experiences_num(0)++;
    }
  }
  FOR1D(responsibilities, i) {
    responsibilities(i) = (1.0 * covered_experiences_num(i)) / p_experiences->N;
  }
}


void RuleSetContainer::getPartitionsForAction(MT::Array< uintA >& partitions, Literal* action) const {
  partitions.clear();
  uint i;
  FOR1D_(rules, i) {
    if (rules.elem(i)->action == action) {
      partitions.append(experiences_per_rule(i));
    }
  }
}


// only for visualisation...
void rule_write__specialVisualisation(Rule* rule, MT::Array< uintA >& outcome_tripletts, bool with_action, ostream& os) {
  CHECK(outcome_tripletts.N = rule->outcomes.N, "wrong size");
  //  os << "r" << endl;
  uint i, k;
  // Default rule does not have an action specified...
  if (with_action) {
    os << "ACTION: ";
    if (rule->action != NULL)
      rule->action->write(os, true);
    else
      os << "no_action";
    os << endl;
  }
  os << "CONTEXT: ";
  FOR1D(rule->context, i) {
    os << *rule->context(i) << " ";
  }
  os << endl;
  os << "OUT:" << endl;
  uint total_num_experiences = 0;
  FOR1D(rule->outcomes, i) {total_num_experiences += outcome_tripletts(i).N;}
  FOR1D(rule->outcomes, i) {
    os.precision(2);
    os << "  " << rule->probs(i) << " ";
    FOR1D(rule->outcomes(i), k) {
      rule->outcomes(i)(k)->write(os);
      //       os<<outcomes(i)(k);
      os << " ";
    }
    if (i==rule->outcomes.N-1)
      //       os<<rule->noise_changes;
      os<<"<noise>";
    os<<"    # [";
    FOR1D(outcome_tripletts(i), k) {
      if (k > 10) {
        os<<"...";
        break;
      }
      os<<outcome_tripletts(i)(k)<<" ";
    }
    os <<"]";
    os << " (" << outcome_tripletts(i).N << "/" << total_num_experiences << " = ";
    if (outcome_tripletts(i).N == total_num_experiences) os<<"100";
    else os << ((uint) 100 * (outcome_tripletts(i).N * 1.0 / total_num_experiences));
    os << "%)" << endl;
  }
  if (Rule::isDefaultRule(rule) &&  outcome_tripletts(0).N > 0) {os<<"ACHTUNG!!! Noisy default rule used to model!!"<<endl;  /*MT_MSG("ACHTUNG!!! Noisy default rule used to model!!");*/}
  //   if (outcome_tripletts(rule->outcomes.N-1).N > 0) {os<<"ACHTUNG!!! Using noise-outcome!"<<endl;  MT_MSG("ACHTUNG!!! Using noise-outcome!");}
}


void RuleSetContainer::write(ostream& out, bool only_action, bool additional_experience_info) const {
  uint i, k;
  out<<"# *** RULES ***"<<endl;
  Literal* last_action = NULL;
  FOR1D_(rules, i) {
    if (last_action != rules.elem(i)->action) {
      last_action = rules.elem(i)->action;
      out<<"##### ";  rules.elem(i)->action->write(out, true);  out<<endl;
    }
    out<<"# "<<i<<"    ("<<(100. * experiences_per_rule(i).N) / nonDefaultRules_per_experience.N
      <<"%, "<<experiences_per_rule(i).N<<"/"<<nonDefaultRules_per_experience.N<<")";
    //     " experiences [";
    //     FOR1D(experiences_per_rule(i), k) {
    //       out<<experiences_per_rule(i)(k)<<" ";
    //       if (k > 10) {
    //         out<<"...";
    //         break;
    //       }
    //     }
    //     out << "]";
    out<<endl;
    if (only_action) {
      rules.elem(i)->action->write(out); out<<endl;
    }
    else {
      rule_write__specialVisualisation(rules.elem(i), experiences_per_ruleOutcome(i), false, out);
    }
  }
  
  if (additional_experience_info && rules.num() > 0) {
    uint num_explanations = p_experiences->N;
    uint num_explanations_by_default_rule = experiences_per_rule(0).N;
    uint num_explanations_by_nonDefault_rule = num_explanations - num_explanations_by_default_rule;
    uint num_explanations_as_noise = 0;
    FOR1D(experiences_per_ruleOutcome, k) {
      num_explanations_as_noise += experiences_per_ruleOutcome(k).last().N;
    }
    uint num_explanations_as_nonNoise = num_explanations - num_explanations_as_noise;
    out<<endl<<"#  *** "<<nonDefaultRules_per_experience.N<<" EXPERIENCES (experience:rule-id) ***"<<endl;
    out<<"#  NON-default rules explain: "<<(100. * num_explanations_by_nonDefault_rule) / num_explanations
      <<"% ("<<num_explanations_by_nonDefault_rule<<"/"<<num_explanations<<")"<<endl;
    out<<"#  NON-noise outcomes explain: "<<(100. * num_explanations_as_nonNoise) / num_explanations
      <<"% ("<<num_explanations_as_nonNoise<<"/"<<num_explanations<<")"<<endl;
    FOR1D(nonDefaultRules_per_experience, i) {
      out<<"# "<<i<<":";
      if (nonDefaultRules_per_experience(i).N > 1)
        out<<nonDefaultRules_per_experience(i)<<"  ";
      else if (nonDefaultRules_per_experience(i).N == 1)
        out<<nonDefaultRules_per_experience(i)(0)<<"  ";
      else
        out<<"default  ";
      out << " (weight=" << experience_weights(i) << ")" << endl;
    }
    out<<endl;
  }
}


void RuleSetContainer::write(const char* filename, bool only_action, bool additional_experience_info) const {
  ofstream out(filename);
  write(out, only_action, additional_experience_info);
}


void RuleSetContainer::write_experiencesWithRules(ostream& os) const {
  uint i, k;
  FOR1D((*p_experiences), i) {
    os<<"--------------"<<endl;
    os<<"EXPERIENCE #"<<i<<":"<<endl;
    ((*p_experiences)(i))->write(os);
    os<<"  ---> COVERING RULES: "<<nonDefaultRules_per_experience(i)<<endl;
    FOR1D(nonDefaultRules_per_experience(i), k) {
      rules.elem(nonDefaultRules_per_experience(i)(k))->write(os);
    }
  }
}


void RuleSetContainer::write_rulesWithExperiences(ostream& os) const {
  uint i, k;
  FOR1D_(rules, i) {
    os<<"--------------"<<endl;
    os<<"RULE #"<<i<<":"<<endl;
    rules.elem(i)->write(os);
    os<<"COVERS "<<experiences_per_rule(i).N<<" experiences."<<endl;
    FOR1D(experiences_per_rule(i), k) {
      os<<"Experience #"<<experiences_per_rule(i)(k)<<endl;
      ((*p_experiences)(experiences_per_rule(i)(k)))->write(os);
    }
  }
}


void RuleSetContainer::sanityCheck(bool ignore_default_rule) const {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"RuleSetContainer::sanityCheck [START]"<<endl;}
  rules.sanityCheck();
  if (rules.num() != experiences_per_rule.N) {
    cout<<"FAILED SANITY CHECK:"<<endl;
    write(cout, true);
    HALT("sanity check failed 0-A");
  }
  if (p_experiences->N != nonDefaultRules_per_experience.N) {
    cout<<"FAILED SANITY CHECK:"<<endl;
    write(cout, true);
    HALT("sanity check failed 0-B");
  }
  if (rules.num() != experiences_per_ruleOutcome.N) {
    cout<<"FAILED SANITY CHECK:"<<endl;
    write(cout, true);
    HALT("sanity check failed 0-C");
  }
  
  uint i, k;

  // Non-default rules
  FOR1D_(rules, i) {
    if (i == 0)
      continue;
    FOR1D((*p_experiences), k) {
      SubstitutionSet subs;
      bool supposed_to_cover__1 = experiences_per_rule(i).findValue(k) >= 0;
      bool supposed_to_cover__2 = nonDefaultRules_per_experience(k).findValue(i) >= 0;
      if (supposed_to_cover__1 != supposed_to_cover__2) {
        cout<<"FAILED SANITY CHECK:"<<endl;
        write(cout, true);
        PRINT(i);
        PRINT(k);
        PRINT(experiences_per_rule(i));
        PRINT(nonDefaultRules_per_experience(k));
        PRINT(supposed_to_cover__1);
        PRINT(supposed_to_cover__2);
        HALT("sanity check failed 1:  supposed_to_cover__1 != supposed_to_cover__2");
      }
      bool covers = reason::calcSubstitutions_rule_groundAction(subs, (*p_experiences)(k)->pre, (*p_experiences)(k)->action, rules.elem(i));
      if (covers != supposed_to_cover__1) {
        cout<<"FAILED SANITY CHECK 2:  covers != supposed_to_cover__1 "<<endl;
        if (covers) {
          cout<<"Rule covers experience although it is supposed to NOT cover experience according to rulesC-information."<<endl;
        }
        else {
          cout<<"Rule does NOT cover experience although it is supposed to cover experience according to rulesC-information."<<endl;
        }
        write(cout, true);
        cout<<"rule #i="<<i<<endl;
        cout<<"experience #k="<<k<<endl;
        PRINT(supposed_to_cover__1);
        PRINT(covers);
        cout<<"Rule:"<<endl;  rules.elem(i)->write();  cout<<endl;
        cout<<"StateTransition:"<<endl;  (*p_experiences)(k)->write(cout,2);  cout<<endl;
        HALT("sanity check failed 2");
      }
    }
    uint total_experiences_in_outcomes = 0;
    FOR1D(experiences_per_ruleOutcome(i), k) {
      total_experiences_in_outcomes += experiences_per_ruleOutcome(i)(k).N;
      uint l;
      FOR1D(experiences_per_ruleOutcome(i)(k), l) {
        if (experiences_per_rule(i).findValue(experiences_per_ruleOutcome(i)(k)(l)) < 0) {
          cout<<"FAILED SANITY CHECK:"<<endl;
          write(cout, false);
          PRINT(i);
          PRINT(k);
          PRINT(l);
          PRINT(experiences_per_ruleOutcome(i)(k)(l));
          PRINT(experiences_per_ruleOutcome(i)(k));
          PRINT(experiences_per_ruleOutcome(i));
          PRINT(experiences_per_rule(i));
          HALT("sanity check failed 3");
        }
      }
    }
    if (total_experiences_in_outcomes < experiences_per_rule(i).N) {
      cout<<"FAILED SANITY CHECK:"<<endl;
      write(cout, false);
      PRINT(i);
      PRINT(total_experiences_in_outcomes);
      PRINT(experiences_per_ruleOutcome(i));
      PRINT(experiences_per_rule(i));
      HALT("sanity check failed 4");
    }
  }
  
  if (!ignore_default_rule) {
    // Check that default rule does not cover experiences which other rules cover
    FOR1D(experiences_per_rule(0), i) {
      FOR1D_(rules, k) {
        if (k == 0)
          continue;
        if (experiences_per_rule(k).findValue(experiences_per_rule(0)(i)) >= 0) {
          cout<<"FAILED SANITY CHECK:  default rule covers experience which also other rules covers"<<endl;
          write(cout, false);
          PRINT(k);
          PRINT(experiences_per_rule(0));
          PRINT(experiences_per_rule(k));
          PRINT(i);
          PRINT(experiences_per_rule(0)(i));
          HALT("sanity check failed 5  --  default rule A");
        }
      }
    }
    
    if (Rule::isDefaultRule(rules.elem(0))) {
      FOR1D(nonDefaultRules_per_experience, i) {
        if (nonDefaultRules_per_experience(i).findValue(0) >= 0) {
          cout<<"FAILED SANITY CHECK:"<<endl;
          write(cout, false);
          PRINT(i);
          PRINT(nonDefaultRules_per_experience);
          HALT("sanity check failed 6  --  default rule B");
        }
      }
    }
  }
  if (DEBUG>0) {cout<<"RuleSetContainer::sanityCheck [END]"<<endl;}
}







/************************************************
 *
 *     learn -- details
 *
 *     (only of interest for programmers)
 *
 ************************************************/

namespace learn {
// for parameter learning
double __pen_sum__base = 20.;
double __pen_sum__final = 20.;
double __pen_pos__base = 20.;
double __pen_pos__final = 20.;
}

namespace learn { namespace CostFunction {

StateTransitionL __cf_experiences_coveredByCurrentRule;
boolA __cf_coverage_outcome_experience;

} }

void learn::calcCoverage(StateTransitionL& covered_experiences, uintA& covered_experiences_ids, const Rule* r, const StateTransitionL& experiences, bool consider_non_unique_covering) {
  arr covered_experience_weights;
  arr experience_weights(experiences.N);
  experience_weights.setUni(1);
  learn::calcCoverage(covered_experiences, covered_experiences_ids, covered_experience_weights, r, experiences, experience_weights, consider_non_unique_covering);
}

void learn::calcCoverage(StateTransitionL& covered_experiences, uintA& covered_experiences_ids, arr& covered_experience_weights, const Rule* r, const StateTransitionL& experiences, const arr& experience_weights, bool consider_non_unique_covering) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"learn::calcCoverage [START]"<<endl;
  if (DEBUG>0) r->write(cout);
  covered_experiences.clear();
  covered_experiences_ids.clear();
  covered_experience_weights.clear();
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
      if (!reason::calcSubstitutions_rule_groundAction(subsNB, experiences(i)->pre, experiences(i)->action, &ruleWithoutNonBinaries, consider_non_unique_covering)) {
        if (DEBUG>0) cout<<" --> 0  (NO_DEICTICREFS_BY_NONBINARY)"<<endl;
        continue;
      }
    }
#endif
    SubstitutionSet subs;
    if (reason::calcSubstitutions_rule_groundAction(subs, experiences(i)->pre, experiences(i)->action, r, consider_non_unique_covering)) {
      covered_experiences.append(experiences(i));
      covered_experiences_ids.append(i);
      covered_experience_weights.append(experience_weights(i));
      if (DEBUG>0) cout<<" --> 1"<<endl;
    }
    else {
      if (DEBUG>0) cout<<" --> 0"<<endl;
    }
  }
  if (DEBUG>0) cout<<"learn::calcCoverage [END]"<<endl;
}


void learn_outcomes__calcSubsumption(boolA& subsumes, const boolA& coverage) {
  CHECK(coverage.nd == 2, "invalid coverage matrix")
      subsumes.resize(coverage.d0, coverage.d0);
  subsumes.setUni(false);
  uint i, j, k;
  FOR1D(subsumes, i) {
    subsumes(i,i) = true;
    for (j=i+1; j<subsumes.d0; j++) {
      subsumes(i,j) = true;
      subsumes(j,i) = true;
      for (k=0; k<coverage.d1; k++) {
        if (coverage(i, k) < coverage(j, k))
          subsumes(i,j) = false;
        if (coverage(j, k) < coverage(i, k))
          subsumes(j,i) = false;
      }
      if (!subsumes(i,j) && !subsumes(j,i))
        break;
    }
  }
}



/************************************************
 *
 *     learn -- details -- outcomes
 *
 ************************************************/



void __calcCoverage_outcomes(boolA& coverage, const MT::Array< LitL >& potential_outcomes, const StateTransitionL& covered_experiences, const Rule* rule) {
  uint DEBUG = 0;
  if (DEBUG>0) cout << "calcOutcomesCoverage [START]" << endl;
  Rule rule_full;
  rule_full.copyBody(*rule);
  rule_full.outcomes = potential_outcomes;
  rule_full.probs.resize(rule_full.outcomes.N);
  if (DEBUG>0) {cout<<rule_full<<endl;}
  coverage.resize(potential_outcomes.N, covered_experiences.N);
  coverage.setUni(false);
  uint i, e;
  FOR1D(covered_experiences, e) {
    if (DEBUG>1) {cout<<"### Ex "<<e<<":"<<endl; covered_experiences(e)->write(cout);}
    uintA covering_outcomes;
    reason::calc_coveringOutcomes(covering_outcomes, &rule_full, covered_experiences(e)->pre, covered_experiences(e)->action, covered_experiences(e)->post);
    FOR1D(covering_outcomes, i) {
      coverage(covering_outcomes(i), e) = true;
    }
    if (covering_outcomes.N > 1) {  // ignore noise outcome if there has been other covering outcome
      coverage(potential_outcomes.N-1, e) = false;
    }
    if (DEBUG>1) {PRINT(covering_outcomes);}
  }
  if (DEBUG>0) {PRINT(coverage);}
  if (DEBUG>0) cout << "calcOutcomesCoverage [END]" << endl;
}


// remove outcomes that (i) do not cover any experience and (ii) have zero-probability  and (iii) sets coverage for cost function
void __get_trimmed_outcomes(MT::Array< LitL >& outcomes, arr& probs, boolA& coverage, const StateTransitionL& coveredExperiences, arr covered_experience_weights, const Rule& rule, uint DEBUG=0) {
  if (DEBUG>0) cout<<"get_trimmed_outcomes [START]"<<endl;
  uint i;
  bool removed, atleastone;
  if (DEBUG>1) {cout<<"Outcomes before:\n"; write(outcomes); PRINT(coveredExperiences);}
  
  // (i) remove outcomes that do not cover any experience
  if (DEBUG>1) cout<<"Removing outcomes that do not cover any experience..."<<endl;
  __calcCoverage_outcomes(coverage, outcomes, coveredExperiences, &rule);
  if (DEBUG>1) {PRINT(coverage);}
  learn::CostFunction::setOutcomesCoverage(coverage);
  if (DEBUG>3) {PRINT(coverage);}
  removed = false;
  MT::Array< LitL > o_help;
  FOR1D(outcomes, i) {
    atleastone = sum(coverage.sub(i,i,0,coverage.d1-1));
    if (i < outcomes.N-1  &&  !atleastone) { // only for non-noise outcomes
      removed = true;
      if (DEBUG>2) cout << "Outcome " << i << " covers 0 experiences and will be removed."<<endl;
    }
    else
      o_help.append(outcomes(i));
  }
  if (DEBUG>1) {cout<<"Outcomes now:\n"; write(o_help);}
  if (removed) {
    __calcCoverage_outcomes(coverage, o_help, coveredExperiences, &rule);
    learn::CostFunction::setOutcomesCoverage(coverage);
    if (DEBUG>3) PRINT(coverage)
  }

  // (ii) remove zero-probability outcomes
  if (DEBUG>1) cout<<"Removing zero-prob outcomes"<<endl;
  // learn params


  learn::PasulaScoreFunction* psf = dynamic_cast<learn::PasulaScoreFunction*>(learn::ScoreFunction::getCurrentScoreFunction());
  if (psf->getPmin() == 0.) {
    // COMPUTE PARAMETERS
    // we do not need to learn parameters; we can just calculate them as the frequencies / weights of the experiences
    double weight_sum=0.;
    probs.resize(o_help.N);

    uint o;
    FOR1D(o_help, o) {
      double o_weight_sum=0.;
      for(uint i=0; i < coverage.dim(1); i++)  {
        if (coverage(o,i) == 1) {
          o_weight_sum += covered_experience_weights(i);
        }
      }
      weight_sum += o_weight_sum;
      probs(o) = o_weight_sum;
    }
    probs /= weight_sum;
  }

  // LEARN PARAMETERS
  else {
    learn::learn_parameters(o_help, probs);
  }

  // remove zero-prob outcomes (except default outcome which can never be deleted)
  removed = false;
  outcomes.clear();
  FOR1D(o_help, i) {
    if (i!=o_help.N-1  &&  TL::isZero(probs(i))) {
      probs.remove(i);
      removed = true;
    }
    else {
      outcomes.append(o_help(i));
    }
  }
  
  if (removed) {
    __calcCoverage_outcomes(coverage, outcomes, coveredExperiences, &rule);
    learn::CostFunction::setOutcomesCoverage(coverage);
    if (DEBUG>3) PRINT(coverage)
  }
  if (DEBUG>1) {cout<<"Outcomes final:\n"; write(outcomes);}
  if (DEBUG>0) cout<<"get_trimmed_outcomes [END]"<<endl;
}

namespace learn {
bool __compareUintAscending(const uint& a, const uint& b) {
  return a > b ? 1 : (a < b ? - 1 : 0);
}
}

// Steps:
// (1) Determine basic outcomes: changes from pre to post in experiences
// (2) Collapse identical outcomes
// (3) Greedily improve outcomes based on add- and remove-operator
//
// When we create new outcomes, we have to check whether there are outcomes with
// (i) no experiences covered or (ii) a probability of 0. In both cases, the corresponding
// outcomes need to be deleted. --> method produceTrimmedOutcomes(...), see above
//
// TODO cope with comparisons (a la Zettlemoyer et al., 2004)
void learn::learn_outcomes(Rule* r, MT::Array< uintA >& coveredExperiences_per_outcome_slim, const StateTransitionL& coveredExperiences_slim,const uintA& covered_experiences_ids_slim, const arr& covered_experience_weights_slim, uint DEBUG) {
  //DEBUG = 10;
  if (DEBUG>0) cout << "induceOutcomes [START]" << endl;
  uint i, k;

  CHECK(covered_experiences_ids_slim.N == coveredExperiences_slim.N, "");
  CHECK(covered_experiences_ids_slim.N == covered_experience_weights_slim.N, "");
  CHECK(coveredExperiences_slim.N>0, "No experiences, bro!");

  // begin expansion
  // in order to use experience_weights we are going to assume that weights are
  // integer and we create artificial duplicates of experiences with weight > 1
  MT::Array< uintA > coveredExperiences_per_outcome;
  StateTransitionL coveredExperiences;
  uintA covered_experiences_ids;

  // normalize weights
  // find minimum weights > 0
  arr covered_experience_weights_slim_nrm;

  arr _non_min_weights;
  static const double min_weight = 1e-3;
  FOR1D(covered_experience_weights_slim, i) { if (covered_experience_weights_slim(i) > min_weight) _non_min_weights.append(covered_experience_weights_slim(i)); }

  if (_non_min_weights.N > 0) {
    double min_w = _non_min_weights.min();

    // normalize weights
    FOR1D(covered_experience_weights_slim,i){
        covered_experience_weights_slim_nrm.append(covered_experience_weights_slim(i) / min_w);
    }
  } else {
    covered_experience_weights_slim_nrm.resize(covered_experience_weights_slim.N);
    covered_experience_weights_slim_nrm.setUni(1.);
  }

//  covered_experience_weights_slim_nrm = covered_experience_weights_slim;

  uintA compressed2Expanded;
  FOR1D(coveredExperiences_slim, i) {
    uint id = covered_experiences_ids_slim(i);
    //int w = round(covered_experience_weights_slim(i));
    int w = round(covered_experience_weights_slim_nrm(i));
    if (w > 1000) {
        if (DEBUG>0)
          cerr << "  WARN clipping learn_outcomes weight from " << w << " to 1000" << endl;
        w = 1000;
    }
    for(uint t = 0; t < w; t++) {
      StateTransition* st = coveredExperiences_slim(i);
      coveredExperiences.append(st);
      covered_experiences_ids.append(coveredExperiences.N-1);
      compressed2Expanded.append(id);
    }
  }

  arr covered_experience_weights(coveredExperiences.N);
  covered_experience_weights.setUni(1);
  // end expansion

  CHECK(covered_experiences_ids.N == coveredExperiences.N, "");
  CHECK(coveredExperiences.N>0, "No experiences, bro!");

  // calc changes first per experience
  if (DEBUG > 0) {r->write(cout);  PRINT(coveredExperiences.N);}
  
  // prepare cost function (stays like this for complete following induceOutcomes-procedure)
  CostFunction::setRuleCoveredExperiences(coveredExperiences);
  __pen_sum__final = __pen_sum__base * coveredExperiences.N;
  __pen_pos__final = __pen_pos__base * coveredExperiences.N;

  MT::Array<Literal*> varComparisons;
  FOR1D(r->context, i) {
    if (r->context(i)->comparison_type == Literal::comparison_variable)
      varComparisons.append(r->context(i));
  }

  // (1) Determine basic outcomes = changes from pre to post = for each covered experience a separate outcome
  MT::Array< LitL > outcomes_basic;
  FOR1D(coveredExperiences, i) {
    if (DEBUG>4) {cout<<"====== Using ex "<<i<<":"<<endl; coveredExperiences(i)->write(cout);}
    SubstitutionSet subs;
    bool covered = reason::calcSubstitutions_rule_groundAction(subs, coveredExperiences(i)->pre, coveredExperiences(i)->action, r);
    CHECK(covered, "An uncovered experience! Be careful when calling methods, noob!");
    if (subs.num() != 1) {
      cout << "FAILING: TOO " << (subs.num() > 1 ? "MANY" : "FEW") << "SUBS"<<endl;
      cout<<"SymbolicState: ";coveredExperiences(i)->pre.write(cout);cout<<endl;
      cout<<"Action: ";coveredExperiences(i)->action->write(cout);cout<<endl;
      cout<<"Rule: "<<endl;r->write(cout);cout<<endl;
      cout << "Substitutions: ";
      uint s;
      FOR1D_(subs, s) {
        subs.elem(s)->write(cout); cout << endl;
      }
    }
    CHECK(subs.num()==1, "Cannot be deictic variant! (Should've been taken care for somewhere else.)")

        LitL nextOutcome;
    // just take first subs --> DEICTIC
    uint subsId = 0;
    
    Substitution invSub;
    subs.elem(subsId)->getInverse(invSub);
    if (DEBUG > 4) {
      cout<<"Substitution: "; subs.elem(subsId)->write(cout); cout << endl;
      cout<<"Inverse Substitution: "; invSub.write(cout);cout << endl;
    }

    //compute substitutions for function value variables seperate, not very elegant but with the existing data structures the easiest solution
    std::multimap<Symbol*, Literal*> functionValueSubs;
    FOR1D(varComparisons, k) {
      uint l;
      FOR1D(coveredExperiences(i)->pre.lits, l) {
        Literal *litInv = invSub.apply(coveredExperiences(i)->pre.lits(l));
        if (varComparisons(k)->s == litInv->s && varComparisons(k)->args == litInv->args)
          functionValueSubs.insert(std::make_pair(varComparisons(k)->s, litInv));
      }
    }

    // insert changed literals
    FOR1D(coveredExperiences(i)->changes, k) {
      if (coveredExperiences(i)->changes(k)->s->symbol_type == Symbol::primitive) {
        Literal *litInv = invSub.apply(coveredExperiences(i)->changes(k));

        std::multimap<Symbol*, Literal*>::iterator find = functionValueSubs.find(coveredExperiences(i)->changes(k)->s);
        for (; find != functionValueSubs.end(); find++) {
          if (find->second->args == litInv->args) {
            //create Literal X + offset
            double offset = litInv->value - find->second->value;
            nextOutcome.append(Literal::get(Symbol::get(litInv->s->name, litInv->s->arity, Symbol::primitive, litInv->s->range_type), litInv->args, offset, Literal::comparison_offset));
            break;
          }
        }

        //default
        if (find == functionValueSubs.end())
          nextOutcome.setAppend(litInv);
      }
    }

    LitL nextOutcome_pureAbstract;
    FOR1D(nextOutcome, k) {
      if (reason::isPurelyAbstract(nextOutcome(k)))
        nextOutcome_pureAbstract.append(nextOutcome(k));

    }
    outcomes_basic.append(nextOutcome_pureAbstract);

    if (DEBUG>3) {
      cout<<"nextOutcome: "; write(nextOutcome); cout<<endl;
      cout<<"nextOutcome_pureAbstract: "; write(nextOutcome_pureAbstract); cout<<endl;
    }
  }
  // add noise outcomes
  LitL noiseOutcome;
  outcomes_basic.append(noiseOutcome);
  
  if (DEBUG > 0) {
    cout << "StateTransition outcomes (incl. noise outcome):" << endl;
    FOR1D(outcomes_basic, i) {
      cout << "(" << i << ") ";
      write(outcomes_basic(i));
      cout << endl;
    }
  }
  
  
  // (2) Collapse identical outcomes
  // ignore noisy outcome
  boolA prune(outcomes_basic.d0);
  prune.setUni(false);
  FOR1D(outcomes_basic, i) {
    if (i==outcomes_basic.N-1)
      break;
    if (prune(i))
      continue;
    for (k=i+1; k<outcomes_basic.d0-1; k++) {
      if (prune(k))
        continue;
      if (Literal::equivalent(outcomes_basic(i), outcomes_basic(k)))
        prune(k) = true;
    }
  }
  MT::Array< LitL > outcomes;
  FOR1D(outcomes_basic, i) {
    if (!prune(i)) {
      outcomes.append(outcomes_basic(i));
    }
  }
  if (DEBUG > 0) {
    cout << "Collapsed experiences outcomes (copies removed)  (incl. noise outcome):" << endl;
    FOR1D(outcomes, i) {cout << "(" << i << ") "; write(outcomes(i)); cout << endl;}
  }

  // trim outcomes
  arr probs;
  boolA coverage;
  __get_trimmed_outcomes(outcomes, probs, coverage, coveredExperiences, covered_experience_weights, *r);
  
  // (3) Greedily improve outcomes
  // score needs to be optimized

  double score, bestScore;
  double loglik;
  
  // evaluate and calc score
  loglik = CostFunction::loglikelihood(probs);
  uint num_outcome_literals = 0;
  FOR1D(outcomes, i) {num_outcome_literals += outcomes(i).N;}
  score = loglik - __alpha_PEN * num_outcome_literals;
  if (DEBUG > 1) {
    cout << "\nBASIC OUTCOMES:" << endl;
    FOR1D(outcomes, i) {
      cout << "("<< i << ") ";
      write(outcomes(i));
      cout << " " << probs(i);
      cout << endl;
    }
    cout << " --> score=" << score << endl;
  }

  bestScore = score;
  bool perform_add = true;  // war urspruenglich false!!
  bool add_possible = true;
  bool remove_possible = true;
  bool change = false;
  boolA subsumes;
  arr probs_new;
  MT::Array< LitL > outcomes_new;
  boolA coverage_new;
  do {
    change = false;
    outcomes_new.clear();
    // ADD [start]
    if (perform_add) {
      if (DEBUG>4) cout<<"Outcome adding"<<endl;
      boolA unifiable(outcomes.N-1, outcomes.N-1);
      uint numUnifiables = 0;
      unifiable.setUni(false);
      // calculate which ones are unifiable
      FOR1D(unifiable, i) {
        for(k=i+1; k<unifiable.d1; k++) {
          unifiable(i,k) = Literal::nonContradicting(outcomes(i), outcomes(k));
          if (unifiable(i,k))
            numUnifiables++;
        }
      }
      // randomly choose two and unify them
      if (numUnifiables>0) {
        change = true;
        int ix = rnd.num(numUnifiables);
        bool found = false;
        for(i=0; i<unifiable.d0; i++) {
          for(k=i+1; k<unifiable.d1; k++) {
            if (unifiable(i,k)) {
              if (ix-- == 0) {
                found = true; // we'll combining outcomes i and k
                break;
              }
            }
          }
          if (found)
            break;
        }
        LitL unifiedOutcome;
        unifiedOutcome.memMove = true;
        unifiedOutcome.append(outcomes(i));
        unifiedOutcome.append(outcomes(k));
        uint o, o2;
        FOR1D(unifiedOutcome, o) {
          uintA indices;
          unifiedOutcome.findValues(indices, unifiedOutcome(o));
          for (o2=indices.N-1; o2>0; o2--) {
            unifiedOutcome.remove(indices(o2));
          }
        }
        outcomes_new.append(unifiedOutcome);
        FOR1D(outcomes, o) {
          if (o==i || o==k)
            continue;
          else
            outcomes_new.append(outcomes(o));
        }
        if (DEBUG>4) {
          PRINT(numUnifiables)
              PRINT(unifiable)
              cout<<"We gonna add:"<<endl;
          cout<<i<<" :";write(outcomes(i));cout<<endl;
          cout<<k<<" :";write(outcomes(k));cout<<endl;
          cout<<"  --> ";write(unifiedOutcome);cout<<endl;
        }
      }
    }
    // ADD [end]
    // REMOVE [start]
    else {
      if (DEBUG>4) cout<<"Outcome deleting"<<endl;
      // remove condition: outcome overlapping with other outcomes on every experience
      // calc which outcomes overlap each other completely
      learn_outcomes__calcSubsumption(subsumes, coverage);
      // Randomly determine overlapped outcome that is to delete.
      // As usual, ignore default outcome.
      boolA overlapped(outcomes.N-1);
      overlapped.setUni(false);
      for (i=0; i<outcomes.N-1; i++) {
        for (k=0; k<outcomes.N-1; k++) {
          if (k==i)
            continue;
          if (subsumes(k,i)) {
            overlapped(i) = true;
            break;
          }
        }
      }
      uint numOverlapped = 0;
      FOR1D(overlapped,i)
          if (overlapped(i))
          numOverlapped++;
      if (DEBUG > 2) {
        PRINT(coverage_new)
            PRINT(subsumes)
            PRINT(overlapped)
            PRINT(numOverlapped)
      }
      if (numOverlapped > 0) {
        uint overlapped2delete = rnd.num(numOverlapped);
        if (DEBUG>1) {cout<<"Delete overlapped outcome #" << overlapped2delete << endl;}
        numOverlapped = 0;
        FOR1D(outcomes, i) {
          // always keep default outcome:
          if (i==outcomes.N-1) {
            outcomes_new.append(outcomes(i));
            break;
          }
          // other outcomes:
          else if (overlapped(i)) {
            if (numOverlapped != overlapped2delete) {
              outcomes_new.append(outcomes(i));
            }
            else {
              CHECK(!change, "Deleting twice is not allowed.")
                  // HERE IS THE ESSENTIAL LINE WHERE WE DELETE!
                  change = true;
            }
            numOverlapped++;
          }
          else
            outcomes_new.append(outcomes(i));
        }
      }
    }
    // REMOVE [end]

    // POSTPROCESSING OF NEW OUTCOME

    // If no new outcome found:
    // (1) Stop, if no more modifications possible (i.e. neither remove nor add)
    // (2) Change to remove or add
    if (!change) {
      bool stop = false;
      if (perform_add) {
        if (!remove_possible)
          stop = true;
        else {
          add_possible = false;
          perform_add = false;
          continue;
        }
      }
      else {
        if (!add_possible)
          stop = true;
        else {
          remove_possible = false;
          perform_add = true;
          continue;
        }
      }
      if (stop) {
        break;
      }
    }
    // Outcomes _have_ changed.
    else {
      if (DEBUG > 2) {
        cout << "\nNEW OUTCOMES:" << endl;
        FOR1D(outcomes_new, i) {
          cout << i << ": ";
          write(outcomes_new(i));
          cout << endl;
        }
      }

      // procduce trimmed outcomes
      __get_trimmed_outcomes(outcomes_new, probs_new, coverage_new, coveredExperiences, covered_experience_weights, *r);

      // evaluate and calc score
      loglik = CostFunction::loglikelihood(probs_new);
      num_outcome_literals = 0;
      FOR1D(outcomes, i) {num_outcome_literals += outcomes(i).N;}
      score = loglik - __alpha_PEN * num_outcome_literals;

      if (DEBUG > 1) {
        cout << "\nNEW OUTCOMES (now with new_probs and score):" << endl;
        FOR1D(outcomes_new, i) {
          cout << i << ": ";
          write(outcomes_new(i));
          cout << " " << probs_new(i);
          cout << endl;
        }
        cout << " --> score=" << score << endl;
      }

      if (score > bestScore  /*||   (!perform_add  &&  score > bestScore)*/) {
        if (DEBUG>1) cout<<"New outcomes accepted."<<endl;
        outcomes = outcomes_new;
        coverage = coverage_new;
        probs = probs_new;
        bestScore = score;
        if (perform_add)
          remove_possible = true; // ensures that we will try the remove-operator again
        else
          add_possible = true;
      }
      else {
        if (DEBUG>1) cout<<"New outcomes rejected."<<endl;
        if (perform_add)
          add_possible = false;
        else
          remove_possible = false;
        if (!add_possible && !remove_possible)
          break;
        perform_add = !perform_add;
      }
    }
  } while (outcomes.N > 1);
  
  //     CHECK(outcomes.N>1, "No non-noise outcome!")

  // set outcomes in rule
  r->outcomes = outcomes;
  r->probs = probs;
  
  // Calc experiences_per_outcome
  coveredExperiences_per_outcome.clear();
  coveredExperiences_per_outcome.resize(coverage.d0);
  uint i_out, i_ex;
  FOR2D(coverage, i_out, i_ex) {
    if (coverage(i_out, i_ex)) {
      coveredExperiences_per_outcome(i_out).append(covered_experiences_ids(i_ex));
    }
  }
  
  
  // sort for probabilities
  MT::Array< LitL > sorted__outcomes(r->outcomes.N);
  arr sorted__probs(r->outcomes.N);
  MT::Array< uintA > sorted__coveredExperiences_per_outcome(r->outcomes.N);
  uintA sorted_indices;
  TL::sort_desc_keys(sorted_indices, probs);
  uint current_i = 0;
  FOR1D(sorted_indices, i) {
    if (sorted_indices(i) == r->outcomes.N-1)   // skip noise outcome
      continue;
    sorted__outcomes(current_i) = r->outcomes(sorted_indices(i));
    sorted__probs(current_i) = r->probs(sorted_indices(i));
    sorted__coveredExperiences_per_outcome(current_i) = coveredExperiences_per_outcome(sorted_indices(i));
    current_i++;
  }
  CHECK(current_i == r->outcomes.N-1, "");
  sorted__probs.last() = r->probs.last();   // noise outcome
  sorted__coveredExperiences_per_outcome.last() = coveredExperiences_per_outcome.last();   // noise outcome
  r->outcomes = sorted__outcomes;
  r->probs = sorted__probs;
  coveredExperiences_per_outcome = sorted__coveredExperiences_per_outcome;
  
  if (DEBUG > 0) {
    cout << "\nFINAL NEW OUTCOMES (now with probs and score):" << endl;
    PRINT(covered_experiences_ids);
    PRINT(coveredExperiences_per_outcome);
    FOR1D(r->outcomes, i) {
      cout << i << ": ";
      write(r->outcomes(i));
      cout << " " << r->probs(i) << " " << coveredExperiences_per_outcome(i);
      cout << endl;
    }
    cout << " --> score=" << bestScore << endl;
  }
  
  // check [START]
  uint used_outcomes = 0;
  FOR1D(coveredExperiences_per_outcome, i) {
    used_outcomes += coveredExperiences_per_outcome(i).N;
  }
  if (used_outcomes < coveredExperiences.N) {
    PRINT(used_outcomes);
    PRINT(covered_experiences_ids);
    PRINT(coveredExperiences_per_outcome);
    FOR1D(r->outcomes, i) {
      cout << i << ": ";
      write(r->outcomes(i));
      cout << " " << r->probs(i) << " " << coveredExperiences_per_outcome(i);
      cout << endl;
    }
    HALT("error in experiences_per_outcome calculation");
  }
  // check [END]
  
  // now map back to compressed version
  coveredExperiences_per_outcome_slim.clear();
  coveredExperiences_per_outcome_slim.resize(coveredExperiences_per_outcome.N);
  FOR1D(coveredExperiences_per_outcome, i) {
    uint j;
    FOR1D(coveredExperiences_per_outcome(i), j) {
      uint slimId = compressed2Expanded(coveredExperiences_per_outcome(i)(j));
      if (!coveredExperiences_per_outcome_slim.elem(i).contains(slimId))
        coveredExperiences_per_outcome_slim.elem(i).insertInSorted(slimId, learn::__compareUintAscending);
    }
  }

  //sanity check & ordering
  FOR1D(coveredExperiences_per_outcome_slim, i) {
    CHECK(!coveredExperiences_per_outcome_slim(i).containsDoubles(), "coveredExperiences contains double");
    coveredExperiences_per_outcome_slim.elem(i).reverse();
  }

  FOR1D(r->probs, i) {
    CHECK(!isnan(r->probs(i)), "outcome " << i << " is NaN!");
  }

  if (DEBUG>0) cout << "induceOutcomes [END]" << endl;
}





/************************************************
 *
 *     learn -- details -- parameters
 *
 ************************************************/


// Return value needs to be MINIMIZEd.
double learn::learn_parameters(const MT::Array< LitL >& outcomes, doubleA& probs, uint DEBUG) {
  
  if (DEBUG > 0) {
    cout << endl;
    cout << "Learning Parameters [START]" << endl;
    PRINT(__pen_sum__base);
    PRINT(__pen_pos__base);
  }

  uint i;
  
  // ATTENTION WE ASSUME COVERAGE HAS BEEN SET CORRECTLY BEFORE
  FOR1D(outcomes, i) {
    if (i == outcomes.N-1)  // omit noise outcome
      continue;
    CHECK(sum(CostFunction::__cf_coverage_outcome_experience.sub(i,i,0,CostFunction::__cf_coverage_outcome_experience.d1-1)), "At least one experience should be covered for outcome i="<<i<<"!");
  }
  
  // ----------------------------
  // INITIALIZATION
  double cost, oldCost, diff_cost=10.0;
  // init probs
  probs.resize(outcomes.N);
  // default prob = 1 / 2N
  if (outcomes.N > 1)
    probs.last() = 1. / (5. * probs.N);
  else {
    probs.last() = 1.0;
    cost = CostFunction::calc(probs);
    return cost;
  }
  // non-default probs = 1/(N-1) (1 - default_prob)
  for (i=0; i<outcomes.N-1; i++) {
    probs(i) = (1. / (probs.N-1.)) * (1. - probs.last());
  }
  CHECK(TL::isZero(sum(probs)-1), "Bad init probs!")

      if (DEBUG > 0) {
    PRINT(CostFunction::__cf_coverage_outcome_experience) // cf_coverage_outcome_experience == coverage
  }
  
  
  
  // ----------------------------
  // OPTIMIZATION
  
  // prepare optimization algorithm
  double (*f)(const arr&);
  f = CostFunction::calc;
  void (*df)(arr&,const arr&);
  df = CostFunction::calc_grad;

  //  MT::checkGradient(f, df, probs, 0.05);
  
  arr gradients;
  cost = CostFunction::calc(probs);
  
  if (DEBUG > 0)
    cout << "init_probs = " << probs << " C=" << cost << endl;

  double STOPPING_THRESHOLD = 0.001;
  uint MAX_STEPS = 1000;
  // RProp
  TL::Rprop rp;
  rp.init(0.025);
  rp.dMin = 1e-9;
  i=0;
  if (DEBUG>0) cout<<"RProp:"<<endl;
  while(fabs(diff_cost) > STOPPING_THRESHOLD) {
    CostFunction::calc_grad(gradients, probs);
    rp.step(probs, gradients);
    oldCost = cost;
    cost = CostFunction::calc(probs);
    diff_cost = cost - oldCost;
    if (DEBUG > 0) {
      if (DEBUG>1)
        cout << i << ": " << probs << " C=" << cost << " diff=" << diff_cost << " sum=" << sum(probs) << endl;
      else {
        if (i>1000) {
          cout << i << ": " << probs << " C=" << cost << " diff=" << diff_cost << " sum=" << sum(probs) << endl;
        }
        else {
          if (i%10==0)
            cout << i << ": " << probs << " C=" << cost << " diff=" << diff_cost << " sum=" << sum(probs) << endl;
        }
      }
    }
    i++;
    if (i>MAX_STEPS) {
      std::cerr <<endl<<endl<<endl<< "probs = " << probs << endl;
      MT_MSG("WARNING!!! Cannot learn rule probabilities! (No convergence.)");
    }
  }

  
  // ----------------------------
  // POST-PROCESSING
  
  // In the end ensure by hand that pi >=0
  double EPSILON__PROB_NEG1 = 0.10;  // 0.01
  double EPSILON__PROB_NEG2 = 0.15;
  double negativstProb = 1.0;
  FOR1D(probs, i) {
    if (probs(i) < -EPSILON__PROB_NEG2) {
      MT::String warning;
      warning << "Param. optimization failed - Significant negative probability: " << probs(i);
      HALT(warning)
    }
    if (probs(i) < -EPSILON__PROB_NEG1) {
      MT::String warning;
      warning << "Param. optimization awkward - Significant negative probability: " << probs(i);
      __pen_pos__base *= 1.3;
      MT_MSG(warning)
    }
    if (probs(i) < 0  &&  probs(i) < negativstProb)
      negativstProb = probs(i);
  }
  if (negativstProb < 0) {
    FOR1D(probs, i) {
      probs(i) -= negativstProb;
    }
  }
  // In the end, ensure by hand that SUM_i pi = 1
  double probSum = sum(probs);
  probs /= probSum;
  double EPSILON__PROB_SUM1 = 0.15;
  double EPSILON__PROB_SUM2 = 0.25;
  if (fabs(probSum - 1) > EPSILON__PROB_SUM2) {
    MT::String warning;
    warning << "Param. optimization failed - sum clearly different from 1: " << probSum << " found=" << (probs*probSum) << " rescaled=" << probs;
    HALT(warning);
  }
  if (fabs(probSum - 1) > EPSILON__PROB_SUM1) {
    MT::String warning;
    warning << "Param. optimization awkward - sum clearly different from 1: " << probSum << " found=" << (probs*probSum) << " rescaled=" << probs;
    MT_MSG(warning);
    __pen_sum__base *= 1.3;
  }
  
  // In the end, ensure that noise outcome (= last outcome) has non-zero prob.
  probs.last() += 1e-5;
  FOR1D(probs, i) {
    if (probs(i) > 1e-5) {
      probs(i) -= 1e-5;
      break;
    }
  }

  if (DEBUG > 0) {
    cout << "Learned params: " << probs << endl;
    cout << "Learning Parameters [END]" << endl;
  }

  // determine final score
  return cost;
}



/************************************************
 *
 *     learn -- details -- cost function
 *
 ************************************************/


void learn::CostFunction::setRuleCoveredExperiences(const StateTransitionL& coveredEx) {
  __cf_experiences_coveredByCurrentRule.clear();
  __cf_experiences_coveredByCurrentRule= coveredEx;
}

void learn::CostFunction::setOutcomesCoverage(const boolA& coverage) {
  __cf_coverage_outcome_experience = coverage;
}


double learn::CostFunction::loglikelihood(const arr& probs) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"loglikelihood [START]"<<endl;
  double loglikelihood = 0.0, innerSum;
  uint i,e;
  FOR1D(__cf_experiences_coveredByCurrentRule, e) {
    innerSum = 0.0;
    if (DEBUG>0) cout<<"ex "<<e<<":"<<endl;
    FOR1D(probs, i) {
      if (__cf_coverage_outcome_experience(i, e)) {
        if (i < probs.N - 1) {
          innerSum += probs(i);
          if (DEBUG>0) cout<<"  o"<<i<<" 1 * "<<probs(i)<<endl;
        }
        else { // noise outcome
          innerSum += __p_min * probs(i);
          if (DEBUG>0) cout<<"  o"<<i<<" "<<__p_min<<" * "<<probs(i)<<endl;
        }
      }
    }
    if (DEBUG>0) cout<<" lik="<<innerSum<<"   log(lik)="<<log(innerSum)<<endl;
    loglikelihood += log(innerSum);
  }
  if (DEBUG>0) PRINT(loglikelihood)
      if (DEBUG>0) cout<<"loglikelihood [END]"<<endl;
  return loglikelihood;
}


// see Latex-paper for mathematical details
// to minimize!!!
double learn::CostFunction::calc(const arr& in) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"CostFunction::calc [START]"<<endl;
  CHECK(in.N == __cf_coverage_outcome_experience.d0, "invalid number of arguments")
      uint i;
  // calc log-likelihood
  double loglik = loglikelihood(in);
  double sumConstraint = __pen_sum__final * pow(sum(in) - 1.0, 2);
  double posConstraint = 0.0;
  FOR1D(in, i) {
    if (in(i) < 0.0) {
      posConstraint += pow(in(i), 2);
    }
  }
  posConstraint *= __pen_pos__final;
  double cost = -loglik + sumConstraint + posConstraint;
  if (DEBUG>0) cout<<"cost="<<cost<<" (-loglik="<<-loglik<<", sumConstraint="<<sumConstraint<<", posConstraint="<<posConstraint<<")"<<endl;
  if (DEBUG>0) cout<<"CostFunction::calc [END]"<<endl;
  return cost;
}


// points into direction of steepest ascent
void learn::CostFunction::calc_grad(arr& out, const arr& in) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"calc_grad [START]"<<endl;
  CHECK(in.N == __cf_coverage_outcome_experience.d0, "invalid number of arguments")
      out.resize(in.N);
  double sumIn = sum(in);
  double loglik_grad, const1_grad, const2_grad, denom_sum;
  uint i, i2, e;
  // calc gradient for each prob
  FOR1D (in, i) {
    // check for all experiences
    loglik_grad = 0.0;
    if (i < in.N - 1) {
      CHECK(sum(__cf_coverage_outcome_experience.sub(i,i,0,__cf_coverage_outcome_experience.d1-1)), "At least one experience should be covered!");
    }
    FOR1D(__cf_experiences_coveredByCurrentRule, e) {
      // I[covers(s', o_i)]
      if (__cf_coverage_outcome_experience(i, e)) {
        // denominator sum
        denom_sum = 0.0;
        FOR1D(in, i2) {
          if (__cf_coverage_outcome_experience(i2, e)) {
            if (i2 < in.N - 1) // non-noise
              denom_sum += in(i2);
            else // noise outcome
              denom_sum += __p_min * in(i2);
          }
        }
        if (i < in.N-1) // non-noise
          loglik_grad += 1. / denom_sum;
        else // noise
          loglik_grad += __p_min / denom_sum;
      }
    }
    
    // constraint 1
    const1_grad = 2. * __pen_sum__final * (sumIn - 1.0) * in(i);
    
    // constraint 2
    if (in(i) < 0.0) {
      const2_grad = __pen_pos__final * 2. * in(i);
    }
    else
      const2_grad = 0.0;
    out(i) = - loglik_grad + const1_grad + const2_grad;

    if (DEBUG>1) {
      cout << i << ": grad=" << out(i) << " -loglik=" << (-loglik_grad) << " sumPen=" << const1_grad << " posPen=" << const2_grad << endl;
    }

  }
  if (DEBUG>0) cout<<"calc_grad [END]"<<endl;
}


}
