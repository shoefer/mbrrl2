/*
 * score_functions.h
 *
 *  Created on: Jul 12, 2013
 *      Author: shoefer
 */

#ifndef SCORE_FUNCTIONS_H_
#define SCORE_FUNCTIONS_H_

#include <relational/reason.h>
#include <string>

namespace relational {

class RuleSetContainer;

namespace learn {

class ScoreFunction {
public:
    static void activateScoreFunction(ScoreFunction* sf, bool delete_old=true);
    static ScoreFunction* getCurrentScoreFunction();

  virtual ~ScoreFunction() {}

	virtual double score(RuleSetContainer& rulesC,
			StateTransitionL& experiences,
			double cutting_threshold=TL::TL_DOUBLE_MIN,
			std::ostream& out=std::cout,
			uint DEBUG=0) {
	  arr experience_weights(experiences.N);
	  experience_weights.setUni(1.0);
	  return this->score(rulesC, experiences, cutting_threshold, experience_weights, out, DEBUG);
	}

	virtual double score(RuleSetContainer& rulesC,
			StateTransitionL& experiences,
			double cutting_threshold,
			arr& experience_weights,
			std::ostream& out=std::cout,
			uint DEBUG=0) = 0;

	virtual std::string name() = 0;
};


class PasulaScoreFunction : public ScoreFunction {
public:
	virtual std::string name() {
		return "PasulaScoreFunction";
	}
	double score(RuleSetContainer& rulesC,
				StateTransitionL& experiences,
				double cutting_threshold,
				arr& experience_weights,
				std::ostream& out=std::cout,
				uint DEBUG=0);
    double getAlpha() const;
    double getPmin() const;
    double getPminNID() const;
};

class WeightVarianceScoreFunction : public PasulaScoreFunction {
    double __beta; // weight for penalizing the reward variance of every rule
    double __beta_negrew; // weight for penalizing the likelihood of negative rewards
    double __nid_penalty; // constant factor for penalizing reward variance of NID rule
    double __beta_outcome; // weight for penalizing the reward variance of every outcome

    // multiply the NID rule by the number of experiences currently covered by NID rule.
    // if an experience is taken out of the NID rule, but results in a rule with
    // the same variance as the NID rule, the experience stays in the NID rule.
    // by setting this to true, the NID rule variance is multiplied by the number
    // of experiences covered by the NID rule, yielding in a score decrease if
    // experience is taken out of the NID rule
    // (~ a rule with variance is better than a NID rule with variance)
    bool __weigh_nid_by_no_exp;
public:
    WeightVarianceScoreFunction(double beta, double beta2=0, double nid_penalty=1,
                                double beta_outcome=0, bool weigh_nid_by_no_exp=false)
    : __beta(beta), __beta_negrew(beta2), __nid_penalty(nid_penalty),
      __beta_outcome(beta_outcome),
      __weigh_nid_by_no_exp(weigh_nid_by_no_exp) {}

	virtual std::string name() {
		return "WeightVarianceScoreFunction";
	}
	double score(RuleSetContainer& rulesC,
				StateTransitionL& experiences,
				double cutting_threshold,
				arr& experience_weights,
				std::ostream& out=std::cout,
                uint DEBUG=1);

    double getBeta() const {
        return __beta;
    }
    double getBeta2() const {
        return __beta_negrew;
    }
    double getNIDPenalty() const {
        return __nid_penalty;
    }
};


}
}

#endif /* SCORE_FUNCTIONS_H_ */
