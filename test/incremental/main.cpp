#define MT_IMPLEMENT_TEMPLATES

#include <relational/learn.h>
#include <relational/score_functions.h>

using namespace relational;

void test_learn() {
  cout<<"********************************"<<endl;
  cout<<" libPRADA learn demo"<<endl;
  cout<<" This demo shows how to learn probabilistic relational rules"<<endl
      <<" with the algorithm by Pasula et al., JAIR (2007)."<<endl
      <<" The demo uses the robot manipulation domain of the experiments in"<<endl
      <<" the paper Lang & Toussaint, JAIR (2010)."<<endl;
  cout<<"********************************"<<endl;
	
  // Rule learning algorithm is heuristic and makes some random choices.
  rnd.seed(12345);
	
	
  // -------------------------------------
  //  PARAMETERS
  // -------------------------------------
  // Regularizer
  double alpha_pen = 1e-10;
  // Lower bounds for probabilities of states in case of noise outcome
  double prob_state_given_NoisyOutcome = 1e-15; // p_min
  // ... same, only for noisy default rule
  double prob_state_given_NoisyOutcome__in_noisyDefaultRule = 1e-15;

  relational::learn::ScoreFunction::activateScoreFunction(
        new relational::learn::WeightVarianceScoreFunction(
          0.,0., 1., 10, true));

  // Log-file
  MT::String logfile("learn.log");
	
  // Symbols
  relational::SymL symbols;
  relational::ArgumentTypeL types;
  relational::readSymbolsAndTypes(symbols, types, "symbols.log");
	
  // Data
  cout<<"Loading first set of experiences from 'transitions.log'."<<endl;
  relational::StateTransitionL transitions = relational::StateTransition::read_SAS_SAS("transitions.log");
  PRINT(transitions.N);
//   write(transitions);

  cout<<"Loading second, additional set of experiences from 'transitions.log'."<<endl;
  relational::StateTransitionL transitions_new = relational::StateTransition::read_SAS_SAS("transitions_new.log");
  PRINT(transitions_new.N);

  //cout << "constants:" << endl;
  //relational::reason::getConstants().write(cout);
 
  // -------------------------------------
  //  LEARN
  // -------------------------------------
	
  relational::learn::set_penalty(alpha_pen);
  relational::learn::set_p_min(prob_state_given_NoisyOutcome, prob_state_given_NoisyOutcome__in_noisyDefaultRule);
  relational::RuleSetContainer rulesC;
  cout<<"Starting rule-learning on first set... (might take quite a while; watch '"<<logfile<<"')"<<endl;
  relational::learn::learn_rules(rulesC, transitions,
                                 false, // incremental
                                 "relearn.log", "relearn_verbose.log", 1);
  
  relational::write(rulesC.rules, "learned_rules.dat");
  rulesC.write("learned_RSC.dat");
  cout<<"Learned rules have been written to 'learned_rules.dat'."<<endl;

  // -------------------------------------
  //  LEARN INCREMENTALLY
  // -------------------------------------
  cout<<"Starting rule-learning on second, additional set... (might take quite a while; watch '"<<logfile<<"')"<<endl;
  relational::learn::learn_rules(rulesC,
                                 transitions_new,
                                 true, // incremental
                                 "relearn.log", "relearn_verbose.log", 1);
  relational::write(rulesC.rules, "relearned_rules.dat");
  rulesC.write("relearned_RSC.dat");
  cout<<"Re-Learned rules have been written to 'relearned_rules.dat'."<<endl;

}



int main(int argc, char** argv){
  cout.precision(2);
  test_learn();
  
  return 0;
}

