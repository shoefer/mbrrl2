#define MT_IMPLEMENT_TEMPLATES

#include <relational/learn.h>
#include <relational/score_functions.h>
#include <relational/rule_learner2.h>

#include <string>
#include <sstream>

#include <relational/logging_macros.h>
#include <boost/filesystem.hpp>

#include <ctime>

void test_learn(std::string sf_type, double alpha_pen = 1e-10, double prob_state=1e-12, double beta = 100., int seed=12345, std::string table_file=std::string(""), std::string symbol_file=std::string("")) {
  cout<<"********************************"<<endl;
  cout<<" mbrrl2 task-sensitive rule learning demo"<<endl;
  cout<<" This demo shows the new rule learning interface for"<<endl
      <<"  the algorithm by Pasula et al., JAIR (2007)."<<endl;
  cout<<"********************************"<<endl;

  boost::filesystem::create_directories("logs");
  
  if (sf_type != "ts" && sf_type != "pasula")
	throw "score function type must be ts or pasula!";
	
  if (table_file == "")
	table_file="data.dat";
  
  if (symbol_file == "") 
	symbol_file = "symbols.dat";
  
  // Rule learning algorithm is heuristic and makes some random choices.
  rnd.seed(12345);
		
  // -------------------------------------
  //  PARAMETERS
  // -------------------------------------
  // Regularizer
  // Lower bounds for probabilities of states in case of noise outcome
  double prob_state_given_NoisyOutcome = prob_state; // p_min
  // ... same, only for noisy default rule
  double prob_state_given_NoisyOutcome__in_noisyDefaultRule = prob_state;
  // Log-file
  std::stringstream stype_prefix;
  stype_prefix << "logs/" << sf_type << "_";
  MT::String logfile((stype_prefix.str() + "learn.log").c_str());
	
  // Symbols
  cout<<"Reading symbols"<<endl;
  relational::SymL symbols;
  relational::ArgumentTypeL types;
  relational::readSymbolsAndTypes(symbols, types, symbol_file.c_str());
	
  // Data
  cout<<"Reading data"<<endl;
  relational::StateTransitionL transitions = relational::StateTransition::read_SAS_SAS(table_file.c_str());
  PRINT(transitions.N);
//   write(transitions);
 
  // -------------------------------------
  //  LEARN
  // -------------------------------------
	
  // Set the hyperparameters for Pasula's score function
  relational::learn::set_penalty(alpha_pen);
  relational::learn::set_p_min(prob_state_given_NoisyOutcome, prob_state_given_NoisyOutcome__in_noisyDefaultRule);
  
  // Instantiate the score functions
  // Pasula's
  relational::learn::ScoreFunction* pasula_sf = new relational::learn::PasulaScoreFunction;
  // Task-sensitive score
  relational::learn::ScoreFunction* ts_sf = new relational::learn::WeightVarianceScoreFunction(beta, 0, 100*beta, 0, 0);
  
  if (sf_type == "ts")
    relational::learn::ScoreFunction::activateScoreFunction(ts_sf);
  else
    relational::learn::ScoreFunction::activateScoreFunction(pasula_sf);
  
  relational::learn::ScoreFunction* sf = relational::learn::ScoreFunction::getCurrentScoreFunction();
  cout << "Score function " << sf->name() << endl;
  cout << " alpha: " << dynamic_cast<relational::learn::PasulaScoreFunction*>(sf)->getAlpha() << endl;
  cout << " pmin: " << dynamic_cast<relational::learn::PasulaScoreFunction*>(sf)->getPmin() << endl;
  cout << " beta: " << dynamic_cast<relational::learn::WeightVarianceScoreFunction*>(ts_sf)->getBeta() << endl;
  cout << " beta_NID: " << dynamic_cast<relational::learn::WeightVarianceScoreFunction*>(ts_sf)->getNIDPenalty() << endl;
 
  /*
  // OLD interface
  relational::RuleSetContainer rulesC;
  cout<<"Starting rule-learning... (might take quite a while; watch 'learn.log' and 'verbose_learn.log')"<<endl;

  clock_t begin = clock();
  relational::learn::learn_rules(rulesC, transitions); 
  clock_t end = clock();
  printf("Learning took %.5f seconds\n", (double(end - begin) / CLOCKS_PER_SEC));
  
  relational::write(rulesC.rules, (stype_prefix.str() + "learned_rules.dat").c_str() );
  rulesC.write( (stype_prefix.str() + "current_learned_rules.dat").c_str() );
  cout<<"Learned rules have been written to '" << stype_prefix.str() << "learned_rules.dat'."<<endl;
  
  // -------------------------------------
  // EVALUATE
  double score = ts_sf->score(rulesC, transitions);
  printf(" learn -  task-sensitive RULE score %10.5f\n", score);

  score = pasula_sf->score(rulesC, transitions);
  printf(" learn -  pasula RULE score %10.5f\n", score);
  */

  // NEW interface
  relational::RuleLearner2 rsc;
  
  // This switch accelerates learning significantly if you have
  // many redundant, identical experiences
  rsc.setCompressBeforeLearning(true);
  
  clock_t begin = clock();
  rsc.learn(transitions, (stype_prefix.str() + "learned_rules.dat").c_str());
  clock_t end = clock();
  
  printf("Learning took %.5f seconds\n", (double(end - begin) / CLOCKS_PER_SEC));
  
  rsc.write((stype_prefix.str() + "current_learned_rules.dat").c_str() );

  // -------------------------------------
  // EVALUATE
  
  // The new interface can also directly call the scoe functions
  double score = rsc.score(ts_sf);
  printf(" learn -  task-sensitive RULE score %10.5f\n", score);

  score = rsc.score(pasula_sf);
  printf(" learn -  pasula RULE score %10.5f\n", score);
}



int main(int argc, char** argv){
	if (argc < 2) {
		cout << "Usage: ./x.exe pasula|ts [alpha_pen] [p_min] [beta] [seed] [table_file] [symbols_file]" << endl;
		return 0;
	}

  cout.precision(2);
  
  if (argc == 2) test_learn(argv[1]);
  if (argc == 3) test_learn(argv[1], atof(argv[2]));
  if (argc == 4) test_learn(argv[1], atof(argv[2]), atof(argv[3]));
  if (argc == 5) test_learn(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]));
  if (argc == 6) test_learn(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]), atoi(argv[5]));
  if (argc == 7) test_learn(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]), atoi(argv[5]), argv[6]);
  if (argc == 8) test_learn(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]), atoi(argv[5]), argv[6], argv[7]);

  return 0;
}

