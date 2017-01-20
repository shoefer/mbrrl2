#define MT_IMPLEMENT_TEMPLATES

#include <ios>
#include <relational/learn.h>
#include <relational/prada.h>
#include <relational/robotManipulationSymbols.h>



// ------------------------------------------------------------------
// ------------------------------------------------------------------
//   learn
// ------------------------------------------------------------------
// ------------------------------------------------------------------


void learn_rules() {
  // Rule learning algorithm is heuristic and makes some random choices.
  rnd.seed(12345);
        
  // -------------------------------------
  //  PARAMETERS
  // -------------------------------------
  // Regularizer
  double alpha_pen = 0.5;
  // Lower bounds for probabilities of states in case of noise outcome
  double prob_state_given_NoisyOutcome = 1e-8; // p_min
  // ... same, only for noisy default rule
  double prob_state_given_NoisyOutcome__in_noisyDefaultRule = 1e-9;
  // Log-file
  MT::String logfile("rulelearn.log");
        
 
  
  // -------------------------------------
  //  SET UP LOGIC
  // -------------------------------------

  // Symbols
  relational::SymL symbols;
  relational::ArgumentTypeL types;
  relational::readSymbolsAndTypes(symbols, types, "scaleSymbols.dat");

  // -------------------------------------
  //  READ DATA
  // -------------------------------------
  
  relational::StateTransitionL transitions = relational::StateTransition::read_SASAS("scaleSamples.dat");

  // -------------------------------------
  //  LEARN
  // -------------------------------------

  time_t timeStart, timeEnd;
  time(&timeStart);
	
  relational::learn::set_penalty(alpha_pen);
  relational::learn::set_p_min(prob_state_given_NoisyOutcome, prob_state_given_NoisyOutcome__in_noisyDefaultRule);
  relational::RuleSetContainer rulesC;
  cout<<"Starting rule-learning... (might take quite a while; watch '"<<logfile<<"')"<<endl;
  relational::learn::learn_rules(rulesC, transitions); 
  
  time(&timeEnd);
  cout << "Rule learning finished after " << std::fixed << difftime(timeEnd, timeStart) << " seconds." << endl;
  
  relational::write(rulesC.rules, "learned_rules.dat");
  cout<<"Learned rules have been written to 'learned_rules.dat'."<<endl;
}










// ------------------------------------------------------------------
// ------------------------------------------------------------------
//   plan
// ------------------------------------------------------------------
// ------------------------------------------------------------------


#define PLAN_TYPE__PRADA 1
#define PLAN_TYPE__A_PRADA 2

void test_plan() {
  MT::String config_file;
  MT::getParameter(config_file, MT::String("confFile"), MT::String("config"));
  cout << "Config-file: " << config_file << endl;
  MT::openConfigFile(config_file);
  
  cout<<"********************************"<<endl;
  cout<<" libPRADA plan demo"<<endl;
  cout<<" Domain: Robot Manipulation"<<endl;
  cout<<" This demo shows how to set up a planning problem and solve it."<<endl
      <<" It uses the robot manipulation domain of the experiments in"<<endl
      <<" the paper Lang & Toussaint, JAIR (2010)."<<endl;
  cout<<"********************************"<<endl;
  
  
  // -------------------------------------
  // READ CONFIG
  // -------------------------------------
  
  cout<<endl<<"READ CONFIG:"<<endl;
  
  uint randSeed;
  MT::getParameter(randSeed, "randSeed");
  rnd.seed(randSeed);
  PRINT(randSeed);
  
  uint plan_type;
  MT::getParameter(plan_type, "plan_type");
  PRINT(plan_type);
  
  double discountFactor;
  MT::getParameter(discountFactor, "discountFactor");
  PRINT(discountFactor);
        
  
  uint PRADA_horizon;
  MT::getParameter(PRADA_horizon, "PRADA_horizon");
  PRINT(PRADA_horizon);
  
  uint PRADA_num_samples;
  MT::getParameter(PRADA_num_samples, "PRADA_num_samples");
  PRINT(PRADA_num_samples);
  
  double PRADA_noise_softener;
  MT::getParameter(PRADA_noise_softener, "PRADA_noise_softener");
  PRINT(PRADA_noise_softener);
    
  uint horizon = PRADA_horizon;
  PRINT(horizon);
  
    
  MT::String rulesFile_name;
  MT::getParameter(rulesFile_name, "file_rules");
  PRINT(rulesFile_name);
  
  MT::String stateFile_name;
  MT::getParameter(stateFile_name, "file_state");
  PRINT(stateFile_name);
  
  MT::String rewardFile_name;
  MT::getParameter(rewardFile_name, "file_reward");
  PRINT(rewardFile_name);
  
  MT::String symbolsFile_name;
  MT::getParameter(symbolsFile_name, "file_symbols");
  PRINT(symbolsFile_name);
  

  // -------------------------------------
  //  SET UP LOGIC
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"SYMBOLS:"<<endl;
  cout<<"Reading symbols from file \""<<symbolsFile_name<<"\"..."<<flush;
  relational::SymL symbols;
  relational::ArgumentTypeL types;
  relational::readSymbolsAndTypes(symbols, types, symbolsFile_name);
  cout<<"done!"<<endl;
  
  relational::writeSymbolsAndTypes("used_symbols.dat");


  // -------------------------------------
  //   STATE
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"STARTING STATE:"<<endl;
  cout<<"Reading state from file \""<<stateFile_name<<"\"... "<<flush;
  relational::SymbolicState s;
  ifstream in_state(stateFile_name);
  s.read(in_state);
  cout<<"done!"<<endl<<endl;
  cout<<"State:"<<endl<<s<<endl<<endl;
  
  // -------------------------------------
  //   TYPES
  // -------------------------------------

  relational::ArgumentType *scale = relational::ArgumentType::get(MT::String("scale"));
  relational::ArgumentType *abstractConcept = relational::ArgumentType::get(MT::String("abstractConcept"));
  relational::ArgumentType *physical_free = relational::ArgumentType::get(MT::String("physical_free"));
  relational::ArgumentType *physical_fixed = relational::ArgumentType::get(MT::String("physical_fixed"));
  relational::ArgumentTypeL physicalTypes;
  physicalTypes.append(physical_free);
  physicalTypes.append(physical_fixed);
  relational::ArgumentType *physical = relational::DisjunctionArgumentType::get(MT::String("physical"), physicalTypes);

  relational::Symbol *grab = relational::Symbol::get("grab");
  relational::Symbol *puton = relational::Symbol::get("puton");
  cout << "grab arity " << grab->arity << endl;
  grab->arg_types.append(physical_free);
  puton->arg_types.append(physical);

  relational::Symbol *scaleBalance = relational::Symbol::get("scaleBalance");
  scaleBalance->arg_types.append(scale);

  uintA objs = s.state_constants;
  relational::ArgumentTypeL objTypes;
  uint i;
  FOR1D(objs, i) {
    if (objs(i) == 71)
      objTypes.append(scale);
    else if (objs(i)-1 == 71 || objs(i)-2 == 71)
      objTypes.append(abstractConcept);
    else if (objs(i)-3 == 71 || objs(i)-4 == 71 || objs(i) == 60)
      objTypes.append(physical_fixed);
    else objTypes.append(physical_free);

  }
  cout << "Set constants " << objs << endl;
  relational::reason::setConstants(objs, objTypes); 
  
  // -------------------------------------
  // REWARD
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"REWARD: "<<endl;
  relational::Reward* reward = NULL;
  // Conjunction of literal list reward and maximize reward
  relational::LitL rewardConds;
  relational::Literal::get(rewardConds, MT::String("scaleBalance(71)=4"));

  relational::MaximizeReward *maxReward = new relational::MaximizeReward(relational::Literal::get("numAbove(74)=1"));
  maxReward->offset = 0;

  reward = new relational::RewardConjunction();
  ((relational::RewardConjunction*)reward)->addReward(new relational::LiteralListReward(rewardConds));
  ((relational::RewardConjunction*)reward)->addReward(maxReward);
  reward->write();
  
  
  // -------------------------------------
  // RULES
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"RULES:"<<endl;
  relational::RuleSet rules;
  relational::RuleSet::read(rulesFile_name, rules);
  cout<<"Rules successfully read from file \""<<rulesFile_name<<"\"."<<endl;
  cout<<"Rules ("<<rules.num()<<"):"<<endl;   rules.write(cout);
  
#if 1
  // Manually create "doNothing"-rule if desired
  relational::Rule* rule_doNothing = relational::Rule::getDoNothingRule();
  rules.append(rule_doNothing);
  rule_doNothing->write();
#endif
  
  // -------------------------------------
  // GROUND THE RULES
  // -------------------------------------
  relational::RuleSet ground_rules;
  relational::RuleSet::ground_with_filtering(ground_rules, rules, relational::reason::getConstants(), s);
  cout<<endl<<endl;
  cout<<"GROUND RULES: (plenty!!)"<<endl;
  cout<<"# = "<<ground_rules.num()<<endl;
  // ground_rules.write();
  
  
  
  // -------------------------------------
  // PLANNERs
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"PLANNER:"<<endl;
  relational::NID_Planner* planner = NULL;
  if (plan_type == PLAN_TYPE__PRADA) {
    planner = new relational::PRADA_Planner();
    ((relational::PRADA_Planner* ) planner)->setNumberOfSamples(PRADA_num_samples);
    ((relational::PRADA_Planner* ) planner)->setNoiseSoftener(PRADA_noise_softener);
    cout<<"Planner: PRADA"<<endl;
  }
  else if (plan_type == PLAN_TYPE__A_PRADA) {
    planner = new relational::A_PRADA();
    ((relational::A_PRADA* ) planner)->setNumberOfSamples(PRADA_num_samples);
    ((relational::A_PRADA* ) planner)->setNoiseSoftener(PRADA_noise_softener);
    cout<<"Planner: A-PRADA"<<endl;
  }
  else {
    HALT("no planner defined");
  }
  planner->setDiscount(discountFactor);
  planner->setHorizon(horizon);
  planner->setGroundRules(ground_rules);
  planner->setReward(reward);
  cout<<"Planner has been fully set up."<<endl;cout<<endl;
  PRINT(planner->ground_actions);
  
  
 
  
  
  // -------------------------------------
  //    PLANNING
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"PLANNING:"<<endl;
    
  if (plan_type == PLAN_TYPE__PRADA || plan_type == PLAN_TYPE__A_PRADA) {
    cout<<endl<<endl<<"*** Planning for a complete plan."<<endl;
    relational::LitL plan;
    double value;
    ((relational::PRADA_Planner*) planner)->plan_full(plan, value, s);
    cout<<endl;
    cout<<"Also, it would like to kindly recommend the following plan with value "<<value<<" to you:"<<endl<<plan<<endl;
  }
}




int main(int argc, char** argv){
  cout.precision(2);
  cout << std::fixed;
  // learn
  learn_rules();
  // plan
  test_plan();
  return 0;
}

