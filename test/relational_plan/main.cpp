#define MT_IMPLEMENT_TEMPLATES

#include <relational/prada.h>
#include <relational/robotManipulationSymbols.h>



#define PLAN_TYPE__PRADA 1
#define PLAN_TYPE__A_PRADA 2
#define PLAN_TYPE__UCT 3
#define PLAN_TYPE__SST 4



void test_plan() {
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
  
  
  double SST_noise_scaling_factor;
  MT::getParameter(SST_noise_scaling_factor, "noise_scaling_factor");
  PRINT(SST_noise_scaling_factor);
  
  uint SST_branch;
  MT::getParameter(SST_branch, "SST_branch");
  PRINT(SST_branch);
  
  uint SST_horizon;
  MT::getParameter(SST_horizon, "SST_horizon");
  PRINT(SST_horizon);
        
  
  uint PRADA_horizon;
  MT::getParameter(PRADA_horizon, "PRADA_horizon");
  PRINT(PRADA_horizon);
  
  uint PRADA_num_samples;
  MT::getParameter(PRADA_num_samples, "PRADA_num_samples");
  PRINT(PRADA_num_samples);
  
  double PRADA_noise_softener;
  MT::getParameter(PRADA_noise_softener, "PRADA_noise_softener");
  PRINT(PRADA_noise_softener);
  
  
  uint UCT_horizon;
  MT::getParameter(UCT_horizon, "UCT_horizon");
  PRINT(UCT_horizon);
  
  uint UCT_c;
  MT::getParameter(UCT_c, "UCT_c");
  PRINT(UCT_c);
  
  uint UCT_numEpisodes;
  MT::getParameter(UCT_numEpisodes, "UCT_numEpisodes");
  PRINT(UCT_numEpisodes);
  
  uint horizon;
  if (plan_type == PLAN_TYPE__SST) horizon = SST_horizon;
  else if (plan_type == PLAN_TYPE__UCT) horizon = UCT_horizon;
  else horizon = PRADA_horizon;
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
  relational::reason::setConstants(s.state_constants);
  cout<<"CONSTANTS:"<<endl;  PRINT(relational::reason::getConstants());
  
  
  
  // -------------------------------------
  // REWARD
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"REWARD: "<<endl;
  relational::Reward* reward = NULL;
  // (1) SIMPLE REWARD FOR TOWER
  //     Create here by hand
  // (1a) LiteralReward
//   relational::Literal* lit = relational::Literal::get(relational::Symbol::get("on"), TUP(66, 67), 1.);
//   reward = new relational::LiteralReward(lit);
  // (1b) LiteralListReward
  relational::LitL lits_reward;
  relational::Literal::get(lits_reward, MT::String("on(66 69) on(69 67)"));
  reward = new relational::LiteralListReward(lits_reward);

  // (2) STACKING REWARD 
  //     Use specification in robotManipulationInterface
//   reward = relational::RobotManipulationSymbols::RewardLibrary::stack();

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

  
  // Example: check which rules cover the current state
  relational::RuleSet coveringGroundRules;
  relational::reason::calc_coveringRules(coveringGroundRules, rules, s);
//   cout<<endl<<endl<<"COVERING RULE GROUNDINGS FOR START STATE (#="<<coveringGroundRules.num()<<"):"<<endl;
//   coveringGroundRules.write();
//   if (coveringGroundRules.num() == 0  
//     ||  (coveringGroundRules.num() == 1 && relational::Rule::isDefaultRule(coveringGroundRules.elem(0)))) {
//     cout<<"No covering rules for current state!"<<endl;
//   }
//   cout<<endl<<endl<<endl;
  cout<<endl<<endl<<"Ground actions with unique covering rules:"<<endl;
  uint i;
  FOR1D_(coveringGroundRules, i) {
    cout<<*coveringGroundRules.elem(i)->action<<" ";
  }
  cout<<endl;
  
  
  // -------------------------------------
  // GROUND THE RULES
  // -------------------------------------
  
  relational::RuleSet ground_rules;
  relational::RuleSet::ground_with_filtering(ground_rules, rules, relational::reason::getConstants(), s);
  cout<<endl<<endl;
  cout<<"GROUND RULES: (plenty!!)"<<endl;
  cout<<"# = "<<ground_rules.num()<<endl;
//   ground_rules.write();
  
  
  
  // -------------------------------------
  // PLANNERs
  // -------------------------------------
  
  cout<<endl<<endl;
  cout<<"PLANNER:"<<endl;
  relational::NID_Planner* planner = NULL;
  if (plan_type == PLAN_TYPE__SST) {
    planner = new relational::NID_SST(SST_branch);
    planner->setNoiseScalingFactor(SST_noise_scaling_factor);
    cout<<"Planner: SST"<<endl;
  }
  else if (plan_type == PLAN_TYPE__PRADA) {
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
  else if (plan_type == PLAN_TYPE__UCT) {
    planner = new relational::NID_UCT();
    ((relational::NID_UCT*) planner)->setNumEpisodes(UCT_numEpisodes);
    ((relational::NID_UCT*) planner)->setC(UCT_c);
    planner->setNoiseScalingFactor(SST_noise_scaling_factor);
    cout<<"Planner: UCT"<<endl;
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
  relational::Literal* action;
  cout<<endl<<endl<<"*** Planning for a single action."<<endl;
  action = planner->plan_action(s);
  
  if (action != NULL)
    cout<<"The planner would like to kindly recommend the following action to you:"<<endl<<*action<<endl;
  else
    cout<<"No action has been found."<<endl;
	 	
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
  MT::initCmdLine(argc, argv);
  cout.precision(3);
  MT::String config_file;
  MT::getParameter(config_file, MT::String("confFile"), MT::String("config"));
  cout << "Config-file: " << config_file << endl;
  MT::openConfigFile(config_file);
  test_plan();
  return 0;
}

