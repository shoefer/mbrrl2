#define MT_IMPLEMENT_TEMPLATES

#include <relational/reason.h>


void test_basics() {
  cout<<"********************************"<<endl;
  cout<<" libPRADA basics demo"<<endl;
  cout<<" This demo shows how to define symbols and rules."<<endl;
  cout<<"********************************"<<endl;
  
  // ----------------------------------
  // SOME SYMBOLS
  relational::Symbol* s_incity = relational::Symbol::get(MT::String("incity"), 1);
  relational::Symbol* s_cruiseto = relational::Symbol::get(MT::String("cruiseto"), 1, relational::Symbol::action);
  relational::Symbol* s_gas = relational::Symbol::get(MT::String("gas"), 0, relational::Symbol::primitive, relational::Symbol::integer_set);
  uintA gas_values;  gas_values.append(0);  gas_values.append(1);  gas_values.append(2);  gas_values.append(3);  gas_values.append(4);
  s_gas->range = gas_values;

  cout<<"SYMBOLS:"<<endl;
  cout<<*s_incity<<endl;
  cout<<*s_cruiseto<<endl;
  cout<<*s_gas<<endl;
  cout<<endl;

  relational::writeSymbolsAndTypes("used_symbols.dat");

  
  // ----------------------------------
  // SOME CONSTANTS
  // (represented by large uints)
  uintA constants;
  constants.append(20); // e.g. city 1 which might be Berlin
  constants.append(21); // e.g. city 2 which might be Frankfurt
  constants.append(22); // e.g. a person which might be Angela Merkel
  
  // Tell the reason-component which uints represent constants
  relational::reason::setConstants(constants);
    
  
  
  // ----------------------------------
  // LITERALS, FUNCTION-VALUES and ATOMS
  uintA arguments(1);
  
  arguments(0) = 0;   // Convention if constants are not set (like above): small uints (<10) denote variables.
  relational::Literal* l1 = relational::Literal::get(s_incity, arguments, 1);
  
  // or write
  relational::Literal* l2 = relational::Literal::get(s_incity, TUP(0), 0);
  
  arguments(0) = constants(2);  // Convention if constants are not set (like above): larger uints (>=10) denote constants.
  relational::Literal* l3 = relational::Literal::get(s_incity, arguments, 1);
  
  arguments(0) = 0;
  relational::Literal* l_action = relational::Literal::get(s_cruiseto, arguments, 1.0);
  
  arguments.clear();  // arguments = {}
  relational::Literal* l4 = relational::Literal::get(s_gas, arguments, 3);
  relational::Literal* l5 = relational::Literal::get(s_gas, arguments, 0);
  relational::Literal* l_comp =  relational::Literal::get(s_gas, arguments, 1., relational::Literal::comparison_greater);
  
  relational::LitL lits_all;   //  list of literals
  lits_all.append(l1);  lits_all.append(l2);  lits_all.append(l3);  lits_all.append(l4);  lits_all.append(l5);  lits_all.append(l_comp);
  
  cout<<"LITERALS: "<<lits_all<<endl<<endl;



  // ----------------------------------
  // SYMBOLIC STATE
  relational::LitL lits_ground;  lits_ground.append(l3);  lits_ground.append(l4);
  relational::SymbolicState state(lits_ground);
  cout<<"STATE:  "<<state<<endl<<endl;



  // ----------------------------------
  // RULES
  
  // RULE 1: noisy default rule
  //  --> Always required!!
  double probability_noise_outcome = 0.7;
  relational::Rule* rule1 = relational::Rule::generateDefaultRule(probability_noise_outcome);
  rule1->noise_changes = 1.7;
  
  // RULE 2
  relational::Rule* rule2 = new relational::Rule;
  rule2->context.append(l_comp);
  rule2->context.append(relational::Literal::get("incity(Y)"));
  rule2->context.append(l2);
  
  rule2->action = l_action;
  
  // Outcome 1:  action succeeded: different city
  relational::LitL outcome1;
  outcome1.append(l1);
  outcome1.append(relational::Literal::get("-incity(Y)"));
  rule2->outcomes.append(outcome1);
  
  // Outcome 2:  action failed: still in same city
  //  (We've driving around, but could not find our way out of the city.)
  relational::LitL outcome2;
  rule2->outcomes.append(outcome2);
  
  // Outcome 3:  noise outcome
  //  (E.g. we ended in a different city, or our car broke on the Autobahn or...)
  relational::LitL outcome3; // empty dummy
  rule2->outcomes.append(outcome3);
  rule2->noise_changes = 2.4; // for PRADA's noise outcome heuristic
  
  arr outcome_probabilities(3);
  outcome_probabilities(0) = 0.7;
  outcome_probabilities(1) = 0.2;
  outcome_probabilities(2) = 0.1;
  rule2->probs = outcome_probabilities;
  
  relational::RuleSet abstract_rules;
  abstract_rules.append(rule1);
  abstract_rules.append(rule2);
  abstract_rules.sanityCheck();
  // The data-structure RuleSet will take care of deleting
  // the rule pointers.
  cout<<"RULES:"<<endl;
  abstract_rules.write();
  
  
  // Rule coverage
  relational::RuleSet coveringGroundRules;
  relational::reason::calc_coveringRules(coveringGroundRules, abstract_rules, state);
  cout<<endl<<endl<<"COVERING RULE GROUNDINGS FOR START STATE (#="<<coveringGroundRules.num()<<"):"<<endl;
  coveringGroundRules.write();
  cout<<endl;


  // Successor state
  cout<<"SAMPLE SUCESSOR STATE"<<endl;
  cout<<"Current state: "<<state<<endl;
  // 	relational::Literal* ground_action = relational::Literal::get("cruiseto(21)");
  MT::String name__ground_action;  name__ground_action<<"cruiseto("<<constants(0)<<")";
  relational::Literal* ground_action = relational::Literal::get(name__ground_action);
  cout<<"Ground action: "<<*ground_action<<endl;
    
  relational::SymbolicState state_successor1;
  relational::reason::sampleSuccessorState_abstractRules(state_successor1, state, abstract_rules, ground_action);
  cout<<"Sampled successor state #1:  "<<state_successor1<<endl;
  
  relational::SymbolicState state_successor2;
  relational::reason::sampleSuccessorState_abstractRules(state_successor2, state, abstract_rules, ground_action);
  cout<<"Sampled successor state #2:  "<<state_successor2<<endl;
}



int main(int argc, char** argv) {
  test_basics();
}

