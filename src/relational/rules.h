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

#ifndef RELATIONAL_rules_h
#define RELATIONAL_rules_h

#include <MT/util.h>
#include <relational/utilTL.h>
#include <relational/symbols.h>
#include <relational/literals.h>

#define FOR1D_(x,i)   for(i=0;i<x.num();i++)
#define FOR1D_DOWN_(x,i)  for(i=x.num();i--;)


namespace relational {
  

/************************************************
 * 
 *     Rule
 * 
 ************************************************/

struct Rule {
  // -----------------------------------
  // Essential data-fields
  Literal* action;
  LitL context;
  // Important convenience for "outcomes" and "probs": the last outcome is the noise outcome.
  MT::Array< LitL > outcomes;
  doubleA probs;
  double noise_changes; // PRADA's noise outcome heuristic: average number of state properties that have changed in case of the noise outcome
  arr outcome_rewards; // optional

  // -----------------------------------
  // Convenience methods
  Rule();
  ~Rule();
  void write(ostream& os = std::cout, bool withAction = true) const;
  void getSymbols(SymL& symbols) const;
  void getArguments(uintA& args) const;
  void getDeicticRefs(uintA& drefs) const;
  void getDeicticRefs(uintA& drefs_pos, uintA& drefs_neg, uintA& drefs_nonBinary) const;
  void getDeicticRefsContains(boolA& containsDR) const;  // calcs which context literals contain deictic refs
  void getAbsentLiterals(LitL& literals, bool positiveOnly = false) const;
  bool usesAtom(Literal* lit) const;
  uint numberLiterals() const;
  uint getNegStartPos();
  uint getNonBinaryStartPos();
  void insertContext(Literal* literal);
  void cleanup();  // order context and outcomes, clean-up DRs
  void copyBody(const Rule& other);
  void sanityCheck();// const;

  ///Checks whether symbol s with the argument args is specified in the outcome
  bool existsInOutcome(Symbol *s, uintA args);
  
  static Rule* read(istream& in);
  
  int numRefs; // used for object management by RuleSet container class
  static uint globalRuleCounter;
  
  // Special rules
  // default rule (which explains everything as noise)
  static Rule* generateDefaultRule(double noiseProb = 0.5, double minProb = 0., double change = 2.0);
  static bool isDefaultRule(const Rule* rule);
  // no-action rule
  static Rule* getDoNothingRule();
  
  private:
    // cannot copy because of numRefs
    DISALLOW_COPY_AND_ASSIGN(Rule);
};







/************************************************
 * 
 *     RuleSet
 * 
 ************************************************/

class RuleSet {
  MT::Array< Rule* > ra;
  
  public:
    RuleSet();
    ~RuleSet();
    void append(Rule* r);
    void remove(uint i);
    Rule* elem(uint i) const;
    int findValue(Rule* r);
    void clear();
    uint num() const;
    void overwrite(uint i, Rule* r);
    RuleSet& operator=(const RuleSet& rs);
    void sort(); // sorts rule-set by actions
    void sort_using_args();
    void sort_using_args(uintA& new_order);
    void write(ostream& os = std::cout) const;

    void sanityCheck() const;
    void changingSymbols(SymL& symbols) const; 
    void removeDoubleLiterals();
    void removeNonChangingSymbols();

    bool hasDefaultRule() const;
  
    static void read(const char* filename, RuleSet& rules);
    static void read(istream& in, RuleSet& rules, bool sort_rules=true);
    
    static void ground(RuleSet& rules_ground, const RuleSet& rules_abstract, const uintA& constants);
    static void ground_with_filtering(RuleSet& rules_ground, const RuleSet& rules_abstract, const uintA& constants, const SymbolicState& s, bool delete_nonchanging_concepts = false);
    static void ground_with_filtering(RuleSet& rules_ground, uintA& abstract_rules_idx, MT::Array<uintA>& abstract2GroundMap, const RuleSet& rules_abstract, const uintA& constants, const SymbolicState& s, bool delete_nonchanging_concepts = false);
    static void groundFunctionVars(const RuleSet& rules, RuleSet& rules_expanded);
};


void write(const RuleSet& rules, ostream& os = cout);
void write(const RuleSet& rules, const char* filename);




/************************************************
 * 
 *     Substitution
 * 
 ************************************************/

class Substitution {
  uintA ins;
  uintA outs;

  public:
  int numRefs;      
  static int globalCounter_Substitution;
      
  Substitution();
  Substitution(const Substitution& s);
  Substitution& operator=(const Substitution& s);
  ~Substitution();
  void apply(LitL& sub_lits, const LitL& unsub_lits) const;
  Literal* apply(Literal* unsub_lit) const;
  SymbolicState* apply(const SymbolicState& state) const;
  Rule* apply(const Rule& Rule) const;
  uint getSubs(uint in) const;
  void getIns(uintA& ids) const;
  void getOuts(uintA& ids) const;
  void getOutsDistinct(uintA& ids) const;
  void addSubs2Variable(uint in);
  void addSubs(uint in, uint out);
  bool hasSubs(uint in) const;
  bool mapsToDistinct() const;
  bool empty() const;
  uint num() const;
  void getInverse(Substitution& invSub) const;
      
  void write(ostream& os = cout) const;
      
  static Substitution* combine(Substitution& sub1, Substitution& sub2);
};



/************************************************
 * 
 *     SubstitutionSet
 * 
 ************************************************/

class SubstitutionSet {
  MT::Array< Substitution* > sa;

public:
  SubstitutionSet();
  ~SubstitutionSet();

  void append(Substitution* s);
  void append(SubstitutionSet& ss);
  void remove(uint i);
  int findValue(Substitution* s);
  void clear();
  uint num() const;
  Substitution* elem(uint i) const;
  Substitution* last() const;
  void overwrite(uint i, Substitution* s);
  
  SubstitutionSet& operator=(const SubstitutionSet& ss);
  
  void write(ostream& os = cout);
  
  static void createAllPossibleSubstitutions(SubstitutionSet& subs, const uintA& vars, const uintA& constants);
};




} // relational namespace



std::ostream& operator<<(std::ostream& os, const relational::Rule& r);
std::ostream& operator<<(std::ostream& os, const relational::Substitution& s);



#endif // RELATIONAL_rules_h
