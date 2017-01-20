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

#ifndef RELATIONAL_literals_h
#define RELATIONAL_literals_h

#include <MT/util.h>
#include <relational/symbols.h>



namespace relational {
  
/************************************************
 * 
 *     Literal
 * 
 ************************************************/


class Literal;
class SymbolicState;
class StateTransition;
typedef MT::Array< Literal* > LitL;
typedef MT::Array< SymbolicState* > SymbolicStateL;
typedef MT::Array< StateTransition* > StateTransitionL;



class Literal {
  protected:
    Literal();
  
  public:
    // -----------------------------------
    // Essential data-fields
    Symbol* s;
    uintA args;
    double value;
    
    enum ComparisonType { comparison_equal, comparison_less, comparison_lessEqual, comparison_greater, comparison_greaterEqual, comparison_variable, comparison_offset };
    ComparisonType comparison_type;
    
    // -----------------------------------
    // Convenience methods
    virtual ~Literal() {}
    static Literal* get(Symbol* s, const uintA& args, double value, ComparisonType comparison_type = comparison_equal);
    static Literal* getVarComparison(Symbol* s, const uintA& args);
    static Literal* get(const char* text);
    static void get(MT::Array< Literal* >& lits, const char* text);
    bool operator==(Literal& lit) const;
    bool operator!=(Literal& lit) const;
    virtual void write(ostream& os = cout, bool withTypes = false) const;
    bool compareValueTo(Literal::ComparisonType compType, double a);
    static bool compare(double a, Literal::ComparisonType compType, double b);

    ///Checks whether all arguments have correct types according to the symbol
    bool typeCheck();
    
    bool isNegated() const;
    
    Literal* getNegated();
    static void getLiterals(LitL& lits, Symbol* s, const uintA& constants, double value, bool withRepeatingArguments = true);
    static void getLiterals_state(LitL& lits, const uintA& constants, double value, bool binaryOnly = false);
    static void getLiterals_state(LitL& lits, const uintA& constants, const uintA& constants_mustBeContained, double value, bool binaryOnly = false);
    static void getLiterals_actions(LitL& lits_actions, const uintA& constants);  // get all action atoms
    static Literal* getLiteral_doNothing();
    static Literal* getLiteral_default_action();
    
    // lists
    static void sort(LitL& lits); // order according to (i) symbol dependencies and (ii) binary value=0 (negative)
    static double getValue(const LitL& literals, const Symbol* s, const uintA& args);
    static uint getArguments(uintA& args, const LitL& lits);
    static bool equivalent(const LitL& lits1, const LitL& lits2);
    static bool nonContradicting(const LitL& l1, const LitL& l2);
    // TODO the below methods seem to be not required anymore
    static int findPattern(const LitL& actions, uint minRepeats);
    static bool negativeBinaryLiteralsLast(const LitL& lits);
};
void write(const LitL& lits, ostream& os = cout);
void write(const MT::Array< LitL >& list_of_lits_list, ostream& os = cout);





/************************************************
 * 
 *     SymbolicState
 * 
 ************************************************/

struct SymbolicState {
  // -----------------------------------
  // Essential data-fields
  MT::Array<Literal*> lits;
  uintA state_constants; // not necessarily set
  bool derived_lits_are_calculated; // Have literals derived symbols been calculated?
  
  // -----------------------------------
  // Convenience methods
  SymbolicState();
  SymbolicState(const MT::Array<Literal*>& lits);
  ~SymbolicState();
  bool operator==(const SymbolicState& lit) const;
  bool operator!=(const SymbolicState& lit) const;
  void write(ostream& os = std::cout, bool primOnly = false) const;
  virtual std::string toString(bool primOnly = false) const;
  void read(istream& in, bool read_constants = false);
  static uint getArgument(const SymbolicState& state, const Symbol& symbol);
  static void getArguments(uintA& args, const SymbolicState& state, const Symbol& symbol);
  static void getArguments(uintA& args, const SymbolicState& state);
  static double getValue(const Symbol* s, const SymbolicState& state);
  static void getValues(arr& values, const SymbolicState& s, const Symbol& f, const uintA& constants);
  static void filterState_full(SymbolicState& state_filtered, const SymbolicState& state_full, const uintA& filter_constants, bool primOnly = true);  // only "filter_constants" as args
  static void filterState_atleastOne(SymbolicState& state_filtered, const SymbolicState& state_full, const uintA& filter_constants, bool primOnly = true); // at least one in "filter_constants" as arg
  static void calcDifferences(LitL& lits_diff_1to2, LitL& lits_diff_2to1, uintA& changedConstants, const SymbolicState& state1, const SymbolicState& state2);
  static void getRelatedConstants(uintA& constants_related, uint id, bool id_covers_first, const Symbol& pred, const SymbolicState& state);  // pred(X,Y), pred(X,Z) --> {Y,Z}
};





/************************************************
 * 
 *     StateTransition
 * 
 ************************************************/

struct StateTransition {
  // -----------------------------------
  // Essential data-fields
  SymbolicState pre, post;
  Literal* action;
  double reward;
  MT::Array< Literal* > changes;
  uintA changedConstants;
  
  // -----------------------------------
  // Convenience methods
  StateTransition(const SymbolicState& pre, Literal* action, const SymbolicState& post, double reward = 0);
  StateTransition();
  // change by sebastian -> make base class polymorphic
  virtual ~StateTransition();
  bool noChange();
  void calcChanges();
  void write(ostream& os = cout, uint detailLevel = 0) const;
  
  void getConstants(uintA& constants);

  static void write(const StateTransitionL& exs, ostream& os = cout, uint detailLevel = 0);
  static StateTransitionL& read_SASAS(const char* filename);  // file: state1 action1 state2 action2 state3...
  static StateTransitionL read_SAS_SAS(const char* filename);  // file: state1 action1 state1' state2 action2 state2'...

  bool operator==(const StateTransition& lit) const;
  bool operator!=(const StateTransition& lit) const;

};


void getConstants(uintA& constants, const StateTransitionL& transitions);
void getStateSymbols(SymL& symbols, const StateTransitionL& transitions);



} // relational namespace



std::ostream& operator<<(std::ostream& os, const relational::SymbolicState& s);
std::ostream& operator<<(std::ostream& os, const relational::Literal& s);
std::ostream& operator<<(std::ostream& os, const relational::StateTransition& e);
std::ostream& operator<<(std::ostream& os, const relational::LitL& lits);
std::ostream& operator<<(std::ostream& os, const relational::StateTransitionL& e);




#endif // RELATIONAL_literals_h
