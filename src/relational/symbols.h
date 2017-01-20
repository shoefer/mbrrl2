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


#ifndef RELATIONAL_symbols_h
#define RELATIONAL_symbols_h

#include <MT/util.h>
#include <relational/utilTL.h>



namespace relational {
  

/************************************************
 * 
 *     Symbol
 * 
 ************************************************/
  

class ArgumentType;
class Symbol;
typedef MT::Array< Symbol* > SymL;
typedef MT::Array< ArgumentType* > ArgumentTypeL;

class Symbol {
  protected:
    Symbol();
    
  public:
    // -----------------------------------
    // Essential data-fields
    MT::String name;
    uint arity;
    
    enum SymbolType { 
                      action,
                      primitive,
                      conjunction,
                      transclosure,
                      count,
                      avg,
                      max,
                      sum,
                      function_change,
                      function_reward,
                      function_difference
    };
    SymbolType symbol_type;
    
    enum RangeType {
        binary,
        integer_set,
        integers,
        reals
    };
    RangeType range_type;
    uintA range;
    
    ArgumentTypeL arg_types;
    
    // -----------------------------------
    // Convenience methods
    static Symbol* get(const MT::String& name, uint arity, SymbolType symbol_type = primitive, RangeType range_type = binary);
    static Symbol* get(const MT::String& name);
    static Symbol* get(const char* name);
    static Symbol* getDefaultAction();
    static void get(MT::Array< Symbol* >& all_symbols);
    static void get_state(MT::Array< Symbol* >& symbols);
    static void get_state_nonBinary(MT::Array< Symbol* >& symbols);
    static void get_action(MT::Array< Symbol* >& symbols);
    virtual ~Symbol() {}
    virtual void write(ostream& os = cout) const;
    virtual bool operator==(const Symbol& s) const;
    virtual bool operator!=(const Symbol& s) const;
    virtual void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
    
    static MT::Array< Symbol* > mem__all_symbols;
    static bool symbols_sorted;

    static void sort(SymL& symbols);
};



class Literal;
class ConjunctionSymbol : public Symbol {
  protected:
    ConjunctionSymbol() {
      symbol_type = conjunction;
      range_type = binary;
    }
  public:
    // -----------------------------------
    // Essential data-fields
    MT::Array<Literal* > base_literals;
    bool free_vars_all_quantified;
    
    // -----------------------------------
    // Convenience methods
    static ConjunctionSymbol* get(const MT::String& name, uint arity, MT::Array<Literal* > base_literals, bool free_vars_all_quantified = true);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
    void getFreeVars(uintA& freeVars) const;
    void getFreeVars(uintA& freeVars_pos, uintA& freeVars_neg) const;
    
    // auxiliarly data-structure
    boolA redundant__base_literal_with_free_vars;
};




// assume d=2 thus far!
class TransClosureSymbol : public Symbol {
  protected:
    TransClosureSymbol();
  public:
    // -----------------------------------
    // Essential data-fields
    Symbol* base_symbol;
    
    // -----------------------------------
    // Convenience methods
    static TransClosureSymbol* get(const MT::String& name, uint arity, Symbol* base_symbol);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};


class Literal;
class CountSymbol : public Symbol {
  protected:
    CountSymbol();
  public:
    // -----------------------------------
    // Essential data-fields
    Literal* base_literal;
    
    // -----------------------------------
    // Convenience methods
    static CountSymbol* get(const MT::String& name, uint arity, Literal* base_literal);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};


class AverageFunction : public Symbol {
  protected:
    AverageFunction();
  public:
    // -----------------------------------
    // Essential data-fields
    Symbol* base_symbol;
    
    // -----------------------------------
    // Convenience methods
    static AverageFunction* get(const MT::String& name, uint arity, Symbol* base_symbol);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};


class SumFunction : public Symbol {
  protected:
    SumFunction();
  public:
    // -----------------------------------
    // Essential data-fields
    Symbol* base_symbol;
    
    // -----------------------------------
    // Convenience methods
    static SumFunction* get(const MT::String& name, uint arity, Symbol* base_symbol);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};


class MaxFunction : public Symbol {
  protected:
    MaxFunction();
  public:
    // -----------------------------------
    // Essential data-fields
    Symbol* base_symbol;
    
    // -----------------------------------
    // Convenience methods
    static MaxFunction* get(const MT::String& name, uint arity, Symbol* base_symbol);    
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};


class ChangeFunction : public Symbol {
  protected:
    ChangeFunction();
  public:
    // -----------------------------------
    // Essential data-fields
    Symbol* base_symbol;
    
    // -----------------------------------
    // Convenience methods
    static ChangeFunction* get(const MT::String& name, uint arity, Symbol* base_symbol);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};


class RewardFunction : public Symbol {
  protected:
    RewardFunction();
  public:
    // -----------------------------------
    // Essential data-fields
    MT::Array< Literal* > base_literals;
    double reward_value;
    
    // -----------------------------------
    // Convenience methods
    static RewardFunction* get(const MT::String& name, uint arity, MT::Array< Literal* > base_literals, double reward_value);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};

class DifferenceFunction : public Symbol {
  protected:
    DifferenceFunction() {
      symbol_type = function_difference;
      range_type = integers;
    }

  public:
    // -----------------------------------
    // Essential data-fields
    Symbol *baseSymbol;
    Literal *baseFunctionLit1, *baseFunctionLit2;
    MT::Array< Literal* > restrictionLits;

     void getFreeVars(uintA& freeVars) const;
    
    // -----------------------------------
    // Convenience methods
    static DifferenceFunction* get(const MT::String& name, uint arity, Symbol* base_symbol, Literal *baseFunctionLit1, Literal *baseFunctionLit2, const MT::Array< Literal* > &restriction_lits);
    void write(ostream& os) const;
    void getDefiningSymbols(MT::Array< Symbol* > & symbols, bool only_direct_predecessors = true) const;
};


/************************************************
 * 
 *     ArgumentType
 * 
 ************************************************/


class ArgumentType {
  protected:
    ArgumentType() {type = simple;}
    
  public:
    // -----------------------------------
    // Essential data-fields
    enum ArgumentType_Type { simple, disjunction }; 
        
    ArgumentType_Type type;
    MT::String name;
    
    // -----------------------------------
    // Convenience methods
    ~ArgumentType() {}
    virtual bool operator==(const ArgumentType& t) const;
    virtual bool operator!=(const ArgumentType& t) const;
    virtual bool subsumes(const ArgumentType& t) const;
    virtual void write(ostream& os = cout) const;
    static ArgumentType* get(const MT::String& name);
    static ArgumentType* read(ifstream& in);
    static void get(ArgumentTypeL& all_argumentTypes);
    
    static ArgumentTypeL mem_all_argument_types;
};


struct DisjunctionArgumentType : public ArgumentType {
  protected:
    DisjunctionArgumentType() {type = disjunction;}
  
  public:
    // -----------------------------------
    // Essential data-fields
    ArgumentTypeL base_types;
  
    // -----------------------------------
    // Convenience methods
    ~DisjunctionArgumentType() {}
    bool subsumes(const ArgumentType& t) const;
    void write(ostream& os = cout) const;
    static DisjunctionArgumentType* get(const MT::String& name, const ArgumentTypeL& base_types);
};




void writeSymbols(const SymL& symbols, ostream& os = cout);
void writeSymbolsAndTypes(const SymL& symbols, const ArgumentTypeL& types, const char* filename);
void writeSymbolsAndTypes(const SymL& symbols, const ArgumentTypeL& types, ostream& out);
void writeSymbolsAndTypes(const char* filename);

void readSymbolsAndTypes(SymL& symbols, ArgumentTypeL& types, const char *filename);
void readSymbolsAndTypes(SymL& symbols, ArgumentTypeL& types, istream& in);




} // relational namespace



std::ostream& operator<<(std::ostream& os, const relational::Symbol& s);
std::ostream& operator<<(std::ostream& os, const relational::ArgumentType& s);
std::ostream& operator<<(std::ostream& os, const relational::SymL& symbols);


#endif // RELATIONAL_symbols_h
