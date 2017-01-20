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

#include "symbols.h"
#include "rules.h"
#include "reason.h"
#include <stdlib.h>

#define FAST__RULES_REASONING 1


namespace relational {



/************************************************
 * 
 *     Rule
 * 
 ************************************************/

uint Rule::globalRuleCounter = 0;


Rule::Rule() {
    numRefs=0;
    globalRuleCounter++;
    noise_changes = 0.0;
//     cout<<"new rule - now "<<globalRuleCounter<<endl;
}


Rule::~Rule() {
      globalRuleCounter--;
//     cout<<"destroyed rule - now "<<globalRuleCounter<<endl;
}


void Rule::getSymbols(SymL& symbols) const {
  uint i, k;
  FOR1D(this->context, i) {
    symbols.setAppend(this->context(i)->s);
    SymL syms_pre;
    this->context(i)->s->getDefiningSymbols(syms_pre, true);
    symbols.setAppend(syms_pre);
  }
  FOR1D(this->outcomes, i) {
    FOR1D(this->outcomes(i), k) {
      symbols.setAppend(this->outcomes(i)(k)->s);
      SymL syms_pre;
      this->outcomes(i)(k)->s->getDefiningSymbols(syms_pre, true);
      symbols.setAppend(syms_pre);
    }
  }
}


void Rule::copyBody(const Rule& other) {
  this->context = other.context;
  this->action = other.action;
  this->outcomes.resize(other.outcomes.N);
  uint i;
  FOR1D(other.outcomes, i) {
    this->outcomes(i) = other.outcomes(i);
  }
  this->probs = other.probs;
  this->noise_changes = other.noise_changes;
  this->outcome_rewards = other.outcome_rewards;
}


void Rule::write(ostream& os, bool withAction) const {
//  os << "r" << endl;
  uint i, j;
  // Default rule does not have an action specified...
  if (withAction) {
    os << "ACTION:"<<endl<<"  ";
    if (action != NULL)
      action->write(os, true);
    else
      os << "default";
    os << endl;
  }
  os << "CONTEXT:"<<endl<<"  ";
  if (context.N == 0)
    os << "--";
  FOR1D(context, i) {
    context(i)->write(os);
//     os<<context(i);
    os << " ";
  }
  os << endl;
  os << "OUTCOMES:" << endl;
  FOR1D(outcomes, i) {
    os.precision(3);
    os << "  ";
    if (probs(i) < 0.001) os << "0";
    else os << probs(i);
    os << " ";
    FOR1D(outcomes(i), j) {
      outcomes(i)(j)->write(os);
//       os<<outcomes(i)(j);
      os << " ";
    }
    if (i!=outcomes.N-1  &&  outcomes(i).N == 0)
      os << "<no-change>";
    if (i==outcomes.N-1)
//       os<<noise_changes;
      os << "<noise>";
    if (outcome_rewards.N > 0) {
      if (!TL::isZero(outcome_rewards(i)))
        os<<"  = REWARD "<<outcome_rewards(i);
    }
    os << endl;
  }
}


void Rule::getArguments(uintA& args) const {
  args.clear();
  uintA drefs;
  getDeicticRefs(drefs);
  args.setAppend(action->args);
  args.setAppend(drefs);
}


void Rule::getDeicticRefs(uintA& drefs) const {
  uintA drefs_pos, drefs_neg, drefs_nonBinary;
  getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
  drefs.clear();
  drefs.append(drefs_pos);  drefs.append(drefs_neg);
}


void Rule::getDeicticRefs(uintA& drefs_pos, uintA& drefs_neg, uintA& drefs_nonBinary) const {
  CHECK(Literal::negativeBinaryLiteralsLast(context), "negative literals should come last");
  drefs_pos.clear();  drefs_neg.clear();
  uint i, k;
  FOR1D(context, i) {
    FOR1D(context(i)->args, k) {
      uint arg = context(i)->args(k);
      if (action->args.findValue(arg) < 0) {
        if (context(i)->s->range_type == Symbol::binary) {
          if (context(i)->isNegated() && drefs_pos.findValue(arg)<0  &&  drefs_nonBinary.findValue(arg)<0)
            drefs_neg.setAppend(arg);
          else
            drefs_pos.setAppend(arg);
        }
        else {
          //cont deictic references in non-binary as neg deictic for the moment
          if (drefs_pos.findValue(arg)<0 && drefs_neg.findValue(arg))
            drefs_nonBinary.setAppend(arg);
        }
      }
    }
  }
}


void Rule::getDeicticRefsContains(boolA& containsDR) const {
  containsDR.clear();
  uintA drefs;
  getDeicticRefs(drefs);
  uint i;
  FOR1D(context, i) {
    if (numberSharedElements(drefs, action->args) > 0)
      containsDR(i) = true;
    else
      containsDR(i) = false;
  }
}


void Rule::getAbsentLiterals(LitL& literals, bool positiveOnly) const {
  literals.clear();
  uintA args;
  getArguments(args);
  LitL literals_cands;
  LitL literals_cands_pos;
  Literal::getLiterals_state(literals_cands_pos, args, 1., true);
  literals_cands.append(literals_cands_pos);
  if (!positiveOnly) {
    LitL literals_cands_neg;
    Literal::getLiterals_state(literals_cands_neg, args, 0., true);
    literals_cands.append(literals_cands_neg);
  }
  uint i, j;
  FOR1D(literals_cands, i) {
    FOR1D(context, j) {
      if (*literals_cands(i) == *context(j))
        break;
    }
    if (j==context.N)
      literals.append(literals_cands(i));
  }
}


bool Rule::usesAtom(Literal* lit) const {
  uint i, k;
  FOR1D(context, i) {
    if (context(i)->s == lit->s  &&  context(i)->args == lit->args)
      return true;
  }
  FOR1D(outcomes, i) {
    FOR1D(outcomes(i), k) {
      if (outcomes(i)(k)->s == lit->s  &&  outcomes(i)(k)->args == lit->args)
      return true;
    }
  }
  return false;
}


uint Rule::numberLiterals() const {
  uint num = 0;
  uint i;
  num += context.N;
  FOR1D(outcomes, i) {
    num += outcomes(i).N;
  }
  return num;
}


uint Rule::getNegStartPos() {
  uint i;
  FOR1D(context, i) {
    if (context(i)->isNegated())
      break;
  }
  return i;
}


uint Rule::getNonBinaryStartPos() {
  uint i;
  FOR1D(context, i) {
    if (context(i)->s->range_type != Symbol::binary || context(i)->isNegated())
      break;
  }
  return i;
}


void Rule::insertContext(Literal* literal) {
  if (context.findValue(literal) >= 0) return;
  uint i;
  uintA drefs;
  getDeicticRefs(drefs);
  uintA literal_deicticRefs;
  setSection(literal_deicticRefs, drefs, literal->args);
  context.memMove = true;
  // if positive always place before neg and non-binaries!
  if (literal->s->range_type == Symbol::binary && !literal->isNegated()) {
    // does not contain deictic refs
    if (literal_deicticRefs.N == 0) {
      context.insert(0, literal);
    }
    // does contain deictic refs
    else {
      uint pos = 0;
      // check deictic refs
      for (i = 0; i < getNonBinaryStartPos(); i++) {
        if (numberSharedElements(context(i)->args, literal_deicticRefs) > 0) {
          pos = i+1;
        }
      }
      context.insert(pos, literal);
    }
  }
  else if (literal->s->range_type == Symbol::binary) {
    uint firstNeg = getNegStartPos();
    // does not contain deictic refs
    if (literal_deicticRefs.N == 0) {
      context.insert(firstNeg, literal);
    }
    // does contain deictic refs
    else {
      uint pos = firstNeg;
      // check deictic refs
      for (i=firstNeg; i < context.N; i++) {
        if (numberSharedElements(context(i)->args, literal_deicticRefs) > 0) {
          pos = i+1;
        }
      }
      context.insert(pos, literal);
    }
  }
  else {    //non-binary
    FOR1D(context, i) {
      if (context(i)->s->range_type != Symbol::binary)
        break;
    }
    uint firstNonBinary = getNonBinaryStartPos();
    // does not contain deictic refs
    if (literal_deicticRefs.N == 0) {
      context.insert(firstNonBinary, literal);
    }
    // does contain deictic refs
    else {
      uint pos = firstNonBinary;
      // check deictic refs
      for (i = firstNonBinary; i < getNegStartPos(); i++) {
        if (numberSharedElements(context(i)->args, literal_deicticRefs) > 0) {
          pos = i+1;
        }
      }
      context.insert(pos, literal);
    }
  }
}


void Rule::cleanup() {
  uint i;
  // clean-up drefs (example: [2,3,5] -->[2,3,4])
  if (!reason::isGround(*this)) {
    uintA drefs;
    getDeicticRefs(drefs);
    if (drefs.N > 0) {
      uint max_dr = drefs.max();
      Substitution sub;
      for (i=action->s->arity; i<max_dr; i++) {
        if (drefs.findValue(i) < 0) {
          sub.addSubs(max_dr-sub.num(), i);
        }
      }
      if (!sub.empty()) {
  //       cout<<endl<<endl<<"Vorher:"<<endl; rule.writeNice();
        LitL new_context;
        sub.apply(new_context, context);
        context = new_context;
        for (i=0; i<outcomes.N-1; i++) {
          LitL new_outcome;
          sub.apply(new_outcome, outcomes(i));
          outcomes(i) = new_outcome;
        }
  //       cout<<endl<<endl<<"Nachher:"<<endl; rule.writeNice();
      }
    }
  }
  // order context
  Literal::sort(context);
  // order outcomes
  for (i=0; i<outcomes.N-1; i++) {
    Literal::sort(outcomes(i));
  }
}


void Rule::sanityCheck() { //const {
  // (1) Check no double in context
  if (context.containsDoubles()) {
    write(cerr);
    HALT("Context contains doubles");
  }
  // (2) Check no double in outcomes
  uint i, k;
  FOR1D(outcomes, i) {
    if (outcomes(i).containsDoubles()) {
      write(cerr);
      HALT("Outcome " << i << " contain doubles");
    }
  }
  // (3) Check probs sum to 1
  //if (fabs(sum(probs) - 1.0) > 0.001) {
  if (fabs(sum(probs) - 1.0) > 0.05) {
    write(cerr);
    printf("%10.10f\n", sum(probs));
    HALT("bad probs for rule");
  }
  // (4) Check no derived symbols in outcomes
  FOR1D(outcomes, i) {
    FOR1D(outcomes(i), k) {
      if (outcomes(i)(k)->s->symbol_type != Symbol::primitive) {
        write(cerr);
        HALT("Outcome " << i << " contains derived symbols.");
      }
    }
  }
  // (5) check no negated free deictic references
  uintA drefs_pos, drefs_neg, drefs_nonBinary;
  getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
  if (drefs_neg.N > 0) {
    cout<<endl<<endl<<endl;
    cout<<"FAILED SANITY CHECK  -  Negated free deictic references:   "<<drefs_neg<<endl;
    cout<<"RULE:  "<<endl<<*this;
    cout<<"Remove this sanity check if you want to allow for negated free deictic references."<<endl;
    cerr<<"Remove this sanity check if you want to allow for negated free deictic references."<<endl;
    HALT("FAILED SANITY CHECK  -  Negated free deictic references:   "<<drefs_neg);
  }
  if (drefs_nonBinary.N > 0) {
    cout<<endl<<endl<<endl;
    cout<<"FAILED SANITY CHECK  -  free non-binary deictic references:   "<<drefs_nonBinary<<endl;
    cout<<"RULE:  "<<endl<<*this;
    cout<<"Remove this sanity check if you want to allow for free non-binary references."<<endl;
    cerr<<"Remove this sanity check if you want to allow for free non-binary references."<<endl;
//    HALT("FAILED SANITY CHECK  -  free non-binary deictic references:   "<< drefs_nonBinary);
  }
}


bool Rule::existsInOutcome(Symbol *s, uintA args) {
  uint i, j;
  FOR1D(outcomes, i) {
    FOR1D(outcomes(i), j) {
    if (outcomes(i)(j)->s == s && outcomes(i)(j)->args == args)
      return true;
    }
  }

  return false;
}


Rule* Rule::read(istream& in) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"Rule::read [START]"<<endl;}
  
  //CHECK(in.is_open(), "Input stream ain't open!");
  Rule* r = new Rule;
  
  MT::String line, untrimmed_line;
  
  // Action
  line.read(in, NULL, "\n"); // ACTION:
  if (DEBUG>1) PRINT(line);
  CHECK(line(0)=='A',"bad action (expecting action first);  line="<<line);   CHECK(line.N < 10, "bad action: too short");
  line.read(in, NULL, "\n");
  if (DEBUG>1) PRINT(line);
  LitL actions_wrapper;
  Literal::get(actions_wrapper, line);
  r->action = actions_wrapper(0);
  if (DEBUG>1) PRINT(*r->action);
  
  // Context
  untrimmed_line.read(in, NULL, "\n"); // CONTEXT:
  line.read(untrimmed_line, NULL, " \n\t");
  if (DEBUG>1) PRINT(line);
  CHECK((line(0)=='C' && line(1)=='O'), "bad context:  line="<<line);
  CHECK(line.N < 10, "bad context due to length (line.N="<<line.N<<"; line=\""<<line<<"\")");
  line.read(in, NULL, "\n");
  if (DEBUG>1) PRINT(line);
  if (line(0) != '-'   ||   (line.N > 1  &&  line(1) != '-'  &&  line(1) != ' '))
    Literal::get(r->context, line);
  
  // Outcomes
  if (DEBUG>0) {cout<<"Reading the outcomes:"<<endl;}
  line.read(in, NULL, "\n"); // OUTCOMES:
  CHECK((line(0)=='O' && line(1)=='U'),"bad outcomes:  "<<line);   CHECK(line.N < 11, "bad outcomes due to length");
  MT::Array< LitL > outcomes;
  bool noise_outcome_has_been_read = false;
  while ( MT::peerNextChar(in) != ' '  &&  MT::peerNextChar(in) != 'P'  &&   MT::peerNextChar(in) != 'C'
                                       &&  MT::peerNextChar(in) != 'R'&&  MT::peerNextChar(in) != 'A') {
    line.read(in, NULL, "\n");
    if (DEBUG>1) PRINT(line);
    if (line.N<2) {
      cout<<"bad line: "<<endl;
      PRINT(line);
      HALT("");
    }
    double prob;
    line >> prob;
    r->probs.append(prob);
    MT::skip(line);
    LitL outcome;
    MT::String rest_line;
    rest_line.read(line);
    if (rest_line(0)=='<' && rest_line(3)=='i') {  // <noise> 
      noise_outcome_has_been_read = true;
      if (DEBUG>1) {cout<<"<noise> outcome detected"<<endl;}
    }
    else if (rest_line(0)=='<' && rest_line(3)=='-') {   // <no-change>
      if (DEBUG>1) {cout<<"<no-change> outcome detected"<<endl;}
    }
    else  // if neither noise outcome nor no-change
      Literal::get(outcome, rest_line);
    r->outcomes.append(outcome);
    if (MT::skip(in) == -1)
      break;
    //  mglw. noch reward fuers outcome einlesen
    if (MT::peerNextChar(in) == 'R') {
      while (r->outcome_rewards.N < r->outcomes.N-1) {  // possibly fill for previous unrewarded outcomes
        r->outcome_rewards.append(0.);
      }
      MT::skipUntil(in, " ");
      double outcome_reward;
      in >> outcome_reward;
      r->outcome_rewards.append(outcome_reward);
      if (MT::skip(in) == -1)
        break;
    }
  }
  if (!noise_outcome_has_been_read) {
    LitL dummy_noise_outcome;
    r->outcomes.append(dummy_noise_outcome);
    r->noise_changes = 0.;
    r->probs.append(1.0 - sum(r->probs));
    if (r->outcome_rewards.N > 0)
      r->outcome_rewards.append(0.);
  }
  if (TL::isZero(r->probs.last())) {
    r->probs(0) -= 0.000001;
    r->probs.last() = 0.000001;  // hack to avoid zero-prob noise-outcome
  }
  
  if (DEBUG>0) {
    cout<<"New rule:"<<endl;
    r->write(cout);
  }
  
  if (DEBUG>0) {cout<<"Rule::read [END]"<<endl;}
  return r;
}


Rule* Rule::generateDefaultRule(double noiseProb, double minProb, double change) {
  CHECK(!isnan(noiseProb), "noiseProb is nan!");
  CHECK(!isnan(minProb), "minProb is nan!");
  CHECK(!isnan(change), "change is nan!");

  Rule* newDefaultRule = new Rule;
  uintA args_empty;
  newDefaultRule->action = Literal::get(Symbol::getDefaultAction(), args_empty, 1.);
  LitL sameOutcome;
  newDefaultRule->outcomes.append(sameOutcome);
  LitL noiseOutcome;
  newDefaultRule->outcomes.append(noiseOutcome);
  if (noiseProb < minProb)
    noiseProb = minProb;
  else if (1.-noiseProb < minProb)
    noiseProb = 1.-minProb;
  newDefaultRule->probs.append(1-noiseProb);
  newDefaultRule->probs.append(noiseProb);
  newDefaultRule->noise_changes = change;
  return newDefaultRule;
}


bool Rule::isDefaultRule(const Rule* rule) {
  return     rule->context.N == 0
          && rule->action == Literal::getLiteral_default_action();
}


Rule* Rule::getDoNothingRule() {
  Rule* r = new Rule;
  r->action = Literal::getLiteral_doNothing();
  r->outcomes.resize(2);
  r->probs.resize(2);
  r->probs(0) = 1.0;
  r->probs(1) = 0.0;
  r->noise_changes = 0.0;
  return r;
}




/************************************************
 * 
 *     RuleSet
 * 
 ************************************************/

RuleSet::RuleSet() {
  ra.memMove = true;
}


RuleSet::~RuleSet() {
  clear();
}


int RuleSet::findValue(Rule* r) {
  return ra.findValue(r);
}


void RuleSet::clear() {
  uint i;
  FOR1D(ra, i) {
    TL::del(ra(i));
  }
  ra.clear();
}


void RuleSet::append(Rule* r) {
  ra.append(TL::getRef(r));
}


uint RuleSet::num() const {return ra.N;}


Rule* RuleSet::elem(uint i) const {return ra(i);}


void RuleSet::remove(uint i) {
  TL::del(ra(i));
  ra.remove(i);
}


void RuleSet::overwrite(uint i, Rule* r) {
  TL::del(ra(i));
  ra(i) = TL::getRef(r);
}

RuleSet& RuleSet::operator=(const RuleSet& rs) {
  this->clear();
  uint i;
  FOR1D(rs.ra, i) {
    this->append(rs.ra(i));
  }
  return *this;
}


void RuleSet::sort() {
  MT::Array< Rule* > ra_sorted;
  uint i, k;
  // default rule to the beginning!
  FOR1D(ra, i) {
    if (Rule::isDefaultRule(ra.elem(i))) {
      ra_sorted.append(ra.elem(i));
      break;
    }
  }
  SymL symbols_action;
  Symbol::get_action(symbols_action);
  FOR1D(symbols_action, i) {
    if (symbols_action(i)->name == "default") continue;
    FOR1D(ra, k) {
      if (ra(k)->action->s == symbols_action(i))
        ra_sorted.append(ra(k));
    }
  }
  ra = ra_sorted;
}

void RuleSet::sort_using_args() {
  uintA new_order;
  sort_using_args(new_order);
}


void RuleSet::sort_using_args(uintA& new_order) {
  new_order.clear();

  uint max_arity = 0;
  uint r2;
  FOR1D(this->ra, r2) {
    if (this->ra(r2)->action->s->arity > max_arity)
      max_arity = this->ra(r2)->action->s->arity;
  }
  
  SymL symbols_action;
  Symbol::get_action(symbols_action);
  
  if (max_arity <= 2) {
    MT::Array< Rule* > ra_sorted;
    uintA action_ids;
    uint r, i;
    uint id;
    
    FOR1D(this->ra, r) {
      int action_idx = symbols_action.findValue(this->ra(r)->action->s);
      CHECK(action_idx >= 0, "action in rule is unknown");
      if (this->ra(r)->action->s->name == "default")
        continue;
      if (this->ra(r)->action->s->arity == 0)
        id = action_idx * 100;
      else if (this->ra(r)->action->s->arity == 1)
        id = action_idx * 100  +  this->ra(r)->action->args(0);
      else if (this->ra(r)->action->s->arity == 2)
        id = action_idx * 10000  +  100 * this->ra(r)->action->args(0) + this->ra(r)->action->args(1);
      else
        HALT("");
      if (action_ids.findValue(id) >= 0)
        continue;
      FOR1D(action_ids, i) {
        if (action_ids(i) > id) {
          action_ids.insert(i, id);
          break;
        }
      }
      if (action_ids.N == i)
        action_ids.append(id);
    }
    
    // DEFAULT RULE always first rule
    int action_idx = symbols_action.findValue(Symbol::getDefaultAction());
    uint DEFAULT_ACTION_ID = action_idx * 100;
    action_ids.insert(0, DEFAULT_ACTION_ID);
    
    FOR1D(action_ids, i) {
      FOR1D(this->ra, r) {
        int action_idx = symbols_action.findValue(this->ra(r)->action->s);
        if (this->ra(r)->action->s->arity == 0) {
          if (action_idx * 100 == action_ids(i)) {
            ra_sorted.append(this->ra(r));
            new_order.append(r);
          }

        }
        else if (this->ra(r)->action->s->arity == 1) {
          if (action_idx * 100  +  this->ra(r)->action->args(0) == action_ids(i)) {
            ra_sorted.append(this->ra(r));
            new_order.append(r);
          }
        }
        else if (this->ra(r)->action->s->arity == 2) {
          if (action_idx * 10000  +  100 * this->ra(r)->action->args(0) + this->ra(r)->action->args(1) == action_ids(i)) {
            ra_sorted.append(this->ra(r));
            new_order.append(r);
          }
        }
        else
          HALT("");
      }
    }
    this->ra = ra_sorted;
  }
  // max_arity > 2
  else {
    relational::write(*this, "ground_rules.dat.unsorted_backup");
    
    MT::Array< Rule* > ra_sorted;
    ra_sorted.memMove = true;
    new_order.memMove = true;
    uint r, r2, d;
    FOR1D(this->ra, r) {
      bool inserted = false;
      FOR1D(ra_sorted, r2) {
        int action_idx1 = symbols_action.findValue(this->ra(r)->action->s);
        CHECK(action_idx1 >= 0, "action in rule is unknown");
        int action_idx2 = symbols_action.findValue(ra_sorted(r2)->action->s);
        CHECK(action_idx2 >= 0, "action in rule is unknown");
        // put before action-predicate with higher id
        if (action_idx1 > action_idx2) {
          ra_sorted.insert(r2, this->ra(r));
          new_order.insert(r2, r);
          inserted = true;
        }
        // same action-predicate
        else if (action_idx1 == action_idx2) {
          FOR1D(this->ra(r)->action->args, d) {
            if (this->ra(r)->action->args(d) < ra_sorted(r2)->action->args(d)) {
              ra_sorted.insert(r2, this->ra(r));
              new_order.insert(r2, r);
              inserted = true;
              break;
            }
          }
        }
        if (inserted)
          break;
      }
      if (!inserted) {
        ra_sorted.append(this->ra(r));
        new_order.append(r);
      }

//       FOR1D(ra_sorted, r3) {
//         ra_sorted(r3)->action->write(); cout<<"  ";
//       }
//       cout<<endl;
    }
    
//     FOR1D(ra_sorted, r3) {
//       ra_sorted(r3)->action->write(); cout<<"  ";
//     }
//     cout<<endl;
    
    CHECK(ra.N == ra_sorted.N, "Some strange rule-sorting mistake.");
    this->ra = ra_sorted;
  }
}


void RuleSet::write(ostream& os) const {
  uint i;
  os<<"#rules = "<<ra.N<<endl;
  FOR1D(ra, i) {
    os << "# Rule #" << i << "  ("<< (i+1) << " out of " << num() << ")" <<endl;
    os << *elem(i) << endl;
  }
}


void RuleSet::sanityCheck() const {
  uint i;
  FOR1D(ra, i) {
    ra(i)->sanityCheck();
  }
}

bool RuleSet::hasDefaultRule() const {
  uint i;
  FOR1D(ra, i) {
    if (Rule::isDefaultRule(ra(i)))
      return true;
  }
  return false;
}

void RuleSet::changingSymbols(SymL& symbols) const {
  symbols.clear();
  uint i, k, l;
  FOR1D(ra, i) {
    FOR1D(ra(i)->outcomes, k) {
      FOR1D(ra(i)->outcomes(k), l) {
        symbols.setAppend(ra(i)->outcomes(k)(l)->s);
      }
    }
  }
}


void RuleSet::removeDoubleLiterals() {
  int DEBUG = 0;
  uint i, k, l, m;
  Rule* r;
  FOR1D(ra, i) {
    r = ra.elem(i);
    if (DEBUG>0) {
      cout<<"BEFORE:";
      r->write(cout);
    }
    // (1) context
    for (k=0; ; k++) {
      if (r->context.N <= k+1)
        break;
      l=k+1;
      while (l < r->context.N) {
        if (r->context(k) == r->context(l)) {
          r->context.memMove = true;
          r->context.remove(l);
          if (DEBUG>0) {cout<<"HUI"<<endl;}
        }
        else {
          l++;
        }
      }
    }
    // (2) Outcomes
    FOR1D(r->outcomes, m) {
      for (k=0; ; k++) {
        if (r->outcomes(m).N <= k+1)
          break;
        l=k+1;
        while (l < r->outcomes(m).N) {
          if (r->outcomes(m)(k) == r->outcomes(m)(l)) {
            r->outcomes(m).memMove = true;
            r->outcomes(m).remove(l);
            if (DEBUG>0) {cout<<"HUI"<<endl;}
          }
          else {
            l++;
          }
        }
      }
    }
    if (DEBUG>0) {
      cout<<"AFTER:";
      r->write(cout);
    }
  }
}


void RuleSet::removeNonChangingSymbols() {
  SymL changing_symbols;
  changingSymbols(changing_symbols);
  uint i, k, l;
  FOR1D(ra, i) {
    Rule* r_ground = ra(i);
    r_ground->context.memMove = true;
    FOR1D_DOWN(r_ground->context, k) {
      NIY;
//       if (changingIds_preds.findValue(r_ground->context(k)->atom->s->id) < 0)
//         r_ground->context.remove(k);
    }
    FOR1D(r_ground->outcomes, l) {
      r_ground->outcomes.memMove = true;
      FOR1D_DOWN(r_ground->outcomes(l), k) {
        NIY;
//         if (changingIds_preds.findValue(r_ground->outcomes(l)(k)->atom->s->id) < 0)
//           r_ground->outcomes(l).remove(k);
      }
    }
  }
}

void RuleSet::read(const char* filename, RuleSet& rules) {
  uint DEBUG = 0;
  ifstream in(filename);
  if (DEBUG>0) {PRINT(filename);}
//   PRINT(filename);
  if (!in.is_open()) {
    cerr<<"Rule-file " << filename << " can't be opened!"<<endl;
    HALT("");
  }
  read(in, rules);
}

void RuleSet::read(istream& in, RuleSet& rules, bool sort_rules) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"read [START]"<<endl;}
  rules.clear();
  
  // read other rules
  while (MT::skip(in) != -1) {
    rules.append(Rule::read(in));
    if (DEBUG>0) {rules.elem(rules.num()-1)->write();}
  }
  // default rule
  if (!Rule::isDefaultRule(rules.elem(0))) {
    rules.append(Rule::generateDefaultRule());
  }
  if (sort_rules)
    rules.sort();
  rules.sanityCheck();
  if (DEBUG>0) {cout<<"read [END]"<<endl;}
}


#ifndef FAST__RULES_REASONING
// Assumption: Rule arguments have to be different vom Deictic References
void RuleSet::ground(RuleSet& rules_ground, const RuleSet& rules_abstract, const uintA constants) {
  rules_ground.clear();
  uint r, c1, c2, i;
  FOR1D_(rules_abstract, r) {
    uintA args, drefs;
    calcTerms(*rules_abstract.elem(r), args);
    boolA inNegatedLiteralsOnly;
    calcDeicticRefs(*r_abs, drefs, inNegatedLiteralsOnly);
    if (inNegatedLiteralsOnly.findValue(true) >= 0) {HALT("disallow negated free DRs");}
    setMinus(args, drefs);
    MT::Array< uintA > combos_args;
    TL::allPermutations(combos_args, constants, args.N, true, true);
    FOR1D(combos_args, c1) {
      uintA dref_constants = constants;
      setMinus(dref_constants, combos_args(c1));
      MT::Array< uintA > combos_drefs;
      TL::allPermutations(combos_drefs, dref_constants, drefs.N, true, true);
      FOR1D(combos_drefs, c2) {
        Substitution sub;
        FOR1D(args, i) {
          sub.addSubs(args(i), combos_args(c1)(i));
        }
        FOR1D(drefs, i) {
          sub.addSubs(drefs(i), combos_drefs(c2)(i));
        }
        Rule* ground_rule = ground(rules_abstract.elem(r), &sub);
        rules_ground.append(ground_rule);
      }
    }
  }
  rules_ground.sort_using_args();
}
#endif
#ifdef FAST__RULES_REASONING
// Assumption: Rule arguments have to be different vom Deictic References
void RuleSet::ground(RuleSet& rules_ground, const RuleSet& rules_abstract, const uintA& constants) {
  uint DEBUG = 0;
  rules_ground.clear();
  uint r, c1, c2, i, k;
  uint changingDR;
  boolA copyTable_context(15); // assume max 15 context
  boolA copyTable_outcomes(10, 15); // assume max 10 outcomes with max 15 changes each
//   uint NUM_REJECTED = 0;
  
  FOR1D_(rules_abstract, r) {
    Rule* r_abs = rules_abstract.elem(r);
    uintA args, drefs_pos, drefs_neg, drefs_nonBinary;
    r_abs->getArguments(args);
    boolA inNegatedLiteralsOnly;
    r_abs->getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
    if (drefs_neg.N > 0 || drefs_nonBinary.N > 0) {HALT("negated or non-Binary free DRs are not allowed");}
    setMinus(args, drefs_pos);
    if (DEBUG>0) {
      cout<<endl<<"******* New abstract rule:"<<endl; r_abs->write(cout);
      PRINT(args);
      PRINT(drefs_pos);
      PRINT(constants);
    }
    MT::Array< uintA > combos_args;
    TL::allPermutations(combos_args, constants, args.N, false, true);
    if (drefs_pos.N > 0) {
      changingDR = drefs_pos(0);
      FOR1D(r_abs->context, i) {
        uintA context_lit_args;
        if (r_abs->context(i)->args.findValue(changingDR) < 0)
          copyTable_context(i) = true;
        else
          copyTable_context(i) = false;
      }
      FOR1D(r_abs->outcomes, i) {
        FOR1D(r_abs->outcomes(i), k) {
          if (r_abs->outcomes(i)(k)->args.findValue(changingDR) < 0)
            copyTable_outcomes(i,k) = true;
          else
            copyTable_outcomes(i,k) = false;
        }
      }
    }
    if (DEBUG>0) {PRINT(combos_args);}
    FOR1D(combos_args, c1) {
      if (DEBUG>1) {PRINT(combos_args(c1));}
      // Check types [START]
      FOR1D(r_abs->action->s->arg_types, i) {
        ArgumentType* arg_type = r_abs->action->s->arg_types(i);
        ArgumentType* object_type = reason::getArgumentTypeOfObject(combos_args(c1)(i));
        if (DEBUG>1) {
          r_abs->action->write(); cout<<endl;
          cout<<"i="<<i<<"  arg "; arg_type->write();
          cout<<"  with   const ";  object_type->write();
          cout<<"  [combos_args(c1)(i)="<<combos_args(c1)(i)<<"]"<<endl;
        }
        if (!arg_type->subsumes(*object_type)) {
//           cout<<"not subsumed!"<<endl;
          break;
        }
      }
      // Abort if types are inadequate
      if (r_abs->action->s->arg_types.N > 0   &&  i < combos_args(c1).N) {
//         cout<<"HUHU KILLING me softly"<<endl;
        if (DEBUG>1) {cout<<"Inadequate types"<<endl;}
        continue;
      }
      // Check types [END]
      
      uintA dref_constants = constants;
      setMinus(dref_constants, combos_args(c1));
      MT::Array< uintA > combos_drefs;
      TL::allPermutations(combos_drefs, dref_constants, drefs_pos.N, true, true);
      if (DEBUG>0) PRINT(combos_drefs);
      Rule* r_last = NULL;
      FOR1D(combos_drefs, c2) {
        Substitution sub;
        FOR1D(args, i) {
          sub.addSubs(args(i), combos_args(c1)(i));
        }
        FOR1D(drefs_pos, i) {
          sub.addSubs(drefs_pos(i), combos_drefs(c2)(i));
        }
        if (DEBUG>0) {cout<<endl<<"Grounding with: ";sub.write(cout);}
        Rule* r_ground;
        // ACHTUNG! combos_drefs change from left to right 
        //          -->  most recent change = combos_drefs(c2)(0)
        // If there are deictic references, we'll copy from last rule all PTs
        // that do not include the latest DR (the one to change in the new rule).
        if (c2>0  &&  r_last != NULL  &&
            ( drefs_pos.N==1 // safe if only 1 DR
            ||
            combos_drefs(c2)(1) == combos_drefs(c2-1)(1))
            // ensure that second-order changing DR is the same
           ) {
          CHECK(combos_drefs(c2)(0) != combos_drefs(c2-1)(0), "");
          if (DEBUG>0) cout<<"PT-copying"<<endl;
          r_ground = new Rule;
          r_ground->copyBody(*r_abs);
          
          FOR1D(r_abs->context, i) {
            if (DEBUG>0) r_ground->context(i)->write(cout);
            if (copyTable_context(i)) {
              r_ground->context(i) = r_last->context(i);
              if (DEBUG>0) cout<<" -> copied   ";
            }
            else {
              r_ground->context(i) = sub.apply(r_abs->context(i));
              if (DEBUG>0) cout<<" -> subsed   ";
            }
          }
          
          r_ground->action = r_last->action;
  
          r_ground->outcomes = r_abs->outcomes.resize(r_abs->outcomes.N);
          FOR1D(r_abs->outcomes, i) {
            FOR1D(r_abs->outcomes(i), k) {
              if (DEBUG>0) r_abs->outcomes(i)(k)->write(cout);
              if (copyTable_outcomes(i,k)) {
                r_ground->outcomes(i)(k) = r_last->outcomes(i)(k);
                if (DEBUG>0) cout<<" -> copied   ";
              }
              else {
                r_ground->outcomes(i)(k) = sub.apply(r_abs->outcomes(i)(k));
                if (DEBUG>0) cout<<" -> subsed   ";
              }
            }
          }
          r_ground->probs = r_abs->probs;
          r_ground->noise_changes = r_abs->noise_changes;
        }
        else {
          r_ground = sub.apply(*rules_abstract.elem(r));
        }
        if (DEBUG>0) {cout<<endl<<"Ground rule:"<<endl;r_ground->write(cout);}        
        rules_ground.append(r_ground);
        r_last = r_ground;
      }
    }
  }
//   PRINT(NUM_REJECTED);
  rules_ground.sort_using_args();
}
#endif

void RuleSet::ground_with_filtering(RuleSet& rules_ground,
                                    const RuleSet& rules_fullyAbstract, const uintA& constants,
                                    const SymbolicState& state, bool delete_nonchanging_concepts) {
  uintA abstract_rules_idx;
  MT::Array<uintA> abstract2GroundMap;
  ground_with_filtering(rules_ground, abstract_rules_idx, abstract2GroundMap, rules_fullyAbstract, constants,
                        state, delete_nonchanging_concepts);

}

void RuleSet::ground_with_filtering(RuleSet& rules_ground,
                                    uintA& abstract_rules_idx,
                                    MT::Array<uintA>& abstract2GroundMap,
                                    const RuleSet& rules_fullyAbstract, const uintA& constants, const SymbolicState& state, bool delete_nonchanging_concepts) {
  uint DEBUG = 0;
  if (DEBUG>0) {cout<<"RuleSet::ground_with_filtering [START]"<<endl;}
  if (DEBUG>0) {
    cout<<endl<<endl<<endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout<< "**************  ATTENTION  ************"<<endl;
    cout << "Assumes that DRs lists change from left to right!!!"<<endl <<endl <<endl <<endl;
  }

  abstract_rules_idx.clear();
  abstract2GroundMap.resize(rules_fullyAbstract.num());

  MT::Array<uintA> abstract2GroundMapUnsorted;
  abstract2GroundMapUnsorted.resize(rules_fullyAbstract.num());
  uintA abstract_rules_idxUnsorted;

  RuleSet rules_abstract;
  groundFunctionVars(rules_fullyAbstract, rules_abstract);
  if (DEBUG > 0) { cout << "Rules with grounded function values:" << endl; rules_abstract.write(); cout << endl; }
  
  // (1) Determine never changing predicate instances and function values
  
  uint i, k;
  SymL symbols_state;
  Symbol::get_state(symbols_state);
  SymL changing_symbols;
  rules_abstract.changingSymbols(changing_symbols);
  SymL nonchanging_symbols;
  FOR1D(symbols_state, i) {
    if (changing_symbols.findValue(symbols_state(i)) < 0  &&  symbols_state(i)->symbol_type == Symbol::primitive)
      nonchanging_symbols.append(symbols_state(i));
  }
  
  LitL nonchanging_lits;
  FOR1D(nonchanging_symbols, i) {
    Symbol* s = nonchanging_symbols(i);
    MT::Array< uintA > combos;
    TL::allPermutations(combos, constants, s->arity, true, true);
    FOR1D(combos, k) {
      if (s->range_type == Symbol::binary) {
        Literal *lit = Literal::get(s, combos(k), 1);
        if (reason::holds(state.lits, lit))
          nonchanging_lits.append(lit);
      }
      else {
        for (uint j = 0; j < s->range.N; j++) {
          Literal *lit = Literal::get(s, combos(k), s->range(j));
        if (reason::holds(state.lits, lit))
          nonchanging_lits.append(lit);
        }
      }

    }
  }
  
  if (DEBUG>0) {
    cout<<state<<endl;
    cout<<"Never changing symbols: "<<nonchanging_symbols<<endl;
    cout<<"Never changing literals: "<<nonchanging_lits<<endl;
    cout<<endl;
  }
  
  
  // (2) Ground rules with pruning
  rules_ground.clear();
  uint r, c1, c2;
  uint latestDR;
  uint MAX_NUM_CONTEXTLITERALS = 20;
  uint MAX_NUM_OUTCOME_PROPERTIES = 20;
  boolA context_containsLatestDR(MAX_NUM_CONTEXTLITERALS); // assume max 15 context
  boolA context_containsDR(MAX_NUM_CONTEXTLITERALS); // assume max 15 context
  boolA outcomes_containsLatestDR(10, MAX_NUM_OUTCOME_PROPERTIES); // assume max 10 outcomes with max 15 changes each
  Rule* r_last = NULL;
  bool valid_action_args;
  bool valid_nonlatest_drefs = true;
  bool valid_all;
  
  uint NUM_REJECTED = 0;
  
  // ABSTRACT RULES
  FOR1D_(rules_abstract, r) {
    Rule* r_abs = rules_abstract.elem(r);
    uintA args, drefs_pos, drefs_neg, drefs_nonBinary;
    r_abs->getArguments(args);
    boolA inNegatedLiteralsOnly;
    r_abs->getDeicticRefs(drefs_pos, drefs_neg, drefs_nonBinary);
    if (drefs_neg.N > 0 || drefs_nonBinary.N > 0) {
      cerr<<endl<<endl<<endl<<"Bad abstract rule:"<<endl;  r_abs->write(cerr);  HALT("disallow negated or non-Binary free DRs!");
    }
    setMinus(args, drefs_pos);
    
    // Prepare copying and pruning:  we copy everything except PIs that contain the first DR
    context_containsDR.setUni(false);
    context_containsLatestDR.setUni(false);
    outcomes_containsLatestDR.setUni(false);
    if (drefs_pos.N > 0) {
      latestDR = drefs_pos(0);
      FOR1D(r_abs->context, i) {
        if (r_abs->context(i)->args.findValue(latestDR) < 0)
          context_containsLatestDR(i) = false;
        else
          context_containsLatestDR(i) = true;
        if (numberSharedElements(r_abs->context(i)->args, drefs_pos) == 0)
          context_containsDR(i) = false;
        else
          context_containsDR(i) = true;
      }
      FOR1D(r_abs->outcomes, i) {
        FOR1D(r_abs->outcomes(i), k) {
          if (r_abs->outcomes(i)(k)->args.findValue(latestDR) < 0)
            outcomes_containsLatestDR(i,k) = false;
          else
            outcomes_containsLatestDR(i,k) = true;
        }
      }
    }
    
    if (DEBUG>1) {
      cout<<endl<<"***************************************************"<<endl;
      cout<<endl<<"New abstract rule:"<<endl; r_abs->write(cout);
      PRINT(args);
      PRINT(drefs_pos);
      PRINT(constants);
      cout<<"context_containsDR:"<<endl;
      for (i=0; i<r_abs->context.N; i++) {
        cout<<context_containsDR(i)<<" ";
      }
      cout<<endl;
      PRINT(context_containsDR);
      
      cout<<"context_containsLatestDR:"<<endl;
      for (i=0; i<r_abs->context.N; i++) {
        cout<<context_containsLatestDR(i)<<" ";
      }
      cout<<endl;
      
      cout<<"outcomes_containsLatestDR:"<<endl;
      FOR1D(r_abs->outcomes, i) {
        cout<<i<<":";
        FOR1D(r_abs->outcomes(i), k) {
          cout<<" "<<outcomes_containsLatestDR(i, k);
        }
        cout<<endl;
      }
      cout<<endl;
    }
    
    MT::Array< uintA > combos_args;
//     TL::allPermutations(combos_args, constants, args.N, true, true);
    TL::allPermutations(combos_args, constants, args.N, false, true);

    // Level 1:  ARGUMENT combos
//     PRINT(combos_args.N);
//     cout<<"Current #ground rules = "<<rules_ground.num()<<endl;
    FOR1D(combos_args, c1) {
      if (c1%1000 == 0) {
        cout<<"."<<std::flush;
      }
      
      // Check types
      FOR1D(r_abs->action->s->arg_types, i) {
        ArgumentType* arg_type = r_abs->action->s->arg_types(i);
        ArgumentType* object_type = reason::getArgumentTypeOfObject(combos_args(c1)(i));
//         cout<<"SDHF: ";
//         r_abs->action->writeNice();  cout<<"  "; arg_type->writeNice();  cout<<endl;
//         PRINT(combos_args(c1)(i));
//         object_type->writeNice();  cout<<endl;
        if (!arg_type->subsumes(*object_type))
          break;
      }
      // Abort if types are inadequate;  but only if types are specified at all
      if (r_abs->action->s->arg_types.N > 0   &&  i < combos_args(c1).N) {
//         cout<<"HUHU KILLING me softly"<<endl;
        continue;
      }
      
      uintA dref_constants = constants;
      setMinus(dref_constants, combos_args(c1)); // Assumption: DRef different from arguments -- but not necessarily in rule-learning!
      MT::Array< uintA > combos_drefs;
      //TL::allPermutations(combos_drefs, dref_constants, drefs_pos.N, true, true);
      try {
        TL::allPermutations(combos_drefs, dref_constants, drefs_pos.N, false, true); // disallow DRs to repeat
      } catch (...) {
        std::cerr << "error computing TL::allPermutations(combos_drefs, dref_constants, drefs_pos.N, false, true)" << std::endl;
        std::cerr << " RULE: " << r << std::endl;
        rules_abstract.elem(r)->write(std::cerr); std::cerr << std::endl;
        PRINTOUT(combos_drefs, cerr);
        PRINTOUT(dref_constants, cerr);
        PRINTOUT(drefs_pos.N, cerr);
      }

      if (DEBUG>2) {cout<<"------- News action args -------"<<endl; PRINT(combos_args);  PRINT(combos_drefs);}
      valid_action_args = true;
      r_last = NULL;
      // Level 2:  DR combos
      FOR1D(combos_drefs, c2) {
        if (c2>0  &&  drefs_pos.N > 1  &&  combos_drefs(c2)(1) == combos_drefs(c2-1)(1)   &&  !valid_nonlatest_drefs) {
          if (DEBUG>2) cout<<"Omitting combos_drefs="<<combos_drefs(c2)<<endl;
          NUM_REJECTED++;
          continue;
        }
        else
          valid_nonlatest_drefs = true;
        // If we have new non_latest_DR --> cannot copy anymore.
        if (c2>0  &&  drefs_pos.N > 1  &&  combos_drefs(c2)(1) != combos_drefs(c2-1)(1))
          r_last = NULL;

        Substitution sub;
        FOR1D(args, i) {
          sub.addSubs(args(i), combos_args(c1)(i));
        }
        FOR1D(drefs_pos, i) {
          sub.addSubs(drefs_pos(i), combos_drefs(c2)(i));
        }
        if (DEBUG>2) {cout<<endl<<"++++++++++ Next try ++++++++++"<<endl<<"Grounding with: ";sub.write(cout);cout<<endl;}
        Rule* r_ground;
        // (I)   RULE CREATION BY COPYING  
        //       If there are deictic references, we'll copy from last rule all PTs
        //       that do not include the latest DR (the one to change in the new rule),
        //       which is combos_drefs(c2)(0).
        //       Attention:  combos_drefs change from left to right 
        //          -->  most recent change = combos_drefs(c2)(0)
        if (r_last != NULL) {
          CHECK(combos_drefs(c2)(0) != combos_drefs(c2-1)(0), "");
          if (DEBUG>2) cout<<"--- Creation by copying"<<endl;
          r_ground = new Rule;
          r_ground->copyBody(*r_abs);
          
          valid_all = true;
          FOR1D(r_abs->context, i) {
            if (DEBUG>3) r_ground->context(i)->write(cout);
            if (context_containsLatestDR(i)) {
              r_ground->context(i) = sub.apply(r_abs->context(i));
              if (DEBUG>3) cout<<" -> subsed   ";
              if (nonchanging_symbols.findValue(r_ground->context(i)->s) >= 0) {
                if (!reason::holds(nonchanging_lits, r_ground->context(i))) {
                  valid_all = false;
                  if (DEBUG>3) {cout<<" invalid! ";}
                }
              }
               if (!r_ground->context(i)->typeCheck()) valid_all = false; //type check
            }
            else {
              r_ground->context(i) = r_last->context(i);
              if (DEBUG>3) cout<<" -> copied   ";
            }
          }
          
          if (!valid_all) {
            delete r_ground;
            NUM_REJECTED++;
            if (DEBUG>2) cout<<" --> REJECTED"<<endl;
          }
          else {
            r_ground->action = r_last->action;
    
            r_ground->outcomes = r_abs->outcomes.resize(r_abs->outcomes.N);
            FOR1D(r_abs->outcomes, i) {
              FOR1D(r_abs->outcomes(i), k) {
                if (DEBUG>3) r_abs->outcomes(i)(k)->write(cout);
                if (outcomes_containsLatestDR(i,k)) {
                  r_ground->outcomes(i)(k) = sub.apply(r_abs->outcomes(i)(k));
                  if (DEBUG>3) cout<<" -> subsed   ";
                }
                else {
                  r_ground->outcomes(i)(k) = r_last->outcomes(i)(k);
                  if (DEBUG>3) cout<<" -> copied   ";
                }
              }
            }
            
            r_ground->noise_changes = r_abs->noise_changes;
            rules_ground.append(r_ground);
            abstract2GroundMapUnsorted(r).append(rules_ground.num()-1);
            abstract_rules_idxUnsorted.append(r);

            r_last = r_ground;
            if (DEBUG>2) {cout<<endl<<"Candidate ground rule:"<<endl;r_ground->write(cout);}
            if (DEBUG>2) cout<<" --> ACCEPTED"<<endl;
          }
        }
        // (II)   RULE CREATION BY FULL GROUNDING
        else {
          if (DEBUG>2) cout<<"--- Creation by full grounding"<<endl;
          r_ground = sub.apply(*rules_abstract.elem(r));
          // FILTERING
          // Check whether context 
          //   (A)  that contain only action-arguments
          //   (B)  that contain not the latest DR
          //  are invalid.
          valid_all = true;
          FOR1D(r_ground->context, i) {
            // (i) Is non-changing predicate?
            if (nonchanging_symbols.findValue(r_ground->context(i)->s) >= 0) {
                // (ii) Does context hold?
                if (!reason::holds(nonchanging_lits, r_ground->context(i))) {
                  valid_all= false;
                  if (DEBUG>3) {cout<<" invalid!  !reason::holds("<<nonchanging_lits<<", "<<*r_ground->context(i)<<")";}
                  // (iii A) Only action arguments?
                  if (!context_containsDR(i)) {
                    valid_action_args = false;
                    if (DEBUG>2) {cout<<"valid_action_args=false:  literal "<<i<<" "<<r_ground->context(i)<<endl;}
                  }
                  // (iii B) Only non_latest_DR arguments?
                  else if (!context_containsLatestDR(i)) {
                    valid_nonlatest_drefs = false;
                    if (DEBUG>2) {cout<<"valid_nonlatest_drefs=false because of literal "<<*r_ground->context(i)<<endl;}
                  }
                }
            } // End of checking r_ground->context(i)
            if (!r_ground->context(i)->typeCheck()) valid_all = false; //type check
          } // End of filtering of r_ground->context
          if (DEBUG>2) {cout<<endl<<"Candidate ground rule:"<<endl;r_ground->write(cout);}
          if (DEBUG>2) {PRINT(valid_all);  PRINT(valid_action_args);  PRINT(valid_nonlatest_drefs);}
          if (valid_all  &&  valid_action_args  &&  valid_nonlatest_drefs) {
            rules_ground.append(r_ground);
            abstract2GroundMapUnsorted(r).append(rules_ground.num()-1);
            abstract_rules_idxUnsorted.append(r);
            r_last = r_ground;
            if (DEBUG>2) cout<<" --> ACCEPTED"<<endl;
          }
          else if (!valid_action_args) {
            delete r_ground;
            NUM_REJECTED += combos_drefs.N - c2;
            if (DEBUG>2) cout<<" --> REJECTED"<<endl;
            break;
          }
          else if (!valid_nonlatest_drefs) { // valid_action_args = true,  valid_nonlatest_drefs = false
            delete r_ground;
            NUM_REJECTED++;
            if (DEBUG>2) cout<<" --> REJECTED"<<endl;
          }
          else {   // valid_all = false  or stupidContext(.)
            delete r_ground;
            NUM_REJECTED++;
            if (DEBUG>2) cout<<" --> REJECTED"<<endl;
          }
        } // End of else for creation of new rule
      } // End of combos_drefs
    } // End of combos_args
  } // End of abstract rules
  
  
  // resort index maps

  uintA new_order;
  rules_ground.sort_using_args(new_order);

  FOR1D(new_order, i) {
    abstract_rules_idx.append( abstract_rules_idxUnsorted(new_order(i)) );
  }
  FOR1D(abstract2GroundMapUnsorted, i) {
    uint j;
    FOR1D(abstract2GroundMapUnsorted(i),j) {
      abstract2GroundMap(i).append( new_order(abstract2GroundMapUnsorted(i)(j)) );
    }
  }

  if (DEBUG > 0) {
    PRINT(rules_ground.num());
    PRINT(abstract_rules_idxUnsorted);
    PRINT(abstract2GroundMapUnsorted);
    PRINT(new_order);
    PRINT(abstract_rules_idx);
    PRINT(abstract2GroundMap);
  }

  if (delete_nonchanging_concepts) {
    rules_ground.removeNonChangingSymbols();
  }

  FOR1D(rules_ground.ra, r) {
    if (rules_ground.ra(r)->context.containsDoubles()) {    //remove doubles if necessary
      LitL newContext;
      for (int ci = 0; ci < rules_ground.ra(r)->context.N; ci++)
        newContext.setAppend(rules_ground.ra(r)->context(ci));
      rules_ground.ra(r)->context = newContext;
    }
  }
  
  rules_ground.sanityCheck();
  
  if (DEBUG>0) {
    cout<<"# Accepted rules: "<<rules_ground.num()<<endl;
    cout<<"# Rejected rules: "<<NUM_REJECTED<<endl;
  }
  if (DEBUG>0) {cout<<"ruleReasoning::ground_with_filtering [END]"<<endl;}
}


void RuleSet::groundFunctionVars(const RuleSet& rules, RuleSet& rules_expanded) {
  uint r, i, j, x, y;
  FOR1D_(rules, r) {
    Rule* rule = rules.elem(r);
    bool addedNewRules = false;
    FOR1D(rule->context, i) {
      if (rule->context(i)->comparison_type == Literal::comparison_variable) {
        FOR1D(rule->context(i)->s->range, j) {
          Rule *newRule = new Rule();
          newRule->action = rule->action;
          newRule->probs = rule->probs;
          double bound = rule->context(i)->s->range(j);
          Literal *newLit = Literal::get(rule->context(i)->s, rule->context(i)->args, bound, Literal::comparison_equal);
          //build new context and outcomes with replaced any/offset literal
          FOR1D(rule->context, x) {
            if (x != i) newRule->insertContext(rule->context(x));
            else newRule->insertContext(newLit);
          }
          //now ground the outcomes
          FOR1D(rule->outcomes, x) {
            MT::Array<Literal*> newOutcome;
            FOR1D(rule->outcomes(x), y) {
              if (rule->outcomes(x)(y)->s == rule->context(i)->s && rule->outcomes(x)(y)->args == rule->context(i)->args && rule->outcomes(x)(y)->comparison_type == Literal::comparison_offset) {
                double finalValue = bound + rule->outcomes(x)(y)->value;
                if (rule->outcomes(x)(y)->s->range_type == Symbol::integer_set) {   //clamp   
                  uint rangeMin = rule->outcomes(x)(y)->s->range.min();
                  uint rangeMax = rule->outcomes(x)(y)->s->range.max();  
                  if (finalValue < rangeMin) finalValue = rangeMin;
                  if (finalValue > rangeMax) finalValue = rangeMax;
                }
                if (finalValue != bound)
                  newOutcome.append(Literal::get(rule->outcomes(x)(y)->s, rule->outcomes(x)(y)->args, finalValue));
              }
              else newOutcome.append(rule->outcomes(x)(y));
            }
            newRule->outcomes.append(newOutcome);
          }
          addedNewRules = true;
          rules_expanded.append(newRule);
        }
        if (addedNewRules) break;        //only one any literal per call
      }
    }
    if (!addedNewRules) rules_expanded.append(rule);
  }
}


void write(const RuleSet& rules, ostream& out) {
  rules.write(out);
}


void write(const RuleSet& rules, const char* filename) {
  ofstream out;
  out.open(filename);
  CHECK(out.is_open(), "Can't write to simulation sequence file.");
  write(rules, out);
  out.close();
}




/************************************************
 * 
 *     Substitution
 * 
 ************************************************/


int Substitution::globalCounter_Substitution = 0;


Substitution::Substitution() {
  numRefs=0;
  globalCounter_Substitution++;
//   PRINT(globalCounter_Substitution);
}


Substitution::Substitution(const Substitution& s) {
  *this = s;
  numRefs = 0;
  globalCounter_Substitution++;
//   PRINT(globalCounter_Substitution);
}


Substitution::~Substitution() {
  globalCounter_Substitution--;
//   PRINT(globalCounter_Substitution);
}


Substitution& Substitution::operator=(const Substitution& s) {
  this->ins = s.ins;
  this->outs = s.outs;
  return *this;
}


void Substitution::apply(LitL& sub_lits, const LitL& unsub_lits) const {
  sub_lits.resize(unsub_lits.N);
  uint i;
  FOR1D(unsub_lits, i) {
    sub_lits(i) = apply(unsub_lits(i));
  }
}


SymbolicState* Substitution::apply(const SymbolicState& state) const {
  SymbolicState* s_new = new SymbolicState;
  apply(s_new->lits, state.lits);
  return s_new;
}


Literal* Substitution::apply(Literal* in_lit) const {
  uintA out_args(in_lit->args.N);
  uint i;
  FOR1D(in_lit->args, i) {
    out_args(i) = getSubs(in_lit->args(i));
  }
  return Literal::get(in_lit->s, out_args, in_lit->value, in_lit->comparison_type);
}


Rule* Substitution::apply(const Rule& r) const {
  Rule* new_r = new Rule;
  uint i;
  new_r->context.resize(r.context.N);
  FOR1D(r.context, i) {
    new_r->context(i) = apply(r.context(i));
  }
  new_r->action = apply(r.action);
  uint j;
  new_r->outcomes.resize(r.outcomes.N);
  FOR1D(r.outcomes, i) {
    new_r->probs.append(r.probs(i));
    FOR1D(r.outcomes(i), j) {
      new_r->outcomes(i).append(apply(r.outcomes(i)(j)));
    }
  }
  new_r->noise_changes = r.noise_changes;
  new_r->outcome_rewards = r.outcome_rewards;
  return new_r;
}


uint Substitution::getSubs(uint in) const {
//   cout<<"getSubs [START]"<<endl;
//   PRINT(in);
//   PRINT(ins);
//   PRINT(outs);
  uint i = ins.rankInSorted(in, TL::uint_compare);
//   PRINT(i);
  if (i < ins.N) {
    if (in == ins(i)) {
//       cout<<"raus = "<<outs(i)<<endl;
//       cout<<"getSubs [END]"<<endl;
      return outs(i);
    }
  }
//   cout<<"raus = "<< in <<endl;
//   cout<<"getSubs [END]"<<endl;
  return in;
}


bool Substitution::empty() const {
  return ins.N == 0;
}


uint Substitution::num() const {
  return this->ins.N;
}


void Substitution::getIns(uintA& ids) const {
  ids.resize(ins.N);
  uint i;
  FOR1D(ins, i) {
    ids(i) = ins(i);
  }
}


void Substitution::getOuts(uintA& ids) const {
  ids.resize(outs.N);
  uint i;
  FOR1D(outs, i) {
    ids(i) = outs(i);
  }
}


void Substitution::getOutsDistinct(uintA& ids) const {
  ids.clear();
  uint i;
  FOR1D(outs, i) {
    ids.setAppend(outs(i));
  }
}


void Substitution::addSubs2Variable(uint in) {
  // FIND FREE VARIABLE
  uint out;
  for (out=0;;out++) {
    if (outs.findValue(out) < 0)
      break;
  }
  addSubs(in, out);
}


void Substitution::addSubs(uint in, uint out) {
  uint i = ins.rankInSorted(in, TL::uint_compare);
  if (i < ins.N) {
    if (in == ins(i)) { // in already in usage
      outs(i) = out;
    }
    else {  // in is new
      ins.insert(i, in);
      outs.insert(i, out);
    }
  }
  else {
    ins.append(in);
    outs.append(out);
  }
}


bool Substitution::hasSubs(uint in) const {
  uint i = ins.rankInSorted(in, TL::uint_compare);
  if (i < ins.N) {
    if (in == ins(i))
      return true;
  }
  return false;
}


void Substitution::getInverse(Substitution& invSub) const {
  CHECK(invSub.empty(), "Substitution already filled!");
  uint i;
  FOR1D(ins, i) {
    invSub.addSubs(outs(i), ins(i));
  }
}


bool Substitution::mapsToDistinct() const {
  return ins.N == outs.N  &&  !outs.containsDoubles();
}


Substitution* Substitution::combine(Substitution& sub1, Substitution& sub2) {
  Substitution* s = new Substitution;
  uint i;
  uintA ids1;
  sub1.getIns(ids1);
  FOR1D(ids1, i) {
    s->addSubs(ids1(i), sub1.getSubs(ids1(i)));
  }
  uintA ids2;
  sub2.getIns(ids2);
  FOR1D(ids2, i) {
    s->addSubs(ids2(i), sub2.getSubs(ids2(i)));
  }
  return s;
}


void Substitution::write(ostream& os) const {
// 	os << "s ";
	uint i;
	FOR1D(ins, i) {
    if (ins(i) == 0) os << "X";
    else if (ins(i) == 1) os << "Y";
    else if (ins(i) == 2) os << "Z";
    else if (ins(i) == 3) os << "V";
    else if (ins(i) == 4) os << "W";
    else if (ins(i) == 5) os << "U";
    else os << ins(i);
		os << "->" << getSubs(ins(i)) << " ";
	}
}




/************************************************
 * 
 *     SubstitutionSet
 * 
 ************************************************/

SubstitutionSet::SubstitutionSet() {
  sa.memMove = true;
}


SubstitutionSet::~SubstitutionSet() {
  clear();
}


int SubstitutionSet::findValue(Substitution* s) {
  return sa.findValue(s);
}


void SubstitutionSet::clear() {
  uint i;
  FOR1D(sa, i) {
    TL::del(sa(i));
  }
  sa.clear();
}


void SubstitutionSet::append(Substitution* s) {
  sa.append(TL::getRef(s));
}


void SubstitutionSet::append(SubstitutionSet& ss) {
  uint i;
  FOR1D_(ss, i) {
    this->append(ss.elem(i));
  }
}


uint SubstitutionSet::num() const {
  return sa.N;
}


Substitution* SubstitutionSet::elem(uint i) const {
  return sa(i);
}


Substitution* SubstitutionSet::last() const {
  return sa(sa.d0-1);
}


void SubstitutionSet::remove(uint i) {
    TL::del(sa(i));
    sa.remove(i);
}


void SubstitutionSet::overwrite(uint i, Substitution* s) {
    TL::del(sa(i));
    sa(i) = TL::getRef(s);
}


SubstitutionSet& SubstitutionSet::operator=(const SubstitutionSet& ss) {
    this->clear();
    uint i;
    FOR1D(ss.sa, i) {
        this->append(ss.sa(i));
    }
    return *this;
}


void SubstitutionSet::write(ostream& os) {
  uint i;
  for (i=0; i<num(); i++) {
    this->elem(i)->write(os); os<<endl;
  }
}


void SubstitutionSet::createAllPossibleSubstitutions(SubstitutionSet& subs, const uintA& vars, const uintA& constants) {
  uint DEBUG = 0;
  if (DEBUG>0) cout<<"createAllPossibleSubstitutions [START]"<<endl;
  if (DEBUG>0) {PRINT(vars); PRINT(constants);}
  if (vars.N > 0 && constants.N == 0) HALT("No constants: constants="<<constants);
  Substitution* sub;
  uint i;
  uintA assignment(vars.N);
  assignment.setUni(0);
  int id = vars.N-1;
  // beim change ganz nach hinten rutschen und wieder alle zuruecksetzen
  while (true) {
    if (DEBUG>0) PRINT(assignment)
    sub = new Substitution;
    FOR1D(vars, i) {sub->addSubs(vars(i), constants(assignment(i)));}
    if (DEBUG>1) {sub->write(cout); cout<<endl;}
    subs.append(sub);
    while (id>=0 && assignment(id) >= constants.N-1) {
      id--;
    }
    if (id<0)
      break;
    assignment(id)++;
    while (++id < (int) assignment.N) {
      assignment(id) = 0;
    }
    id=assignment.N-1;
  }
  CHECK(subs.num() == pow(constants.N, vars.N), "Oh no, we forgot some subs!")
  if (DEBUG>0) cout<<"createAllPossibleSubstitutions [END]"<<endl;
}



}  // namespace relational


std::ostream& operator<<(std::ostream& os, const relational::Rule& r) {
  r.write(os); return os;
}


std::ostream& operator<<(std::ostream& os, const relational::Substitution& s) {
  s.write(os); return os;
}
