#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

/**
 * A literal is a boolean variable that may be either negated or non-negated.
 */
class Literal {
private:
    // The name of the boolean variable.
    std::string _name;
    // Is this a negated variable symbol.
    bool _neg;

public:
    /**
     * Construct a new Literal, with specified negation.
     */
    Literal(std::string name, bool neg) {
        _name = name;
        _neg = neg;
    }

    /**
     * Construct a literal given as a string [<neg>]<name>, where <neg> is
     * an optional "~" character. For example "~x" constructs the negation
     * of variable x, and "x" constructs the non-negated variable x.
     */
    Literal(std::string name) {
        if (name.at(0) == '~') {
            _name = name.substr(1);
            _neg = true;
        } else {
            _name = name;
            _neg = false;
        }
    }

    std::string toString() {
        auto prefix = _neg ? "~" : "";
        return prefix + _name;
    }

    std::string getVarName() {
        return _name;
    }

    /**
     * Evaluate this literal for a given truth value of its variable.
     */
    bool eval(bool val) {
        return _neg ? !val : val;
    }
};


/**
 * A clause is a set of literals L, representing the logical disjunction of the
 * literals in L.
 */
class Clause {
private:
    std::vector<Literal> _elems;

public:
    Clause(std::vector<Literal> elems) {
        _elems = elems;
    }

    std::string toString() {
        std::string outStr;
        for (auto e : _elems) {
            outStr += e.toString();
            outStr += " | ";
        }
        // Remove the extra disjunction symbol.
        return outStr.substr(0, outStr.length() - 3);
    }

    /**
     * Return the set of variables that appear in this clause.
     */
    std::vector<Literal> getLiterals() {
        return _elems;
    }

    /**
     * Return the set of variables that appear in this clause.
     */
    std::set<std::string> getVariableSet() {
        std::set<std::string> outSet;
        for (auto e : _elems) {
            outSet.insert(e.getVarName());
        }
        return outSet;
    }
};

/**
 * A boolean formula in conjunctive normal form (CNF) is a conjunction of clauses.
 */
class CNF {
private:
    std::vector<Clause> _clauses;

public:
    CNF() {}

    CNF(std::vector<Clause> clauses) {
        _clauses = clauses;
    }

    /**
     * Returns whether this CNF contains no clauses.
     */
    bool isEmpty() {
        return _clauses.size() == 0;
    }

    /**
     * Returns whether this CNF contains some empty clause.
     */
    bool hasEmptyClause() {
        if (isEmpty()) {
            return false;
        }

        for (auto c : _clauses) {
            if (c.getLiterals().size() == 0) {
                return true;
            }
        }
        return false;
    }


    /**
     * Returns the set of variables that appear in the CNF formula.
     */
    std::set<std::string> getVariableSet() {
        std::set<std::string> outSet;
        for (auto c : _clauses) {
            auto clauseVars = c.getVariableSet();
            for (auto v : clauseVars) {
                outSet.insert(v);
            }
        }
        return outSet;
    }

    std::string toString() {
        std::string outStr;
        for (auto e : _clauses) {
            outStr += "(" + e.toString() + ")";
            outStr += " & ";
        }
        // Remove the extra conjunction symbol.
        return outStr.substr(0, outStr.length() - 2);
    }

    /**
     * Assign values to the specified variables and return a simplified CNF formula.
     *
     * An assignment is given as a map from variable names to boolean values.
     */
    CNF assign(std::map<std::string, bool> assgmt) {
        // Record the set of variables that have been assigned values.
        std::set<std::string> assignedVars;
        for (auto it = assgmt.begin(); it != assgmt.end(); it++) {
            assignedVars.insert(it->first);
        }

        // Initialize the set of clauses that will appear in the resulting CNF after assignment.
        std::vector<Clause> outClauses;
        for (auto c : _clauses) {
            std::vector<Literal> outClauseLits;
            bool clauseIsTrue = false;
            // For each literal in this clause, check if it is assigned a value.
            for (auto lit : c.getLiterals()) {
                if (clauseIsTrue) {
                    continue;
                }
                if (assignedVars.find(lit.getVarName()) != assignedVars.end()) {
                    // If this literal evaluates to TRUE under the
                    // assignment, then the whole clause is satisfied, so we
                    // remove it i.e. don't add it to the output CNF.
                    bool assignedVal = assgmt[lit.getVarName()];
                    if (lit.eval(assignedVal)) {
                        clauseIsTrue = true;
                    } else {
                        // If the literal evaluates to FALSE, then we simply
                        // remove it from the output clause i.e. we don't
                        // add it.
                        continue;
                    }

                } else {
                    // If this literal is not assigned a value, then simply add it the
                    // current clause as is.
                    outClauseLits.push_back(lit);
                }
            }

            // Add this clause if it was not necessarily satisfied under the
            // given assignment.
            if (!clauseIsTrue) {
                Clause newClause = Clause(outClauseLits);
                outClauses.push_back(newClause);
            }
        }

        return CNF(outClauses);
    }
};

class Context {
public:
    // Current variable node.
    std::string _currVar;
    // Index of the variable in a given variable ordering.
    int _currVarInd;
    // Current partial assignment.
    std::map<std::string, bool> _assmt;

    Context(std::string currVar, int currVarInd, std::map<std::string, bool> assmt) {
        _assmt = assmt;
        _currVar = currVar;
        _currVarInd = currVarInd;
    }
};

// TODO: Move this into a dedicated class representing a partial variable assignment.
std::string assignmentToString(std::map<std::string, bool> a) {
    if (a.size() == 0) {
        return "{}";
    }
    std::string outStr = "{ ";
    for (auto it = a.begin(); it != a.end(); it++) {
        outStr += it->first;
        outStr += "=";
        outStr += it->second ? "1" : "0";
        outStr += " ";
    }
    return outStr + "}";
}

/**
 * Satisfiability solver.
 */
class Solver {
private:
    std::map<std::string, bool> _currAssignment;

public:
    Solver() {}

    std::map<std::string, bool> getAssignment() {
        return _currAssignment;
    }

    bool isSat(CNF f) {
        std::cout << "# checking isSat # " << std::endl;

        auto varSet = f.getVariableSet();
        std::vector<std::string> varList;
        // Initialize list of variables in some fixed, arbitrary order.
        for (auto v : varSet) {
            varList.push_back(v);
        }


        // The set of tree nodes in the frontier i.e. discovered but not explored yet.
        std::vector<Context> frontier;

        std::map<std::string, bool> initVals;
        frontier.push_back(Context(varList.at(0), 0, initVals));

        // Explore all possible assignments in a depth first manner.
        while (frontier.size() > 0) {
            Context currNode = frontier.back();
            frontier.pop_back();

            std::cout << "currVar: '" << currNode._currVar << "', varInd=" << currNode._currVarInd
                      << std::endl;
            std::cout << "curr assignment: " << assignmentToString(currNode._assmt) << std::endl;

            // Reduce based on current assignment.
            auto currAssmt = currNode._assmt;
            CNF fassigned = f.assign(currAssmt);
            std::cout << "fassigned: " << fassigned.toString() << std::endl;
            if (fassigned.isEmpty()) {
                _currAssignment = currAssmt;
                return true;
            }
            if (fassigned.hasEmptyClause()) {
                // Current assignment is necessarily UNSAT, so no need to
                // explore further down this branch.
                continue;
            }

            int varNextInd = currNode._currVarInd + 1;

            // Explore each possible true/false assignment for this variable.
            std::map<std::string, bool> tAssign(currNode._assmt);
            std::map<std::string, bool> fAssign(currNode._assmt);

            tAssign.insert(std::make_pair(currNode._currVar, true));
            fAssign.insert(std::make_pair(currNode._currVar, false));

            // Reached the last variable i.e. a leaf.
            if (varNextInd >= varList.size()) {
                continue;
            }

            Context tctx(varList.at(varNextInd), varNextInd, tAssign);
            Context fctx(varList.at(varNextInd), varNextInd, fAssign);

            frontier.push_back(tctx);
            frontier.push_back(fctx);
        }

        return false;
    }
};


int main(int argc, char const* argv[]) {
    auto l1 = Literal("~x");
    auto l2 = Literal("y");
    std::cout << l1.toString() << std::endl;
    std::cout << l2.toString() << std::endl;

    auto c1 = Clause({l1, l2});
    auto c2 = Clause({l1});
    std::cout << c1.toString() << std::endl;
    std::cout << c2.toString() << std::endl;

    // (x \/ y) /\ (~y \/ ~z)
    auto lx = Literal("x");
    auto ly = Literal("y");
    auto lny = Literal("~y");
    auto lnz = Literal("~z");

    auto cxy = Clause({lx, ly});
    auto cynz = Clause({lny, lnz});
    auto fa = CNF({cxy, cynz});

    std::map<std::string, bool> ax1 = {{"x", true}};
    std::map<std::string, bool> ax0 = {{"x", false}};
    std::cout << "fa: " << fa.toString() << std::endl;
    std::cout << "fa (x=1): " << fa.assign(ax1).toString() << std::endl;
    std::cout << "fa (x=0): " << fa.assign(ax0).toString() << std::endl;


    Solver solver = Solver();
    auto ret = solver.isSat(fa);
    std::cout << "fa isSat: " << ret << std::endl;
    std::cout << "fa assignment: " << assignmentToString(solver.getAssignment()) << std::endl;

    return 0;
}
