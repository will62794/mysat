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

    bool isNegated() {
        return _neg;
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

    Clause(std::vector<std::string> elems) {
        for (auto ls : elems) {
            _elems.push_back(Literal(ls));
        }
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
     * Returns the number of literals in this clause.
     */
    int size() {
        return _elems.size();
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

    CNF(std::vector<std::vector<std::string>> clauses) {
        for (auto c : clauses) {
            _clauses.push_back(Clause(c));
        }
    }

    static CNF randomCNF(int nclauses, int nvars, int clause_size) {
        std::vector<std::string> allVars;
        for (int i = 0; i < nvars; i++) {
            std::string baseVarName = "x";
            std::string varIndStr = std::to_string(i);
            allVars.push_back(baseVarName + varIndStr);
        }

        std::vector<Clause> allClauses;
        for (int i = 0; i < nclauses; i++) {
            std::vector<Literal> currClauseLits;
            for (int j = 0; j < clause_size; j++) {
                int randVarInd = rand() % nvars;
                auto randVar = allVars[randVarInd];
                auto negation = rand() % 2 ? "~" : "";
                currClauseLits.push_back(Literal(negation + randVar));
            }
            allClauses.push_back(Clause(currClauseLits));
        }
        return CNF(allClauses);
    }

    /**
     * Is this CNF formula equivalent to 'true' i.e. is it an empty set of
     * clauses.
     */
    bool isTrue() {
        return isEmpty();
    }

    /**
     * Is this CNF formula equivalent to 'true' i.e. does it contain a single
     * empty clause.
     */
    bool isFalse() {
        return _clauses.size() == 1 && _clauses.at(0).getLiterals().size() == 0;
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

    /**
     * Returns the set of variables that appear in the CNF formula in a fixed, arbitrary order.
     */
    std::vector<std::string> getVariableList() {
        auto varSet = getVariableSet();
        std::vector<std::string> varList;
        for (auto v : varSet) {
            varList.push_back(v);
        }
        return varList;
    }

    std::string toString() {
        if (_clauses.size() == 0) {
            return "{}";
        }
        std::string outStr = "{ ";
        for (auto e : _clauses) {
            outStr += "(" + e.toString() + ")";
            outStr += " & ";
        }
        // Remove the extra conjunction symbol.
        return outStr.substr(0, outStr.length() - 2) + "}";
    }

    /**
     * Generates an encoding of this CNF formula in the DIMACS CNF file format.
     *
     * See https://people.sc.fsu.edu/~jburkardt/data/cnf/cnf.html for a
     * description of the file format.
     */
    std::string toDIMACS() {
        std::string outStr;

        // Write the 'problem' line.
        // p cnf <num_vars> <num_clauses>
        std::string numVarsStr = std::to_string(getVariableList().size());
        outStr += "p cnf " + numVarsStr + " " + std::to_string(_clauses.size());
        outStr += "\n";

        // Map from varialbe names to numeric identifiers.
        std::map<std::string, int> varNumMap;
        for (int ind = 0; ind < getVariableList().size(); ind++) {
            auto v = getVariableList()[ind];
            // Variable identifiers in DIMACS are 1-indexed.
            varNumMap.insert(std::make_pair(v, ind + 1));
        }

        // Write each clause.
        for (auto c : _clauses) {
            for (auto l : c.getLiterals()) {
                int varNum = varNumMap[l.getVarName()];
                std::string neg = l.isNegated() ? "-" : "";
                outStr += (neg + std::to_string(varNum));
                outStr += " ";
            }
            // The definition of each clause is terminated by a final value of "0".
            outStr += "0\n";
        }
        return outStr;
    }

    /**
     * Evaluate a CNF given a full (non partial) assignment.
     */
    bool eval(std::map<std::string, bool> assgmt) {
        bool allClausesTrue = true;
        for (Clause c : _clauses) {
            bool clauseTrue = false;
            for (auto l : c.getLiterals()) {
                clauseTrue = clauseTrue || l.eval(assgmt[l.getVarName()]);
            }
            allClausesTrue = allClausesTrue && clauseTrue;
        }
        return allClausesTrue;
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

    /**
     * Does there exist a unit clause in this CNF.
     */
    bool hasUnitClause() {
        for (auto c : _clauses) {
            if (c.size() == 1) {
                return true;
            }
        }
        return false;
    }

    /**
     * Get some literal which has a unit clause in this CNF.
     *
     * Assumes there exists some unit clause.
     */
    Literal getUnitLiteral() {
        for (auto c : _clauses) {
            if (c.size() == 1) {
                return c.getLiterals()[0];
            }
        }
        // Should never reach.
        assert(false);
        return Literal("");
    }

    /**
     * Apply unit resolution to the current CNF and return a new formula with the result.
     *
     * Assumes that the given literal is a unit clause in this CNF.
     */
    CNF unitPropagate(Literal l) {
        // std::vector<Clause> outClauses;
        return this->assign({{l.getVarName(), !l.isNegated()}});
    }
};

class Context {
public:
    CNF _f;
    // Current variable node.
    std::string _currVar;
    // Index of the variable in a given variable ordering.
    int _currVarInd;
    // Current partial assignment.
    std::map<std::string, bool> _assmt;

    Context(CNF f, std::string currVar, int currVarInd, std::map<std::string, bool> assmt) {
        _f = f;
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

    bool _isSatBruteForceRec(CNF f,
                             std::vector<std::string> varList,
                             int varInd,
                             std::map<std::string, bool> currAssign) {
        if (varInd == varList.size()) {
            return f.eval(currAssign);
        } else {
            std::map<std::string, bool> tAssign(currAssign);
            tAssign.insert(std::make_pair(varList.at(varInd), true));

            std::map<std::string, bool> fAssign(currAssign);
            fAssign.insert(std::make_pair(varList.at(varInd), false));

            bool tbranch = _isSatBruteForceRec(f, varList, varInd + 1, tAssign);
            bool fbranch = _isSatBruteForceRec(f, varList, varInd + 1, fAssign);
            return tbranch || fbranch;
        }
    }

    /**
     * Check satisfiability by brute force exploration of all possible assignments.
     *
     * Return true if satisfiable and false if unsatisfiable.
     */
    bool isSatBruteForce(CNF f) {
        auto varList = f.getVariableList();
        std::map<std::string, bool> initAssign;
        return _isSatBruteForceRec(f, varList, 0, initAssign);
    }

    bool isSat(CNF f) {
        std::cout << "# checking isSat # " << std::endl;

        std::vector<std::string> varList = f.getVariableList();

        // The set of tree nodes in the frontier i.e. discovered but not explored yet.
        std::vector<Context> frontier;

        std::map<std::string, bool> initVals;
        frontier.push_back(Context(f, varList.at(0), 0, initVals));

        // Explore all possible assignments in a depth first manner.
        while (frontier.size() > 0) {
            Context currNode = frontier.back();
            frontier.pop_back();

            std::cout << "currVar: '" << currNode._currVar << "', varInd=" << currNode._currVarInd
                      << std::endl;
            // std::cout << "curr assignment: " << assignmentToString(currNode._assmt) << std::endl;

            // Reduce based on current assignment.
            auto currAssmt = currNode._assmt;
            CNF currF = currNode._f;
            CNF fassigned = currF.assign(currAssmt);

            // Close the formula under unit resolution.
            while (fassigned.hasUnitClause()) {
                auto unitLit = fassigned.getUnitLiteral();
                fassigned = fassigned.unitPropagate(unitLit);
            }

            std::cout << "fassigned: " << fassigned.toString() << std::endl;
            if (fassigned.isEmpty()) {
                _currAssignment = currAssmt;
                std::cout << "fassigned is empty." << std::endl;
                return true;
            }
            if (fassigned.hasEmptyClause()) {
                // Current assignment is necessarily UNSAT, so no need to
                // explore further down this branch.
                std::cout << "fassigned has empty clause." << std::endl;
                continue;
            }

            int varNextInd = currNode._currVarInd + 1;
            std::cout << "varNextInd: " << varNextInd << std::endl;

            // Check if we reached the last variable i.e. a leaf.
            if (varNextInd <= varList.size()) {

                // Explore each possible true/false assignment for this variable.
                std::map<std::string, bool> tAssign(currNode._assmt);
                std::map<std::string, bool> fAssign(currNode._assmt);

                tAssign.insert(std::make_pair(currNode._currVar, true));
                fAssign.insert(std::make_pair(currNode._currVar, false));

                // If we have reached the last variable, then pass in a dummy
                // 'leaf' variable node for the next level in the search tree.
                // This value isn't explicitly used, but it makes it convenient
                // so that we can traverse one more level deep in the search
                // tree.
                auto nextVar = varNextInd == varList.size() ? "LEAF" : varList.at(varNextInd);

                Context tctx(fassigned, nextVar, varNextInd, tAssign);
                Context fctx(fassigned, nextVar, varNextInd, fAssign);

                frontier.push_back(tctx);
                frontier.push_back(fctx);
            }
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

    std::map<std::string, bool> ax1 = {{"x", true}, {"y", false}};
    std::map<std::string, bool> ax0y1z1 = {{"x", false}, {"y", true}, {"z", true}};
    std::cout << "fa: " << fa.toString() << std::endl;
    std::cout << "fa (x=1): " << fa.assign(ax1).toString() << std::endl;
    std::cout << "fa (x0y1z1): " << fa.assign(ax0y1z1).toString() << std::endl;
    std::cout << "fa (x0y1z1) false: " << fa.assign(ax0y1z1).isFalse() << std::endl;


    Solver solver = Solver();
    auto ret = solver.isSat(fa);
    std::cout << "fa isSat: " << ret << std::endl;
    std::cout << "fa assignment: " << assignmentToString(solver.getAssignment()) << std::endl;


    auto a1 = solver.getAssignment();
    auto ret1 = fa.assign(a1);
    std::cout << "fa:" << ret1.toString() << std::endl;
    std::cout << "fa true:" << ret1.isTrue() << std::endl;

    //
    // Preliminary test cases.
    //

    CNF ct1 = CNF({{"x", "y"}, {"~y"}});
    CNF ct2 = CNF({{"x", "y"}, {"~x"}});
    CNF ct3 = CNF({{"x"}, {"~x"}});
    CNF ct4 = CNF({{"x", "y"}, {"x", "~y"}, {"~x", "y"}, {"~x", "~y"}});

    assert(solver.isSatBruteForce(ct1));
    assert(solver.isSatBruteForce(ct2));
    assert(!solver.isSatBruteForce(ct3));
    assert(!solver.isSatBruteForce(ct4));

    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    assert(solver.isSatBruteForce(ct2) == solver.isSat(ct2));
    assert(solver.isSatBruteForce(ct3) == solver.isSat(ct3));
    assert(solver.isSatBruteForce(ct4) == solver.isSat(ct4));

    auto d1 = ct4.toDIMACS();
    std::cout << d1;

    // Generate random CNF SAT formulas.
    srand(time(NULL));

    int niters = 5;
    for (int k = 0; k < niters; k++) {
        int nclauses = 40;
        int nvars = 8;
        int clause_size = 3;
        auto cf = CNF::randomCNF(nclauses, nvars, clause_size);
        std::cout << cf.toString() << std::endl;
        auto retOracle = solver.isSatBruteForce(cf);
        assert(retOracle == solver.isSat(cf));
        std::cout << "Result: " << (retOracle ? "SAT" : "UNSAT") << std::endl;
        std::cout << cf.toDIMACS() << std::endl;
    }

    CNF t1 = CNF({{"x", "y", "z"}, {"z"}});
    std::cout << "t1: " << t1.toString() << std::endl;
    auto t1p = t1.unitPropagate(Literal("z"));
    std::cout << "t1p: " << t1p.toString() << std::endl;

    return 0;
}
