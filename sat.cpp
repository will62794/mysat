#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

/**
 * A (potentially partial) assignment of values to variables of a CNF formula.
 */
class Assignment {
private:
    std::map<std::string, bool> _vals;

public:
    Assignment() {}

    Assignment(std::map<std::string, bool> a) {
        _vals = a;
    }

    std::map<std::string, bool> getVals() {
        return _vals;
    }

    std::set<std::string> getAssignedVars() {
        std::set<std::string> assignedVars;
        for (auto it = _vals.begin(); it != _vals.end(); it++) {
            assignedVars.insert(it->first);
        }
        return assignedVars;
    }

    bool getVal(std::string varName) {
        return _vals[varName];
    }

    bool hasVar(std::string varName) {
        return _vals.find(varName) != _vals.end();
    }

    void set(std::string varName, bool val) {
        _vals.insert(std::make_pair(varName, val));
    }

    void unset(std::string varName) {
        auto it = _vals.find(varName);
        _vals.erase(it);
    }

    std::string toString() {
        if (_vals.size() == 0) {
            return "{}";
        }
        std::string outStr = "{ ";
        for (auto it = _vals.begin(); it != _vals.end(); it++) {
            outStr += it->first;
            outStr += "=";
            outStr += it->second ? "1" : "0";
            outStr += " ";
        }
        return outStr + "}";
    }

    std::string toStringCompact() {
        if (_vals.size() == 0) {
            return "";
        }
        std::string outStr = "";
        for (auto it = _vals.begin(); it != _vals.end(); it++) {
            // outStr += it->first;
            // outStr += "_";
            outStr += (it->second ? "1" : "0");
            outStr += "";
        }
        return outStr;
    }
};

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

    static CNF fromDIMACS(){
        // TODO: Parse DIMACS CNF instances.
    }

    /**
     * Evaluate a CNF given a full (non partial) assignment.
     */
    bool eval(Assignment assgmt) {
        bool allClausesTrue = true;
        for (Clause c : _clauses) {
            bool clauseTrue = false;
            for (auto l : c.getLiterals()) {
                clauseTrue = clauseTrue || l.eval(assgmt.getVal(l.getVarName()));
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
    CNF assign(Assignment assgmt) {
        // Record the set of variables that have been assigned values.
        std::set<std::string> assignedVars = assgmt.getAssignedVars();

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
                    bool assignedVal = assgmt.getVal(lit.getVarName());
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
        return this->assign(Assignment({{l.getVarName(), !l.isNegated()}}));
    }
};

class Context {
public:
    CNF _f;
    // Current variable node.
    std::string _currVar;
    std::string _parentVar;
    // Index of the variable in a given variable ordering.
    int _currVarInd;
    // Current partial assignment.
    Assignment _assmt;
    Assignment _parentAssmt;

    Context() {
        _currVar = "";
        _parentVar = "";
    }

    Context(CNF f,
            std::string currVar,
            std::string parentVar,
            int currVarInd,
            Assignment assmt,
            Assignment parentAssmt) {
        _f = f;
        _currVar = currVar;
        _parentVar = parentVar;
        _currVarInd = currVarInd;
        _assmt = assmt;
        _parentAssmt = parentAssmt;
    }
};

// Termination tree edge.
class TreeEdge {
public:
    std::string from;
    std::string to;
    std::string fromid;
    std::string toid;

    TreeEdge(std::string f, std::string t, std::string ida, std::string idb) {
        from = f;
        to = t;
        fromid = ida;
        toid = idb;
    }

    bool operator<(const TreeEdge& other) const {
        return std::make_tuple(other.from, other.to, other.fromid, other.toid) <
            std::make_tuple(from, to, fromid, toid);
    }
};

/**
 * Satisfiability solver.
 */
class Solver {
private:
    Assignment _currAssignment;

    // Records the termination tree of the SAT solver execution i.e.
    // a record of the path in the tree it visited during solving.
    std::set<TreeEdge> terminationTree;

    bool enableUnitPropagation = true;

public:
    Solver() {}

    Assignment getAssignment() {
        return _currAssignment;
    }

    void printTerminationTree() {
        std::cout << "termination tree:" << std::endl;
        for (auto e : terminationTree) {
            std::cout << "\"" << e.from + " | " + e.fromid << "\" -> \"" << e.to + " | " + e.toid
                      << "\"" << std::endl;
        }
        std::cout << "===" << std::endl;
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
        frontier.push_back(Context(f, varList.at(0), "root", 0, initVals, {}));

        Context lastNode = Context();

        // Explore all possible assignments in a depth first manner.
        while (frontier.size() > 0) {
            Context currNode = frontier.back();
            frontier.pop_back();

            std::cout << "* currVar: '" << currNode._currVar << "', varInd=" << currNode._currVarInd
                      << ", parentVar: " << currNode._parentVar << std::endl;
            std::cout << "parentVar: " << currNode._parentVar << std::endl;
            std::cout << "curr assignment: " << currNode._assmt.toString() << std::endl;
            std::cout << "par assignment: " << currNode._parentAssmt.toString() << std::endl;

            // Record information for termination tree.
            auto localCurrAssmt = currNode._assmt;
            TreeEdge e = TreeEdge(currNode._parentVar,
                                  currNode._currVar,
                                  currNode._parentAssmt.toStringCompact(),
                                  currNode._assmt.toStringCompact());
            terminationTree.insert(e);

            // Reduce based on current assignment.
            auto currAssmt = currNode._assmt;
            CNF currF = currNode._f;
            CNF fassigned = currF.assign(currAssmt);
            std::cout << "f after assignment: " << fassigned.toString() << std::endl;

            // Close the formula under unit resolution, if enabled.
            if (enableUnitPropagation) {
                while (fassigned.hasUnitClause()) {
                    auto unitLit = fassigned.getUnitLiteral();
                    // Assign the unit literal's variable to truthify the unit clause
                    // and simplify the formula based on this assignment.
                    currAssmt.set(unitLit.getVarName(), !unitLit.isNegated());
                    fassigned = fassigned.unitPropagate(unitLit);
                }
                std::cout << "f after unit prop: " << fassigned.toString() << std::endl;
                std::cout << "curr assignment after unit prop: " << currAssmt.toString()
                          << std::endl;
            }

            if (fassigned.isEmpty()) {
                // If we determined that the formula is necessarily satisfiable
                // under the current partial assignment, then we can fill in
                // arbitrary values for the remaining, unassigned variables.
                for (auto v : varList) {
                    // Assign a value to the currently unassigned variable.
                    if (!currAssmt.hasVar(v)) {
                        bool arbitraryVal = true;
                        currAssmt.set(v, arbitraryVal);
                    }
                }
                _currAssignment = currAssmt;
                return true;
            }
            if (fassigned.hasEmptyClause()) {
                // Current assignment is necessarily UNSAT, so no need to
                // explore further down this branch.
                std::cout << "fassigned has empty clause." << std::endl;
                lastNode = currNode;
                continue;
            }

            int varNextInd = currNode._currVarInd + 1;
            std::cout << "varNextInd: " << varNextInd << std::endl;

            // Check if we reached the last variable i.e. a leaf.
            if (varNextInd <= varList.size()) {

                // Explore each possible true/false assignment for this variable.
                Assignment tAssign(currAssmt);
                Assignment fAssign(currAssmt);

                tAssign.set(currNode._currVar, true);
                fAssign.set(currNode._currVar, false);

                // If we have reached the last variable, then pass in a dummy
                // 'leaf' variable node for the next level in the search tree.
                // This value isn't explicitly used, but it makes it convenient
                // so that we can traverse one more level deep in the search
                // tree.
                auto nextVar = varNextInd == varList.size() ? "LEAF" : varList.at(varNextInd);

                std::string parentVar = currNode._currVar;
                Context tctx(fassigned, nextVar, parentVar, varNextInd, tAssign, currNode._assmt);
                Context fctx(fassigned, nextVar, parentVar, varNextInd, fAssign, currNode._assmt);

                frontier.push_back(tctx);
                frontier.push_back(fctx);
            }

            lastNode = currNode;
        }

        return false;
    }
};
