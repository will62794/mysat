#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP

std::vector<std::string> split(std::string const& str, std::string delim) {
    std::vector<std::string> out;
    size_t start;
    size_t end = 0;

    while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
        end = str.find(delim, start);
        out.push_back(str.substr(start, end - start));
    }
    return out;
}

/**
 * A (potentially partial) assignment of values to variables of a CNF formula.
 */
class Assignment {
private:
    std::map<std::string, bool> _vals;
    // Allows for recording of decision levels for each variable assignment.
    std::map<std::string, int> _levels;

public:
    Assignment() {}

    Assignment(std::map<std::string, bool> a) {
        _vals = a;
    }

    Assignment(std::map<std::string, bool> a, std::map<std::string, int> levels) {
        _vals = a;
        _levels = levels;
    }

    std::map<std::string, bool> getVals() {
        return _vals;
    }

    int getDecisionLevel(std::string varname) {
        return _levels.at(varname);
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

    void set(std::string varName, bool val, int level) {
        _vals.insert(std::make_pair(varName, val));
        _levels.insert(std::make_pair(varName, level));
    }

    void unset(std::string varName) {
        auto it = _vals.find(varName);
        _vals.erase(it);
        auto it2 = _levels.find(varName);
        if (it2 != _levels.end()) {
            _levels.erase(it2);
        }
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
            if (_levels.find(it->first) != _levels.end()) {
                outStr += "(l" + std::to_string(_levels.at(it->first)) + ")";
            }
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

    std::string getVarName() const {
        return _name;
    }

    bool isNegated() const {
        return _neg;
    }

    /**
     * Evaluate this literal for a given truth value of its variable.
     */
    bool eval(bool val) {
        return _neg ? !val : val;
    }

    bool operator<(const Literal& other) const {
        return std::make_tuple(_neg, _name) <
            std::make_tuple(other.isNegated(), other.getVarName());
    }

    /**
     * Return the negated version of this literal.
     */
    Literal negate() {
        return Literal(getVarName(), !_neg);
    }
};


/**
 * A clause is a set of literals L, representing the logical disjunction of the
 * literals in L.
 */
class Clause {
private:
    std::vector<Literal> _elems;

    // Represents a clause that has one of its literals evaluate to true, so
    // trivially evaluates to true.
    bool _isSat = false;

public:
    Clause(std::vector<Literal> elems) {
        _elems = elems;
    }

    Clause(std::vector<std::string> elems) {
        for (auto ls : elems) {
            _elems.push_back(Literal(ls));
        }
    }

    static Clause makeTrueClause() {
        std::vector<Literal> empty;
        Clause trueClause(empty);
        trueClause.setTrue();
        return trueClause;
    }

    void setTrue() {
        _isSat = true;
    }

    bool isTrue() {
        return _isSat;
    }

    /**
     * Evaluate truth value of a clause under a given assignment.
     */
    bool eval(Assignment a) {
        if (_isSat) {
            return true;
        }
        bool clauseTrue = false;
        for (auto l : getLiterals()) {
            clauseTrue = clauseTrue || l.eval(a.getVal(l.getVarName()));
        }
        return clauseTrue;
    }

    std::string toString() {
        std::string outStr;
        if (_isSat) {
            return "T";
        }
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

    /**
     * Does this clause contain a literal with the given variable.
     */
    bool hasVariable(std::string varname) {
        for (auto l : _elems) {
            if (l.getVarName() == varname) {
                return true;
            }
        }
        return false;
    }

    /**
     * Perform resolution on the two given clauses and return the resolvent.
     *
     * Assumes that each clause contains a literal with opposite polarity i.e.
     *
     * C1 = ~L \/ a1 \/ ... \/ ak
     * C2 =  L \/ b1 \/ ... \/ bk
     *
     * where L is given as 'varToResolveOn'.
     */
    static Clause resolve(Clause c1, Clause c2, std::string varToResolveOn) {
        std::set<Literal> resolvedClauseLits;
        for (auto l : c1.getLiterals()) {
            if (l.getVarName() != varToResolveOn) {
                resolvedClauseLits.insert(l);
            }
        }
        for (auto l : c2.getLiterals()) {
            if (l.getVarName() != varToResolveOn) {
                resolvedClauseLits.insert(l);
            }
        }
        std::vector<Literal> litVec;
        for (auto l : resolvedClauseLits) {
            litVec.push_back(l);
        }

        return Clause(litVec);
    }

    /**
     * Return the given clause with all literals negated.
     */
    Clause negate() {
        std::vector<Literal> outLits;
        for (auto l : getLiterals()) {
            outLits.push_back(l.negate());
        }
        return Clause(outLits);
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
            std::set<Literal> currClauseLits;
            for (int j = 0; j < clause_size; j++) {
                int randVarInd = rand() % nvars;
                auto randVar = allVars[randVarInd];
                auto negation = rand() % 2 ? "~" : "";
                currClauseLits.insert(Literal(negation + randVar));
            }
            std::vector<Literal> litVec(currClauseLits.begin(), currClauseLits.end());
            allClauses.push_back(Clause(litVec));
        }
        return CNF(allClauses);
    }

    /**
     * Appends the given clause to this CNF formula, modifying it in place.
     */
    void appendClause(Clause c) {
        _clauses.push_back(c);
    }

    /**
     * Is this CNF formula equivalent to 'true' i.e. are all of its clauses satisfied.
     */
    bool isTrue() {
        bool allSat = true;
        for (auto c : _clauses) {
            allSat = allSat && c.isTrue();
        }
        return allSat;
    }

    /**
     * Is this CNF formula equivalent to 'false' i.e. are all of its clauses unsatisfied.
     */
    bool isFalse() {
        bool allUnsat = true;
        for (auto c : _clauses) {
            allUnsat = allUnsat && c.size() == 1 && c.getLiterals().size() == 0;
        }
        return allUnsat;

        // return _clauses.size() == 1 && _clauses.at(0).getLiterals().size() == 0;
    }

    /**
     * Returns whether this CNF contains no clauses.
     */
    bool isEmpty() {
        return _clauses.size() == 0;
    }

    // Is the current CNF satisfied i.e. all clauses trivially evaluate to true.
    bool isSatisfied() {
        if (_clauses.empty()) {
            return true;
        }
        bool allSat = true;
        for (auto c : _clauses) {
            allSat = allSat && c.isTrue();
        }
        return allSat;
    }

    /**
     * Returns whether this CNF contains some empty clause.
     */
    bool hasEmptyClause() {
        if (isEmpty()) {
            return false;
        }

        for (auto c : _clauses) {
            // Check for empty clause that is not trivially true.
            if (c.getLiterals().size() == 0 && !c.isTrue()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Return index of the first empty clause in this CNF.
     */
    int getEmptyClause() {
        for (int ci = 0; ci < _clauses.size(); ci++) {
            auto c = _clauses.at(ci);
            // Check for empty clause that is not trivially true.
            if (c.getLiterals().size() == 0 && !c.isTrue()) {
                return ci;
            }
        }
        return -1;
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

    static CNF fromDIMACS(std::string dimacsStr) {
        int nvars;
        int nclauses;
        std::vector<Clause> allClauses;
        bool parsedProblemLine = false;

        auto lines = split(dimacsStr, "\n");
        for (auto line : lines) {
            // Comment line.
            if (line.at(0) == 'c') {
                continue;
            }
            // Problem line.
            // p cnf <num_vars> <num_clauses>
            else if (line.at(0) == 'p') {
                auto args = split(line, " ");
                assert(args[0] == "p");
                assert(args[1] == "cnf");  // Must be CNF type.
                nvars = stoi(args[2]);
                nclauses = stoi(args[3]);
                parsedProblemLine = true;
            }
            // Read another clause line.
            else if (parsedProblemLine) {
                // Assume all clause vars are on same line for now.
                std::vector<std::string> clauseVars = split(line, " ");
                std::vector<Literal> clauseLits;
                for (auto cvar : clauseVars) {
                    if (cvar == "0") {
                        // Clauses are terminated by a zero, which is not
                        // considered as a variable in the clause.
                        break;
                    }

                    // Minus sign is used to indicate a negated variable.
                    std::string varPrefix = "x";
                    if (cvar.at(0) == '-') {
                        bool neg = true;
                        clauseLits.push_back(Literal(varPrefix + cvar.substr(1), neg));
                    } else {
                        clauseLits.push_back(Literal(varPrefix + cvar));
                    }
                }
                Clause newClause = Clause(clauseLits);
                allClauses.push_back(newClause);
            }
        }
        assert(allClauses.size() == nclauses);
        return CNF(allClauses);
    }

    /**
     * Evaluate a CNF given a full (non partial) assignment.
     */
    bool eval(Assignment assgmt) {
        bool allClausesTrue = true;
        for (Clause c : _clauses) {
            allClausesTrue = allClausesTrue && c.eval(assgmt);
            // if(c.setTrue)
            // bool clauseTrue = false;
            // for (auto l : c.getLiterals()) {
            // clauseTrue = clauseTrue || l.eval(assgmt.getVal(l.getVarName()));
            // }
            // clauseTrue;
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
            if (c.isTrue()) {
                outClauses.push_back(Clause::makeTrueClause());
                continue;
            }

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
                Clause newClause(outClauseLits);
                outClauses.push_back(newClause);
            } else {
                // Represent a satisfied clause as trivially true, but don't
                // remove it from the formula.
                std::vector<Literal> emptyEls;
                Clause newClause(emptyEls);
                newClause.setTrue();
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
     * Return index of a unit clause, if one exists.
     */
    int getUnitClause() {
        for (int ind = 0; ind < _clauses.size(); ind++) {
            auto c = _clauses[ind];
            if (c.size() == 1) {
                return ind;
            }
        }
        // No unit clause exists.
        return -1;
    }

    std::vector<Clause> getClauses() {
        return _clauses;
    }

    /**
     * Return the nth clause from this CNF, 0-indexed.
     */
    Clause getClause(int ind) {
        return _clauses.at(ind);
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

    // Decision level that can be stored during SAT tree exploration.
    int _decisionLevel = -1;

    Context() {
        _currVar = "";
        _parentVar = "";
    }

    Context(CNF f,
            std::string currVar,
            std::string parentVar,
            int currVarInd,
            Assignment assmt,
            Assignment parentAssmt,
            int decisionLevel) {
        _f = f;
        _currVar = currVar;
        _parentVar = parentVar;
        _currVarInd = currVarInd;
        _assmt = assmt;
        _parentAssmt = parentAssmt;
        _decisionLevel = decisionLevel;
    }

    int getDecisionLevel() {
        return _decisionLevel;
    }
};

class TreeNode {
public:
    std::string name;
    std::string id;
    std::string currF;
    TreeNode() {}
    TreeNode(std::string f, std::string i, std::string cf) {
        name = f;
        id = i;
        currF = cf;
    }
    bool operator<(const TreeNode& other) const {
        return std::make_tuple(other.name, other.id, other.currF) <
            std::make_tuple(name, id, currF);
    }

    std::string getDOTId() {
        return name + " | " + id + " | " + currF;
    }
};

// Termination tree edge.
class TreeEdge {
public:
    TreeNode from;
    TreeNode to;

    TreeEdge(TreeNode f, TreeNode t) {
        from = f;
        to = t;
    }

    TreeEdge(std::string f, std::string t, std::string ida, std::string idb) {
        from = TreeNode(f, ida, "");
        to = TreeNode(t, idb, "");
    }

    bool operator<(const TreeEdge& other) const {
        return std::make_tuple(from, to) < std::make_tuple(other.from, other.to);
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
    bool recordTerminationTree = true;

    // Current set of clauses learned by CDCL.
    std::vector<Clause> learnedClauses;

    // Number of UNSAT conflicts encountered during search.
    int numConflicts = 0;

    //
    // CDCL data structures.
    //

    // Trail of assigned literals along with their decision levels.
    std::vector<std::pair<Literal, int>> currTrail;
    // Current assignment, incorporating all assignments made on the current trail.
    Assignment currTrailAssmt;

public:
    Solver() {}

    Assignment getAssignment() {
        return _currAssignment;
    }

    void addTerminationEdge(TreeEdge e) {
        terminationTree.insert(e);
    }

    int getNumConflicts() {
        return numConflicts;
    }

    void printTerminationTree() {
        std::cout << "termination tree:" << std::endl;
        for (auto e : terminationTree) {
            std::cout << "\"" << e.from.getDOTId() << "\" -> \"" << e.to.getDOTId() << "\""
                      << std::endl;
        }
        std::cout << "===" << std::endl;
    }

    std::string getTerminationTreeDOT() {
        std::string outStr = "digraph G {\n";
        for (auto e : terminationTree) {
            std::string edgeStr = "";
            edgeStr += ("\"" + e.from.getDOTId() + "\"");
            edgeStr += (" -> ");
            edgeStr += ("\"" + e.to.getDOTId() + "\"");
            edgeStr += ";";
            outStr += (edgeStr + "\n");
        }
        outStr += "}\n";
        return outStr;
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

    CNF cdclUnitPropagation(CNF f,
                            int currDecisionLevel,
                            std::map<std::string, int>& antecedents,
                            std::map<std::string, int>& varDecisionLevels,
                            std::set<std::string>& varsAssignedAtCurrLevel) {
        while (f.hasUnitClause()) {
            // Also return the index of the clause that this variable is now unit in.
            // A literal l is unit in a clause (l1 \/ ... \/ lk \/ l) if all l1..lk have
            // been assigned false.
            Literal unitLit = f.getUnitLiteral();
            int unitClauseInd = f.getUnitClause();

            // Assign the unit literal's variable to truthify the unit clause
            // and simplify the formula based on this assignment.
            // currAssmt.set(unitLit.getVarName(), !unitLit.isNegated(),
            // currNode.getDecisionLevel());
            f = f.unitPropagate(unitLit);
            LOG(DEBUG) << "f assigned after unit prop round: " << f.toString();
            antecedents[unitLit.getVarName()] = unitClauseInd;

            // Save assignment to the trail.
            currTrail.push_back({unitLit, currDecisionLevel});

            // Mark as assigned at this decision level.
            varsAssignedAtCurrLevel.insert(unitLit.getVarName());
            varDecisionLevels[unitLit.getVarName()] = currDecisionLevel;

            currTrailAssmt.set(unitLit.getVarName(), !unitLit.isNegated());
        }
        return f;
    }

    /**
     * Analyze conflict, return a learned clause and a backjump level.
     */
    std::pair<Clause, int> cdclAnalyzeConflict(CNF f,
                                               Clause conflictClause,
                                               std::map<std::string, int> antecedents,
                                               std::map<std::string, int> varDecisionLevels,
                                               std::set<std::string> varsAssignedAtCurrLevel) {

        // If we have now encountered a conflict, we can now start
        // with the conflicting clause C, and pick the most recent
        // literal in it that was assigned, L, Then take the
        // antecedent of L, and compute resolve(C,L). This gives a
        // new clause C', which we can then continute the process
        // with. That is, find the most recent assigned literal
        // appearing in C', and then find its antecedent and do
        // resolution again. We should be able to find the most
        // recent assigned literal in a clause by walking backwards
        // through the trail (??)

        // Start with the conflicting clause.
        // auto emptyCInd = fassigned.getEmptyClause();
        LOG(DEBUG) << "cdcl analyzing conflict ";

        int backjumpLevel = -1;

        // Clause currClause = currCNF.getClause(emptyCInd);
        Clause currClause = conflictClause;
        // currCNF.getClause(emptyCInd);

        // Walk backwards through the trail.
        int ti = currTrail.size() - 1;
        while (ti > 0) {
            auto lit = currTrail[ti].first;
            auto litLevel = currTrail[ti].second;
            if (currClause.hasVariable(lit.getVarName())) {
                // This is the most recent var in trail that was assigned in this
                // clause.
                LOG(DEBUG) << "first assigned: " << lit.toString();

                // Then, figure out the antecedent clause for this variable.
                int anteInd = antecedents[lit.getVarName()];

                // Have reached a decision (non-forced) variable assignment in the
                // trail.
                if (anteInd < 0) {
                    LOG(DEBUG) << "Reached decision variable.";
                    break;
                }

                Clause ante = f.getClause(anteInd);
                LOG(DEBUG) << "curr conflict clause: " << currClause.toString();
                LOG(DEBUG) << "antecedent: " << ante.toString() << " (c" << anteInd << ")";

                // Now, resolve the current clause with the antecedent.
                Clause resolvent = Clause::resolve(ante, currClause, lit.getVarName());
                LOG(DEBUG) << "resolvent: " << resolvent.toString();
                currClause = resolvent;

                // TODO (4/25/22): Once we're done resolving, we need to add the
                // discovered conflict clause as a learned clause, and then backjump
                // appropriately.

                // Termination condition for resolution: stop when
                // the current resolvent contains exactly one
                // literal whose value was assigned at the current
                // decision level.

                // Check for the number of variables assigned at this decision level
                // that appear in the current conflict clause.
                int numVarsFromCurrLevel = 0;
                for (auto v : currClause.getVariableSet()) {
                    if (varsAssignedAtCurrLevel.find(v) != varsAssignedAtCurrLevel.end()) {
                        numVarsFromCurrLevel += 1;
                    }
                }
                LOG(DEBUG) << "num vars from curr decision level: " << numVarsFromCurrLevel;

                if (numVarsFromCurrLevel == 1) {
                    // Terminate.
                    Clause learnedClause = currClause;
                    learnedClauses.push_back(learnedClause);

                    // We will backtrack to the "assertion level",
                    // which is the second highest (i.e. second
                    // deepest) level of any variable that appears
                    // in the learned conflict clause.
                    int levelHighest = -1;
                    int levelSecondHighest = -1;
                    for (auto l : learnedClause.getLiterals()) {
                        // Determine the decision level of this variable.
                        auto level = varDecisionLevels[l.getVarName()];
                        // LOG(DEBUG) << "level: " << level;
                        if (level > levelHighest) {
                            levelSecondHighest = levelHighest;
                            levelHighest = level;
                        } else if (level > levelSecondHighest && level != levelHighest) {
                            levelSecondHighest = level;
                        }
                    }

                    if (levelSecondHighest == -1) {
                        levelSecondHighest = levelHighest;
                    }
                    backjumpLevel = levelSecondHighest;
                    // currCNF.appendClause(learnedClause);
                    return {learnedClause, backjumpLevel};
                }
            }
            ti--;
        }
        // Should never reach this return.
        return {Clause::makeTrueClause(), -1};
    }

    bool isSatCDCL(CNF f) {
        LOG(DEBUG) << "checking isSat CDCL:" << f.toString();

        std::vector<std::string> varList = f.getVariableList();

        // The set of tree nodes in the frontier i.e. discovered but not explored yet.
        std::vector<Context> frontier;

        // Sequence of currently made assignments, along with decision level.
        // std::vector<std::pair<Literal, int>> currTrail;
        // Assignment currTrailAssmt;

        std::map<std::string, bool> initVals;
        int initDecisionLevel = -1;
        // frontier.push_back(Context(f, varList.at(0), "root", 0, initVals, {},
        // initDecisionLevel));

        int currDecisionLevel = 0;

        // Current formula, with potentially learned clauses added over time.
        CNF currCNF = f;

        CNF currF = currCNF;
        CNF fassigned = currF;
        // currF.assign(currAssmt);
        LOG(DEBUG) << "f after assignment: " << fassigned.toString();

        std::set<std::string> unchosenVars = f.getVariableSet();

        bool backjumped = false;

        // while (true) {
        for (int j = 0; j < 8; j++) {

            // TODO: We need to run unit propagation after we've backjumped and start exploring
            // again.

            // Store a mapping from each assigned variable back to the
            // index of the clause in which it became unit i.e. store a
            // pointer back to its "antecedent" clause. That is, the
            // clause that "caused" it be assigned via unit resolution.
            std::map<std::string, int> antecedents;
            // std::vector<Literal> trail;
            // The decision levels of each assigned variable.
            std::map<std::string, int> varDecisionLevels;

            // Set of variables that are assigned at the current decision
            // level, either by an explicit decision, or as an assignment
            // derived via unit propagation.
            std::set<std::string> varsAssignedAtCurrLevel;
            // The parent var is the one that has been assigned at this decision level.
            // varsAssignedAtCurrLevel.insert(currNode._parentVar);

            // Extend the trail with a new variable assignment, and update the current
            // assignment.

            // If we just backjumped here, don't assign a new variable yet. First, apply unit
            // propagation.
            LOG(DEBUG) << "---> current f: " << currCNF.toString();
            LOG(DEBUG) << "current decision level: " << currDecisionLevel;
            if (!backjumped) {
                LOG(DEBUG) << "## choosing new variable assignment";
                auto it = unchosenVars.begin();
                std::string chosenVar = *it;
                unchosenVars.erase(it);
                auto newLit = Literal(chosenVar);
                currTrail.push_back({newLit, currDecisionLevel});
                currTrailAssmt.set(chosenVar, !newLit.isNegated());
                varDecisionLevels[chosenVar] = currDecisionLevel;

                LOG(DEBUG) << "extending trail with " << newLit.toString();
                LOG(DEBUG) << "current trail:";
                for (auto l : currTrail) {
                    LOG(DEBUG) << l.first.toString() << " (l" << l.second << ")";
                }
            } else {
                LOG(DEBUG) << "just backjumped, so not choosing new variable yet.";
            }

            // Update the current set of variables assigned at this decision level.
            for (auto l : currTrail) {
                if (l.second == currDecisionLevel) {
                    varsAssignedAtCurrLevel.insert(l.first.getVarName());
                }
            }

            CNF currF = currCNF;
            CNF fassigned = currF.assign(currTrailAssmt);
            LOG(DEBUG) << "f after assignment: " << fassigned.toString();
            if (fassigned.isSatisfied()) {
                // If we determined that the formula is necessarily satisfiable
                // under the current partial assignment, then we can fill in
                // arbitrary values for the remaining, unassigned variables.
                for (auto v : varList) {
                    // Assign a value to the currently unassigned variable.
                    if (!currTrailAssmt.hasVar(v)) {
                        bool arbitraryVal = true;
                        currTrailAssmt.set(v, arbitraryVal);
                    }
                }
                _currAssignment = currTrailAssmt;
                return true;
            }

            // Initialize the antecedent graph with the existing variable assignments.
            auto avars = currTrailAssmt.getVals();
            for (auto it = avars.begin(); it != avars.end(); it++) {
                antecedents[it->first] = -1;  // no predecessor.
                // trail.push_back(Literal(it->first, !it->second));
                LOG(DEBUG) << "trail element " << it->first << "=" << it->second;
                // varDecisionLevels[it->first] = currAssmt.getDecisionLevel(it->first);
            }

            // Apply unit propagation.
            LOG(DEBUG) << "running unit propagation for CDCL";
            fassigned = cdclUnitPropagation(fassigned,
                                            currDecisionLevel,
                                            antecedents,
                                            varDecisionLevels,
                                            varsAssignedAtCurrLevel);

            LOG(DEBUG) << "f after unit prop: " << fassigned.toString();
            LOG(DEBUG) << "curr assignment after unit prop: " << currTrailAssmt.toString();

            // Conflict analysis. Learn a clause and determine backjump level.
            int backjumpLevel = -1;
            if (fassigned.hasEmptyClause()) {

                numConflicts++;

                int emptyCInd = fassigned.getEmptyClause();
                Clause currClause = currCNF.getClause(emptyCInd);
                LOG(DEBUG) << "!! conflict encountered, conflicting clause index: " << emptyCInd;

                auto ret = cdclAnalyzeConflict(
                    currCNF, currClause, antecedents, varDecisionLevels, varsAssignedAtCurrLevel);
                Clause learnedClause = ret.first;
                backjumpLevel = ret.second;
                if (backjumpLevel >= 0) {
                    currCNF.appendClause(learnedClause);
                }

                LOG(DEBUG) << "**** learned clause: " << learnedClause.toString();
                LOG(DEBUG) << "**** backjump level: " << backjumpLevel;

                if (backjumpLevel == -1) {
                    // Pop all items off the trail and backjump to restart the whole search.
                    while (currTrail.size() > 0) {
                        currTrail.pop_back();
                    }
                    backjumped = true;
                    currDecisionLevel = backjumpLevel;
                    currTrailAssmt = Assignment();
                    unchosenVars = f.getVariableSet();
                    continue;
                }
            }

            // Backjump to lower decision level if necessary.
            if (backjumpLevel >= 0) {
                LOG(DEBUG) << "Backjumping to level " << backjumpLevel;
                LOG(DEBUG) << "current trail stack:";
                for (auto l : currTrail) {
                    LOG(DEBUG) << l.first.toString() << " (l" << l.second << ")";
                }

                LOG(DEBUG) << "popping off items on trail stack. ";

                // TODO: I don't think this backjump level calculation is correct.
                // int frontierBackjumpLevel = (backjumpLevel + 1);
                int frontierBackjumpLevel = backjumpLevel;
                while (currTrail.back().second > frontierBackjumpLevel && currTrail.size() > 0) {
                    // while (frontier.back().getDecisionLevel() > frontierBackjumpLevel &&
                    //    frontier.size() > 0) {
                    // LOG(DEBUG) << "popping " << frontier.back()._assmt.toString();

                    LOG(DEBUG) << "popping " << currTrail.back().first.toString() << ", (l"
                               << currTrail.back().second << ")";

                    // Clear this assignment.
                    // assmtToRestartFrom.unset(frontier.back()._currVar);
                    currTrailAssmt.unset(currTrail.back().first.getVarName());

                    // Add this variable back into the set of unchosen vars.
                    unchosenVars.insert(currTrail.back().first.getVarName());

                    // frontier.pop_back();
                    currTrail.pop_back();
                    LOG(DEBUG) << "newAssmt: " << currTrailAssmt.toString();
                }
                backjumped = true;
                currDecisionLevel = backjumpLevel;

            } else {
                // If we didn't backjump, proceed to the next decision level.
                currDecisionLevel++;
                backjumped = false;
            }
        }

        return false;
    }

    // // Explore all possible assignments in a depth first manner.
    // while (frontier.size() > 0) {
    //     Context currNode = frontier.back();
    //     frontier.pop_back();

    //     LOG(DEBUG) << "### currVar to assign next: '" << currNode._currVar
    //                << "', varInd=" << currNode._currVarInd
    //                << ", parentVar: " << currNode._parentVar << " ###";
    //     LOG(DEBUG) << "curr decision level: " << currNode.getDecisionLevel();
    //     LOG(DEBUG) << "parentVar: " << currNode._parentVar;
    //     LOG(DEBUG) << "curr assignment: " << currNode._assmt.toString();
    //     // LOG(DEBUG) << "par assignment: " << currNode._parentAssmt.toString();

    //     // Reduce based on current assignment.
    //     auto currAssmt = currNode._assmt;
    //     CNF currF = currCNF;
    //     CNF fassigned = currF.assign(currAssmt);
    //     LOG(DEBUG) << "f after assignment: " << fassigned.toString();

    //     // Pure literal elimination.
    //     // TODO: Not implemented currently but may add it.

    //     if (fassigned.isSatisfied()) {
    //         // If we determined that the formula is necessarily satisfiable
    //         // under the current partial assignment, then we can fill in
    //         // arbitrary values for the remaining, unassigned variables.
    //         for (auto v : varList) {
    //             // Assign a value to the currently unassigned variable.
    //             if (!currAssmt.hasVar(v)) {
    //                 bool arbitraryVal = true;
    //                 currAssmt.set(v, arbitraryVal);
    //             }
    //         }
    //         _currAssignment = currAssmt;
    //         return true;
    //     }

    //     // Close the formula under unit resolution, if enabled.
    //     int backjumpLevel = -1;
    //     if (enableUnitPropagation) {
    //         LOG(DEBUG) << "running unit propagation";

    //         // Store a mapping from each assigned variable back to the
    //         // index of the clause in which it became unit i.e. store a
    //         // pointer back to its "antecedent" clause. That is, the
    //         // clause that "caused" it be assigned via unit resolution.
    //         std::map<std::string, int> antecedents;
    //         std::vector<Literal> trail;
    //         // The decision levels of each assigned variable.
    //         std::map<std::string, int> varDecisionLevels;

    //         // Set of variables that are assigned at the current decision
    //         // level, either by an explicit decision, or as an assignment
    //         // derived via unit propagation.
    //         std::set<std::string> varsAssignedAtCurrLevel;
    //         // The parent var is the one that has been assigned at this decision level.
    //         varsAssignedAtCurrLevel.insert(currNode._parentVar);

    //         // Initialize the antecedent graph with the existing variable assignments.
    //         auto avars = currAssmt.getVals();
    //         for (auto it = avars.begin(); it != avars.end(); it++) {
    //             antecedents[it->first] = -1;  // no predecessor.
    //             trail.push_back(Literal(it->first, !it->second));
    //             LOG(DEBUG) << "trail element " << it->first << "=" << it->second;

    //             varDecisionLevels[it->first] = currAssmt.getDecisionLevel(it->first);
    //         }

    //         while (fassigned.hasUnitClause()) {
    //             // Also return the index of the clause that this variable is now unit in.
    //             // A literal l is unit in a clause (l1 \/ ... \/ lk \/ l) if all l1..lk have
    //             // been assigned false.
    //             auto unitLit = fassigned.getUnitLiteral();
    //             auto unitClauseInd = fassigned.getUnitClause();

    //             // Assign the unit literal's variable to truthify the unit clause
    //             // and simplify the formula based on this assignment.
    //             currAssmt.set(
    //                 unitLit.getVarName(), !unitLit.isNegated(), currNode.getDecisionLevel());
    //             fassigned = fassigned.unitPropagate(unitLit);
    //             LOG(DEBUG) << "f assigned after unit prop round: " << fassigned.toString();
    //             antecedents[unitLit.getVarName()] = unitClauseInd;

    //             // Save assignment to the trail.
    //             trail.push_back(unitLit);

    //             // Mark as assigned at this decision level.
    //             varsAssignedAtCurrLevel.insert(unitLit.getVarName());
    //             varDecisionLevels[unitLit.getVarName()] = currNode.getDecisionLevel();
    //         }
    //         LOG(DEBUG) << "current decision level: " << currNode.getDecisionLevel();
    //         LOG(DEBUG) << "f after unit prop: " << fassigned.toString();
    //         LOG(DEBUG) << "curr assignment after unit prop: " << currAssmt.toString();

    //         for (auto v : varsAssignedAtCurrLevel) {
    //             LOG(DEBUG) << "after unit prop, var assigned at curr level: " << v;
    //         }

    //         // TODO: Finish fleshing out antecedent graph and its use for deriving conflict
    //         // clauses.
    //         // std::cout << "antecedent graph:" << std::endl;
    //         for (auto it = antecedents.begin(); it != antecedents.end(); it++) {
    //             // std::cout << it->first << " -> " << it->second << std::endl;

    //             std::string litVarName = it->first;
    //             int clauseInd = it->second;

    //             // Record predecessor variable assignments for each node.
    //             // TODO: Decide how to store this.
    //             if (clauseInd >= 0) {
    //                 auto unitClause = currCNF.getClause(clauseInd);
    //                 std::vector<std::string> anteVars;
    //                 for (auto l : unitClause.getLiterals()) {
    //                     if (l.getVarName() != litVarName) {
    //                         anteVars.push_back(l.getVarName());
    //                     }
    //                 };
    //             }
    //         }

    //         LOG(DEBUG) << "trail:";
    //         for (auto l : trail) {
    //             LOG(DEBUG) << l.toString();
    //         }

    //         // If a contradiction has been derived, apply conflict analysis.
    //         if (fassigned.hasEmptyClause()) {

    //             numConflicts++;

    //             // If we have now encountered a conflict, we can now start
    //             // with the conflicting clause C, and pick the most recent
    //             // literal in it that was assigned, L, Then take the
    //             // antecedent of L, and compute resolve(C,L). This gives a
    //             // new clause C', which we can then continute the process
    //             // with. That is, find the most recent assigned literal
    //             // appearing in C', and then find its antecedent and do
    //             // resolution again. We should be able to find the most
    //             // recent assigned literal in a clause by walking backwards
    //             // through the trail (??)

    //             // Start with the conflicting clause.
    //             auto emptyCInd = fassigned.getEmptyClause();
    //             Clause currClause = currCNF.getClause(emptyCInd);
    //             LOG(DEBUG) << "!! conflict encountered, conflicting clause index: "
    //                        << emptyCInd;

    //             // Walk backwards through the trail.
    //             int ti = trail.size() - 1;
    //             while (ti > 0) {
    //                 auto lit = trail[ti];
    //                 if (currClause.hasVariable(lit.getVarName())) {
    //                     // This is the most recent var in trail that was assigned in this
    //                     // clause.
    //                     LOG(DEBUG) << "first assigned: " << lit.toString();

    //                     // Then, figure out the antecedent clause for this variable.
    //                     int anteInd = antecedents[lit.getVarName()];

    //                     // Have reached a decision (non-forced) variable assignment in the
    //                     // trail.
    //                     if (anteInd < 0) {
    //                         LOG(DEBUG) << "Reached decision variable.";
    //                         break;
    //                     }

    //                     Clause ante = currCNF.getClause(anteInd);
    //                     LOG(DEBUG) << "curr conflict clause: " << currClause.toString();
    //                     LOG(DEBUG)
    //                         << "antecedent: " << ante.toString() << " (c" << anteInd << ")";

    //                     // Now, resolve the current clause with the antecedent.
    //                     Clause resolvent = Clause::resolve(ante, currClause,
    //                     lit.getVarName()); LOG(DEBUG) << "resolvent: " <<
    //                     resolvent.toString(); currClause = resolvent;

    //                     // TODO (4/25/22): Once we're done resolving, we need to add the
    //                     // discovered conflict clause as a learned clause, and then backjump
    //                     // appropriately.

    //                     // Termination condition for resolution: stop when
    //                     // the current resolvent contains exactly one
    //                     // literal whose value was assigned at the current
    //                     // decision level.

    //                     // Check for the number of variables assigned at this decision level
    //                     // that appear in the current conflict clause.
    //                     int numVarsFromCurrLevel = 0;
    //                     for (auto v : currClause.getVariableSet()) {
    //                         if (varsAssignedAtCurrLevel.find(v) !=
    //                             varsAssignedAtCurrLevel.end()) {
    //                             numVarsFromCurrLevel += 1;
    //                         }
    //                     }
    //                     LOG(DEBUG)
    //                         << "num vars from curr decision level: " << numVarsFromCurrLevel;

    //                     if (numVarsFromCurrLevel == 1) {
    //                         // Terminate.
    //                         Clause learnedClause = currClause;
    //                         learnedClauses.push_back(learnedClause);
    //                         LOG(DEBUG) << "-> learned clause: " << learnedClause.toString();

    //                         // We will backtrack to the "assertion level",
    //                         // which is the second highest (i.e. second
    //                         // deepest) level of any variable that appears
    //                         // in the learned conflict clause.
    //                         int levelHighest = -1;
    //                         int levelSecondHighest = -1;
    //                         for (auto l : learnedClause.getLiterals()) {
    //                             // Determine the decision level of this variable.
    //                             auto level = varDecisionLevels[l.getVarName()];
    //                             // LOG(DEBUG) << "level: " << level;
    //                             if (level > levelHighest) {
    //                                 levelSecondHighest = levelHighest;
    //                                 levelHighest = level;
    //                             } else if (level > levelSecondHighest &&
    //                                        level != levelHighest) {
    //                                 levelSecondHighest = level;
    //                             }
    //                         }

    //                         if (levelSecondHighest == -1) {
    //                             levelSecondHighest = levelHighest;
    //                         }
    //                         backjumpLevel = levelSecondHighest;
    //                         currCNF.appendClause(learnedClause);
    //                         break;
    //                     }
    //                 }
    //                 ti--;
    //             }
    //         }

    //         LOG(DEBUG) << "backjump level: " << backjumpLevel;
    //     }


    //     // Record information for termination tree if enabled.
    //     // Mostly for debugging/visualization.
    //     if (recordTerminationTree) {
    //         auto localCurrAssmt = currNode._assmt;
    //         auto from = TreeNode(currNode._parentVar,
    //                              currNode._parentAssmt.toStringCompact(),
    //                              currNode._f.toString());
    //         auto to = TreeNode(
    //             currNode._currVar, currNode._assmt.toStringCompact(), fassigned.toString());
    //         TreeEdge e = TreeEdge(from, to);
    //         terminationTree.insert(e);
    //     }


    //     //
    //     // Backjump to lower decision level if we ran conflict analysis!
    //     //
    //     if (backjumpLevel >= 0) {
    //         LOG(DEBUG) << "current frontier stack:";
    //         for (Context e : frontier) {
    //             LOG(DEBUG) << "  " << e._assmt.toString();
    //         }

    //         LOG(DEBUG) << "popping current levels on stack. ";

    //         // Context with the assignments we backtrack past removed, that
    //         // we will restart our search from.
    //         Assignment assmtToRestartFrom = currAssmt;

    //         // The stack records 'Context' objects, so we offset the backjump level by 1 to
    //         // account for this.

    //         // TODO: I don't think this backjump level calculation is correct.
    //         int frontierBackjumpLevel = (backjumpLevel + 1);
    //         while (frontier.back().getDecisionLevel() > frontierBackjumpLevel &&
    //                frontier.size() > 0) {
    //             LOG(DEBUG) << "popping " << frontier.back()._assmt.toString();

    //             // Clear this assignment.
    //             assmtToRestartFrom.unset(frontier.back()._currVar);

    //             frontier.pop_back();
    //         }
    //         // Restart the search process after backjumping.
    //         LOG(DEBUG) << "Backjumping in search. ";

    //         //
    //         // TODO: Need to think more carefully about what needs to be done here.
    //         //
    //         // auto restartVar = varList.at(backjumpLevel);
    //         // LOG(DEBUG) << "Backjumping to variable: " << restartVar << ", assignment: " <<
    //         // assmtToRestartFrom.toString(); Context ctxToRestartFrom(currCNF, restartVar,
    //         "",
    //         // backjumpLevel, assmtToRestartFrom, {}, backjumpLevel);
    //         // frontier.push_back(ctxToRestartFrom);

    //         continue;
    //     }

    //     // if (fassigned.isEmpty()) {
    //     if (fassigned.isSatisfied()) {
    //         // If we determined that the formula is necessarily satisfiable
    //         // under the current partial assignment, then we can fill in
    //         // arbitrary values for the remaining, unassigned variables.
    //         for (auto v : varList) {
    //             // Assign a value to the currently unassigned variable.
    //             if (!currAssmt.hasVar(v)) {
    //                 bool arbitraryVal = true;
    //                 currAssmt.set(v, arbitraryVal);
    //             }
    //         }
    //         _currAssignment = currAssmt;
    //         return true;
    //     }
    //     if (fassigned.hasEmptyClause()) {
    //         // Current assignment is necessarily UNSAT, so no need to
    //         // explore further down this branch.
    //         LOG(DEBUG) << "fassigned has empty clause. (** CONFLICT **)";
    //         lastNode = currNode;
    //         continue;
    //     }

    //     int varNextInd = currNode._currVarInd + 1;
    //     LOG(DEBUG) << "varNextInd: " << varNextInd;

    //     // Check if we reached the last variable i.e. a leaf.
    //     if (varNextInd <= varList.size()) {

    //         // Explore each possible true/false assignment for this variable.
    //         Assignment tAssign(currAssmt);
    //         Assignment fAssign(currAssmt);

    //         tAssign.set(currNode._currVar, true, currNode.getDecisionLevel() + 1);
    //         fAssign.set(currNode._currVar, false, currNode.getDecisionLevel() + 1);

    //         // If we have reached the last variable, then pass in a dummy
    //         // 'leaf' variable node for the next level in the search tree.
    //         // This value isn't explicitly used, but it makes it convenient
    //         // so that we can traverse one more level deep in the search
    //         // tree.
    //         auto nextVar = varNextInd == varList.size() ? "LEAF" : varList.at(varNextInd);

    //         std::string parentVar = currNode._currVar;
    //         Context tctx(fassigned,
    //                      nextVar,
    //                      parentVar,
    //                      varNextInd,
    //                      tAssign,
    //                      currNode._assmt,
    //                      currNode.getDecisionLevel() + 1);
    //         Context fctx(fassigned,
    //                      nextVar,
    //                      parentVar,
    //                      varNextInd,
    //                      fAssign,
    //                      currNode._assmt,
    //                      currNode.getDecisionLevel() + 1);

    //         // TODO: Order of true/false exploration could be chosen differently or
    //         dynamically. LOG(DEBUG) << "pushing " << fAssign.toString() << " nextvar=" <<
    //         nextVar
    //                    << " onto search frontier";
    //         frontier.push_back(fctx);
    //         LOG(DEBUG) << "pushing " << tAssign.toString() << " nextvar=" << nextVar
    //                    << " onto search frontier";
    //         frontier.push_back(tctx);
    //     }
    // }


    bool isSat(CNF f) {
        LOG(DEBUG) << "checking isSat:" << f.toString();

        std::vector<std::string> varList = f.getVariableList();

        // The set of tree nodes in the frontier i.e. discovered but not explored yet.
        std::vector<Context> frontier;

        std::map<std::string, bool> initVals;
        int initDecisionLevel = -1;
        frontier.push_back(Context(f, varList.at(0), "root", 0, initVals, {}, initDecisionLevel));

        Context lastNode = Context();

        // Current formula, with potentially learned clauses added over time.
        CNF currCNF = f;

        // Explore all possible assignments in a depth first manner.
        while (frontier.size() > 0) {
            Context currNode = frontier.back();
            frontier.pop_back();

            LOG(DEBUG) << "### currVar to assign next: '" << currNode._currVar
                       << "', varInd=" << currNode._currVarInd
                       << ", parentVar: " << currNode._parentVar << " ###";
            LOG(DEBUG) << "curr decision level: " << currNode.getDecisionLevel();
            LOG(DEBUG) << "parentVar: " << currNode._parentVar;
            LOG(DEBUG) << "curr assignment: " << currNode._assmt.toString();
            // LOG(DEBUG) << "par assignment: " << currNode._parentAssmt.toString();

            // Reduce based on current assignment.
            auto currAssmt = currNode._assmt;
            CNF currF = currCNF;
            CNF fassigned = currF.assign(currAssmt);
            LOG(DEBUG) << "f after assignment: " << fassigned.toString();

            // Pure literal elimination.
            // TODO: Not implemented currently but may add it.

            if (fassigned.isSatisfied()) {
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

            // Close the formula under unit resolution, if enabled.
            int backjumpLevel = -1;
            if (enableUnitPropagation) {
                LOG(DEBUG) << "running unit propagation";

                // Store a mapping from each assigned variable back to the
                // index of the clause in which it became unit i.e. store a
                // pointer back to its "antecedent" clause. That is, the
                // clause that "caused" it be assigned via unit resolution.
                std::map<std::string, int> antecedents;
                std::vector<Literal> trail;
                // The decision levels of each assigned variable.
                std::map<std::string, int> varDecisionLevels;

                // Set of variables that are assigned at the current decision
                // level, either by an explicit decision, or as an assignment
                // derived via unit propagation.
                std::set<std::string> varsAssignedAtCurrLevel;
                // The parent var is the one that has been assigned at this decision level.
                varsAssignedAtCurrLevel.insert(currNode._parentVar);

                // Initialize the antecedent graph with the existing variable assignments.
                auto avars = currAssmt.getVals();
                for (auto it = avars.begin(); it != avars.end(); it++) {
                    antecedents[it->first] = -1;  // no predecessor.
                    trail.push_back(Literal(it->first, !it->second));
                    LOG(DEBUG) << "trail element " << it->first << "=" << it->second;

                    varDecisionLevels[it->first] = currAssmt.getDecisionLevel(it->first);
                }

                while (fassigned.hasUnitClause()) {
                    // Also return the index of the clause that this variable is now unit in.
                    // A literal l is unit in a clause (l1 \/ ... \/ lk \/ l) if all l1..lk have
                    // been assigned false.
                    auto unitLit = fassigned.getUnitLiteral();
                    auto unitClauseInd = fassigned.getUnitClause();

                    // Assign the unit literal's variable to truthify the unit clause
                    // and simplify the formula based on this assignment.
                    currAssmt.set(
                        unitLit.getVarName(), !unitLit.isNegated(), currNode.getDecisionLevel());
                    fassigned = fassigned.unitPropagate(unitLit);
                    LOG(DEBUG) << "f assigned after unit prop round: " << fassigned.toString();
                    antecedents[unitLit.getVarName()] = unitClauseInd;

                    // Save assignment to the trail.
                    trail.push_back(unitLit);

                    // Mark as assigned at this decision level.
                    varsAssignedAtCurrLevel.insert(unitLit.getVarName());
                    varDecisionLevels[unitLit.getVarName()] = currNode.getDecisionLevel();
                }
                LOG(DEBUG) << "current decision level: " << currNode.getDecisionLevel();
                LOG(DEBUG) << "f after unit prop: " << fassigned.toString();
                LOG(DEBUG) << "curr assignment after unit prop: " << currAssmt.toString();

                for (auto v : varsAssignedAtCurrLevel) {
                    LOG(DEBUG) << "after unit prop, var assigned at curr level: " << v;
                }

                // TODO: Finish fleshing out antecedent graph and its use for deriving conflict
                // clauses.
                // std::cout << "antecedent graph:" << std::endl;
                for (auto it = antecedents.begin(); it != antecedents.end(); it++) {
                    // std::cout << it->first << " -> " << it->second << std::endl;

                    std::string litVarName = it->first;
                    int clauseInd = it->second;

                    // Record predecessor variable assignments for each node.
                    // TODO: Decide how to store this.
                    if (clauseInd >= 0) {
                        auto unitClause = currCNF.getClause(clauseInd);
                        std::vector<std::string> anteVars;
                        for (auto l : unitClause.getLiterals()) {
                            if (l.getVarName() != litVarName) {
                                anteVars.push_back(l.getVarName());
                            }
                        };
                    }
                }

                LOG(DEBUG) << "trail:";
                for (auto l : trail) {
                    LOG(DEBUG) << l.toString();
                }

                // If a contradiction has been derived, apply conflict analysis.
                if (fassigned.hasEmptyClause()) {

                    numConflicts++;

                    // If we have now encountered a conflict, we can now start
                    // with the conflicting clause C, and pick the most recent
                    // literal in it that was assigned, L, Then take the
                    // antecedent of L, and compute resolve(C,L). This gives a
                    // new clause C', which we can then continute the process
                    // with. That is, find the most recent assigned literal
                    // appearing in C', and then find its antecedent and do
                    // resolution again. We should be able to find the most
                    // recent assigned literal in a clause by walking backwards
                    // through the trail (??)

                    // Start with the conflicting clause.
                    auto emptyCInd = fassigned.getEmptyClause();
                    Clause currClause = currCNF.getClause(emptyCInd);
                    LOG(DEBUG) << "!! conflict encountered, conflicting clause index: "
                               << emptyCInd;

                    // Walk backwards through the trail.
                    int ti = trail.size() - 1;
                    while (ti > 0) {
                        auto lit = trail[ti];
                        if (currClause.hasVariable(lit.getVarName())) {
                            // This is the most recent var in trail that was assigned in this
                            // clause.
                            LOG(DEBUG) << "first assigned: " << lit.toString();

                            // Then, figure out the antecedent clause for this variable.
                            int anteInd = antecedents[lit.getVarName()];

                            // Have reached a decision (non-forced) variable assignment in the
                            // trail.
                            if (anteInd < 0) {
                                LOG(DEBUG) << "Reached decision variable.";
                                break;
                            }

                            Clause ante = currCNF.getClause(anteInd);
                            LOG(DEBUG) << "curr conflict clause: " << currClause.toString();
                            LOG(DEBUG)
                                << "antecedent: " << ante.toString() << " (c" << anteInd << ")";

                            // Now, resolve the current clause with the antecedent.
                            Clause resolvent = Clause::resolve(ante, currClause, lit.getVarName());
                            LOG(DEBUG) << "resolvent: " << resolvent.toString();
                            currClause = resolvent;

                            // TODO (4/25/22): Once we're done resolving, we need to add the
                            // discovered conflict clause as a learned clause, and then backjump
                            // appropriately.

                            // Termination condition for resolution: stop when
                            // the current resolvent contains exactly one
                            // literal whose value was assigned at the current
                            // decision level.

                            // Check for the number of variables assigned at this decision level
                            // that appear in the current conflict clause.
                            int numVarsFromCurrLevel = 0;
                            for (auto v : currClause.getVariableSet()) {
                                if (varsAssignedAtCurrLevel.find(v) !=
                                    varsAssignedAtCurrLevel.end()) {
                                    numVarsFromCurrLevel += 1;
                                }
                            }
                            LOG(DEBUG)
                                << "num vars from curr decision level: " << numVarsFromCurrLevel;

                            if (numVarsFromCurrLevel == 1) {
                                // Terminate.
                                Clause learnedClause = currClause;
                                learnedClauses.push_back(learnedClause);
                                LOG(DEBUG) << "-> learned clause: " << learnedClause.toString();

                                // We will backtrack to the "assertion level",
                                // which is the second highest (i.e. second
                                // deepest) level of any variable that appears
                                // in the learned conflict clause.
                                int levelHighest = -1;
                                int levelSecondHighest = -1;
                                for (auto l : learnedClause.getLiterals()) {
                                    // Determine the decision level of this variable.
                                    auto level = varDecisionLevels[l.getVarName()];
                                    // LOG(DEBUG) << "level: " << level;
                                    if (level > levelHighest) {
                                        levelSecondHighest = levelHighest;
                                        levelHighest = level;
                                    } else if (level > levelSecondHighest &&
                                               level != levelHighest) {
                                        levelSecondHighest = level;
                                    }
                                }

                                if (levelSecondHighest == -1) {
                                    levelSecondHighest = levelHighest;
                                }
                                backjumpLevel = levelSecondHighest;
                                currCNF.appendClause(learnedClause);
                                break;
                            }
                        }
                        ti--;
                    }
                }

                LOG(DEBUG) << "backjump level: " << backjumpLevel;
            }


            // Record information for termination tree if enabled.
            // Mostly for debugging/visualization.
            if (recordTerminationTree) {
                auto localCurrAssmt = currNode._assmt;
                auto from = TreeNode(currNode._parentVar,
                                     currNode._parentAssmt.toStringCompact(),
                                     currNode._f.toString());
                auto to = TreeNode(
                    currNode._currVar, currNode._assmt.toStringCompact(), fassigned.toString());
                TreeEdge e = TreeEdge(from, to);
                terminationTree.insert(e);
            }


            //
            // Backjump to lower decision level if we ran conflict analysis!
            //
            if (backjumpLevel >= 0) {
                LOG(DEBUG) << "current frontier stack:";
                for (Context e : frontier) {
                    LOG(DEBUG) << "  " << e._assmt.toString();
                }

                LOG(DEBUG) << "popping current levels on stack. ";

                // Context with the assignments we backtrack past removed, that
                // we will restart our search from.
                Assignment assmtToRestartFrom = currAssmt;

                // The stack records 'Context' objects, so we offset the backjump level by 1 to
                // account for this.

                // TODO: I don't think this backjump level calculation is correct.
                int frontierBackjumpLevel = (backjumpLevel + 1);
                while (frontier.back().getDecisionLevel() > frontierBackjumpLevel &&
                       frontier.size() > 0) {
                    LOG(DEBUG) << "popping " << frontier.back()._assmt.toString();

                    // Clear this assignment.
                    assmtToRestartFrom.unset(frontier.back()._currVar);

                    frontier.pop_back();
                }
                // Restart the search process after backjumping.
                LOG(DEBUG) << "Backjumping in search. ";

                //
                // TODO: Need to think more carefully about what needs to be done here.
                //
                // auto restartVar = varList.at(backjumpLevel);
                // LOG(DEBUG) << "Backjumping to variable: " << restartVar << ", assignment: "
                // << assmtToRestartFrom.toString(); Context ctxToRestartFrom(currCNF,
                // restartVar, "", backjumpLevel, assmtToRestartFrom, {}, backjumpLevel);
                // frontier.push_back(ctxToRestartFrom);

                continue;
            }

            // if (fassigned.isEmpty()) {
            if (fassigned.isSatisfied()) {
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
                LOG(DEBUG) << "fassigned has empty clause. (** CONFLICT **)";
                lastNode = currNode;
                continue;
            }

            int varNextInd = currNode._currVarInd + 1;
            LOG(DEBUG) << "varNextInd: " << varNextInd;

            // Check if we reached the last variable i.e. a leaf.
            if (varNextInd <= varList.size()) {

                // Explore each possible true/false assignment for this variable.
                Assignment tAssign(currAssmt);
                Assignment fAssign(currAssmt);

                tAssign.set(currNode._currVar, true, currNode.getDecisionLevel() + 1);
                fAssign.set(currNode._currVar, false, currNode.getDecisionLevel() + 1);

                // If we have reached the last variable, then pass in a dummy
                // 'leaf' variable node for the next level in the search tree.
                // This value isn't explicitly used, but it makes it convenient
                // so that we can traverse one more level deep in the search
                // tree.
                auto nextVar = varNextInd == varList.size() ? "LEAF" : varList.at(varNextInd);

                std::string parentVar = currNode._currVar;
                Context tctx(fassigned,
                             nextVar,
                             parentVar,
                             varNextInd,
                             tAssign,
                             currNode._assmt,
                             currNode.getDecisionLevel() + 1);
                Context fctx(fassigned,
                             nextVar,
                             parentVar,
                             varNextInd,
                             fAssign,
                             currNode._assmt,
                             currNode.getDecisionLevel() + 1);

                // TODO: Order of true/false exploration could be chosen differently or
                // dynamically.
                LOG(DEBUG) << "pushing " << fAssign.toString() << " nextvar=" << nextVar
                           << " onto search frontier";
                frontier.push_back(fctx);
                LOG(DEBUG) << "pushing " << tAssign.toString() << " nextvar=" << nextVar
                           << " onto search frontier";
                frontier.push_back(tctx);
            }

            lastNode = currNode;
        }

        return false;
    }
};
