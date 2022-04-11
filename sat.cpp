#include <iostream>
#include <vector>

/**
 * A literal is a boolean variable that may be either negated or non-negated.
 */
class Literal{
    private:
        // The name of the boolean variable.
        std::string _name;
        // Is this a negated variable symbol.
        bool _neg;
    public:
        /**
         * Construct a new Literal, with specified negation.
         */
        Literal(std::string name, bool neg){
            _name = name;
            _neg = neg;
        }

        /**
         * Construct a literal given as a string [<neg>]<name>, where <neg> is
         * an optional "~" character. For example "~x" constructs the negation
         * of variable x, and "x" constructs the non-negated variable x.
         */
        Literal(std::string name){
            if(name.at(0) == '~'){
                _name = name.substr(1);
                _neg = true;
            } else{
                _name = name;
                _neg = false;
            }
        }

        std::string toString(){
            auto prefix = _neg ? "~" : "";
            return prefix + _name;
        }
};


/**
 * A clause is a set of literals L, representing the logical disjunction of the
 * literals in L.
 */
class Clause{
    private:
        std::vector<Literal> _elems;
    public:
        Clause(std::vector<Literal> elems){
            _elems = elems;
        }

        std::string toString(){
            std::string outStr;
            for(auto e : _elems){
                outStr += e.toString();
                outStr += " | ";
            }
            // Remove the extra disjunction symbol.
            return outStr.substr(0,outStr.length()-3);
        }
};

/**
 * A boolean formula in conjunctive normal form (CNF) is a conjunction of clauses.
 */
class CNF{
    private:
        std::vector<Clause> _clauses;
    public:
        CNF(std::vector<Clause> clauses){
            _clauses = clauses;
        }

        std::string toString(){
            std::string outStr;
            for(auto e : _clauses){
                outStr += "(" + e.toString() + ")";
                outStr += " & ";
            }
            // Remove the extra conjunction symbol.
            return outStr.substr(0,outStr.length()-2);            
        }
};


int main(int argc, char const *argv[])
{
    auto l1 = Literal("~x");
    auto l2 = Literal("y");
    std::cout << l1.toString() << std::endl;
    std::cout << l2.toString() << std::endl;

    auto c1 = Clause({l1,l2});
    auto c2 = Clause({l1});
    std::cout << c1.toString() << std::endl;
    std::cout << c2.toString() << std::endl;

    auto f1 = CNF({c1,c2});
    std::cout << f1.toString() << std::endl;

    return 0;
}
