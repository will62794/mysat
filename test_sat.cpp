#include "sat.cpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
using namespace std::chrono;

/**
 * Test that the output of the fast/optimized SAT solver implementation matches
 * that of the naive, brute force method for randomly generated CNF
 * formulas.
 */
void testConformanceRandomCNF(int niters, int nclauses, int nvars, int clause_size) {
    srand(33);
    for (int k = 0; k < niters; k++) {
        Solver solver = Solver();
        auto randf = CNF::randomCNF(nclauses, nvars, clause_size);
        std::cout << "- random CNF: " << randf.toString() << std::endl;
        auto retOracle = solver.isSatBruteForce(randf);
        auto retImpl = solver.isSat(randf);
        assert(retOracle == retImpl);

        // If the formula was satisfiable, test that the discovered satisfying
        // assignment is correct.
        if (retImpl) {
            assert(randf.eval(solver.getAssignment()));
        }

        // solver.printTerminationTree();

        // std::cout << "Result: " << (retOracle ? "SAT" : "UNSAT") << std::endl;
        // std::cout << randf.toDIMACS() << std::endl;
    }
}

void scratchTest() {
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

    auto ax1 = Assignment({{"x", true}, {"y", false}});
    auto ax0y1z1 = Assignment({{"x", false}, {"y", true}, {"z", true}});
    std::cout << "fa: " << fa.toString() << std::endl;
    std::cout << "fa (x=1): " << fa.assign(ax1).toString() << std::endl;
    std::cout << "fa (x0y1z1): " << fa.assign(ax0y1z1).toString() << std::endl;
    std::cout << "fa (x0y1z1) false: " << fa.assign(ax0y1z1).isFalse() << std::endl;


    Solver solver = Solver();
    auto ret = solver.isSat(fa);
    std::cout << "fa isSat: " << ret << std::endl;
    std::cout << "fa assignment: " << solver.getAssignment().toString() << std::endl;


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
}

void testSimple1() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"~x", "y"}, {"~y", "z"}});
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple2() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"~a", "b"}, {"~b", "~c"}, {"c", "~d"}});
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple3() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"a", "b"},
                   {"b", "c"},
                   {"~a", "~x", "y"},
                   {"~a", "x", "z"},
                   {"~a", "~y", "z"},
                   {"~a", "x", "~z"},
                   {"~a", "~y", "~z"}});
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

//
// Randomized conformance checking, with checking correctness
// of smaller formulas first.
//
// testConformanceRandomCNF(int niters, int nclauses, int nvars, int clause_size)
//
void testConformance() {
    auto start = high_resolution_clock::now();
    testConformanceRandomCNF(50, 1, 2, 2);
    testConformanceRandomCNF(50, 2, 2, 2);
    testConformanceRandomCNF(50, 2, 3, 2);
    testConformanceRandomCNF(50, 4, 4, 2);
    testConformanceRandomCNF(50, 8, 4, 2);
    testConformanceRandomCNF(50, 16, 8, 4);
    testConformanceRandomCNF(100, 50, 10, 4);

    auto stop = high_resolution_clock::now();
    auto durationMS = duration_cast<milliseconds>(stop - start);
    std::cout << "ran conformance checking in " << durationMS.count() << "ms" << std::endl;
}

void testDIMACSParse() {
    std::string dimacsCNFStr;
    std::string fname = "benchmarks/cnf_samples/simple_v3_c2.cnf";
    std::ifstream t(fname);
    std::stringstream buffer;
    buffer << t.rdbuf();
    dimacsCNFStr = buffer.str();
    CNF f = CNF::fromDIMACS(dimacsCNFStr);
    assert(f.toString() == "{ (x1 | ~x3) & (x2 | x3 | ~x1) }");
}

void testCNF(std::string cnfFile, bool expectSat) {
    std::ifstream t(cnfFile);
    std::stringstream buffer;
    buffer << t.rdbuf();
    CNF f = CNF::fromDIMACS(buffer.str());

    Solver solver = Solver();
    auto start = high_resolution_clock::now();
    auto ret = solver.isSat(f);
    auto stop = high_resolution_clock::now();

    assert(ret == expectSat);

    auto durationMS = duration_cast<milliseconds>(stop - start);
    std::cout << "Checked CNF file '" << cnfFile << "' in " << durationMS.count() << "ms"
              << std::endl;

    // If the formula was satisfiable, test that the discovered satisfying
    // assignment is correct.
    if (ret) {
        assert(f.eval(solver.getAssignment()));
    }
}

int main(int argc, char const* argv[]) {
    // Turn off debug level here.
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Debug, el::ConfigurationType::Format, "%datetime %level %msg");
    el::Loggers::reconfigureLogger("default", defaultConf);

    // testDIMACSParse();

    // testSimple1();
    // testSimple2();
    testSimple3();

    return 0;

    // Turn off debug level here.
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", defaultConf);

    testConformance();

    testCNF("benchmarks/cnf_samples/aim-50-1_6-yes1-4.cnf", true);

    //
    // Graph coloring benchmarks.
    //

    // Easier.
    testCNF("benchmarks/flat30-60/flat30-1.cnf", true);
    testCNF("benchmarks/flat30-60/flat30-2.cnf", true);
    testCNF("benchmarks/flat30-60/flat30-3.cnf", true);
    testCNF("benchmarks/flat30-60/flat30-4.cnf", true);
    testCNF("benchmarks/flat30-60/flat30-5.cnf", true);
    testCNF("benchmarks/flat30-60/flat30-6.cnf", true);

    // Harder.
    testCNF("benchmarks/flat50-115/flat50-1.cnf", true);


    // Still too slow with current implementation.

    // testCNF("benchmarks/flat50-115/flat50-2.cnf", true);
    // testCNF("benchmarks/flat50-115/flat50-3.cnf", true);
    // testCNF("benchmarks/cnf_samples/par8-1-c.cnf", true);
    // testCNF("benchmarks/cnf_samples/aim-100-1_6-no-1.cnf", false);
    // testCNF("benchmarks/cnf_samples/dubois20.cnf", false);

    // CNF t1 = CNF({{"x", "y", "z"}, {"z"}});
    // std::cout << "t1: " << t1.toString() << std::endl;
    // auto t1p = t1.unitPropagate(Literal("z"));
    // std::cout << "t1p: " << t1p.toString() << std::endl;

    return 0;
}
