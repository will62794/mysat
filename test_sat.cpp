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
        CNF randf = CNF::randomCNF(nclauses, nvars, clause_size);
        std::cout << "- random CNF: " << randf.toString() << std::endl;
        auto retOracle = solver.isSatBruteForce(randf);
        auto retImpl = solver.isSat(randf);
        // std::cout << randf.toDIMACS() << std::endl;
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

void testSimple4() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"x0", "~x0"}, {"~x1"}});
    std::cout << ct1.toDIMACS();
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple5() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"x1"}, {"x2", "~x0"}});
    std::cout << ct1.toDIMACS();
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple6() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"~x1", "~x3"}, {"~x0"}, {"x0", "x2"}, {"~x2"}});
    std::cout << ct1.toDIMACS();
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple7() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"x3", "~x0"}, {"~x2", "~x3"}, {"x2", "~x3"}, {"x3", "~x3"}});
    std::cout << ct1.toDIMACS();
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple8() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"~x2", "~x3"}, {"~x1"}, {"x3", "~x2"}, {"x3", "~x3"}});
    std::cout << ct1.toDIMACS();
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple9() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"~x2", "~x3"}, {"x0", "~x1"}, {"~x1", "~x2"}, {"x2", "~x1"}});
    std::cout << ct1.toDIMACS();
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple10() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"x3", "~x2"}, {"~x2", "~x3"}, {"x0", "x1"}, {"x3", "~x3"}});
    std::cout << ct1.toDIMACS();
    assert(solver.isSatBruteForce(ct1) == solver.isSat(ct1));
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;
    std::cout << solver.getTerminationTreeDOT();
}

void testSimple11() {
    Solver solver = Solver();
    CNF ct1 = CNF({{"x0", "~x1", "~x2"},
                   {"~x0", "~x1", "~x2"},
                   {"x2", "~x2"},
                   {"x1", "~x0", "~x2"},
                   {"x1", "x2", "~x0"},
                   {"x2", "~x1"}});
    std::cout << ct1.toDIMACS();
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

    // Small clause size.
    testConformanceRandomCNF(500, 1, 2, 2);
    testConformanceRandomCNF(500, 2, 2, 2);
    testConformanceRandomCNF(500, 2, 3, 2);
    testConformanceRandomCNF(500, 4, 4, 2);
    testConformanceRandomCNF(500, 6, 4, 2);
    testConformanceRandomCNF(500, 8, 4, 2);

    // Larger clause size.
    testConformanceRandomCNF(1000, 4, 3, 3);
    testConformanceRandomCNF(1000, 6, 3, 3);
    testConformanceRandomCNF(1000, 8, 3, 3);
    testConformanceRandomCNF(1000, 10, 3, 3);
    testConformanceRandomCNF(1000, 12, 3, 3);
    testConformanceRandomCNF(1000, 16, 3, 3);

    // Larger clause size.
    testConformanceRandomCNF(1000, 4, 3, 4);
    testConformanceRandomCNF(1000, 4, 3, 4);
    testConformanceRandomCNF(1000, 6, 3, 4);
    testConformanceRandomCNF(1000, 10, 4, 4);
    testConformanceRandomCNF(1000, 12, 3, 4);
    testConformanceRandomCNF(1000, 32, 8, 4);
    testConformanceRandomCNF(1000, 50, 10, 4);

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
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;


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

    testSimple1();
    testSimple2();
    testSimple3();
    testSimple4();
    testSimple5();
    testSimple6();
    testSimple7();
    testSimple8();
    testSimple9();
    testSimple10();
    testSimple11();

    testCNF("benchmarks/random/random0.cnf", true);
    testCNF("benchmarks/random/random1.cnf", true);
    testCNF("benchmarks/random/random2.cnf", true);
    testCNF("benchmarks/random/random3.cnf", true);
    testCNF("benchmarks/random/random4.cnf", true);
    testCNF("benchmarks/random/random5.cnf", true);

    // return 0;

    // Turn off debug level here.
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", defaultConf);

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
    testCNF("benchmarks/flat30-60/flat30-7.cnf", true);
    testCNF("benchmarks/flat30-60/flat30-8.cnf", true);

    // Harder.
    testCNF("benchmarks/flat50-115/flat50-1.cnf", true);
    testCNF("benchmarks/flat50-115/flat50-2.cnf", true);
    testCNF("benchmarks/flat50-115/flat50-3.cnf", true);
    testCNF("benchmarks/flat50-115/flat50-4.cnf", true);
    testCNF("benchmarks/flat50-115/flat50-5.cnf", true);
    testCNF("benchmarks/flat50-115/flat50-6.cnf", true);
    testCNF("benchmarks/flat50-115/flat50-7.cnf", true);
    testCNF("benchmarks/flat50-115/flat50-8.cnf", true);

    testCNF("benchmarks/cnf_samples/par8-1-c.cnf", true);
    testCNF("benchmarks/cnf_samples/aim-100-1_6-no-1.cnf", false);
    testCNF("benchmarks/cnf_samples/dubois20.cnf", false);

    //
    // Pigeonhole benchmark (UNSAT).
    //
    testCNF("benchmarks/pigeon-hole/hole6.cnf", false);

    // Randomized conformance testing.
    testConformance();

    return 0;
}
