#include "sat.cpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
using namespace std::chrono;

struct Result {
    int durationMS;
    bool isSat;
};

Result runCNF(std::string cnfFile) {
    std::ifstream t(cnfFile);
    std::stringstream buffer;
    buffer << t.rdbuf();
    CNF f = CNF::fromDIMACS(buffer.str());

    Solver solver = Solver();
    auto start = high_resolution_clock::now();
    auto ret = solver.isSat(f);
    auto stop = high_resolution_clock::now();

    auto durationMS = duration_cast<milliseconds>(stop - start);
    std::cout << "Checked CNF file '" << cnfFile << "' in " << durationMS.count() << "ms"
              << std::endl;
    std::cout << "num conflicts: " << solver.getNumConflicts() << std::endl;


    // If the formula was satisfiable, test that the discovered satisfying
    // assignment is correct.
    if (ret) {
        assert(f.eval(solver.getAssignment()));
    }

    Result res;
    res.durationMS = durationMS.count();
    res.isSat = ret;
    return res;
}

int main(int argc, char const* argv[]) {
    // Turn off debug level here.
    el::Configurations defaultConf;
    defaultConf.setToDefault();
    defaultConf.set(el::Level::Debug, el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureLogger("default", defaultConf);

    if (argc < 3) {
        std::cout << "Must supply CSV out file and CNF files as arguments.";
        return 0;
    }

    std::string outcsvFile(argv[1]);

    std::stringstream outCSV;
    outCSV << "duration_ms,is_sat\n";
    for (int i = 2; i < argc; i++) {
        std::string argCNF(argv[i]);
        std::cout << "Running : " << argCNF << std::endl;
        Result res = runCNF(argCNF);
        outCSV << res.durationMS << "," << res.isSat << "\n";
    }

    // Write the results into CSV.
    std::ofstream csvfile;
    csvfile.open(outcsvFile);
    csvfile << outCSV.str();
    csvfile.close();

    return 0;
}
