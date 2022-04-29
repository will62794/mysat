from pysat.solvers import Minisat22
from pysat.formula import CNF
import time
import subprocess
import csv

bms = [
    "benchmarks/cnf_samples/aim-50-1_6-yes1-4.cnf",
    "benchmarks/cnf_samples/aim-100-1_6-no-1.cnf",
    # Random-3-SAT instances with controlled backbone size.
    # "benchmarks/CBS_k3_n100_m403_b10/CBS_k3_n100_m403_b10_0.cnf",
    # "benchmarks/CBS_k3_n100_m403_b10/CBS_k3_n100_m403_b10_1.cnf",
    # "benchmarks/CBS_k3_n100_m403_b10/CBS_k3_n100_m403_b10_2.cnf",
    # Graph coloring.
    "benchmarks/flat50-115/flat50-1.cnf",
    "benchmarks/flat50-115/flat50-2.cnf",
    "benchmarks/flat50-115/flat50-3.cnf",
    "benchmarks/flat50-115/flat50-4.cnf",
    "benchmarks/flat50-115/flat50-5.cnf",
    "benchmarks/flat50-115/flat50-6.cnf",
    "benchmarks/flat50-115/flat50-7.cnf",
    "benchmarks/flat50-115/flat50-8.cnf",
    "benchmarks/flat75-180/flat75-1.cnf",
    "benchmarks/flat75-180/flat75-2.cnf",
    "benchmarks/flat75-180/flat75-3.cnf",
    "benchmarks/flat100-239/flat100-1.cnf",
    # Pigeonhole.
    "benchmarks/pigeon-hole/hole6.cnf",
    # "benchmarks/pigeon-hole/hole7.cnf",
    # "benchmarks/pigeon-hole/hole8.cnf"
]

def run_minisat22(bms_to_run):
    for bm in bms_to_run:
        f = CNF(from_file=bm)
        m = Minisat22(bootstrap_with=f.clauses)
        start = time.time()
        ret = m.solve()
        dur = time.time() - start
        is_sat = "SAT" if ret else "UNSAT"
        print(f"{bm}: {is_sat}, Solved in {round(dur*1000, 2)}ms")

def run_mysat(bms_to_run):
    # Run my implementation on benchmarks.
    for bm in bms_to_run:
        csv_result_file = "results/mysat_results.csv"
        args = " ".join([csv_result_file] + [bm])
        cmd = "./main " + args
        # print(cmd)
        res = subprocess.run(cmd, shell=True, capture_output=True)
        csvfile = open(csv_result_file)
        reader = csv.DictReader(csvfile, delimiter=',')
        row = list(reader)[0]
        is_sat = "SAT" if row['is_sat'] else "UNSAT"
        print(f"{bm}: {is_sat}, Solved in {row['duration_ms']}ms")

# Run MiniSAT on benchmarks.
print("# Running MiniSAT22")
run_minisat22(bms)

# Run mysat on benchmarks.
print("# Running mysat")
run_mysat(bms)