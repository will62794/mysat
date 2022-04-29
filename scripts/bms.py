from pysat.solvers import Minisat22
from pysat.formula import CNF
import time
import subprocess
import csv
import matplotlib.pyplot as plt
import numpy as np

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
    # "benchmarks/pigeon-hole/hole6.cnf",
    # "benchmarks/pigeon-hole/hole7.cnf",
    # "benchmarks/pigeon-hole/hole8.cnf"
]

def run_minisat22(bms_to_run):
    results = []
    for bm in bms_to_run:
        f = CNF(from_file=bm)
        m = Minisat22(bootstrap_with=f.clauses)
        start = time.time()
        ret = m.solve()
        dur = time.time() - start
        is_sat = "SAT" if ret else "UNSAT"
        print(f"{bm}: {is_sat}, Solved in {round(dur*1000, 2)}ms")
        results.append({"bm": bm, "duration_ms": round(dur*1000, 2), "is_sat":is_sat})
    return results

def run_mysat(bms_to_run):
    results = []
    # Run my implementation on benchmarks.
    for bm in bms_to_run:
        csv_result_file = "results/result.csv"
        args = " ".join([csv_result_file] + [bm])
        cmd = "./main " + args
        # print(cmd)
        res = subprocess.run(cmd, shell=True, capture_output=True)
        csvfile = open(csv_result_file)
        reader = csv.DictReader(csvfile, delimiter=',')
        row = list(reader)[0]
        is_sat = "SAT" if row['is_sat'] else "UNSAT"
        print(f"{bm}: {is_sat}, Solved in {row['duration_ms']}ms")
        csvfile.close()
        results.append({"bm": bm, "duration_ms": float(row['duration_ms']), "is_sat": is_sat})
    return results

def save_results(results, outfilename):
    outfile = open(outfilename, 'w')
    writer = csv.DictWriter(outfile, fieldnames=["bm", "duration_ms", "is_sat"])
    writer.writeheader()
    for r in results:
        writer.writerow(r)
    outfile.close()   

# Run MiniSAT on benchmarks.
print("# Running MiniSAT22")
minisat_results = run_minisat22(bms)
save_results(minisat_results, "results/minisat22_results.csv")

# Run mysat on benchmarks.
print("# Running mysat")
mysat_results = run_mysat(bms)
save_results(mysat_results, "results/mysat_results.csv")



fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
mysat_durations = [r["duration_ms"] for r in mysat_results]
minisat_durations = [r["duration_ms"] for r in minisat_results]
results = [mysat_durations, minisat_durations]
print(results)
ax.plot(bms,results[0],color="blue")
ax.plot(bms,results[1],color="orange")
# plt.show()
plt.savefig("results/compare")