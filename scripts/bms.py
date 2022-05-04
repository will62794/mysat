from pysat.solvers import Minisat22
from pysat.formula import CNF
import time
import subprocess
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys

bms = [
    "benchmarks/cnf_samples/aim-50-1_6-yes1-4.cnf",
    "benchmarks/cnf_samples/aim-100-1_6-no-1.cnf",
    # Random-3-SAT instances with controlled backbone size.
    # "benchmarks/CBS_k3_n100_m403_b10/CBS_k3_n100_m403_b10_0.cnf",
    # "benchmarks/CBS_k3_n100_m403_b10/CBS_k3_n100_m403_b10_1.cnf",
    # "benchmarks/CBS_k3_n100_m403_b10/CBS_k3_n100_m403_b10_2.cnf",
   
    # Graph coloring.
    "benchmarks/flat30-60/flat30-1.cnf",
    "benchmarks/flat30-60/flat30-2.cnf",
    "benchmarks/flat30-60/flat30-3.cnf",
    "benchmarks/flat30-60/flat30-4.cnf",
    "benchmarks/flat50-115/flat50-1.cnf",
    "benchmarks/flat50-115/flat50-2.cnf",
    "benchmarks/flat50-115/flat50-3.cnf",
    "benchmarks/flat50-115/flat50-4.cnf",
    "benchmarks/flat50-115/flat50-5.cnf",
    "benchmarks/flat50-115/flat50-6.cnf",
    "benchmarks/flat50-115/flat50-7.cnf",
    "benchmarks/flat50-115/flat50-8.cnf",

    # "benchmarks/flat75-180/flat75-1.cnf",
    # "benchmarks/flat75-180/flat75-2.cnf",
    # Pigeonhole.
    "benchmarks/pigeon-hole/hole6.cnf",
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

def run_mysat(bms_to_run, cdcl):
    results = []
    # Run my implementation on benchmarks.
    for bm in bms_to_run:
        csv_result_file = "results/result.csv"
        args = " ".join([csv_result_file, cdcl] + [bm])
        cmd = "./main " + args
        # print(cmd)
        maxtime_secs = 25
        try:
            res = subprocess.run(cmd, shell=True, capture_output=True, timeout=maxtime_secs)
        except subprocess.TimeoutExpired as e:
            print(e)
            results.append({"bm": bm, "duration_ms": maxtime_secs * 1000, "is_sat": False, "timeout": True})
            continue
        csvfile = open(csv_result_file)
        reader = csv.DictReader(csvfile, delimiter=',')
        row = list(reader)[0]
        is_sat = "SAT" if row['is_sat'] else "UNSAT"
        print(f"{bm}: {is_sat}, Solved in {row['duration_ms']}ms")
        csvfile.close()
        results.append({"bm": bm, "duration_ms": float(row['duration_ms']), "is_sat": is_sat, "timeout": False})
    return results

def save_results(results, outfilename):
    outfile = open(outfilename, 'w')
    writer = csv.DictWriter(outfile, fieldnames=["bm", "duration_ms", "is_sat", "timeout"])
    writer.writeheader()
    for r in results:
        writer.writerow(r)
    outfile.close()   

if len(sys.argv) == 1:
    # Run MiniSAT on benchmarks.
    print("# Running MiniSAT22")
    minisat_results = run_minisat22(bms)
    save_results(minisat_results, "results/minisat22_results.csv")

    # Run mysat on benchmarks.
    print("# Running mysat without CDCL")
    mysat_dpll_results = run_mysat(bms, "false")
    save_results(mysat_dpll_results, "results/mysat_results_dpll.csv")

    # Run mysat on benchmarks.
    print("# Running mysat with CDCL")
    mysat_cdcl_results = run_mysat(bms, "true")
    save_results(mysat_cdcl_results, "results/mysat_results_cdcl.csv")


fig = plt.figure()
fig, ax = plt.subplots()
# fig.set_figwidth(20)

reader = csv.DictReader(open("results/minisat22_results.csv"))
minisat_results = list(reader)

reader = csv.DictReader(open("results/mysat_results_dpll.csv"))
mysat_dpll_results = list(reader)

reader = csv.DictReader(open("results/mysat_results_cdcl.csv"))
mysat_cdcl_results = list(reader)

print(mysat_cdcl_results)

# ax = fig.add_axes([0,0,1,1])
mysat_dpll_durations = np.array([float(r["duration_ms"]) for r in mysat_dpll_results])
mysat_cdcl_durations = np.array([float(r["duration_ms"]) for r in mysat_cdcl_results])
minisat_durations = np.array([float(r["duration_ms"]) for r in minisat_results])
results = [
    mysat_dpll_durations, 
    mysat_cdcl_durations, 
    minisat_durations
]
print(results)
width = 0.2
x = np.arange(len(bms))  # the label locations

cs = ["red" if r["timeout"]=="True" else "black" for r in mysat_dpll_results]

# TODO: Consider distinguishing bars by texture so that color can be used to indicate timeouts as well.
r1 = ax.barh(x,results[0],0.2, width, color=cs, label="mysat_dpll")
r1 = ax.barh(x + 0.2,results[1],0.2, width, color="blue", label="mysat_cdcl")
r2 = ax.barh(x + 0.4,results[2], width, color="orange", label="minisat2.2")
# plt.plot(bms,results[0],color="blue")
# plt.plot(bms,results[1],color="orange")
# plt.set_xticks(bms)
# plt.show()
ax.set_yticks(x, [bm[bm.index("/")+1:] for bm in bms])
ax.invert_yaxis()
ax.set_xlabel("time to solve (ms)")
ax.set_xscale('log')
ax.legend()
# ax.bar_label(r1)
# ax.bar_label(r2, padding=3)
plt.tight_layout()
plt.savefig("results/compare.pdf")