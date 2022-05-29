import os
import sys
import subprocess
import time
import cpuinfo

######## CONFIGURATION ########

NUM_CORES = 1
NUM_SAMPLES = 3
TIMEOUT = 3.0  # 3 seconds

REDUCE = r"bootstrapreduce"
BIOMDIR = "~/signature-based-gb/f5/benchmarks/BIOMD"
ARG2 = '-w ~/signature-based-gb/f5/correctness/output.rlg'

########## TEMPLATES ##########

FILE_DIR = f"{BIOMDIR}/"

load_f5 = "load_package f5$"
load_groebner = "load_package groebner$ torder({}, revgradlex)$"

run_f5 = "f5(odes)"
run_groebner = "groebner(odes)"

FILE_TEMPLATE = """{0}\n
                   {1}\n
                   parameters := {3}$\n
                   operator diff$\n
                   odes := {2}$
                   odes := for each o in odes collect part(o, 2)$ \n
                   {4};\n
                   end;"""

########## MARKDOWN ##########

MD_HEADLINES = ["f5", "groebner"]

########### TIMINGS ###########

EXCEPT = [
    "BIOMD0000000002",
    "BIOMD0000000011",
    "BIOMD0000000014",
    "BIOMD0000000038",
]

runtimes = dict()

############ MAIN #############


def runsystem(fname):
    ans = 60*10  # 10 minutes

    for _ in range(NUM_SAMPLES):
        starttime = time.time()
        try:
            p = subprocess.run([REDUCE, ARG2, fname],
                                # shell=True,
                                capture_output=True,
                                timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            print(f"Not finished in {TIMEOUT} seconds")
            return -1

        ans = min(ans, time.time() - starttime)

    return ans

def create_file(filename, odes, parameters):
    tmpname = filename + ".red"
    with open(tmpname, "w") as tmpfile:
        tmpfile.write(FILE_TEMPLATE.format(load_f5, load_groebner, odes,
                                        parameters, run_f5))
    return tmpname

def cleanup():
    _, dirs, files = next(os.walk("."))
    for filename in files:
        if filename.endswith(".red"):
            os.remove(f"{filename}")

def markdown():
    md = ""
    md += "|Model|" + "|".join(MD_HEADLINES) + "|\n"
    md += "|-----|" +  "|".join(["---" for _ in MD_HEADLINES]) + "|\n"
    for name in runtimes.keys():
        md += f"|{name}|"
        md += str(runtimes[name]) + "|" + " - " + "|"
        md += "\n"

    md += "\n*Benchmarking environment:*\n\n"
    cpu = cpuinfo.get_cpu_info()
    cpu = [f"{cpu['arch']}",
            f"{cpu['brand_raw']}"]
    md += '\n\n'.join(cpu) + "\n"

    with open("benchmark_result.md", "w") as iomd:
        iomd.write(md)

def runall(f5=False, groebner=False):
    # for filename in BENCHMARKS:
    _, dirs, _ = next(os.walk("."))
    for filename in dirs:
        print(filename)
        if filename in EXCEPT:
            print("skipping..")
            continue
        with open(f"{filename}/odes.txt", "r") as fodes:
            odes = fodes.read()
            with open(f"{filename}/parameters.txt", "r") as fparameters:
                parameters = fparameters.read()
                fname = create_file(filename, odes, parameters)
                elapsed = runsystem(fname)
                runtimes[filename] = elapsed
    cleanup()

def main(argv):
    f5 = "-f5" in argv
    groebner = "-groebner" in argv
    runall(f5=f5, groebner=groebner)
    markdown()

if __name__ == '__main__':
    main(sys.argv)

