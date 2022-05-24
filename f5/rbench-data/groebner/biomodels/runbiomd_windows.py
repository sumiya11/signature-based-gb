
######## DISCLAIMER #########
# works only on windows #

import os
import sys
import subprocess
import time
import cpuinfo

######## CONFIGURATION ########

NUM_CORES = 1
NUM_SAMPLES = 3
TIMEOUT = 180.0  # 180 seconds

REDUCE = r"C:\Program Files\Reduce\lib\csl\reduce.exe"
BIOMDIR = "C:\\data\\projects\\mpi\\try2\\signature-based-gb\\f5\\rbench-data\\groebner\\biomodels"
ARG2 = '-w C:\\data\\projects\\mpi\\signature-based-gb\\f5\\correctness\\output.rlg'

########## TEMPLATES ##########

FILE_DIR = f"{BIOMDIR}\\"

load_f5 = "in \"C:\\data\\projects\\mpi\\signature-based-gb\\f5\\f5.red\"$"
load_groebner = "load_package groebner$"

torder = "torder({}, revgradlex)$"

run_f5 = "f5(odes)"
run_groebner = "groebner(odes)"

FILE_TEMPLATE = "{0}\n" + "{1}\n\n" + "{3}\n" + "operator diff$\n" +\
                   "odes := {2}$\n" +\
                   "odes := for each o in odes collect part(o, 2)$\n\n" +\
                   "gb := {4}$\n\n" +\
                   "end; % of file"

########## MARKDOWN ##########

MD_HEADLINES = ["f5", "groebner"]

########### TIMINGS ###########

EXCEPT = [
    # "BIOMD0000000002",
    # "BIOMD0000000011",
    # "BIOMD0000000014",
    # "BIOMD0000000038",
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
    tmpname = FILE_DIR + filename + "\\" + filename + ".red"
    with open(tmpname, "w") as tmpfile:
        tmpfile.write(FILE_TEMPLATE.format(load_groebner, torder, odes,
                                        parameters, run_groebner))
    return tmpname

def cleanup():
    _, dirs, _ = next(os.walk("."))
    print(dirs)
    for dir in dirs:
        _, _, files = next(os.walk(f".\\{dir}"))
        for f in files:
            if f == (dir + ".red"):
                os.remove(F"{BIOMDIR}\\{dir}\\{f}")

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

    with open("benchmark_result", "w") as iomd:
        iomd.write(md)

def runall(f5=False, groebner=False):
    # for filename in BENCHMARKS:
    _, dirs, _ = next(os.walk("."))
    for filename in dirs:
        print(filename)
        if filename in EXCEPT:
            print("skipping..")
            continue
        with open(f"{BIOMDIR}\\{filename}\\odes.txt", "r") as fodes:
            odes = fodes.read()
            with open(f"{BIOMDIR}\\{filename}\\parameters.txt", "r") as fparameters:
                parameters = fparameters.read()
                parameters = parameters.replace("=", ":=")
                fname = create_file(filename, odes, parameters)
                # elapsed = runsystem(fname)
                # runtimes[filename] = elapsed
    # cleanup()

def main(argv):
    f5 = "-f5" in argv
    groebner = "-groebner" in argv
    runall(f5=f5, groebner=groebner)
    # markdown()

if __name__ == '__main__':
    main(sys.argv)
