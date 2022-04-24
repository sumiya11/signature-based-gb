import os
import sys
import subprocess

REDUCE = r"C:\Program Files\Reduce\lib\csl\reduce.exe"
ARG2 = '-w C:\\data\\projects\\mpi\\signature-based-gb\\f5\\correctness\\output.rlg'
ARG1 = 'C:\\data\\projects\\mpi\\signature-based-gb\\f5\\benchmarks\{0}\{1}{0}.red'

TINY = ["cyclic6", "cyclic7", "cyclic8",
        "noon4", "noon5", "ojika4", "ku10", "eco5"]
LARGE = []

class Bench:
    def __init__(self, dir, time, best_time, samples):
        self.dir = dir
        self.time = time
        self.best_time = best_time
        self.samples = samples

    def __str__(self):
        msg = f"{self.dir}: {self.time}s from {self.samples} samples. Old best: {self.best_time}s"
        return msg

def compare_message(bench1, bench2):
    return f"{bench1.dir}: f5 {bench1.time}s vs. groebner {bench2.time}s in {bench1.samples} samples"

def extract_time(output):
    tidx = output.index("TIME")
    fromstart = output[tidx+5:]
    fromstart = fromstart[0:fromstart.index(")")]
    return round(float(fromstart) / 1e3, 2)

def runfile(dir, samples, groebner=False):
    prefix = "gb" if groebner else ""
    arg1 = f'C:\\data\\projects\\mpi\\signature-based-gb\\f5\\benchmarks\{dir}\{prefix}{dir}.red'

    min_time = float("inf")
    with open(f"{dir}\\{prefix}res", 'r') as stored:
        best_time = stored.read()
        if not best_time:
            best_time = "inf"
        best_time = float(best_time)

    for i in range(samples):
        p = subprocess.run([REDUCE, ARG2, ARG1.format(dir, prefix)],
                            capture_output=True,
                            shell=True)

        output = p.stdout.decode()

        time = extract_time(output)
        if time < min_time:
            min_time = time

        if min_time < best_time:
            with open(f"{dir}\\{prefix}res", 'w') as stored:
                stored.write(str(min_time))

    bench = Bench(dir, min_time, best_time, samples)

    return bench

def runbench_small(cmp):
    samples = 3
    totaltime = 0
    for root, dirs, files in os.walk("."):
        for dir in dirs:
            if dir not in TINY:
                continue
            if not cmp:
                bench = runfile(dir, samples)
                print(bench)
                totaltime += bench.time
            else:
                bench1 = runfile(dir, samples, groebner=False)
                totaltime += bench1.time
                bench2 = runfile(dir, 1, groebner=True)
                msg = compare_message(bench1, bench2)
                print(msg)

    return totaltime

def main(argv):
    cmp = "-cmp" in argv

    if "-tiny" in argv:
        print("Running tiny benchmarks..\n-------------------")
        report = runbench_small(cmp)
        print(f"-------------------\nTotal time: {report}s")

    if "-large" in argv:
        print("Running LARGE benchmarks..\n-------------------")
        report = runbench_large(cmp)
        print(f"-------------------\nTotal time: {report}s")

if __name__ == '__main__':
    main(sys.argv)
