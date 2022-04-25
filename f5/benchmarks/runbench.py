import os
import sys
import subprocess

REDUCE = r"redcsl"
ARG1 = '~/signature-based-gb/f5/benchmarks/{0}/{1}{0}.red'
ARG2 = '-w ~/signature-based-gb/f5/correctness/output.rlg'


UNDETECTABLE = ["cyclic6", "cyclic7", "cyclic8", "eco5", "noon4",
                "ojika4", "ku10", "trinks"]

TINY = ["noon5", "noon6", 
        "cyclic9", "cyclic10"]

LARGE = ["eco7", "katsura6", "nbody4sym", "nbody4", 
        "cassou", "hairer2"]

ENORMOUS = ["henrion5", "nbody5"]


class Bench:
    def __init__(self, dir, time, best_time, samples):
        self.dir = dir
        self.time = time
        self.best_time = best_time
        self.samples = samples

    def __str__(self):
        msg = f"{self.time}s from {self.samples} samples. Old best: {self.best_time}s"
        return msg

def compare_message(bench1, bench2):
    return f"f5 {bench1.time}s vs. groebner {bench2.time}s in {bench1.samples} samples"

def extract_time(output):
    tidx = output.index("TIME")
    fromstart = output[tidx+5:]
    fromstart = fromstart[0:fromstart.index(")")]
    return round(float(fromstart) / 1e3, 3)

def runfile(dir, samples, groebner=False):
    prefix = "gb" if groebner else ""
    
    if not groebner:
        print(f"{dir}: ", end='')
        sys.stdout.flush()

    min_time = float("inf")
    with open(f"{dir}/{prefix}res", 'r') as stored:
        best_time = stored.read()
        if not best_time or best_time == '\n':
            best_time = "inf"
        best_time = float(best_time)

    for i in range(samples):
        p = subprocess.run([REDUCE, ARG2, ARG1.format(dir, prefix)],
                            capture_output=True,
                            shell=False)

        output = p.stdout.decode()
    
        time = extract_time(output)
        if time < min_time:
            min_time = time

        if min_time < best_time:
            with open(f"{dir}/{prefix}res", 'w') as stored:
                stored.write(str(min_time))

    bench = Bench(dir, min_time, best_time, samples)

    return bench

def runbench_small(cmp):
    samples = 1
    totaltime = 0
    print(f"Systems: {TINY}\n----------------------")
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

def runbench_large(cmp):
    samples = 1
    totaltime = 0
    print(f"Systems: {LARGE}\n----------------------")
    for root, dirs, files in os.walk("."):
        for dir in dirs:
            if dir not in LARGE:
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

