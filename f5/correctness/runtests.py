import os
import sys
import subprocess

REDUCE = r"reduce"
ARG1 = '~/signature-based-gb/f5/correctness/{0}.red'

KEYWORDS = ['Wrong Answer', 'Regression', 'Modular Error']

class Report:
    def __init__(self, msg, ok):
        self.msg = msg
        self.ok = ok
    def __str__(self):
        if self.ok:
            return "Tests Passed!"
        else:
            return f"Tests Failed:\n{self.msg}"

def slice_message(msg, i):
    return msg[max(0, i - 20):min(len(msg), i + 40)]

def compare_output(test, true):

    for kw in KEYWORDS:
        if kw in test:
            i = test.index(kw)
            one = slice_message(test, i)
            msg = f"----------- {kw} -----------\n{one}\n----------------------------"
            return Report(msg, False)

    ltest = list(test)
    ltrue = list(true)

    nonsign = ['\r', '\n']

    i, j = 0, 0
    while i < len(ltest) and j < len(ltrue):
        while i < len(ltest) and ltest[i] in nonsign:
            i += 1
        while j < len(ltrue) and ltrue[j] in nonsign:
            j += 1

        if i < len(ltest) and j < len(ltrue) and ltest[i] != ltrue[j]:
            one = slice_message(test, i)
            two = slice_message(true, j)
            msg = f"----------- Expected -----------\n{two}\n----------- Found -----------\n{one}\n----------------------------"
            return Report(msg, False)
        i += 1
        j += 1

    return Report("", True)

def runtests_f5(upd):
    p = subprocess.run([REDUCE, ARG1.format("general_tests")],
                        capture_output=True,
                        shell=True)

    output = p.stdout.decode()

    if upd:
        with open("f5correct.rlg", 'w') as correctfile:
            correctfile.write(output)

    with open("f5correct.rlg", 'r') as correctfile:
        correct = correctfile.read()

    report = compare_output(output, correct)

    return report


def runtests_mod(upd):
    p = subprocess.run([REDUCE, ARG1.format("modular_tests")],
                        capture_output=True,
                        shell=True)

    output = p.stdout.decode()

    if upd:
        with open("modularcorrect.rlg", 'w') as correctfile:
            correctfile.write(output)

    with open("modularcorrect.rlg", 'r') as correctfile:
        correct = correctfile.read()

    report = compare_output(output, correct)

    return report

def runtests_regression(upd):
    p = subprocess.run([REDUCE, ARG2, ARG1.format("regression_tests")],
                        capture_output=True,
                        shell=True)
    print(p)
    output = p.stdout.decode()

    if upd:
        with open("regrcorrect.rlg", 'w') as correctfile:
            correctfile.write(output)

    with open("regrcorrect.rlg", 'r') as correctfile:
        correct = correctfile.read()

    report = compare_output(output, correct)

    return report

def main(argv):

    upd = "-upd" in argv

    if "-f5" in argv:
        print("Running F5 tests..")
        report = runtests_f5(upd)
        print(report)
        print("------------------------")

    if "-mod" in argv:
        print("Running modular tests..")
        report = runtests_mod(upd)
        print(report)
        print("------------------------")

    print("Running regression tests..")
    report = runtests_regression(upd)
    print(report)

if __name__ == '__main__':
    main(sys.argv)