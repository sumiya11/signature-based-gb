
import subprocess

p = subprocess.run(["reduce", "test.red"], capture_output=True)

print(p)

