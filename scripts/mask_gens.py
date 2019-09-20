import sys
import argparse

parser = argparse.ArgumentParser(
        description="Reads lines from stdin in format 'time\\tnum_lineages and " +\
                "strips those with time less than step_size apart from previous" +\
                "line. If provided, recent_step_size applies until switch_time" +\
                "generations in the past."
        )
parser.add_argument("--global_step_size", type=float, default=1)
parser.add_argument("--recent_step_size", type=float, default=0.001)
parser.add_argument("--switch_time", type=float, default=20)
args = parser.parse_args()

time = 0
while True:
    line = sys.stdin.readline().strip()
    if line == '':
        break

    try:
        t, num_lineages = line.strip().split('\t')
    except ValueError:
        print("Error parsing line:", file=sys.stderr)
        print(list(line), file=sys.stderr)
        continue

    if float(num_lineages) > 1e6 or float(num_lineages) != int(float(num_lineages)):
        print("Error parsing line:", file=sys.stderr)
        print(list(line), file=sys.stderr)
        continue

    if float(t) >= time:
        print(t + '\t' + num_lineages)
        while float(t) >= time:
            if time < args.switch_time:
                time += args.recent_step_size
            else:
                time += args.global_step_size
