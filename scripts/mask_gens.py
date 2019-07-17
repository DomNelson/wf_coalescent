import sys

step_size = 1
if len(sys.argv) > 1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("Reads lines from stdin in format 'time\\tnum_lineages and " +\
                "strips those with time less than step_size apart from previous line")
        print("Usage: cat num_lineages_file.txt | python mask_gens.py step_size > masked_num_lineages_file.txt")
        sys.exit()
    else:
        step_size = float(sys.argv[1])

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
        time += step_size
