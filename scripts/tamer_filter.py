#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys

"""
To execute the command in Python:
>>> from tamer_filter import filter_table
>>> filter_table('myfile.csv', 0.8, 0.1, upper_count=2, lines_to_skip=1, delimiter=',', outfile='myfile.filtered.csv')

To execute from the shell:
$ ./tamer_filter.py --upper_threshold=0.8 --lower_threshold=0.1 --upper_count=2, --header_lines=1 --delimiter=, myfile.csv > myfile.filtered.csv
"""

def filter_table(filename, upper_threshold, lower_threshold, upper_count,
                 lines_to_skip, delimiter, outfile=sys.stdout):
    with open(filename, 'r') as infile:
        if lines_to_skip > 0:
            for _ in range(lines_to_skip):
                header = next(infile)
                print(header.rstrip())
        for line in infile:
            values = line.rstrip().split(delimiter)
            while True:
                try:
                    values.remove('NA')
                except:
                    break 

            #label = values[0]
            values.pop(0)
            nvalues = [float(v) for v in values]
            nupper = sum([1 for v in nvalues if v >= upper_threshold])
            nlower = sum([1 for v in nvalues if v <= lower_threshold])
            #if nupper >= upper_count and nlower + nupper + 1 == len(values):
            if nupper >= upper_count and nlower + nupper == len(values) and nlower > 1:
                print(line.rstrip(), file=outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--upper_threshold', metavar='UPR', type=float, default=0.8, help='minimum upper threshold (default: 0.8)')
    parser.add_argument('--lower_threshold', metavar='LWR', type=float, default=0.1, help='maximum lower threshold (default: 0.1)')
    parser.add_argument('--upper_count', metavar='CNT', type=int, default=2, help='at least CNT samples must have value UPR or more; remaining samples must have value LWR or less')
    parser.add_argument('--header_lines', metavar='HDR', type=int, default=1, help='number of header lines that should not be parsed (default: 1)')
    parser.add_argument('--delimiter', metavar='DLM', default='\t', help='character separating values (default: tab)')
    parser.add_argument('infile', help='input file')

    args = parser.parse_args()
    filter_table(args.infile, args.upper_threshold, args.lower_threshold,
                 args.upper_count, args.header_lines, args.delimiter)
