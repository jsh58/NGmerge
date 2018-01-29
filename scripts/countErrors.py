#!/usr/bin/python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# July 2017

# Count errors (substitions, insertions, and Ns)
#   in a SAM, broken down by quality scores.
#   - errors are determined from each alignment's
#     CIGAR and MD string; no reference genome
#     sequence is required
#   - for merged reads, four separate counts tables
#     are produced: one each for unmerged ends,
#     merge matches, merge mismatches, and merge
#     mismatches due to Ns
#     - distinguishing unmerged from merged is done
#       by the length argument (-l)
#     - distinguishing matches/mismatches/Ns is done
#       by the mismatch file (-m), which is usually
#       produced by NGmerge (-j <file> -b) or by
#       findDiffs.py (or findDiffs_SP.py)

import sys
import gzip
import re
import argparse

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

def parseCigar(cigar, diff):
  '''
  Save positions of inserted bases.
  '''
  ops = re.findall(r'(\d+)([IM])', cigar)
  pos = 0
  for op in ops:
    if op[1] == 'I':
      for i in range(pos, pos + int(op[0])):
        diff[i] = 2
    pos += int(op[0])

def getTag(lis, tag):
  '''
  Get optional tag from a SAM record.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == tag:
      return spl[-1]
  sys.stderr.write('Error! Cannot find %s in SAM record\n' % tag)
  sys.exit(-1)

def loadMismatch(filename, d):
  '''
  Load stitching mismatches to given dict.
  '''
  f = openRead(filename)
  for line in f:
    spl = line.rstrip().split('\t')
    if len(spl) < 6:
      sys.stderr.write('Error! Poorly formatted mismatch record:\n' + line)
      sys.exit(-1)
    if spl[0] in d:
      d[spl[0]].append(spl[1:])
    else:
      d[spl[0]] = [spl[1:]]

  if f != sys.stdin:
    f.close()

def findDiffs(diff, md):
  '''
  Find positions of substitutions using MD.
  '''
  loc = 0  # location on read
  parts = re.findall(r'(\d+|\D+|\^\D+)', md)
  for part in parts:
    try:
      # sequence match
      val = int(part)

      # adjust for insertions (which do not consume MD parts)
      while val:
        if diff[loc] != 2:
          val -= 1
        loc += 1

    except ValueError:

      # skip deletion
      if part[0] == '^':
        pass

      else:

        # substitution
        for i in range(len(part)):
          # skip inserted bases
          while diff[loc] == 2:
            loc += 1
          diff[loc + i] = 1
          loc += 1

def countBases(res, res2, res3, res4, diff,
    qual, start, end, pos, stitch, maxQual):
  '''
  Count matches/mismatches/insertions/Ns.
    For stitched reads, separate counts into
      - unstitched (res)
      - stitch matches (res2)
      - stitch mismatches (res3) -- 'pos' dict has
          list of positions of stitch mismatches
      - stitch mismatches due to Ns (res4)
  '''
  if stitch:
    for i in range(len(diff)):
      q = ord(qual[i]) - 33  # assume Sanger scale
      if q > maxQual or q < 0:
        sys.stderr.write('Error! Quality score outside of range [0, %d]\n' % maxQual)
        sys.exit(-1)
      if i >= start and i < end:
        if i in pos:
          if pos[i][0] == 'N' or pos[i][2] == 'N':
            res4[ q ][ diff[i] ] += 1
          else:
            res3[ q ][ diff[i] ] += 1
        else:
          res2[ q ][ diff[i] ] += 1
      else:
        res[ q ][ diff[i] ] += 1

  else:
    for i in range(len(diff)):
      q = ord(qual[i]) - 33  # assume Sanger scale
      if q > maxQual or q < 0:
        sys.stderr.write('Error! Quality score outside of Sanger range\n')
        sys.exit(-1)
      res[ q ][ diff[i] ] += 1

def printOutput(fOut, res):
  '''
  Produce output.
  '''
  fOut.write('\t'.join(['qual', 'match', 'sub', 'ins', 'N',
    'subRate', 'insRate', 'NRate']) + '\n')
  for i in range(len(res)):
    fOut.write('\t'.join(map(str, [i] + res[i])))
    total = float(sum(res[i]))
    for j in range(1, 4):
      if total:
        fOut.write('\t%.9f' % (res[i][j] / total))
      #else:
      #  fOut.write('\tNA')
    fOut.write('\n')

def processSAM(fIn, fOut, len_R1, len_R2, mismatch,
    stitch, maxQual):
  '''
  Process the SAM file. Count errors.
  '''
  res = [[0, 0, 0, 0] for i in range(maxQual + 1)]  # for collecting results
  res2 = [[0, 0, 0, 0] for i in range(maxQual + 1)]  # for results of stitch matches
  res3 = [[0, 0, 0, 0] for i in range(maxQual + 1)]  # for results of stitch mismatches
  res4 = [[0, 0, 0, 0] for i in range(maxQual + 1)]  # for results of mismatches due to Ns
  d = dict()  # for read headers (checking for duplicates)
  count = 0
  for line in fIn:
    if line[0] == '@': continue
    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record:\n' + line)
      sys.exit(-1)
    flag = int(spl[1])
    if flag & 0x904: continue  # skip unmapped, sec/supp
    if (spl[0], flag & 0xC0) in d:
      sys.stderr.write('Warning! Skipping duplicate for read %s, %d\n' \
        % (spl[0], flag & 0xC0))
      continue
    d[(spl[0], flag & 0xC0)] = 1

    # determine positions of differences using CIGAR and MD
    diff = [0] * len(spl[9])  # assume matches at every position (value=0)
    parseCigar(spl[5], diff)  # add insertions (value=2)
    findDiffs(diff, getTag(spl[11:], 'MD'))  # add substitutions (value=1)
    for i in range(len(spl[9])):  # add Ns (value=3)
      if spl[9][i] == 'N':
        diff[i] = 3

    # determine positions to analyze
    start = 0
    end = len(spl[10])
    if len_R1:
      if flag & 0x10:
        start = max(len(spl[10]) - len_R1, 0)
        end = min(len_R2, len(spl[10]))
      else:
        start = max(len(spl[10]) - len_R2, 0)
        end = min(len_R1, len(spl[10]))

    # find positions of stitch mismatches
    pos = dict()
    if spl[0] in mismatch:
      for t in mismatch[spl[0]]:
        if flag & 0x10:
          pos[len(spl[10])-int(t[0])-1] = t[1:] # adjust if rc
        else:
          pos[int(t[0])] = t[1:]

    # count bases by quality score
    countBases(res, res2, res3, res4, diff,
      spl[10], start, end, pos, stitch, maxQual)
    count += 1

  if stitch:
    fOut.write('Unstitched ends:\n')
  printOutput(fOut, res)
  if stitch:
    fOut.write('\nStitch matches:\n')
    printOutput(fOut, res2)
    fOut.write('\nStitch mismatches:\n')
    printOutput(fOut, res3)
    fOut.write('\nStitch mismatches due to Ns:\n')
    printOutput(fOut, res4)
  sys.stderr.write('Reads analyzed: %d\n' % count)

def main():
  '''Main.'''
  # Set command-line arguments
  parser = argparse.ArgumentParser(prog='python ' + sys.argv[0],
    add_help=False)
  parser._action_groups.pop()

  required = parser.add_argument_group('Required arguments')
  required.add_argument('-i', dest='infile', required=True,
    metavar='<file>', help='SAM file (use \'-\' for stdin; ' +
    'must have \'MD\' optional fields)')
  required.add_argument('-o', dest='outfile', required=True,
    metavar='<file>', help='Output file of base counts (matches and errors ' +
    '[substitutions, insertions, and Ns]), broken down by quality score')

  optional = parser.add_argument_group('Arguments for analyzing merged reads')
  optional.add_argument('-l', dest='length', type=str,
    metavar='<int>', help='Lengths of original PE reads (for asymmetric ' +
    'reads, specify comma-separated values: <int>,<int>)')
  optional.add_argument('-m', dest='misfile', metavar='<file>',
    help='File listing merge mismatches (produced by NGmerge ' +
    '[-j <file> -b] or by findDiffs.py')

  other = parser.add_argument_group('Other options')
  other.add_argument('-q', dest='maxQual', type=int, default=40,
    metavar='<int>', help='Maximum quality score (0-based; def. 40)')
  other.add_argument('-h', '--help', dest='help', action='help',
    help='Show help message and exit')

  # parse req. arguments, open files
  args = parser.parse_args()
  fIn = openRead(args.infile)
  fOut = openWrite(args.outfile)

  # parse merging arguments
  stitch = False   # analyzing stitched reads?
  len_R1 = len_R2 = 0   # length of original reads
  if args.length:
    spl = args.length.split(',')
    len_R1 = int(spl[0])
    if len(spl) > 1:
      len_R2 = int(spl[1])
    else:
      len_R2 = len_R1
    stitch = True
  mismatch = dict()
  if args.misfile:
    loadMismatch(args.misfile, mismatch)

  # process SAM file
  processSAM(fIn, fOut, len_R1, len_R2, mismatch,
    stitch, args.maxQual)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
