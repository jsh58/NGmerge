#!/usr/bin/python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Aug. 2017

# Determining mismatches of merged PE reads
#   from a SeqPrep alignment file (-E).
#   - Output is a tab-delimited file listing the
#     following for each mismatch: read header,
#     position, and the bases and quality scores
#     of the R1 and R2 reads at that position.
#     Note that since the quality scores are not
#     in the SeqPrep alignment file, they will
#     be arbitrarily assigned '!'.

import sys
import gzip
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

def printDiff(f, read, seq1, seq2):
  '''
  Print differences between two seqs.
  '''
  count = 0
  # determine length of leading/trailing gaps ('---')
  leadGap = 0
  while seq1[leadGap] == '-' or seq2[leadGap] == '-':
    leadGap += 1
  tailGap = min(len(seq1), len(seq2)) - 1
  while seq1[tailGap] == '-' or seq2[tailGap] == '-':
    tailGap -= 1
  for j in xrange(leadGap, tailGap + 1):
    if (seq1[j] != seq2[j] or seq1[j] == 'N' or seq2[j] == 'N') \
        and seq1[j] != ' ' and seq2[j] != ' ':
      f.write('\t'.join([read, str(j - leadGap), seq1[j], '!', \
        seq2[j], '!']) + '\n')
      count += 1
  return count

def parseFile(fIn, fOut):
  '''
  Produce list of mismatches from SeqPrep alignment file.
  '''
  line = fIn.readline()
  read = ''
  subj = quer = ''
  count = diffs = 0
  while line:
    spl = line.rstrip()
    if spl[0:3] == 'ID:':
      read = spl[3:].lstrip().split(' ')[0]
    elif spl[0:5] == 'SUBJ:':
      subj = spl[6:]
    elif spl[0:6] == 'READ1:':
      subj = spl[7:]
    elif spl[0:5] == 'QUER:':
      quer = spl[6:]
      diffs += printDiff(fOut, read, subj, quer)
      read = subj = quer = ''
      count += 1
    elif spl[0:6] == 'READ2:':
      quer = spl[7:]
      diffs += printDiff(fOut, read, subj, quer)
      read = subj = quer = ''
      count += 1
    line = fIn.readline()

  return count, diffs

def main():
  '''Main.'''
  # Set command-line arguments
  parser = argparse.ArgumentParser(prog='python ' + sys.argv[0])

  required = parser.add_argument_group('Required arguments')
  required.add_argument('-i', dest='alnfile', required=True,
    metavar='<file>', help='Alignment file produced by SeqPrep (-E argument)')
  required.add_argument('-o', dest='outfile', required=True,
    metavar='<file>', help='Output file of sequence mismatches (listing ' +
    'header, position, bases, and quality scores [\'!\'], tab-delimited)')

  # parse arguments
  args = parser.parse_args()

  # process file
  fIn = openRead(args.alnfile)
  fOut = openWrite(args.outfile)
  count, diffs = parseFile(fIn, fOut)
  sys.stderr.write('Alignments analyzed: %d\n' % count)
  sys.stderr.write('Sequence diffs/Ns printed: %d\n' % diffs)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
