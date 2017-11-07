#!/usr/bin/python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Aug. 2017

# Reconstructing alignments and determining
#   mismatches of merged PE reads.
#   - Alignments are assumed not to have indels.
#   - Output is a tab-delimited file listing the
#     following for each mismatch: read header,
#     position, and the bases and quality scores
#     of the R1 and R2 reads at that position.

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

def revComp(dna):
  '''
  Reverse-complements the given DNA sequence.
  '''
  rc = ''
  for nuc in dna[::-1]:
    comp = ''
    if nuc == 'A': comp = 'T'
    elif nuc == 'C': comp = 'G'
    elif nuc == 'G': comp = 'C'
    elif nuc == 'T': comp = 'A'
    elif nuc == 'N': comp = 'N'
    else:
      sys.stderr.write('Error! Unknown nucleotide: %s\n' % nuc)
      sys.exit(-1)
    rc += comp
  return rc

def loadFastq(f):
  '''Return a dict of fastq read lengths.'''
  d = dict()
  count = 0
  line = f.readline()
  while line:
    head = line.rstrip().split(' ')[0]
    if head[0] != '@':
      sys.stderr.write('Error! Not FASTQ format\n')
      sys.exit(-1)
    for i in xrange(3):
      line = f.readline()
      if not line:
        sys.stderr.write('Error! Not FASTQ format\n')
        sys.exit(-1)
      if i == 0:
        seq = line.rstrip()
    d[head] = len(seq)

    count += 1
    line = f.readline()

  if f != sys.stdin:
    f.close()
  return d, count

def retrieveReads(r1, r2, d):
  '''Retrieve reads from r1/r2 based on headers in d.'''
  raw = dict()
  count = len(d)
  line1 = r1.readline().rstrip().split(' ')[0]
  line2 = r2.readline().rstrip().split(' ')[0]
  while line1 and line2 and count:
    match = False
    head = line1
    if head[0] != '@':
      sys.stderr.write('Error! Not FASTQ format\n')
      sys.exit(-1)
    if line1 != line2:
      sys.stderr.write('Error! R1/R2 files do not match:\n' + line1 + line2)
      sys.exit(-1)
    if head in d:
      match = True
    for i in xrange(3):
      line1 = r1.readline()
      line2 = r2.readline()
      if not line1 or not line2:
        sys.stderr.write('Error! Not FASTQ format\n')
        sys.exit(-1)
      if match:
        if i == 0:
          seq1 = line1.rstrip()
          seq2 = line2.rstrip()
        elif i == 2:
          qual1 = line1.rstrip()
          qual2 = line2.rstrip()
    if match:
      raw[head] = [seq1, qual1, revComp(seq2), qual2[::-1]]
      count -= 1
    line1 = r1.readline().rstrip().split(' ')[0]
    line2 = r2.readline().rstrip().split(' ')[0]

  if r1 != sys.stdin:
    r1.close()
  if r1 != sys.stdin:
    r2.close()
  return raw

def printDiffs(fOut, d, raw):
  '''Print alignment diffs and Ns.'''
  count = 0
  for head in raw:

    # position alignment according to length in d dict
    offset = d[head] - len(raw[head][2])
    if offset < 0:
      raw[head][2] = raw[head][2][-offset:]
      raw[head][3] = raw[head][3][-offset:]
    if offset > 0:
      raw[head][2] = ' ' * (offset) + raw[head][2]
      raw[head][3] = ' ' * (offset) + raw[head][3]

    # find and print diffs, Ns
    for i in range(min(len(raw[head][0]), len(raw[head][2]))):
      if (raw[head][0][i] != raw[head][2][i] and raw[head][2][i] != ' ') \
          or (raw[head][0][i] == 'N' and raw[head][2][i] != ' ') \
          or (raw[head][2][i] == 'N' and raw[head][0][i] != ' '):
        fOut.write('\t'.join([head[1:], str(i), \
          raw[head][0][i], raw[head][1][i], \
          raw[head][2][i], raw[head][3][i]]) + '\n')
        count += 1

  return count

def main():
  '''Main.'''
  # Set command-line arguments
  parser = argparse.ArgumentParser(prog='python ' + sys.argv[0])

  required = parser.add_argument_group('Required arguments')
  required.add_argument('-m', dest='mergefile', required=True,
    metavar='<file>', help='FASTQ file of merged reads')
  required.add_argument('-1', dest='r1file', required=True,
    metavar='<file>', help='FASTQ file of original R1 reads')
  required.add_argument('-2', dest='r2file', required=True,
    metavar='<file>', help='FASTQ file of original R2 reads')
  required.add_argument('-o', dest='outfile', required=True,
    metavar='<file>', help='Output file of sequence mismatches (listing ' +
    'header, position, bases, and quality scores, tab-delimited)')

  # parse arguments
  args = parser.parse_args()

  # load merged read lengths
  f = openRead(args.mergefile)
  d, count = loadFastq(f)
  sys.stderr.write('Reads in %s: %d\n' % (args.mergefile, count))

  # retrieve original reads
  r1 = openRead(args.r1file)
  r2 = openRead(args.r2file)
  raw = retrieveReads(r1, r2, d)
  sys.stderr.write('Reads retrieved: %d\n' % (len(raw)))

  # print differences
  fOut = openWrite(args.outfile)
  count = printDiffs(fOut, d, raw)
  sys.stderr.write('Sequence diffs/Ns printed: %d\n' % (count))

  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
