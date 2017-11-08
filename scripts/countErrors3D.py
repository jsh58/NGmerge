#!/usr/bin/python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Aug. 2017

# Count substitution errors of merged reads in a
#   SAM, broken down by merge matches/mismatches/Ns
#   and both original reads' quality scores.
#   - errors are determined from each alignment's
#     CIGAR and MD string; no reference genome
#     sequence is required
#   - the merged alignments are reconstructed from
#     the original reads; they are assumed not to
#     have indels
#   - four separate counts tables are produced:
#     one for unmerged ends (broken down by quality
#     score), and one each for merge matches, merge
#     mismatches, and merge mismatches due to Ns
#     (which are broken down by both reads' quality
#     scores)

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
    seq1, qual1, seq2, qual2, maxQual):
  '''
  Count matches/mismatches/insertions/Ns.
    For stitched reads, separate counts into
      - unstitched (res)
      - stitch matches (res2)
      - stitch mismatches (res3)
      - stitch mismatches due to Ns (res4)
  '''
  for i in range(len(diff)):
    if i < len(qual1):
      q1 = ord(qual1[i]) - 33  # assume Sanger scale
    if i < len(qual2):
      q2 = ord(qual2[i]) - 33  # assume Sanger scale
    if (q1 > maxQual or q1 < 0) and qual1[i] != ' ':
      sys.stderr.write('Error! Quality score \'%s\' outside of range [0, %d]\n' % (qual1[i], maxQual))
      sys.exit(-1)
    if (q2 > maxQual or q2 < 0) and qual2[i] != ' ':
      sys.stderr.write('Error! Quality score \'%s\' outside of range [0, %d]\n' % (qual2[i], maxQual))
      sys.exit(-1)
    if i >= len(seq1) or seq1[i] == ' ':
      res[ q2 ][ diff[i] ] += 1
    elif i >= len(seq2) or seq2[i] == ' ':
      res[ q1 ][ diff[i] ] += 1
    elif seq1[i] == 'N' or seq2[i] == 'N':
      res4[ q1 ][ q2 ][ diff[i] ] += 1
    elif seq1[i] != seq2[i]:
      res3[ q1 ][ q2 ][ diff[i] ] += 1
    else:
      res2[ q1 ][ q2 ][ diff[i] ] += 1

def printOutput(fOut, res):
  '''
  Produce output.
  '''
  count = 0
  fOut.write('\t'.join(['qual', 'match', 'sub', 'ins', 'N',
    'subRate', 'insRate', 'NRate']) + '\n')
  for i in range(len(res)):
    fOut.write('\t'.join(map(str, [i] + res[i])))
    total = float(sum(res[i]))
    count += sum(res[i])
    for j in range(1, 4):
      if total:
        fOut.write('\t%.9f' % (res[i][j] / total))
      #else:
      #  fOut.write('\tNA')
    fOut.write('\n')
  return count

def print3d(fOut, res):
  '''Print 3d results. Ignore insertions, Ns'''
  count = 0
  fOut.write('\t' + '\t'.join(map(str, range(len(res)))) + '\n')
  for i in range(len(res)):
    fOut.write(str(i))
    for j in range(len(res[i])):
      fOut.write('\t%d/%d' % (res[i][j][1], sum(res[i][j][0:2])))
      count += sum(res[i][j])  # include insertions and Ns in count
    fOut.write('\n')
  return count

def alignReads(head, flag, raw, length):
  '''Align paired reads. Assume no indels.'''
  if head not in raw:
    sys.stderr.write('Error! Cannot find %s in raw reads\n' % head)
    sys.exit(-1)
  if flag & 0x10:
    offset = length - len(raw[head][0])
    if offset < 0:
      seq1 = revComp(raw[head][0])[-offset:]
      qual1 = raw[head][1][offset-1::-1]
    elif offset >= 0:
      seq1 = ' ' * (offset) + revComp(raw[head][0])
      qual1 = ' ' * (offset) + raw[head][1][::-1]
    seq2 = raw[head][2]
    qual2 = raw[head][3]
  else:
    offset = length - len(raw[head][2])
    if offset < 0:
      seq2 = revComp(raw[head][2])[-offset:]
      qual2 = raw[head][3][offset-1::-1]
    elif offset >= 0:
      seq2 = ' ' * (offset) + revComp(raw[head][2])
      qual2 = ' ' * (offset) + raw[head][3][::-1]
    seq1 = raw[head][0]
    qual1 = raw[head][1]
  return seq1, qual1, seq2, qual2

def processSAM(fIn, fOut, raw, maxQual):
  '''
  Process the SAM file. Count errors.
  '''
  res = [[0, 0, 0, 0] for i in range(maxQual + 1)]  # for results of unstitched ends
  res2 = [[[0, 0, 0, 0] for i in range(maxQual + 1)] for j in range(maxQual + 1)]  # for results of stitch matches
  res3 = [[[0, 0, 0, 0] for i in range(maxQual + 1)] for j in range(maxQual + 1)]  # for results of stitch mismatches
  res4 = [[[0, 0, 0, 0] for i in range(maxQual + 1)] for j in range(maxQual + 1)]  # for results of mismatches due to Ns
  d = dict()  # for read headers (checking for duplicates)
  count = 0
  bases = 0
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

    # determine alignment using raw dict and length of seq
    seq1, qual1, seq2, qual2 = alignReads(spl[0], flag, raw, len(spl[9]))

    # determine positions of differences using CIGAR and MD
    diff = [0] * len(spl[9])  # assume matches at every position (value=0)
    parseCigar(spl[5], diff)  # add insertions (value=2)
    findDiffs(diff, getTag(spl[11:], 'MD'))  # add substitutions (value=1)
    for i in range(len(spl[9])):  # add Ns (value=3)
      if spl[9][i] == 'N':
        diff[i] = 3

    # count bases by quality score
    countBases(res, res2, res3, res4, diff,
      seq1, qual1, seq2, qual2, maxQual)
    count += 1

  # print output
  sys.stderr.write('\t' + str(count))
  fOut.write('Unstitched ends:\n')
  tally = printOutput(fOut, res)
  sys.stderr.write('\t' + str(tally))  # unstitched
  fOut.write('\nStitch matches:\n')
  tally = print3d(fOut, res2)
  sys.stderr.write('\t' + str(tally))  # stitch matches
  fOut.write('\nStitch mismatches:\n')
  tally = print3d(fOut, res3)
  sys.stderr.write('\t' + str(tally))  # stitch mismatches
  fOut.write('\nStitch mismatches due to Ns:\n')
  tally = print3d(fOut, res4)
  sys.stderr.write('\t' + str(tally))  # mismatches due to Ns
  sys.stderr.write('\n')

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

def loadReads(r1, r2):
  '''Return a dict of fastq reads loaded from r1/r2.'''
  raw = dict()
  line1 = r1.readline()
  line2 = r2.readline()
  while line1 and line2:
    if line1[0] != '@':
      sys.stderr.write('Error! Not FASTQ format\n')
      sys.exit(-1)
    head = line1.rstrip().split(' ')[0][1:]
    if head != line2.rstrip().split(' ')[0][1:]:
      sys.stderr.write('Error! R1/R2 files do not match\n')
      sys.exit(-1)
    for i in xrange(3):
      line1 = r1.readline()
      line2 = r2.readline()
      if i == 0:
        seq1 = line1.rstrip()
        seq2 = line2.rstrip()
      elif i == 2:
        qual1 = line1.rstrip()
        qual2 = line2.rstrip()
    raw[head] = (seq1, qual1, seq2, qual2)
    line1 = r1.readline()
    line2 = r2.readline()

  if r1 != sys.stdin:
    r1.close()
  if r2 != sys.stdin:
    r2.close()
  return raw

def main():
  '''Main.'''
  # Set command-line arguments
  parser = argparse.ArgumentParser(prog='python ' + sys.argv[0])

  required = parser.add_argument_group('Required arguments')
  required.add_argument('-i', dest='infile', required=True,
    metavar='<file>', help='SAM file (use \'-\' for stdin; ' +
    'must have \'MD\' optional fields)')
  required.add_argument('-1', dest='r1file', required=True,
    metavar='<file>', help='FASTQ file of original R1 reads')
  required.add_argument('-2', dest='r2file', required=True,
    metavar='<file>', help='FASTQ file of original R2 reads')
  required.add_argument('-o', dest='outfile', required=True,
    metavar='<file>', help='Output file of substitution error rates, ' +
    'broken down by overlap (matches/mismatches/Ns) and by quality ' +
    'scores of both original reads')

  # parse arguments
  args = parser.parse_args()

  # open input files
  fIn = openRead(args.infile)
  fR1 = openRead(args.r1file)
  fR2 = openRead(args.r2file)

  # load original fastq reads
  d = loadReads(fR1, fR2)

  # process SAM file
  fOut = openWrite(args.outfile)
  maxQual = 40
  processSAM(fIn, fOut, d, maxQual)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
