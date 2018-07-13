#!/usr/bin/env python3

import sys
import textwrap
import statistics
from Bio.Blast import NCBIWWW

def seq_list(fasta_infile):
  infile = open(fasta_infile, "r")
  seqs = ""

  for line in infile:
    if len(line) == 0:
      continue
    elif line[0] == ">":
      seqs += line[1:]
    else:
      continue

  print(seqs)
  return True

def gc(fasta_infile):
  infile = open(fasta_infile, "r")
  scaffolds = {}

  for line in infile:
    line = line.strip()
    if len(line) == 0:
      continue
    elif line[0] == ">":
      scaffolds[line[1:]] = [0, 0]
      key = line[1:]
    else:
      for char in line:
        if char == "G" or char == "C" or char == "g" or char == "c":
          scaffolds[key][0] += 1
        else:
          scaffolds[key][1] += 1
  infile.close()

  for key in scaffolds:
    scaffold_length = sum(scaffolds[key])
    gc_prop = scaffolds[key][0] / scaffold_length
    print(key + "\t" + str(scaffold_length) + "\t" + str(gc_prop))

  return True

def gc_nr(fasta_infile):
  infile = open(fasta_infile, "r")
  scaffolds = {}

  for line in infile:
    line = line.strip()
    if len(line) == 0:
      continue
    elif line[0] == ">":
      scaffolds[line[1:]] = [0, 0, 0]
      key = line[1:]
    else:
      for char in line:
        if char == "A" or char == "T":
          scaffolds[key][0] += 1
        elif char == "G" or char == "C":
          scaffolds[key][1] += 1
        else:
          scaffolds[key][2] += 1
  infile.close()

  for key in scaffolds:
    scaffold_length = sum(scaffolds[key])
    gc_prop = scaffolds[key][1] / (scaffolds[key][0] + scaffolds[key][1])
    print(key + "\t" + str(scaffold_length) + "\t" + str(gc_prop))
  return True

def perc_repeat(fasta_infile):
  infile = open(fasta_infile, "r")
  scaffolds = {}

  for line in infile:
    line = line.strip()
    if len(line) == 0:
      continue
    elif line[0] == ">":
      scaffolds[line[1:]] = [0, 0]
      key = line[1:]
    else:
      for char in line:
        # repeats are marked as lowercases
        if char == "a" or char == "t" or char == "g" or char == "c":
          scaffolds[key][0] += 1
        else:
          scaffolds[key][1] += 1
  infile.close()

  for key in scaffolds:
    scaffold_length = sum(scaffolds[key])
    repeat_prop = scaffolds[key][0] / scaffold_length
    print(key + "\t" + str(scaffold_length) + "\t" + str(repeat_prop))

  return True

def len_sd(fasta_infile):
  infile = open(fasta_infile, "r")
  scaffolds = {}
  tot_length = 0

  for line in infile:
    line = line.strip()
    if len(line) == 0:
      continue
    elif line[0] == ">":
      scaffolds[line[1:]] = 0
      key = line[1:]
    else:
      scaffolds[key] += len(line)
      tot_length += len(line)
  infile.close()

  numseq = len(scaffolds)
  mean = tot_length / numseq
  median = statistics.median(scaffolds.values())
  sample_sd = statistics.stdev(scaffolds.values(), mean)

  print("#MEAN=" + str(mean) + "\tMEDIAN=" + str(median) + "\tNUMSEQ=" + str(numseq))
  print("#You may want to check for outlier(s) if mean and median values greatly differ.")
  for key in scaffolds:
    sd = (scaffolds[key] - mean) / sample_sd
    print(key + "\t" + str(scaffolds[key]) + "\t" + str(sd))
  return True

def fasta_separate(fasta_infile):
  infile = open(fasta_infile, "r")
  seq_processed = 0
  out = ""

  for line in infile:
    if len(line) == 0:
      continue
    elif line[0] == ">":
      # Case: the first sequence
      if out == "":
        out += line
        seq_name = line[1:].strip()
      else:
        outfile = open((seq_name + ".fasta"), "w")
        outfile.write(out)
        outfile.close()
        seq_processed += 1
        if seq_processed % 1000 == 0:
          print(str(seq_processed) + " sequences processed")
        out = line
        seq_name = line[1:].strip()
    else:
      out += line

  # Writes out the very last sequence
  outfile = open((seq_name + ".fasta"), "w")
  outfile.write(out)
  outfile.close()
  return True

def blast(fasta_infile):
  try:
    to_blast = open(fasta_infile).read()
    result = NCBIWWW.qblast("blastn", "nt", to_blast)
    print(result.read())
    return True
  except ModuleNotFoundError:
    print("Please install biopython")
    print("Recommended installation command: python3 -m pip install biopython")
    print("Download BioPython: https://biopython.org/wiki/Download")

def versionInfo():
  print("fastaStats 1.01 :: Hyunjin Park 2018")
  print("hyunjin.park@utexas.edu")
  print("my_github")
  return

def fieldInfo():
  print("seq_list: sequence name")
  print("gc: sequence name\tsequence length\tgc proportion")
  print("gc_nr: sequence name\tsequence length\tgc proportion in region excluding repeats and egenerate bases")
  print("perc_repeat: sequence name\tsequence length\tproportion of repeated region")
  print("len_sd: sequence name\tsequence length\tstandard deviation")
  print("fasta_separate: N./A.")
  return

def copyRight():
  statement = """
    Copyright (c) 2018 Hyunjin Park

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    """
  print(statement)
  return

def printHelp():
  print("\nfastaStats provides utilities to calculate elementary statistics of your FASTA file.")
  print("Usage: python3 fastaStats.py [util] [infile]")
  print("type \"python3 fastaStats.py -m\" for descriptions of output fields")
  print("type \"python3 fastaStats.py -v\" for version information")
  print("type \"python3 fastaStats.py -c\" for copyright statement")
  print("type \"python3 fastaStats.py -h\" to see this message")
  print()
  print("Util\tDescription")
  print("seq_list\tprints the list of sequences in the file")
  print("gc\tcalculates gc propertion")
  print("gc_nr\tcalculates gc proportion, excluding repeats and degenerate bases")
  print("perc_repeat\tcalculates the proportion of repeated region")
  print("len_sd\tcalculates the mean of lengths of given sequences and the standard deviation for each sequence")
  print("fasta_separate\tseparates each sequence record into individual file")
  print("blast\tperforms NCBI blast")
  print()
  print("Check out gffStats!\n")
  return

def main(argv):
  _utils = ["seq_list", "gc", "gc_nr", "perc_repeat", "len_sd", "fasta_separate", "blast", "-v", "-m", "-c"]
  if len(argv) != 2 or argv[0] not in _utils:
    printHelp()
    exit(1)
  elif argv[0] == "-v":
    versionInfo()
    exit()
  elif argv[0] == "-m":
    fieldInfo()
    exit()
  elif argv[0] == "-c":
    copyRight()
    exit()
  else:
    if (not globals()[argv[0]](argv[1])):
      print("An error occurred. The program is terminated.")
      exit(1)

if __name__ == "__main__":
  # Cheers
  main(sys.argv[1:])

