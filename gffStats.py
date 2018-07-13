#!/usr/bin/env python3

import sys
import textwrap
import statistics
import re

def unique_genes(gff_infile):
  infile = open(gff_infile, "r")
  pattern = re.compile("gene=[a-zA-Z0-9_]+;")

  total_matches = 0
  unique_genes = []
  for line in infile:
    match = pattern.search(line)
    if match != None:
      total_matches += 1
      if match.group(0) not in unique_genes:
        print(match.group(0)[5:-1])
        unique_genes.append(match.group(0))

  infile.close()
  return True

def unique_coding_genes(gff_infile):
  infile = open(gff_infile, "r")
  pattern = re.compile("gene=[a-zA-Z0-9_]+;")

  unique_genes = []
  for line in infile:
    match = pattern.search(line)
    if match != None and "protein_coding" in line:
      if match.group(0) not in unique_genes:
        print(match.group(0)[5:-1])
        unique_genes.append(match.group(0))

  infile.close()
  return True

def rrna(gff_infile):
  infile = open(gff_infile, "r")
  pattern = re.compile("gene=[a-zA-Z0-9_]+;")

  unique_genes = []
  for line in infile:
    match = pattern.search(line)
    if match != None and "rRNA" in line:
      if match.group(0) not in unique_genes:
        print(match.group(0)[5:-1])
        unique_genes.append(match.group(0))

  infile.close()
  return True

def avgExon(gff_infile):
  infile = open(gff_infile, "r")
  pattern = re.compile("gene=[a-zA-Z0-9_]+;")

  match = None
  genes = {}

  for line in infile:
    if len(line) > 0:
      if line[0] != "#":
        match = pattern.search(line)
        # 3rd field of GFF file contains feature type name, e.g. gene
        # See GFF specification at: https://useast.ensembl.org/info/website/upload/gff.html
        feature = line.split("\t")[2]
        if match != None and feature == "gene":
          gene = match.group(0)[5:-1]
          genes[gene] = []
        if match != None and feature == "exon":
          # 4th and 5th fields of GFF file contains start and end positions of the feature
          # GFF position numbering starts at 1; nonetheless, since we are calculating lengths, it is irrelevant
          start, end = line.split("\t")[3], line.split("\t")[4]
          length = int(end) - int(start) + 1
          genes[gene].append(length)

  infile.close()

  for gene in genes:
    try:
      mean = statistics.mean(genes[gene])
      print(gene + "\t" + str(len(genes[gene])) + "\t" + str(mean))
    except statistics.StatisticsError:
      continue
  return True

def versionInfo():
  print("gffStats 1.01 :: Hyunjin Park 2018")
  print("hyunjin.park@utexas.edu")
  print("my_github")
  return

def fieldInfo():
  print("avgExon: gene\tnumber of exons\taverage exon length")
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
  print("\ngffStats provides utilities to calculate elementary statistics of your FASTA file.")
  print("Usage: python3 gffStats.py [util] [infile]")
  print("type \"python3 gffStats.py -m\" for descriptions of output fields")
  print("type \"python3 gffStats.py -v\" for version information")
  print("type \"python3 gffStats.py -c\" for copyright statement")
  print("type \"python3 gffStats.py -h\" to see this message")
  print()
  print("unique_genes\treturns the list of unique genes")
  print("unique_coding_genes\treturns the list of unique coding genes")
  print("rrna\treturns the list of rRNA genes")
  print("avgExon\treturns the average number of exons and the average exon lengths for each coding gene")
  print()
  print("Check out fastaStats!\n")
  return

def main(argv):
  _utils = ["unique_genes", "unique_coding_genes", "rrna", "avgExon", "-v", "-m", "-c"]
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

  