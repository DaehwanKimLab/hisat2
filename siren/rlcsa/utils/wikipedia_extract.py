#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-

import getopt, sys, struct


def main():
  if len(sys.argv) < 4:
    return

  sequences = 0
  in_sequence = False
  part_size = 1048576 * int(sys.argv[2])
  print "Extract size:", part_size

  if sys.argv[1] == "-":
    infile = sys.stdin
  else:
    infile = open(sys.argv[1], "r")
  partname = "extract"

  start_tag = "<" + sys.argv[3] + ">"
  end_tag = "</" + sys.argv[3] + ">"
  print "Tag:", sys.argv[3]

  endmarker = 0
  if len(sys.argv) >= 5:
    endmarker = int(sys.argv[4])
  print "End marker:", endmarker
  marker = struct.pack("B", endmarker)
  print

  output = open(partname, "wb")

  for line in infile:
    if in_sequence:
      if line.find(end_tag) >= 0:
        output.write(marker)
        in_sequence = False
      else:
        output.write(line)
    else:
      if line.find(start_tag) >= 0:
        if output.tell() >= part_size:
          output.close()
          break
        in_sequence = True
        sequences += 1

  infile.close()
  output.close()
  print "Sequences: ", sequences
  print

if __name__ == "__main__":
    main()

