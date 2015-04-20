#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-

import getopt, sys


def main():
  if len(sys.argv) < 4:
    return

  sequences = 0
  in_sequence = False
  part_size = 1048576 * int(sys.argv[2])
  current_file = 1
  print "Part size:", part_size
  print

  if sys.argv[1] == "-":
    infile = sys.stdin
  else:
    infile = open(sys.argv[1], "r")
  partname = "part"

  start_tag = "<" + sys.argv[3] + ">"
  end_tag = "</" + sys.argv[3] + ">"

  output = open(partname + "." + str(current_file), "wb")
  print "Writing part", output.name, "..."

  for line in infile:
    if in_sequence:
      if line.find(end_tag) >= 0:
        output.write("\0")
        in_sequence = False
      else:
        output.write(line)
    else:
      if line.find(start_tag) >= 0:
        if output.tell() >= part_size:
          output.close()
          current_file += 1
          output = open(partname + "." + str(current_file), "wb")
          print "Writing part", output.name, "..."
        in_sequence = True
        sequences += 1

  infile.close()
  output.close()
  print
  print "Sequences: ", sequences

if __name__ == "__main__":
    main()
