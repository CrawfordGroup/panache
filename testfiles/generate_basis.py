#!/usr/bin/env python3

import argparse
import re


parser = argparse.ArgumentParser(description='Generate a function that will contain basis set information for a libpanache test')
parser.add_argument("-n", required=True, metavar="func_name", help="The name of the function to create")
parser.add_argument("file", help="The file to use")
args = parser.parse_args()


# Read lines into a single list, stripping from
# both ends, and splitting into components
# I'm using regular expressions so that I can merge delimiters
lines = [re.split(" +", line.strip()) for line in open(args.file)]

nshells = int(lines[0][0])
ncenters = int(lines[0][1])

nshellspercenter = [int(i) for i in lines[1]]
shells = []

curline = 2
for l in range(0,nshells):
  # nprim has to be an int because we will use it later
  #  to calculate line numbers, etc
  shell = {}
  shell['nprim'] = int(lines[curline][1])
  shell['am'] = lines[curline][2]
  shell['ispure'] = lines[curline][3]
  shell['coef'] = []
  shell['exp'] = []

  for (c,e) in lines[curline+1:curline+shell['nprim']+1]:
    shell['coef'].append(c)
    shell['exp'].append(e)

  shells.append(shell)
  curline = curline + shell['nprim'] + 1

if nshells != len(shells):
  raise Exception("Error: nshellss is {0} but I have {1}".format(nshells, len(shells)))

print("//                  From file: {0}".format(args.file))
print("//          Number of centers: {0}".format(ncenters))
print("//           Number of shells: {0}".format(nshells))
print("//     # of shells per center: {0}".format(nshellspercenter))
print("")
print("int {0}(int * & nshellspercenter, C_ShellInfo * &shells)".format(args.n))
print("{")
print("int ncenters = {0};\n"
      "int nshells = {1};\n\n"
      "nshellspercenter = new int[ncenters]\n"
      "shells = new C_ShellInfo[nshells];".format(ncenters, nshells))
print("")

for i in range(0, ncenters):
  print("nshellspercenter[{0}] = {1};".format(i,nshellspercenter[i]))
print("\n")

for i in range(0,nshells):
  nprim = shells[i]['nprim']
  print("shells[{0}].nprim = {1};".format(i, nprim))
  print("shells[{0}].am = {1};".format(i, shells[i]['am']))
  print("shells[{0}].ispure = {1};".format(i, shells[i]['ispure']))
  print("shells[{0}].coef = new double[{1}];".format(i, nprim))
  print("shells[{0}].exp = new double[{1}];".format(i, nprim))

  for j in range(0, nprim):
    print("shells[{0}].exp[{1}] = {2};".format(i, j, shells[i]['exp'][j]))
  for j in range(0, nprim):
    print("shells[{0}].coef[{1}] = {2};".format(i, j, shells[i]['coef'][j]))

  print("")

print("return ncenters;\n")
print("}")

