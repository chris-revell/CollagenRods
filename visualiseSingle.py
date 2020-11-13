#! /usr/local/bin/python3
#
#  visualise.py
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

import numpy as np
import os
from sys import argv

conditions = {}
with open("{}/conditions.txt".format(argv[1])) as f:
    for line in f:
        (key, val) = line.split(",")
        conditions[key] = float(val)

boxSize = float(conditions["boxSize"])
σ = float(conditions["σ"])
N = int(conditions["N"])
L = float(conditions["L"])

data = np.genfromtxt("{}/output.txt".format(argv[1]),delimiter=", ")

outfile = open("{}/povrayTmp{}.pov".format(argv[1],"RealTime"),"w")
outfile.write("#include \"colors.inc\"\n")
outfile.write("camera {\n" )
outfile.write("  sky <0,0,1>           \n")
outfile.write("  direction <-1,0,0>      \n")
outfile.write("  right <-4/3,0,0>      \n")
outfile.write("  location <{},{},{}> \n".format(boxSize*10,boxSize*2,boxSize*2))
outfile.write("  look_at <0.,0.,0.>     \n" )
outfile.write("  angle 15      \n")
outfile.write("}\n")
outfile.write("global_settings { ambient_light White }\n")
outfile.write("light_source {\n" )
outfile.write("  <{},{},{}>   \n".format(boxSize*15,-boxSize*2,boxSize*2))
outfile.write("  color White \n")
outfile.write("}\n")
outfile.write("background { color White }\n" )

r = data[-N*2:-N,:]
Ω = data[-N:,:]

r1 = r - L*Ω/2.0
r2 = r + L*Ω/2.0

for j in range(N):
  outfile.write("cylinder{{<{},{},{}>,<{},{},{}>,{} texture{{pigment{{color Green}}}}}}\n".format(r1[j,0],r1[j,1],r1[j,2],r2[j,0],r2[j,1],r2[j,2],σ/2.0))
  outfile.write("sphere{{<{},{},{}>,{} texture{{pigment{{color Blue}}}}}}\n".format(r1[j,0],r1[j,1],r1[j,2],σ))
  outfile.write("sphere{{<{},{},{}>,{} texture{{pigment{{color Red}}}}}}\n".format(r2[j,0],r2[j,1],r2[j,2],σ))

outfile.close()
os.system("povray {}/povrayTmp{}.pov +W1600 +H1600 > /dev/null 2>&1".format(argv[1],"RealTime"))
os.system("rm {}/povrayTmp{}.pov".format(argv[1],"RealTime"))
