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

data = np.genfromtxt("{}/output.txt".format(argv[1]),delimiter=", ")

conditions = {}
with open("{}/conditions.txt".format(argv[1])) as f:
    for line in f:
        (key, val) = line.split(",")
        conditions[key] = float(val)

boxSize = 2.0*float(conditions["boxSize"])
σ = float(conditions["σ"])
N = int(conditions["N"])
L = float(conditions["L"])

nImages = int(data.shape[0]/(2*N))

for i in range(nImages):
    #os.system("clear")
    print("Rendering {}".format(i))
    outfile = open("{}/povrayTmp{:03d}.pov".format(argv[1],i),"w")
    outfile.write("#include \"colors.inc\"\n")
    outfile.write("camera {\n" )
    outfile.write("  sky <0,0,1>           \n")
    outfile.write("  direction <-1,0,0>      \n")
    outfile.write("  right <-4/3,0,0>      \n")
    outfile.write("  location <{},{},{}> \n".format(boxSize*5,boxSize*2,boxSize*2))
    outfile.write("  look_at <0.,0.,0.>     \n" )
    outfile.write("  angle 15      \n")
    outfile.write("}\n")
    outfile.write("global_settings { ambient_light White }\n")
    outfile.write("light_source {\n" )
    outfile.write("  <{},{},{}>   \n".format(boxSize*15,-boxSize*2,boxSize*2))
    outfile.write("  color White \n")
    outfile.write("}\n")
    outfile.write("background { color White }\n" )

    #outfile.write("box{{<{},{},{}>,<{},{},{}> texture {{pigment{{color White filter 0.8}} finish {{phong 1.0}} }} }}".format(boxSize/2.0,boxSize/2.0,boxSize/2.0,-boxSize/2.0,-boxSize/2.0,-boxSize/2.0))

    r = data[i*N*2:i*N*2+N,:]
    Ω = data[i*N*2+N:(i+1)*N*2,:]

    r1 = r - L*Ω/2.0
    r2 = r + L*Ω/2.0

    for j in range(N):
      outfile.write("cylinder{{<{},{},{}>,<{},{},{}>,{} texture{{pigment{{color Green}}}}}}\n".format(r1[j,0],r1[j,1],r1[j,2],r2[j,0],r2[j,1],r2[j,2],σ/2.0))
      outfile.write("sphere{{<{},{},{}>,{} texture{{pigment{{color Blue}}}}}}\n".format(r1[j,0],r1[j,1],r1[j,2],σ/2.0))
      outfile.write("sphere{{<{},{},{}>,{} texture{{pigment{{color Red}}}}}}\n".format(r2[j,0],r2[j,1],r2[j,2],σ/2.0))

    outfile.close()
    os.system("povray {}/povrayTmp{:03d}.pov +W800 +H800 > /dev/null 2>&1".format(argv[1],i))
    os.system("rm {}/povrayTmp{:03d}.pov".format(argv[1],i))

os.system("convert -delay 10 -loop 0 {}/povrayTmp*.png {}/animated.gif".format(argv[1],argv[1]))
#os.system("rm {}/*.png".format(argv[1]))
