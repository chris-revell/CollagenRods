using DelimitedFiles


a = readdlm("output.txt",',',Float64)

r = a[4041:60,:]
Ω = a[4061:end,:]

r1 = r .- Ω/2.0
r2 = r .+ Ω/2.0

f = open("test.pov","w")

write(f,"#include \"colors.inc\"\n")
write(f,"camera {\n" )
write(f,"  sky <0,0,1>           \n")
write(f,"  direction <-1,0,0>      \n")
write(f,"  right <-4/3,0,0>      \n")
write(f,"  location <$(boxSize*6),$(boxSize*2),$(boxSize*2)> \n")
write(f,"  look_at <0.,0.,0.>     \n" )
write(f,"  angle 15      \n")
write(f,"}\n")
write(f,"global_settings { ambient_light White }\n")
write(f,"light_source {\n" )
write(f,"  <$(boxSize*5),$(boxSize*2),$(boxSize*2)>   \n")
write(f,"  color White \n")
write(f,"}\n")
write(f,"background { color White }\n" )

for i=1:20
   write(f,"cylinder{<$(r1[i,1]),$(r1[i,2]),$(r1[i,3])>,<$(r2[i,1]),$(r2[i,2]),$(r2[i,3])>,$(σ/2.0) texture{pigment{color Green}}}\n")
end
