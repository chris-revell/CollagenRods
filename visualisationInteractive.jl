#
#  Simulate.jl
#  collagen-rods
#
#  Created by Christopher Revell on 23/11/2020.
#
#
# Script to produce 3D visualisation of a single time step that can be rotated and zoomed interactively.

using DelimitedFiles
using Makie
using LinearAlgebra
using AbstractPlotting
using GeometryBasics
using Colors

foldername = ARGS[1]
step = parse(Int64,ARGS[2])

conditions = readdlm(foldername*"/conditions.txt",',')
conditionsDict = Dict(conditions[:,1] .=> conditions[:,2])
nMonomers = conditionsDict["N"]
L = conditionsDict["L"]
σ = conditionsDict["σ"]
containerHeight=conditionsDict["containerHeight"]

readdlm(foldername*"/output.txt",',',Float64)
data = readdlm(foldername*"/output.txt",',',Float64)

r = zeros(Float64,nMonomers,3)
Ω = zeros(Float64,nMonomers,3)

set_theme!(show_axis=false,scale_plot=false,resolution=(1600,1600))
#lims = max(containerHeight,L)
#lim = FRect3D((-lims/2.0,-lims/2.0,-lims/2.0),(lims,lims,lims))
scene = Scene()

r .= data[step*2*nMonomers+1:step*2*nMonomers+nMonomers,:]
Ω .= data[step*2*nMonomers+nMonomers+1:(step+1)*2*nMonomers,:]
for j=1:nMonomers
    r₁ = Point3(r[j,:].-(L/2.0).*Ω[j,:])
    r₂ = Point3(r[j,:].+(L/2.0).*Ω[j,:])
    mesh!(Cylinder(r₁,r₂,σ/2.0),color=:green)
    mesh!(Sphere(r₁,σ),color=:red)
    mesh!(Sphere(r₂,σ),color=:blue)
end
#gui()
display(scene)
readline()
