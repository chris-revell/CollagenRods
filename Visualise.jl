#
#  Visualise.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 22/11/2020.
#
#

module Visualise

using DelimitedFiles
using Makie
using LinearAlgebra
using AbstractPlotting
using GeometryBasics
using Colors
using Printf


function visualise(foldername,nTrimers,L,σ)

    nImages = 100
    data = readdlm(foldername*"/output.txt",',',Float64)
    nImages = floor(Int64,size(data)[1]/(2.0*nParticles))
    r = zeros(Float64,nTrimers,3)
    Ω = zeros(Float64,nTrimers,3)
    r₁ = fill(Point3(0.0,0.0,0.0),nTrimers)
    r₂ = fill(Point3(0.0,0.0,0.0),nTrimers)
    objects = fill(Cylinder(Point3(0.0,0.0,0.0),Point3(1.0,1.0,1.0),1.0),nTrimers)

    for i in 0:nImages-1
        println("Rendering $(i+1)/$nImages")
        scene = Scene()
        set_theme!(show_axis = false, scale_plot = false, resolution = (1600, 1600))
        r .= data[i*2*nTrimers+1:i*2*nTrimers+nTrimers,:]
        Ω .= data[i*2*nTrimers+nTrimers+1:(i+1)*2*nTrimers,:]
        for j=1:nTrimers
            r₁[j] = Point3(r[j,:].-(L/2.0).*Ω[j,:])
            r₂[j] = Point3(r[j,:].+(L/2.0).*Ω[j,:])
            objects[j] = Cylinder(r₁[j],r₂[j],σ/2.0)
            mesh!(objects[j],color=:green)
        end
        save("$foldername/$(@sprintf("%03d",i)).png",scene)

        #meshscatter(objects,color=:green)

    end

end

export visualise

end
