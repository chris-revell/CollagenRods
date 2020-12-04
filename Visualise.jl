#
#  Visualise.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 22/11/2020.
#
#
# Functions to visualise results using Makie.

module Visualise

using DelimitedFiles
using Makie
using LinearAlgebra
using AbstractPlotting
using GeometryBasics
using Colors
using Printf


function visualise(foldername,nTrimers,L,σ,boxSize)

    data = readdlm(foldername*"/output.txt",',',Float64)
    nImages = floor(Int64,size(data)[1]/(2.0*nTrimers))
    r = zeros(Float64,nTrimers,3)
    Ω = zeros(Float64,nTrimers,3)

    set_theme!(show_axis=false,scale_plot=false,resolution=(1600,1600))
    lims = max(boxSize,L)
    lim = FRect3D((-lims/2.0,-lims/2.0,-lims/2.0),(lims,lims,lims))
    scene = Scene(limits=lim)
    mesh!(Sphere(Point3(zeros(3)),1.0))
    save("$foldername/Tmp001.png",scene)

    for i in 0:nImages-1
        run(`clear`)
        println("Rendering $(i+1)/$nImages")
        scene = Scene(limits=lim)

        r .= data[i*2*nTrimers+1:i*2*nTrimers+nTrimers,:]
        Ω .= data[i*2*nTrimers+nTrimers+1:(i+1)*2*nTrimers,:]
        for j=1:nTrimers
            r₁ = Point3(r[j,:].-(L/2.0).*Ω[j,:])
            r₂ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(L/3.0).*Ω[j,:])
            mesh!(Cylinder(r₁,r₂,σ/2.0),color=:red)

            r₁ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(L/3.0).*Ω[j,:])
            r₂ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(2.0*L/3.0).*Ω[j,:])
            mesh!(Cylinder(r₁,r₂,σ/2.0),color=:green)

            r₁ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(2.0*L/3.0).*Ω[j,:])
            r₂ = Point3(r[j,:].+(L/2.0).*Ω[j,:])
            mesh!(Cylinder(r₁,r₂,σ/2.0),color=:blue)
        end
        save("$foldername/Tmp$(@sprintf("%03d",i+1)).png",scene)
    end
    run(`clear`)
    println("Animating...")
    run(`convert -delay 10 -loop 0 $foldername/Tmp"*".png $foldername/animated.gif`)
end

function visualise(foldername)

    conditions = readdlm(foldername*"/conditions.txt",',')
    conditionsDict = Dict(conditions[:,1] .=> conditions[:,2])
    nTrimers = conditionsDict["N"]
    L = conditionsDict["L"]
    σ = conditionsDict["σ"]
    boxSize=conditionsDict["boxSize"]

    readdlm(foldername*"/output.txt",',',Float64)
    data = readdlm(foldername*"/output.txt",',',Float64)
    nImages = floor(Int64,size(data)[1]/(2.0*nTrimers))

    r = zeros(Float64,nTrimers,3)
    Ω = zeros(Float64,nTrimers,3)

    set_theme!(show_axis=false,scale_plot=false,resolution=(1600,1600))
    lims = max(boxSize,L)
    lim = FRect3D((-lims/2.0,-lims/2.0,-lims/2.0),(lims,lims,lims))
    scene = Scene(limits=lim)
    mesh!(Sphere(Point3(zeros(3)),1.0))
    save("$foldername/Tmp001.png",scene)

    for i in 0:nImages-1
        run(`clear`)
        println("Rendering $(i+1)/$nImages")
        scene = Scene(limits=lim)

        r .= data[i*2*nTrimers+1:i*2*nTrimers+nTrimers,:]
        Ω .= data[i*2*nTrimers+nTrimers+1:(i+1)*2*nTrimers,:]
        for j=1:nTrimers
            r₁ = Point3(r[j,:].-(L/2.0).*Ω[j,:])
            r₂ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(L/3.0).*Ω[j,:])
            mesh!(Cylinder(r₁,r₂,σ/2.0),color=:red)

            r₁ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(L/3.0).*Ω[j,:])
            r₂ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(2.0*L/3.0).*Ω[j,:])
            mesh!(Cylinder(r₁,r₂,σ/2.0),color=:green)

            r₁ = Point3(r[j,:].-(L/2.0).*Ω[j,:].+(2.0*L/3.0).*Ω[j,:])
            r₂ = Point3(r[j,:].+(L/2.0).*Ω[j,:])
            mesh!(Cylinder(r₁,r₂,σ/2.0),color=:blue)
        end
        save("$foldername/Tmp$(@sprintf("%03d",i+1)).png",scene)
    end
    run(`clear`)
    println("Animating...")
    run(`convert -delay 10 -loop 0 $foldername/Tmp"*".png $foldername/animated.gif`)
end

export visualise

end
