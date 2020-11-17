#
#  CreateRunDirectory.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 15/10/2020.
#
#

module CreateRunDirectory

using Dates
using Base.Filesystem

@views function createRunDirectory(N,L,σ,ϵ,p,η,kT,tMax,boxSize,D₀,DParallel,DPerpendicular,DRotation,interactionThresh)

    # Create directory for run data labelled with current time.
    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkpath("output/$(foldername)")

    # Store system parameters.
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile,"N,                $N")
        println(conditionsfile,"L,                $L")
        println(conditionsfile,"σ,                $σ")
        println(conditionsfile,"ϵ,                $ϵ")
        println(conditionsfile,"p,                $p")
        println(conditionsfile,"η,                $η")
        println(conditionsfile,"kT,               $kT")
        println(conditionsfile,"tMax,             $tMax")
        println(conditionsfile,"boxSize,          $boxSize")
        println(conditionsfile,"D₀,               $D₀")
        println(conditionsfile,"DParallel,        $DParallel")
        println(conditionsfile,"DPerpendicular,   $DPerpendicular")
        println(conditionsfile,"DRotation,        $DRotation")
        println(conditionsfile,"interactionThresh,$interactionThresh")
    end

    return foldername
end

export createRunDirectory

end
