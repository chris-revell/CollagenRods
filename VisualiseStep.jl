#
#  VisualiseStep.jl
#  collagen-diffusive-rods
#
#  Created by Christopher Revell on 07/10/2020.
#
#

module VisualiseStep

using Plots

# Plot all rods
@inline @views function visualiseStep(N,r,Ω,L)
    plot()
    for i=1:N
        plot!([r[i,1],r[i,1]+L.*Ω[i,1]./2.0],[r[i,2],r[i,2]+L.*Ω[i,2]./2.0],[r[i,3],r[i,3]+L.*Ω[i,3]./2.0],xlims=[-1,1],ylims=[-1,1],zlims=[-1,1],linecolor=:red,aspect_ratio=:equal,legend=:false)
        #plot!([r[i,1]-L.*Ω[i,1]./2.0,r[i,1]],[r[i,2]-L.*Ω[i,2]./2.0,r[i,2]],[r[i,3]-L.*Ω[i,3]./2.0,r[i,3]],xlims=[-1,1],ylims=[-1,1],zlims=[-1,1],linecolor=:red,aspect_ratio=:equal,legend=:false)
        #plot!([r[i,1],r[i,1]+L.*Ω[i,1]./2.0],[r[i,2],r[i,2]+L.*Ω[i,2]./2.0],[r[i,3],r[i,3]+L.*Ω[i,3]./2.0],xlims=[-1,1],ylims=[-1,1],zlims=[-1,1],linecolor=:green,aspect_ratio=:equal,legend=:false)
    end
    return nothing
end

export visualiseStep

end
