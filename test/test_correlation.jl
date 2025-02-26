using Test
using ITensorInfiniteMPS
using ITensors
using ITensorIMPSTools

using Random: seed!
using HDF5
seed!(42)

#base_path = joinpath(pkgdir(ITensorInfiniteMPS), "examples", "vumps", "src")
#src_files = ["vumps_subspace_expansion.jl", "entropy.jl"]
#for f in src_files
#  include(joinpath(base_path, f))
#end
#
###############################################################################
## VUMPS parameters
##
#
#maxdim = 50 # Maximum bond dimension
#cutoff = 1e-6 # Singular value cutoff when increasing the bond dimension
#max_vumps_iters = 200 # Maximum number of iterations of the VUMPS algorithm at each bond dimension
#vumps_tol = 1e-14
#outer_iters = 5 # Number of times to increase the bond dimension
#localham_type = MPO # or ITensor
#conserve_qns = true
#eager = true
#
#model_params = (t=1.0, U=10.0, V=0.0)
#
###############################################################################
## CODE BELOW HERE DOES NOT NEED TO BE MODIFIED
##
#
#N = 2 # Unit cell size
#
#@show N
#@show localham_type
#
#initstate(n) = isodd(n) ? "↑" : "↓"
#s = infsiteinds("Electron", N; initstate, conserve_qns)
#ψ = InfMPS(s, initstate)
#
#model = Model"hubbard"()
#@show model, model_params
#
## Form the Hamiltonian
#H = InfiniteSum{localham_type}(model, s; model_params...)
#
## Check translational invariance
#println("\nCheck translational invariance of initial infinite MPS")
#@show norm(contract(ψ.AL[1:N]..., ψ.C[N]) - contract(ψ.C[0], ψ.AR[1:N]...))
#
#outputlevel = 1
#vumps_kwargs = (tol=vumps_tol, maxiter=max_vumps_iters, outputlevel, eager)
#subspace_expansion_kwargs = (cutoff=cutoff, maxdim=maxdim)
#
## For now, to increase the bond dimension you must alternate
## between steps of VUMPS and subspace expansion (which outputs
## a new state that is equal to the original state but with
## a larger bond dimension)
#
#println("\nRun VUMPS on initial product state, unit cell size $N")
#ψ = vumps_subspace_expansion(H, ψ; outer_iters, subspace_expansion_kwargs, vumps_kwargs)
#
## Check translational invariance
#println("\nCheck translational invariance of optimized infinite MPS")

@testset "Correlation Checks" begin
  ψ = h5open(joinpath(@__DIR__, "hubb_test.h5"), "r") do fi
    read(fi, "psi", InfiniteCanonicalMPS)
  end

  N = length(ψ)
  # should be 1e-14, but tests are still very conservative 
  @show norm(contract(ψ.AL[1:N]..., ψ.C[N]) - contract(ψ.C[0], ψ.AR[1:N]...))

  @testset "Commuting Obs" begin
    @test round.(correlation_approx(ψ, "Sz", "Sz", 10); digits=8) ≈
      round.(correlation_exact(ψ, "Sz", "Sz", 10); digits=8)

    @test round.(correlation_approx(ψ, "S+", "S-", 10); digits=8) ≈
      round.(correlation_exact(ψ, "S+", "S-", 10); digits=8)
  end
  @testset "Non-commuting Obs" begin
    @test round.(correlation_approx(ψ, "ntot * Sz", "Nup", 10); digits=8) ≈
      round.(correlation_exact(ψ, "ntot * Sz", "Nup", 10); digits=8)
  end
end
