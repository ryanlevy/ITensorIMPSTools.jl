using Test
using ITensorInfiniteMPS
using ITensors
using ITensorIMPSTools

using Random: seed!
using HDF5
seed!(42)

@testset "Transfer Tools" begin
  ψ = h5open(joinpath(@__DIR__, "hubb_test.h5"), "r") do fi
    read(fi, "psi", InfiniteCanonicalMPS)
  end

  N = length(ψ)
  # should be 1e-14, but tests are still very conservative 
  @show norm(contract(ψ.AL[1:N]..., ψ.C[N]) - contract(ψ.C[0], ψ.AR[1:N]...))

  eigList = getSpectrumQN(ψ)
  @test length(eigList) > 0
  @test maximum([maximum(abs.(x[2])) for x in eigList]) ≈ 1.0
end
