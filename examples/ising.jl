using ITensors, ITensorMPS
using ITensorInfiniteMPS

include(
  joinpath(
    pkgdir(ITensorInfiniteMPS), "examples", "vumps", "src", "vumps_subspace_expansion.jl"
  ),
)

##############################################################################
# VUMPS parameters
#

maxdim = 10 # Maximum bond dimension
cutoff = 1e-6 # Singular value cutoff when increasing the bond dimension
max_vumps_iters = 20 # Maximum number of iterations of the VUMPS algorithm at a fixed bond dimension
tol = 1e-5 # Precision error tolerance for outer loop of VUMPS or TDVP
outer_iters = 5 # Number of times to increase the bond dimension
time_step = -Inf # -Inf corresponds to VUMPS
solver_tol = (x -> x / 100) # Tolerance for the local solver (eigsolve in VUMPS and exponentiate in TDVP)
multisite_update_alg = "parallel" # Choose between ["sequential", "parallel"]
conserve_qns = true
N = 2 # Number of sites in the unit cell (1 site unit cell is currently broken)

# Parameters of the transverse field Ising model
model_params = (J=1.0, h=0.9)

converge_well = false # do another 1e-14 loop

##############################################################################
# CODE BELOW HERE DOES NOT NEED TO BE MODIFIED
#

model = Model("ising")
initstate(n) = "↑"
s = infsiteinds("S=1/2", N; initstate, conserve_szparity=conserve_qns)
ψ = InfMPS(s, initstate)

# Form the Hamiltonian
H = InfiniteSum{MPO}(model, s; model_params...)

# Check translational invariance
@show norm(contract(ψ.AL[1:N]..., ψ.C[N]) - contract(ψ.C[0], ψ.AR[1:N]...))

vumps_kwargs = (
  tol=tol,
  maxiter=max_vumps_iters,
  solver_tol=solver_tol,
  multisite_update_alg=multisite_update_alg,
  outputlevel=0,
)
subspace_expansion_kwargs = (cutoff=cutoff, maxdim=maxdim)

ψ = tdvp_subspace_expansion(
  H, ψ; time_step, outer_iters, subspace_expansion_kwargs, vumps_kwargs
)

# Check translational invariance
@show norm(contract(ψ.AL[1:N]..., ψ.C[N]) - contract(ψ.C[0], ψ.AR[1:N]...))

stop = 100
d_exact = [expect(ψ, "Sz", i) for i in 1:stop]
using ITensorIMPSTools

d = ITensorIMPSTools.finite_onsite(ψ, "Sz", stop)

SzSz_approx = correlation_approx(ψ, "Sz", "Sz", stop) - d * d'
SzSz_exact = correlation_exact(ψ, "Sz", "Sz", stop) - d_exact * d_exact'

@show SzSz_approx ≈ SzSz_exact

@show sum(SzSz_approx .< 0)
@show sum(SzSz_exact .< 0)

vumps_kwargs = (
  tol=1e-13,
  maxiter=400,
  solver_tol=solver_tol,
  multisite_update_alg=multisite_update_alg,
  outputlevel=0,
)
ψ2 = tdvp_subspace_expansion(
  H, ψ; time_step, outer_iters=1, subspace_expansion_kwargs, vumps_kwargs
)

# Check translational invariance
@show norm(contract(ψ2.AL[1:N]..., ψ2.C[N]) - contract(ψ2.C[0], ψ2.AR[1:N]...))

stop = 100
Op1 = "Sz"
Op2 = "Sz"
d1_exact = [expect(ψ2, Op1, i) for i in 1:stop]
d2_exact = [expect(ψ2, Op2, i) for i in 1:stop]

# the diagonals will be measured differently than the Celled object
d1 = ITensorIMPSTools.finite_onsite(ψ2, Op1, stop; redo_gauge=false)
d2 = ITensorIMPSTools.finite_onsite(ψ2, Op2, stop; redo_gauge=false)

d1_g = ITensorIMPSTools.finite_onsite(ψ2, Op1, stop)
d2_g = ITensorIMPSTools.finite_onsite(ψ2, Op2, stop)

SzSz_approx = correlation_approx(ψ2, Op1, Op2, stop; redo_gauge=false) - d1 * d2'
SzSz_approx_gauged = correlation_approx(ψ2, Op1, Op2, stop; redo_gauge=true) - d1_g * d2_g'
SzSz_exact = correlation_exact(ψ2, Op1, Op2, stop; tol=1e-8) - d1_exact * d2_exact'
@show SzSz_approx ≈ SzSz_exact
@show SzSz_approx_gauged ≈ SzSz_exact

@show sum(SzSz_approx .< 0)
@show sum(SzSz_approx_gauged .< 0)
@show sum(SzSz_exact .< 0)
nothing
