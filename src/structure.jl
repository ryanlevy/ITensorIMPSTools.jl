using ITensors: dag, δ, denseblocks, dag, prime, commoninds, noprime, order
using ITensorMPS: linkinds, siteinds
using ITensorInfiniteMPS: InfiniteCanonicalMPS, translatecell, translator, expect, nsites
using KrylovKit: linsolve

# One should really use this, which makes an ITensorMap
# But the main code uses direct methods instead 
# (probably to include the C*C term)

#TransferMatrixMixed(ψ::InfiniteMPS, ψh::InfiniteMPS) = TransferMatrixMixed(ψ, ψh, Cell(1))
#
#function TransferMatrixMixed(ψ::InfiniteMPS, ψh::InfiniteMPS, c::Cell)
#  N = nsites(ψ)
#  ψᴴ = prime(linkinds, dag(ψh))
#  ψᶜ = ψ[c]
#  ψᶜᴴ = ψᴴ[c]
#  r = unioninds(linkinds(ψ, N => N + 1), linkinds(ψᴴ, N => N + 1))
#  l = unioninds(linkinds(ψ, 1 => 0), linkinds(ψᴴ, 1 => 0))
#  return ITensorMap(ψᶜ, ψᶜᴴ; input_inds=r, output_inds=l)
#end

# This code uses notation from https://scipost.org/SciPostPhysLectNotes.7/pdf
# Tangent-space methods for uniform matrix product states

struct Aᴿₗ
  ψ::InfiniteCanonicalMPS
  n::Int
end

function (A::Aᴿₗ)(x)
  ψ = A.ψ
  ψᴴ = dag(ψ)
  ψ′ = ψᴴ'
  ψ̃ = prime(linkinds, ψᴴ)
  n = A.n

  N = length(ψ)
  #@assert n == N

  l = linkinds(only, ψ.AL)
  l′ = linkinds(only, ψ′.AL)
  r = linkinds(only, ψ.AR)
  r′ = linkinds(only, ψ′.AR)

  xT = (translatecell(translator(ψ), x, -1))
  for k in (n - N + 1):n
    xT = xT * ψ.AR[k] * ψ̃.AL[k]
  end
  δˡ = δ(Bool, r[n], l′[n])
  δʳ = δ(Bool, l[n], r′[n])
  xR = x * ψᴴ.C[n] * ψ′.C[n] * dag(δʳ) * denseblocks(δˡ)
  return xT - xR
end

struct Aᴸᵣ
  ψ::InfiniteCanonicalMPS
  n::Int
end

function (A::Aᴸᵣ)(x)
  ψ = A.ψ
  ψᴴ = dag(ψ)
  ψ′ = ψᴴ'
  ψ̃ = prime(linkinds, ψᴴ)
  n = A.n

  N = length(ψ)
  #@assert n == N

  l = linkinds(only, ψ.AL)
  l′ = linkinds(only, ψ′.AL)
  r = linkinds(only, ψ.AR)
  r′ = linkinds(only, ψ′.AR)

  xT = (translatecell(translator(ψ), x, 1))
  for j in reverse(1:N)
    xT = xT * ψ.AL[j] * ψ̃.AR[j]
  end
  δˡ(n) = δ(Bool, r[n], l′[n])
  δʳ(n) = δ(Bool, dag(l[n]), prime(r[n]))
  xR = x * ψᴴ.C[0] * (ψ′.C[0] * δˡ(0)) * denseblocks(δʳ(0))
  return xT - xR
end

function overlapTerm(
  ψ::InfiniteCanonicalMPS, op1::String, s1::Int, op2::String, s2::Int, const_terms::Vector
)
  s = siteinds(ψ)
  N = length(ψ.AL)

  lp = (commoninds(ψ.AL[0], ψ.C[0]))
  ϕ = 1.0
  for j in 1:N
    ϕ *= prime(ψ.AL[j], lp)
    ϕ = if (j == s1)
      noprime(ϕ * (op(s[s1], op1) - const_terms[1] * op(s[s1], "I")), "Site")
    else
      ϕ
    end
    ϕ = if (j == s2)
      noprime(ϕ * (op(s[s2], op2) - const_terms[2] * op(s[s2], "I")), "Site")
    else
      ϕ
    end
    ϕ *= prime(dag(ψ.AL[j]), "Link")
    @assert order(ϕ) == 2
  end
  lpp = commonind(ψ.C[N], prime(ψ.C[N], lp))
  ϕ *= prime(dag(ψ.C[N]), "Link")
  ϕ = noprime(ϕ) * ψ.C[N]
  ans = ϕ[]

  return ans
end

function leftTerm(
  q::Number,
  ψ::InfiniteCanonicalMPS,
  op1::String,
  s1::Int,
  op2::String,
  s2::Int,
  const_terms::Vector;
  exact_fixedpoint=false,
  tol=1e-10,
  debug=false,
  lin_args=(;),
)
  s = siteinds(ψ)
  N = nsites(ψ)

  ψᴴ = prime(dag(ψ), "Link")

  # Form the vector ĀL*Oα*Ac (connected on left)
  lp = (commoninds(ψ.AL[0], ψ.C[0]))
  ϕ = 1.0
  for j in 1:N
    ϕ *= prime(ψ.AL[j], lp)
    ϕ = if (j == s1)
      noprime(ϕ * (op(s[s1], op1) - const_terms[1] * op(s[s1], "I")), "Site")
    else
      ϕ
    end
    ϕ *= ψᴴ.AL[j]
    @assert order(ϕ) == 2
  end
  ϕ *= ψ.C[N]

  A = Aᴿₗ(ψ, N)
  factor = exp(-1im * q)
  factor = imag(factor) == 0 ? real(factor) : factor

  Lα, info = linsolve(A, ϕ, 1, -factor; tol=tol, lin_args...)
  (debug) && @show info.converged == 1
  if info.converged != 1
    @warn "Left term did not converge!"
  end

  # Form the vector Āc*Oβ*AR (connected on right)
  lp = (commoninds(ψ.AR[N], ψ.C[N]))
  ϕ = prime(dag(ψ).C[N], "Link")
  for j in reverse(1:N)
    ϕ *= ψᴴ.AL[j]
    ϕ = if (j == s2)
      noprime(prime(ϕ, "Site") * (op(s[s2], op2) - const_terms[2] * op(s[s2], "I")), "Site")
    else
      ϕ
    end
    ϕ *= prime(ψ.AR[j], lp)
    @assert order(ϕ) == 2
  end

  ϕ = (translatecell(translator(ψ), ϕ, 1))
  total = factor * (ϕ * Lα)[]
  (debug) && @show total
  return total
end

function rightTerm(
  q::Number,
  ψ::InfiniteCanonicalMPS,
  op1::String,
  s1::Int,
  op2::String,
  s2::Int,
  const_terms::Vector;
  exact_fixedpoint=false,
  tol=1e-10,
  debug=false,
  lin_args=(;),
)
  s = siteinds(ψ)
  N = nsites(ψ)

  ψᴴ = prime(dag(ψ), "Link")

  # Form the vector ĀR*Oα*Ac (connected on right)
  lp = linkinds(ψ.AR, (N, N + 1))
  ϕ = prime(ψ.C[N], lp)
  for j in reverse(1:N)
    ϕ *= prime(ψ.AL[j], lp)
    ϕ = if (j == s1)
      noprime(ϕ * (op(s[s1], op1) - const_terms[1] * op(s[s1], "I")), "Site")
    else
      ϕ
    end
    ϕ *= ψᴴ.AR[j]
    @assert order(ϕ) == 2
  end

  A = Aᴸᵣ(ψ, N)
  factor = exp(1im * q)
  factor = imag(factor) == 0 ? real(factor) : factor

  Lα, info = linsolve(A, ϕ, 1, -factor; tol=tol, lin_args...)
  (debug) && @show info.converged == 1
  if info.converged != 1
    @warn "Right term did not converge!"
  end

  # Form the vector Āc*Oβ*AL (connected on left)
  lp = (commoninds(ψ.AL[0], ψ.C[0]))
  ϕ = 1.0
  for j in 1:N
    ϕ *= prime(ψ.AL[j], lp)
    ϕ = if (j == s2)
      noprime(ϕ * (op(s[s2], op2) - const_terms[2] * op(s[s2], "I")), "Site")
    else
      ϕ
    end
    ϕ *= ψᴴ.AL[j]
    @assert order(ϕ) == 2
  end
  ϕ *= prime(dag(ψ.C[N]), "Link")

  ϕ = (translatecell(translator(ψ), ϕ, -1))

  total = factor * (ϕ * Lα)[]
  (debug) && @show total
  return total
end
function static_sq(
  q::Number,
  ψ::InfiniteCanonicalMPS,
  op1::String,
  s1::Int,
  op2::String,
  s2::Int;
  tol=1e-10, # also a lin_arg
  debug=false,
  lin_args=(;),
)
  @assert (s1 <= length(ψ)) && (s2 ≤ length(ψ)) # eventually fix
  const_terms = [expect(ψ, op1, s1), expect(ψ, op2, s2)]
  debug && @show const_terms
  sq = overlapTerm(ψ, op1, s1, op2, s2, const_terms)
  sq += leftTerm(q, ψ, op1, s1, op2, s2, const_terms; tol, debug, lin_args)
  sq += rightTerm(q, ψ, op1, s1, op2, s2, const_terms; tol, debug, lin_args)
  return sq
end
