using ITensorInfiniteMPS:
  TransferMatrix,
  input_inds,
  InfiniteCanonicalMPS,
  siteinds,
  finite_mps,
  translator,
  set_ortho_lims!
using ITensors: random_itensor, dag
using ITensorMPS: correlation_matrix
using KrylovKit: eigsolve, schursolve, Arnoldi
using ITensors:
  dag, swapprime, ITensor, commonind, adapt, datatype, using_auto_fermion, scalartype, inds
using ITensorMPS:
  _op_prod,
  siteinds,
  op,
  orthogonalize,
  MPS,
  has_fermion_string,
  norm,
  delta,
  prime,
  siteind,
  scalar
using ITensorInfiniteMPS: TransferMatrix, translator

function getLR(T; tol=1e-14, neigs=1)
  vR = getLM(T; tol, neigs)
  vL = getLM(transpose(T); tol, neigs)
  return vL, vR
  #  #T = TransferMatrix(ψ.AL)
  #  Tᵀ = transpose(T)
  #  vⁱᴿ = random_itensor(dag(input_inds(T)))
  #  vⁱᴸ = random_itensor(dag(input_inds(Tᵀ)))
  #
  #  λR, vR, right_info = eigsolve(T, vⁱᴿ, neigs, :LM; tol=tol)
  #  λL, vL, left_info = eigsolve(Tᵀ, vⁱᴸ, neigs, :LM; tol=tol)
  #  return vL[1], vR[1]
end

function getLM(T; tol=1e-14, neigs=1)
  vi = random_itensor(dag(input_inds(T)))
  TT, vecs, vals, info = schursolve(T, vi, neigs, :LM, Arnoldi(; tol))
  if size(TT, 2) > 1 && TT[2, 1] != 0
    @warn("something wrong with LM vectors")
  end
  return vecs[1]
end

function finishRightMult(psi, N_, L_, currSpot, vR, normN)
  maxUnitCell = ceil(Int, currSpot / N_)
  L = L_
  for i in (currSpot + 1):(maxUnitCell * N_)
    L *= psi[i]
    L *= dag(prime(psi[i], !siteind(psi, i)))
  end
  #vR = tlator(vR, maxUnitCell - 1)
  return scalar(L * vR) / normN
end

function getAllNorms(psi, maxUnitCell, vL, vR)
  norms = Vector{Number}(undef, maxUnitCell)
  L = vL
  T = transpose(TransferMatrix(psi.AL))
  vR = translator(psi)(vR, -1)
  for i in 1:maxUnitCell
    L = T(L)
    norms[i] = scalar(L * vR)
  end
  return norms
end

# correlation_matrix but always dots <L| [op strings] |R>
# as well as continuing the unit cell
# TODO: norms is nearly the same but maybe it isn't for low convergence?
function correlation_matrix_env(psi, _Op1, _Op2, lEnv, rEnv, stop; ishermitian=nothing)
  sites = 1:stop
  site_range = nothing
  if !isnothing(site_range)
    @warn "The `site_range` keyword arg. to `correlation_matrix` is deprecated: use the keyword `sites` instead"
    sites = site_range
  end
  if !(sites isa AbstractRange)
    sites = collect(sites)
  end

  start_site = first(sites)
  end_site = last(sites)

  N = stop
  N_ = length(psi)
  sToC = x -> ceil(Int, x / N_)
  ElT = scalartype(psi)
  s = siteinds(psi)

  tlator = translator(psi)
  norms = getAllNorms(psi, sToC(stop), lEnv, rEnv)

  Op1 = _Op1 #make copies into which we can insert "F" string operators, and then restore.
  Op2 = _Op2
  onsiteOp = _op_prod(Op1, Op2)
  fermionic1 = has_fermion_string(Op1, s[start_site])
  fermionic2 = has_fermion_string(Op2, s[end_site])
  if fermionic1 != fermionic2
    error(
      "correlation_matrix: Mixed fermionic and bosonic operators are not supported yet."
    )
  end

  # Decide if we need to calculate a non-hermitian corr. matrix, which is roughly double the work.
  is_cm_hermitian = ishermitian
  if isnothing(is_cm_hermitian)
    # Assume correlation matrix is non-hermitian
    is_cm_hermitian = false
    O1 = op(Op1, s[start_site])
    O2 = op(Op2, s[start_site])
    O1 /= norm(O1)
    O2 /= norm(O2)
    #We need to decide if O1 ∝ O2 or O1 ∝ O2^dagger allowing for some round off errors.
    eps = 1e-10
    is_op_proportional = norm(O1 - O2) < eps
    is_op_hermitian = norm(O1 - dag(swapprime(O2, 0, 1))) < eps
    if is_op_proportional || is_op_hermitian
      is_cm_hermitian = true
    end
    # finally if they are both fermionic and proportional then the corr matrix will
    # be anti symmetric insterad of Hermitian. Handle things like <C_i*C_j>
    # at this point we know fermionic2=fermionic1, but we put them both in the if
    # to clarify the meaning of what we are doing.
    if is_op_proportional && fermionic1 && fermionic2
      is_cm_hermitian = false
    end
  end

  #psi = orthogonalize(psi, start_site)
  #norm2_psi = norm(psi[start_site])^2

  # Nb = size of block of correlation matrix
  Nb = length(sites)

  C = zeros(ElT, Nb, Nb)

  #if start_site == 1
  #  L = ITensor(1.0)
  #else
  #  lind = commonind(psi[start_site], psi[start_site - 1])
  #  L = delta(dag(lind), lind')
  #end
  L = lEnv
  pL = start_site - 1

  for (ni, i) in enumerate(sites[1:(end - 1)])
    while pL < i - 1
      pL += 1
      sᵢ = siteind(psi, pL)
      L = (L * psi[pL]) * prime(dag(psi[pL]), !sᵢ)
    end

    Li = L * psi[i]

    # Get j == i diagonal correlations
    rind = commonind(psi[i], psi[i + 1])
    oᵢ = adapt(datatype(Li), op(onsiteOp, s[i]))
    #C[ni, ni] = ((Li * oᵢ) * prime(dag(psi[i]), !rind))[] / norm2_psi
    tempL = (Li * oᵢ) * dag(psi[i])'
    C[ni, ni] = finishRightMult(
      psi, N_, tempL, i, tlator(rEnv, sToC(i) - 1), norms[sToC(i)]
    )

    # Get j > i correlations
    if !using_auto_fermion() && fermionic2
      Op1 = "$Op1 * F"
    end

    oᵢ = adapt(datatype(Li), op(Op1, s[i]))

    Li12 = (dag(psi[i])' * oᵢ) * Li
    pL12 = i

    for (n, j) in enumerate(sites[(ni + 1):end])
      nj = ni + n

      while pL12 < j - 1
        pL12 += 1
        if !using_auto_fermion() && fermionic2
          oᵢ = adapt(datatype(psi[pL12]), op("F", s[pL12]))
          Li12 *= (oᵢ * dag(psi[pL12])')
        else
          sᵢ = siteind(psi, pL12)
          Li12 *= prime(dag(psi[pL12]), !sᵢ)
        end
        Li12 *= psi[pL12]
      end

      lind = commonind(psi[j], Li12)
      Li12 *= psi[j]

      oⱼ = adapt(datatype(Li12), op(Op2, s[j]))
      sⱼ = siteind(psi, j)
      #val = (Li12 * oⱼ) * prime(dag(psi[j]), (sⱼ, lind))

      # XXX: This gives a different fermion sign with
      # ITensors.enable_auto_fermion()
      # val = prime(dag(psi[j]), (sⱼ, lind)) * (oⱼ * Li12)
      tempL = (Li12 * oⱼ) * dag(psi[j])'

      C[ni, nj] = finishRightMult(
        psi, N_, tempL, j, tlator(rEnv, sToC(j) - 1), norms[sToC(j)]
      )
      #C[ni, nj] = scalar(val) / norm2_psi
      if is_cm_hermitian
        C[nj, ni] = conj(C[ni, nj])
      end

      pL12 += 1
      if !using_auto_fermion() && fermionic2
        oᵢ = adapt(datatype(psi[pL12]), op("F", s[pL12]))
        Li12 *= (oᵢ * dag(psi[pL12])')
      else
        sᵢ = siteind(psi, pL12)
        Li12 *= prime(dag(psi[pL12]), !sᵢ)
      end
      @assert pL12 == j
    end #for j
    Op1 = _Op1 #"Restore Op1 with no Fs"

    if !is_cm_hermitian #If isHermitian=false the we must calculate the below diag elements explicitly.

      #  Get j < i correlations by swapping the operators
      if !using_auto_fermion() && fermionic1
        Op2 = "$Op2 * F"
      end
      oᵢ = adapt(datatype(psi[i]), op(Op2, s[i]))
      Li21 = (Li * oᵢ) * dag(psi[i])'
      pL21 = i
      if !using_auto_fermion() && fermionic1
        Li21 = -Li21 #Required because we swapped fermionic ops, instead of sweeping right to left.
      end

      for (n, j) in enumerate(sites[(ni + 1):end])
        nj = ni + n

        while pL21 < j - 1
          pL21 += 1
          if !using_auto_fermion() && fermionic1
            oᵢ = adapt(datatype(psi[pL21]), op("F", s[pL21]))
            Li21 *= oᵢ * dag(psi[pL21])'
          else
            sᵢ = siteind(psi, pL21)
            Li21 *= prime(dag(si[pL21]), !sᵢ)
          end
          Li21 *= psi[pL21]
        end

        lind = commonind(psi[j], Li21)
        Li21 *= psi[j]

        oⱼ = adapt(datatype(psi[j]), op(Op1, s[j]))
        sⱼ = siteind(psi, j)
        #val = (prime(dag(psi[j]), (sⱼ, lind)) * (oⱼ * Li21))[]
        #C[nj, ni] = val / norm2_psi
        #tempL = (L * (oⱼ * psi[j])) * dag(psi[j])'
        tempL = dag(psi[j])' * (oⱼ * Li21)
        C[nj, ni] = finishRightMult(
          psi, N_, tempL, j, tlator(rEnv, sToC(j) - 1), norms[sToC(j)]
        )

        pL21 += 1
        if !using_auto_fermion() && fermionic1
          oᵢ = adapt(datatype(psi[pL21]), op("F", s[pL21]))
          Li21 *= (oᵢ * dag(psi[pL21])')
        else
          sᵢ = siteind(psi, pL21)
          Li21 *= prime(dag(psi[pL21]), !sᵢ)
        end
        @assert pL21 == j
      end #for j
      Op2 = _Op2 #"Restore Op2 with no Fs"
    end #if is_cm_hermitian

    # use translational invariance
    (i == N_) && break

    pL += 1
    sᵢ = siteind(psi, i)
    L = Li * prime(dag(psi[i]), !sᵢ)
  end #for i

  ## Get last diagonal element of C
  #i = end_site
  #while pL < i - 1
  #  pL += 1
  #  sᵢ = siteind(psi, pL)
  #  L = L * psi[pL] * prime(dag(psi[pL]), !sᵢ)
  #end
  #lind = commonind(psi[i], psi[i - 1])
  #oᵢ = adapt(datatype(psi[i]), op(onsiteOp, s[i]))
  #sᵢ = siteind(psi, i)
  ##val = (L * (oᵢ * psi[i]) * prime(dag(psi[i]), (sᵢ, lind)))[]
  #tempL = (L * (oᵢ * psi[i])) * dag(psi[i])'
  #C[Nb, Nb] = finishRightMult(
  #  psi, tlator, N_, tempL, i, tlator(rEnv, sToC(i) - 1), norms[sToC(i)]
  #)
  #C[Nb, Nb] = val / norm2_psi

  # everything should get mapped back into the appropriate unit cell and its distance
  # careful to handle non-commuting ops
  for i in (N_ + 1):stop
    for j in (N_ + 1):stop
      a = (i > j) ? (i - 1) % N_ + 1 : (j - 1) % N_ + 1
      b = a + abs(i - j)
      C[i, j] = (i > j) ? C[a, b] : C[b, a]
    end
  end

  return C
end

function correlation_slow(
  ψ::InfiniteCanonicalMPS, op1::String, op2::String, stop::Int; tol=1e-14, kwargs...
)
  T = TransferMatrix(ψ.AL)
  vL, vR = getLR(T; tol)
  return correlation_matrix_env(ψ, op1, op2, vL, vR, stop; kwargs...)
end

function correlation_fast(ψ::InfiniteCanonicalMPS, op1::String, op2::String, stop::Int)
  corr_infinite = correlation_matrix(
    finite_mps(ψ, 1:(stop + 1)), op1, op2; sites=2:(stop + 1)
  )
  return corr_infinite
end

function finite_onsite(ψ::InfiniteCanonicalMPS, op::String, stop::Int)
  return expect(finite_mps(ψ, 1:(stop + 1)), op; sites=2:(stop + 1))
end
