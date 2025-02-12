
T = TransferMatrix(ψ.AL)
vL, vR = ITensorIMPSTools.getLR(T)

spot1 = 1
spot2 = 100
maxUnitCell = ceil(Int, spot2 / N)
L = vL
LN = vL
for i in 1:(maxUnitCell * N)
  global L, LN
  if i == spot1 || i == spot2
    L *= ψ.AL[i] * op("Sz", s[i])
    L *= dag(prime(ψ.AL[i]))
  else
    L *= ψ.AL[i]
    L *= dag(prime(ψ.AL[i], !s[i]))
  end
  LN *= ψ.AL[i]
  LN *= dag(prime(ψ.AL[i], !s[i]))
  @assert order(L) == 2
end
pn = scalar(LN * (translator(ψ))(vR, maxUnitCell - 1))

@show scalar(L * translator(ψ)(vR, maxUnitCell - 1)) / pn
