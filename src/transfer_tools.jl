using ITensorInfiniteMPS: TransferMatrix, combiner, combinedind, dag, input_inds
using ITensors: space, Index
using KrylovKit: eigsolve

function getSpectrumQN(ψ; neigs=4, tol=1e-10)
  eigList = []

  T = TransferMatrix(ψ.AL)
  ii = combinedind(combiner(dag(input_inds(T))))
  qs = space(ii)
  for q in qs
    vⁱᴿ = random_itensor(dag(input_inds(T)), Index(q[1] => 1))
    neigs_ = min(q[2], neigs)
    λR, vvR, right_info = eigsolve(T, vⁱᴿ, neigs_, :LM; tol=tol)
    push!(eigList, (q, λR[1:neigs_]))
  end
  return eigList
end

""" Take list of (QN,size), eigs to a text file """
function write_spec_to_file(eigList, fname)
  neigs = maximum([length(p[2]) for p in eigList])
  open(fname, "w") do fo
    for ((q, qN), eigs) in eigList
      line = "$q $qN " #"$(q[1]) $(q[2])"
      for i in 1:neigs
        if i > length(eigs)
          line *= "0.0 0.0 "
        else
          line *= "$(real(eigs[i])) $(imag(eigs[i])) "
        end
      end
      write(fo, line * "\n")
    end
  end
end
