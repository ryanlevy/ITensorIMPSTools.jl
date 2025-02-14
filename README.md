# ITensorIMPSTools.jl
Tools for ITensorInfiniteMPS

## Features
* Static Structure Factor code
* Exact (and fast) correlator computation
* Transfer Matrix eigenvalues (sorted by QN/charge sectors)

## Installation
Because this package isn't registered please use
```julia
julia> using Pkg; Pkg.add(url="https://github.com/ryanlevy/ITensorIMPSTools.jl")
```

## Structure Factor Code
Structure Factor for VUMPS


The command `static_sq(q,ψ,A,α,B,β)`Computes 

$$
S^{\alpha\beta}(q) \sim \sum_{m,n} e^{iq(m-n)} \langle \psi|A^\alpha(n)B^\beta(m)|\psi\rangle
$$

where $\alpha, \beta$ are unit cell positions and the local part is removed, i.e. $A \to A - \langle A \rangle$ (as well as $B$)

### Usage

Obtain a wave function from somewhere (InfiniteCanonicalMPS)

Default usage:
```julia
julia> using ITensorIMPSTools: static_sq
       static_sq(π , ψ, "Ntot", 1, "Ntot", 1)
0.2507357963758881 - 1.1720606343985638e-27im
```
The number correspond to where in the unit cell to place the operator (essentially how to get bands)

`debug=True` can reveal convergence info from the linsolve and the total of each of the terms

```julia
julia> using ITensorIMPSTools: static_sq
       static_sq(π , ψ, "Ntot", 1, "Ntot", 1, debug=true,)
const_terms = [0.9997998645953666, 0.9997998645953666]
info.converged == 1 = true
total = -0.0020754459740376266 - 3.555571131735114e-19im
info.converged == 1 = true
total = -0.0020754459729724166 + 3.555571120014508e-19im
0.2507357963758881 - 1.1720606343985638e-27im
```

## Exact Correlators
While the structure factor code is nearly exact, sometimes its useful to visualize a real space computation of $\langle O_i O_j \rangle$. The examples use a `finite_mps` to measure this,  however, the inexact gauge of the IMPS can mean that at large distances the `finite_mps` can become non-translationally invariant. This code uses a slightly slower method to exactly compute the correlator with no dependence on gauge condition. See [Future Blog Post] for a demonstration. 

See [the ising example](https://github.com/ryanlevy/ITensorIMPSTools.jl/blob/main/examples/ising.jl) for a comparison of the fast and exact way of calculating this correlation: 
```julia
# makes a finite_mps and then puts that into correlation_matrix
O1O2_approx = correlation_fast(ψ2, Op1, Op2, stop) 

# takes advantage of translation symmetry
# and assumes there's no gauge 
O1O2_exact = correlation_slow(ψ2, Op1, Op2, stop) 
```

## Transfer Matrix Eigenvalues

There is [an example](https://github.com/ITensor/ITensorInfiniteMPS.jl/blob/main/examples/vumps/transfer_matrix_spectrum.jl) on how to obtain the transfer matrix spectrum, but this doesn't show how to obtain the charges (or QN blocks) that each eigenvalue comes from. The `getSpectrumQN` function will find the first `neigs` eigenvalues in each block, for all possible blocks. 

Default usage:
```julia
julia> using ITensorIMPSTools: getSpectrumQN
julia> eigList = getSpectrumQN(p; neigs=4)
# ...
 (QN(("Nf",3,-1),("Sz",5)) => 268, ComplexF64[-0.0005918088527009802 + 0.024003123783658695im, -0.0005918088527009802 - 0.024003123783658695im, 0.0005253780255565451 + 0.023955607507603508im, 0.0005253780255565451 - 0.023955607507603508im])
 (QN(("Nf",3,-1),("Sz",7)) => 6, ComplexF64[5.626383910108263e-7 + 0.0im, -2.710505431213761e-20 + 0.0im])
 (QN(("Nf",4,-1),("Sz",-6)) => 6, ComplexF64[-0.00012168676378200158 + 0.0im, 0.0 + 0.0im])
 (QN(("Nf",4,-1),("Sz",-4)) => 179, ComplexF64[0.00724475732832309 + 0.0im, -0.005071919371676677 + 0.0im, -0.004725498423612593 + 0.0im, 0.0034402159774847875 + 0.00014587167219837032im])
```

There is also a tool for writing this to a text file `write_spec_to_file(eigList, fname)`. 

### TODO
- Structure Factor
    - Switch to output level style reporting
    - Multi-site operators
- Transfer Matrix
    - Target specific QNs by hand  
