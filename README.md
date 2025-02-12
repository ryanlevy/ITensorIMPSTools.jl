# ITensorIMPSTools.jl
Tools for ITensorInfiniteMPS

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

### TODO
- Switch to output level style reporting
- Multi-site operators
