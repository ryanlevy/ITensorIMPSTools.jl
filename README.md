# ITensorIMPSTools.jl
Tools for ITensorInfiniteMPS

## Structure Factor Code
Structure Factor for VUMPS
Computes 
$$
S(q) = \sum_{m,n} e^{iq(m-n)} \langle \psi|O_nO_m|\psi\rangle
$$

with the local part removed, i.e. $O \to O - \langle O \rangle$

### Usage

Obtain a wave function from somewhere (InfiniteCanonicalMPS)

Default usage:
```julia
julia> include("structure.jl");
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

## TODO
- Switch to output level style reporting
- Multi-site operators
