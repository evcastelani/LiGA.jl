```@setup liga
using Liga
layout(3)
```
# Projective Space (\$\mathbb{G}_{k+1}\$)

The projective space is a subspace of \$\mathbb{G}_{k+1}\$, which is nothing more than \$\mathbb{G}_k\$ space with an extra dimension, so there is no need to create a new type for this.

This space can be used through embedding called homogenization, it takes a vector in \$ \mathbb{R}^n\$into its embedding in the affine space.

    H(X::Vector{Float64})
    H(X::kmultvec)

Since homogenization takes an element of \$\mathbb{R}^n\$ to its affine space, it is possible to enter an element of the type `kmultvec` as input parameter for the function `H`, however this must be in the Euclidean subspace of \$\mathbb{G}_k\$, that is, its components must be both of grade 1.

## Operations and Functions

Since types do not change in this space, all operations and functions previously mentioned can be used in the same way. For example, considering *p-basis* elements we have
```@repl liga

v=pbasis(e12,false,2.0)

u=pbasis(e123,true,1.0)

u ∘ v

u⋅v

u^v

```
or considering multivectors in this space,

```@repl liga

u=pbasis(e123,true,1.0)+pbasis(e12,true,-2.5)

v=pbasis(e1,true,1.0)+pbasis(e12,true,5.0)

u ∘ v

u⋅v

u^v

```

The function `iH` can also be used, which is the inverse of the `H` function, that is, it takes a 1-vector of \$\mathbb{G}_{k+1}\$ with a component in the direction `e(k+1)` different from zero and projects it into the k-dimensional Euclidean subspace of \$\mathbb{G}_{k}\$.
