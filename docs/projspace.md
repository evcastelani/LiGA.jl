---
layout: page
title: Projective Space (G(k+1))
navigation: 5
---

## Projective Space

The projective space is a subspace of G(k+1), which is nothing more than Gₖ space with an extra dimension, so there is no need to create a new type for this.

This space can be used through embedding called homogenization, it takes a vector in Rⁿ into its embedding in the affine space.

    H(X::Vector{Float64})
    H(X::kmultvec)

Since homogenization takes an element of Rⁿ to its affine space, it is possible to enter an element of the type `kmultvec` as input parameter for the function `H`, however this must be in the Euclidean subspace of Gₖ, that is, its components must be both of grade 1.

### Operations and Functions

Since types do not change in this space, all operations and functions previously mentioned can be used in the same way.

The function `iH` can also be used, which is the inverse of the `H` function, that is, it takes a 1-vector of G(k+1) with a component in the direction `e(k+1)` different from zero and projects it into the k-dimensional Euclidean subspace of Gₖ.
