---
layout: page
title: Conformal Space (G(k+1,1))
navigation: 6
---

## Conformal Space

In space Gₖ we can use the _stereographic embedding_, taking an element of Gₖ to G(k+1) using the function `S` that represents the _stereographic embedding_, so to take this element to the _conformal space_ is used the function `Hm`, which represents the homogenization in the projective space with negative signature (Minkowski space), that is, instead of a new basis element that squares to `+1`, another basis element is added which is the `e-` that squares to `-1`.

Since these functions depend on `e+` and `e-`, the types based on `pbasis` were created, however, one can use `e∞` and `e∘` as basis elements for _conformal space_ which is a subspace of G(k+1,1), then the types based on `cbasis` were created for this purpose.

### `pbasis` and `cbasis`

The main difference between these two forms is in these types, the components `e∞` and `e∘` are generated from `e+` and `e-` through the following expressions.

    e∞ = e₋ + e₊
    e∘ = 0.5(e₋ - e₊)

thus,

    (e∞)² = (e₋ + e₊)² = 0
    (e∘)² = 0.5(e₋ - e₊)² = 0
    e∞ ⋅ e∘ = -1

Another difference is in the creation of a `cbasis` element, while creating an element of type `pbasis` requires a logical array and an additional logical element, to create an element of the `cbasis` type we need two more logical elements, that is,

    pbasis(e1, true, 1.0)

creates the element

    1.0e1-

of type `pbasis` and

    cbasis(e1, true, false, 1.0)

creates the element

    e1∞

which is of the type `cbasis`.

### Multivectors and blades  

The only difference between `pmultvec` and `cmultvec` as well as for `pblade` and `cblade`, is the type of basis that compose them.

For multivectors

    pmultvec(Vector{pbasis})
    cmultvec(Vector{cbasis})

and for blades

    pblade(Vector{pmultvec})
    cblade(Vector{cmultvec})

### Operations and functions

The operations and functions for these types are the same as for the Gₖ space, only the conjugate function is added, which returns the conjugate of an element.

The main difference here is that all functions depend on the types `pbasis`, `pmultvec` and `pblade`, that is, when elements of the `cbasis`, `cmultvec` or `cblade` types are operated they are automatically converted to one of the first three types and them reconverted back.
