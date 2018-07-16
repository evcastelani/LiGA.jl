```@setup liga
using Liga
layout(3)
```

# Conformal Space (\$\mathbb{G}_{k+1,1}\$)

In space \$\mathbb{G}_k\$ we can use the *stereographic embedding*, taking an element of 
\$\mathbb{G}_k\$ to \$\mathbb{G}_{k+1}\$ using the function `S` that represents the *stereographic embedding*, so to take this element to the *conformal space* is used the function ```Hm```, which represents the homogenization in the projective space with negative signature (Minkowski space), that is, instead of a new basis element that squares to `+1`, another basis element is added which is the `e-` that squares to `-1`.

Since these functions depend on `e+` and `e-`, the types based on `pbasis` were created, however, one can use `e∞` and `e∘` as basis elements for *conformal space* which is a subspace of \$\mathbb{G}_k\$ to \$\mathbb{G}_{k+1}\$, then the types based on `cbasis` were created for this purpose.

## ```pbasis``` and ```cbasis```

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

which is of the type `cbasis`. Let us try to explicit more these concepts.
We can create e∞ and e∘ directly using the ```cbasis```type.

```@repl liga
einf=cb(id, true, false, 1.0)

typeof(einf)

e0=cb(id, false, true, 1.0)

typeof(e0)

einf ∘ e0

typeof(einf ∘ e0)

```


## Multivectors and blades  

The only difference between `pmultvec` and `cmultvec` as well as for `pblade` and `cblade`, is the type of basis that compose them.

For multivectors

    pmultvec(Vector{pbasis})
    cmultvec(Vector{cbasis})

and for blades

    pblade(Vector{pmultvec})
    cblade(Vector{cmultvec})

Let us to see an example to define a *c-multivector*.

```@repl liga
u = cmultvec([cb(id, true, false, 1.0),cb(e12, true, false, 2.0)])

typeof(u)

v = cmultvec([cb(e1, true, false, 3.0),cb(e123, true, false, -2.0)])

typeof(v)
```

## Operations and functions

The operations and functions for these types are the same as for the Gₖ space, only the conjugate function is added, which returns the conjugate of an element.

The main difference here is that all functions depend on the types `pbasis`, `pmultvec` and `pblade`, that is, when elements of the `cbasis`, `cmultvec` or `cblade` types are operated they are automatically converted to one of the first three types and them reconverted back.
Following previuos idea, let us see  *cmultivectors* in action.


```@repl liga
u = cmultvec([cb(id, true, false, 3.0),cb(e13, false, false, -7.0)])

v = cmultvec([cb(e1, true, false, 3.0),cb(e123, true, false, -2.0)])

u+v #sum

u ∘ v #geometric product

u⋅v #inner product

typeof(u⋅v)

#and others...

```
