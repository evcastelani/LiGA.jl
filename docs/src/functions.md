## Main functions

There are several functions programmed in this library, however many others can be programmed from the existing ones. We will therefore highlight the main functions that form the *core* of LIGA. 
Operations such as sum and difference are intuitive and applicable to all structures (`BasisBlade, MutiVector and Blades`). These operations are pretty simple so, we can ommit here. 

### Geometric, Inner and Outer products

```@docs
geometric

inner 

outer
```

### Grade functions 
In most cases, this function handles with BasisBlade and returns a structure of the same type. The exception is the function `gradeprojection`.

```@docs
grade

gradeplus 

grademinus

gradeprojection
```

### Involutions

The first involution that we highlight is the reverse. Although this function can be used by any structure, it is important to note that it returns an object of the same structure (not a general multivector as in other functions).

```@docs

reverse(b::BasisBlade)

reverse(b::MultiVector)

reverse(A::Blade)

```

Another important involution is the `conjugate` which is avaliable just for `BasisBlade` and `MultiVector`. (**Need to extend to Blade**).

```@docs

conjugate(b::BasisBlade)

conjugate(b::MultiVector)

conjugate(A::Blade)

```

### Duality

Duality is an important concept within GA and depends explicitly on the definition of pseudo-scalar which is a variable given in Liga module by `pseudo`.
COMPLEMENTAR COM ALGO AQUI!

With the definition of pseudo scalar in mind, we can define the duality function by the geometric product of the entity (`BasisBlade`, `MultiVector` or `Blade`) with this pseudo scalar.

```@docs

dual(a::GAAbstractType)

```

### Magnitude

Is the equivalent way to define a norm in ``\mathbb{G}_{p,q}``.

```@docs

magnitude(A::GAAbstractType)

```

### Inversion

One of the main functions in this library is the inversion of the entities (`BasisBlade`, `MultiVectors` and `Blades`). The need to invert such structures often arise in a number of applications. For this reason, was defined a function with specific calculations for each structure, trying to optimize its practical use. The inversion of `MultiVector` is the most expensive operation because, it depends on tensor representation and solve a linear system. 

## Convertions
In Liga is possible convert a `Number`, `Vector`, `BasisBlade` or `Blade` into `MultiVector` type. 
```@docs
multivector
```

## Projective and Conformal Spaces
