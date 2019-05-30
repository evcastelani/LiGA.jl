## Installation

This package is supported just for Julia version 1.0. Consequently, 
it uses package 3.0. Currently Liga is in [Metadata.jl](https://github.com/JuliaLang/METADATA.jl), so the package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Liga
```

## Layouts

To start to use Liga, just type:

```@repl 1
using Liga
```

The second step to use Liga is define an environment. The environment is defined through the ```layout function```.


```@repl 1
layout(3,1,"GA")
```

In this case, were created the space ``\mathbb{G}_{3,1}``. If you want to see more information about the created space, just type

```@repl 1
Liga.layoutinfo()
```

Essentially, the layout function has 4 arguments (one is an optional kwarg):
```layout(p::Int,q::Int,algebra::String; custom_repl="disable")```. 

The `p,q, algebra` arguments are mandatory and they are used to define the space with signature `p,q`. For example, if you want to create the conformalized space ``\mathbb{G}_{4,1}``, you need to type:

```@repl 1
layout(4,1,"Conformal")
```

wich is equivalent to 

```@repl 1
layout(4,1,"GA")
```

but with differents notations, see the `layout_info()` for each layout tested. 

The `custom_repl` keyword argument is optional. By default this is `disable` and changes the input prompt format (we are using [OhMyREPL package](https://github.com/KristofferC/OhMyREPL.jl) here) in order to identify the space. If you want enable this option just type  something like
```@repl 1
layout(4,1,"GA",custom_repl="enable")
```
## Basis Blade
All types introduced in Liga are subtype of  `GAAbstractType`. `BasisBlade` type is the most primitive subtype wich is used to define more complete structures, like MultiVectors. Once we define a layout, a set of `BasisBlade` is created. For example, let us create the space ``\mathbb{G}_{3}``.
 ```@repl 1
layout(3,0,"GA")
```
Now, the elements `id,e1,e2,e3,e12,e13,e23,e123` were created and they are associated to integer. 

```@repl 1
typeof(e1)
e1.index
```
In addition, operation tables related to the geometric, internal and outer products are created too. In order to acess these tables you just need to type
```@repl 1
BBgeoprodTable
BBinnerprodTable
BBouterprodTable
```
The `id` element is a BasisBlade used in place of the number 1. By default, if (for example)`e12` exists does not exist `e21`. This does not mean that it can not be set. I you type
```@repl 1
e21=-e12
typeof(e21)
e21.index
```
With this command we create a new `BasisBlade` element but, in general, it isn't a necessary definition. `BasisBlade` elements are immutable types so, alternatively, you can create a new one by 

```@repl 1
bb=BasisBlade(-5)
bb.index
```

## MultiVectors

`MultiVectors` are mutable structures associated to sparse vectors. To ilustrate, consider the following examples:

```@repl 1
using SparseArrays

MultiVector(sparsevec([2,4,6],[1.5,-2.0,3.0]))

MultiVector(sparsevec([1],[1.5]))

MultiVector(sparsevec([1,8],[-2.0,3.0]))
```
Of course, this is not a very didactic way of defining a multivector.  A simpler way to create them is through `BasisBlade`. A better (**and recommended**) way is ilustrated below. 

```@repl 1
v=1.5*e1 - 2.0*e3 + 3.0*e13

1.5*id

2.0*id+3.0*e123

```

Note that there when is defined a scalar product or sum or diference between an real scalar and BasisBlade element a new type (MultiVector) is returned. A MultiVector struct is showed just in function and order of the Basis Blade set that define the space.  

## Blades

Blades is an important structure in context of GA. In Liga, a Blade is a subtype of GAAbstract type given by a set of LI vectors. For define a Blade element you need a set of vectors, the number of vectors in this set and a bool variable. There are some functions in Liga that are devotated to work just with Blade. In other cases, a Blade is converted to MultiVector type and all functions can be used too.  

```@repl 1
Blade([[1.0,2.0],[1.0,0.0]],2,true)
```

This way is the most primitive way to define a Blade. If you need to check if the set of vectors is a LI set, you need to setup the last parameter like `false`. In order to understand the properties of a Blade, see the following example.

```@repl 1
v = Blade([[1.0,2.0],[1.0,0.0]],2,true)
v.vectors
v.grade
v.checked
```
If a Blade is defined with `false` value in `checked` parameter, the Blade object will check and if the vector set is LI a Blade with `true` value is returned.   

```@repl 1
v = Blade([[1.0,2.0],[1.0,0.0]],2,false)
v.vectors
v.grade
v.checked
```

In the case of LD vector set an error message is showed.

```@repl 1
w = Blade([[1.0,2.0],[1.0,2.0]],2,false)
```

The main reason for the introduction of this Boolean variable (`checked`) is to avoid, in situations that require further processing, unnecessary checks.
