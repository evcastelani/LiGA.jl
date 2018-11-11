# This file contains all objects created and
# used in Liga

"""
```
GAAbstractType (Liga type)
```
GAAbstractType is the most primitive type in Liga.
This type is used to indentify others types like BasisBlade
and MultiVector.
"""
abstract type GAAbstractType end

"""
```
BasisBlade (Liga type)
```
Basis Blade is important to define several table operations in Liga.
## Example
```julia-repl
julia> BasisBlade(1) #in layout(3,0,"GA")
```
returns e1.

```julia-repl
julia> e1.index #in layout(3,0,"GA")
```
returns 1 (::Int).
"""
struct BasisBlade <: GAAbstractType
	index::Int
end

"""
```
MultiVector (Liga type)
```
MultiVector is the most used object of Liga. It is defined
by BasisBlade operations like bbscalar (*). Essentially, this
object is a SparseVector.
## Example
```julia-repl
julia> v=1.0*e1  #in layout(3,0,"GA")
```
returns 1.0e1 (MultiVector not BasisBlade).
```julia-repl
julia> v.comp
```
returns
```
8-element SparseVector{Float64,Int64} with 1 stored entry:
[2]  =  1.0
```
A more directed way to create a MultiVector is to define
a vector. But in this case, an order must be respected.
```
julia>MultiVector([1.0,0.0,2.0]) #in layout(3,0,"GA")
```
returns 1.0id+2.0e2
"""
mutable struct MultiVector <: GAAbstractType
	comp::SparseVector{Float64}
end

"""
```
Blade (Liga type)
```
In progress!
"""
mutable struct Blade <: GAAbstractType
	mvector::Vector{MultiVector}
end
