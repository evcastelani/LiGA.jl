# This file contains all objects created and
# used in Liga

"""
```
GAAbstractType (Liga type)
```
GAAbstractType is the most primitive type in Liga.
This type is used to identify others types like BasisBlade
MultiVectors and Blade.
"""
abstract type GAAbstractType end

"""
```
BasisBlade (Liga type)
```
Basis Blade is an imutable type. 
It is important to define several table operations in Liga.


## Example (in layout(3,0,"GA")=𝔾₃)

```julia-repl
julia> BasisBlade(1) #
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
MultiVector is a mutable subtype of GAAbstract type. It is the 
most used object of Liga. It is defined by BasisBlade operations 
like bbscalar (*). Essentially, this object is a SparseVector.

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
#	function MultiVector(u)
#		for i=1:length(u)
#			if abs(u[i])<LigaPrecision
#				u[i]=0.0
#			end
#		end
#		new(u)
#	end
end

"""
```
Blade (Liga type)
```
A Blade is a mutable subtype of GAAbstractType.
In order to define a Blade, you need to setup tree arguments:
vectors, grade and checked. The first one is related to the set 
of linear independent vectors wich define the Blade. The grade 
is the number of vectors in this set. If you trust that the vectors 
are L. I. so you can set the argument `checked` with `true` value. 
If not, you can setup like `false`. In last case, is tested if the 
vector's set is really a L. I. set. 

```
julia> Blade([[1.0,2.0,4.0]],1,false)
```

returns [1.0,2.0,3.0]. 

```
julia> Blade([[1.0,0.0,0.0],[0.0,1.0,0.0]],2,true)
```
returns [1.0, 0.0, 0.0]∧[0.0, 1.0, 0.0].


"""
mutable struct Blade <: GAAbstractType
	scalar::Float64
	vectors::Array{Array{Float64,1},1}
	grade::Int
	checked::Bool 
	function Blade(s,v,gr,li)
		if li == false 
			lenv=length(v)
			lenvi=length(v[1])
			A=zeros(lenvi,lenv)
			for i=1:lenv
				if lenvi!=length(v[i])
					error("Some vector component with wrong dimension")
				else
					A[:,i]=v[i]
				end
			end
			if rank(A)<lenv
				error("Not LI vectors components, impossible to define a Blade")
			else 
				new(s,v,lenv,true)
			end
		else
			new(s,v,gr,true)
		end 
	end
end
