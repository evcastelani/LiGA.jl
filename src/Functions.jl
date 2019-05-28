# TODO: INNER AND GEOMETRIC PRODUCT DEDICATED TO BLADES 
# TODO: NEW STRUTUCT VERSOR
##############################################################################
# Multivector convertion                                                     #
##############################################################################
"""
```
multivector(b::BasisBlade) (Liga function)
multivector(a::Number) (Liga function)
multivector(v::Vector{Float64}) (Liga function)
multivector(A::Blade) (Liga function)
```
This function is used to create a MultiVector type 
from numbers, vectors, basis blade or Blades. 

## Example (in 𝔾₃₁)

```julia-repl
julia> multivector(e1)
```
returns the multivector `(1.0e1)`.

## Example (in 𝔾₃₁)

```julia-repl
julia> multivector([1.0,2.0,3.0])
```
returns the multivector `(1.0 e1)+(2.0 e2)+(3.0 e3)`.

## Example (in 𝔾₃₁)

```julia-repl
julia> multivector(Blade([[1.0,2.0,3.0]],1,true))
```
returns the multivector `(1.0 e1)+(2.0 e2)+(3.0 e3)`.
"""
function multivector(a::BasisBlade)
	ms=spzeros(GIdLayout[4])
	if a.index<0
		ms[-a.index]=-1.0
	else
		ms[a.index]=1.0
	end
	return MultiVector(ms)
end
function multivector(a::Number)
	ms=spzeros(GIdLayout[4])
	ms[1]=Float64(a)
	return MultiVector(ms)
end

# rever estas funções
function multivector(v::Vector{Float64})
	ms=spzeros(GIdLayout[4])
	for i=1:length(v)
		ms[i+1]=v[i]
	end
	return MultiVector(ms)
end
function multivector(A::Blade)
	ms=multivector(A.vectors[1]) 
	for i=2:A.grade 
		ms=outer(ms,multivector(A.vectors[i]))
	end
	if A.scalar≈1.0
		return ms
	else
		return A.scalar*ms
	end
end

##############################################################################
#  Blade convertion and constructor                                          #                                
##############################################################################
"""
```
blade(scalar,args::Vector{Float64}...) (Liga Function)
```
This function is useful for two purposes: 
 1. Is a simple to create a Blade type element;
 2. Convert a Float64 vector into Blade type element.

 ## Example (in 𝔾₃₁)

 ```julia-repl 
julia> blade(1.0,[1.0,2.0],[1.0,-1.0])
```
returns the Blade type element `1.0*[1.0,2.0]∧[1.0,-1.0]`.
"""
function blade(scalar,args::Vector{Float64}...)
	nargs=length(args)
	v=Vector{Vector{Float64}}(undef,nargs)
	for i=1:nargs
		v[i]=args[i] 
	end
	return Blade(scalar,v,nargs,true)
end

##############################################################################
# Grade like functions                                                       #
##############################################################################
"""
```
grade(b::BasisBlade) (Liga function)
```
This function is used to determine the grade of a basis blade element

## Example (in 𝔾₃₁)

```julia-repl 
julia> grade(e123)
```
return the value (Int) `3`.
"""
function grade(b::BasisBlade)
	if b.index==0
		return 0
	else
		return length(findall(LogicalBasisBlade[abs(b.index)]))
	end
end

"""
```
gradeplus(b::BasisBlade) (Liga function)
```
This function is used to determine the grade plus a Basis Blade element.
It returns an integer.
"""
function gradeplus(b::BasisBlade)
	if b.index==0
		return 0
	else
		return length(findall(LogicalBasisBlade[abs(b.index)][1:GIdLayout[1]]))
	end
end

"""
```
grademinus(b::BasisBlade) (Liga function)
```
This function is used to determine the grade minus of a basis blade element
It returns an integer.
"""
function grademinus(b::BasisBlade)
	if b.index==0
		return 0
	else
		return length(findall(LogicalBasisBlade[abs(b.index)][GIdLayout[1]+1:GIdLayout[1]+GIdLayout[2]]))
	end
end

"""
```
gradeprojection(b::BasisBlade,k::Int) (Liga function)
```
This function is used to determine the grade projection of a basis blade.
This function returns a Basis Blade element. 

## Example (in 𝔾₃₁)

```julia-repl
julia> gradeprojection(-e123,3)
```
returns the Basis Blade `-e123`
"""
function gradeprojection(b::BasisBlade,k::Int)
	if grade(b)==k
		return b
	else
		return BasisBlade(0)
	end
end

"""
```
gradeprojection(b::MultiVector,k::Int) (Liga function)
```
This function is used to determine the grade projection of a multivector.
This function returns a MultiVector element.

## Example (in 𝔾₃₁)

```julia-repl
julia> gradeprojection(1.0*id+e2-2.0*e12p,3)
```
returns the multivector  `- ( 2.0 e12+ )`.
"""
function gradeprojection(b::MultiVector,k::Int)
	ind=findnz(b.comp)
	nnv=nnz(b.comp)
	ms=spzeros(GIdLayout[4])
	for i=1:nnv
		bb=gradeprojection(BasisBlade(ind[1][i]),k)
		if bb.index!=0
			ms[bb.index]=ind[2][i]
		end
	end
	return MultiVector(ms)
end

############################################################################
# Geometric product functions                                              #
############################################################################

"""
```
geometric(a::GAAbstractType,b::GAAbstractType) (Liga function)
```
This function is used to calculate the geometric product of
GAAbstractType elements (Basis Blade, MultiVector or Blades).
It returns a MultiVector type.

## Example (in 𝔾₃)

```julia-repl
julia> geometric(e12,e12)
```
returns a multivector  `-1.0id`.

Note that we can use other structs in arguments. 

## Example (in 𝔾₃)

```julia-repl
julia> geometric(1.0+3.0*e12,e1)
```
returns `( 1.0 e1 ) - ( 3.0 e2)`.

## Example (in 𝔾₄₁ conformal)
```julia-repl
julia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0

julia> v=3.0ep

julia> geometric(u,v)
``` 
returns `( 16.5 id ) + ( 6.0 e+ ) + ( 4.5 e1+ ) - ( 9.0 e2+ ) - ( 7.5 e+- )`.

You can use the infix operator ∘ instead of geometric command.

## Example (in 𝔾₃)

```julia-repl
julia> u = 1.0+3.0*e12

julia> v = e1

julia> u ∘ v
```
"""
function geometric(a::GAAbstractType,b::GAAbstractType)
	if !isa(a,MultiVector)
		a=multivector(a)
	end
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	ind1=findnz(a.comp)
	ind2=findnz(b.comp)
	nnv1=nnz(a.comp)
	nnv2=nnz(b.comp)
	ms=spzeros(GIdLayout[4])
	for i=1:nnv1
		vl=ind1[2][i]
		for j=1:nnv2
			ind=BBgeoprodTable[ind1[1][i],ind2[1][j]].index
			if ind>0
				ms[ind]=ms[ind]+vl*ind2[2][j]
			end
			if ind<0
				ms[-ind]=ms[-ind]-vl*ind2[2][j]
			end
		end
	end
	return MultiVector(ms)
end

"""
```
∘(a::GAAbstractType,b::GAAbstractType) (Liga function)
```
It is an infix operator for geometric product. This function is 
used to calculate the geometric product of GAAbstractType elements 
(Basis Blade, MultiVector or Blades). It returns a MultiVector type.

## Example (in 𝔾₃)

```julia-repl
julia> u = 1.0+3.0*e12

julia> v = e1

julia> u ∘ v
```
"""
function Base.:∘(a::GAAbstractType,b::GAAbstractType)
	return geometric(a,b)
end

################################################################################
# inner product                                                                #
################################################################################

"""
```
inner(a::GAAbstractType,b::GAAbstractType) (Liga function)
```
This function is used to determine the inner product
of GAAbstractType elements (Basis Blade, MultiVector or Blades).

## Example (in 𝔾₃)

```julia-repl
julia> inner(1.0+3.0*e12,e1)
```
returns `- ( 3.0 e2 )`.

## Example (in 𝔾₄₁ conformal)

```julia-repl
julia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0

julia> v=3.0*ep

julia> inner(u,v)
``` 
returns `( 16.5 id )`.

An alternative way to run this functions is using the infix 
operator (⋅).
"""
function inner(a::GAAbstractType,b::GAAbstractType)
	if !isa(a,MultiVector)
		a=multivector(a)
	end
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	ind1=findnz(a.comp)
	ind2=findnz(b.comp)
	nnv1=nnz(a.comp)
	nnv2=nnz(b.comp)
    ms=spzeros(GIdLayout[4])
	for i=1:nnv1
		vl=ind1[2][i]
		for j=1:nnv2
			ind=BBinnerprodTable[ind1[1][i],ind2[1][j]].index
			if ind>0
				ms[ind]=ms[ind]+vl*ind2[2][j]
			end
			if ind<0
				ms[-ind]=ms[-ind]-vl*ind2[2][j]
			end
		end
    end
    return MultiVector(ms)
end

"""
```
⋅(a::GAAbstractType,b::GAAbstractType) (Liga function)
```
This function is used to determine the inner product
of GAAbstractType elements (Basis Blade or MultiVector).

## Example (in 𝔾₄₁ conformal)

```julia-repl
julia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0

julia> v=3.0*ep

julia> u⋅v
``` 
returns `( 16.5 id )`.
"""
function ⋅(a::GAAbstractType,b::GAAbstractType)
	return inner(a,b)
end

##############################################################################
# outer product                                                              #
##############################################################################

"""
```
outer(a::GAAbstractType,b::GAAbstractType) (Liga function)
```
This function is used to determine the outer product
of GAAbstractType elements (Basis Blade, MultiVector or Blades).

## Example (in 𝔾₃)

```julia-repl
julia> outer(1.0+3.0*e12,e1)
```
returns `( 1.0 e1 )`.

## Example (in 𝔾₄₁ conformal)

```julia-repl
julia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0

julia> v=3.0*ep

julia> outer(u,v)
``` 
returns `( 6.0 e+ ) + ( 4.5 e1+ ) - ( 9.0 e2+ ) - ( 7.5 e+- )`.

An alternative way to run this function is using the infix 
operator (^).

```julia-repl
julia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0

julia> v=3.0*ep

julia> u^v
``` 
"""
function outer(a::GAAbstractType,b::GAAbstractType)
	if isa(a,Blade) && isa(b,Blade)
		len = length(a.vectors[1])
		if a.grade+b.grade>GIdLayout[1]+GIdLayout[2]
			return Blade(0.0,[zeros(len)],0,true)
		else
			A=zeros(len,a.grade+b.grade)
			for i=1:a.grade
				A[:,i]=a.vectors[i]
			end
			for i=a.grade+1:b.grade 
				A[:,i]=b.vectors[i]
			end
			if rank(A)<a.grade+b.grade
				return Blade(0.0,[zeros(len)],0,true)
			else
				return Blade(a.scalar*b.scalar,[a.vectors;b.vectors],a.grade+b.grade,true)
			end
		end
	else
		if !isa(a,MultiVector)
			a=multivector(a)
		end
		if !isa(b,MultiVector)
			b=multivector(b)
		end
		ind1=findnz(a.comp)
		ind2=findnz(b.comp)
		nnv1=nnz(a.comp)
		nnv2=nnz(b.comp)
    	ms=spzeros(GIdLayout[4])
		for i=1:nnv1
			vl=ind1[2][i]
			for j=1:nnv2
				ind=BBouterprodTable[ind1[1][i],ind2[1][j]].index
				if ind>0
					ms[ind]=ms[ind]+vl*ind2[2][j]
				end
				if ind<0
					ms[-ind]=ms[-ind]-vl*ind2[2][j]
				end
			end
    	end
		return MultiVector(ms)
	end
end

"""
```
^(a::GAAbstractType,b::GAAbstractType) (Liga function)
```
In an infix operator to outer product. 
This function is used to determine the outer product
of GAAbstractType elements (Basis Blade, MultiVector or Blades).
	
## Example (in 𝔾₃)

```julia-repl

julia> 1.0+3.0*e12^e1

```
returns `( 1.0 e1 )`.

## Example (in 𝔾₄₁ conformal)	

```julia-repl
julia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0
	
julia> v=3.0*ep
	
julia> u^v
``` 
returns `( 6.0 e+ ) + ( 4.5 e1+ ) - ( 9.0 e2+ ) - ( 7.5 e+- )`.
	
"""
function Base.:^(a::GAAbstractType,b::GAAbstractType)
	return outer(a,b)
end

#############################################################################
# sum and difference operations                                             #
#############################################################################
function Base.:+(a::Number,b::GAAbstractType)
	a=multivector(a)
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	return MultiVector(a.comp+b.comp)
end
function Base.:+(b::GAAbstractType,a::Number)
	a=multivector(a)
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	return MultiVector(a.comp+b.comp)
end
function Base.:+(a::GAAbstractType,b::GAAbstractType)
	if !isa(a,MultiVector)
		a=multivector(a)
	end
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	return MultiVector(a.comp+b.comp)
end
function Base.:-(a::Number,b::GAAbstractType)
	a=multivector(a)
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	return MultiVector(a.comp-b.comp)
end
function Base.:-(b::GAAbstractType,a::Number)
	a=multivector(a)
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	return MultiVector(b.comp-a.comp)
end
function Base.:-(a::GAAbstractType,b::GAAbstractType)
	if !isa(a,MultiVector)
		a=multivector(a)
	end
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	return MultiVector(a.comp-b.comp)
end
#Basis blade signal
function Base.:-(a::BasisBlade)
	if a.index>0
		return BasisBlade(-a.index)
	else
		return BasisBlade(a.index)
	end
end


##############################################################################
# scalar product                                                             #
##############################################################################
"""
```
scalar(a::Number,b::GAAbstractType) (Liga function)
```
This function is used to define the most simple multivector structure.

## Example (in 𝔾₃₁)
```julia-repl
julia> u = scalar(0.13,e1)
```
returns a multivector defined by [0.0,0.13,0.0,...,0.0] (see more about
MultiVectors). Consequently it is a constructor if b is a blade basis. If b is
a MultiVector, it return the scalar product. For example,
```julia-repl
julia> u = scalar(0.13,e2-e1+2.0)
```
returns the multivector ( 0.26 id ) - (0.13 e1 ) + ( 0.13 e2 ).
An alternative way to use this command is by the simbol ``*``. For example,
```julia-repl
julia> u = 0.13*(e2-e1+2.0)
```
"""
function scalar(a::Number,b::GAAbstractType)
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	return MultiVector(a.*(b.comp))
end
function scalar(b::GAAbstractType,a::Number)
	return scalar(a,b)
end
"""
```
scalar(A::Blade,B::Blade) (Liga function)
```
This function returns a scalar product of 
two blades
"""
function scalar(A::Blade,B::Blade)
	if A.grade==B.grade
		return inner(A,B)
	else
		return 0.0
	end
end
"""
```
*(a::Number,b::GAAbstractType) (Liga function)
```
This function is used to define the most simple multivector structure.

## Example (in 𝔾₃₁)

```julia-repl
julia> u = scalar(0.13,e1)
```
returns a multivector defined by [0.0,0.13,0.0,...,0.0] (see more about
MultiVectors). Consequently it is a constructor if b is a blade basis. If b is
a MultiVector, it return the scalar product. For example,
```julia-repl
julia> u = scalar(0.13,e2-e1+2.0)
```
returns the multivector ( 0.26 id ) - (0.13 e1 ) + ( 0.13 e2 ).
An alternative way to use this command is by the simbol ``*``. For example,
```julia-repl
julia> u = 0.13*(e2-e1+2.0)
```
"""
function Base.:*(a::Number,b::GAAbstractType)
	return scalar(a,b)
end
function Base.:*(b::GAAbstractType,a::Number)
	return scalar(a,b)
end
function Base.:*(A::Blade,B::Blade)
	return scalar(A,B)
end
##############################################################################
# Involutions                                                                #
##############################################################################

"""
```
reverse(b::BasisBlade) (Liga function)
```
This function calculates the reverse of a Basis Blade.
For example,
```julia-repl
julia> u = reverse(e12)
```
returns the BasisBlade `-e12`.
"""
function reverse(b::BasisBlade)
	k=grade(b)
	if (-1.0)^(k*(k-1)/2.0)>0.0
		return b
	else
		return BasisBlade(-b.index)
	end
end

"""
```
reverse(V::MultiVector) (Liga function)
```
This function calculates the reverse of a MultiVector.
It returns a new MultiVector.

## Example (in 𝔾₃₁)

```julia-repl
julia> u=1.0*id -2.0*e1 -5.0*e∞ +3.0*e∞0 
julia> reverse(u)
```
returns the multivector `(1.0 id) - ( 2.0 e1) - ( 5.0 e∞ ) - ( 3.0 e∞0 )`.
"""
function reverse(b::MultiVector)
	ind=findnz(b.comp)
	nnv=nnz(b.comp)
	ms=spzeros(GIdLayout[4])
	for i=1:nnv
		a=reverse(BasisBlade(ind[1][i]))
		if a.index>0
			ms[a.index]=ind[2][i]
		else
			ms[-a.index]=-ind[2][i]
		end
	end
	return MultiVector(ms)
end

"""
```
reverse(B::Blade) (Liga function)
```
This function calculates the reverse of a Blade.
It returns a new Blade. 

## Example (in 𝔾₃₁)
```julia-repl
julia> B=Blade([[1.0,0.0,0.0],[0.0,1.0,0.0]],2,true)
julia> reverse(B)
```
returns `-1.0*[1.0, 0.0, 0.0]∧[0.0, 1.0, 0.0]`.
"""
function reverse(A::Blade)
	#tmp=Vector{Float64}(undef,GIdLayout[])
	A.scalar=(-1.0)^((A.grade*(A.grade-1))/2.0)
	return A   
end

"""
```
conjugate(b::BasisBlade) (Liga function)
```
This function calculates the conjugate of a Basis Blade.
	
## Example (in 𝔾₃₁)

```julia-repl
julia> u = conjugate(e12)
```
returns the BasisBlade `-e12`.
"""
function conjugate(b::BasisBlade)
	#pergunta: como fica no conformal?
	r=grademinus(b)
	if (-1.0)^r>0.0
		return reverse(b)
	else
		return -reverse(b)
	end
end

"""
```
conjugate(V::MultiVector) (Liga function)
```
This function calculates the conjugate of a MultiVector
and return a multivector.

## Example (in 𝔾₃₁)

```julia-repl
julia> u = -1.0+3.0*e1-1.0*e12

julia> conjugate(u)
```
returns the MultiVector `- ( 1.0 id ) + ( 3.0 e1 ) + ( 1.0 e12 )`.
"""
function conjugate(b::MultiVector)
	ind=findnz(b.comp)
	nnv=nnz(b.comp)
	ms=spzeros(GIdLayout[4])
	for i=1:nnv
		a=conjugate(BasisBlade(ind[1][i]))
		if a.index>0
			ms[a.index]=ind[2][i]
		else
			ms[-a.index]=-ind[2][i]
		end
	end
	return MultiVector(ms)
end

function conjugatevector(v::Vector{Float64})
	for i=2:length(v)+1
		if grademinus(BasisBlade(i))==1
			v[i-1]=-v[i-1]
		end
	end
	return v 
end

"""
```
conjugate(B::Blade) (Liga function)
```
This function calculates the conjugate of a Blade
and return a new Blade.

## Example (in 𝔾₃₂)

```julia-repl
julia> u = blade([1.0,2.0,3.0,4.0,5.0],[-1.0,-1.0,-1.0,-1.0,3.0])

julia> conjugate(u)
```
returns the Blade `[-1.0, -1.0, -1.0, 1.0, -3.0]∧[1.0, 2.0, 3.0, -4.0, -5.0]`.
"""
function conjugate(A::Blade)
	v=Vector{Vector{Float64}}(undef,A.grade)
	for i=1:A.grade 
		v[i]=conjugatevector(A.vectors[A.grade+1-i])	
	end
	return Blade(v,A.grade,true)
end

###############################################################################
#               Duality                                                       #
###############################################################################
"""
```
dual(b::GAAbstractType) (Liga function)
```
This function calculates the dual of a GAAbstractType
and returns a multivector element.

## Example (in 𝔾₃₁)

```julia-repl
julia> u = dual(e12)
```
returns the multivector ( 1.0 e+- ).

## Example (in 𝔾₃₂)

```julia-repl
julia> b=blade([1.0,0.0,3.0,0.0,5.0],[0.0,-1.0,-1.0,-1.0,0.0])
```

returns the multivector  `- ( 5.0 e123 ) - ( 5.0 e124 ) + ( 5.0 e134 ) - ( 3.0 e125 ) - ( 1.0 e235 ) - ( 3.0 e145 ) - ( 1.0 e245 ) + ( 1.0 e345 )`.
"""
function dual(b::GAAbstractType)
	if !isa(b,MultiVector)
		b=multivector(b)
	end
	#b=multivector(a)
	ind=findnz(b.comp)
	nnv=nnz(b.comp)
	ms=ind[2][1]*bbdual(BasisBlade(ind[1][1]))
	for i=2:nnv
		ms=ms+ind[2][i]*bbdual(BasisBlade(ind[1][i]))
	end
	return ms
end

"""
```
bbdual(b::BasisBlade) (Liga function)
```
This function calculates the dual of a BasisBlade type
and returns a multivector element. It is an auxiliary function.
"""
function bbdual(a::BasisBlade)
	#cj=conjugate(I)
	#if cj.index<0
	#	return -bbgeometric(a,BasisBlade(-cj.index))
	#else
	#	return  bbgeometric(a,cj)
	#end
	return geometric(a,conjugate(I))
end
#function dual(b::MultiVector)
#	ind=findnz(b.comp)
#	nnv=nnz(b.comp)
#	ms=ind[2][1]*bbdual(BasisBlade(ind[1][1]))
#	for i=2:nnv
#		ms=ms+ind[2][i]*bbdual(BasisBlade(ind[1][i]))
#	end
#	return ms
#end
#function dual(A::Blade)
#	return dual(multivector(A))
#end

###############################################################################
#               Magnitude                                                     #
###############################################################################
"""
```
magnitude(A::GAAbstractType) (Liga function)
```
This function calculates the magnitude of a general GAAbstractType
and returns a scalar. 

## Example in 𝔾₄₁ (layout(4,1,"GA"))

```julia-repl
julia> magnitude(2.0*id -2.0e3)
```
returns `2.8284271247461903`.
"""
function magnitude(A::GAAbstractType)
	if !isa(A,MultiVector)
		A=multivector(A)
	end
	return norm(A.comp,2)
end

###############################################################################
#              Inverse                                                        #
###############################################################################
"""
```
inverse(b::BasisBlade) (Liga function)
```
This function calculates the inverse of a Basis Blade.
For example, in layout(3,1,"GA")
```julia-repl
julia> u = inverse(e12)
```
returns the BasisBlade -e12. Note that conjugate (for basis blade) is equal to
inverse command.
"""
function inverse(a::BasisBlade)
	return conjugate(a)
end
"""
```
inverse(b::MultiVector) (Liga function)
```
This function calculates the inverse of a multivector element.
For example, in layout(3,1,"GA")
```julia-repl
julia>  u = inverse(e1+e12-(1.5*en))
```
returns the multivector ( 0.5714 e1 ) - ( 0.0952 e- ) - ( 0.5714 e12 ) - ( 0.7619 e2- ).

"""
function inverse(a::MultiVector)
	len=GIdLayout[4]
	A=spzeros(len,len)
	for i=1:len
		A[i,:]=a.comp'*Γ[i]
	end
	b=zeros(len)
	b[1]=1.0
	#display(det(A))
	return MultiVector(sparse(A\b))
end

#need to check! 
function inverse(A::Blade)
	B=reverse(A)
	B.vectors[1]=(1.0/geometric(A,B).comp[1])*B.vectors[1]
	return Blade(B.vectors,A.grade,true)
end
###############################################################################
# embedding                                                                   #
###############################################################################
"""
```
embedding(v::Vector{Float64}) (Liga function)
```
This function makes sense just in Conformal or Projective layout. 
It was define to embedding a real vector in corresponding enviromment. 
For exemplo, let us suppose that layout(4,1,"Conformal") was defined. 
So, it is possible embedding a vector (like [1.0,2.0,3.0]) 	in conformal 
space as the following.

```julia-repl
julia>  embedding([1.0,2.0,3.0])
```
returns the multivector ( 1.0 e1 ) + ( 2.0 e2 ) + ( 3.0 e3 ) + ( 6.5 e+ ) + ( 7.5 e- ).
If projective space was defined (like layout(4,0,"Projective")) the same vector 
becomes ( 1.0 e1 ) + ( 2.0 e2 ) + ( 3.0 e3 ) + ( 1.0 e4 ).
"""
function embedding(x::Vector{Float64})
	if GIdLayout[3]=="GA"
		error("function defined but not enable for this layout.")
	else
		ms=spzeros(GIdLayout[4])
		n=length(x)
		if n!=GIdLayout[1]-1
			error("something wrong with layout dimension and vector dimension")
		else
			if GIdLayout[3]=="Projective"
				ms[2:n+1]=x
				ms[n+2]=1.0
				return MultiVector(ms)
			end
			if GIdLayout[3]=="Conformal"
				ms[2:n+1]=x
				a=0.5*(x'*x)
				ms[n+2]=a-0.5
				ms[n+3]=a +0.5
				return MultiVector(ms)
			end
		end
	end
end
###############################################################################
# convertions                                                                 #
###############################################################################
function Base.convert(::Type{Float64}, x::MultiVector)
	return x.comp[1]
end
function Base.Float64(x::MultiVector)
	return convert(Float64,x)
end

##############################################################################
# copy function                                                              #
# ############################################################################
function Base.copy(A::MultiVector)
	return MultiVector(copy(A.comp))
end
function Base.copy(A::Blade)
	return Blade(copy(A.vectors),copy(A.grade),copy(A.checked))
end

##############################################################################
# equal function                                                             #
# ############################################################################
function Base.isequal(A::MultiVector,B::MultiVector)
	return isequal(A.comp,B.comp)
end

function Base.:(==)(A::MultiVector,B::MultiVector)
	return A.comp==B.comp
end
##############################################################################
# approx function                                                            #
##############################################################################
function Base.:≈(A::MultiVector,B::MultiVector)
	return isapprox(A.comp,B.comp,atol=LigaPrecision)
end