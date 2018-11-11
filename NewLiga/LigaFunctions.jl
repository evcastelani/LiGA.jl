##############################################################################
# Multivector convertion                                                     #
##############################################################################
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
##############################################################################
# Grade like functions                                                       #
##############################################################################
"""
```
grade(b::BasisBlade) (Liga function)
```
This function is used to determine the grade of a basis blade element
## Example
```julia-repl
julia> grade(e123)
```
return the value (Int) 3.
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
## Example
```julia-repl
julia> gradeprojection(-e123,3)
```
returns the Basis Blade -e123
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
## Example
```julia-repl
julia> gradeprojection(1.0*id+e2-2.0*e12p,3)
```
returns the multivector  - ( 2.0 e12+ ).
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
GAAbstractType elements (Basis Blade or MultiVector).
It returns a multivector.

## Example
```julia-repl
julia> geometric(e12,e12)
```
returns a multivector  -1.0id.
You can use ∘ operator too.
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
This function is used to calculate the geometric product of
GAAbstractType elements (Basis Blade or MultiVector). It returns a multivector.
## Example
```julia-repl
julia> u = 0.13*id + (0.12*e1) + (-0.15*e3) + (0.18*e12) + (0.1*e23) + (-0.29*e123)

julia> v = 0.23*id + (0.32*e1) + (-0.11*e3) + (0.128*e12) + (0.4*e23) + (-1.29*e123)

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
of GAAbstractType elements (Basis Blade or MultiVector).
## Example
```julia-repl
julia> u = 0.13*id + (0.12*e1) + (-0.15*e3) + (0.18*e12) + (0.1*e23) + (-0.29*e123)
julia> v = 0.23*id + (0.32*e1) + (-0.11*e3) + (0.128*e12) + (0.4*e23) + (-1.29*e123)
julia> inner(u,v) #layout(3,0,"GA")
```
returns -0.38224id+0.245e1+0.006759999999999997e2+0.26932e3+0.2254e12-0.2476e23.
An alternative way to run this functions is using (⋅) operator.
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
## Example
```julia-repl
julia> u = 0.13*id + (0.12*e1) + (-0.15*e3) + (0.18*e12) + (0.1*e23) + (-0.29*e123)
julia> v = 0.23*id + (0.32*e1) + (-0.11*e3) + (0.128*e12) + (0.4*e23) + (-1.29*e123)
julia> u⋅v #layout(3,0,"GA")
```
returns -0.38224id+0.245e1+0.006759999999999997e2+0.26932e3+0.2254e12-0.2476e23.
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
of GAAbstractType elements (Basis Blade or MultiVector).
## Example
```julia-repl
julia> u = 0.13*id + (0.12*e1) + (-0.15*e3) + (0.18*e12) + (0.1*e23) + (-0.29*e123)
julia> v = 0.23*id + (0.32*e1) + (-0.11*e3) + (0.128*e12) + (0.4*e23) + (-1.29*e123)
julia> outer(u,v) #layout(3,0,"GA")
```
returns 0.03id+0.0692e1-0.0488e3+0.05804e12+0.0348e13+0.075e23-0.1934e123 (like multivector). An alternative way to run this
function is using (^) operator.
"""
function outer(a::GAAbstractType,b::GAAbstractType)
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

"""
```
outer(a::GAAbstractType,b::GAAbstractType) (Liga function)
```
This function is used to determine the outer product
of GAAbstractType elements (Basis Blade or MultiVector).
## Example
```julia-repl
julia> u = 0.13*id + (0.12*e1) + (-0.15*e3) + (0.18*e12) + (0.1*e23) + (-0.29*e123)
julia> v = 0.23*id + (0.32*e1) + (-0.11*e3) + (0.128*e12) + (0.4*e23) + (-1.29*e123)
julia> u^v #layout(3,0,"GA")
```
returns 0.03id+0.0692e1-0.0488e3+0.05804e12+0.0348e13+0.075e23-0.1934e123 (like multivector).
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
For example,
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
*(a::Number,b::GAAbstractType) (Liga function)
```
This function is used to define the most simple multivector structure.
For example,
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
returns the BasisBlade -e12.
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
reverse(b::MultiVector) (Liga function)
```
This function calculates the reverse of a Multi Vector.
For example,
```julia-repl
julia> u=1.0*id -2.0*e1 -5.0*e∞ +3.0*e∞0 #in layout(3,1,"CGA")
julia> reverse(u)
```
returns the multivector 1.0id-2.0e1-5.0e∞-3.0e∞0.
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
conjugate(b::BasisBlade) (Liga function)
```
This function calculates the conjugate of a Basis Blade.
For example, in layout(3,1,"GA")
```julia-repl
julia> u = conjugate(e12)
```
returns the BasisBlade -e12.
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

###############################################################################
#               Duality                                                       #
###############################################################################
"""
```
dual(b::GAAbstractType) (Liga function)
```
This function calculates the dual of a basis blade or multivector element and returns
a multivector element.
For example, in layout(3,1,"GA")
```julia-repl
julia> u = dual(e12)
```
returns the multivector ( 1.0 e+- ).
"""
function dual(a::BasisBlade)
	#cj=conjugate(I)
	#if cj.index<0
	#	return -bbgeometric(a,BasisBlade(-cj.index))
	#else
	#	return  bbgeometric(a,cj)
	#end
	return geometric(a,conjugate(I))
end
function dual(b::MultiVector)
	ind=findnz(b.comp)
	nnv=nnz(b.comp)
	ms=ind[2][1]*dual(BasisBlade(ind[1][1]))
	for i=2:nnv
		ms=ms+ind[2][i]*dual(BasisBlade(ind[1][i]))
	end
	return ms
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

###############################################################################
# embedding                                                                   #
###############################################################################
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
				return MultiVector(ms)+0.5*(x'*x)*e∞+e0
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
