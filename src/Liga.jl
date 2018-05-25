__precompile__()

module Liga
import Base.show
import Base.copy
import Base.subtypes
importall Base.Operators
#########################################################
export kbasis, kmultvec, kblade,
pbasis, pmultvec, pblade,
cbasis, cmultvec, cblade,
kb, grade, mvectovec, copy,
geoprod, inner, outer, scalar,
bscalar, mvsum, mvreverse, dual,
magnitude, remove, reduct, bltomv,
inverse, projection, rejection,
pb, kbtopb, conjugate, prtore,
retopr, plen, mvtopmv, cb,
cbtopb, pbtocb, cbtore, retocb,
retoaffin, affine, iretoaffin,
euctoga, S, Hm, H, ovector, iH, irH,
pconformal, ipconformal, conformal,
iconformal, cbltopbl, pbltocbl, ==,
conftore,layout,tree,buildtree


#########################################################

#layout

function tree(v,niv,pos,lista,n)
    if niv<=n
        if pos==0
            aux=copy(v);
            aux[niv]=true
            niv=niv+1
            lista=push!(lista,aux)
            lista=tree(aux,niv,0,lista,n)
            lista=tree(aux,niv,1,lista,n)

        else
            niv=niv+1
            aux=copy(v);
            lista=tree(aux,niv,0,lista,n)
            lista=tree(aux,niv,1,lista,n)

        end
    end
    return lista
end

function buildtree(n)
   v=falses(n)
   lista=[v]
   for i in tree(v,1,0,[v],n)[2:length(tree(v,1,0,[v],n))]
        lista=push!(lista,i)
   end
   for i in tree(v,1,1,[v],n)[2:length(tree(v,1,1,[v],n))]
        lista=push!(lista,i)
   end
   return lista
end
""" This function is used to define the space to work
## Example
```julia-repl
julia> layout(3)
```
generates a G3 space with base 1,e1,e2,e3,e12,e13,e23,e123

"""
function layout(dim::Int)
	bn=buildtree(dim);
	abn=BitArray{1}(dim)
	abnval=0.0
	#sorting bn
	bnval=zeros(2^dim)
	for i=1:2^dim
		bnval[i]=1
		for j=1:length(bn[i])
			if bn[i][j]==true
				bnval[i]+=2^(j-1)
			end			
		end
	end
	for i=1:2^dim
		for j=i+1:2^dim
			if bnval[i]>bnval[j]
				abnval=bnval[j]
				abn=bn[j]
				bnval[j]=bnval[i]
				bn[j]=bn[i]
				bnval[i]=abnval
				bn[i]=abn
			end
		end
	end
	#display(bnval)
	#display(bn)
    bnew=Array{Array{Bool,1},1}(2^dim)
    ind=1
    for v in bn
        s=find(x->x==true,v)
        if isempty(s)
			eval(parse("const id = $v ;"));
			eval(parse("export id;"));
            bnew[ind]=eval(parse("$v "));
            ind+=1
        else
            conc=string(s[1])
            for k=2:length(s)
                conc=string(conc,s[k])
            end
			eval(parse("const e$(conc) = $v ;"));
			eval(parse("export e$(conc);"));
            bnew[ind]=eval(parse("$v ;"));
            ind+=1
        end
	end
	eval(parse("export setupbase  ;"));
	eval(parse("const setupbase = $(bnew);"));
	return ;
end



########################################################
"""
```
c(x::Vector{Bool}, a::Number)
```
Creates a scaled basis element of geometric algebra space from x and a.

```julia-repl
julia> kbasis(e13, 2.3)
2.3e13
julia> kbasis(e23)
1.0e23
```
"""
type kbasis
    e::Vector{Bool}
	scl::Number

	kbasis(a::Vector{Bool}) =   new(a, 1.0)
	kbasis(a::Vector{Bool}, b::Number) =   new(a, b)

end
#########################################################
function kb(a::Vector{Bool})
	return kbasis(a, 1.0)
end
#########################################################
function kb(a::Vector{Bool}, b::Number)
	return kbasis(a, b)
end
#########################################################
function Base.:*(a::Number, b::Vector{Bool})
	return kbasis(b,a)
end
#########################################################
"""
```
kmultvec(X::Vector{kbasis})
```

Creates a multivector in the Euclidean geometric algebra (Gₚ) composed by the sum of the coordinates of X.

```julia-repl
julia> kmultvec([kbasis(e1, 2.0), kbasis(e12, 3.3), kbasis(e123)])
2.0e1 + 3.3e12 + 1.0e123
```
"""
type kmultvec
    comp::Vector{kbasis}

	kmultvec(a::kbasis) = new([a])
	kmultvec(A::Vector{kbasis}) = new(A)
	kmultvec(a::Number) = new([kbasis(id, a)])
end
#########################################################
"""
```
grade(a::kbasis)
grade(a::pbasis)
grade(a::cbasis)
```
Returns the grade of the basis element x

```julia-repl
julia> grade(kbasis(e13, 1.5))
2
```
"""
function grade(a::kbasis)
	acu = 0
	for i=1:length(a.e)
		if a.e[i] == true
			acu = acu + 1
		end
	end
	return acu
end
#########################################################
"""
```
mvectovec(X::Vector{kmultvec})
```
Auxiliary function.

"""
function mvectovec(X::Vector{kmultvec})
	l = length(X)
	flag = 0
	for i=1:l
		m = length(X[i].comp)
		for j=1:m
			if grade(X[i].comp[j]) != 1
				flag = 1
			end
		end
	end
	if flag == 1
		error("Grade error!")
		return 0
	else
		n = length(X[1].comp[1].e)
		A = Matrix{Number}(l, n)
		fill!(A, 0.0)
		for i=1:l
			m = length(X[i].comp)
			for j=1:m
				for k=1:n
					if X[i].comp[j].e[k] == true
						A[i,k] = X[i].comp[j].scl
					end
				end
			end
		end
		return A
	end
end
#########################################################
"""
```
kblade(X::Vector{kmultvec})
```
Creates a blade in the Euclidean geometric algebra (Gp) composed by the outer product of the elements of X.

The blade is created only if the multivectors that compound X are both 1-vectors and X is a L.i set

```julia-repl
julia> a = kmultvec([kbasis(e1,2.0),kbasis(e3,4.0),kbasis(e3,4.0)])
2.0e1 + 4.0e3 + 4.0e3
julia> b = kmultvec([kbasis(e1, 1.0), kbasis(e2, 1.0), kbasis(e3,2.0)])
1.0e1 + 1.0e2 + 2.0e3
julia> A = kblade([a,b])
(2.0e1 + 4.0e3 + 4.0e3)∧(1.0e1 + 1.0e2 + 2.0e3)
```

"""
type kblade
	conj::Vector{kmultvec}
	function kblade(X::Vector{kmultvec})
		T = mvectovec(X)
		if det(T*transpose(T)) != 0
			new(X)
		else
			error("L.d")
		end
	end
end
#########################################################
function copy(a::kbasis)
	e = copy(a.e)
	esc = copy(a.scl)
	result = kbasis(e, esc)
	return result
end
#########################################################
"""
```
geoprod(a::kbasis, b::kbasis)
geoprod(X::kmultvec, Y::kmultvec)
geoprod(A::kblade, B::kblade)
geoprod(a::pbasis, b::pbasis)
geoprod(X::pmultvec, Y::pmultvec)
geoprod(A::pblade, B::pblade)
geoprod(a::cbasis, b::cbasis)
geoprod(X::cmultvec, Y::cmultvec)
geoprod(A::cblade, B::cblade)
geoprod(A::Array{kbasis,1},B::Array{kbasis,1})
geoprod(A::Array{kbasis,2},B::Array{kbasis,2})
geoprod(A::Array{kbasis,1},B::Array{kbasis,2})
geoprod(A::Array{kbasis,2},B::Array{kbasis,1})

```
Receives as input parameters elements of the same type and returns the geometric product between them.

It can be used the operator "∘ (circ)" instead of call the function.
"""
function geoprod(a::kbasis, b::kbasis)
	l = length(a.e)
	Ea = copy(a.e)
	Eb = copy(b.e)
	ab = kbasis(copy(Ea), 0.0)
	scl = 1
	eq = 0
	cont = 0
	for i=1:l-1
		if Eb[i] == true
			for j=i+1:l
				if Ea[j] == true
					cont = cont + 1
				end
			end
		end
	end
	for i=1:l
		if Ea[i] != Eb[i]
			ab.e[i] = true
		else
			ab.e[i] = false
		end
	end
	ab.scl = a.scl*b.scl*(-1)^cont
	return ab
end
#########################################################
function geoprod(A::Array{kbasis,2},B::Array{kbasis,2})
    (ma,na)=size(A)
    (mb,nb)=size(B)
    if mb != na
        error("DimensionMismatch")
    else
        C=Array{kmultvec,2}(ma,nb)
        for i=1:ma,j=1:nb
            C[i,j]=kmultvec([kbasis(id,0.0)])
            for k=1:na
                C[i,j]=C[i,j]+geoprod(A[i,k],B[k,j])
            end
        end
    end
    return C
end
function geoprod(A::Array{kbasis,1},B::Array{kbasis,2})
    ma=length(A)
    (mb,na)=size(B)
    if ma != mb
        error("DimensionMismatch")
    else
        C=Array{kmultvec,1}(na)
        for i=1:na
            C[i]=kmultvec([kb(id,0.0)])
            for j=1:ma
                C[i]=C[i]+geoprod(A[j],B[j])
            end
        end
    end
    return C
end
function geoprod(A::Array{kbasis,2},B::Array{kbasis,1})
    geoprod(A::Array{kbasis,1},B::Array{kbasis,2})
end
function geoprod(A::Array{kbasis,1},B::Array{kbasis,1})
    m=length(A)
    n=length(B)
    if m != n
        error("DimensionMismatch")
    else
        C=kmultvec([kbasis(id,0.0)])
        for i=1:m
            C=C+geoprod(A[i],B[i])
        end
    end
    return C
end
#########################################################
function Base.:∘(A::Array{kbasis},B::Array{kbasis})
    return geoprod(A,B)
end
function Base.:∘(a::kbasis,b::kbasis)
    return geoprod(a,b)
end
function Base.:∘(a::Number, b::kbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
function Base.:*(a::Number, b::kbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
"""
```
inner(a::kbasis, b::kbasis)
inner(X::kmultvec, Y::kmultvec)
inner(A::kblade, B::kblade)
inner(a::pbasis, b::pbasis)
inner(X::pmultvec, Y::pmultvec)
inner(A::pblade, B::pblade)
inner(a::cbasis, b::cbasis)
inner(X::cmultvec, Y::cmultvec)
inner(A::cblade, B::cblade)

```
Receives as input parameters elements of the same type and returns the inner product between them.

It can be used the operator "⋅ (cdot)" instead of call the function.
"""
function inner(a::kbasis, b::kbasis)
	if grade(geoprod(a,b)) == abs(grade(a) - grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kbasis(x, 0.0)
	end
end
#########################################################
function Base.:⋅(a::kbasis,b::kbasis)
    return inner(a,b)
end
#########################################################
"""
```
outer(a::kbasis, b::kbasis)
outer(X::kmultvec, Y::kmultvec)
outer(A::kblade, B::kblade)
outer(a::pbasis, b::pbasis)
outer(X::pmultvec, Y::pmultvec)
outer(A::pblade, B::pblade)
outer(a::cbasis, b::cbasis)
outer(X::cmultvec, Y::cmultvec)
outer(A::cblade, B::cblade)

```
Receives as input parameters elements of the same type and returns the outer product between them.

It can be used the operator " ^ " instead of call the function.
"""
function outer(a::kbasis, b::kbasis)
	if grade(geoprod(a,b)) == (grade(a) + grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kbasis(x, 0.0)
	end
end
#########################################################
function Base.:^(a::kbasis,b::kbasis)
    return outer(a,b)
end
#########################################################
"""
```
scalar(a::kbasis, b::kbasis)
scalar(X::kmultvec, Y::kmultvec)
scalar(A::kblade, B::kblade)
scalar(a::pbasis, b::pbasis)
scalar(X::pmultvec, Y::pmultvec)
scalar(A::pblade, B::pblade)
scalar(a::cbasis, b::cbasis)
scalar(X::cmultvec, Y::cmultvec)
scalar(A::cblade, B::cblade)

```
Receives as input parameters elements of the same type and returns the scalar product between them.

It can be used the operator " * " instead of call the function.
"""
function scalar(a::kbasis, b::kbasis)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b).scl
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kbasis(x, 0.0)
	end
end
#########################################################
function Base.:*(a::kbasis,b::kbasis)
    return scalar(a,b)
end
#########################################################
function bscalar(a::kbasis, b::kbasis)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.e))
		fill!(x, false)
		return kbasis(x, 0.0)
	end
end
#########################################################
"""
```
mvsum(a::kbasis, b::kbasis)
mvsum(X::kMulvec, Y::kmultvec)
mvsum(a::pbasis, b::pbasis)
mvsum(X::pMulvec, Y::pmultvec)
mvsum(a::cbasis, b::cbasis)
mvsum(X::cMulvec, Y::cmultvec)

```
Receives as input parameters elements of the same type and returns the sum of them.

It can be used the operators "+" or "-" instead of call the function.
"""

function mvsum(a::kbasis, b::kbasis)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, false)
	sume = kbasis(e, 0.0)
	if a.e == b.e
		sume.e = a.e
		sume.scl = a.scl + b.scl
	elseif a.e != b.e
		sume = kmultvec([a,b])
	end
	return sume
end
#########################################################
function Base.:+(a::kbasis,b::kbasis)
    return mvsum(copy(a),copy(b))
end
function Base.:-(a::kbasis,b::kbasis)
	c = copy(b)
	c.scl = -c.scl
	return mvsum(copy(a),c)
end
#########################################################
"""
```
mvreverse(a::kbasis)
mvreverse(X::kmultvec)
mvreverse(A::kblade)
mvreverse(a::pbasis)
mvreverse(A::pblade)
mvreverse(a::cbasis)
mvreverse(A::cblade)

```
Returns the reverse of the element.
"""
function mvreverse(a::kbasis)
	k = grade(a)
	e = copy(a.e)
	b = kbasis(e)
	b.scl = (-1)^(k*(k-1)/2)*a.scl
	return b
end
#########################################################
"""
```
dual(a::kbasis)
dual(X::kmultvec)
dual(A::kblade)
dual(X::pmultvec)
dual(a::pbasis)
dual(A::pblade)

```
Returns the dual element of the basis blade a, the multivetor X or the blade A.
"""
function dual(a::kbasis)
	l = length(a.e)
	e = Vector{Bool}(l)
	fill!(e, true)
	I = kbasis(e)
	Inv = mvreverse(I)
	dual = geoprod(a, Inv)
	return dual
end
#########################################################
"""
```
magnitude(a::kbasis)
magnitude(a::kmultvec)
magnitude(A::kblade)
magnitude(a::pbasis)
magnitude(a::pmultvec)
magnitude(A::pblade)
magnitude(a::cbasis)
magnitude(A::cblade)

```
Returns the magnitude of the input element.
"""
function magnitude(a::kbasis)
	b = mvreverse(a)
	return sqrt(scalar(a,b))
end
#########################################################
function ovector(X::kmultvec)
	l = length(X.comp)
	f = 1
	for i=1:l
		if grade(X.comp[i]) != 1
			f = 0
		end
	end
	if f == 1
		return true
	else
		return false
	end
end
#########################################################
#we are starting the inverse compute
"""
```
inverse(a::kmultvec)
```
Returns the inverse of a multi vector input. Be 
careful with this function because in current version
it isn't optimized!
"""
function inverse(mvet::kmultvec)
    n=length(Liga.setupbase)
    C=Array{kbasis,2}(n,n)
	D=zeros(n,n,n)
	E=zeros(n,n)
	mvetscl=zeros(n)
	for i=1:n
		for j=1:length(mvet.comp)
			if mvet.comp[j].e == Liga.setupbase[i]
				mvetscl[i]=mvet.comp[j].scl
			end
		end
	end
        for i=1:n
            for j=1:n
                C[i,j]=geoprod(kb(Liga.setupbase[i],1.0),kb(Liga.setupbase[j],1.0))
                for k=1:n
                    if C[i,j]==kb(Liga.setupbase[k],1.0)
						D[i,j,k]=1.0
						#E[i,j]=D[i,j,k]*mvetscl[k]
                    end
                    if C[i,j]==kb(Liga.setupbase[k],-1.0)
						D[i,j,k]=-1.0
						#E[i,j]=D[i,j,k]*mvetscl[k]
                    end
				end
            end
		end
		for k=1:n
			E[k,:]=mvetscl'*D[:,:,k]
		end
    #display(C)
	#display(D)
	#display(E)
	#display(inv(E))
	invE=inv(E)[:,1]
	ivect=kb(Liga.setupbase[1],invE[1])
	for i=2:n
		ivect=ivect+kb(Liga.setupbase[i],invE[i])
	end
	return ivect
	#dmvet=length(mvet.comp)
	#for i=1:2^dim
	#	for j=1:dmvet
	#		for k=1:n 
	#			#if mvet.comp.e == 
	#		end
	#	end
	#	#if mvet.comp 
	#end
end
#########################################################
"""
```
remove(A::kmultvec, b::Number)
remove(A::pmultvec, b::Number)
remove(A::cmultvec, b::Number)

```
Auxiliary function.
"""
function remove(A::kmultvec, b::Number)
	l = length(A.comp)
	aux = 0
	for i=1:l
		if A.comp[i].scl == b
			aux = aux+1
		end
	end
	B = Vector{kbasis}(l-aux)
	for i=1:l
		for j=1:l-1
			if A.comp[j].scl == b
				aux2 = A.comp[j+1]
				A.comp[j+1] = A.comp[j]
				A.comp[j] = aux2
			end
		end
	end
	for i=1:(l-aux)
		B[i] = A.comp[i]
	end
	return kmultvec(B)
end
#########################################################
function mvsum(A::kmultvec, B::kmultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Vector{kbasis}(l+m)
	for i=1:l
		AB[i] = copy(A.comp[i])
	end
	for i=l+1:l+m
		AB[i] = copy(B.comp[i-l])
	end
	for j=1:(m+l-1)
		aux = AB[j]
		for i=j+1:m+l
			if AB[i].e == aux.e && AB[i].scl != 0.0
				AB[j] = mvsum(aux, AB[i])
				aux = AB[j]
				AB[i].scl = 0.0
			end
		end
	end
	result = remove(kmultvec(AB), 0.0)
	if length(result.comp) != 0
		return result
	else
		return kmultvec(kbasis(id, 0.0))
	end
end
#########################################################
function Base.:+(A::kmultvec,B::kmultvec)
    return mvsum(copy(A),copy(B))
end
function Base.:-(A::kmultvec,B::kmultvec)
	l=length(B.comp)
	C = copy(B)
	if l != 0
		for i=1:l
			C.comp[i].scl = -C.comp[i].scl
		end
	else
		C = kmultvec([kbasis(id, false, 0.0)])
	end
	return mvsum(copy(A),C)
end
#########################################################
"""
```
reduct(A::kmultvec)
reduct(A::pmultvec)
reduct(A::cmultvec)

```
Auxiliary function.
"""
function reduct(A::kmultvec)
	l = length(A.comp)
	if l != 0
		a = kmultvec([A.comp[1]])
		B = Vector{kbasis}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, kmultvec(B))
		return C
	else
		return kmultvec([kbasis(id, 0.0)])
	end
end
#########################################################
function geoprod(A::kmultvec, B::kmultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kbasis}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = geoprod(A.comp[i], B.comp[j])
		end
	end
	p = kmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0.0)
	if length(result.comp) != 0
		return result
	else
		return kmultvec([kbasis(id, 0.0)])
	end
end
#########################################################
function Base.:∘(A::kmultvec,B::kmultvec)
    return geoprod(A,B)
end
function Base.:∘(a::Number,B::kmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
function Base.:*(a::Number,B::kmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function inner(A::kmultvec, B::kmultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kbasis}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = inner(A.comp[i], B.comp[j])
		end
	end
	p = kmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return kmultvec([kbasis(id, 0.0)])
	end
end
#########################################################
function Base.:⋅(A::kmultvec,B::kmultvec)
    return inner(A,B)
end
#########################################################
function outer(A::kmultvec, B::kmultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kbasis}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = outer(A.comp[i], B.comp[j])
		end
	end
	p = kmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return kmultvec([kbasis(id, 0.0)])
	end
end
#########################################################
function Base.:^(A::kmultvec,B::kmultvec)
    return outer(A,B)
end
#########################################################
function copy(A::kmultvec)
	l = length(A.comp)
	X = Vector{kbasis}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return kmultvec(X)
end
#########################################################
function dual(A::kmultvec)
	l = length(A.comp)
	d = Vector{kbasis}(l)
	for i=1:l
		d[i] = dual(A.comp[i])
	end
	return kmultvec(d)
end
#########################################################
function scalar(A::kmultvec, B::kmultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{kbasis}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bscalar(A.comp[i], B.comp[j])
		end
	end
	p = kmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result.comp[1].scl
	else
		return 0.0
	end
end
#########################################################
function Base.:*(A::kmultvec,B::kmultvec)
    return scalar(A,B)
end
#########################################################
function mvreverse(A::kmultvec)
	l = length(A.comp)
	B = copy(A)
	for i=1:l
		B.comp[i] = mvreverse(B.comp[i])
	end
	return B
end
#########################################################
function magnitude(A::kmultvec)
	B = mvreverse(A)
	return sqrt(scalar(A, B))
end
#########################################################
"""
```
bltomv(A::kblade)
bltomv(A::pblade)
bltomv(A::cblade)

```
Auxiliary function.
"""
function bltomv(A::kblade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################################
function geoprod(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = geoprod(X, Y)
	return result
end
#########################################################
function Base.:∘(A::kblade,B::kblade)
    return geoprod(A,B)
end
function Base.:∘(a::Number,B::kblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
function Base.:*(a::Number,B::kblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function outer(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = outer(X, Y)
	return result
end
#########################################################
function Base.:^(A::kblade,B::kblade)
    return outer(A,B)
end
#########################################################
function inner(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = inner(X, Y)
	return result
end
#########################################################
function Base.:⋅(A::kblade,B::kblade)
    return inner(A,B)
end
#########################################################
function scalar(A::kblade, B::kblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = scalar(X, Y)
	return result
end
#########################################################
function Base.:*(A::kblade,B::kblade)
    return scalar(A,B)
end
#########################################################
function mvsum(A::kblade, B::kblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return mvsum(AA,BB)
end
#########################################################
function Base.:+(A::kblade,B::kblade)
    return mvsum(copy(A),copy(B))
end
function Base.:-(A::kblade,B::kblade)
	AA = bltomv(copy(A))
	BB = bltomv(copy(B))
	l = length(BB.comp)
	for i=1:l
		BB.comp[i].scl = -BB.comp[i].scl
	end
    return mvsum(AA,BB)
end
#########################################################
function mvreverse(A::kblade)
	l = length(A.conj)
	T = Vector{kmultvec}(l)
	fill!(T, kmultvec(kbasis(id, 0.0)))
	for i=1:l
		T[i] = A.conj[l-i+1]
	end
	return kblade(T)
end
#########################################################
function magnitude(A::kblade)
	X = scalar(A, mvreverse(A))
	return sqrt(X)
end
#########################################################
function dual(A::kblade)
	X = bltomv(A)
	return dual(X)
end
#########################################################
function copy(A::kblade)
	l = length(A.conj)
	X = Vector{kmultvec}(l)
	for i=1:l
		X[i] = copy(A.conj[i])
	end
	return kblade(X)
end
#########################################################
"""
```
inverse(A::kblade)
inverse(A::pblade)
inverse(A::cblade)

```
Returns the inverse of A if A is a non-null-blade or the pseudoinverse A if A is a null-blade.
"""
function inverse(A::kblade)
	div = scalar(A, mvreverse(A))
	l = length(A.conj)
	aux = mvreverse(A)
	X = copy(aux)
	m = length(aux.conj[1].comp)
	for j=1:m
		X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
	end
	return X
end
#########################################################
"""
```
projection(A::kblade, N::kblade)
projection(A::pblade, N::pblade)

```
Returns the projection of the blade A onto the blade N.
"""
function projection(A::kblade, N::kblade)
	N2 = inverse(N)
	result = geoprod(inner(A, N2), bltomv(N))
	return result
end
#########################################################
"""
```
rejection(A::kblade, N::kblade)
rejection(A::pblade, N::pblade)

```
Returns the rejection of the blade A from the blade N.
"""
function rejection(A::kblade, N::kblade)
	aux = projection(A, N)
	l = length(aux.comp)
	for i=1:l
		aux.comp[i].scl = -(aux.comp[i].scl)
	end
	result = mvsum(bltomv(A), aux)
	return result
end
#########################################################
function show(io::IO, a::kbasis)
	if (a.scl == -1) && (grade(a) == 0)
		print(io, "-1")
	elseif (a.scl == -1) && (grade(a) != 0)
		print(io, "-")
	elseif (a.scl != 1) || (grade(a) == 0)
		print(io, "$(a.scl)")
	end
	if (a.scl == 1) && (grade(a) != 0)
		print(io, "$(a.scl)")
	end
	flag = 0
	for i=1:length(a.e)
		if a.e[i] == true
			flag = 1
		end
	end
	if flag == 1 && a.scl != 0
		print(io, "e")
	end
	for i=1:length(a.e)
		if a.e[i] == true && a.scl != 0
			print(io, "$i")
		end
	end
end
#########################################################
function show(io::IO, A::kmultvec)
	l = length(A.comp)
	if l != 0
		for i=1:(l-1)
			if A.comp[i+1].scl>=0
				print(io, "$(A.comp[i])+")
			else
				print(io, "$(A.comp[i])")
			end
		end
		print(io, "$(A.comp[l])")
	else
		print(io, 0)
	end
end
#########################################################
function show(io::IO, X::kblade)
	l = length(X.conj)
	for i=1:l-1
		print(io, "($(X.conj[i]))∧")
	end
	print(io, "($(X.conj[l]))")
end
#########################################################
function H(X::Vector{Float64})
	l = length(X)
	Y = Vector{Float64}(l+1)
	for i=1:l
		Y[i] = X[i]
	end
	Y[l+1] = 1.0
	return euctoga(Y)
end
#########################################################
function H(X::kmultvec)
	Y = mvectovec([X])
	l = length(Y)
	Z = Vector{Float64}(l)
	for i=1:l
		Z[i] = Y[i]
	end
	return H(Z)
end
#########################################################
function iH(X::kmultvec)
	Y = reduct(X)
	l = length(Y.comp)
	f = 0
	if l != 0
		m = length(Y.comp[1].e)
		en = fill!(Vector{Bool}(m), false)
		en[m] = true
		E = kbasis(en, 1.0)
		Z = Vector{kbasis}(l)
		ide = fill!(Vector{Bool}(m-1), false)
		fill!(Z, kbasis(ide, 0.0))
	end
	for i=1:l
		if Y.comp[i].e[m] == true && Y.comp[i].scl != 0.0
			f = 1
		end
	end
	if ovector(Y) == true && f == 1
		sl = ((Y ⋅ E).comp[1].scl)^(-1)
		for i=1:l
			Y.comp[i].scl = sl * Y.comp[i].scl
		end
		for i=1:l
			e = fill!(Vector{Bool}(m-1), false)
			if Y.comp[i].e[m] == false
				for j=1:m-1
					if Y.comp[i].e[j] == true
						e[j] = true
					end
				end
				Z[i] = kbasis(e, Y.comp[i].scl)
			end
		end
		return reduct(kmultvec(Z))
	elseif ovector(Y) == true && f == 0
		error("Null e(k+1) component error")
	else
		error("Grade error")
	end
end
#########################################################
function irH(X::kmultvec)
	Y = iH(X)
	l = length(Y.comp)
	Z = mvectovec([Y])
	Z2 = Vector{Float64}(l)
	for i=1:l
		Z2[i] = Z[i]
	end
	return Z2
end
#########################################################
"""
```
pbasis(a::Vector{Bool}, b::Bool, c::Number)
```
Creates a scaled basis element of geometric algebra space from a, b and c.
The b::Bool parameter represent the conponet of the basis that squares to -1.

"""
type pbasis
	eb::Vector{Bool}
	ep::Bool
	scl::Number

	pbasis(a::Vector{Bool}, b::Bool, c::Number) = new(a, b, c)
	pbasis(a::Vector{Bool}, b::Bool) = new(a, b, 1.0)
	pbasis(b::Number) = new(id, false, b)
end
#########################################################
function pb(a::Vector{Bool}, b::Bool, c::Number)
	return pbasis(a, b, c)
end
#########################################################
function pb(a::Vector{Bool}, b::Bool)
	return pbasis(a, b, 1.0)
end
#########################################################
function pb(b::Number)
	return pb(id, false, b)
end
#########################################################
"""
```
pmultvec(X::Vector{pbasis})
```
Creates a multivector in the geometric algebra (Gp,1) composed by the sum of the coordinates of X.

"""
type pmultvec
	comp::Vector{pbasis}

	pmultvec(a::pbasis) = new([a])
	pmultvec(X::Vector{pbasis}) = new(X)
end
#########################################################
function grade(a::pbasis)
	l = length(a.eb)
	cont = 0
	for i=1:l
		if a.eb[i] == true
			cont = cont + 1
		end
	end
	if a.ep == true
		cont = cont + 1
	end
	return cont
end
#########################################################
function show(io::IO, a::pbasis)
	l = length(a.eb)
	print(io, "$(a.scl)")
	if grade(a) != 0 && a.scl != 0
		print(io, "e")
		for i=1:l-1
			if a.eb[i] == true
				print(io, "$i")
			end
		end
		if a.eb[l] == true
			print(io, "₊")
		end
		if a.ep == true
			print(io, "₋")
		end
	end
end
#########################################################
function show(io::IO, A::pmultvec)
	l = length(A.comp)
	if l == 0
		print(io, "0")
	else
		print(io, "$(A.comp[1])")
		for i=2:l
			if A.comp[i].scl != 0
				if A.comp[i].scl > 0
					print(io, " +$(A.comp[i])")
				else
					print(io, " $(A.comp[i])")	
				end
			else
				print(io, " +0.0")
			end
		end
	end
end
#########################################################
"""
```
kbtopb(a::kbasis)
```
Auxiliary function.
"""
function kbtopb(a::kbasis)
	b = pbasis(a.e, false, a.scl)
	return b
end
#########################################################
function copy(a::pbasis)
	eb = copy(a.eb)
	ep = copy(a.ep)
	esc = copy(a.scl)
	result = pbasis(eb,ep, esc)
	return result
end
#########################################################
function geoprod(a::pbasis, b::pbasis)
	l = length(a.eb)
	if a.ep == b.ep == true
		aux = -1
		ep = false
	elseif a.ep == b.ep == false
		aux = 1
		ep = false
	else
		aux = 1
		ep = true
	end
	c = Vector{Bool}(l+1)
	d = Vector{Bool}(l+1)
	for i=1:l
		c[i] = a.eb[i]
		d[i] = b.eb[i]
	end
	c[l+1] = a.ep
	d[l+1] = b.ep
	f = geoprod(kbasis(c, a.scl), kbasis(d, b.scl))
	scl = aux*f.scl
	g = Vector{Bool}(l)
	for i=1:l
		g[i] = f.e[i]
	end
	result = pbasis(g, ep, scl)
	return result
end
#########################################################
function Base.:∘(a::pbasis,b::pbasis)
    return geoprod(a,b)
end
function Base.:∘(a::Number, b::pbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
function Base.:*(a::Number, b::pbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function inner(a::pbasis, b::pbasis)
	if grade(geoprod(a,b)) == abs(grade(a) - grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pbasis(x, false, 0.0)
	end
end
#########################################################
function Base.:⋅(a::pbasis,b::pbasis)
    return inner(a,b)
end
#########################################################
function outer(a::pbasis, b::pbasis)
	if grade(geoprod(a,b)) == (grade(a) + grade(b))
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pbasis(x, false, 0.0)
	end
end
#########################################################
function Base.:^(a::pbasis,b::pbasis)
    return outer(a,b)
end
#########################################################
function scalar(a::pbasis, b::pbasis)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b).scl
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pbasis(x, false, 0.0)
	end
end
#########################################################
function bscalar(a::pbasis, b::pbasis)
	if grade(geoprod(a,b)) == 0
		return geoprod(a,b)
	else
		x = Vector{Bool}(length(a.eb))
		fill!(x, false)
		return pbasis(x, false, 0.0)
	end
end
#########################################################
function Base.:*(a::pbasis,b::pbasis)
    return scalar(a,b)
end
#########################################################
function mvsum(a::pbasis, b::pbasis)
	l = length(a.eb)
	eb = Vector{Bool}(l)
	fill!(eb, false)
	if a.eb == b.eb && a.ep == b.ep
		scl = a.scl + b.scl
	 	return pbasis(a.eb, a.ep, scl)
	else
		return pmultvec([a,b])
	end
end
#########################################################
function Base.:+(a::pbasis,b::pbasis)
    return mvsum(copy(a),copy(b))
end
function Base.:-(a::pbasis,b::pbasis)
	c = copy(b)
	c.scl = -c.scl
	return mvsum(copy(a),c)
end
#########################################################
function mvreverse(a::pbasis)
	k = grade(a)
	eb = copy(a.eb)
	ep = copy(a.ep)
	b = pbasis(eb,ep, 1)
	b.scl = (-1)^(k*(k-1)/2)*a.scl
	return b
end
#########################################################
"""
```
conjugate(a::pbasis)
conjugate(X::pmultvec)
conjugate(A::pblade)
conjugate(A::cblade)
```
Returns the basis element a, multivector X or blade A.
"""
function conjugate(a::pbasis)
	grm = 0
	if a.ep == true
		grm = 1
	end
	b = mvreverse(a)
	b.scl = b.scl*(-1)^grm
	return b
end
#########################################################
function dual(a::pbasis)
	l = length(a.eb)
	eb = Vector{Bool}(l)
	fill!(eb, true)
	I = pbasis(eb,true,1.0)
	Inv = mvreverse(I)
	dual = geoprod(a, Inv)
	return dual
end
#########################################################
"""
```
prtore(a::pbasis)
```
Auxiliary function.
"""
function prtore(a::pbasis)
	l = length(a.eb)
	e = Vector{Bool}(l+1)
	for i=1:l
		e[i] = a.eb[i]
	end
	e[l+1] = a.ep
	b = kbasis(e, a.scl)
	return b
end
#########################################################
"""
```
retopr(a::kbasis)
```
Auxiliary function.
"""
function retopr(a::kbasis)
	l = length(a.e)
	eb = Vector{Bool}(l-1)
	ep = a.e[l]
	for i=1:l-1
		eb[i] = a.e[i]
	end
	b = pbasis(eb,ep,a.scl)
end
#########################################################
"""
```
plen(A::pmultvec)
plen(A::pblade)
```
Auxiliary function.
"""
function plen(A::pmultvec)
	l = length(A.comp)
	if l != 0
		r = length(A.comp[1].eb)
		eb = fill!(Vector{Bool}(r), false)
		return eb
	else
		return id
	end
end
#########################################################
"""
```
mvtopmv(A::kmultvec)
```
Auxiliary function.
"""
function mvtopmv(A::kmultvec)
	l = length(A.comp)
	X = Vector{pbasis}(l)
	for i=1:l
		X[i] = kbtopb(A.comp[i])
	end
	return pmultvec(X)
end
#########################################################
function ovector(X::pmultvec)
	l = length(X.comp)
	f = 1
	for i=1:l
		if grade(X.comp[i]) != 1
			f = 0
		end
	end
	if f == 1
		return true
	else
		return false
	end
end
#########################################################
function remove(A::pmultvec, b::Number)
	l = length(A.comp)
	aux = 0
	for i=1:l
		if A.comp[i].scl == b
			aux = aux+1
		end
	end
	B = Vector{pbasis}(l-aux)
	for i=1:l
		for j=1:l-1
			if A.comp[j].scl == b
				aux2 = A.comp[j+1]
				A.comp[j+1] = A.comp[j]
				A.comp[j] = aux2
			end
		end
	end
	for i=1:(l-aux)
		B[i] = A.comp[i]
	end
	return pmultvec(B)
end
#########################################################
function copy(A::pmultvec)
	l = length(A.comp)
	X = Vector{pbasis}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return pmultvec(X)
end
#########################################################
function mvsum(A::pmultvec, B::pmultvec)
	l = length(A.comp)
	m = length(B.comp)
	A2 = Vector{kbasis}(l)
	B2 = Vector{kbasis}(m)
	for i=1:l
		A2[i] = prtore(A.comp[i])
	end
	for i=1:m
		B2[i] = prtore(B.comp[i])
	end
	A2 = kmultvec(A2)
	B2 = kmultvec(B2)
	AB2 = mvsum(A2,B2)
	n = length(AB2.comp)
	AB = Vector{pbasis}(n)
	for i=1:n
		AB[i] = retopr(AB2.comp[i])
	end
	AB = pmultvec(AB)
	result = remove(AB, 0.0)
	if length(result.comp) != 0
		return result
	else
		return pmultvec(pbasis(id,false, 0.0))
	end
end
#########################################################
function Base.:+(A::pmultvec,B::pmultvec)
    return mvsum(copy(A),copy(B))
end
function Base.:-(A::pmultvec,B::pmultvec)
	l=length(B.comp)
	C = copy(B)
	if l != 0
		for i=1:l
			C.comp[i].scl = -C.comp[i].scl
		end
	else
		eb = plen(C)
		C = pmultvec([pbasis(eb, false, 0.0)])
	end
	return mvsum(copy(A),C)
end
#########################################################
function reduct(A::pmultvec)
	l = length(A.comp)
	eb = plen(A)
	if l != 0
		a = pmultvec([A.comp[1]])
		B = Vector{pbasis}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, pmultvec(B))
		return C
	else
		return pmultvec([pbasis(eb, false, 0.0)])
	end
end
#########################################################
function geoprod(A::pmultvec, B::pmultvec)
	l = length(A.comp)
	m = length(B.comp)
	eb = plen(A)
	AB = Matrix{pbasis}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = geoprod(A.comp[i], B.comp[j])
		end
	end
	p = pmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pmultvec([pbasis(eb, false, 0.0)])
	end
end
#########################################################
function Base.:∘(A::pmultvec,B::pmultvec)
    return geoprod(A,B)
end
function Base.:∘(a::Number,B::pmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
function Base.:*(a::Number,B::pmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function inner(A::pmultvec, B::pmultvec)
	l = length(A.comp)
	m = length(B.comp)
	eb = plen(A)
	AB = Matrix{pbasis}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = inner(A.comp[i], B.comp[j])
		end
	end
	p = pmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pmultvec([pbasis(eb, false, 0.0)])
	end
end
#########################################################
function Base.:⋅(A::pmultvec,B::pmultvec)
    return inner(A,B)
end
#########################################################
function outer(A::pmultvec, B::pmultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{pbasis}(l,m)
	eb = plen(A)
	for i=1:l
		for j=1:m
			AB[i,j] = outer(A.comp[i], B.comp[j])
		end
	end
	p = pmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result
	else
		return pmultvec([pbasis(eb, false, 0.0)])
	end
end
#########################################################
function Base.:^(A::pmultvec,B::pmultvec)
    return outer(A,B)
end
#########################################################
function scalar(A::pmultvec, B::pmultvec)
	l = length(A.comp)
	m = length(B.comp)
	AB = Matrix{pbasis}(l,m)
	for i=1:l
		for j=1:m
			AB[i,j] = bscalar(A.comp[i], B.comp[j])
		end
	end
	p = pmultvec(reshape(AB,l*m))
	p = reduct(p)
	result = remove(p, 0)
	if length(result.comp) != 0
		return result.comp[1].scl
	else
		return 0.0
	end
end
#########################################################
function Base.:*(A::pmultvec,B::pmultvec)
    return scalar(A,B)
end
#########################################################
function conjugate(A::pmultvec)
	l = length(A.comp)
	B = Vector{pbasis}(l)
	for i=1:l
		B[i] = conjugate(A.comp[i])
	end
	return pmultvec(B)
end
#########################################################
function magnitude(A::pmultvec)
	scl = scalar(A,conjugate(A))
	result = sqrt(scl)
	return result
end
#########################################################
function dual(A::pmultvec)
	k = length(A.comp)
	d = Vector{pbasis}(k)
	for i=1:k
		d[i] = dual(A.comp[i])
	end
	return pmultvec(d)
end
#########################################################
"""
```
mvectovec(X::Vector{pmultvec})
```
Auxiliary function.
"""
function mvectovec(X::Vector{pmultvec})
	l = length(X)
	flag = 0
	for i=1:l
		m = length(X[i].comp)
		for j=1:m
			if grade(X[i].comp[j]) != 1
				flag = 1
			end
		end
	end
	if flag == 1
		error("Grade error (mvectovec function)")
		return 0
	else
		n = length(X[1].comp[1].eb)
		A = Matrix{Number}(l, n+1)
		fill!(A, 0.0)
		for i=1:l
			m = length(X[i].comp)
			for j=1:m
				for k=1:n
					if X[i].comp[j].eb[k] == true
						A[i,k] = X[i].comp[j].scl
					end
				end
				if X[i].comp[j].ep == true
					A[i,n+1] = X[i].comp[j].scl
				end
			end
		end
		return A
	end
end
#########################################################
"""
```
pblade(X::Vector{pmultvec})
```
Creates a blade in the geometric algebra (Gp,1) composed by the outer product of the elements of X.

The blade is created only if the multivectors that compound X are both 1-vectors and X is a L.i set

"""
type pblade
	conj::Vector{pmultvec}
	function pblade(X::Vector{pmultvec})
		T = mvectovec(X)
		if det(T*transpose(T)) != 0
			new(X)
		else
			error("L.d(pblade conversion)")
		end
	end
end
#########################################################
function show(io::IO, A::pblade)
	if length(A.conj) != 0
		print(io, "($(A.conj[1]))")
		for i=2:length(A.conj)
			print(io, "∧($(A.conj[i]))")
		end
	else
		print(io, "0")
	end
end
#########################################################
function bltomv(A::pblade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################################
function geoprod(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = geoprod(X, Y)
	return result
end
#########################################################
function Base.:∘(A::pblade, B::pblade)
    return geoprod(A,B)
end
function Base.:∘(a::Number,B::pblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
function Base.:*(a::Number,B::pblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function outer(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = outer(X, Y)
	return result
end
#########################################################
function Base.:^(A::pblade, B::pblade)
    return outer(A,B)
end
#########################################################
function inner(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = inner(X, Y)
	return result
end
#########################################################
function Base.:⋅(A::pblade, B::pblade)
    return inner(A,B)
end
#########################################################
function scalar(A::pblade, B::pblade)
	X = bltomv(A)
	Y = bltomv(B)
	result = scalar(X, Y)
	return result
end
#########################################################
function Base.:*(A::pblade, B::pblade)
    return scalar(A,B)
end
#########################################################
function plen(A::pblade)
	B = bltomv(A)
	l = length(B.comp)
	if l != 0
		r = length(B.comp[1].eb)
		eb = fill!(Vector{Bool}(r), false)
		return eb
	else
		return id
	end
end
#########################################################
function mvreverse(A::pblade)
	l = length(A.conj)
	eb = plen(A)
	T = Vector{pmultvec}(l)
	fill!(T, pmultvec([pbasis(eb, false, 0.0)]))
	for i=1:l
		T[i] = A.conj[l-i+1]
	end
	return pblade(T)
end
#########################################################
function conjugate(A::pblade)
	l = length(A.conj)
	T = Vector{pmultvec}(l)
	eb = plen(A)
	fill!(T, pmultvec([pbasis(eb, false, 0.0)]))
	for i=1:l
		T[i] = conjugate(A.conj[l-i+1])
	end
	return pblade(T)
end
#########################################################
function magnitude(A::pblade)
	X = scalar(A, conjugate(A))
	return sqrt(X)
end
#########################################################
function dual(A::pblade)
	X = bltomv(A)
	return dual(X)
end
#########################################################
function copy(A::pblade)
	l = length(A.conj)
	X = Vector{pmultvec}(l)
	for i=1:l
		X[i] = copy(A.conj[i])
	end
	return pblade(X)
end
#########################################################
function inverse(A::pblade)
	if geoprod(A,A).comp[1].scl != 0
		div = scalar(A, mvreverse(A))
		l = length(A.conj)
		aux = mvreverse(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	else
		div = inner(A, conjugate(A)).comp[1].scl
		l = length(A.conj)
		aux = conjugate(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	end
	return X
end
#########################################################
function projection(A::pblade, N::pblade)
	X = geoprod(inner(A, inverse(N)), bltomv(N))
	return X
end
#########################################################
function rejection(A::pblade, N::pblade)
	X = projection(A, N)
	l = length(X.comp)
	for i=1:l
		X.comp[i].scl = (-1)*(X.comp[i].scl)
	end
	result = mvsum(bltomv(A), X)
	return result
end
#########################################################
"""
```
cbasis(a::Vector{Bool},b::Bool,c::Bool,d::Number)
```
Creates a scaled basis element of geometric algebra space (Gp+1,1) from a, b, c and d using e∞ and e∘ as basis elements instead of the regular basis.

"""
type cbasis
	er::Vector{Bool}
	ei::Bool
	eo::Bool
	scl::Number

	function cbasis(a::Vector{Bool},b::Bool,c::Bool,d::Number)
		if d != 0.0
			return new(a,b,c,d)
		else
			l = length(a)
			x = Vector{Bool}(l)
			fill!(x, false)
			return new(x, false, false, 0.0)
		end
	end
end
#########################################################
function cb(a::Vector{Bool},b::Bool,c::Bool,d::Number)
	return cbasis(a,b,c,d)
end
#########################################################
"""
```
cmultvec(X::Vector{cbasis})
```
Creates a multivector in the geometric algebra (Gp+1,1) composed by the sum of the coordinates of X.

"""
type cmultvec
	comp::Vector{cbasis}

	cmultvec(a::cbasis) = new([a])
	cmultvec(X::Vector{cbasis}) = new(X)
end
#########################################################
function grade(a::cbasis)
	l = length(a.er)
	cont = 0
	for i=1:l
		if a.er[i] == true
			cont = cont + 1
		end
	end
	if a.ei == true
		cont = cont + 1
	end
	if a.eo == true
		cont = cont + 1
	end
	return cont
end
#########################################################
function show(io::IO, a::cbasis)
	l = length(a.er)
	print(io, "$(a.scl)")
	if grade(a) != 0
		print(io, "e")
		for i=1:l
			if a.er[i] == true
				print(io, "$i")
			end
		end
		if a.ei == true
			print(io, "\∞")
		end
		if a.eo == true
			print(io, "∘")
		end
	end
end
#########################################################
function show(io::IO, A::cmultvec)
	l = length(A.comp)
	if l == 0
		print(io, "0")
	else
		print(io, "$(A.comp[1])")
		for i=2:l
			if A.comp[i].scl != 0
				print(io, " + $(A.comp[i])")
			else
				print(io, " + 0.0")
			end
		end
	end
end
#########################################################
function copy(a::cbasis)
	er = copy(a.er)
	ei = copy(a.ei)
	eo = copy(a.eo)
	scl = copy(a.scl)
	return cbasis(er,ei,eo,scl)
end
#########################################################
function copy(A::cmultvec)
	l = length(A.comp)
	X = Vector{cbasis}(l)
	for i=1:l
		X[i] = copy(A.comp[i])
	end
	return cmultvec(X)
end
#########################################################
function ovector(X::cmultvec)
	l = length(X.comp)
	f = 1
	for i=1:l
		if grade(X.comp[i]) != 1
			f = 0
		end
	end
	if f == 1
		return true
	else
		return false
	end
end
#########################################################
function mvsum(a::cbasis, b::cbasis)
	l = length(a.er)
	er = Vector{Bool}(l)
	fill!(er, false)
	if a.er == b.er && a.ei == b.ei && a.eo == b.eo
		scl = a.scl + b.scl
	 	return cmultvec(cbasis(a.er, a.ei,a.eo, scl))
	else
		return cmultvec([a,b])
	end
end
#########################################################
function Base.:+(a::cbasis, b::cbasis)
    return mvsum(copy(a),copy(b))
end
function Base.:-(a::cbasis, b::cbasis)
	c = copy(b)
	c.scl = - c.scl
    return mvsum(copy(a),b)
end
#########################################################
function remove(A::cmultvec, b::Number)
	l = length(A.comp)
	aux = 0
	for i=1:l
		if A.comp[i].scl == b
			aux = aux+1
		end
	end
	B = Vector{cbasis}(l-aux)
	if aux != 0
		for i=1:l
			for j=1:l-1
				if A.comp[j].scl == b
					aux2 = A.comp[j+1]
					A.comp[j+1] = A.comp[j]
					A.comp[j] = aux2
				end
			end
		end
		for i=1:(l-aux)
			B[i] = A.comp[i]
		end
		if length(B) != 0
			return cmultvec(B)
		else
			return cmultvec(cbasis(id, false, false, 0.0))
		end
	else
		return A
	end
end
#########################################################
function mvsum(A::cmultvec, B::cmultvec)
	l = length(A.comp)
	m = length(B.comp)
	A2 = Vector{kbasis}(l)
	B2 = Vector{kbasis}(m)
	for i=1:l
		A2[i] = cbtore(A.comp[i])
	end
	for i=1:m
		B2[i] = cbtore(B.comp[i])
	end
	A2 = kmultvec(A2)
	B2 = kmultvec(B2)
	AB2 = mvsum(A2,B2)
	n = length(AB2.comp)
	AB = Vector{cbasis}(n)
	for i=1:n
		AB[i] = retocb(AB2.comp[i])
	end
	AB = cmultvec(AB)
	result = remove(AB, 0.0)
	if length(result.comp) != 0
		return result
	else
		return cmultvec(cbasis(id,false,false, 0.0))
	end
end
#########################################################
function Base.:+(A::cmultvec,B::cmultvec)
    return mvsum(copy(A),copy(B))
end
function Base.:-(A::cmultvec,B::cmultvec)
	l=length(B.comp)
	C = copy(B)
	if l != 0
		for i=1:l
			C.comp[i].scl = -C.comp[i].scl
		end
	else
		C = cmultvec([cbasis(id,false, false, 0.0)])
	end
	return mvsum(copy(A),C)
end
#########################################################
function reduct(A::cmultvec)
	l = length(A.comp)
	if l != 0
		a = cmultvec([A.comp[1]])
		B = Vector{cbasis}(l-1)
		for i=2:l
			B[i-1] = A.comp[i]
		end
		C = mvsum(a, cmultvec(B))
		return C
	else
		return cmultvec([cbasis(id,false, false, 0.0)])
	end
end
#########################################################
"""
```
cbtopb(a::cbasis)
```
Auxiliary function.
"""
function cbtopb(a::cbasis)
	scl = a.scl
	l = length(a.er)
	eb1 = vcat(copy(a.er), true)
	eb2 = vcat(copy(a.er), false)
	if a.ei == a.eo == true
		result = mvsum(pbasis(eb1, true, 1.0), pbasis(eb2, false, -1.0))
	elseif a.ei == true && a.eo == false
		result = mvsum(pbasis(eb1, false, 1.0), pbasis(eb2, true, 1.0))
	elseif a.ei == false && a.eo == true
		result = mvsum(pbasis(eb1, false, -0.5), pbasis(eb2, true, 0.5))
	elseif a.ei == a.eo == false
		result = pmultvec(pbasis(eb2, false, 1.0))
	end
	l = length(result.comp)
	for i=1:l
		result.comp[i].scl = scl*(result.comp[i].scl)
	end
	return result
end
#########################################################
"""
```
pbtocb(a::pbasis)
```
Auxiliary function.
"""
function pbtocb(a::pbasis)
	l = length(a.eb)
	er = Vector{Bool}(l-1)
	for i=1:l-1
		er[i] = a.eb[i]
	end
	scl = a.scl
	if a.eb[l] == true && a.ep == false
		result = mvsum(cbasis(er,true,false,0.5*scl),cbasis(er,false,true,-scl))
	elseif a.eb[l] == false && a.ep == true
		result = mvsum(cbasis(er,true,false,0.5*scl),cbasis(er,false,true,scl))
	elseif a.eb[l] == a.ep == true
		result = mvsum(cbasis(er,true,true,scl),cbasis(er,false,false,scl))
	else
		result = cmultvec(cbasis(er, false, false, scl))
	end
	return result
end
#########################################################
function cbtopb(A::cmultvec)
	l = length(A.comp)
	X = Vector{pbasis}(2*l)
	m = length(A.comp[1].er)
	eb = Vector{Bool}(m+1)
	fill!(eb, false)
	for i=1:l
		if length(cbtopb(A.comp[i]).comp) == 2
			X[2*i-1] = cbtopb(A.comp[i]).comp[1]
			X[2*i] = cbtopb(A.comp[i]).comp[2]
		elseif length(cbtopb(A.comp[i]).comp) == 1
			X[2*i-1] = cbtopb(A.comp[i]).comp[1]
			X[2*i] = pbasis(eb,false,0.0)
		else
			error("Unexpected error")
		end
	end
	return reduct(pmultvec(X))
end
#########################################################
function pbtocb(A::pmultvec)
	l = length(A.comp)
	X = Vector{cbasis}(2*l)
	m = length(A.comp[1].eb)

	eb = Vector{Bool}(m+1)
	fill!(eb, false)
	for i=1:l
		if length(pbtocb(A.comp[i]).comp) == 2
			X[2*i-1] = pbtocb(A.comp[i]).comp[1]
			X[2*i] = pbtocb(A.comp[i]).comp[2]
		elseif length(pbtocb(A.comp[i]).comp) == 1
			X[2*i-1] = pbtocb(A.comp[i]).comp[1]
			X[2*i] = cbasis(id,false,false,0.0)
		else
			error("Unexpected error")
		end
	end
	return reduct(cmultvec(X))
end
#########################################################
function geoprod(a::cbasis, b::cbasis)
	A = cbtopb(a)
	B = cbtopb(b)
	AB = geoprod(A,B)
	return pbtocb(AB)
end
#########################################################
function Base.:∘(a::cbasis, b::cbasis)
    return geoprod(a,b)
end
function Base.:∘(a::Number, b::cbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
function Base.:*(a::Number, b::cbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function inner(a::cbasis, b::cbasis)
	A = cbtopb(a)
	B = cbtopb(b)
	AB = inner(A,B)
	return pbtocb(AB)
end
#########################################################
function Base.:⋅(a::cbasis, b::cbasis)
    return inner(a,b)
end
#########################################################
function outer(a::cbasis, b::cbasis)
	A = cbtopb(a)
	B = cbtopb(b)
	AB = outer(A,B)
	return pbtocb(AB)
end
#########################################################
function Base.:^(a::cbasis, b::cbasis)
    return outer(a,b)
end
#########################################################
function scalar(a::cbasis, b::cbasis)
	A = cbtopb(a)
	B = cbtopb(b)
	AB = scalar(A,B)
	return pbtocb(AB)
end
#########################################################
function Base.:*(a::cbasis, b::cbasis)
    return scalar(a,b)
end
#########################################################
"""
```
cbtore(a::cbasis)
```
Auxiliary function.
"""
function cbtore(a::cbasis)
	l = length(a.er)
	e = Vector{Bool}(l+2)
	for i=1:l
		e[i] = a.er[i]
	end
	e[l+1] = a.ei
	e[l+2] = a.eo
	b = kbasis(e, a.scl)
	return b
end
#########################################################
"""
```
retocb(a::kbasis)
```
Auxiliary function.
"""
function retocb(a::kbasis)
	l = length(a.e)
	er = Vector{Bool}(l-2)
	for i=1:l-2
		er[i] = a.e[i]
	end
	ei = a.e[l-1]
	eo = a.e[l]
	b = cbasis(er, ei, eo, a.scl)
	return b
end
#########################################################
function geoprod(A::cmultvec,B::cmultvec)
	X = cbtopb(A)
	Y = cbtopb(B)
	XY = geoprod(X,Y)
	return pbtocb(XY)
end
#########################################################
function Base.:∘(A::cmultvec,B::cmultvec)
    return geoprod(A,B)
end
function Base.:∘(a::Number,B::cmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
function Base.:*(a::Number,B::cmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function inner(A::cmultvec,B::cmultvec)
	X = cbtopb(A)
	Y = cbtopb(B)
	XY = inner(X,Y)
	return pbtocb(XY)
end
#########################################################
function Base.:⋅(A::cmultvec,B::cmultvec)
    return inner(A,B)
end
#########################################################
function outer(A::cmultvec,B::cmultvec)
	X = cbtopb(A)
	Y = cbtopb(B)
	XY = outer(X,Y)
	return pbtocb(XY)
end
#########################################################
function Base.:^(A::cmultvec,B::cmultvec)
    return outer(A,B)
end
#########################################################
function scalar(A::cmultvec,B::cmultvec)
	X = cbtopb(A)
	Y = cbtopb(B)
	XY = scalar(X,Y)
	return pbtocb(XY)
end
#########################################################
function Base.:*(A::cmultvec,B::cmultvec)
    return scalar(A,B)
end
#########################################################
"""
```
retoaffin(a::Vector{Float64})
```
Takes a vector of Euclidean space into a vector of affine space.
"""
function retoaffin(a::Vector{Float64})
	l = length(a)
	A = Vector{kbasis}(l+1)
	aux = Vector{Bool}(l+1)
	fill!(aux, false)
	for i=1:l
		aux[i] = true
		A[i] = copy(kbasis(aux, a[i]))
		aux[i] = false
	end
	aux[l+1] = true
	A[l+1] = copy(kbasis(aux, 1.0))
	return kmultvec(A)
end
#########################################################
"""
```
affine(A::kmultvec)
```
Takes a vector of Euclidean subspace of Gp into a vector of affine space.
"""
function affine(A::kmultvec)
	X = copy(reduct(A))
	l = length(X.comp)
	m = length(X.comp[1].e)
	f1 = 1
	f2 = 0
	for i=1:l
		if grade(A.comp[i]) != 1
			f1 = 0
		end
		if length(A.comp[i].e[m]) == true
			f2 = 1
		end
	end
	if f1 == f2 == 1
		aux = Vector{Bool}(m)
		fill!(aux, false)
		aux[m] = true
		scl = inner(A, kmultvec(kbasis(aux, 1.0))).comp[1].scl
		a = copy(X)
		for i=1:l
			X.comp[i].scl = (X.comp[i].scl)/scl
		end
		return X
	elseif f1 == 1 && f2 == 0
		error("Image error")
	elseif f1 == 0 && f2 == 1
		error("Grade error")
	end
end
#########################################################
"""
```
iretoaffin(A::kmultvec)
```
The inverse transformation of "affine" function.

Takes a vector in affine space back into a vector of Euclidean space subspace of Gp.
"""
function iretoaffin(A::kmultvec)
	l = length(A.comp)
	m =	length(A.comp[1].e)
	f1 = 1
	f2 = 0
	aux = 0
	X = reduct(copy(A))
	for i=1:l
		if grade(X.comp[i]) != 1
			f1 = 0
		end
		if length(X.comp[i].e[m]) == true
			f2 = 1
			aux = i
		end
	end
	if f1 == f2 == 1
		X = affine(X)
		e = fill!(Vector{Bool}(m),false)
		X.comp[aux] = kbasis(e, 0.0)
		return reduct(X)
	elseif f1 == 1 && f2 == 0
		error("Image error")
	elseif f1 == 0 && f2 == 1
		error("Grade error")
	end
end
#########################################################
"""
```
euctoga(x::Vector{Float64})
```
Auxiliary function.
"""
function euctoga(x::Vector{Float64})
	l = length(x)
	X = Vector{kbasis}(l)
	e = fill!(Vector{Bool}(l), false)
	for i=1:l
		e[i] = true
		X[i] = kbasis(copy(e), x[i])
		e[i] = false
	end
	return kmultvec(X)
end
#########################################################
"""
```
S(x::Vector{Float64})
S(X::kmultvec)
```
The stereographic embedding of Euclidean space.

Takes a vector of Euclidean space in its stereographic embedding.
"""
function S(x::Vector{Float64})
	l = length(x)
	x2 = vecdot(x,x)
	aux = 2/(x2 + 1)
	X = Vector{Float64}(l+1)
	for i=1:l
		X[i] = copy(aux*x[i])
	end
	X[l+1] = (x2 - 1)/(x2 + 1)
	return euctoga(X)
end
#########################################################
function S(X::kmultvec)
	x = mvectovec([X])
	l = length(x)
	y = Vector{Float64}(l)
	for i=1:l
		y[i] = x[i]
	end
	return S(y)
end
#########################################################
"""
```
Hm(x::Vector{Float64})
```
Auxiliary function.
"""
function Hm(x::Vector{Float64})
	l = length(x)
	y = euctoga(x)
	X = Vector{pbasis}(l+1)
	for i=1:l
		X[i] = kbtopb(y.comp[i])
	end
	m = length(id)
	ide = Vector{Bool}(m+1)
	fill!(ide, false)
	X[l+1] = pbasis(ide, true, 1.0)
	return pmultvec(X)
end
#########################################################
function Hm(S::kmultvec)
	l = length(S.comp)
	X = Vector{pbasis}(l+1)
	for i=1:l
		X[i] = kbtopb(S.comp[i])
	end
	m = length(id)
	ide = Vector{Bool}(m+1)
	fill!(ide, false)
	X[l+1] = pbasis(ide, true, 1.0)
	return pmultvec(X)
end
#########################################################
"""
```
pconformal(x::Vector{Float64})
```
Conformal embedding of a Euclidean vector into a vector of Gp+1,1 with new basis elements being e+ and e-.
"""
function pconformal(x::Vector{Float64})
	x2 = vecdot(x,x)
	a = 0.5*(x2 + 1)
	X = Hm(S(x))
	l = length(X.comp)
	for i=1:l
		X.comp[i].scl = (X.comp[i].scl)*a
	end
	return X
end
function pconformal(x::kmultvec)
	y = mvectovec(x)
	return pconformal(x)
end
#########################################################
"""
```
ipconformal(X::pmultvec)
```
The inverse of "pconformal" function.
"""
function ipconformal(X::pmultvec)
	l=length(X.comp)
	Y = copy(X)
	if l != 0
		m = length(X.comp[1].eb)
	end
	eb = Vector{Bool}(m)
	fill!(eb,false)
	eb2 = copy(eb)
	eb2[m] = true
	ei = pmultvec([pbasis(eb, true, 1.0), pbasis(eb2, false, 1.0)])
	eo = pmultvec([pbasis(eb, true, 0.5), pbasis(eb2, false, -0.5)])
	d = inner(X, ei)
	n = length(d.comp)
	if n != 0
		div = -d.comp[1].scl
	end
	for i=1:l
		Y.comp[i].scl = Y.comp[i].scl/(div)
	end
	B = pblade([Y])
	result = rejection(B, pblade([ei, eo]))
	return result
end
#########################################################
"""
```
conformal(x::Vector{Float64})
```
Conformal embedding of a Euclidean vector into a vector of Gp+1,1 with new basis elements being e∞ and e∘.
"""
function conformal(x::Vector{Float64})
	x2 = vecdot(x,x)
	a = 0.5*(x2 + 1)
	X = Hm(S(x))
	l = length(X.comp)
	for i=1:l
		X.comp[i].scl = (X.comp[i].scl)*a
	end
	return pbtocb(X)
end
#########################################################
"""
```
iconformal(X::cmultvec)
```
The inverse of "conformal" function.
"""
function iconformal(X::cmultvec)
	Y = cbtopb(X)
	Y = ipconformal(Y)
	return pbtocb(Y)
end
#########################################################
"""
```
mvectovec(X::Vector{cmultvec})
```
Auxiliary function.
"""
function mvectovec(X::Vector{cmultvec})
	l = length(X)
	Y = Vector{pmultvec}(l)
	for i=1:l
		Y[i] = cbtopb(X[i])
	end
	return mvectovec(Y)
end
#########################################################
"""
```
cblade(X::Vector{cmultvec})
```
Creates a blade in the geometric algebra \$ \\mathbb{G}_{k+1,1} \$ composed by the 
outer product of the elements of X.

"""
type cblade
	conj::Vector{cmultvec}
	function cblade(X::Vector{cmultvec})
		T = mvectovec(X)
		if det(T*transpose(T)) != 0
			new(X)
		else
			error("L.d (cblade conversion)")
		end
	end
end
#########################################################
function show(io::IO, A::cblade)
	if length(A.conj) != 0
		print(io, "($(A.conj[1]))")
		for i=2:length(A.conj)
			print(io, "∧($(A.conj[i]))")
		end
	else
		print(io, "0")
	end
end
#########################################################
function bltomv(A::cblade)
	l = length(A.conj)
	X = A.conj[1]
	for i=2:l
		X = outer(X, A.conj[i])
	end
	return X
end
#########################################################
function geoprod(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return geoprod(AA,BB)
end
#########################################################
function Base.:∘(A::cblade,B::cblade)
	return geoprod(A,B)
end
function Base.:∘(a::Number,B::cblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
function Base.:*(a::Number,B::cblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function inner(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return inner(AA,BB)
end
#########################################################
function Base.:⋅(A::cblade,B::cblade)
	return inner(A,B)
end
#########################################################
function outer(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return outer(AA,BB)
end
#########################################################
function Base.:^(A::cblade,B::cblade)
	return outer(A,B)
end
#########################################################
function scalar(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return scalar(AA,BB)
end
#########################################################
function Base.:*(A::cblade,B::cblade)
	return scalar(A,B)
end
#########################################################
function mvsum(A::cblade,B::cblade)
	AA = bltomv(A)
	BB = bltomv(B)
	return mvsum(AA,BB)
end
#########################################################
function Base.:+(A::cblade,B::cblade)
	return mvsum(copy(A),copy(B))
end
function Base.:-(A::cblade,B::cblade)
	AA = bltomv(copy(A))
	BB = bltomv(copy(B))
	return AA-BB
end
#########################################################
"""
```
cbltopbl(A::cblade)
```
Auxiliary function.
"""
function cbltopbl(A::cblade)
	l = length(A.conj)
	V = Vector{pmultvec}(l)
	for i=1:l
		V[i] = cbtopb(A.conj[i])
	end
	B = pblade(V)
	return B
end
#########################################################
"""
```
pbltocbl(A::pblade)
```
Auxiliary function.
"""
function pbltocbl(A::pblade)
	l = length(A.conj)
	V = Vector{cmultvec}(l)
	for i=1:l
		V[i] = pbtocb(A.conj[i])
	end
	B = cblade(V)
	return B
end
#########################################################
function mvreverse(A::cblade)
	B = cbltopbl(A)
	B = mvreverse(B)
	return pbltocbl(B)
end
#########################################################
function conjugate(A::cblade)
	B = cbltopbl(A)
	B = conjugate(B)
	return pbltocbl(B)
end
#########################################################
function magnitude(A::cblade)
	X = A * conjugate(A)
	return sqrt(X)
end
#########################################################
function copy(A::cblade)
	l = length(A.conj)
	X = Vector{cmultvec}(l)
	for i=1:l
		X[i] = copy(A.conj[i])
	end
	return cblade(X)
end
#########################################################
function inverse(A::cblade)
	if (A∘A).comp[1].scl != 0
		div = (A * mvreverse(A))
		l = length(A.conj)
		aux = mvreverse(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	else
		div = (A ⋅ conjugate(A)).comp[1].scl
		l = length(A.conj)
		aux = conjugate(A)
		X = copy(aux)
		m = length(aux.conj[1].comp)
		for j=1:m
			X.conj[1].comp[j].scl = (aux.conj[1].comp[j].scl)*(div)^(-1)
		end
	end
	return X
end
#########################################################
function Base.:∘(a::kbasis,B::kmultvec)
	A = kmultvec(a)
	return A ∘ B
end
#########################################################
function Base.:∘(A::kmultvec,b::kbasis)
	B = kmultvec(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::kbasis,B::kmultvec)
	A = kmultvec(a)
	return A ⋅ B
end
#########################################################
function Base.:⋅(A::kmultvec,b::kbasis)
	B = kmultvec(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::kbasis,B::kmultvec)
	A = kmultvec(a)
	return A ^ B
end
#########################################################
function Base.:^(A::kmultvec,b::kbasis)
	B = kmultvec(b)
	return A ^ B
end
#########################################################
function Base.:*(a::kbasis,B::kmultvec)
	A = kmultvec(a)
	return A * B
end
#########################################################
function Base.:*(A::kmultvec,b::kbasis)
	B = kmultvec(b)
	return A * B
end
#########################################################
function Base.:+(a::kbasis,B::kmultvec)
	A = kmultvec(a)
	return A + B
end
#########################################################
function Base.:+(A::kmultvec,b::kbasis)
	B = kmultvec(b)
	return A + B
end
#########################################################
function Base.:-(a::kbasis,B::kmultvec)
	A = kmultvec(a)
	return A - B
end
#########################################################
function Base.:-(A::kmultvec,b::kbasis)
	B = kmultvec(b)
	return A - B
end
#########################################################
function Base.:∘(a::kblade,B::kmultvec)
	A = bltomv(a)
	return A ∘ B
end
#########################################################
function Base.:∘(A::kmultvec,b::kblade)
	B = bltomv(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::kblade,B::kmultvec)
	A = bltomv(a)
	return A ⋅ B
end
#########################################################
function Base.:⋅(A::kmultvec,b::kblade)
	B = bltomv(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::kblade,B::kmultvec)
	A = bltomv(a)
	return A ^ B
end
#########################################################
function Base.:^(A::kmultvec,b::kblade)
	B = bltomv(b)
	return A ^ B
end
#########################################################
function Base.:*(a::kblade,B::kmultvec)
	A = bltomv(a)
	return A * B
end
#########################################################
function Base.:*(A::kmultvec,b::kblade)
	B = bltomv(b)
	return A * B
end
#########################################################
function Base.:+(a::kblade,B::kmultvec)
	A = bltomv(a)
	return A + B
end
#########################################################
function Base.:+(A::kmultvec,b::kblade)
	B = bltomv(b)
	return A + B
end
#########################################################
function Base.:-(a::kblade,B::kmultvec)
	A = bltomv(a)
	return A - B
end
#########################################################
function Base.:-(A::kmultvec,b::kblade)
	B = bltomv(b)
	return A - B
end
#########################################################
function Base.:∘(a::kblade,b::kbasis)
	A = bltomv(a)
	B = kmultvec(b)
	return A ∘ B
end
#########################################################
function Base.:∘(a::kbasis,b::kblade)
	A = kmultvec(a)
	B = bltomv(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::kblade,b::kbasis)
	A = bltomv(a)
	B = kmultvec(b)
	return A ⋅ B
end
#########################################################
function Base.:⋅(a::kbasis,b::kblade)
	A = kmultvec(a)
	B = bltomv(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::kblade,b::kbasis)
	A = bltomv(a)
	B = kmultvec(b)
	return A ^ B
end
#########################################################
function Base.:^(a::kbasis,b::kblade)
	A = kmultvec(a)
	B = bltomv(b)
	return A ^ B
end
#########################################################
function Base.:*(a::kblade,b::kbasis)
	A = bltomv(a)
	B = kmultvec(b)
	return A * B
end
#########################################################
function Base.:*(a::kbasis,b::kblade)
	A = kmultvec(a)
	B = bltomv(b)
	return A * B
end
#########################################################
function Base.:+(a::kblade,b::kbasis)
	A = bltomv(a)
	B = kmultvec(b)
	return A + B
end
#########################################################
function Base.:+(a::kbasis,b::kblade)
	A = kmultvec(a)
	B = bltomv(b)
	return A + B
end
#########################################################
function Base.:-(a::kblade,b::kbasis)
	A = bltomv(a)
	B = kmultvec(b)
	return A - B
end
#########################################################
function Base.:-(a::kbasis,b::kblade)
	A = kmultvec(a)
	B = bltomv(b)
	return A - B
end
#########################################################
function Base.:∘(a::pbasis,B::pmultvec)
	A = pmultvec(a)
	return A ∘ B
end
#########################################################
function Base.:∘(A::pmultvec,b::pbasis)
	B = pmultvec(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::pbasis,B::pmultvec)
	A = pmultvec(a)
	return A ⋅ B
end
#########################################################
function Base.:⋅(A::pmultvec,b::pbasis)
	B = pmultvec(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::pbasis,B::pmultvec)
	A = pmultvec(a)
	return A ^ B
end
#########################################################
function Base.:^(A::pmultvec,b::pbasis)
	B = pmultvec(b)
	return A ^ B
end
#########################################################
function Base.:*(a::pbasis,B::pmultvec)
	A = pmultvec(a)
	return A * B
end
#########################################################
function Base.:*(A::pmultvec,b::pbasis)
	B = pmultvec(b)
	return A * B
end
#########################################################
function Base.:+(a::pbasis,B::pmultvec)
	A = pmultvec(a)
	return A + B
end
#########################################################
function Base.:+(A::pmultvec,b::pbasis)
	B = pmultvec(b)
	return A + B
end
#########################################################
function Base.:-(a::pbasis,B::pmultvec)
	A = pmultvec(a)
	return A - B
end
#########################################################
function Base.:-(A::pmultvec,b::pbasis)
	B = pmultvec(b)
	return A - B
end
#########################################################
function Base.:∘(a::pblade,B::pmultvec)
	A = bltomv(a)
	return A ∘ B
end
#########################################################
function Base.:∘(A::pmultvec,b::pblade)
	B = bltomv(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::pblade,B::pmultvec)
	A = bltomv(a)
	return A ⋅ B
end
#########################################################
function Base.:⋅(A::pmultvec,b::pblade)
	B = bltomv(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::pblade,B::pmultvec)
	A = bltomv(a)
	return A ^ B
end
#########################################################
function Base.:^(A::pmultvec,b::pblade)
	B = bltomv(b)
	return A ^ B
end
#########################################################
function Base.:*(a::pblade,B::pmultvec)
	A = bltomv(a)
	return A * B
end
#########################################################
function Base.:*(A::pmultvec,b::pblade)
	B = bltomv(b)
	return A * B
end
#########################################################
function Base.:+(a::pblade,B::pmultvec)
	A = bltomv(a)
	return A + B
end
#########################################################
function Base.:+(A::pmultvec,b::pblade)
	B = bltomv(b)
	return A + B
end
#########################################################
function Base.:-(a::pblade,B::pmultvec)
	A = bltomv(a)
	return A - B
end
#########################################################
function Base.:-(A::pmultvec,b::pblade)
	B = bltomv(b)
	return A - B
end
#########################################################
function Base.:∘(a::pblade,b::pbasis)
	A = bltomv(a)
	B = pmultvec(b)
	return A ∘ B
end
#########################################################
function Base.:∘(a::pbasis,b::pblade)
	A = pmultvec(a)
	B = bltomv(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::pblade,b::pbasis)
	A = bltomv(a)
	B = pmultvec(b)
	return A ⋅ B
end
#########################################################
function Base.:⋅(a::pbasis,b::pblade)
	A = pmultvec(a)
	B = bltomv(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::pblade,b::pbasis)
	A = bltomv(a)
	B = pmultvec(b)
	return A ^ B
end
#########################################################
function Base.:^(a::pbasis,b::pblade)
	A = pmultvec(a)
	B = bltomv(b)
	return A ^ B
end
#########################################################
function Base.:*(a::pblade,b::pbasis)
	A = bltomv(a)
	B = pmultvec(b)
	return A * B
end
#########################################################
function Base.:*(a::pbasis,b::pblade)
	A = pmultvec(a)
	B = bltomv(b)
	return A * B
end
#########################################################
function Base.:+(a::pblade,b::pbasis)
	A = bltomv(a)
	B = pmultvec(b)
	return A + B
end
#########################################################
function Base.:+(a::pbasis,b::pblade)
	A = pmultvec(a)
	B = bltomv(b)
	return A + B
end
#########################################################
function Base.:-(a::pblade,b::pbasis)
	A = bltomv(a)
	B = pmultvec(b)
	return A - B
end
#########################################################
function Base.:-(a::pbasis,b::pblade)
	A = pmultvec(a)
	B = bltomv(b)
	return A - B
end
#########################################################
function Base.:∘(a::cbasis,B::cmultvec)
	A = cmultvec(a)
	return A ∘ B
end
#########################################################
function Base.:∘(A::cmultvec,b::cbasis)
	B = cmultvec(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::cbasis,B::cmultvec)
	A = cmultvec(a)
	return A ⋅ B
end
#########################################################
function Base.:⋅(A::cmultvec,b::cbasis)
	B = cmultvec(b)
	return A ⋅ B
end
################c#########################################
function Base.:^(a::cbasis,B::cmultvec)
	A = cmultvec(a)
	return A ^ B
end
#########################################################
function Base.:^(A::cmultvec,b::cbasis)
	B = cmultvec(b)
	return A ^ B
end
#########################################################
function Base.:*(a::cbasis,B::cmultvec)
	A = cmultvec(a)
	return A * B
end
#########################################################
function Base.:*(A::cmultvec,b::cbasis)
	B = cmultvec(b)
	return A * B
end
###################################################c#####
function Base.:+(a::cbasis,B::cmultvec)
	A = cmultvec(a)
	return A + B
end
#########################################################
function Base.:+(A::cmultvec,b::cbasis)
	B = cmultvec(b)
	return A + B
end
#########################################################
function Base.:-(a::cbasis,B::cmultvec)
	A = cmultvec(a)
	return A - B
end
#########################################################
function Base.:-(A::cmultvec,b::cbasis)
	B = cmultvec(b)
	return A - B
end
#########################################################
function Base.:∘(a::cblade,B::cmultvec)
	A = bltomv(a)
	return A ∘ B
end
#########################################################
function Base.:∘(A::cmultvec,b::cblade)
	B = bltomv(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::cblade,B::cmultvec)
	A = bltomv(a)
	return A ⋅ B
end
#########################################################
function Base.:⋅(A::cmultvec,b::cblade)
	B = bltomv(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::cblade,B::cmultvec)
	A = bltomv(a)
	return A ^ B
end
#########################################################
function Base.:^(A::cmultvec,b::cblade)
	B = bltomv(b)
	return A ^ B
end
#########################################################
function Base.:*(a::cblade,B::cmultvec)
	A = bltomv(a)
	return A * B
end
#########################################################
function Base.:*(A::cmultvec,b::cblade)
	B = bltomv(b)
	return A * B
end
#########################################################
function Base.:+(a::cblade,B::cmultvec)
	A = bltomv(a)
	return A + B
end
#########################################################
function Base.:+(A::cmultvec,b::cblade)
	B = bltomv(b)
	return A + B
end
#########################################################
function Base.:-(a::cblade,B::cmultvec)
	A = bltomv(a)
	return A - B
end
#########################################################
function Base.:-(A::cmultvec,b::cblade)
	B = bltomv(b)
	return A - B
end
#########################################################
function Base.:∘(a::cblade,b::cbasis)
	A = bltomv(a)
	B = cmultvec(b)
	return A ∘ B
end
#########################################################
function Base.:∘(a::cbasis,b::cblade)
	A = cmultvec(a)
	B = bltomv(b)
	return A ∘ B
end
#########################################################
function Base.:⋅(a::cblade,b::cbasis)
	A = bltomv(a)
	B = cmultvec(b)
	return A ⋅ B
end
#########################################################
function Base.:⋅(a::cbasis,b::cblade)
	A = cmultvec(a)
	B = bltomv(b)
	return A ⋅ B
end
#########################################################
function Base.:^(a::cblade,b::cbasis)
	A = bltomv(a)
	B = cmultvec(b)
	return A ^ B
end
#########################################################
function Base.:^(a::cbasis,b::cblade)
	A = cmultvec(a)
	B = bltomv(b)
	return A ^ B
end
#########################################################
function Base.:*(a::cblade,b::cbasis)
	A = bltomv(a)
	B = cmultvec(b)
	return A * B
end
#########################################################
function Base.:*(a::cbasis,b::cblade)
	A = cmultvec(a)
	B = bltomv(b)
	return A * B
end
#########################################################
function Base.:+(a::cblade,b::cbasis)
	A = bltomv(a)
	B = cmultvec(b)
	return A + B
end
#########################################################
function Base.:+(a::cbasis,b::cblade)
	A = cmultvec(a)
	B = bltomv(b)
	return A + B
end
#########################################################
function Base.:-(a::cblade,b::cbasis)
	A = bltomv(a)
	B = cmultvec(b)
	return A - B
end
#########################################################
function Base.:-(a::cbasis,b::cblade)
	A = cmultvec(a)
	B = bltomv(b)
	return A - B
end
#########################################################
function Base.:*(a::Number,b::Tuple{Vector{Bool},Bool})
	return pb(b[1], b[2],a)
end
#########################################################
function Base.:*(a::Number,b::Tuple{Vector{Bool},Bool, Bool})
	return cb(b[1], b[2], b[3], a)
end
#########################################################
function mvreverse(a::cbasis)
	b = cbtopb(a)
	B = mvreverse(b)
	return pbtocb(B)
end
#########################################################
function mvreverse(a::cmultvec)
	b = cbtopb(a)
	B = mvreverse(b)
	return pbtocb(B)
end
#########################################################
function dual(a::cbasis)
    b = cbtopb(a)
    c = dual(b)
    return pbtocb(c)
end
#########################################################
function dual(X::cmultvec)
    Y = cbtopb(X)
    Z = dual(Y)
    return pbtocb(Z)
end
#########################################################
function dual(A::cblade)
    B = cbltopbl(A)
    C = dual(B)
    return pbltocbl(C)
end
#########################################################
function kblade(a::kbasis)
	l = length(a.e)
	m = grade(a)
	e = Vector{Bool}(l)
	fill!(e,false)
	A = Vector{kmultvec}(l)
	for i=1:l
		if a.e[i] == true
			e[i] = true
			A[i] = copy(kmultvec([kbasis(e,1.0)]))
		else
			A[i] = copy(kmultvec([kbasis(e,0.0)]))
		end
		e[i] = false
	end
	B = Vector{kmultvec}(m)
	aux = 1
	for i=1:m
		for j=l:(-1):(i)
			if A[j].comp[1].scl != 0.0
				B[i] = A[j]
				aux = j
			end
		end
		A[aux] = copy(kmultvec([kbasis(e,0.0)]))
	end
	B[1].comp[1].scl = a.scl
	return kblade(B)
end
#########################################################
function pblade(a::pbasis)
	l = length(a.eb)
	m = grade(a)
	e = Vector{Bool}(l)
	fill!(e,false)
	A = Vector{pmultvec}(l+1)
	for i=1:l
		if a.eb[i] == true
			e[i] = true
			A[i] = copy(pmultvec([pbasis(e,false,1.0)]))
		else
			A[i] = copy(pmultvec([pbasis(e,false,0.0)]))
		end
		e[i] = false
	end
	if a.ep == true
		A[l+1] = copy(pmultvec([pbasis(e,true,1.0)]))
	else
		A[l+1] = copy(pmultvec([pbasis(e,false,0.0)]))
	end
	B = Vector{pmultvec}(m)
	aux = 1
	for i=1:m
		for j=l+1:(-1):(i)
			if A[j].comp[1].scl != 0.0
				B[i] = A[j]
				aux = j
			end
		end
		A[aux] = copy(pmultvec([pbasis(e,false,0.0)]))
	end
	B[1].comp[1].scl = a.scl
	return pblade(B)
end
#########################################################
function projection(a::kbasis, N::kblade)
    A = kblade(a)
    return (A ⋅ inverse(N)) ∘ N
end
#########################################################
function projection(a::kbasis, n::kbasis)
    A = kblade(a)
    N = kblade(n)
    return (A ⋅ inverse(N)) ∘ N
end
#########################################################
function projection(X::kmultvec, N::kblade)
    l = length(X.comp)
    A = Vector{kbasis}(l)
    for i=1:l
        A[i] = projection(X.comp[i], N).comp[1]
    end
    A = kmultvec(A)
    return reduct(A)
end
#########################################################
function projection(X::kmultvec, n::kbasis)
    N = kblade(n)
    return projection(X,N)
end
#########################################################
function projection(a::pbasis, N::pblade)
	A = pblade(a)
    return (A ⋅ inverse(N)) ∘ N
end
#########################################################
function projection(a::pbasis, n::pbasis)
    A = pblade(a)
    N = pblade(n)
    return (A ⋅ inverse(N)) ∘ N
end
#########################################################
function projection(X::pmultvec, N::pblade)
    l = length(X.comp)
    A = Vector{pbasis}(l)
    for i=1:l
	        A[i] = projection(X.comp[i], N).comp[1]
    end
    A = pmultvec(A)
    return reduct(A)
end
#########################################################
function projection(X::pmultvec, n::pbasis)
    N = pblade(n)
    return projection(X,N)
end
#########################################################
function projection(X::pmultvec, n::pmultvec)
    N = pblade(n)
    return projection(X,N)
end
#########################################################
function projection(X::pbasis, n::pmultvec)
    N = pblade([n])
    return projection(X,N)
end
#########################################################
function projection(X::pblade, n::pmultvec)
    N = pblade([n])
    return projection(X,N)
end
#########################################################
function projection(a::cbasis, n::cbasis)
    x = cbtopb(a)
    y = cbtopb(n)
    z = projection(x,y)
    return pbtocb(z)
end
#########################################################
function projection(a::cbasis, N::cblade)
    x = cbtopb(a)
    y = cbltopbl(N)
    z = projection(x,y)
    return pbtocb(z)
end
#########################################################
function projection(X::cmultvec, N::cblade)
    x = cbtopb(X)
    y = cbltopbl(N)
    z = projection(x,y)
    return pbtocb(z)
end
#########################################################
function projection(X::cmultvec, n::cbasis)
    x = cbtopb(X)
    y = cbtopb(n)
    z = projection(x,y)
    return pbtocb(z)
end
#########################################################
function projection(X::cblade, N::cblade)
    A = cbltopbl(X)
    M = cbltopbl(N)
    Z = projection(A,M)
    return pbtocb(Z)
end
#########################################################
function Base.:+(A::pblade,B::pblade)
	return mvsum(copy(A),copy(B))
end
function Base.:-(A::pblade,B::pblade)
	AA = bltomv(copy(A))
	BB = bltomv(copy(B))
	return AA-BB
end
#########################################################
function Base.:∘(a::kbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:∘(A::kmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:∘(A::kblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:∘(a::pbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:∘(A::pmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:∘(A::pblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:∘(a::cbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:∘(A::cmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:∘(A::cblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:⋅(a::kbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:⋅(A::kmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:⋅(A::kblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:⋅(a::pbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:⋅(A::pmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:⋅(A::pblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:⋅(a::cbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:⋅(A::cmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:⋅(A::cblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:⋅(a::Number, b::kbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function Base.:⋅(a::Number,B::kmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:⋅(a::Number,B::kblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:⋅(a::Number, b::pbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function Base.:⋅(a::Number,B::pmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:⋅(a::Number,B::pblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:⋅(a::Number, b::cbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function Base.:⋅(a::Number,B::cmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:⋅(a::Number,B::cblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:*(a::kbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:*(A::kmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:*(A::kblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:*(a::pbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:*(A::pmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:*(A::pblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:*(a::cbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:*(A::cmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:*(A::cblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:^(a::kbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:^(A::kmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:^(A::kblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:^(a::pbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:^(A::pmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:^(A::pblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:^(a::cbasis, b::Number)
	c = copy(a)
	c.scl = b*c.scl
	return c
end
#########################################################
function Base.:^(A::cmultvec,b::Number)
    C = copy(A)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = b*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:^(A::cblade,b::Number)
	C = copy(A)
	if length(C.conj) != 0
		C.conj[1] = b * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:^(a::Number, b::kbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function Base.:^(a::Number,B::kmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:^(a::Number,B::kblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:^(a::Number, b::pbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function Base.:^(a::Number,B::pmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:^(a::Number,B::pblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function Base.:^(a::Number, b::cbasis)
	c = copy(b)
	c.scl = a*c.scl
	return c
end
#########################################################
function Base.:^(a::Number,B::cmultvec)
    C = copy(B)
	l = length(C.comp)
	for i=1:l
		C.comp[i].scl = a*C.comp[i].scl
	end
	return C
end
#########################################################
function Base.:^(a::Number,B::cblade)
	C = copy(B)
	if length(C.conj) != 0
		C.conj[1] = a * C.conj[1]
		return C
	else
		return C
	end
end
#########################################################
function ==(a::kbasis,b::kbasis)
	if a.e == b.e && a.scl == b.scl
		return true
	else
		return false
	end
end
#########################################################
function ==(a::kmultvec,b::kmultvec)
	A = a - b
	if A.comp[1].scl == 0.0
		return true
	else
		return false
	end
end
#########################################################
function ==(a::kblade,b::kblade)
	A = a - b
	if A.comp[1].scl == 0.0
		return true
	else
		return false
	end
end
#########################################################
function ==(a::pbasis,b::pbasis)
	if a.eb == b.eb && a.ep == b.ep && a.scl == b.scl
		return true
	else
		return false
	end
end
#########################################################
function ==(a::pmultvec,b::pmultvec)
	A = a - b
	if A.comp[1].scl == 0.0
		return true
	else
		return false
	end
end
#########################################################
function ==(a::pblade,b::pblade)
	A = a - b
	if A.comp[1].scl == 0.0
		return true
	else
		return false
	end
end
#########################################################
function ==(a::cbasis,b::cbasis)
	if a.er == b.er && a.ei == b.ei && a.eo == a.eo && a.scl == b.scl
		return true
	else
		return false
	end
end
#########################################################
function ==(a::cmultvec,b::cmultvec)
	A = a - b
	if A.comp[1].scl == 0.0
		return true
	else
		return false
	end
end
#########################################################
function ==(a::cblade,b::cblade)
	A = a - b
	if A.comp[1].scl == 0.0
		return true
	else
		return false
	end
end
#########################################################
function pblade(a::pmultvec)
	l = length(a.comp)
	if l != 1
		return error("Invalid!")
	else
		return pblade(a.comp[1])
	end
end


function conftore(A::cmultvec)
	l = length(A.comp)
	f1 = 0
	f2 = 0
	for i=1:l
		if grade(A.comp[i]) != 1
			f1 = 1
		end
		if A.comp[i].ei == true || A.comp[i].eo == true
			f2 = 1
		end 
	end
	if f1 == 0 && f2 == 0
		re = Vector{kbasis}(l)
		for i=1:l
			re[i] = kbasis(A.comp[i].er, A.comp[i].scl)
		end
		m = length(re[1].e)
		re = kmultvec(re)
		a = mvectovec([re])
		b = Vector{Float64}(m)
		for i=1:m
			b[i] = a[i]
		end
		return b
	else
		return error("Grade error")
	end
end


#########################################################
function Base.:/(a::cmultvec, b::cmultvec)
	l = length(b.comp)
	m = length(a.comp)
	if l!=1
		return error("Can not divide by multivectors!")
	else
		c = copy(b.comp[1])
		if grade(c)!=0
			return error("Can not divide by multivectors!")
		else
			k = c.scl
			A = copy(a)
			for i=1:m
				A.comp[i].scl = A.comp[i].scl/k
			end
			return A
		end
	end
end
#########################################################
end




#function arvore(v,niv,pos,lista,n)
#    if niv<=n
#        if pos==0
#            aux=copy(v)
#            aux[niv]=true
#            niv=niv+1
#            lista=push!(lista,aux)
#            lista=arvore(aux,niv,0,lista,n)
#            lista=arvore(aux,niv,1,lista,n)
#
#        else
#            niv=niv+1
#            aux=copy(v)
#            lista=arvore(aux,niv,0,lista,n)
#            lista=arvore(aux,niv,1,lista,n)
#
#        end
#    end
#    return lista
#end
#
#function buildtree(n)
#   v=falses(n)
#   lista=[v]
#   for i in arvore(v,1,0,[v],n)[2:length(arvore(v,1,0,[v],n))]
#        lista=push!(lista,i)
#   end
#   for i in arvore(v,1,1,[v],n)[2:length(arvore(v,1,1,[v],n))]
#        lista=push!(lista,i)
#   end
#   return lista
#end
#""" This function is used to define the space to work
### Example
#```julia-repl
#julia> layout(3)
#```
#generates a G3 space with base 1,e1,e2,e3,e12,e13,e23,e123
#
#"""
#function layout(dim::Int)
#    bn=buildtree(dim)
#    bnew=Array{Array{Bool,1},1}(2^dim)
#    ind=1
#    for v in bn
#        s=find(x->x==true,v)
#        if isempty(s)
#            eval(parse("const id = $v ;"))
#            bnew[ind]=eval(parse("$v"))
#            ind+=1
#        else
#            conc=string(s[1])
#            for k=2:length(s)
#                conc=string(conc,s[k])
#            end
#            eval(parse("const e$(conc) = $v ;"))
#            bnew[ind]=eval(parse("$v"))
#            ind+=1
#        end
#    end
#    eval(parse("const setupbase = $(bnew) ;"))
#Liga.layout(dim)
#end
#