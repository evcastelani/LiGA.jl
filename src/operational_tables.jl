# This file contains essentials functions to define
# tables to geometric, inner and outer products.

"""
```
bbgeometric(a::BasisBlade, b::BasisBlade) (Liga function)
```
This function is used to determine the geometric product
of basis blade. Essentially, this function is used when layout
function is trigged, but it can be used separated too.
## Example
```julia-repl
julia> bbgeometric(e1,e12)
```
returns e2.
"""
function bbgeometric(a::BasisBlade, b::BasisBlade)
		len = GIdLayout[1]+GIdLayout[2] #p+q
		Ea = copy(LogicalBasisBlade[a.index])
		Eb = copy(LogicalBasisBlade[b.index])
		abscalar = 1
		abindex = 1
		ablogical = Vector{Bool}(undef,len)
		cont = 0
		for i=1:len-1
			if Eb[i] == true
				for j=i+1:len
					if Ea[j] == true
						cont = cont + 1
					end
				end
			end
		end
		for i=1:len
			if Ea[i] != Eb[i]
				ablogical[i] = true
			else
					ablogical[i] = false
					if Ea[i] == true && i>GIdLayout[1]
						cont =cont +1
					end
			end
		end
		abscalar = ((-1)^cont)
		ind=1
		while ablogical != LogicalBasisBlade[ind]
			ind = ind+1
		end
		abindex = ind
		ab=BasisBlade(abscalar*abindex)
		return ab
end


"""
```
bbinner(a::BasisBlade,b::BasisBlade) (Liga function)
```
This function is used to determine the inner product
of basis blades. Essentially, this function is used when layout
function is trigged, but it can be used separated too.
## Example
```julia-repl
julia> bbinner(e12,e2)
```
returns the Basis Blade e1.
"""
function bbinner(a::BasisBlade,b::BasisBlade)
		k=grade(a)
		l=grade(b)
		if k != 0 && l !=0
			return gradeprojection(bbgeometric(a,b),abs(k-l))
		else
			return BasisBlade(0)
		end
end

"""
```
bbouter(a::BasisBlade,b::BasisBlade) (Liga function)
```
This function is used to determine the outer product
of basis blades. Essentially, this function is used when layout
function is trigged, but it can be used separated too.
## Example
```julia-repl
julia> bbouter(e12,e3)
```
returns the Basis Blade e123.
"""
function bbouter(a::BasisBlade,b::BasisBlade)
	k=grade(a)
	l=grade(b)
	return gradeprojection(bbgeometric(a,b),k+l)
end

"""
```
OperationTable() (Liga function)
```
OperationTable is function trigged by layout function and
it is used to define the operational table related to geometric,
inner and outer product. As result, constants array are created
to optimize these products to multivectors.
## Example
```julia-repl
julia> A,B,C,D=OperationTable()
```
returns A a table for geometric product, B a table for inner product,
C a table for outer product and D a tensor.
"""
function OperationTable()
	len=2^(GIdLayout[1]+GIdLayout[2])
	OT=Array{Any,2}(undef,len,len)
	IOT=Array{Any,2}(undef,len,len)
	OOT=Array{Any,2}(undef,len,len)
	Tensor=Vector{SparseMatrixCSC}(undef,len)
	for i=1:len
        for j=1:len
			OT[i,j]=bbgeometric(BasisBlade(i),BasisBlade(j))
			IOT[i,j]=bbinner(BasisBlade(i),BasisBlade(j))
			OOT[i,j]=bbouter(BasisBlade(i),BasisBlade(j))
		end
	end
	for k=1:len
		Tensor[k]=spzeros(len,len)
		for i=1:len
			for j=1:len
				if OT[i,j]==BasisBlade(k)
					Tensor[k][i,j]=1.0
				end
				if OT[i,j]==-BasisBlade(k)
				   Tensor[k][i,j]=-1.0
			    end
			end
		end
	end
	#display(Tensor)
	return OT,IOT,OOT,Tensor
end
