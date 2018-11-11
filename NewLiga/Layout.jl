#layout
#auxiliary function to binary association
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

""" This function is used to define the algebra  to work
## Example
```julia-repl
julia> layout(3,0,"GA")
```
generates a G3 algebra  with base 1,e1,e2,e3,e12,e13,e23,e123

"""
function layout(p::Int,q::Int,algebra::String)
	eval(:(const GIdLayout=($p,$q,$(algebra),$(2^(p+q))))) #Global Indentify Layout Parameters
	eval(:(const LigaPrecision=10.0^(-10)))
	dim=p+q
	bn=buildtree(dim);
    #abn=BitArray{1}(dim)
    abn=BitArray(undef,dim,1)
	abnval=0.0
	#sorting bn
	bngr=zeros(2^dim)
	abngr=0
	#sorting bn
	bnval=zeros(2^dim)
	for i=1:2^dim
		bnval[i]=0.0
		bngr[i]=0
		for j=1:length(bn[i])
			if bn[i][j]==true
				bnval[i]+=2^(j-1)
				bngr[i]+=1
			end
		end
	end
	#sorting by grade
	for i=1:2^dim #sorting by grade
		for j=i+1:2^dim
			if bngr[i]>bngr[j]
				abnval=bnval[j]
				abn=bn[j]
				abngr=bngr[j]
				bnval[j]=bnval[i]
				bn[j]=bn[i]
				bngr[j]=bngr[i]
				bnval[i]=abnval
				bn[i]=abn
				bngr[i]=abngr
			end
		end
	end
	#sorting by value
	for i=1:2^dim
		for j=i+1:2^dim
			if bnval[i]>bnval[j] && bngr[i]==bngr[j]
				abnval=bnval[j]
				abn=bn[j]
			#	abngr=bngr[j]
				bnval[j]=bnval[i]
				bn[j]=bn[i]
			#	bngr[j]=bngr[i]
				bnval[i]=abnval
				bn[i]=abn
			#	bngr[i]=abngr
			end
		end
	end

	if algebra  == "GA"
    	bnew=Array{Vector{Bool},1}(undef,2^dim)
    	ind=1
    	for v in bn
    	    s=findall(x->x==true,v)
    	    if isempty(s)
				eval(Meta.parse("const id = BasisBlade(1) ;"));
				eval(Meta.parse("export id;"));
    	        bnew[ind]=v
    	        ind+=1
    	    else
    	        conc=string(s[1])
    	        for k=2:length(s)
    	            conc=string(conc,s[k])
    	        end
				eval(Meta.parse("const e$(conc) = BasisBlade($(ind)) ;"));
				eval(Meta.parse("export e$(conc);"));
    	        bnew[ind]=v
    	        ind+=1
    	    end
		end
		eval(Meta.parse("export LogicalBasisBlade  ;"));
		eval(Meta.parse("const LogicalBasisBlade = $(bnew);"));
		#Run operation table
		eval(:(const BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ=OperationTable()))

	else
		if algebra == "Projective"
			if q != 0
				error("Something wrong with the second argument of layout(arg1,arg2,arg3), probably you need to setup equal 0.")
            else
				printstyled("\n An apropriated embedding for Projective Space was created \n",color=:blue)
				bnew=Array{Vector{Bool},1}(undef,2^dim)
    			ind=1
    			for v in bn
    	    		s=findall(x->x==true,v)
    	    		if isempty(s)
						eval(Meta.parse("const id = BasisBlade(1) ;"));
						eval(Meta.parse("export id;"));
    	        		bnew[ind]=v
    	        		ind+=1
    	    		else
    	        		conc=string(s[1])
    	        		for k=2:length(s)
    	            		conc=string(conc,s[k])
    	        		end
						eval(Meta.parse("const e$(conc) = BasisBlade($(ind)) ;"));
						eval(Meta.parse("export e$(conc);"));
    	        		bnew[ind]=v
    	        		ind+=1
    	    		end
				end
				eval(Meta.parse("export LogicalBasisBlade  ;"));
				eval(Meta.parse("const LogicalBasisBlade = $(bnew);"));
				#Run operation table
				eval(:(const BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ=OperationTable()))
			end
		end
		if algebra =="Conformal"
			if q != 1
				error("Something wrong with the second argument of layout(arg1,arg2,arg3), probably you need to setup equal 1.")
            else
				printstyled("\n An apropriated embedding for Conformal Space was created \n",color=:blue)
				#auxmask=zeros(2^dim)
    	        bnew=Array{Vector{Bool},1}(undef,2^dim)
    	        ind=1
    	        for v in bn
    	            s=findall(x->x==true,v)
    	            if isempty(s)
		    		    eval(Meta.parse("const id = BasisBlade(1) ;"));
		    		    eval(Meta.parse("export id;"));
    	                bnew[ind]=v
    	                ind+=1
    	            else
                        if s[1]==dim-1
                            conc="p"
							#auxmask[ind]+=1
                        end
                        if s[1]==dim
                            conc="n"
							#auxmask[ind]+=2
                        end
                        if s[1]<dim-1
                            conc=string(s[1])
                        end
                        for k=2:length(s)
                            if s[k]==dim-1
                                conc=string(conc,"p")
							#	auxmask[ind]+=1
                            end
                            if s[k]==dim
                                conc=string(conc,"n")
							#	auxmask[ind]+=2
                            end
                            if s[k]<dim-1
                                conc=string(conc,s[k])
                            end
    	                end
		        		eval(Meta.parse("const e$(conc) = BasisBlade($(ind)) ;"));
		        		eval(Meta.parse("export e$(conc);"));
    	                bnew[ind]=v
    	                ind+=1
    	            end
		        end

				#eval(Meta.parse("const vmask=$(auxmask)"))
		        eval(Meta.parse("export LogicalBasisBlade  ;"));
		        eval(Meta.parse("const LogicalBasisBlade = $(bnew);"));
				#Run operation table

				eval(:(const BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ=OperationTable()))
				eval(Meta.parse("const e∞ = 1.0*en+1.0*ep ;"));
				eval(Meta.parse("export e∞ "));
				eval(Meta.parse("const e0 = 0.5*en-0.5*ep ;"));
				eval(Meta.parse("export e0 "));
			end

		end
		if algebra!="Projective" && algebra!="Conformal"
			error("Not indentified algebra- Necessary redefine the layout!")
		end
	end

	#eval(:(const Layout=LayoutInfo([$p,$q],$(2^(p+q)),$(algebra ))))
	eval(Meta.parse("const I=BasisBlade(2^($(dim)))"));
	return ;
end

function layout_info()
	printstyled(" ( p, q ) signature  = ( $(GIdLayout[1]),$(GIdLayout[2]) )\n",color=:blue)
	printstyled(" algebra   = $(GIdLayout[3])\n",color=:blue)
	printstyled(" dimension = $(GIdLayout[4])\n",color=:blue)
	printstyled(" identification of basis blade elements \n",color=:blue)
	for i=1:2^(GIdLayout[1]+GIdLayout[2])
		printstyled(" $(i) -> $(BasisBlade(i))\n",color=:blue )
	end
end
