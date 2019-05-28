#layout
#auxiliary function to binary association
""" This function is used to define a list of logical elements 
for use in order to define an algebra of layout function. It is 
an internal function.
"""
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

"""
This function compose the list given in tree function. This is 
an internal function.
"""
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

"""
This is an internal function used to print in REPL an UTF8 char 
when the layout pretty option is trigged.
"""
function prettylayout(p,q)
	D=['₁','₂','₃','₄','₅','₆','₇','₈','₉']
	if 1<=p<10 && 1<=q<10
		return [D[p],D[q]]
	else
		if q==0
			return [D[p]]
		end
	end
end

""" This function is used to define the algebra  to work. If you
want use Liga, is mandatory to use this function. 
## Example
```julia-repl
julia> layout(3,0,"GA")
```
generates a G3 algebra  with base 1,e1,e2,e3,e12,e13,e23,e123.
If you want to redefine the space, is highly recommended to restart 
Julia. 
"""
function layout(p::Int,q::Int,algebra::String; custom_repl="enable")
	# changing mode prompt
	
	eval(:( GIdLayout=($p,$q,$(algebra),$(2^(p+q))))) #Global Indentify Layout Parameters
	eval(:( LigaPrecision=10.0^(-10)))
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
				#abngr=bngr[j]
				bnval[j]=bnval[i]
				bn[j]=bn[i]
				#bngr[j]=bngr[i]
				bnval[i]=abnval
				bn[i]=abn
				#bngr[i]=abngr
			end
		end
	end

	if algebra  == "GA"
		if custom_repl=="enable"
			#OhMyREPL.colorscheme!("JuliaDefault")
			OhMyREPL.input_prompt!("𝔾$(String(prettylayout(p,q))) >", :green)
			#OhMyREPL.output_prompt!("𝔾$(String(prettylayout(p,q))) >", :yellow)
			#Base.active_repl.interface.modes[1].prompt=" 𝔾$(String(prettylayout(p,q))) > "
		else
			Base.active_repl.interface.modes[1].prompt="julia> "
		end
		bnew=Array{Vector{Bool},1}(undef,2^dim)
		ind=1
		for v in bn
			s=findall(x->x==true,v)
			if isempty(s)
				eval(Meta.parse(" id = BasisBlade(1) ;"));
				eval(Meta.parse("export id;"));
				bnew[ind]=v
				ind+=1
			else
				conc=string(s[1])
				for k=2:length(s)
					conc=string(conc,s[k])
				end
				eval(Meta.parse(" e$(conc) = BasisBlade($(ind)) ;"));
				eval(Meta.parse("export e$(conc);"));
				bnew[ind]=v
				ind+=1
			end
		end
		eval(Meta.parse("export LogicalBasisBlade  ;"));
		eval(Meta.parse(" LogicalBasisBlade = $(bnew);"));
		#Run operation table
		eval(Meta.parse("export BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ"));
		eval(Meta.parse("BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ=OperationTable()"));

	else
		if algebra == "Projective"
			if q != 0
				error("Something wrong with the second argument of layout(arg1,arg2,arg3), probably you need to setup equal 0.")
			else
				if custom_repl=="enable"
					#OhMyREPL.colorscheme!("JuliaDefault")
					OhMyREPL.input_prompt!("𝔾$(String(prettylayout(p,q))) [Proj]>", :green)
					#OhMyREPL.output_prompt!("𝔾$(String(prettylayout(p,q))) [Proj]>", :yellow)
					#Base.active_repl.interface.modes[1].prompt=" 𝔾$(String(prettylayout(p,q))) > "
				else
					Base.active_repl.interface.modes[1].prompt="julia> "
				end
				#printstyled("\n An apropriated embedding for Projective Space was created \n",color=:blue)
				bnew=Array{Vector{Bool},1}(undef,2^dim)
				ind=1
				for v in bn
					s=findall(x->x==true,v)
					if isempty(s)
						eval(Meta.parse(" id = BasisBlade(1) ;"));
						eval(Meta.parse("export id;"));
						bnew[ind]=v
						ind+=1
					else
						conc=string(s[1])
						for k=2:length(s)
							conc=string(conc,s[k])
						end
						eval(Meta.parse(" e$(conc) = BasisBlade($(ind)) ;"));
						eval(Meta.parse("export e$(conc);"));
						bnew[ind]=v
						ind+=1
					end
				end
				eval(Meta.parse("export LogicalBasisBlade  ;"));
				eval(Meta.parse(" LogicalBasisBlade = $(bnew);"));
				#Run operation table
				eval(Meta.parse("export BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ"));
				eval(Meta.parse("BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ=OperationTable()"));
			end
		end
		if algebra =="Conformal"
			if q != 1
				error("Something wrong with the second argument of layout(arg1,arg2,arg3), probably you need to setup equal 1.")
			else
				if custom_repl=="enable"
					OhMyREPL.input_prompt!("𝔾$(String(prettylayout(p,q)))[Conf]>", :green)
			#		OhMyREPL.output_prompt!("𝔾$(String(prettylayout(p,q)))[Conf]>", :yellow)
					#Base.active_repl.interface.modes[1].prompt=" 𝔾$(String(prettylayout(p,q))) > "
				else
					Base.active_repl.interface.modes[1].prompt="julia> "
				end
				#printstyled("\n An apropriated embedding for Conformal Space was created \n",color=:blue)
				#auxmask=zeros(2^dim)
				bnew=Array{Vector{Bool},1}(undef,2^dim)
				ind=1
				for v in bn
					s=findall(x->x==true,v)
					if isempty(s)
						eval(Meta.parse(" id = BasisBlade(1) ;"));
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
								#auxmask[ind]+=1
							end
							if s[k]==dim
								conc=string(conc,"n")
								#auxmask[ind]+=2
							end
							if s[k]<dim-1
								conc=string(conc,s[k])
							end
						end
						eval(Meta.parse(" e$(conc) = BasisBlade($(ind)) ;"));
						eval(Meta.parse("export e$(conc);"));
						bnew[ind]=v
						ind+=1
					end
				end
				#eval(Meta.parse("const vmask=$(auxmask)"))
				eval(Meta.parse("export LogicalBasisBlade  ;"));
				eval(Meta.parse(" LogicalBasisBlade = $(bnew);"));
				#Run operation table
				eval(Meta.parse("export BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ"));
				eval(Meta.parse("BBgeoprodTable,BBinnerprodTable,BBouterprodTable,Γ=OperationTable()"));
				eval(Meta.parse(" e∞ = 1.0*en+1.0*ep ;"));
				eval(Meta.parse("export e∞ "));
				eval(Meta.parse(" e0 = 0.5*en-0.5*ep ;"));
				eval(Meta.parse("export e0 "));
			end
		end
		if algebra!="Projective" && algebra!="Conformal"
			error("Not indentified algebra- Necessary redefine the layout!")
		end
	end

	#eval(:(const Layout=LayoutInfo([$p,$q],$(2^(p+q)),$(algebra ))))
	eval(Meta.parse(" I=BasisBlade(2^($(dim)))"));
	return ;
end

"""
This function show basic information after run layout function. For example if you 
want do see basis elements and quantities related to created space is enough to follow
## Example
```julia-repl
julia> layoutinfo()
```
"""
function layoutinfo()
	printstyled(" ( p, q ) signature  = ( $(GIdLayout[1]),$(GIdLayout[2]) )\n",color=:blue)
	printstyled(" algebra   = $(GIdLayout[3])\n",color=:blue)
	printstyled(" dimension = $(GIdLayout[4])\n",color=:blue)
	printstyled(" identification of basis blade elements \n",color=:blue)
	for i=1:2^(GIdLayout[1]+GIdLayout[2])
		printstyled(" $(i) -> $(BasisBlade(i))\n",color=:blue )
	end
end
