function arvore(v,niv,pos,lista,n)
    if niv<=n
        if pos==0
            aux=copy(v)
            aux[niv]=true
            niv=niv+1
            lista=push!(lista,aux)
            lista=arvore(aux,niv,0,lista,n)
            lista=arvore(aux,niv,1,lista,n)

        else
            niv=niv+1
            aux=copy(v)
            lista=arvore(aux,niv,0,lista,n)
            lista=arvore(aux,niv,1,lista,n)

        end
    end
    return lista
end

function buildtree(n)
   v=falses(n)
   lista=[v]
   for i in arvore(v,1,0,[v],n)[2:length(arvore(v,1,0,[v],n))]
        lista=push!(lista,i)
   end
   for i in arvore(v,1,1,[v],n)[2:length(arvore(v,1,1,[v],n))]
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
bn=buildtree(dim)
bnew=Array{Array{Bool,1},1}(2^dim)
ind=1
for v in bn
    s=find(x->x==true,v)
    if isempty(s)
        eval(parse("const id = $v"))
        bnew[ind]=eval(parse("$v"))
        ind+=1
    else
        conc=string(s[1])
        for k=2:length(s)
            conc=string(conc,s[k])
        end
        eval(parse("const e$(conc) = $v"))
        bnew[ind]=eval(parse("$v"))
        ind+=1
    end
end
eval(parse("const setupbase = $(bnew)"))
#Liga.layout(dim)
end
