#printIO
function show(io::IO, b::BasisBlade)
	if GIdLayout[3]=="GA" || GIdLayout[3]=="Projective"
		if b.index>0
			if b.index == 1
				print(io, "id")
			else
				print(io, "e")
				for i=1:GIdLayout[1]+GIdLayout[2]
					if LogicalBasisBlade[b.index][i]==true
						print(io, "$i" )
					end
				end
			end
		else
			if b.index == -1
				print(io, "-id")
			else
				if b.index==0
					print(io,"0")
				else
					print(io, "-e")
					for i=1:GIdLayout[1]+GIdLayout[2]
						if LogicalBasisBlade[-b.index][i]==true
							print(io, "$i" )
						end
					end
				end
			end
		end
	end
	if GIdLayout[3]=="Conformal"
		if b.index>0
			if b.index == 1
				print(io, "id")
			else
				print(io, "e")
				for i=1:GIdLayout[1]+GIdLayout[2]
					if LogicalBasisBlade[b.index][i]==true
						if i<GIdLayout[1]+GIdLayout[2]-1
							print(io, "$i" )
						else
							if i==GIdLayout[1]+GIdLayout[2]-1
								print(io,"+")
							else
								print(io,"-")
							end
						end
					end
				end

			end
		else
			if b.index == -1
				print(io, "-id")
			else
				if b.index==0
					print(io,"0")
				else
					print(io, "-e")
					for i=1:GIdLayout[1]+GIdLayout[2]
						if LogicalBasisBlade[-b.index][i]==true
							if i<GIdLayout[1]+GIdLayout[2]-1
								print(io, "$i" )
							else
								if i==GIdLayout[1]+GIdLayout[2]-1
									print(io,"+")
								else
									print(io,"-")
								end
							end
						end
					end
				end
			end
		end

	end
end
function show(io::IO,M::MultiVector)
	#len=length(M.comp)
	#if maximum(abs.(M.comp))<=LigaPrecision
	#	print(io,"0.0")
	#else
	#	if M.comp[1]!=0.0
	#		if M.comp[1]>0.0
	#			print(io,"( $(M.comp[1])$(BasisBlade(1)) )")
	#		else
	#	end
	#	for i=2:len
	#		if M.comp[i]!=0.0
	#			if M.comp[i]>0.0
	#				print(io," + ( $(M.comp[i])$(BasisBlade(i)) )")
	#			else
	#				print(io," - ( $(abs(M.comp[i]))$(BasisBlade(i)) )")
	#			end
	#		end
	#	end
	#end

		ind,val=findnz(M.comp)
		if maximum(abs.(M.comp))<=LigaPrecision
			print(io,"0.0")
		else
			k=true
			for i in ind
				if abs(M.comp[i])>LigaPrecision
					if k==true
						if M.comp[i]>0.0
							print(io,"( $(round(M.comp[i],digits=4)) $(BasisBlade(i)) )")
						else
							print(io," - ( $(round(abs(M.comp[i]),digits=4)) $(BasisBlade(i)) )")
						end
						k=false
					else
						if M.comp[i]>0.0
							print(io," + ( $(round(M.comp[i],digits=4)) $(BasisBlade(i)) )")
						else
							print(io," - ( $(round(abs(M.comp[i]),digits=4)) $(BasisBlade(i)) )")
						end
					end
				end
			end
		end
end
