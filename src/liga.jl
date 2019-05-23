# This version is focused in Julia 1.0 version
module Liga
	# Dependencies
	import Base.show
	import Base.reverse
	import Base.convert
	#   import OhMyREPL
	using SparseArrays,LinearAlgebra
	using OhMyREPL
	# New definitions
	export  BasisBlade,MultiVector,
	    GAAbstractType,Blade,
	    layout,layout_info,
	    tree,buildtree,
	    grade,grademinus,
	    gradeplus,gradeprojection,
	    bbgeometric,bbinner,bbouter,
	    geometric,scalar,inner,
	    outer,mvsum,mvdiff,
	    reverse,conjugate,⋅,
	    dual,inverse,OperationTable,
	    embedding,copy,isequal,==,
	    multivector,≈,magnitude,blade

	# Main files
	# Objects.jl produces two objects: BasisBlade and MultiVectors
	include("liga_objects.jl")
	# Layout.jl produces the enviroment to work. In current version
	# three enviroments are suported: GA, Conformal and Projetive
	include("layout.jl")
	# OperationalTables.jl produces arrays of BasisBlade geometric,
	# inner and outer products. This is an important function because
	# this is the optimization of engine of Liga
	include("operational_tables.jl")
	# LigaFunctions.jl contains all implemented functions like MultiVector
	# geometric product.
	include("liga_functions.jl")
	# LigaShow.jl contains some modifications to print BasisBlade and
	# MultiVectors.
	include("liga_show.jl")
end
