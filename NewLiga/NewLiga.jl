# This version is focused in Julia 1.0 version
module NewLiga
    # Dependencies
    import Base.show
    import Base.reverse
    import Base.convert
    using SparseArrays,LinearAlgebra
    # New definitions
    export  BasisBlade,MultiVector,
            GAAbstractType,LayoutInfo,
            layout,layout_info,
            tree,buildtree,
            grade,grademinus,
            gradeplus,gradeprojection,
            bbgeometric,bbinner,bboute,
            geometric,scalar,inner,
            outer,mvsum,mvdiff,I,
            reverse,conjugate,⋅,
            dual,inverse,Γ,
            embedding

    # Main files

    # Objects.jl produces two objects: BasisBlade and MultiVectors
    include("Objects.jl")

    # Layout.jl produces the enviroment to work. In current version
    # two enviroment are suported: GA and CGA (Geometric Algebra and
    # Conformal Geometric Algebra)
    include("Layout.jl")

    # OperationalTables.jl produces arrays of BasisBlade geometric,
    # inner and outer products. This is an important function because
    # this is the optimization of engine of Liga
    include("OperationalTables.jl")

    # LigaFunctions.jl contains all implemented functions like MultiVector
    # geometric product.
    include("LigaFunctions.jl")

    # LigaShow.jl contains some modifications to print BasisBlade and
    # MultiVectors.
    include("LigaShow.jl")
end
