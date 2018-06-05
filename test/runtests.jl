using Liga
using Base.Test

# write your own tests here
layout(3)

@testset "layout_tests" begin
    @test id == [false,false,false]
    @test e1 == [true,false,false]
    @test e2 == [false,true,false]
    @test e123 == [true,true,true]
end

@testset "geoprod_tests" begin
    @test geoprod(kb(e1,1.0),kb(e12,2.0)) == kb(e2,2.0)
    v=3.0*e23-1.0*e13+3.0*e1+1.0*e2
    @test geoprod(3.0*e1+1.0*e2,1.0*e123+1.0*id) == v
end

