@testset "Testing objects and Layout" begin
    layout(3,1,"GA")
    obj1=BasisBlade(2)
    answer1=e1
    @test obj1==answer1 
    obj2=[BasisBlade(1),BasisBlade(2),BasisBlade(3),BasisBlade(4),BasisBlade(5)] 
    answer2=[id,e1,e2,e3,e4]
    @test obj2==answer2
    obj3=BasisBlade(16)
    answer3=e1234
    @test obj3==answer3

    layout(5,0,"Projective")
    obj1=BasisBlade(2)
    answer1=e1
    @test obj1==answer1 
    obj2=[BasisBlade(1),BasisBlade(2),BasisBlade(3),BasisBlade(4),BasisBlade(5),BasisBlade(6) ]
    answer2=[id,e1,e2,e3,e4,e5]
    @test obj2==answer2
    obj3=BasisBlade(32)
    answer3=e12345
    @test obj3==answer3

    layout(3,1,"Conformal")
    obj1=BasisBlade(2)
    answer1=e1
    @test obj1==answer1 
    obj2=[BasisBlade(1),BasisBlade(2),BasisBlade(3),BasisBlade(4),BasisBlade(5)] 
    answer2=[id,e1,e2,ep,en]
    @test obj2==answer2
    obj3=BasisBlade(16)
    answer3=e12pn
    @test obj3==answer3
    obj4=e0
    obj5=e∞
    @test obj4==-0.5*ep+0.5*en
    @test obj5==1.0*ep+1.0*en

end