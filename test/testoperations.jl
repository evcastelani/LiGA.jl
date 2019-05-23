@testset "Testing operations and Functions" begin

    layout(3,1,"Conformal")
    u=1.0*id+2.0*e1-e0
    v=2.0*id+1.0*e1+e∞
    # sum
    @test u+v == 3.0*id + 3.0*e1 + 1.5*ep + 0.5*en

    # difference
    @test u-v ==-1.0*id + 1.0*e1 - 0.5*ep - 1.5*en

    # scalar
    @test 1.5*v ==3.0*id + 1.5*e1 + 1.5*ep + 1.5*en
    @test scalar(1.5,v) == 3.0*id + 1.5*e1 + 1.5*ep + 1.5*en

    # reverse
    @test reverse(v) == 2.0*id + 1.0*e1 + 1.0*ep + 1.0*en

    # geometric
    @test u∘v == 5.0*id + 5.0*e1 + 2.0*ep + 1.5*e1p + 2.5*e1n + 1.0*epn    
    @test geometric(u,v) == 5.0*id + 5.0*e1 + 2.0*ep + 1.5*e1p + 2.5*e1n + 1.0*epn    

    # inner
    @test u⋅v == 3.0*id 
    @test inner(u,v) == 3.0*id

    # outer 
    @test u^v == 2.0*id + 5.0*e1 + 2.0*ep + 1.5*e1p + 2.5*e1n + 1.0*epn
    @test outer(u,v) == 2.0*id + 5.0*e1 + 2.0*ep + 1.5*e1p + 2.5*e1n + 1.0*epn

    # gradeprojection
    @test gradeprojection(u,1) == 2.0*e1 + 0.5*ep - 0.5*en
    @test gradeprojection(v,1) == 1.0*e1 + 1.0*ep + 1.0*en

    # reverse
    @test reverse(u) == 1.0*id + 2.0*e1 + 0.5*ep - 0.5*en
    @test reverse(v) == 2.0*id + 1.0*e1 + 1.0*ep + 1.0*en

    # dual 
    @test dual(u) ==  0.5*e12p - 0.5*e12n - 2.0*e2pn - 1.0*e12pn 
    @test dual(v) == - 1.0*e12p - 1.0*e12n - 1.0*e2pn - 2.0*e12pn

    # inverse
    @test inverse(u)∘u ≈ 1.0*id
    @test inverse(v)∘v ≈ 1.0*id

end