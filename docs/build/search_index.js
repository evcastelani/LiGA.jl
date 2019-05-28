var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": "This page is devoted to document and make easier the use of Liga- Library for  Geometric Algebra. Our intent is provide a good experience of use for students  and reseachers interesting in use (and programming) a Geometric Algebra  framework. All scripts were made in Julia language,  version 1.0. There is a miscellany of libraries for GA study. Our intent here is not to provide software to compete with such libraries in performance or generality. Our goal is to provide a relatively simple, open, easy to install and focused environment for Julia users. This library was built  just reading basic operations in books, that is, no worries about analyzing the best way to program all functions based on other libraries.We highlight that the version 0.5 of Liga is a considerable faster than previous versions with more optimized calculus. Of course, this work is all the time in progress  and several things are changing in order to get better results. It\'s important a previous contact with Geometric Algebra in theoretical aspects. In this sense, we recommend the following references:[1] Christian Perwass, Geometric Algebra with Applications in Engineering,  Springer Series in Geometry and Computing, 2009.[2] L.Dorst, D.Fontijne and S.Mann, Geometric Algebra for Computer Science: an object oriented approach to geometry, Morgan Kaufman publishers, 2007.We really like the book given by [1], so we tried to provide an experience focused on that notation."
},

{
    "location": "index.html#Current-Status-1",
    "page": "Overview",
    "title": "Current Status",
    "category": "section",
    "text": "The current status of this project is beta quality, don\'t use for anything important. We can provide some multi-dimensional spaces but without visualization of elements (it is a goal feature for future)."
},

{
    "location": "index.html#Developed-by-1",
    "page": "Overview",
    "title": "Developed by",
    "category": "section",
    "text": "This project was developed as part of master degree dissertation of Vinicius Riter. The main contributors of this work areVinicius Riter\nEmerson Vitor Castelani\nJair da Silva\nWesley Shirabayashi\nValter Soares Camargo"
},

{
    "location": "gettingstarted.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "gettingstarted.html#Installation-1",
    "page": "Getting Started",
    "title": "Installation",
    "category": "section",
    "text": "This package is supported just for Julia version 1.0. Consequently,  it uses package 3.0. Currently Liga is in Metadata.jl, so the package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:pkg> add Liga"
},

{
    "location": "gettingstarted.html#Layouts-1",
    "page": "Getting Started",
    "title": "Layouts",
    "category": "section",
    "text": "To start to use Liga, just type:push!(LOAD_PATH,\"/Users/emersonvitor/Google Drive/NewLiga/NewLiga_16\")\nusing NewLigaThe second step to use Liga is define an environment. The environment is defined through the layout function.layout(3,1,\"GA\")In this case, were created the space mathbbG_31. If you want to see more information about the created space, just typelayout_info()Essentially, the layout function has 4 arguments (one is an optional kwarg): layout(p::Int,q::Int,algebra::String; custom_repl=\"disable\"). The p,q, algebra arguments are mandatory and they are used to define the space with signature p,q. For example, if you want to create the conformalized space mathbbG_41, you need to type:layout(4,1,\"Conformal\")wich is equivalent to layout(4,1,\"GA\")but with differents notations, see the layout_info() for each layout tested. The custom_repl keyword argument is optional. By default this is disable and changes the input prompt format (we are using OhMyREPL package here) in order to identify the space. If you want enable this option just type  something likelayout(4,1,\"GA\",custom_repl=\"enable\")"
},

{
    "location": "gettingstarted.html#Basis-Blade-1",
    "page": "Getting Started",
    "title": "Basis Blade",
    "category": "section",
    "text": "All types introduced in Liga are subtype of  GAAbstractType. BasisBlade type is the most primitive subtype wich is used to define more complete structures, like MultiVectors. Once we define a layout, a set of BasisBlade is created. For example, let us create the space mathbbG_3.layout(3,0,\"GA\")Now, the elements id,e1,e2,e3,e12,e13,e23,e123 were created and they are associated to integer. typeof(e1)\ne1.indexIn addition, operation tables related to the geometric, internal and outer products are created too. In order to acess these tables you just need to typeBBgeoprodTable\nBBinnerprodTable\nBBouterprodTableThe id element is a BasisBlade used in place of the number 1. By default, if (for example)e12 exists does not exist e21. This does not mean that it can not be set. I you typee21=-e12\ntypeof(e21)\ne21.indexWith this command we create a new BasisBlade element but, in general, it isn\'t a necessary definition. BasisBlade elements are immutable types so, alternatively, you can create a new one by bb=BasisBlade(-5)\nbb.index"
},

{
    "location": "gettingstarted.html#MultiVectors-1",
    "page": "Getting Started",
    "title": "MultiVectors",
    "category": "section",
    "text": "MultiVectors are mutable structures associated to sparse vectors. To ilustrate, consider the following examples:using SparseArrays\n\nMultiVector(sparsevec([2,4,6],[1.5,-2.0,3.0]))\n\nMultiVector(sparsevec([1],[1.5]))\n\nMultiVector(sparsevec([1,8],[-2.0,3.0]))Of course, this is not a very didactic way of defining a multivector.  A simpler way to create them is through BasisBlade. A better (and recommended) way is ilustrated below. v=1.5*e1 - 2.0*e3 + 3.0*e13\n\n1.5*id\n\n2.0*id+3.0*e123\nNote that there when is defined a scalar product or sum or diference between an real scalar and BasisBlade element a new type (MultiVector) is returned. A MultiVector struct is showed just in function and order of the Basis Blade set that define the space.  "
},

{
    "location": "gettingstarted.html#Blades-1",
    "page": "Getting Started",
    "title": "Blades",
    "category": "section",
    "text": "Blades is an important structure in context of GA. In Liga, a Blade is a subtype of GAAbstract type given by a set of LI vectors. For define a Blade element you need a set of vectors, the number of vectors in this set and a bool variable. There are some functions in Liga that are devotated to work just with Blade. In other cases, a Blade is converted to MultiVector type and all functions can be used too.  Blade([[1.0,2.0],[1.0,0.0]],2,true)This way is the most primitive way to define a Blade. If you need to check if the set of vectors is a LI set, you need to setup the last parameter like false. In order to understand the properties of a Blade, see the following example.v = Blade([[1.0,2.0],[1.0,0.0]],2,true)\nv.vectors\nv.grade\nv.checkedIf a Blade is defined with false value in checked parameter, the Blade object will check and if the vector set is LI a Blade with true value is returned.   v = Blade([[1.0,2.0],[1.0,0.0]],2,false)\nv.vectors\nv.grade\nv.checkedIn the case of LD vector set an error message is showed.w = Blade([[1.0,2.0],[1.0,2.0]],2,false)The main reason for the introduction of this Boolean variable (checked) is to avoid, in situations that require further processing, unnecessary checks."
},

{
    "location": "functions.html#",
    "page": "Guide",
    "title": "Guide",
    "category": "page",
    "text": ""
},

{
    "location": "functions.html#Main-functions-1",
    "page": "Guide",
    "title": "Main functions",
    "category": "section",
    "text": "There are several functions programmed in this library, however many others can be programmed from the existing ones. We will therefore highlight the main functions that form the core of LIGA.  Operations such as sum and difference are intuitive and applicable to all structures (BasisBlade, MutiVector and Blades). These operations are pretty simple so, we can ommit here. "
},

{
    "location": "functions.html#NewLiga.geometric",
    "page": "Guide",
    "title": "NewLiga.geometric",
    "category": "function",
    "text": "geometric(a::GAAbstractType,b::GAAbstractType) (Liga function)\n\nThis function is used to calculate the geometric product of GAAbstractType elements (Basis Blade, MultiVector or Blades). It returns a MultiVector type.\n\nExample (in 𝔾₃)\n\njulia> geometric(e12,e12)\n\nreturns a multivector  -1.0id.\n\nNote that we can use other structs in arguments. \n\nExample (in 𝔾₃)\n\njulia> geometric(1.0+3.0*e12,e1)\n\nreturns ( 1.0 e1 ) - ( 3.0 e2).\n\nExample (in 𝔾₄₁ conformal)\n\njulia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0\n\njulia> v=3.0ep\n\njulia> geometric(u,v)\n\nreturns ( 16.5 id ) + ( 6.0 e+ ) + ( 4.5 e1+ ) - ( 9.0 e2+ ) - ( 7.5 e+- ).\n\nYou can use the infix operator ∘ instead of geometric command.\n\nExample (in 𝔾₃)\n\njulia> u = 1.0+3.0*e12\n\njulia> v = e1\n\njulia> u ∘ v\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.inner",
    "page": "Guide",
    "title": "NewLiga.inner",
    "category": "function",
    "text": "inner(a::GAAbstractType,b::GAAbstractType) (Liga function)\n\nThis function is used to determine the inner product of GAAbstractType elements (Basis Blade, MultiVector or Blades).\n\nExample (in 𝔾₃)\n\njulia> inner(1.0+3.0*e12,e1)\n\nreturns - ( 3.0 e2 ).\n\nExample (in 𝔾₄₁ conformal)\n\njulia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0\n\njulia> v=3.0*ep\n\njulia> inner(u,v)\n\nreturns ( 16.5 id ).\n\nAn alternative way to run this functions is using the infix  operator (⋅).\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.outer",
    "page": "Guide",
    "title": "NewLiga.outer",
    "category": "function",
    "text": "outer(a::GAAbstractType,b::GAAbstractType) (Liga function)\n\nThis function is used to determine the outer product of GAAbstractType elements (Basis Blade, MultiVector or Blades).\n\nExample (in 𝔾₃)\n\njulia> outer(1.0+3.0*e12,e1)\n\nreturns ( 1.0 e1 ).\n\nExample (in 𝔾₄₁ conformal)\n\njulia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0\n\njulia> v=3.0*ep\n\njulia> outer(u,v)\n\nreturns ( 6.0 e+ ) + ( 4.5 e1+ ) - ( 9.0 e2+ ) - ( 7.5 e+- ).\n\nAn alternative way to run this function is using the infix  operator (^).\n\njulia> u=2.0*id + 1.5*e1 -3.0*e2 +4.0*e∞ -3.0*e0\n\njulia> v=3.0*ep\n\njulia> u^v\n\n\n\n\n\n"
},

{
    "location": "functions.html#Geometric,-Inner-and-Outer-products-1",
    "page": "Guide",
    "title": "Geometric, Inner and Outer products",
    "category": "section",
    "text": "geometric\n\ninner \n\nouter"
},

{
    "location": "functions.html#NewLiga.grade",
    "page": "Guide",
    "title": "NewLiga.grade",
    "category": "function",
    "text": "grade(b::BasisBlade) (Liga function)\n\nThis function is used to determine the grade of a basis blade element\n\nExample (in 𝔾₃₁)\n\njulia> grade(e123)\n\nreturn the value (Int) 3.\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.gradeplus",
    "page": "Guide",
    "title": "NewLiga.gradeplus",
    "category": "function",
    "text": "gradeplus(b::BasisBlade) (Liga function)\n\nThis function is used to determine the grade plus a Basis Blade element. It returns an integer.\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.grademinus",
    "page": "Guide",
    "title": "NewLiga.grademinus",
    "category": "function",
    "text": "grademinus(b::BasisBlade) (Liga function)\n\nThis function is used to determine the grade minus of a basis blade element It returns an integer.\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.gradeprojection",
    "page": "Guide",
    "title": "NewLiga.gradeprojection",
    "category": "function",
    "text": "gradeprojection(b::BasisBlade,k::Int) (Liga function)\n\nThis function is used to determine the grade projection of a basis blade. This function returns a Basis Blade element. \n\nExample (in 𝔾₃₁)\n\njulia> gradeprojection(-e123,3)\n\nreturns the Basis Blade -e123\n\n\n\n\n\ngradeprojection(b::MultiVector,k::Int) (Liga function)\n\nThis function is used to determine the grade projection of a multivector. This function returns a MultiVector element.\n\nExample (in 𝔾₃₁)\n\njulia> gradeprojection(1.0*id+e2-2.0*e12p,3)\n\nreturns the multivector  - ( 2.0 e12+ ).\n\n\n\n\n\n"
},

{
    "location": "functions.html#Grade-functions-1",
    "page": "Guide",
    "title": "Grade functions",
    "category": "section",
    "text": "In most cases, this function handles with BasisBlade and returns a structure of the same type. The exception is the function gradeprojection.grade\n\ngradeplus \n\ngrademinus\n\ngradeprojection"
},

{
    "location": "functions.html#Base.reverse-Tuple{BasisBlade}",
    "page": "Guide",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(b::BasisBlade) (Liga function)\n\nThis function calculates the reverse of a Basis Blade. For example,\n\njulia> u = reverse(e12)\n\nreturns the BasisBlade -e12.\n\n\n\n\n\n"
},

{
    "location": "functions.html#Base.reverse-Tuple{MultiVector}",
    "page": "Guide",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(V::MultiVector) (Liga function)\n\nThis function calculates the reverse of a MultiVector. It returns a new MultiVector.\n\nExample (in 𝔾₃₁)\n\njulia> u=1.0*id -2.0*e1 -5.0*e∞ +3.0*e∞0 \njulia> reverse(u)\n\nreturns the multivector (1.0 id) - ( 2.0 e1) - ( 5.0 e∞ ) - ( 3.0 e∞0 ).\n\n\n\n\n\n"
},

{
    "location": "functions.html#Base.reverse-Tuple{Blade}",
    "page": "Guide",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(B::Blade) (Liga function)\n\nThis function calculates the reverse of a Blade. It returns a new Blade. \n\nExample (in 𝔾₃₁)\n\njulia> B=Blade([[1.0,0.0,0.0],[0.0,1.0,0.0]],2,true)\njulia> reverse(B)\n\nreturns [0.0, 1.0, 0.0]∧[1.0, 0.0, 0.0].\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.conjugate-Tuple{BasisBlade}",
    "page": "Guide",
    "title": "NewLiga.conjugate",
    "category": "method",
    "text": "conjugate(b::BasisBlade) (Liga function)\n\nThis function calculates the conjugate of a Basis Blade.\n\nExample (in 𝔾₃₁)\n\njulia> u = conjugate(e12)\n\nreturns the BasisBlade -e12.\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.conjugate-Tuple{MultiVector}",
    "page": "Guide",
    "title": "NewLiga.conjugate",
    "category": "method",
    "text": "conjugate(V::MultiVector) (Liga function)\n\nThis function calculates the conjugate of a MultiVector and return a multivector.\n\nExample (in 𝔾₃₁)\n\njulia> u = -1.0+3.0*e1-1.0*e12\n\njulia> conjugate(u)\n\nreturns the MultiVector - ( 1.0 id ) + ( 3.0 e1 ) + ( 1.0 e12 ).\n\n\n\n\n\n"
},

{
    "location": "functions.html#NewLiga.conjugate-Tuple{Blade}",
    "page": "Guide",
    "title": "NewLiga.conjugate",
    "category": "method",
    "text": "conjugate(B::Blade) (Liga function)\n\nThis function calculates the conjugate of a Blade and return a new Blade.\n\nExample (in 𝔾₃₂)\n\njulia> u = blade([1.0,2.0,3.0,4.0,5.0],[-1.0,-1.0,-1.0,-1.0,3.0])\n\njulia> conjugate(u)\n\nreturns the Blade [-1.0, -1.0, -1.0, 1.0, -3.0]∧[1.0, 2.0, 3.0, -4.0, -5.0].\n\n\n\n\n\n"
},

{
    "location": "functions.html#Involutions-1",
    "page": "Guide",
    "title": "Involutions",
    "category": "section",
    "text": "The first involution that we highlight is the reverse. Although this function can be used by any structure, it is important to note that it returns an object of the same structure (not a general multivector as in other functions).\nreverse(b::BasisBlade)\n\nreverse(b::MultiVector)\n\nreverse(A::Blade)\nAnother important involution is the conjugate which is avaliable just for BasisBlade and MultiVector. (Need to extend to Blade).\nconjugate(b::BasisBlade)\n\nconjugate(b::MultiVector)\n\nconjugate(A::Blade)\n"
},

{
    "location": "functions.html#NewLiga.dual-Tuple{GAAbstractType}",
    "page": "Guide",
    "title": "NewLiga.dual",
    "category": "method",
    "text": "dual(b::GAAbstractType) (Liga function)\n\nThis function calculates the dual of a GAAbstractType and returns a multivector element.\n\nExample (in 𝔾₃₁)\n\njulia> u = dual(e12)\n\nreturns the multivector ( 1.0 e+- ).\n\nExample (in 𝔾₃₂)\n\njulia> b=blade([1.0,0.0,3.0,0.0,5.0],[0.0,-1.0,-1.0,-1.0,0.0])\n\nreturns the multivector  - ( 5.0 e123 ) - ( 5.0 e124 ) + ( 5.0 e134 ) - ( 3.0 e125 ) - ( 1.0 e235 ) - ( 3.0 e145 ) - ( 1.0 e245 ) + ( 1.0 e345 ).\n\n\n\n\n\n"
},

{
    "location": "functions.html#Duality-1",
    "page": "Guide",
    "title": "Duality",
    "category": "section",
    "text": "Duality is an important concept within GA and depends explicitly on the definition of pseudo-scalar which is a variable given in Liga module by pseudo. COMPLEMENTAR COM ALGO AQUI!With the definition of pseudo scalar in mind, we can define the duality function by the geometric product of the entity (BasisBlade, MultiVector or Blade) with this pseudo scalar.\ndual(a::GAAbstractType)\n"
},

{
    "location": "functions.html#NewLiga.multivector",
    "page": "Guide",
    "title": "NewLiga.multivector",
    "category": "function",
    "text": "multivector(b::BasisBlade) (Liga function)\nmultivector(a::Number) (Liga function)\nmultivector(v::Vector{Float64}) (Liga function)\nmultivector(A::Blade) (Liga function)\n\nThis function is used to create a MultiVector type  from numbers, vectors, basis blade or Blades. \n\nExample (in 𝔾₃₁)\n\njulia> multivector(e1)\n\nreturns the multivector (1.0e1).\n\nExample (in 𝔾₃₁)\n\njulia> multivector([1.0,2.0,3.0])\n\nreturns the multivector (1.0 e1)+(2.0 e2)+(3.0 e3).\n\nExample (in 𝔾₃₁)\n\njulia> multivector(Blade([[1.0,2.0,3.0]],1,true))\n\nreturns the multivector (1.0 e1)+(2.0 e2)+(3.0 e3).\n\n\n\n\n\n"
},

{
    "location": "functions.html#Convertions-1",
    "page": "Guide",
    "title": "Convertions",
    "category": "section",
    "text": "In Liga is possible convert a Number, Vector, BasisBlade or Blade into MultiVector type. multivector"
},

{
    "location": "functions.html#Projective-and-Conformal-Spaces-1",
    "page": "Guide",
    "title": "Projective and Conformal Spaces",
    "category": "section",
    "text": ""
},

]}
