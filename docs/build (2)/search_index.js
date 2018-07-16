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
    "text": "This page is devoted to document and make easier the use of LIGA- Library for  Geometric Algebra. Our intent is provide a good experience of use for students  and reseachers interesting in use (and programming) a Geometric Algebra  framework. All scripts were made in Julia language,  version 0.6. There is a miscellany of libraries for GA study. Our intent here is not to provide software to compete with such libraries in performance or generality. Our goal is to provide a relatively simple, open, easy to install and focused environment for Julia users.It\'s important a previous contact with Geometric Algebra in theoretical aspects. In this sense, we recommend the following references:-Christian Perwass, Geometric Algebra with Applications in Engineering,  Springer Series in Geometry and Computing, 2009.-Valter Soares Camargo, Conformal Geometric Algebra and Distance Geometry, (Doctoral Thesis). IMECC-UNICAMP, 2015. "
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
    "location": "index.html#To-do-things-1",
    "page": "Overview",
    "title": "To do things",
    "category": "section",
    "text": "Extend Inversion to more types (cmultvec and pmultivec)\nExtend Projection to more types (kmultvec,cmultvec,pmultivec)\nExtend Rejection to more types (kmultvec,cmultvec,pmultivec)\nImprovements in matrix products.\nPut more examples in Conformal Space Section about: S, Hm, conformal"
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
    "text": "(In the future) Open Julia REPL and type Pkg.add(\"Liga\")Or (by now) install it yourself as: git clone https://github.com/evcastelani/Liga.jl"
},

{
    "location": "gettingstarted.html#Layouts-1",
    "page": "Getting Started",
    "title": "Layouts",
    "category": "section",
    "text": "As cited, just the second option is avaiable. So, you can type using Liga In order to start to use Liga, the first point that you must be in mind is setup the layout. The layout function defines the space and all logical elements for build $ \\mathbb{G}_k $. For example, if we want create the $\\mathbb{G}_3$ space, we need just typelayout(3)In this case, we are generating a logical structure for create basis of vectors and multivectors for $\\mathbb{G}_3$. Naturally, we can work (with appropriated setup of elements), with $\\mathbb{G}_4$ and $\\mathbb{G}_{4,1}$. Let us started with more basic example. After type layout(3)  was created base logical elements. Type, for examplee1 as we can see e1 is a logical array with value  [false,true,true]All base elements were created using this idea. Try type others logical basis elementse123\n\ne12       An  observation in this case is the order of the index of these elements are important. If you typee21an error mensage is showed. This does not mean that this element can not be built. We will see this later. There is a subtle difference between logical elements and k-base elements. The k-basis element can provide more information, like a logical base element and a scalar associated. To create a k-basis elements we need to use the command kb. For example, typekb(e1)and aftertypeof(kb(e1))Naturally, kb(e1) is a k-basis element. In this case, in order to acess the information related to kb(e1) just ypekb(e1).e\nkb(e1).sclSo, we can see that an kb element needs a logical basis element and a scalar element. We can generate others kb elements by typekb(e1,-2.0)There are other syntaxes functions that are equivalent to command kb. For example,kb([true,false,false])or kbasis(e1,1.0)It is obvious that kb is abbreviation of k-basis element and if we think in pratical purposes, it is just a logical vector. In current version, kb supports a more direct definition. For example, if we want to define the element kb(e12,2.3) we can just type in Julia REPL, the command-2.3*e12wich is equivalent tokb([true,true,false],-2.3)Concerning the diference between e12 and e21, if by some reason, we need to define the e21 element, we can just define the element kb(e12,-1.0).Lets to consider a more complete example. The most import function in GA context is the geometric product. So, we can show an example related to this function. For this, just type these commands in Julia REPL \nu=5.0*e123\n\nv=1.0*id+3.0*e12\n\nv∘u #(v \\circ+TAB u)and a gemetrical product between the multivector u and the multivector v is showed. Note that, we didn\'t define a multivector type, but it is very intuitive that they are composed of k-basis elements.All these commands provide the basics for starting more elaborate constructions in the space that has been defined.   Consequently, we recommend that the reader calmly review documentation on the types defined and functions supported. Naturally, we introduce the gemetrical product in $G_3$ but Liga provides more functions and other kind of elements."
},

{
    "location": "types.html#",
    "page": "Types",
    "title": "Types",
    "category": "page",
    "text": "using Liga\nlayout(3)"
},

{
    "location": "types.html#Types-1",
    "page": "Types",
    "title": "Types",
    "category": "section",
    "text": "For each space we have three types. In all spaces, geometric algebra  elements can be represented divided in basis, multivectors and blades. All operations and how each type works will be explained later in this guide."
},

{
    "location": "types.html#Basis-types-1",
    "page": "Types",
    "title": "Basis types",
    "category": "section",
    "text": ""
},

{
    "location": "types.html#k-basis-(kbasis)-1",
    "page": "Types",
    "title": "k-basis (kbasis)",
    "category": "section",
    "text": "The kbasis type is the first of all types. All of the other types depend on this in some level.This type represents the basis elements of the geometric algebra of $\\mathbb{R}^k$, it depends of a k-dimensional logical array and a scalar, generallykbasis(x::Vector{Bool}, a::Number)will create a basis element from x and a. For example,kbasis([false, true, false], 3.5)return the element 3.5e2 wich is in $\\mathbb{G}_3$ space.As commented in Layout Section it\'s possible to create a kbasis using the base operator * between a Number and a Vector{Bool}, as follows2.0*[false, true, true]will return2.0e23and, using the layout elements, it can be used2.0*e23to create the same element above."
},

{
    "location": "types.html#p-basis-(pbasis)-1",
    "page": "Types",
    "title": "p-basis (pbasis)",
    "category": "section",
    "text": "To create an element of this type is required, in addition to the logical array and scalar, one more logical element, that is,pbasis(x::Vector{Bool}, y::Bool, a::Number)where the second parameter represents the basis element that squares to -1. For example,pbasis([true, false, true, true], true, 1.0)returns the element of $ \\mathbb{G}_{4,1} $1.0e13+-where the basis element e3 is replaced by e+ with the same properties and is added the basis element e- that squares to -1."
},

{
    "location": "types.html#c-basis-(cbasis)-1",
    "page": "Types",
    "title": "c-basis (cbasis)",
    "category": "section",
    "text": "This _type_ instead of use e- as addition basis element, the elements e∞ and e∘ are used. To call an element of this type it is necessary, besides the logical array and scale, two logical elements, as followscbasis(x::Vector{Bool}, y::Bool, z::Bool, a::Number)the second parameter represents the basis element e and the third represents e∘. As in the previous cases, instead of use the Vector{Bool} can be used the elements generated by the layout function, that is, if was called the layout(3) before, thencbasis(e12, true, false, 2.0)will return2.0e12∞Shortly e∞ = e- + e+ and e∘ = 0.5(e- - e+) hence e∞² = e∘² = 0 and e∞e∘ = -1."
},

{
    "location": "types.html#Multivector-types-1",
    "page": "Types",
    "title": "Multivector types",
    "category": "section",
    "text": "The multivectors types are composed of a vector of basis elements, each coordinate of the vector represents a portion of the sum that will generate the multivector itself."
},

{
    "location": "types.html#k-multivector-(kmultvec)-1",
    "page": "Types",
    "title": "k-multivector (kmultvec)",
    "category": "section",
    "text": "This type is composed of a kbasis vector and represents an arbitrary element of the geometric algebra $\\mathbb{G}_k$.kmultvec(X::Vector{kbasis})once created the kmultvec its coordinates elements will be arranged as a sum,kmultvec([kbasis(e1, 1.0), kbasis(id, 2.4), kbasis(e123, -1.0)])will create the elmemnt1.0e1 + 2.4 + -1.0it can be created by using the sum directly as followskbasis(e1, 1.0) + kbasis(id, 2.4) + kbasis(e123, -1.0)since the sum of kbasis, as will be explained later, generates an element of type kmultvec in most cases."
},

{
    "location": "types.html#p-multivector-(pmultvec)-1",
    "page": "Types",
    "title": "p-multivector (pmultvec)",
    "category": "section",
    "text": "The only difference here is that instead of being composed of a kbasis vector the pmultvec is composed of a vector pbasis,pmultvec(X::Vector{pbasis})but, in the same way as the previous one, its coordinates elements  will be arranged as a sum."
},

{
    "location": "types.html#c-multivector-(cmultvec)-1",
    "page": "Types",
    "title": "c-multivector (cmultvec)",
    "category": "section",
    "text": "Again, the difference remains in the vector type, which in this case will be a vector composed of cbasis elements,cmultvec(X::Vector{cbasis})Other than that, everything follows in the same way."
},

{
    "location": "types.html#Blade-*types*-1",
    "page": "Types",
    "title": "Blade types",
    "category": "section",
    "text": "These types represent elements of the geometric algebra called blades, for this we need a set 1.i. of 1-vectors."
},

{
    "location": "types.html#k-Blades-(kblade)-1",
    "page": "Types",
    "title": "k-Blades (kblade)",
    "category": "section",
    "text": "The kblade type represents a blade in $\\mathbb{G}_k$ space and is composed of a kmultvec vectorkblade(A::Vector{kmultvec})and, if the introduced vector is composed of l.i 1-vectors, this will generate a blade formed by the outer product of the coordinates of this vector. For example,kblade([kmultvec(kbasis(e1)), kmultvec(kbasis(e2)), kmultvec(kbasis(e3))])generates the kblade(e1)∧(e2)∧(e3)"
},

{
    "location": "types.html#p-Blades-(pblade)-1",
    "page": "Types",
    "title": "p-Blades (pblade)",
    "category": "section",
    "text": "This type represents a blade in (k,1) geometric algebra space the only difference between this type and the previous one is that it is composed of a pmultvec vector.pblade(A::Vector{pmultvec})however, as in the previous case, this will generate a blade formed by the outer product of the coordinates of the vector."
},

{
    "location": "types.html#c-Blades-(cblade)-1",
    "page": "Types",
    "title": "c-Blades (cblade)",
    "category": "section",
    "text": "As for kblade and pblade this type represents a blade in (k+1,1) geometric algebra space and the difference is in the type that composes it, in this case cmultvec.cblade(A::Vector{cmultvec})and as in previous cases, this will generate the blade composed by the outer product of its coordinates.All basis types can be called using the functions kb for kbasis, pb to pbasis and cb for cbasis"
},

{
    "location": "Gkspace.html#",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "$\\mathbb{G}_k$ Space",
    "category": "page",
    "text": "using Liga\nlayout(3)"
},

{
    "location": "Gkspace.html#\\mathbb{G}_k-Space-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "$\\mathbb{G}_k$ Space",
    "category": "section",
    "text": "The Euclidean geometric algebra $G_k$ space is the basic space of this project. The main operations of the space, that is, the Clifford or geometric product, the inner, outer and scalar products, the reverse, dual, projection and rejection of elements of this space are programmed and ready for use.To create an element of this space is necessary to call the types which have only elements of type kbasis in their composition, that is, the elements of the type kbasis itself, the kmultvec and the kblade types. Functions and operations related to this space  can be used in projetive and conformal space, but we explain in Projetive Space Section this.If the reader has not yet read about the types supported by this library, reading in Section Types is highly recommended."
},

{
    "location": "Gkspace.html#Functions-and-Operations-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Functions and Operations",
    "category": "section",
    "text": "Here we will give an explanation of how each function works, in terms of what is taken as the input parameter and what the function will return, so that some detailed examples are given later in order to simplify the understanding of the program and avoid possible mistakes."
},

{
    "location": "Gkspace.html#Operations-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Operations",
    "category": "section",
    "text": ""
},

{
    "location": "Gkspace.html#Geometric-Product-geoprod(a,b)-or-a-b-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Geometric Product - geoprod(a,b) or a ∘ b",
    "category": "section",
    "text": "As in the literature, the Clifford or geometric product is the main operation of the space, and is the basis for most of the other operations. It takes two elements of the type kbasis, kmultvec or kblade and return their geometric product.geoprod(a::kbasis, b::kbasis)\ngeoprod(X::kmultvec, Y::kmultvec)\ngeoprod(A::kblade, B::kblade)To simplify the operation input it can be used the base operator ∘ instead of call the function geoprod. Base.:∘(a::kbasis, b::kbasis)\n Base.:∘(X::kmultvec, Y::kmultvec)\n Base.:∘(A::kblade, B::kblade)In order to obtain the geometric product of elements of different types,as long as the operator is used, it is not necessary to convert them before . As explained above, since a generic element of geometric algebra is a multivector, the elements of type kbasis and kblade can be converted to the type kmultvec and thus can be apply the geometric product operation. These conversions are done automatically when two parameters of different types are used as input.See Examples section below to better understand the input and output of operations."
},

{
    "location": "Gkspace.html#Inner-Product-inner(a,b)-or-a-b-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Inner Product - inner(a,b) or a ⋅ b",
    "category": "section",
    "text": "The input parameters for the inner product are of the same types as for the geometric product, so to obtain the inner product between two elements of type kbasis, kmultvec or kblade just use the inner function. inner(a::kbasis, b::kbasis)\n inner(X::kmultvec, Y::kmultvec)\n inner(A::kblade, B::kblade)As for the geometric product, to obtain the inner product is sufficient to use a base operator, in this case, the operator ⋅ instead of call the inner function. Base.:⋅(a::kbasis, b::kbasis)\n Base.:⋅(X::kmultvec, Y::kmultvec)\n Base.:⋅(A::kblade, B::kblade)"
},

{
    "location": "Gkspace.html#Outer-Product-outer(a,b)-or-a-b-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Outer Product - outer(a,b) or a ^ b",
    "category": "section",
    "text": "As for the geometric and inner products, one can use the function outer to obtain the outer product between two elements of type kbasis, kmultvec or kblade. outer(a::kbasis, b::kbasis)\n outer(X::kmultvec, Y::kmultvec)\n outer(A::kblade, B::kblade)The base operator that can be used instead of call the function in this case is the ^ symbol. Base.:^(a::kbasis, b::kbasis)\n Base.:^(X::kmultvec, Y::kmultvec)\n Base.:^(A::kblade, B::kblade)"
},

{
    "location": "Gkspace.html#Scalar-Product-scalar(a,b)-or-a-*-b-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Scalar Product - scalar(a,b) or a * b",
    "category": "section",
    "text": "The scalar function is not different from the other three mentioned above, in order to obtain the scalar product between two elements of type kbasis, kmultvec or kblade just use the function inner. scalar(a::kbasis, b::kbasis)\n scalar(X::kmultvec, Y::kmultvec)\n scalar(A::kblade, B::kblade)Or using the base operator *. Base.:*(a::kbasis, b::kbasis)\n Base.:*(X::kmultvec, Y::kmultvec)\n Base.:* (A::kblade, B::kblade)In this operation, unlike the previous ones, the output parameter will always be of type Number."
},

{
    "location": "Gkspace.html#Sum-and-difference-mvsum(a,b)-or-a-b-and-a-b-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Sum and difference - mvsum(a,b) or a + b and a - b",
    "category": "section",
    "text": "This function has a particularity in the type of output, since the sum of two elements of type kbasis will only be of type kbasis if they are l.d, otherwise the output will be of type kmultvec.To obtain the sum of two elements of types kbasis,  kmultvec or kblade, one can call the function mvsum, as well as the previous functions, mvsum(a::kbasis, b::kbasis)\n mvsum(X::kmultvec, Y::kmultvec)\n mvsum(A::kblade, B::kblade)or using the base operator + Base.:+(a::kbasis, b::kbasis)\n Base.:+(X::kmultvec, Y::kmultvec)\n Base.:+(A::kblade, B::kblade)but, to make it easier, in addition to the base operator +, it can be used the base operator - which convert the second element to its additive inverse. Base.:-(a::kbasis, b::kbasis)\n Base.:-(X::kmultvec, Y::kmultvec)\n Base.:-(A::kblade, B::kblade)thus obtaining the difference between the two elements."
},

{
    "location": "Gkspace.html#Functions-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Functions",
    "category": "section",
    "text": ""
},

{
    "location": "Gkspace.html#Grade-grade(a)-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Grade - grade(a)",
    "category": "section",
    "text": "The grade function take as input an element of the type kbasis and return its grade which is an element of the type Int, it\'s used more as an auxiliary function. grade(a::kbasis)"
},

{
    "location": "Gkspace.html#Dual-dual(a)-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Dual - dual(a)",
    "category": "section",
    "text": "The dual function, as the name suggests, returns its dual element, taking as input elements of types kbasis,  kmultvec or kblade. dual(a::kbasis)\n dual(X::kmultvec)\n dual(A::kblade)As for the products and sum, in case the input of the kblade the output will be of the type kmultvec."
},

{
    "location": "Gkspace.html#Reverse-mvreverse(a)-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Reverse - mvreverse(a)",
    "category": "section",
    "text": "Returns the reverse of the input element, but, unlike other functions, this will always return an element of the same input type, being it kbasis,  kmultvec or kblade. mvreverse(a::kbasis)\n mvreverse(X::kmultvec)\n mvreverse(A::kblade)The conjugation function isn\'t necessary in $G_k$ space, since its depends of negative grade."
},

{
    "location": "Gkspace.html#Magnitude-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Magnitude",
    "category": "section",
    "text": "The magnitude function receives as input an element of type kbasis,  kmultvec or kblade and returns its magnitude of the type Ǹumber. magnitude(a::kbasis)\n magnitude(X::kmultvec)\n magnitude(A::kblade)"
},

{
    "location": "Gkspace.html#Special-functions-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Special functions",
    "category": "section",
    "text": "In this space, some functions will work only with elements of the type kblade, one being the inversion function and the other two depending on the first one, since not every multivector of Gₖ has inverse with respect to the geometric product."
},

{
    "location": "Gkspace.html#Inversion-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Inversion",
    "category": "section",
    "text": "This function returns the inverse element of the input with respect to the geometric product. inverse(A::kblade)The output type will also be kblade."
},

{
    "location": "Gkspace.html#Projection-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Projection",
    "category": "section",
    "text": "The projection is a simple function that depends on two elements of the type kblade projection(A::kblade, N::kblade)and it returns, in this case, the projection of A onto N, which will be an element of the type kmultvec."
},

{
    "location": "Gkspace.html#Rejection-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Rejection",
    "category": "section",
    "text": "The rejection function also receives two elements of the type kblade as input rejection(A::kblade, N::kblade)and, as for the projection, it returns an element of the type kmultvec which will be the rejection of the kblade A from the kblade  N."
},

{
    "location": "Gkspace.html#Examples-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Examples",
    "category": "section",
    "text": "All of the following examples are based on layout(3), that is, the space generated will be $\\mathbb{G}_3$ depending on the objects id, e1, e2, e3, e12, e13, e23 and e123, however everything follows in the same way for other dimensions."
},

{
    "location": "Gkspace.html#Grade-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Grade",
    "category": "section",
    "text": "The grade function is the most basic of them, as mentioned above, it\'s used more like an auxiliary function for the other and it returns the grade of a kbasis element as follows:grade(kbasis(e12, 1.0))which is the grade of e12.Another way to do this is by calling the element before the function\na = kbasis(e123, -2.3)\nand then use the grade functiongrade(a)which is the grade of -2.3e123."
},

{
    "location": "Gkspace.html#Geometric,-inner,-outer-and-scalar-products-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Geometric, inner, outer and scalar products",
    "category": "section",
    "text": "The geometric, inner and outer products in case the input elements are of type kd or kmultvec, will always return the same type, even if the result is a scalar or a basis blade, that is, if the following elements are called:\na = kbasis(e3, 2.0)\n\nb = kbasis(e12, -1.0)\n\nc = kbasis(id, 0.0)\n\nd = kbasis(e123)when operated, for examplegeoprod(a, b)which is a kbasis element, although when operatedinner(a, b)this element is seen as an object of type kbasis. Note that if a second parameter is not placed when called the kbasis type it will be generated with the scalar 1.0.As mentioned before, we can use the base operators ∘, ⋅ and ^ to calculate these products instead of call the functionsa ^ c\n\nb ⋅ d\nThe same is true if the input is of the type kmultvec as follows:X = kmultvec([kbasis(e12,2.0),kbasis(e123,-1.5),kbasis(e3,4.0)])\n\nY = kmultvec([a, b, kbasis(id, 2.0)])\n\nZ = kmultvec(d)and these can be operated in the same way as the previous onesX ∘ Y\nY ⋅ Z\nY ^ ZNow, in the case of kblade type, everything follows like the previous cases, the difference is in the output of the function. For example, if we call the following kblade elementsA = kblade([kmultvec([kb(e1)]), kmultvec([kb(e2)]), kmultvec([kb(e3)])])\nB = kblade([kmultvec([kb(e1, 2.0), a]), kmultvec([kb(e1, 3.0)])])\nC = kblade([kmultvec(a)])\n\nA ⋅ B\nreturns an element of the type  kmultvec. In fact, type\ntypeof(A ⋅ B)the same will happen for the other operations B ^ Creturns 0.0and B ∘ Creturns -12.0e1The scalar product follows a different rule, in all cases the output will always be of Number. For example, if we want to calculate the scalar product between the kbasis element d and the kblade A previously called simply do as follows:d=kbasis(id,-2.0)\nd * Aremember that this is an element of type Number.typeof(ans)You can operate any other elements of types kbasis, kmultvec and kblade,   a * a\n b * X\n B * B\n X * A"
},

{
    "location": "Gkspace.html#Sum-and-difference-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Sum and difference",
    "category": "section",
    "text": "The mvsum operation works in the same way as the previous ones with respect to the input parameters, the difference is in the output parameters, which in the case of the sum of elements of the type kbasis the output there are two possible outputsmvsum(kbasis(e23, 2.3), kbasis(e23, 2.7))in this case, it may be noted that the parameter e23 is used in the construction of both inputs, then the it returns 5.0e23which is of the type kbasis,typeof(ans)However, if we calculate the following sum a + bit returns 2.0e3 + -e12which is of the type kmultvec, remembering that a and b are the elements kbasis(e3, 2.0) and kbasis(e12, -1.0) called above.Other than that, the rest follows as before\nA + B\n\ntypeof(ans)\n\nX - C\n\ntypeof(ans)\n\nd + Y\n\ntypeof(ans)\n\na + d\n\ntypeof(ans)\nNote that here we can use either the base operator + or -.Again, the output when applied the sum and the difference  between elements of the type kblade will be of the type kmultvec."
},

{
    "location": "Gkspace.html#Dual-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Dual",
    "category": "section",
    "text": "The dual function returns the dual of the input element as followsdual(a)\n\ntypeof(dual(a))\n\ndual(X)\n\ntypeof(dual(X))\n\ndual(A)\n\ntypeof(dual(A))As for the previous operations, except for the scalar product, in case the input is of type kblade the output will be of type kmultvec."
},

{
    "location": "Gkspace.html#Reverse-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Reverse",
    "category": "section",
    "text": "Unlike other functions and operations the mvreverse function returns an element of the same type as the input, be it kbasis, kmultvec or kblademvreverse(a)\n\ntypeof(mvreverse(a))\n\nmvreverse(X)\n\ntypeof(mvreverse(X))\n\nmvreverse(A)\n\ntypeof(mvreverse(A))"
},

{
    "location": "Gkspace.html#Magnitude-2",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Magnitude",
    "category": "section",
    "text": "The magnitude in $\\mathbb{G}_k$ is nothing more than the square root of the scalar product of an element by its reverse, so it is, like the scalar product, always returns an element of type Float. For example,magnitude(a)\n\ntypeof(magnitude(a))\n\nmagnitude(X)\n\ntypeof(magnitude(X))\n\nmagnitude(A)\n\ntypeof(magnitude(A))"
},

{
    "location": "Gkspace.html#Inversion-2",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Inversion",
    "category": "section",
    "text": "In current version of Liga, this important functions works just for kblade and  kmultvec types. We are improve to more types. The return is an equivalent type to  input.\ninverse(A)\n\ntypeof(inverse(A))\n\ninverse(B)\n\ntypeof(inverse(B))\n\ninverse(C)\n\ntypeof(inverse(C))\nit returns the inverse elements of A,B and C.We can test the function using the geometric product of a kblade by its inverse, and it should return 1.0\nA ∘ inverse(A)\n\nB ∘ inverse(B)\n\nC ∘ inverse(C)\nNote that as mentioned when dealing with the geometric product, this is a element of type kmultvectypeof(B ∘ inverse(B))Let us see an example related to k-multivector. \ninverse(X)\n\ntypeof(inverse(X))\n\nX ∘ inverse(X)\n\ntypeof(X ∘ inverse(X))\nNote that, in X ∘ inverse(X) the result was 2.7755575615628914e-17e12-1.3877787807814457e-17e3+1.0 wich is approximated 1.0 (and kmultvec);"
},

{
    "location": "Gkspace.html#Projection-2",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Projection",
    "category": "section",
    "text": "This function depends on the inversion and by now can be used only with elements of type kblades as input. For example, projection(A, B)\n projection(B, A)\n projection(A, C)will return e123\n 6.0e13\n e123these output elements are of type kmultvec"
},

{
    "location": "Gkspace.html#Rejection-2",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Rejection",
    "category": "section",
    "text": "Like the above function Projection, the rejection receives an element of type kblade and returns an element of type kmultvec which is the rejection of the first input parameter from the second. rejection(B, kblade([kmultvec(kbasis(e2))]))will return 6.0e13which is the difference between the kblade B and the projection rejection(B, kblade([kmultvec(kbasis(e2))]))."
},

{
    "location": "Gkspace.html#Base.transpose-Tuple{Array{Liga.kbasis,1}}",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Base.transpose",
    "category": "method",
    "text": "transpose(A::Array{kbasis,1})\ntranspose(A::Array{kbasis,2})\ntranspose(A::Array{kmultvec,1})\ntranspose(A::Array{kmultvec,2})\n\nReturns the transpose of an Array defined by kbasis or kmultvec elements.\n\n\n\n"
},

{
    "location": "Gkspace.html#Matrix-operations-1",
    "page": "$\\mathbb{G}_k$ Space",
    "title": "Matrix operations",
    "category": "section",
    "text": "The main operations can be extended to array like the classic matrix product. For example, considering the geometric productu=[1.0*e12;2.0*e123]\n\ntranspose(u)∘u\n\nA=[1.0*e1+2.0*e123 2.0*e1+1.0*id; 0.5*e123+1.0*e3 1.0*id+1.0*e1]\n\nB=[5.0*e12+2.0*e123 2.0*e1-1.0*id; 1.0*e123-2.0*e2 1.0*id+1.0*e1]\n\nA∘BNote in this example that additionally, a transpose function for array defined by kbasis or kmultvec types was implemented.transpose(A::Array{kbasis,1})"
},

{
    "location": "projspace.html#",
    "page": "Projective Space ($\\mathbb{G}_{k+1}$)",
    "title": "Projective Space ($\\mathbb{G}_{k+1}$)",
    "category": "page",
    "text": "using Liga\nlayout(3)"
},

{
    "location": "projspace.html#Projective-Space-(\\mathbb{G}_{k1})-1",
    "page": "Projective Space ($\\mathbb{G}_{k+1}$)",
    "title": "Projective Space ($\\mathbb{G}_{k+1}$)",
    "category": "section",
    "text": "The projective space is a subspace of $\\mathbb{G}_{k+1}$, which is nothing more than $\\mathbb{G}_k$ space with an extra dimension, so there is no need to create a new type for this.This space can be used through embedding called homogenization, it takes a vector in $ \\mathbb{R}^n$into its embedding in the affine space.H(X::Vector{Float64})\nH(X::kmultvec)Since homogenization takes an element of $\\mathbb{R}^n$ to its affine space, it is possible to enter an element of the type kmultvec as input parameter for the function H, however this must be in the Euclidean subspace of $\\mathbb{G}_k$, that is, its components must be both of grade 1."
},

{
    "location": "projspace.html#Operations-and-Functions-1",
    "page": "Projective Space ($\\mathbb{G}_{k+1}$)",
    "title": "Operations and Functions",
    "category": "section",
    "text": "Since types do not change in this space, all operations and functions previously mentioned can be used in the same way. For example, considering p-basis elements we have\nv=pbasis(e12,false,2.0)\n\nu=pbasis(e123,true,1.0)\n\nu ∘ v\n\nu⋅v\n\nu^v\nor considering multivectors in this space,\nu=pbasis(e123,true,1.0)+pbasis(e12,true,-2.5)\n\nv=pbasis(e1,true,1.0)+pbasis(e12,true,5.0)\n\nu ∘ v\n\nu⋅v\n\nu^v\nThe function iH can also be used, which is the inverse of the H function, that is, it takes a 1-vector of $\\mathbb{G}_{k+1}$ with a component in the direction e(k+1) different from zero and projects it into the k-dimensional Euclidean subspace of $\\mathbb{G}_{k}$."
},

{
    "location": "confspace.html#",
    "page": "Conformal Space ($\\mathbb{G}_{k+1,1}$)",
    "title": "Conformal Space ($\\mathbb{G}_{k+1,1}$)",
    "category": "page",
    "text": "using Liga\nlayout(3)"
},

{
    "location": "confspace.html#Conformal-Space-(\\mathbb{G}_{k1,1})-1",
    "page": "Conformal Space ($\\mathbb{G}_{k+1,1}$)",
    "title": "Conformal Space ($\\mathbb{G}_{k+1,1}$)",
    "category": "section",
    "text": "In space $\\mathbb{G}_k$ we can use the stereographic embedding, taking an element of  $\\mathbb{G}_k$ to $\\mathbb{G}_{k+1}$ using the function S that represents the stereographic embedding, so to take this element to the conformal space is used the function Hm, which represents the homogenization in the projective space with negative signature (Minkowski space), that is, instead of a new basis element that squares to +1, another basis element is added which is the e- that squares to -1.Since these functions depend on e+ and e-, the types based on pbasis were created, however, one can use e∞ and e∘ as basis elements for conformal space which is a subspace of $\\mathbb{G}_k$ to $\\mathbb{G}_{k+1}$, then the types based on cbasis were created for this purpose."
},

{
    "location": "confspace.html#pbasis-and-cbasis-1",
    "page": "Conformal Space ($\\mathbb{G}_{k+1,1}$)",
    "title": "pbasis and cbasis",
    "category": "section",
    "text": "The main difference between these two forms is in these types, the components e∞ and e∘ are generated from e+ and e- through the following expressions.e∞ = e₋ + e₊\ne∘ = 0.5(e₋ - e₊)thus,(e∞)² = (e₋ + e₊)² = 0\n(e∘)² = 0.5(e₋ - e₊)² = 0\ne∞ ⋅ e∘ = -1Another difference is in the creation of a cbasis element, while creating an element of type pbasis requires a logical array and an additional logical element, to create an element of the cbasis type we need two more logical elements, that is,pbasis(e1, true, 1.0)creates the element1.0e1-of type pbasis andcbasis(e1, true, false, 1.0)creates the elemente1∞which is of the type cbasis. Let us try to explicit more these concepts. We can create e∞ and e∘ directly using the cbasistype.einf=cb(id, true, false, 1.0)\n\ntypeof(einf)\n\ne0=cb(id, false, true, 1.0)\n\ntypeof(e0)\n\neinf ∘ e0\n\ntypeof(einf ∘ e0)\n"
},

{
    "location": "confspace.html#Multivectors-and-blades-1",
    "page": "Conformal Space ($\\mathbb{G}_{k+1,1}$)",
    "title": "Multivectors and blades",
    "category": "section",
    "text": "The only difference between pmultvec and cmultvec as well as for pblade and cblade, is the type of basis that compose them.For multivectorspmultvec(Vector{pbasis})\ncmultvec(Vector{cbasis})and for bladespblade(Vector{pmultvec})\ncblade(Vector{cmultvec})Let us to see an example to define a c-multivector.u = cmultvec([cb(id, true, false, 1.0),cb(e12, true, false, 2.0)])\n\ntypeof(u)\n\nv = cmultvec([cb(e1, true, false, 3.0),cb(e123, true, false, -2.0)])\n\ntypeof(v)"
},

{
    "location": "confspace.html#Operations-and-functions-1",
    "page": "Conformal Space ($\\mathbb{G}_{k+1,1}$)",
    "title": "Operations and functions",
    "category": "section",
    "text": "The operations and functions for these types are the same as for the Gₖ space, only the conjugate function is added, which returns the conjugate of an element.The main difference here is that all functions depend on the types pbasis, pmultvec and pblade, that is, when elements of the cbasis, cmultvec or cblade types are operated they are automatically converted to one of the first three types and them reconverted back. Following previuos idea, let us see  cmultivectors in action.u = cmultvec([cb(id, true, false, 3.0),cb(e13, false, false, -7.0)])\n\nv = cmultvec([cb(e1, true, false, 3.0),cb(e123, true, false, -2.0)])\n\nu+v #sum\n\nu ∘ v #geometric product\n\nu⋅v #inner product\n\ntypeof(u⋅v)\n\n#and others...\n"
},

{
    "location": "advancedex.html#",
    "page": "Advanced Examples",
    "title": "Advanced Examples",
    "category": "page",
    "text": ""
},

{
    "location": "advancedex.html#Advanced-Examples-1",
    "page": "Advanced Examples",
    "title": "Advanced Examples",
    "category": "section",
    "text": "Here we are considering an example of Liga in circle detection.  This problem is an important problem with several applications in Computational Vision. Roughly speaking, this problem consists in, given a binarized image defined by  A={(x1,y1),...,(xn,yn)} points try to get a circle with center (a,b) and radious r that contain p points of A. The more classical way to solve this problem is using Hough Transform.  For a good review of this subject we recommend [1] and [2]. The idea of Hough  Transform for circle detection is use a voting array and discretizations of the  desired parameters (center and radius). Some improvements depending on coordinates changing. Using GA we can get an elegant way for solving this problem without discretizations (for example, without discretization for radius r). We do not expect to propose a faster way to solve but an elegant way. Besides Liga is not optimized for specific calculations. So, be careful. In sequence we have the code. This code can be downloaded here.#load basic packages\nusing PyPlot, Liga\nlayout(2)\n\n#Auxiliary function	\nfunction findind(tensor::Array{Float64,3})\n  	(m,n,p)=size(tensor)\n  	vlmax=maximum(tensor)\n  	for i=1:m,j=1:n,k=1:p\n          if tensor[i,j,k]==vlmax\n              return i,j,k,vlmax\n          end\n	end\nend\n\n\nfunction main(filename)\n	#reading the example file\n	A=readdlm(filename)\n	#setting the size of the image and defining main parameters\n	m=300\n	n=300\n	maxrad=150	\n	npun=length(A[:,1])	\n	Vot=zeros(m,n,maxrad)\n	#embeding A in conformal space and stored in EA\n	EA=Array{cmultvec,2}(npun,1)\n	einf=cb(id, true, false, 1.0)\n	ea12=cb(e12, false, false, 1.0)\n	eaid=cb(id, true, false, -1.0)\n	SS=cmultvec(1.0*(e1,false,true))\n	raio=1.0\n	centro=[1.0,1.0]\n	for i=1:npun\n		EA[i]=conformal([A[i,1],A[i,2]])\n	end\n	#display(EA)\n	#the main block of the idea\n	#for i=1:npun,j=i+1:npun,k=j+1:npun #conservative method\n	for i=1:2:npun,j=i+1:3:npun,k=j+1:5:npun #speed up				\n		X=EA[i]^EA[j]^EA[k]\n		SS=dual(X)\n		raio = sqrt(((SS⋅SS)/((SS⋅einf)⋅(SS⋅einf))).comp[1].scl) \n		centro = conftore(projection(SS, ea12) / (SS ⋅ eaid))\n    	if 1.0<=centro[1]<=m && 1.0<=centro[2]<=n\n			if 1.0<=abs(raio)<=maxrad\n				centro=round.(Int,centro)\n				raio=abs.(round(Int,raio))\n				Vot[centro[1],centro[2],raio]=Vot[centro[1],centro[2],raio]+1\n			end\n		end		\n	end		\n	sol=zeros(3)\n	valsol=0.0\n	sol[1],sol[2],sol[3],valsol=findind(Vot)\n	display(sol)\n	display(valsol)\n	draw_solution(A,sol[1],sol[2],sol[3])\nend\n\nfunction draw_solution(A,cfx,cfy,rf)\n	angulo=[0:((2*pi)/360):2*pi+1;]\n	x=zeros(length(angulo))\n	y=zeros(length(angulo))\n	for i=1:length(angulo)\n    	x[i]=(abs(rf)*cos(angulo[i])+cfx)\n   		y[i]=(abs(rf)*sin(angulo[i])+cfy)\n	end\n	plot(A[:,1],A[:,2],\".\")\n	plot(x,y,color=\"red\",\"-\")\n	ax=gca()\n	ax[:axis](\"equal\")\nendAs result of example datafile-50-300-300-10.txt we have(Image: figure_1.jpeg)[1] Duda, Richard O., and Peter E. Hart. \"Use of the Hough transformation to  detect lines and curves in pictures.\" Communications of the ACM 15.1 (1972): 11-15.[2] Illingworth, John, and Josef Kittler. \"A survey of the Hough transform. \"Computer vision, graphics, and image processing 44.1 (1988): 87-116."
},

{
    "location": "summary.html#",
    "page": "Summaries",
    "title": "Summaries",
    "category": "page",
    "text": ""
},

{
    "location": "summary.html#Summaries-1",
    "page": "Summaries",
    "title": "Summaries",
    "category": "section",
    "text": "Here, we summarize all types definitions, operations  and functions supported by Liga package. "
},

{
    "location": "summary.html#Liga.kbasis",
    "page": "Summaries",
    "title": "Liga.kbasis",
    "category": "type",
    "text": "c(x::Vector{Bool}, a::Number)\n\nCreates a scaled basis element of geometric algebra space from x and a.\n\njulia> kbasis(e13, 2.3)\n2.3e13\njulia> kbasis(e23)\n1.0e23\n\n\n\n"
},

{
    "location": "summary.html#Liga.kmultvec",
    "page": "Summaries",
    "title": "Liga.kmultvec",
    "category": "type",
    "text": "kmultvec(X::Vector{kbasis})\n\nCreates a multivector in the Euclidean geometric algebra (Gₚ) composed by the sum of the coordinates of X.\n\njulia> kmultvec([kbasis(e1, 2.0), kbasis(e12, 3.3), kbasis(e123)])\n2.0e1 + 3.3e12 + 1.0e123\n\n\n\n"
},

{
    "location": "summary.html#Liga.kblade",
    "page": "Summaries",
    "title": "Liga.kblade",
    "category": "type",
    "text": "kblade(X::Vector{kmultvec})\n\nCreates a blade in the Euclidean geometric algebra (Gp) composed by the outer product of the elements of X.\n\nThe blade is created only if the multivectors that compound X are both 1-vectors and X is a L.i set\n\njulia> a = kmultvec([kbasis(e1,2.0),kbasis(e3,4.0),kbasis(e3,4.0)])\n2.0e1 + 4.0e3 + 4.0e3\njulia> b = kmultvec([kbasis(e1, 1.0), kbasis(e2, 1.0), kbasis(e3,2.0)])\n1.0e1 + 1.0e2 + 2.0e3\njulia> A = kblade([a,b])\n(2.0e1 + 4.0e3 + 4.0e3)∧(1.0e1 + 1.0e2 + 2.0e3)\n\n\n\n"
},

{
    "location": "summary.html#Liga.cbasis",
    "page": "Summaries",
    "title": "Liga.cbasis",
    "category": "type",
    "text": "cbasis(a::Vector{Bool},b::Bool,c::Bool,d::Number)\n\nCreates a scaled basis element of geometric algebra space (Gp+1,1) from a, b, c and d using e∞ and e∘ as basis elements instead of the regular basis.\n\n\n\n"
},

{
    "location": "summary.html#Liga.cmultvec",
    "page": "Summaries",
    "title": "Liga.cmultvec",
    "category": "type",
    "text": "cmultvec(X::Vector{cbasis})\n\nCreates a multivector in the geometric algebra (Gp+1,1) composed by the sum of the coordinates of X.\n\n\n\n"
},

{
    "location": "summary.html#Liga.cblade",
    "page": "Summaries",
    "title": "Liga.cblade",
    "category": "type",
    "text": "cblade(X::Vector{cmultvec})\n\nCreates a blade in the geometric algebra $ \\mathbb{G}_{k+1,1} $ composed by the  outer product of the elements of X.\n\n\n\n"
},

{
    "location": "summary.html#Summary-of-types-1",
    "page": "Summaries",
    "title": "Summary of types",
    "category": "section",
    "text": "kbasis\nkmultvec\nkblade\ncbasis\ncmultvec\ncblade"
},

{
    "location": "summary.html#Liga.geoprod",
    "page": "Summaries",
    "title": "Liga.geoprod",
    "category": "function",
    "text": "geoprod(a::kbasis, b::kbasis)\ngeoprod(X::kmultvec, Y::kmultvec)\ngeoprod(A::kblade, B::kblade)\ngeoprod(a::pbasis, b::pbasis)\ngeoprod(X::pmultvec, Y::pmultvec)\ngeoprod(A::pblade, B::pblade)\ngeoprod(a::cbasis, b::cbasis)\ngeoprod(X::cmultvec, Y::cmultvec)\ngeoprod(A::cblade, B::cblade)\ngeoprod(A::Array{kbasis,2},B::Array{kbasis,1})\ngeoprod(A::Array{kbasis,2},B::Array{kbasis,2})\ngeoprod(A::Array{kmultvec,2},B::Array{kmultvec,1})\ngeoprod(A::Array{kmultvec,2},B::Array{kmultvec,2})\n\nReceives as input parameters elements of the same type and returns the geometric product between them. This function is extended to array like a traditional product. It can be used the operator \"∘ (circ)\" instead of call the function.\n\n\n\n"
},

{
    "location": "summary.html#Liga.inner",
    "page": "Summaries",
    "title": "Liga.inner",
    "category": "function",
    "text": "inner(a::kbasis, b::kbasis)\ninner(X::kmultvec, Y::kmultvec)\ninner(A::kblade, B::kblade)\ninner(a::pbasis, b::pbasis)\ninner(X::pmultvec, Y::pmultvec)\ninner(A::pblade, B::pblade)\ninner(a::cbasis, b::cbasis)\ninner(X::cmultvec, Y::cmultvec)\ninner(A::cblade, B::cblade)\n\n\nReceives as input parameters elements of the same type and returns the inner product between them.\n\nIt can be used the operator \"⋅ (cdot)\" instead of call the function.\n\n\n\n"
},

{
    "location": "summary.html#Liga.outer",
    "page": "Summaries",
    "title": "Liga.outer",
    "category": "function",
    "text": "outer(a::kbasis, b::kbasis)\nouter(X::kmultvec, Y::kmultvec)\nouter(A::kblade, B::kblade)\nouter(a::pbasis, b::pbasis)\nouter(X::pmultvec, Y::pmultvec)\nouter(A::pblade, B::pblade)\nouter(a::cbasis, b::cbasis)\nouter(X::cmultvec, Y::cmultvec)\nouter(A::cblade, B::cblade)\n\n\nReceives as input parameters elements of the same type and returns the outer product between them.\n\nIt can be used the operator \" ^ \" instead of call the function.\n\n\n\n"
},

{
    "location": "summary.html#Liga.scalar",
    "page": "Summaries",
    "title": "Liga.scalar",
    "category": "function",
    "text": "scalar(a::kbasis, b::kbasis)\nscalar(X::kmultvec, Y::kmultvec)\nscalar(A::kblade, B::kblade)\nscalar(a::pbasis, b::pbasis)\nscalar(X::pmultvec, Y::pmultvec)\nscalar(A::pblade, B::pblade)\nscalar(a::cbasis, b::cbasis)\nscalar(X::cmultvec, Y::cmultvec)\nscalar(A::cblade, B::cblade)\n\n\nReceives as input parameters elements of the same type and returns the scalar product between them.\n\nIt can be used the operator \" * \" instead of call the function.\n\n\n\n"
},

{
    "location": "summary.html#Liga.dual",
    "page": "Summaries",
    "title": "Liga.dual",
    "category": "function",
    "text": "dual(a::kbasis)\ndual(X::kmultvec)\ndual(A::kblade)\ndual(X::pmultvec)\ndual(a::pbasis)\ndual(A::pblade)\n\n\nReturns the dual element of the basis blade a, the multivetor X or the blade A.\n\n\n\n"
},

{
    "location": "summary.html#Liga.grade",
    "page": "Summaries",
    "title": "Liga.grade",
    "category": "function",
    "text": "grade(a::kbasis)\ngrade(a::pbasis)\ngrade(a::cbasis)\n\nReturns the grade of the basis element x\n\njulia> grade(kbasis(e13, 1.5))\n2\n\n\n\n"
},

{
    "location": "summary.html#Liga.inverse",
    "page": "Summaries",
    "title": "Liga.inverse",
    "category": "function",
    "text": "inverse(a::kmultvec)\n\nReturns the inverse of a multi vector input. Be  careful with this function because in current version it isn\'t optimized!\n\n\n\ninverse(A::kblade)\ninverse(A::pblade)\ninverse(A::cblade)\n\n\nReturns the inverse of A if A is a non-null-blade or the pseudoinverse A if A is a null-blade.\n\n\n\n"
},

{
    "location": "summary.html#Liga.mvsum",
    "page": "Summaries",
    "title": "Liga.mvsum",
    "category": "function",
    "text": "mvsum(a::kbasis, b::kbasis)\nmvsum(X::kMulvec, Y::kmultvec)\nmvsum(a::pbasis, b::pbasis)\nmvsum(X::pMulvec, Y::pmultvec)\nmvsum(a::cbasis, b::cbasis)\nmvsum(X::cMulvec, Y::cmultvec)\n\n\nReceives as input parameters elements of the same type and returns the sum of them.\n\nIt can be used the operators \"+\" or \"-\" instead of call the function.\n\n\n\n"
},

{
    "location": "summary.html#Liga.mvreverse",
    "page": "Summaries",
    "title": "Liga.mvreverse",
    "category": "function",
    "text": "mvreverse(a::kbasis)\nmvreverse(X::kmultvec)\nmvreverse(A::kblade)\nmvreverse(a::pbasis)\nmvreverse(A::pblade)\nmvreverse(a::cbasis)\nmvreverse(A::cblade)\n\n\nReturns the reverse of the element.\n\n\n\n"
},

{
    "location": "summary.html#Liga.magnitude",
    "page": "Summaries",
    "title": "Liga.magnitude",
    "category": "function",
    "text": "magnitude(a::kbasis)\nmagnitude(a::kmultvec)\nmagnitude(A::kblade)\nmagnitude(a::pbasis)\nmagnitude(a::pmultvec)\nmagnitude(A::pblade)\nmagnitude(a::cbasis)\nmagnitude(A::cblade)\n\n\nReturns the magnitude of the input element.\n\n\n\n"
},

{
    "location": "summary.html#Liga.projection",
    "page": "Summaries",
    "title": "Liga.projection",
    "category": "function",
    "text": "projection(A::kblade, N::kblade)\nprojection(A::pblade, N::pblade)\n\n\nReturns the projection of the blade A onto the blade N.\n\n\n\n"
},

{
    "location": "summary.html#Liga.rejection",
    "page": "Summaries",
    "title": "Liga.rejection",
    "category": "function",
    "text": "rejection(A::kblade, N::kblade)\nrejection(A::pblade, N::pblade)\n\n\nReturns the rejection of the blade A from the blade N.\n\n\n\n"
},

{
    "location": "summary.html#Liga.S",
    "page": "Summaries",
    "title": "Liga.S",
    "category": "function",
    "text": "S(x::Vector{Float64})\nS(X::kmultvec)\n\nThe stereographic embedding of Euclidean space.\n\nTakes a vector of Euclidean space in its stereographic embedding.\n\n\n\n"
},

{
    "location": "summary.html#Liga.Hm",
    "page": "Summaries",
    "title": "Liga.Hm",
    "category": "function",
    "text": "Hm(x::Vector{Float64})\n\nAuxiliary function.\n\n\n\n"
},

{
    "location": "summary.html#Liga.H",
    "page": "Summaries",
    "title": "Liga.H",
    "category": "function",
    "text": "H(X::Vector{Float64})\nH(X::kmultvec)\n\nHomegenization. It takes a vector in mathbbR^n into its embedding in the affine space.\n\n\n\n"
},

{
    "location": "summary.html#Liga.iH",
    "page": "Summaries",
    "title": "Liga.iH",
    "category": "function",
    "text": "iH(X::kmultvec)\n\nInversion for homonezation function H.\n\n\n\n"
},

{
    "location": "summary.html#Liga.pconformal",
    "page": "Summaries",
    "title": "Liga.pconformal",
    "category": "function",
    "text": "pconformal(x::Vector{Float64})\n\nConformal embedding of a Euclidean vector into a vector of Gp+1,1 with new basis elements being e+ and e-.\n\n\n\n"
},

{
    "location": "summary.html#Liga.ipconformal",
    "page": "Summaries",
    "title": "Liga.ipconformal",
    "category": "function",
    "text": "ipconformal(X::pmultvec)\n\nThe inverse of \"pconformal\" function.\n\n\n\n"
},

{
    "location": "summary.html#Liga.conformal",
    "page": "Summaries",
    "title": "Liga.conformal",
    "category": "function",
    "text": "conformal(x::Vector{Float64})\n\nConformal embedding of a Euclidean vector into a vector of Gp+1,1 with new basis elements being e∞ and e∘.\n\n\n\n"
},

{
    "location": "summary.html#Liga.iconformal",
    "page": "Summaries",
    "title": "Liga.iconformal",
    "category": "function",
    "text": "iconformal(X::cmultvec)\n\nThe inverse of \"conformal\" function.\n\n\n\n"
},

{
    "location": "summary.html#Summary-of-functions-1",
    "page": "Summaries",
    "title": "Summary of functions",
    "category": "section",
    "text": "geoprod\ninner\nouter\nscalar\ndual\ngrade\ninverse\nmvsum\nmvreverse\nmagnitude\nprojection\nrejection\nS\nHm\nH\niH\npconformal\nipconformal\nconformal\niconformal"
},

]}
