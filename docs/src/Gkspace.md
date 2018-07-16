```@setup liga
using Liga
layout(3)
```

# \$\mathbb{G}_k\$ Space


The Euclidean geometric algebra \$G_k\$ space is the basic space of this project. The main operations of the space, that is, the *Clifford* or *geometric product*, the *inner*, *outer* and *scalar* products, the *reverse*, *dual*, *projection* and *rejection* of elements of this space are programmed and ready for use.

To create an element of this space is necessary to call the **types** which have only elements of type ```kbasis``` in their composition, that is, the elements of the *type* ```kbasis``` itself, the ```kmultvec``` and the ```kblade``` *types*. Functions and operations related to this space 
can be used in projetive and conformal space, but we explain in [Projetive Space Section](projspace.md) this.

**If the reader has not yet read about the types supported by this library, reading in [Section Types](types.md) is highly recommended.**


## Functions and Operations

Here we will give an explanation of how each *function* works, in terms of what is taken as the *input* parameter and what the *function* will return, so that some detailed examples are given later in order to simplify the understanding of the program and avoid possible mistakes.

### Operations

#### Geometric Product - ```geoprod(a,b)``` or ```a ∘ b```

As in the literature, the *Clifford* or *geometric product* is the main operation of the space, and is the basis for most of the other operations. It takes two elements of the *type* ```kbasis```, ```kmultvec``` or ```kblade``` and return their *geometric product*.

```
geoprod(a::kbasis, b::kbasis)
geoprod(X::kmultvec, Y::kmultvec)
geoprod(A::kblade, B::kblade)
```

To simplify the operation *input* it can be used the base operator ```∘``` instead of call the *function* ```geoprod```.

     Base.:∘(a::kbasis, b::kbasis)
     Base.:∘(X::kmultvec, Y::kmultvec)
     Base.:∘(A::kblade, B::kblade)

In order to obtain the *geometric product* of elements of different *types*,as long as the operator is used, it is not necessary to convert them before . As explained above, since a generic element of *geometric algebra* is a *multivector*, the elements of type ```kbasis``` and ```kblade``` can be converted to the *type* ```kmultvec``` and thus can be apply the *geometric product* operation. These conversions are done automatically when two parameters of different *types* are used as input.

See **Examples** section below to better understand the *input* and *output* of operations.

#### Inner Product - ```inner(a,b)``` or ```a ⋅ b```

The *input* parameters for the *inner product* are of the same *types* as for the *geometric product*, so to obtain the *inner product* between two elements of *type* ```kbasis```, ```kmultvec``` or ```kblade``` just use the ```inner``` *function*.

     inner(a::kbasis, b::kbasis)
     inner(X::kmultvec, Y::kmultvec)
     inner(A::kblade, B::kblade)

As for the *geometric product*, to obtain the *inner product* is sufficient to use a base operator, in this case, the operator ```⋅``` instead of call the ```inner``` *function*.

     Base.:⋅(a::kbasis, b::kbasis)
     Base.:⋅(X::kmultvec, Y::kmultvec)
     Base.:⋅(A::kblade, B::kblade)


#### Outer Product - ```outer(a,b)``` or ```a ^ b```

As for the *geometric* and *inner* products, one can use the *function* ```outer``` to obtain the *outer product* between two elements of *type* ```kbasis```, ```kmultvec``` or ```kblade```.

     outer(a::kbasis, b::kbasis)
     outer(X::kmultvec, Y::kmultvec)
     outer(A::kblade, B::kblade)

The base operator that can be used instead of call the *function* in this case is the ```^``` symbol.

     Base.:^(a::kbasis, b::kbasis)
     Base.:^(X::kmultvec, Y::kmultvec)
     Base.:^(A::kblade, B::kblade)


#### Scalar Product - ```scalar(a,b)``` or ```a * b```

The ```scalar``` *function* is not different from the other three mentioned above, in order to obtain the *scalar product* between two elements of *type* ```kbasis```, ```kmultvec``` or ```kblade``` just use the *function* ```inner```.

     scalar(a::kbasis, b::kbasis)
     scalar(X::kmultvec, Y::kmultvec)
     scalar(A::kblade, B::kblade)

Or using the base operator ```*```.

     Base.:*(a::kbasis, b::kbasis)
     Base.:*(X::kmultvec, Y::kmultvec)
     Base.:* (A::kblade, B::kblade)

**In this operation, unlike the previous ones, the output parameter will always be of *type* ```Number```.**

#### Sum and difference - ```mvsum(a,b)``` or ```a + b``` and ```a - b```

This *function* has a particularity in the *type* of *output*, since the *sum* of two elements of *type* ```kbasis``` will only be of *type* ```kbasis``` if they are l.d, otherwise the output will be of *type* ```kmultvec```.

To obtain the *sum* of two elements of *types* ```kbasis```,  ```kmultvec``` or ```kblade```, one can call the *function* ```mvsum```, as well as the previous functions,

     mvsum(a::kbasis, b::kbasis)
     mvsum(X::kmultvec, Y::kmultvec)
     mvsum(A::kblade, B::kblade)

or using the base operator ```+```

     Base.:+(a::kbasis, b::kbasis)
     Base.:+(X::kmultvec, Y::kmultvec)
     Base.:+(A::kblade, B::kblade)

but, to make it easier, in addition to the base operator ```+```, it can be used the base operator ```-``` which convert the second element to its *additive inverse*.

     Base.:-(a::kbasis, b::kbasis)
     Base.:-(X::kmultvec, Y::kmultvec)
     Base.:-(A::kblade, B::kblade)

thus obtaining the *difference* between the two elements.

### Functions

#### Grade - ```grade(a)```

The *grade function* take as *input* an element of the *type* ```kbasis``` and return its grade which is an element of the *type* ```Int```, it's used more as an auxiliary *function*.

     grade(a::kbasis)


#### Dual - ```dual(a)```

The *dual function*, as the name suggests, returns its dual element, taking as *input* elements of *types* ```kbasis```,  ```kmultvec``` or ```kblade```.

     dual(a::kbasis)
     dual(X::kmultvec)
     dual(A::kblade)

As for the *products* and *sum*, in case the *input* of the ```kblade``` the *output* will be of the *type* ```kmultvec```.


#### Reverse - ```mvreverse(a)```

Returns the *reverse* of the *input* element, but, unlike other *functions*, this will always return an element of the same *input type*, being it ```kbasis```,  ```kmultvec``` or ```kblade```.

     mvreverse(a::kbasis)
     mvreverse(X::kmultvec)
     mvreverse(A::kblade)

The conjugation *function* isn't necessary in \$G_k\$ space, since its depends of negative grade.

#### Magnitude

The *magnitude function* receives as *input* an element of type ```kbasis```,  ```kmultvec``` or ```kblade``` and returns its magnitude of the *type* ```Ǹumber```.

     magnitude(a::kbasis)
     magnitude(X::kmultvec)
     magnitude(A::kblade)

### Special functions

In this space, some *functions* will work only with elements of the *type* ```kblade```, one being the *inversion function* and the other two depending on the first one, since not every multivector of Gₖ has inverse with respect to the *geometric product*.

#### Inversion

This *function* returns the *inverse* element of the *input* with respect to the *geometric product*.

     inverse(A::kblade)

The *output type* will also be ```kblade```.

#### Projection

The *projection* is a simple *function* that depends on two elements of the *type* ```kblade```

     projection(A::kblade, N::kblade)

and it returns, in this case, the *projection* of ```A``` onto ```N```, which will be an element of the *type* ```kmultvec```.


#### Rejection

The *rejection function* also receives two elements of the *type* ```kblade``` as input

     rejection(A::kblade, N::kblade)

and, as for the *projection*, it returns an element of the *type* ```kmultvec``` which will be the *rejection* of the ```kblade``` ```A``` from the ```kblade```  ```N```.


## Examples

All of the following examples are based on ```layout(3)```, that is, the space generated will be \$\mathbb{G}_3\$ depending on the objects ```id```, ```e1```, ```e2```, ```e3```, ```e12```, ```e13```, ```e23``` and ```e123```, however everything follows in the same way for other dimensions.

#### Grade

The grade *function* is the most basic of them, as mentioned above, it's used more like an auxiliary *function* for the other and it returns the grade of a ```kbasis``` element as follows:

```@repl liga
grade(kbasis(e12, 1.0))
```

which is the grade of ```e12```.

Another way to do this is by calling the element before the *function*

```@repl liga

a = kbasis(e123, -2.3)

```
and then use the grade *function*

```@repl liga
grade(a)
```
which is the grade of ```-2.3e123```.

#### Geometric, inner, outer and scalar products

The *geometric, inner and outer products* in case the *input* elements are of *type* ```kd``` or ```kmultvec```, will always return the same type, even if the result is a *scalar* or a *basis blade*, that is, if the following elements are called:

```@repl liga

a = kbasis(e3, 2.0)

b = kbasis(e12, -1.0)

c = kbasis(id, 0.0)

d = kbasis(e123)
```
when operated, for example

```@repl liga
geoprod(a, b)
```
which is a ```kbasis``` element, although when operated

```@repl liga
inner(a, b)
```
this element is seen as an object of *type* ```kbasis```.
Note that if a second parameter is not placed when called the ```kbasis``` type it will be generated with the scalar 1.0.

As mentioned before, we can use the base operators ```∘```, ```⋅``` and ```^``` to calculate these products instead of call the functions

```@repl
a ^ c

b ⋅ d

```

The same is true if the input is of the type kmultvec as follows:

```@repl liga
X = kmultvec([kbasis(e12,2.0),kbasis(e123,-1.5),kbasis(e3,4.0)])

Y = kmultvec([a, b, kbasis(id, 2.0)])

Z = kmultvec(d)
```

and these can be operated in the same way as the previous ones

```@repl liga
X ∘ Y
Y ⋅ Z
Y ^ Z
```

Now, in the case of ```kblade``` *type*, everything follows like the previous cases, the difference is in the *output* of the *function*. For example, if we call the following ```kblade``` elements

```@repl liga
A = kblade([kmultvec([kb(e1)]), kmultvec([kb(e2)]), kmultvec([kb(e3)])])
B = kblade([kmultvec([kb(e1, 2.0), a]), kmultvec([kb(e1, 3.0)])])
C = kblade([kmultvec(a)])

A ⋅ B

```
returns an element of the *type*  ```kmultvec```. In fact, type

```@repl liga

typeof(A ⋅ B)
```

the same will happen for the other operations

     B ^ C

returns

     0.0

and

     B ∘ C

returns

     -12.0e1

The *scalar product* follows a different rule, in all cases the *output* will always be of *Number*. For example, if we want to calculate the *scalar product* between the ```kbasis``` element ```d``` and the ```kblade``` ```A``` previously called simply do as follows:

```@repl liga
d=kbasis(id,-2.0)
d * A
```
remember that this is an element of *type* ```Number```.
```@repl liga
typeof(ans)
```

You can operate any other elements of *types* ```kbasis```, ```kmultvec``` and ```kblade```,  
```@repl liga
 a * a
 b * X
 B * B
 X * A
```


#### Sum and difference

The *mvsum* operation works in the same way as the previous ones with respect to the *input* parameters, the difference is in the *output* parameters, which in the case of the *sum* of elements of the *type* ```kbasis``` the *output* there are two possible *outputs*

```@repl liga
mvsum(kbasis(e23, 2.3), kbasis(e23, 2.7))
```
in this case, it may be noted that the parameter ```e23``` is used in the construction of both *inputs*, then the it returns

     5.0e23

which is of the *type* ```kbasis```,
```@repl liga
typeof(ans)
```

However, if we calculate the following sum

     a + b

it returns

     2.0e3 + -e12

which is of the *type* ```kmultvec```, remembering that ```a``` and ```b``` are the elements ```kbasis(e3, 2.0)``` and ```kbasis(e12, -1.0)``` called above.

Other than that, the rest follows as before

```@repl liga

A + B

typeof(ans)

X - C

typeof(ans)

d + Y

typeof(ans)

a + d

typeof(ans)

```     


Note that here we can use either the base operator ```+``` or ```-```.

Again, the *output* when applied the *sum* and the *difference*  between elements of the *type* ```kblade``` will be of the *type* ```kmultvec```.

#### Dual

The *dual function* returns the dual of the *input* element as follows

```@repl liga
dual(a)

typeof(dual(a))

dual(X)

typeof(dual(X))

dual(A)

typeof(dual(A))
```
As for the previous operations, except for the *scalar product*, in case the *input* is of *type* ```kblade``` the *output* will be of *type* ```kmultvec```.

#### Reverse

Unlike other *functions* and *operations* the *mvreverse* function returns an element of the same *type* as the *input*, be it ```kbasis```, ```kmultvec``` or ```kblade```

```@repl liga
mvreverse(a)

typeof(mvreverse(a))

mvreverse(X)

typeof(mvreverse(X))

mvreverse(A)

typeof(mvreverse(A))
```
#### Magnitude

The *magnitude* in \$\mathbb{G}_k\$ is nothing more than the square root of the *scalar product* of an element by its reverse, so it is, like the scalar product, always returns an element of *type* ```Float```. For example,

```@repl liga
magnitude(a)

typeof(magnitude(a))

magnitude(X)

typeof(magnitude(X))

magnitude(A)

typeof(magnitude(A))
```


#### Inversion

In current version of Liga, this important functions works just for ```kblade``` and 
```kmultvec``` types. We are improve to more types. The return is an equivalent type to 
input.

``` @repl liga

inverse(A)

typeof(inverse(A))

inverse(B)

typeof(inverse(B))

inverse(C)

typeof(inverse(C))

```

it returns the inverse elements of ```A```,```B``` and ```C```.

We can test the function using the *geometric product* of a ```kblade``` by its *inverse*, and it should return 1.0

```@repl liga

A ∘ inverse(A)

B ∘ inverse(B)

C ∘ inverse(C)

```

Note that as mentioned when dealing with the *geometric product*, this is a element of *type* ```kmultvec```

```@repl liga
typeof(B ∘ inverse(B))
```

Let us see an example related to *k-multivector*. 

```@repl liga

inverse(X)

typeof(inverse(X))

X ∘ inverse(X)

typeof(X ∘ inverse(X))

```

Note that, in ```X ∘ inverse(X)``` the result was 2.7755575615628914e-17e12-1.3877787807814457e-17e3+1.0
wich is approximated 1.0 (and ```kmultvec```);

#### Projection

This function depends on the *inversion* and by now can be used only with elements of *type* ```kblades``` as *input*. For example,

     projection(A, B)
     projection(B, A)
     projection(A, C)

will return

     e123
     6.0e13
     e123

these *output* elements are of *type* ```kmultvec```

#### Rejection

Like the above function **Projection**, the rejection receives an element of *type* ```kblade``` and returns an element of *type* ```kmultvec``` which is the *rejection* of the first *input* parameter from the second.

     rejection(B, kblade([kmultvec(kbasis(e2))]))

will return

     6.0e13

which is the difference between the ```kblade``` ```B``` and the projection ```rejection(B, kblade([kmultvec(kbasis(e2))]))```.


## Matrix operations

The main operations can be extended to array like the classic matrix product. For example, considering the geometric product

```@repl liga
u=[1.0*e12;2.0*e123]

transpose(u)∘u

A=[1.0*e1+2.0*e123 2.0*e1+1.0*id; 0.5*e123+1.0*e3 1.0*id+1.0*e1]

B=[5.0*e12+2.0*e123 2.0*e1-1.0*id; 1.0*e123-2.0*e2 1.0*id+1.0*e1]

A∘B
```

Note in this example, additionally, a transpose function for array defined by kbasis or kmultvec types was implemented. Be carefull to use this function in keeping the same type, that is, the elements can be kmultvec or kbasis, but not both.

```@docs 
transpose(A::Array{kbasis,1})
```