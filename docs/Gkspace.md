---
layout: page
title: G(k) Space
navigation: 4
---

## Basics

The Euclidean geometric algebra Gₖ space is the basic space of this project. The main operations of the space, that is, the _Clifford_ or _geometric_ product, the _inner_, _outer_ and _scalar_ products, the _reverse_, _dual_, _projection_ and _rejection_ of elements of this space are programmed and ready for use.

To create an element of this space is necessary to call the ___types___ which have only elements of type ```kb``` in their composition, that is, the elements of the _type_ ```kb``` itself, the ```kMultvec``` and the ```kblade``` _types_.

For now, only elements of the same _type_ can be operated, however, it is possible to convert elements of the _type_ ```kb``` and ```kblade``` into elements of the _type_ ```kMultvec``` and thus apply the operations that will be described, these operations and necessary conversions will be explained below.


## Functions and Operations

Here we will give an explanation of how each _function_ works, in terms of what is taken as the _input_ parameter and what the _function_ will return, so that some detailed examples are given later in order to simplify the understanding of the program and avoid possible mistakes.

### Conversion functions

#### Convert ```kb``` to ```kMultvec```

The ```kMultvec``` _function_ will convert an element of the _type_ ```kb``` into an element of the _type_ ```kMultvec```. It can be called as follows

     kMultvec(a::kb)

and this will return an element of the _type_ ```kMultvec```.

#### Convert ```kblade``` to ```kMultvec```

The ```bltomv``` _function_ will convert an element of the _type_ ```kblade``` into an element of the _type_ ```kMultvec``` simply by calculating the outer product of the elements that compose it. It can be called as follows

     bltomv(A::kblade)

and this will also return an element of the _type_ ```kMultvec```.

### Operations

#### Grade - ```grade(a)```
The _grade_ _function_ take as _input_ an element of the _type_ ```kb``` and return its grade which is an element of the _type_ ```Int```, it's used more as an auxiliary _function_.

     grade(a::kb)

#### Geometric Product - ```geoprod(a,b)``` or ```a ∘ b```

Before start the explanation of the following operations it's important to note that the _output_ parameters of these _functions_ are of the same _type_ as _input_, except in the case of ```kblade```, which in both case will return an element of _type_ ```kMultvec```.

As in the literature, the _Clifford_ or _geometric product_ is the main operation of the space, and is the basis for most of the other operations. It takes two elements of the _type_ ```kb```, ```kMultvec``` or ```kblade``` and return their _geometric product_.

     geoprod(a::kb, b::kb)
     geoprod(X::kMultvec, Y::kMultvec)
     geoprod(A::kblade, B::kblade)

To simplify the operation _input_ it can be used the base operator ```∘``` instead of call the _function_ ```geoprod```.

     Base.:∘(a::kb, b::kb)
     Base.:∘(X::kMultvec, Y::kMultvec)
     Base.:∘(A::kblade, B::kblade)

In order to obtain the _geometric product_ of elements of different _types_ it is __necessary__ to convert them before. As explained above, since a generic element of _geometric algebra_ is a _multivector_, the elements of type ```kb``` and ```kblade``` can be converted to the _type_ ```kMultvec``` and thus can be apply the _geometric product_ operation.

See __Examples__ section below to better understand the _input_ of operations.

#### Inner Product - ```inner(a,b)``` or ```a ⋅ b```

The _input_ parameters for the _inner product_ are of the same _types_ as for the _geometric product_, so to obtain the _inner product_ between two elements of _type_ ```kb```, ```kMultvec``` or ```kblade``` just use the ```inner``` _function_.

     inner(a::kb, b::kb)
     inner(X::kMultvec, Y::kMultvec)
     inner(A::kblade, B::kblade)

As for the _geometric product_, to obtain the _inner product_ is sufficient to use a base operator, in this case, the operator ```⋅``` instead of call the ```inner``` _function_.

     Base.:⋅(a::kb, b::kb)
     Base.:⋅(X::kMultvec, Y::kMultvec)
     Base.:⋅(A::kblade, B::kblade)

The same rule is applied to operate elements of different _types_ i.e. it's necessary to convert the elements with the _conversion functions_ before operate them.

#### Outer Product - ```outer(a,b)``` or ```a ^ b```

As for the _geometric_ and _inner_ products, one can use the _function_ ```outer``` to obtain the _outer product_ between two elements of _type_ ```kb```, ```kMultvec``` or ```kblade```.

     outer(a::kb, b::kb)
     outer(X::kMultvec, Y::kMultvec)
     outer(A::kblade, B::kblade)

The base operator that can be used instead of call the _function_ in this case is the ```^``` symbol.

     Base.:^(a::kb, b::kb)
     Base.:^(X::kMultvec, Y::kMultvec)
     Base.:^(A::kblade, B::kblade)

The conversion rules are also applied to the _outer product_.

#### Scalar Product - ```scalar(a,b)``` or ```a * b```

The ```scalar``` _function_ is not different from the other three mentioned above, in order to obtain the _scalar product_ between two elements of _type_ ```kb```, ```kMultvec``` or ```kblade``` just use the _function_ ```inner```.

     outer(a::kb, b::kb)
     outer(X::kMultvec, Y::kMultvec)
     outer(A::kblade, B::kblade)

Or using the base operator ```*```.

     Base.:*(a::kb, b::kb)
     Base.:*(X::kMultvec, Y::kMultvec)
     Base.:* (A::kblade, B::kblade)

__In this operation, unlike the previous ones, the output parameter will always be of _type_ ```Number```.__

Again, the conversion rules are also applied to the _scalar product_.

#### Sum and difference - ```mvsum(a,b)``` or ```a + b``` and ```a - b```

This _function_ has a particularity in the _type_ of _output_, since the _sum_ of two elements of _type_ ```kb``` will only be of _type_ ```kb``` if they are l.d, otherwise the output will be of _type_ ```kMultvec```.

To obtain the _sum_ of two elements of _types_ ```kb```,  ```kMultvec``` or ```kblade```, one can call the _function_ ```mvsum```, as well as the previous functions,

     mvsum(a::kb, b::kb)
     mvsum(X::kMultvec, Y::kMultvec)
     mvsum(A::kblade, B::kblade)

or using the base operator ```+```

     Base.:+(a::kb, b::kb)
     Base.:+(X::kMultvec, Y::kMultvec)
     Base.:+(A::kblade, B::kblade)

but, to make it easier, in addition to the base operator ```+```, it can be used the base operator ```-``` which convert the second element to its _additive inverse_.

     Base.:-(a::kb, b::kb)
     Base.:-(X::kMultvec, Y::kMultvec)
     Base.:-(A::kblade, B::kblade)

thus obtaining the _difference_ between the two elements.

### Other functions

#### Dual - ```dual(a)```

The _dual_ _function_, as the name suggests, returns its dual element, taking as _input_ elements of _types_ ```kb```,  ```kMultvec``` or ```kblade```.

     dual(a::kb)
     dual(X::kMultvec)
     dual(A::kblade)

As for the _products_ and _sum_, in case the _input_ of the ```kblade``` the _output_ will be of the _type_ ```kMultvec```.


#### Reverse - ```mvreverse(a)```

Returns the _reverse_ of the _input_ element, but, unlike other _functions_, this will always return an element of the same _input type_, being it ```kb```,  ```kMultvec``` or ```kblade```.

     mvreverse(a::kb)
     mvreverse(X::kMultvec)
     mvreverse(A::kblade)

The conjugation _function_ isn't necessary in Gₖ space, since its depends of negative grade.

#### Magnitude

The _magnitude function_ receives as _input_ an element of type ```kb```,  ```kMultvec``` or ```kblade``` and returns its magnitude of the _type_ ```Ǹumber```.

     magnitude(a::kb)
     magnitude(X::kMultvec)
     magnitude(A::kblade)

### Special functions

In this space, some _functions_ will work only with elements of the _type_ ```kblade```, one being the _inversion function_ and the other two depending on the first one, since not every multivector of Gₖ has inverse with respect to the _geometric product_.

#### Inversion

This _function_ returns the _inverse_ element of the _input_ with respect to the _geometric product_.

     inverse(A::kblade)

The _output type_ will also be ```kblade```.

#### Projection

The _projection_ is a simple _function_ that depends on two elements of the _type_ ```kblade```

     projection(A::kblade, N::kblade)

and it returns, in this case, the _projection_ of ```A``` onto ```N```, which will be an element of the _type_ ```kMultvec```.


#### Rejection

The _rejection function_ also receives two elements of the _type_ ```kblade``` as input

     rejection(A::kblade, N::kblade)

and, as for the _projection_, it returns an element of the _type_ ```kMultvec``` which will be the _rejection_ of the ```kblade``` ```A``` from the ```kblade```  ```N```.


## Examples

All of the following examples are based on ```layout(3)```, that is, the space generated will be G₃ depending on the objects ```id```, ```e1```, ```e2```, ```e3```, ```e12```, ```e13```, ```e23``` and ```e123```, however everything follows in the same way for other dimensions.

#### Grade

The grade _function_ is the most basic of them, as mentioned above, it's used more like an auxiliary _function_ for the other and it returns the grade of a ```kb``` element as follows:

     grade(kb(e12, 1.0))

it returns

     2

which is the grade of ```e12```.

Another way to do this is by calling the element before the _function_

     a = kb(e123, -2.3)

and then use the grade _function_

     grade(a)

and it returns

     3

which is the grade of ```-2.3e123```.

#### Geometric, inner, outer and scalar products

The _geometric, inner and outer products_ in case the _input_ elements are of _type_ ```kd``` or ```kMultvec```, will always return the same type, even if the result is a _scalar_ or a _basis blade_, that is, if the following elements are called:

     a = kb(e3, 2.0)
     b = kb(e12, -1.0)
     c = kb(id, 0.0)
     d = kb(e123)

when operated, for example

     geoprod(a, b)

it returns

     -2.0e123

which is a ```kb``` element, although when operated

     inner(a, b)

it returns

     0.0

this element is seen as an object of _type_ ```kb```.
Note that if a second parameter is not placed when called the ```kb``` type it will be generated with the scalar 1.0.

As mentioned before, we can use the base operators ```∘```, ```⋅``` and ```^``` to calculate these products instead of call the functions

     a ^ c

will return

     0.0

and

     b ⋅ d

will return

     e3

The same is true if the input is of the type kMultvec as follows:

     X = kMultvec([kb(e12,2.0),kb(e123,-1.5),kb(e3,4.0)])
     Y = kMultvec([a, b, kb(id, 2.0)])
     Z = kMultvec(d)

which will generate the respective elements

     2.0e12 + -1.5e123 + 4.0e3
     2.0e3 + -1.0e12 + 2.0
     e123

and these can be operated in the same way as the previous ones

     X ∘ Y
     Y ⋅ Z
     Y ^ Z

will return respectively

     -3.0e123 + e12 + 10.0 + 6.5e3
     2.0e12 + e3 + 2.0e123
     2.0e123

Note that the ```Z```  has received the previously called ```e123``` converted from ```kb``` to ```kMulvec```.

Now, in the case of ```kblade``` _type_, everything follows like the previous cases, the difference is in the _output_ of the _function_. For example, if we call the following ```kblade``` elements

     A = kblade([kMultvec([kb(e1)]), kMultvec([kb(e2)]), kMultvec([kb(e3)])])
     B = kblade([kMultvec([kb(e1, 2.0), a]), kMultvec([kb(e1, 3.0)])])
     C = kblade([kMultvec(a)])

which created the respective elements

     (e1)∧(e2)∧(e3)
     (2.0e1 + 2.0e3)∧(3.0e1)
     (2.0e3)

and then

     A ⋅ B

returns an element of the _type_  ```kMultvec```, which in this case is

     -6.0e2

the same will happen for the other operations

     B ^ C

returns

     0.0

and

     B ∘ C

returns

     -12.0e1

The _scalar product_ follows a different rule, in all cases the _output_ will always be of _Number_. For example, if we want to calculate the _scalar product_ between the ```kb``` element ```d``` and the ```kblade``` ```A``` previously called can be used the conversion functions mentioned above as follows:

     kMultvec(d) * bltomv(A)

and it returns

     -1.00

remember that this is an element of _type_ ```Number```.

You can operate any other elements of _types_ ```kb```, ```kMultvec``` and ```kblade```,  provided that the conversion rules are respected

     a * a
     kMultvec(b) * X
     B * B
     X * bltomv(A)

it will return respectively

     4.00
     2.00
     -36.0
     1.50

___Remember, before performing a sequence of operations, make sure that all output elements are compatible for the input of the next operation.___

#### Sum and difference

The _mvsum_ operation works in the same way as the previous ones with respect to the _input_ parameters, the difference is in the _output_ parameters, which in the case of the _sum_ of elements of the _type_ ```kb``` the _output_ there are two possible _outputs_

     mvsum(kb(e23, 2.3), kb(e23, 2.7))

in this case, it may be noted that the parameter ```e23``` is used in the construction of both _inputs_, then the it returns

     5.0e23

which is of the _type_ ```kb```, however, if we calculate the following sum

     a + b

it returns

     2.0e3 + -e12

which is of the _type_ ```kMultvec```, remembering that ```a``` and ```b``` are the elements ```kb(e3, 2.0)``` and ```kb(e12, -1.0)``` called above.

Other than that, the rest follows as before

     A + B
     X - bltomv(C)
     kMultvec(d) + Y
     a + d

will return the respective _output_

     e123 + -6.0e13
     2.0e12 + -1.5e123 + 2.0e3
     e123 + -2.0e3 + -e12 + 2.0
     -2.0e3 + e123

Note that here we can use either the base operator ```+``` or ```-```.

Again, the _output_ when applied the _sum_ and the _difference_  between elements of the _type_ ```kblade``` will be of the _type_ ```kMultvec```.

#### Dual

The _dual function_ returns the dual of the _input_ element as follows

     dual(a)
     dual(X)
     dual(A)

will return

     2.0e12
     2.0e3 + -1.5 + -4.0e12
     1.0

As for the previous operations, except for the _scalar product_, in case the _input_ is of _type_ ```kblade``` the _output_ will be of _type_ ```kMultvec```.

#### Reverse

Unlike other _functions_ and _operations_ the _mvreverse_ function returns an element of the same _type_ as the _input_, be it ```kb```, ```kMultvec``` or ```kblade```

     mvreverse(a)
     mvreverse(X)
     mvreverse(A)

will return

     -2.0e3
     -2.0e12 + 1.5e123 + 4.0e3
     (e3)∧(e2)∧(e1)

respectively.

#### Magnitude

The _magnitude_ in Gₖ is nothing more than the square root of the _scalar product_ of an element by its reverse, so it is, like the scalar product, always returns an element of _type_ ```Number```. For example,

     magnitude(a)
     magnitude(X)
     magnitude(A)

returns the elements

     2.00
     4.71…
     1.00

which are both of _type_ ```Number```.

#### Inversion

This is one of the three special functions, since it can only be used with ```kblade``` as input parameter,

     inverse(A)
     inverse(B)
     inverse(C)

it returns the inverse elements of ```A```,```B``` and ```C```,

     (e3)∧(e2)∧(e1)
     (0.08333333333333333e1)∧(2.0e1 + -2.0e3)
     (-0.5e3)

we can test the function using the _geometric product_ of a ```kblade``` by its _inverse_, and it should return 1.0

     B ∘ inverse(B)

will return as expected

     1.0

Note that as mentioned when dealing with the _geometric product_, this is a element of _type_ ```kMultvec```

#### Projection

This function depends on the _inversion_, so it can also only be used with elements of _type_ ```kblades``` as _input_. For example,

     projection(A, B)
     projection(B, A)
     projection(A, C)

will return

     e123
     6.0e13
     e123

these _output_ elements are of _type_ ```kMultvec```

#### Rejection

Like the above functions, the rejection receives an element of _type_ ```kblade``` and returns an element of _type_ ```kMultvec``` which is the _rejection_ of the first _input_ parameter from the second.

     rejection(B, kblade([kMultvec(kb(e2))]))

will return

     6.0e13

which is the difference between the ```kblade``` ```B``` and the projection ```rejection(B, kblade([kMultvec(kb(e2))]))```.
