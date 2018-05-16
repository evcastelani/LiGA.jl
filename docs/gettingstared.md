---
layout: page
title: Getting Started
navigation: 2
---

## Installation


(In the future) Open Julia REPL and type

     Pkg.add("Liga")

Or (by now) install it yourself as:

     git clone https://github.com/evcastelani/Liga.jl

## Layouts

As cited, just the second option is avaiable. So, you can type 

	 using "Liga"

This file contains all functions of the package. The first point that you must
be in mind is *setup* the layout. The ```layout``` function define the space and
all logical elements for build G(k). For example, if we want create the G(3) space, just type

	 layout(3)

In this case, we are generating a logical structure for create basis of vectors and multivectors for G(3). 
Naturally, we can work (with appropriated setup of elements),with G(4) and G(4,1). 

Let us started with more basic example. After type ```layout(3)```  was created base logical elements. Type, for example

     e1

This is a logical array with value 

     [false,true,true]

All base elements were created using this idea. Try type others logical basis elements
    
     e123

     e12       

        
An  observation in this case is the order of the index of these elements are important. If you type

     e21

an error mensage is showed. This does not mean that this element can not be built. We will see this later. There is a subtle difference between logical elements and *k-base elements*. The *k-basis element* can provide more information, like a logical base element and a scalar associated. To create a *k-basis elements* we need to use the command *kb*. For example, type

     kb(e1)

and after

     typeof(e1)

Naturally, ```kb(e1)``` is a *k-basis element*. In this case, in order to acess the information related to ```kb(e1)``` just ype

     kb(e1).e
     kb(e1).scl

So, we can see that an *kb* element needs a logical basis element and a scalar element. We can generate others *kb* elements by type

     kb(e1,-2.0)
     
Other syntax functions that are equivalent to command ```kb```. For example,

     kb([true,false,false])

or 

     kbasis(e1,1.0)


```kb``` is abbreviation of *k-basis element* and if we think in pratical purposes, it is just a logical vector. 

In current version, ```kb``` supports a more direct definition. For example, if we want to define the element ```kb(e12,2.3)``` we can just type in Julia 
REPL, the command

    -2.3*e12

wich is equivalent to

    kb([true,true,false],-2.3)

Return to diference between *e12* and *e21*, if by some reason, we need to define the *e21* element, we can just define the element ```kb(e12,-1.0)```.

All these commands provide the basics for starting more elaborate constructions in the space that has been defined. 

Consequently, we recommend that the reader calmly review documentation on the types defined.        
