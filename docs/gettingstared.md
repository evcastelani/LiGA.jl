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
be in mind is *setup* the layout. The ```layout.jl``` is a file that write
an ```objects.jl```. The ```objects.jl``` generates a basis for some space that we want work. For example, if we type

	 layout(3)

we are generating  the G3 space as default, that is, we are generating a basis of vectors and multivectors for G3. But, we can work (with appropriated setup of elements),with G3,1 and G4,1. 

Let us started with more basic example. After type ```layout(3)```  was created bases elements. Type, for example

     e1

This is a logical array with value 

     [false,true,true]

All bases elements were created using this idea. Try type others basis elements
    
     e123

     e12       

        
An  observation in this case is the order of these elements are important. If you type

     e21

an error mensage is showed. This does not mean that this element can not be built. We will see this later. Reciprocally, type

     kb([true,true,false])

This command show the multivector ```e12```. So ```kb``` is abbreviation of *k-basis element* and if we think in pratical purposes, it is just a logical vector. 
Another important property is ```kb``` supports a scalar argument with the basis element. For example,

    -2.3e12

is equal to type

    kb([true,true,false],-2.3)

All these commands provide the basics for starting more elaborate constructions in the space that has been defined. 

Consequently, we recommend that the reader calmly review documentation on the types defined.        
