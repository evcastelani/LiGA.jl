## Installation


(In the future) Open Julia REPL and type

     Pkg.add("Liga")

Or (by now) install it yourself as:

     git clone https://github.com/evcastelani/Liga.jl

## Layouts

As cited, just the second option is avaiable. So, you can type 

```@repl 1
using Liga 
```

In order to start to use Liga, the first point that you must
be in mind is setup the layout. The ```layout``` function defines the space and
all logical elements for build \$ \mathbb{G}_k \$. For example, if we want create the \$\mathbb{G}_3\$ space, we need just type

```@repl 1
layout(3)
```

In this case, we are generating a logical structure for create basis of vectors and multivectors for \$\mathbb{G}_3\$. Naturally, we can work (with appropriated setup of elements), with \$\mathbb{G}_4\$ and \$\mathbb{G}_{4,1}\$. 

Let us started with more basic example. After type ```layout(3)```  was created base logical elements. Type, for example

```@repl 1
e1 
```

as we can see ```e1``` is a logical array with value 
```[false,true,true]```


All base elements were created using this idea. Try type others logical basis elements

```@repl 1
e123

e12       
```
        
An  observation in this case is the order of the index of these elements are important. If you type

```@repl1
e21
```

an error mensage is showed. This does not mean that this element can not be built. We will see this later. There is a subtle difference between logical elements and *k-base elements*. The *k-basis element* can provide more information, like a logical base element and a scalar associated. To create a *k-basis elements* we need to use the command *kb*. For example, type

```@repl 1
kb(e1)
```

and after
```@repl 1
typeof(kb(e1))
```

Naturally, ```kb(e1)``` is a *k-basis element*. In this case, in order to acess the information related to ```kb(e1)``` just ype

```@repl 1
kb(e1).e
kb(e1).scl
```

So, we can see that an *kb* element needs a logical basis element and a scalar element. We can generate others *kb* elements by type

```@repl 1
kb(e1,-2.0)
```

There are other syntaxes functions that are equivalent to command ```kb```. For example,

```@repl 1
kb([true,false,false])
```
or 
```@repl 1
kbasis(e1,1.0)
```

It is obvious that ```kb``` is abbreviation of *k-basis element* and if we think in pratical purposes, it is just a logical vector. 

In current version, ```kb``` supports a more direct definition. For example, if we want to define the element ```kb(e12,2.3)``` we can just type in Julia REPL, the command

```@repl 1
-2.3*e12
```

wich is equivalent to

```@repl 1
kb([true,true,false],-2.3)
```

Concerning the diference between *e12* and *e21*, if by some reason, we need to define the *e21* element, we can just define the element ```kb(e12,-1.0)```.

Lets to consider a more complete example. The most import function in GA context is the *geometric product*. So, we can show an example related to this function.
For this, just type these commands in Julia REPL 

```@repl 1

u=5.0*e123

v=1.0*id+3.0*e12

v∘u #(v \circ+TAB u)
```
and a gemetrical product between the multivector u and the multivector v is showed. Note that, we didn't define a *multivector type*, but it is very intuitive that they are composed of *k-basis elements*.

All these commands provide the basics for starting more elaborate constructions in the space that has been defined.  
Consequently, we recommend that the reader calmly review documentation on the types defined and functions supported. Naturally, we introduce the gemetrical product in \$G_3\$ but Liga provides more functions and other kind of elements.

  
