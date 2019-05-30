# LIGA - Library for Geometric Algebra

[![Build Status](https://travis-ci.org/evcastelani/Liga.jl.svg?branch=master)](https://travis-ci.org/evcastelani/Liga.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://evcastelani.github.io/Liga.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://evcastelani.github.io/Liga.jl/dev)
[![codecov.io](http://codecov.io/github/evcastelani/Liga.jl/coverage.svg?branch=master)](http://codecov.io/github/evcastelani/Liga.jl?branch=master)

It is new version of Liga. 

### Basic usage

Once the files have been downloaded, you need to import `NewLiga`. In order to do that, type:

``` 
push!(LOAD_PATH,"path to file/NewLiga/")
using NewLiga
```
After that, Liga was imported. To use, just need to type:

```
layout(3,1,"Conformal")
```

Consequently, the Conformal enviroment G(3,1) is created. If you need to setup other enviroment, you need to restart Julia and import Liga again.

