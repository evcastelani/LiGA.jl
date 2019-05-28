# LIGA - Library for Geometric Algebra

[![Build Status](https://travis-ci.org/evcastelani/Liga.jl.svg?branch=master)](https://travis-ci.org/evcastelani/Liga.jl)

[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://evcastelani.github.io/Liga.jl/)

[![Coverage Status](https://coveralls.io/repos/evcastelani/Liga.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/evcastelani/Liga.jl?branch=master)

[![codecov.io](http://codecov.io/github/evcastelani/Liga.jl/coverage.svg?branch=master)](http://codecov.io/github/evcastelani/Liga.jl?branch=master)

It is dev version of Liga. 

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

