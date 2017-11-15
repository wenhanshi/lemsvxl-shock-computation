# dbsk2d-xshock-computation
This is the package for extrinsic shock computation, extracted from Maruthi's version in lemsvxl.
We remove the dependencies on other parts of lemsvxl, and make VXL be its only dependency.

__Input:__ image file (`.png`/`.jpg` ...) and its contour (`.cem`/`.cemv`)  
__Output:__ extrinsic shock graph (`.esf`)

![xshock](https://github.com/wenhanshi/markdown-img-link/blob/master/xshock.png)

xshock computation package can be viewed as a refactor of Maruthi's whole shock computation package.

## Structure of Shock computation

__[ishock computation](https://github.com/wenhanshi/lemsvxl-shock-computation/tree/master/dbsk2d-ishock-computation)__:  
`.png`/`.jpg` + `.cem`/`.cemv` --> `.isf`/`.osf`  
It's a stand alone package for generating __intrinsic shock graph__ from __image contours__.

__[xshock computation](https://github.com/wenhanshi/lemsvxl-shock-computation/tree/master/dbsk2d-xshock-computation) (we are here!)__:  
`.png`/`.jpg` + `.cem`/`.cemv` --> `.esf`  
It's a stand alone package for generating __extrinsic shock graph__ from __image contours__.

__[osf_to_esf](https://github.com/wenhanshi/lemsvxl-shock-computation/tree/master/osf-to-esf)__:  
`.osf` --> `.esf`  
It's a stand alone package for fastly generating __extrinsic shock graph__ from __ishock file__.

## User Guide

### 0. Preparation

- Original image: It can be `.jpg` or `.png`.
- Contour file: Use [MSEL_contour_extraction](https://github.com/yg13/MSEL_contour_extraction_cxx) to get `.cem` (version 2) from image.  
Then use converter to get `.cem`/`.cemv` (version 1) from `.cem` (version 2).
- CMake ( >= v3.5 )
- gcc 4.x

_Note: you may need gcc-4.x to build this package .  
e.g.
Use `$ gcc --version` to check your gcc version.  
Download gcc-4.8_
```commandline
$ sudo apt-get install gcc-4.8
```
_Use gcc-4.8 in your environment._
```commandline
$ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 100
```
_Then the same as g++-4.x._

### 1. Configure VXL

Download the latest version (v1.17.0) of [VXL](https://github.com/vxl/vxl).

Build the VXL with CMake  
e.g. under VXL root dir:
```commandline
$ mkdir vxl-bin
$ cd vxl-bin
$ cmake .. -G "Unix Makefiles"
$ make
```

Set VXL_DIR to your build directory  
e.g.
```commandline
$ export VXL_DIR="~/work/vxl/src/vxl-bin"
```

### 2. Configure dbsk2d-xshock-computation

Build the package with CMake
e.g. under dbsk2d-xshock-computation root dir:
```commandline
$ mkdir test_build
$ test_build
$ cmake .. -G "Unix Makefiles"
$ make 
```

### 3. Run

To use it,
```commandline
$ ./dbsk2d-xshock-computation [-x input.xml] [-print-def-xml] [-?]
```
__Print default parameters XML file: -print-def-xml__  
e.g.
```commandline
$ ./dbsk2d-xshock-computation -print-def-xml
```
You will get `input-defaults.xml` with the content - __input parameters__:  
    `input_contour_extention=".cem"`: specific extension of input contour file (.cem/.cemv)  
    `input_image_extention=".png"`: specific extension of input image file  
    `input_object_dir="/YOUR_PATH"`: input directory  
    `input_object_name="test"`: file name, e.g. `test.cemv` and `test.png`    
    `output_extension=".esf"`: output file extension  
    `output_shock_folder="/YOUR_PATH"`: path to save `.esf`  
    ...  
__Use specific configuration to run__  
e.g.  
```commandline
$ ./dbsk2d-xshock-computation -x my_input.xml
```
or you can use the default input
```commandline
$ ./dbsk2d-xshock-computation -x input_defaults.xml
```

## Features

- Retain the previous structure of directory, such as `/dbgl/algo/dbgl_biarc.cxx`.
- Remove CMakeList in sub dirs, use global one to control.
- Get the 'closure' of functional package of ishock computation (smaller black box).
- Turn LEMSVXL dependencies like `dbgl`, `bpro1`, `vidpro1` ... much lighter.
- Retain the dependencies on VXL
- Use XML to control the parameters
- Use `.esf` as the output of xshock computation, easy to be visualized

## Guidance for Package Extraction

- Use __doxygen__ to find out main process.

- Figure out the start of INPUT and the end of OUTPUT, get initial files/functions.

- Calculate the closure based on these files, get __minimum functional set__.  
_Pay attention to `.txx` files. Do NOT copy content of `.txx` to `.cxx`,
since there may be some templates including these `.txx` files_

- Rebuild dir structure, make it independent package view.  
To be clearer, turn   
`#include <dbsk2d/dbsk2d_ishock_bcurve.h>`  
to  
`#include "dbsk2d_ishock_bcurve.h"`

- Rewrite __CMakeList__.

- Build your new package under the new CMakeList

## Some possible build errors

- __Can not find file or dir ...__  
Maybe your closure is not complete, or `#include ...` is not changed correctly

- __invalid use of incomplete type class xxx, forward declaration of xxx__  
This may caused by the usage of smart pointer.
That is, you only include `class_foo_sptr.h`, but not `class_foo.h`.

- __xxx function is used but not implemented.__  
This should appear during linking process. Find out the unimplemented functions and search them in 
original package. Some `.cxx` or `.txx` may be missed during the closure.
These missing files can also be in external libraries, check `TARGET_LINK_LIBRARIES` in CMakeList, make sure you
have included all libraries.  

    Note: Some missing files can be difficult to target due to bad naming habit.
    You should not use `foo.h` to control `foo.cxx`, `foo_bar1.cxx` and `foo_bar2.cxx`.


## Update: Compatible to Latest VXL Version (1.17.0)

Compared to the latest version of VXL, the previous one has different `vnl_math.h`.

In old version (<1.14.0), vnl_math functions are implemented as a namespace called `vnl_math`.
There are some math functions such as `max`, `min` and `hypot`. But in latest version, `vnl_math`
is redesigned as a class. Those functions are redesigned as `inline` ones. The function name
is changed from `foo()` to `vnl_math_foo()`.

That is, if there is `vnl_math::foo()` in old version, use `foo()` instead. If there is
variables like `vnl_math::pi`, just keep them the same.


_Written and updated by [Wenhan](mailto:shiwenhan@bupt.edu.cn)_
