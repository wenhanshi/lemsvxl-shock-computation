# dbsk2d-ishock-computation
This is the package for intrinsic shock computation, extracted from Maruthi's version in lemsvxl.
We remove the dependencies on other parts of lemsvxl, and make VXL be its only dependency. We redefine
the output of ishock conputation package as .isf/.osf file, and design a proper format to 
save and load.

__Input:__ image file (`.png`/`.jpg` ...) and its contour (`.cem`/`.cemv`)  
__Output:__ intrinsic shock graph (`.isf`) / output shock (`.osf`)

![ishock](https://github.com/wenhanshi/markdown-img-link/blob/master/ishock.png)

## Structure of Shock computation

__[ishock computation](https://github.com/wenhanshi/dbsk2d-ishock-computation) (we are here!):__  
`.png`/`.jpg` + `.cem`/`.cemv` --> `.isf`/`.osf`  
It's a stand alone package for generating __intrinsic shock graph__ from __image contours__.

__[xshock computation](https://github.com/wenhanshi/dbsk2d-xshock-computation):__  
`.png`/`.jpg` + `.cem`/`.cemv` --> `.esf`  
It's a stand alone package for generating __extrinsic shock graph__ from __image contours__.

__[osf_to_esf](https://github.com/wenhanshi/osf-to-esf):__  
`.osf` --> `.esf`  
It's a stand alone package for fastly generating __extrinsic shock graph__ from __ishock file__.


## User Guide

### 0. Preparation

- Original image: It can be `.jpg` or `.png`.
- Contour file: If you do not have contour file of the image, use [MSEL_contour_extraction](https://github.com/yg13/MSEL_contour_extraction_cxx) to get `.cem` (version 2) from image.  
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
$ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 99
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

### 2. Configure dbsk2d-ishock-computation

Build the package with CMake
e.g. under dbsk2d-ishock-computation root dir:
```commandline
$ mkdir test_build
$ cd test_build
$ cmake .. -G "Unix Makefiles"
$ make 
```

### 3. Run

To use it,
```commandline
$ ./dbsk2d-ishock-computation [-x input.xml] [-print-def-xml] [-?]
```
__Print default parameters XML file: -print-def-xml__  
e.g.
```commandline
$ ./dbsk2d-ishock-computation -print-def-xml
```
You will get `input-defaults.xml` with the content - __input parameters__:  
    `input_contour_extention=".cem"`: specific extension of input contour file (.cem/.cemv)  
    `input_image_extention=".png"`: specific extension of input image file  
    `input_object_dir="/YOUR_PATH"`: input directory  
    `input_object_name="test"`: file name, e.g. `test.cemv` and `test.png`  
    `Save_ISF-isfoutput="ishock.isf"`: the output file name for ishock saver, i.e. .isf  
    `Load_ISF-isfinput="ishock.isf"`: the input file name for ishock loader (only for test)   
    ...  
__Use specific configuration to run__  
e.g.  
```commandline
$ ./dbsk2d-ishock-computation -x my_input.xml
```
or you can use the default input
```commandline
$ ./dbsk2d-ishock-computation -x input_defaults.xml
```

## IShock Graph Computation Process

![process](https://github.com/wenhanshi/markdown-img-link/blob/master/stream.png)

## .isf File

We use `.isf` file as the output of the package. It is text file, including some essential
ishock info for the following procedures, such as pruning, sampling and composite fragment computation.

Importantly, we do not know about the specific following procedures using ishock results.
As a result, we decide to store only key information on ishock graph. Fortunately, the .isf
is easy to save and load with the saver and loader in the package. Developer may modify 
the saver and loader to add more information on ishock graph you need.

Currently, we only store boundaries, ishock nodes and ishock edges in `.isf`:

### Boundaries

Format:  
```
# ==============BOUNDARY=================
# Boundary Elements:
# [ID] [TYPE] [x1, y1] [x2, y2] ... [xn, yn]
[BOUNDARY_BEGIN]
497 1 243.915 189.065 243.399 188.706 
498 1 243.399 188.706 243.915 189.065
# ...
```

The boundary information includes:  

- ID: unique id for every boundary element
- TYPE: due to __line-fitting__, only two possible types here: 0(bpoint) and 1(bline)
- <X, Y>: one or more points of the boundary

### IShock Node

Format:
```
# ==============ISHOCK NODE==============
# IShock Nodes:
# [ID] [x, y] [CHILD_ISHOCK_NODE_ID] [BND_1_ID] [BND_2_ID] ...
[ISHOCK_NODE_BEGIN]
6 243.399 188.706 116614 501 498 
10 241.606 188.299 118705 503 506 
# ...
```

The ishock node information includes:

- ID: unique id for every ishock element (node and edge)
- <X, Y>: point extrinsic coordinate for this ishock node
- CHILD_ISHOCK_NODE_ID: index for child ishock node, i.e. parent_node ---> child_node
- BND_1_ID: group of boundaries generating this ishock node

_Note: CHILD_ISHOCK_NODE_ID can be -1, which means no child ishock node_

### IShock Edge

Format:
```
# ==============ISHOCK EDGE==============
# IShock Edges:
# [ID] [TYPE] [FROM_NODE_ID] [TO_NODE_ID] [LEFT_BND_ID] [RIGHT_BND_ID]
[ISHOCK_EDGE_BEGIN]
1 8 -1 122330 497 3
2 8 -1 122540 3 500
# ...
```

The ishock edge information includes:

- ID: unique id for every ishock element (node and edge)
- TYPE: due to __line-fitting__, only 1(point-point), 2(point-line), 4(line-line), 8(contact)
- FROM_NODE_ID & TO_NODE_ID: from_node ---> to_node (or, parent and child)
- LEFT_BND_ID & RIGHT_BND_ID: index of two boundaries which generate this ishock edge

_Note: FROM_NODE_ID/TO_NODE_ID can be -1, which means no from/to ishock node_

## .osf File

The same as `.isf` file, `.osf` file can be the output of ishock computation. However, .isf 
only contains necessary information for creating ishock graph. Developer should design his
own loader for `.isf` to make it work in the following functional packages.

`.osf` is the __super set__ of `.isf`. It contains necessary as well as supplementary information
on ishock graph output. Without changing any details of algorithms in the following packages, 
e.g. xshock computation or composite fragment computation, they will get completely the same results.

Based on `.isf`, we add more additional members to `.osf`, especially math parameters.

![osf_define](https://github.com/wenhanshi/markdown-img-link/blob/master/osf_define.png)

### Boundary Point

Format:
```
# ==============ISHOCK FILE==============
# ==============BOUNDARY POINT=================
# Boundary points:
# [ID] [TYPE] [x] [y]
[BOUNDARY_POINT_BEGIN]
```

### Boundary Line

Format:
```
# ==============BOUNDARY LINE=================
# Boundary lines:
# [ID] [TYPE] [START_POINT] [END_POINT] [U] [N] [L]
[BOUNDARY_LINE_BEGIN]
```

### Ishock Node

Format:
```
# ==============ISHOCK NODE==============
# IShock Nodes:
# [ID] [x] [y] [START_TIME] [END_TIME] [CHILD_SHOCK_LINK_ID_1] [CHILD_SHOCK_LINK_ID_2] [NUM_BND] [BND_1_ID] [BND_2_ID] ...
[ISHOCK_NODE_BEGIN]
```

### Ishock Edge

Format:
```
# ==============ISHOCK EDGE==============
# IShock Edges:
# [ID] [TYPE] [START_TIME] [END_TIME] [LS_ETA] [RS_ETA] [LS_TAU] [LE_TAU] [RS_TAU] [RE_TAU] [H] [FROM_NODE_ID] [TO_NODE_ID] [LEFT_BND_ID] [RIGHT_BND_ID] [OTHER_INFO]
# [OTHER_INFO]: 
# for point-point: [N] [U]
# for point-line: [NU] [U] [N] [L_DELTA] [R_DELTA] [L]
# for line-line: [UR] [PHI] [UL] [SIGMA] [THETA_L] [THETA_R] [LL] [LR]
# for contact: [N]
[ISHOCK_EDGE_BEGIN]
```

### Other Infomation

Format:
```
# ==============OTHER INFO==============
# Other info (pshocks in ishock node):
# [ISHOCK_NODE_ID] [NUM_PSHOCK] [PSHOCK_ID_1] ...
```

## Features

- Retain the previous structure of directory, such as `/dbgl/algo/dbgl_biarc.cxx`.
- Remove CMakeList in sub dirs, use global one to control.
- Get the 'closure' of functional package of ishock computation (smaller black box).
- Turn LEMSVXL dependencies like `dbgl`, `bpro1`, `vidpro1` ... much lighter.
- Retain the dependencies on VXL
- Use XML to control the parameters
- Use `.isf` file as the output
- Use `.osf` file as the output for making sure the whole shock package can work without changing any details

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
