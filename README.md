# lemsvxl-shock-computation

The shock computation packages are built for dealing with the streamline of generating Intrinsic Shock Graph and Extrinsic Shock Graph
from the source image and its contours.

See more details about examples, tests and performance in report.

## Streamline

![streamline](https://github.com/wenhanshi/markdown-img-link/blob/master/two%20pipelines.png)

## Structure of Shock Computation Packages

__[ishock computation](https://github.com/wenhanshi/dbsk2d-ishock-computation):__  
`.png`/`.jpg` + `.cem`/`.cemv` --> `.isf`/`.osf`  
It's a stand alone package for generating __intrinsic shock graph__ from __image contours__.

__[xshock computation](https://github.com/wenhanshi/dbsk2d-xshock-computation):__  
`.png`/`.jpg` + `.cem`/`.cemv` --> `.esf`  
It's a stand alone package for generating __extrinsic shock graph__ from __image contours__.

__[osf_to_esf](https://github.com/wenhanshi/osf-to-esf):__  
`.osf` --> `.esf`  
It's a stand alone package for fastly generating __extrinsic shock graph__ from __ishock file__.


