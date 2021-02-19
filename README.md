# Software for LL1-based hyperspectral super-resolution and blind unmixing with variable images
This repository contains the software for joint HSR and unmixing using LL1 decomposition.

Copyright (c) 2021 Clémence Prévost, Ricardo Borsoi, Konstantin Usevich, José Bermudez, David Brie, Cédric Richard <br>
Contact: ```clemence.prevost@univ-lorraine.fr```

This MATLAB software reproduces the results from the following papers:

```
Incoming
```

## Acknowledgements

CP-based STEREO and Blind-STEREO are courtesy of C. Kanatsoulis. If you want to use this software for publication, please cite the following papers:

```
@article{kanatsoulis2018hyperspectral,
  title={Hyperspectral super-resolution: A coupled tensor factorization approach},
  author={Kanatsoulis, Charilaos I and Fu, Xiao and Sidiropoulos, Nicholas D and Ma, Wing-Kin},
  journal={IEEE Transactions on Signal Processing},
  volume={66},
  number={24},
  pages={6503--6517},
  year={2018},
  publisher={IEEE}
}
``` 

## Content

 - demo.m : demo file with minimal requirements.
 
 - /data : contains synthetic datasets and initial factors.
 
 - /demos : contains demo files that produce tables and figures (including ```demo.m```).

 - /figures : where the figures are saved.

 - main.m : codes that allow to run desired demos.
 
 - /metrics : contains the comparison metrics used to assess performance.
 
 - /src : contains helpful files to run the demos.

## Minimal requirements

In order to run the demo file and reproduce the figures, you will need to:
- Download and install Tensorlab 3.0: https://www.tensorlab.net

## Demo file
 
 A demo with minimal requirements is available. To proceed, please run the ```demo.m``` file.
 
  ### Load data

  ### Perform fusion
  
  ### Perform blind unmixing
  

  
  ## Reproduce figures from the paper
  
  To do so, you need to run the ```main.m``` file. Here, a menu is available and allows you to choose which figure or table you want to generate. Each number in the table below corresponds to a set of figures.

| Number | Content                                        |
|--------|------------------------------------------------|
| 1      | produces Figure 1 from CAMSAP paper            |
| 2      | produces Fig. 1 and 4 from preprint            |
| 3      | produces Fig. 2 and 5 from preprint            |
| 4      | produces Fig. 3 from preprint                  |
| 5      | produces Fig. 6 from preprint                  |
| 6      | produces Fig. 8, 9 and 10 from preprint        |
| 7      | produces Fig. 11 and 12 from preprint          |



