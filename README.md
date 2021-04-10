# Software for LL1-based hyperspectral super-resolution and blind unmixing with variable images
This repository contains the software for joint HSR and unmixing using LL1 decomposition.

Copyright (c) 2021 Clémence Prévost, Ricardo Borsoi, Konstantin Usevich, José Bermudez, David Brie, Cédric Richard <br>
Contact: ```clemence.prevost@univ-lorraine.fr```

This MATLAB software reproduces the results from the following papers:

```
Incoming
```

## Acknowledgements

The baseline algorithms used in the manuscript are courtesy of their respective authors.


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
  
  The synthetic datasets are available in the ```\data``` folders. 
  The Ivanpah Playa and Lake Tahoe datasets are available online https://landsat.gsfc.nasa.gov.
  
  ## Reproduce figures from the paper
  
  To do so, you need to run the ```main.m``` file. Here, a menu is available and allows you to choose which figure or table you want to generate. Each number in the table below corresponds to a set of figures.

| Number | Content                                        |
|--------|------------------------------------------------|
| 1      | Image fusion - Lake Tahoe dataset              |
| 2      | Image fusion - Ivanpah Playa dataset           |
| 3      | Image fusion - Indian Pines dataset            |
| 4      | Unmixing - field-like synthetic dataset        |
| 5      | Unmixing - synthetic dataset, non-id. NMF      |
| 6      | Unmixing - Lake Tahoe dataset                  |
| 7      | Unmixing - Ivanpah Playa dataset               |



