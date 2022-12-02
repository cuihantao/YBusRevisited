# YBusRevisited

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://YBusRevisited.github.io/YBusRevisited.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://YBusRevisited.github.io/YBusRevisited.jl/dev)
[![Build Status](https://github.com/YBusRevisited/YBusRevisited.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/YBusRevisited/YBusRevisited.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/YBusRevisited/YBusRevisited.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/YBusRevisited/YBusRevisited.jl)


This repository contains the code for a preprint.

## Installation

Clone the repository, start `julia` in the root directory of the repository,
and run the following commands:

```julia
]activate .
]instantiate

# Replace the following path to the location of MATPOWER's `data` folder on your machine
julia> MATPOWER_DATA_PATH="path/to/matpower/data/folder"
julia> include("src/timing_launcher.jl")
```

The program will benchmark the Ybus method and the element-wise method for nine
test systems. It may take several minutes to complete.
