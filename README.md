# TLmodels

TLModels.jl (TLM) is a package extension for PowerSystems.jl (PSY). It extends PSY's modeling capabilities in two primary ways:
1. Allows the use of measured real life transmission line impedance data.
2. Allows the use of four transmission line models with these data, including two new high-fidelity transmission line models. 

This package was presented in the following paper: ***G. E. Colon-Reyes, R. Kravis, S. Sharma, and D. Callaway,
“Transmission line dynamics on inverter-dominated grids: analysis and
simulations,” 2023. [Online]. Available: https://arxiv.org/abs/2310.08553***. Please cite this paper if you will be using this package for your work.

TLM offers four TL models:
1. ***Statpi*** - $\pi$ line topology with algebraic equations
2. ***Dynpi*** - $\pi$ line topology with dynamic equations
3. ***MSSB*** - Series-connected multi-segment single-branch $\pi$ line topology
4. ***MSMB*** - Series-connected multi-segment multi-branch $\pi$ line topology


A tutorial of how to use this package can be found in A tutorial of how to use this package can be found in [test\TLmodels_tutorial.ipynb.](https://github.com/gecolonr/TLModels.jl/blob/main/test/TLmodels_tutorial.ipynb)


[![Build Status](https://github.com/gecolonr/TLmodels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gecolonr/TLmodels.jl/actions/workflows/CI.yml?query=branch%3Amain)
