# Jessamine.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wgm-applied-math.github.io/Jessamine.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wgm-applied-math.github.io/Jessamine.jl/dev/)
[![Build Status](https://github.com/wgm-applied-math/Jessamine.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wgmitchener/Jessamine.jl/actions/workflows/CI.yml?query=branch%3Amain)


## About

Jessamine is a collection of [Julia](https://www.julialang.org) packages for machine learning, specifically, evolutionary symbolic regression and classification using static-single-assignment-form expressions.
It is a research project under development and not (yet) easy to use.
Expect ongoing improvements and breaking changes.

An article explaining the design and operation of Jessamine is available as [LaTeX source](https://github.com/wgm-applied-math/Jessamine-Development-article).
Once complete, it will be posted on arXiv and submitted for publication in an academic journal.

## Installation

As of 2026-08-18, I have not registered this package.
To use this package within a Julia project, use the [Pkg.jl](https://pkgdocs.julialang.org/v1/) command line,
```julia-repl
pkg> add https://github.com/wgm-applied-math/Jessamine.jl#main
```
You can use a tag or branch name in place of `main`.

## Further details

This package is the core functionality: genome, evolution, ridge regression.
For additional functionality, see the following:

- [JessamineSymbolics.jl](https://github.com/wgm-applied-math/JessamineSymbolics.jl):
  Core package for using [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) with Jessamine.

- [JessamineCLI.jl](https://github.com/wgm-applied-math/JessamineCLI.jl):
  Command-line interface for some essential functionality.
  This is probably the most useful package to actually try to use right now.

- [JessamineBenchmark.jl](https://github.com/wgm-applied-math/JessamineBenchmark.jl):
  Collection of symbolic-regression benchmarks on standard data sets.
  Includes scripts and configuration files useful for running symbolic regression on a Slurm cluster.

## Still under construction

- [JessamineSciKitLearn.jl](https://github.com/wgm-applied-math/JessamineSciKitLearn.jl)
  and
  [jessaminescikitlearn](https://github.com/wgm-applied-math/jessaminescikitlearn):
  Julia and Python packages for using Jessamine within the Python package [scikit-learn](https://scikit-learn.org/).

- [JessamineMLJ.jl](https://github.com/wgm-applied-math/JessamineMLJ.jl):
  Package for stacking algorithms from [MLJ.jl](https://github.com/JuliaAI/MLJ.jl) on top of Jessamine.
  
  
