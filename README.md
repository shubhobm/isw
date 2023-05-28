# Intrinsic Sliced Wasserstein Distances for Comparing Collections of Probability Distributions on Manifolds and Graphs

This repository contains code and data for the experiments in the paper

> [Intrinsic Sliced Wasserstein Distances for Comparing Collections of Probability Distributions on Manifolds and Graphs](https://arxiv.org/abs/2010.15285)
>
> Raif Rustamov and Subhabrata Majumdar
>
> International Conference on Machine Learning (ICML), 2023.

## Setup

You should have the [R statistical software](https://www.r-project.org) installed in your system. To install the dependencies, run the following in bash:
```
Rscript install_dependencies.R
```

## Contents
Below are the the contents of this repository.

| Subdirectory | Method |
|---|---|
| `scripts/sim_line.R` | ISW implementation on real line, synthetic experiments |
| `scripts/sim_circle.R` | ISW implementation on a circle, synthetic experiments |
| `scripts/sim_cyl.R` | ISW implementation on a cylinder, synthetic experiments |
| `scripts/chicago_crime.R` | ISW implementation on graph data, Chicago Crime |
| `data/chicago_theft2018.RData` | The Chicago Crime dataset |

## Contact
Please feel free to [email](mailto:zoom.subha@gmail.com) with questions of comments.
