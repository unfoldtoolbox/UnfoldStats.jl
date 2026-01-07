# UnfoldStats

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://unfoldtoolbox.github.io/UnfoldStats.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://unfoldtoolbox.github.io/UnfoldStats.jl/dev/)
[![Build Status](https://github.com/unfoldtoolbox/UnfoldStats.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/UnfoldStats.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/unfoldtoolbox/UnfoldStats.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/unfoldtoolbox/UnfoldStats.jl)

|rERP|EEG visualisation|EEG Simulations|BIDS pipeline|Decode EEG data|Statistical testing|
|---|---|---|---|---|---|
| <a href="https://github.com/unfoldtoolbox/Unfold.jl/tree/main"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623787-757575d0-aeb9-4d94-a5f8-832f13dcd2dd.png"></a> | <a href="https://github.com/unfoldtoolbox/UnfoldMakie.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623793-37af35a0-c99c-4374-827b-40fc37de7c2b.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldSim.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623795-328a4ccd-8860-4b13-9fb6-64d3df9e2091.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldBIDS.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622460-2956ca20-9c48-4066-9e50-c5d25c50f0d1.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldDecode.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622487-802002c0-a1f2-4236-9123-562684d39dcf.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldStats.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623799-4c8f2b5a-ea84-4ee3-82f9-01ef05b4f4c6.png"></a>|

A toolbox for statistical testing of EEG/ERP data and Unfold.jl models.
Build on the [Unfold](https://github.com/unfoldtoolbox/unfold.jl/).

> [!CAUTION]
> This package is work in progress. We have not (yet) extensively checked the validity and correctness of the statistical tests. While we think these are valid approaches, we have ongoing studies investigating whether this is actually the case. If you notice some unexpected behavior please open an issue!

## Install

<SUPPORT LIST TO BE ADDED>

### Installing Julia

<details>
<summary>Click to expand</summary>

The recommended way to install julia is [juliaup](https://github.com/JuliaLang/juliaup).
It allows you to, e.g., easily update Julia at a later point, but also test out alpha/beta versions etc.

TL:DR; If you dont want to read the explicit instructions, just copy the following command

#### Windows

AppStore -> JuliaUp,  or `winget install julia -s msstore` in CMD

#### Mac & Linux

`curl -fsSL https://install.julialang.org | sh` in any shell
</details>

### Installing Unfold

```julia
using Pkg
Pkg.add("UnfoldStats")
```

### Contributions

Contributions are very welcome. These could be typos, bugreports, feature-requests, speed-optimization, new solvers, better code, better documentation.

### How-to Contribute

You are very welcome to raise issues and start pull requests!

### Adding Documentation

1. We recommend to write a Literate.jl document and place it in `docs/literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `docs/src/generated/FOLDER/FILENAME.md`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/make.jl) with a reference to `docs/src/generated/FOLDER/FILENAME.md`.


## Citation

If you make use of these statistics, please cite:

<CITATION>

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/jschepers"><img src="https://avatars.githubusercontent.com/u/22366977?v=4?s=100" width="100px;" alt="Judith Schepers"/><br /><sub><b>Judith Schepers</b></sub></a><br /><a href="#infra-jschepers" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#test-jschepers" title="Tests">⚠️</a> <a href="#code-jschepers" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://www.benediktehinger.de"><img src="https://avatars.githubusercontent.com/u/10183650?v=4?s=100" width="100px;" alt="Benedikt Ehinger"/><br /><sub><b>Benedikt Ehinger</b></sub></a><br /><a href="#code-behinger" title="Code">💻</a> <a href="#bug-behinger" title="Bug reports">🐛</a> <a href="#doc-behinger" title="Documentation">📖</a> <a href="#infra-behinger" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/rosagross"><img src="https://avatars.githubusercontent.com/u/43537429?v=4?s=100" width="100px;" alt="rosagross"/><br /><sub><b>rosagross</b></sub></a><br /><a href="#bug-rosagross" title="Bug reports">🐛</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->



This project follows the [all-contributors](https://allcontributors.org/docs/en/specification) specification. 

Contributions of any kind welcome!
You can find the emoji key for the contributors [here](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/contrib-emoji.md)


## Acknowledgements



## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
