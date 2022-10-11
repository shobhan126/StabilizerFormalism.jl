# StabilizerFormalism

This small repository is intended to setup a DSL for Stabilizer Formalism which popular among the quantum error correction community.

WIP at the moment as I am writing this while taking a Quantum Fault Tolerance course and use it for my fun projects.

# Cloning the package

In your terminal, navigate to where you want to store the package type:
```bash
git clone git@github.com:shobhan126/StabilizerFormalism.jl.git
```
(This requires you have `ssh` key setup for your device)
Alternatively this might work:
```bash
git clone https://github.com/shobhan126/StabilizerFormalism.jl.git
```


# Adding the package to Julia
Currently this package is not registered in the julia registry; you may add the package in your environment using the falling method

Open your `julia` REPL and type the following

```sh
julia>] # ] in repl will take you to package maode
pkg> add git@github.com:shobhan126/StabilizerFormalism.jl.git
```

Alternatively you may open a jupyter notebook and evaluate the following cell:


```julia

using Pkg
Pkg.as
Pkg.add("url=https://github.com/shobhan126/StabilizerFormalism.jl.git", rev='main")
```

# Quick Usage
```julia
# creating a Stabilizer from macro
julia> s1 = @pauli XIYZ # stabilizer 1
XIYZ
# creating a stabilizer using a symbol
julia> s2 = Pauli(:IXYZ)
IXYZ

# creating a stabilizer using string macro
julia> s3 = p"XIXI"
XIXI

# creating a stabilizer using string
julia> s4 = Pauli("YZYZ")

# product of two stabilizers
julia> s5 = s1 * s3
XXII

# product of stabilizers that dont commute
julia> p"XIXI" * p"YIII"
imZIXI

# negating a stabilizer
julia> -s1
-XIYZ

# multiplying by imaginary
julia> s6 = -im * s1
-imXIYZ

# absolute valus
julia> abs(s6)
XIYZ

# checking equality
julia> s7 = Pauli("xiyz")

julia> s1 == s7
true

julia> s1 == s2
false

# getting the symplecting representation
julia> symplectic(s1) |> print
Bool[1, 0, 1, 0, 0, 0, 1, 1]


# check matrix for a set of Stabilizers
julia> checkmatrix([s1, s2, s3])
3×8 BitMatrix:
 1  0  1  0  0  0  1  1
 0  1  1  0  0  0  1  1
 1  1  0  0  0  0  0  0

 # check matrix will discard the the imaginary unit / negative sign
```
