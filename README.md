# PhysicalTensors.jl

This module provides a type for physical tensors, overloading both math operators and standard functions.

To use this module you will need to add the following Julia packages to yours:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/MutableTypes.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalSystemsOfUnits.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalScalars.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalVectors.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalTensors.jl")
```

## Re-exported from PhysicalFields.jl

Package [Reexport.jl](https://github.com/simonster/Reexport.jl), or [here](https://juliapackages.com/p/reexport), is used to export symbols from an imported package. Here `Reexport` is used to re-export part of the interface of package [PhsysicalFields.jl](https://github.com/AlanFreed/PhysicalFields.jl), effectively linking these exports so they appear to originate from within this module. Specifically, using the alias

```
const PhysicalUnits = PhysicalSystemsOfUnits.PhysicalSystemOfUnits
```

the following type definition for a physical tensor field is exported as

```
struct PhysicalTensor <: PhysicalField
    r::UInt16           # rows in the matrix representation of the tensor
    c::UInt16           # columns in the matrix representation of the tensor
    m::StaticMatrix     # values of the tensor in its specified system of units
    u::PhysicalUnits    # the tensor's physical units
end
```

which is a sub-type to the abstract type

```
abstract type PhysicalField end
```

Also, a type definition for an array of such tensor fields is exported as

```
struct ArrayOfPhysicalTensors
    e::UInt16           # number of entries or elements held in the array
    r::UInt16           # rows in each physical tensor held in the array
    c::UInt16           # columns in each physical tensor held in the array
    a::Array            # array of matrices holding values of a physical tensor
    u::PhysicalUnits    # units of this physical tensor
end
```

where all entries in the array have the same physical units.

Constructors for these types are also re-exported here, they being

```
function newPhysicalTensor(rows::Integer, cols::Integer, units::PhysicalUnits)::PhysicalTensor
```

which supplies a new tensor of dimension `rows` by `cols` whose values are set to `0` and whose physical units are those supplied by the argument `units`, and

```
function newArrayOfPhysicalTensors(len::Integer, t₁::PhysicalTensor)::ArrayOfPhysicalTensors
```

where `t₁` is the first entry in a new array of tensors whose length is `len`.

To retrieve and assign array values, functions `Base.:(getindex)` and `Base.:(setindex!)` have been overloaded so that the bracket notation `[]` can be used to: *i)* retrieve and assign scalar fields belonging to an instance of `PhysicalTensor`, and *ii)* retrieve and assign tensor fields belonging to an instance of `ArrayOfPhysicalTensors`.

## Type conversions

Conversion to a string is provided for instances of `PhysicalTensor` by the function

```
function toString(t::PhysicalTensor; format::Char='E')::String
```

where the keyword `format` is a character that, whenever its value is 'E' or 'e', represents the tensor components in a scientific notation; otherwise, they will be represented in a fixed-point notation.

Conversion to a matrix of the values held by a tensor is provided by
```
function toMatrix(t::PhysicalTensor)::StaticMatrix
```

## Operators

The following operators have been overloaded so that they can handle objects of type PhysicalTensor, whenever such operations make sense, e.g., one cannot add two tensors with different units or different dimensions. The overloaded logical operators include: `==`, `≠` and `≈`. The overloaded unary operators include: `+` and `-`. And the overloaded binary operators include: `+`, `-`, `*`, `/` and `\`, where the latter solves a linear system of equations, e.g., **Ax** = **b** solved for **x**, which is written in code as `x = A\b`.

## Methods for both PhysicalTensor and ArrayOfPhysicalTensors

The following methods can accept arguments that are objects of either type, viz., PhysicalTensor or ArrayOfPhysicalTensors. They are self explanatory: `copy`, `deepcopy`, `isDimensionless`, `isCGS`, `isSI`, `toCGS` and `toSI`.

## Math functions for PhysicalTensor

To construct a tensor product between two vectors, typically written as `t = y⊗z` or in component notation as `tᵢⱼ = yᵢzⱼ`:
```
function tensorProduct(y::PhysicalVector, z::PhysicalVector)::PhysicalTensor
```

To get the transpose of a tensor:
```
function Base.:(transpose)(t::PhysicalTensor)::PhysicalTensor
```

To get the trace of a tensor:
```
function LinearAlgebra.:(tr)(t::PhysicalTensor)::PhysicalScalar
```

To get the determinant of a tensor:
```
function LinearAlgebra.:(det)(t::PhysicalTensor)::PhysicalScalar
```

To get the inverse of a tensor, when it exists:
```
function Base.:(inv)(t::PhysicalTensor)::PhysicalTensor
```

To get the **QR** or Gram-Schmidt decomposition of a tensor, when it exists:
```
function LinearAlgebra.:(qr)(t::PhysicalTensor)::Tuple
```
e.g., written in code as
```
(Q, R) = qr(t)
```
where **Q** is an orthogonal Gram rotation matrix, and **R** is an upper-triangular (right) matrix.

To get the **LQ** decomposition of a tensor, when it exists:
```
function LinearAlgebra.:(lq)(y::PhysicalMatrix)::Tuple
```
e.g., written in code as
```
(L, Q) = lq(t)
```
where **L** is a lower-triangular (left) matrix, and **Q** is an orthogonal rotation matrix.