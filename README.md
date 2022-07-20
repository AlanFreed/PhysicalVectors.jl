# PhysicalVectors.jl

This module provides a type for physical vectors, overloading both math operators and standard functions.

To use this module you will need to add the following Julia packages to yours:

```
using Pkg
Pkg.add(url = "https://github.com/AlanFreed/MutableTypes.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalSystemsOfUnits.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalFields.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalScalars.jl")
Pkg.add(url = "https://github.com/AlanFreed/PhysicalVectors.jl")
```

## Re-exported from PhysicalFields.jl

Package [Reexport.jl](https://github.com/simonster/Reexport.jl), or [here](https://juliapackages.com/p/reexport), is used to export symbols from an imported package. Here `Reexport` is used to re-export part of the interface of package [PhsysicalFields.jl](https://github.com/AlanFreed/PhysicalFields.jl), effectively linking these exports so they appear to originate from within this module. Specifically, using the alias

```
const PhysicalUnits = PhysicalSystemsOfUnits.PhysicalSystemOfUnits
```

the following type definition for a physical vector field is exported as

```
struct PhysicalVector <: PhysicalField
    l::UInt16           # length of the vector
    v::StaticVector     # values of the vector in its specified system of units
    u::PhysicalUnits    # the vector's physical units
end
```

which is a sub-type to the abstract type

```
abstract type PhysicalField end
```

Also, a type definition for an array of such vector fields is exported as

```
struct ArrayOfPhysicalVectors
    e::UInt16           # number of entries or elements held in the array
    l::UInt16           # length of each physical vector held in the array
    a::Array            # array of vectors holding values of a physical vector
    u::PhysicalUnits    # units of this physical vector
end
```

where all entries in the array have the same physical units.

Constructors for these types are also re-exported here, they being

```
function newPhysicalVector(len::Integer, units::PhysicalUnits)::PhysicalVector
```

which supplies a new vector of dimension `len` whose values are set to `0` and whose physical units are those supplied by the argument `units`, and

```
function newArrayOfPhysicalVectors(len::Integer, v₁::PhysicalVector)::ArrayOfPhysicalVectors
```

where `v₁` is the first entry in a new array of vectors whose length is `len`.

To retrieve and assign array values, functions `Base.:(getindex)` and `Base.:(setindex!)` have been overloaded so that the bracket notation `[]` can be used to: *i)* retrieve and assign scalar fields belonging to an instance of `PhysicalVector`, and *ii)* retrieve and assign vector fields belonging to an instance of `ArrayOfPhysicalVectors`. 

Also, conversion to a string is provided for instances of `PhysicalVector` by the re-exported method

```
function toString(v::PhysicalVector; format::Char='E')::String
```

where the keyword `format` is a character that, whenever its value is 'E' or 'e', represents the vector components in a scientific notation; otherwise, they will be represented in a fixed-point notation.

## Operators

The following operators have been overloaded so that they can handle objects of type PhysicalVector, whenever such operations make sense, e.g., one cannot add two vectors with different units or different dimensions. The overloaded logical operators include: `==`, `≠` and `≈`. The overloaded unary operators include: `+` and `-`. And the overloaded binary operators include: `+`, `-`, `*` and `/`.

## Methods for both PhysicalVector and ArrayOfPhysicalVectors

The following methods can accept arguments that are objects of either type, viz., PhysicalVector or ArrayOfPhysicalVectors. They are self explanatory: `copy`, `deepcopy`, `isDimensionless`, `isCGS`, `isSI`, `toCGS` and `toSI`.

## Math functions for PhysicalVector

Function `norm(v,p)` returns the p-norm of vector v, with the default being the Euclidean norm, i.e., p = 2. Function `unitVector(v)` returns a dimensionless vector of length 1 pointing in the direction of vector v. And function `cross(y,z)` returns the cross product y × z, provided the vectors have dimensions of 3.
