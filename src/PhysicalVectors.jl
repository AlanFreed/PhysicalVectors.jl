module PhysicalVectors

using
    LinearAlgebra,
    MutableTypes,
    PhysicalScalars,
    PhysicalSystemsOfUnits,
    Reexport,
    StaticArrays

@reexport using PhysicalFields:
    PhysicalVector,
    ArrayOfPhysicalVectors,
    newPhysicalVector,
    newArrayOfPhysicalVectors,
    getindex,
    setindex!,
    VtoString

export
    # For both PhysicalVector and ArrayOfPhysicalVectors objects:
    isDimensionless,
    isCGS,
    isSI,
    toCGS,
    toSI,
    copy,
    deepcopy,
    # For PhysicalVector objects:
    toString,
    toArray,
    # math functions
    norm,
    unitVector,
    cross

#=
--------------------------------------------------------------------------------
=#

const MNumber = MutableTypes.MNumber

#=
--------------------------------------------------------------------------------
=#

# Methods for testing the kind of system of units.

function isDimensionless(v::PhysicalVector)::Bool
    return PhysicalSystemsOfUnits.isDimensionless(v.u)
end

function isDimensionless(av::ArrayOfPhysicalVectors)::Bool
    return PhysicalSystemsOfUnits.isDimensionless(av.u)
end

function isCGS(v::PhysicalVector)::Bool
    return PhysicalSystemsOfUnits.isCGS(v.u)
end

function isCGS(av::ArrayOfPhysicalVectors)::Bool
    return PhysicalSystemsOfUnits.isCGS(av.u)
end

function isSI(v::PhysicalVector)::Bool
    return PhysicalSystemsOfUnits.isSI(v.u)
end

function isSI(av::ArrayOfPhysicalVectors)::Bool
    return PhysicalSystemsOfUnits.isSI(av.u)
end

# Methods for converting between systems of units.

function toCGS(v::PhysicalVector)::PhysicalVector
    if isCGS(v)
        return v
    elseif isSI(v)
        units = PhysicalSystemsOfUnits.CGS(v.u.m, v.u.kg, v.u.s, v.u.K)
        vector = newPhysicalVector(v.l, units)
        if (v.u == PhysicalSystemsOfUnits.KELVIN)
            for i in 1:v.l
                vector.v[i] = v.v[i] - 273.15
            end
        else
            for i in 1:v.l
                vector.v[i] = (100.0^v.u.m * 1000.0^v.u.kg) * v.v[i]
            end
        end
        return vector
    else
        msg = "Units for vector 'v' must be either CGS or SI."
        throw(ErrorException(msg))
    end
end

function toCGS(av::ArrayOfPhysicalVectors)::ArrayOfPhysicalVectors
    if isCGS(av)
        return av
    elseif isSI(av)
        units = PhysicalSystemsOfUnits.CGS(av.u.m, av.u.kg, av.u.s, av.u.K)
        pv??? = newPhysicalVector(av.l, units)
        if (av.u == PhysicalSystemsOfUnits.KELVIN)
            for i in 1:av.l
                pv???.v[i] = av.a[1,i] - 273.15
            end
        else
            for i in 1:av.l
                pv???.v[i] = (100.0^av.u.m * 1000.0^av.u.kg) * av.a[1,i]
            end
        end
        vecArr = newArrayOfPhysicalVectors(av.e, pv???)
        for i in 2:av.e
            if (av.u == PhysicalSystemsOfUnits.KELVIN)
                for j in 1:av.l
                    vecArr.a[i,j] = av.a[i,j] - 273.15
                end
            else
                for j in 1:av.l
                    vecArr.a[i,j] = (100.0^av.u.m * 1000.0^av.u.kg) * av.a[i,j]
                end
            end
        end
        return vecArr
    else
        msg = "Units for the array of vectors 'av' must be either CGS or SI."
        throw(ErrorException(msg))
    end
end

function toSI(v::PhysicalVector)::PhysicalVector
    if isSI(v)
        return v
    elseif isCGS(v)
        units = PhysicalSystemsOfUnits.SI(v.u.cm, v.u.g, v.u.s, v.u.C)
        vector = newPhysicalVector(v.l, units)
        if (v.u == PhysicalSystemsOfUnits.CENTIGRADE)
            for i in 1:v.l
                vector.v[i] = 273.15 + v.v[i]
            end
        else
            for i in 1:v.l
                vector.v[i] = (100.0^(-v.u.cm) * 1000.0^(-v.u.g)) * v.v[i]
            end
        end
        return vector
    else
        msg = "Units for vector 'v' must be either CGS or SI."
        throw(ErrorException(msg))
    end
end

function toSI(av::ArrayOfPhysicalVectors)::ArrayOfPhysicalVectors
    if isSI(av)
        return av
    elseif isCGS(av)
        units = PhysicalSystemsOfUnits.SI(av.u.cm, av.u.g, av.u.s, av.u.C)
        pv??? = newPhysicalVector(av.l, units)
        if (av.u == PhysicalSystemsOfUnits.CENTIGRADE)
            for i in 1:av.l
                pv???.v[i] = 273.15 + av.a[1,i]
            end
        else
            for i in 1:av.l
                pv???.v[i] = (100.0^(-av.u.cm) * 1000.0^(-av.u.g)) * av.a[1,i]
            end
        end
        vecArr = newArrayOfPhysicalVectors(av.e, pv???)
        for i in 2:av.e
            if (av.u == PhysicalSystemsOfUnits.CENTIGRADE)
                for j in 1:av.l
                    vecArr.a[i,j] = 273.15 + av.a[i,j]
                end
            else
                for j in 1:av.l
                    vecArr.a[i,j] = (100.0^(-av.u.cm) * 1000.0^(-av.u.g)) * av.a[i,j]
                end
            end
        end
        return vecArr
    else
        msg = "Units for the array of vectors 'av' must be either CGS or SI."
        throw(ErrorException(msg))
    end
end

#=
--------------------------------------------------------------------------------
=#

function Base.:(copy)(y::PhysicalVector)::PhysicalVector
    length = copy(y.l)
    vector = copy(y.v)
    units  = copy(y.u)
    return PhysicalVector(UInt16(length), vector, units)
end

function Base.:(copy)(y::ArrayOfPhysicalVectors)::ArrayOfPhysicalVectors
    entries = copy(y.e)
    length  = copy(y.l)
    array   = copy(y.a)
    units   = copy(y.u)
    return ArrayOfPhysicalVectors(UInt16(entries), UInt16(length), array, units)
end

function Base.:(deepcopy)(y::PhysicalVector)::PhysicalVector
    length = deepcopy(y.l)
    vector = deepcopy(y.v)
    units  = deepcopy(y.u)
    return PhysicalVector(UInt16(length), vector, units)
end

function Base.:(deepcopy)(y::ArrayOfPhysicalVectors)::ArrayOfPhysicalVectors
    entries = deepcopy(y.e)
    length  = deepcopy(y.l)
    array   = deepcopy(y.a)
    units   = deepcopy(y.u)
    return ArrayOfPhysicalVectors(UInt16(entries), UInt16(length), array, units)
end

#=
--------------------------------------------------------------------------------
=#

# Overloaded these operators: logical: ==, ???, ???
#                             unary:   +, -
#                             binary:  +, -, *, /

function Base.:???(y::PhysicalVector, z::PhysicalVector)::Bool
    if y.l ??? z.l
        return true
    end
    if PhysicalSystemsOfUnits.isEquivalent(y.u, z.u)
        if (isCGS(y) && isCGS(z)) || (isSI(y) && isSI(z))
            for i in 1:y.l
                if y.v[i] ??? z.v[i]
                    return true
                end
            end
        elseif (isCGS(y) && isSI(z))
            w = toCGS(z)
            for i in 1:y.l
                if y.v[i] ??? w.v[i]
                    return true
                end
            end
        elseif (isSI(y) && isCGS(z))
            w = toCGS(y)
            for i in 1:w.l
                if w.v[i] ??? z.v[i]
                    return true
                end
            end
        else
            msg = "Vectors must have units of either CGS or SI."
            throw(ErrorException(msg))
        end
    else
        return true
    end
    return false
end

function Base.:(==)(y::PhysicalVector, z::PhysicalVector)::Bool
    notEqual = (y ??? z)
    return !notEqual
end

function Base.:???(y::PhysicalVector, z::PhysicalVector)::Bool
    if y.l ??? z.l
        return false
    end
    if PhysicalSystemsOfUnits.isEquivalent(y.u, z.u)
        if (isCGS(y) && isCGS(z)) || (isSI(y) && isSI(z))
            for i in 1:y.l
                if !(y.v[i] ??? z.v[i])
                    return false
                end
            end
        elseif (isCGS(y) && isSI(z))
            w = toCGS(z)
            for i in 1:y.l
                if !(y.v[i] ??? w.v[i])
                    return false
                end
            end
        elseif (isSI(y) && isCGS(z))
            w = toCGS(y)
            for i in 1:w.l
                if !(w.v[i] ??? z.v[i])
                    return false
                end
            end
        else
            msg = "Vectors must have units of either CGS or SI."
            throw(ErrorException(msg))
        end
    else
        return false
    end
    return true
end

function Base.:+(y::PhysicalVector)::PhysicalVector
    vector = newPhysicalVector(y.l, y.u)
    for i in 1:y.l
        vector[i] = +y[i]
    end
    return vector
end

function Base.:-(y::PhysicalVector)::PhysicalVector
    vector = newPhysicalVector(y.l, y.u)
    for i in 1:y.l
        vector[i] = -y[i]
    end
    return vector
end

function Base.:+(y::PhysicalVector, z::PhysicalVector)::PhysicalVector
    if y.l ??? z.l
        msg = "Vector addition requires vectors to have the same length."
        throw(DimensionMismatch(msg))
    end
    if PhysicalSystemsOfUnits.isEquivalent(y.u, z.u)
        if (isCGS(y) && isCGS(z)) || (isSI(y) && isSI(z))
            vector = newPhysicalVector(y.l, y.u)
            for i in 1:y.l
                vector[i] = y[i] + z[i]
            end
        elseif (isCGS(y) && isSI(z))
            w = toCGS(z)
            vector = newPhysicalVector(y.l, y.u)
            for i in 1:y.l
                vector[i] = y[i] + w[i]
            end
        elseif (isSI(y) && isCGS(z))
            w = toCGS(y)
            vector = newPhysicalVector(z.l, z.u)
            for i in 1:w.l
                vector[i] = w[i] + z[i]
            end
        else
            msg = "Vectors must have either CGS or SI units."
            throw(ErrorException(msg))
        end
    else
        msg = "Vector addition requires vectors to have equivalent units."
        throw(ErrorException(msg))
    end
    return vector
end

function Base.:-(y::PhysicalVector, z::PhysicalVector)::PhysicalVector
    if y.l ??? z.l
        msg = "Vector subtraction requires vectors to have the same length."
        throw(DimensionMismatch(msg))
    end
    if PhysicalSystemsOfUnits.isEquivalent(y.u, z.u)
        if (isCGS(y) && isCGS(z)) || (isSI(y) && isSI(z))
            vector = newPhysicalVector(y.l, y.u)
            for i in 1:y.l
                vector[i] = y[i] - z[i]
            end
        elseif (isCGS(y) && isSI(z))
            w = toCGS(z)
            vector = newPhysicalVector(y.l, y.u)
            for i in 1:y.l
                vector[i] = y[i] - w[i]
            end
        elseif (isSI(y) && isCGS(z))
            w = toCGS(y)
            vector = newPhysicalVector(z.l, z.u)
            for i in 1:w.l
                vector[i] = w[i] - z[i]
            end
        else
            msg = "Vectors must have either CGS or SI units."
            throw(ErrorException(msg))
        end
    else
        msg = "Vector subtraction requires vectors to have equivalent units."
        throw(ErrorException(msg))
    end
    return vector
end

function Base.:*(y::PhysicalVector, z::PhysicalVector)::PhysicalScalar
    if y.l ??? z.l
        msg = "A vector dot product requires vectors to have the same length."
        throw(DimensionMismatch(msg))
    end
    if (isCGS(y) && isCGS(z)) || (isSI(y) && isSI(z))
        units = y.u + z.u
        value = LinearAlgebra.dot(y.v, z.v)
    elseif (isCGS(y) && isSI(z))
        w = toCGS(z)
        units = y.u + w.u
        value = LinearAlgebra.dot(y.v, w.v)
    elseif (isSI(y) && isCGS(z))
        w = toCGS(y)
        units = w.u + z.u
        value = LinearAlgebra.dot(w.v, z.v)
    else
        msg = "Vectors must have either CGS or SI units."
        throw(ErrorException(msg))
    end
    dotProduct = newPhysicalScalar(units)
    PhysicalScalars.set!(dotProduct, value)
    return dotProduct
end

function Base.:*(y::PhysicalScalar, z::PhysicalVector)::PhysicalVector
    if (PhysicalScalars.isCGS(y) && isCGS(z)) || (PhysicalScalars.isSI(y) && isSI(z))
        units = y.u + z.u
        scalarProduct = newPhysicalVector(z.l, units)
        for i in 1:z.l
            scalarProduct[i] = y * z[i]
        end
    elseif (PhysicalScalars.isCGS(y) && isSI(z))
        w = toCGS(z)
        units = y.u + w.u
        scalarProduct = newPhysicalVector(w.l, units)
        for i in 1:w.l
            scalarProduct[i] = y * w[i]
        end
    elseif (PhysicalScalars.isSI(y) && isCGS(z))
        w = PhysicalScalars.toCGS(y)
        units = w.u + z.u
        scalarProduct = newPhysicalVector(z.l, units)
        for i in 1:z.l
            scalarProduct[i] = w * z[i]
        end
    else
        msg = "Scalars and vectors must have either CGS or SI units."
        throw(ErrorException(msg))
    end
    return scalarProduct
end

function Base.:*(y::Union{Real,MNumber}, z::PhysicalVector)::PhysicalVector
    scalarProduct = newPhysicalVector(z.l, z.u)
    for i in 1:z.l
        scalarProduct[i] = y * z[i]
    end
    return scalarProduct
end

function Base.:/(y::PhysicalVector, z::PhysicalScalar)::PhysicalVector
    if (isCGS(y) && PhysicalScalars.isCGS(z)) || (isSI(y) && PhysicalScalars.isSI(z))
        units = y.u - z.u
        scalarDivision = newPhysicalVector(y.l, units)
        for i in 1:y.l
            scalarDivision[i] = y[i] / z
        end
    elseif (isCGS(y) && PhysicalScalars.isSI(z))
        w = PhysicalScalars.toCGS(z)
        units = y.u - w.u
        scalarDivision = newPhysicalVector(y.l, units)
        for i in 1:y.l
            scalarDivision[i] = y[i] / w
        end
    elseif (isSI(y) && PhysicalScalars.isCGS(z))
        w = toCGS(y)
        units = w.u - z.u
        scalarDivision = newPhysicalVector(w.l, units)
        for i in 1:w.l
            scalarDivision[i] = w[i] / z
        end
    else
        msg = "Scalars and vectors must have either CGS or SI units."
        throw(ErrorException(msg))
    end
    return scalarDivision
end

function Base.:/(y::PhysicalVector, z::Union{Real,MNumber})::PhysicalVector
    scalarDivision = newPhysicalVector(y.l, y.u)
    for i in 1:y.l
        scalarDivision[i] = y[i] / z
    end
    return scalarDivision
end

#=
--------------------------------------------------------------------------------
=#

# Functions of type PhysicalVector:

function toString(v::PhysicalVector; format::Char='E')::String
    return VtoString(v; format)
end

function toArray(v::PhysicalVector)::StaticVector
    return deepcopy(v.v)
end

function LinearAlgebra.:(norm)(y::PhysicalVector, p::Real=2)::PhysicalScalar
    value = norm(y.v, p)
    units = y.u
    magnitude = newPhysicalScalar(units)
    PhysicalScalars.set!(magnitude, value)
    return magnitude
end

function unitVector(y::PhysicalVector)::PhysicalVector
    unitVec = y / norm(y,2)
    return unitVec
end

function LinearAlgebra.:(cross)(y::PhysicalVector, z::PhysicalVector)::PhysicalVector
    if (y.l ??? 3) || (z.l ??? 3)
        msg = "Vector cross product is only defined for 3 dimensional vectors."
        throw(DimensionMismatch(msg))
    end
    if (isCGS(y) && isCGS(z)) || (isSI(y) && isSI(z))
        units = y.u + z.u
        value = cross(y.v, z.v)
    elseif (isSI(y) && isCGS(z))
        x = toCGS(y)
        units = x.u + z.u
        value = cross(x.v, z.v)
    elseif (isCGS(y) && isSI(z))
        x = toCGS(z)
        units = y.u + x.u
        value = cross(y.v, x.v)
    else
        msg = "Vectors must have either CGS or SI units."
        throw(ErrorException(msg))
    end
    crossProduct = newPhysicalVector(3, units)
    for i in 1:3
        scalar = newPhysicalScalar(units)
        PhysicalScalars.set!(scalar, value[i])
        crossProduct[i] = scalar
    end
    return crossProduct
end

end # module PhysicalVectors
