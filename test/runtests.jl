module runtests

using
    PhysicalSystemsOfUnits,
    PhysicalScalars,
    PhysicalVectors,
    Test

# Create three scalars with units of force.
s1 = newPhysicalScalar(DYNE)
s2 = newPhysicalScalar(DYNE)
s3 = newPhysicalScalar(DYNE)
set!(s1, 1.0)
set!(s2, 2.0)
set!(s3, 3.0)

# Create two vectors and poplulate them with components whose units are force.
a = newPhysicalVector(3, DYNE)
a[1] = s1
a[2] = s2
a[3] = s3
b = newPhysicalVector(3, DYNE)
b[1] = s1 * s2 / s3
b[2] = s2 * s3 / s1
b[3] = s3 * s1 / s2

@testset "[]" begin
    @test b[1] == a[1] * a[2] / a[3]
    @test b[2] == a[2] * a[3] / a[1]
    @test b[3] == a[3] * a[1] / a[2]
end

@testset "toString" begin
    @test PhysicalVectors.toString(a) == "{ 1.0000E+00  2.0000E+00  3.0000E+00}ᵀ g⋅cm/s²"
    @test PhysicalVectors.toString(b) == "{ 6.6667E-01  6.0000E+00  1.5000E+00}ᵀ g⋅cm/s²"
end

@testset "isDimensionless" begin
    @test PhysicalVectors.isDimensionless(a) == false
    @test PhysicalVectors.isDimensionless(newPhysicalVector(2, SI_DIMENSIONLESS)) == true
end

@testset "isCGS" begin
    @test PhysicalVectors.isCGS(newPhysicalVector(2, DYNE)) == true
    @test PhysicalVectors.isCGS(newPhysicalVector(2, PASCAL)) == false
    @test PhysicalVectors.isCGS(newPhysicalVector(2, CGS_DIMENSIONLESS)) == true
    @test PhysicalVectors.isCGS(newPhysicalVector(2, SI_DIMENSIONLESS)) == false
end

@testset "isSI" begin
    @test PhysicalVectors.isSI(newPhysicalVector(2, DYNE)) == false
    @test PhysicalVectors.isSI(newPhysicalVector(2, PASCAL)) == true
    @test PhysicalVectors.isSI(newPhysicalVector(2, CGS_DIMENSIONLESS)) == false
    @test PhysicalVectors.isSI(newPhysicalVector(2, SI_DIMENSIONLESS)) == true
end

@testset "toCGS & toSI" begin
    @test PhysicalVectors.toCGS(PhysicalVectors.toSI(a)) ≈ a
    c = PhysicalVectors.toSI(b)
    @test PhysicalVectors.toSI(PhysicalVectors.toCGS(c)) ≈ c
end

@testset "norm" begin
    s = newPhysicalScalar(DYNE)
    set!(s, sqrt(14))
    @test norm(a) ≈ s
end

@testset "unitVector" begin
    s = newPhysicalScalar(CGS_DIMENSIONLESS)
    set!(s, 1)
    u = newPhysicalVector(3, CGS_DIMENSIONLESS)
    u[1] = (1 / sqrt(14)) * s
    u[2] = (2 / sqrt(14)) * s
    u[3] = (3 / sqrt(14)) * s
    @test unitVector(a) ≈ u
end

@testset "cross" begin
    one = newPhysicalScalar(SI_LENGTH)
    set!(one, 1)
    u = newPhysicalVector(3, SI_LENGTH)
    v = newPhysicalVector(3, SI_LENGTH)
    w = newPhysicalVector(3, SI_AREA)
    u[1] = one
    v[2] = one
    w[3] = one * one
    @test cross(u, v) ≈ w
end

end # module runtests
