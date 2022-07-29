module runtests

using
    MutableTypes,
    PhysicalSystemsOfUnits,
    PhysicalScalars,
    PhysicalVectors,
    PhysicalTensors,
    Test

# Create a force vector.
f = newPhysicalVector(3, NEWTON)

# Create and populate a normal vector with units of area.
n = newPhysicalVector(3, SI_AREA)
n[1] = PhysicalScalar(MReal(0.6), SI_AREA)
n[1] = PhysicalScalar(MReal(0.8), SI_AREA)
n[1] = PhysicalScalar(MReal(0.0), SI_AREA)

# Create and populate a tensor with units of stress.
σ = newPhysicalTensor(3, 3, PASCAL)
σ[1,1] = PhysicalScalar(MReal(1), PASCAL)
σ[1,2] = PhysicalScalar(MReal(2), PASCAL)
σ[1,3] = PhysicalScalar(MReal(3), PASCAL)
σ[2,1] = PhysicalScalar(MReal(4), PASCAL)
σ[2,2] = PhysicalScalar(MReal(5), PASCAL)
σ[2,3] = PhysicalScalar(MReal(6), PASCAL)
σ[3,1] = PhysicalScalar(MReal(7), PASCAL)
σ[3,2] = PhysicalScalar(MReal(8), PASCAL)
σ[3,3] = PhysicalScalar(MReal(9), PASCAL)

# Create and populate force vector.
f = newPhysicalVector(3, NEWTON)
for i in 1:3
    for j in 1:3
        f[i] = σ[i,j] * n[j]
    end
end

@testset "[]" begin
    force = σ * n  
    @test force[1] == f[1]
    @test force[2] == f[2]
    @test force[3] == f[3]
end

@testset "toString" begin
    @test PhysicalTensors.toString(σ) == "⌈ 1.0000E+00  2.0000E+00  3.0000E+00⌉
| 4.0000E+00  5.0000E+00  6.0000E+00| kg/(m⋅s²)
⌊ 7.0000E+00  8.0000E+00  9.0000E+00⌋"
end

@testset "isDimensionless" begin
    @test PhysicalTensors.isDimensionless(σ) == false
    @test PhysicalTensors.isDimensionless(newPhysicalTensor(2, 2, SI_DIMENSIONLESS)) == true
end

@testset "isCGS" begin
    @test PhysicalTensors.isCGS(newPhysicalTensor(2, 2, DYNE)) == true
    @test PhysicalTensors.isCGS(newPhysicalTensor(2, 2, PASCAL)) == false
    @test PhysicalTensors.isCGS(newPhysicalTensor(2, 2, CGS_DIMENSIONLESS)) == true
    @test PhysicalTensors.isCGS(newPhysicalTensor(2, 2, SI_DIMENSIONLESS)) == false
end

@testset "isSI" begin
    @test PhysicalTensors.isSI(newPhysicalTensor(2, 2, DYNE)) == false
    @test PhysicalTensors.isSI(newPhysicalTensor(2, 2, PASCAL)) == true
    @test PhysicalTensors.isSI(newPhysicalTensor(2, 2, CGS_DIMENSIONLESS)) == false
    @test PhysicalTensors.isSI(newPhysicalTensor(2, 2, SI_DIMENSIONLESS)) == true
end

@testset "toCGS & toSI" begin
    @test PhysicalTensors.toSI(PhysicalTensors.toCGS(σ)) ≈ σ
end

@testset "copy & deepcopy" begin
    @test PhysicalTensors.copy(σ) == σ
    @test PhysicalTensors.deepcopy(σ) == σ
end

@testset "norm" begin
    @test norm(σ) ≈ PhysicalScalar(MReal(sqrt(285)), PASCAL)
end

@testset "tensorProduct" begin
    fn = newPhysicalTensor(3, 3, f.u+n.u)
    for i in 1:3
        for j in 1:3
            fn[i,j] = f[i] * n[j]
        end
    end
    @test tensorProduct(f,n) ≈ fn
end

end # module runtests
