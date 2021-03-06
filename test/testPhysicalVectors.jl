
module testPhysicalVectors

using
    Base,
    LinearAlgebra,
    PhysicalSystemsOfUnits,
    PhysicalScalars,
    PhysicalVectors

export
    run

const
    UtoString = PhysicalSystemsOfUnits.toString
    StoString = PhysicalScalars.toString
    VtoString = PhysicalVectors.toString

function run()
    format = 'F'
    precision = 5
    aligned = true
    println("Test the get and set operators via []:")
    a = newPhysicalVector(3, DYNE)
    println("A new vector: ", VtoString(a; format))
    a1 = newPhysicalScalar(a.u)
    set!(a1, 1.0)
    a2 = newPhysicalScalar(a.u)
    set!(a2, 2.0)
    a3 = newPhysicalScalar(a.u)
    set!(a3, 3.0)
    a[1] = a1
    a[2] = a2
    a[3] = a3
    println("reassigned:   ", VtoString(a; format))
    b     = newPhysicalVector(15, CENTIMETER)
    b1  = newPhysicalScalar(CENTIMETER)
    set!(b1, 1.0)
    b2  = newPhysicalScalar(CENTIMETER)
    set!(b2, 0.9)
    b3  = newPhysicalScalar(CENTIMETER)
    set!(b3, 0.8)
    b4  = newPhysicalScalar(CENTIMETER)
    set!(b4, 0.7)
    b5  = newPhysicalScalar(CENTIMETER)
    set!(b5, 0.6)
    b6  = newPhysicalScalar(CENTIMETER)
    set!(b6, 0.5)
    b7  = newPhysicalScalar(CENTIMETER)
    set!(b7, 0.4)
    b8  = newPhysicalScalar(CENTIMETER)
    set!(b8, 0.3)
    b9  = newPhysicalScalar(CENTIMETER)
    set!(b9, 0.2)
    b10 = newPhysicalScalar(CENTIMETER)
    set!(b10, 0.1)
    b11 = newPhysicalScalar(CENTIMETER)
    b12 = newPhysicalScalar(CENTIMETER)
    set!(b12, -0.1)
    b13 = newPhysicalScalar(CENTIMETER)
    set!(b13, -0.2)
    b14 = newPhysicalScalar(CENTIMETER)
    set!(b14, -0.3)
    b15 = newPhysicalScalar(CENTIMETER)
    set!(b15, -0.4)
    b[1] = b1
    b[2] = b2
    b[3] = b3
    b[4] = b4
    b[5] = b5
    b[6] = b6
    b[7] = b7
    b[8] = b8
    b[9] = b9
    b[10] = b10
    b[11] = b11
    b[12] = b12
    b[13] = b13
    b[14] = b14
    b[15] = b15
    println("Check printing of long vectors:")
    println(VtoString(b; format))
    format = 'E'
    println(VtoString(b; format))
    println("Check printing of intermediate length vectors:")
    c = newPhysicalVector(9, CENTIMETER)
    for i in 1:c.l
        c[i] = b[i]
    end
    format = 'F'
    println(VtoString(c; format))
    d = newPhysicalVector(6, CENTIMETER)
    for i in 1:d.l
        d[i] = b[i]
    end
    format = 'E'
    println(VtoString(d; format))
    println("Check printing of short vectors:")
    println(VtoString(a; format))
    format = 'F'
    println(VtoString(a; format))
    println("Testing vector arithmetic in 3 space:")
    y = newPhysicalScalar(CENTIMETER)
    set!(y, ??)
    b = newPhysicalVector(3, DYNE)
    b1 = newPhysicalScalar(DYNE)
    set!(b1, -3.0)
    b2 = newPhysicalScalar(DYNE)
    set!(b2, -2.0)
    b3 = newPhysicalScalar(DYNE)
    set!(b3, -1.0)
    b[1] = b1
    b[2] = b2
    b[3] = b3
    format = 'E'
    println("    y = ", StoString(y; format, precision, aligned))
    println("    a = ", VtoString(a; format))
    println("    b = ", VtoString(b; format))
    println("   -b = ", VtoString(-b; format))
    println("a + b = ", VtoString(a+b; format))
    println("a - b = ", VtoString(a-b; format))
    println("a * b = ", StoString(a*b; format, precision, aligned))
    println("y * a = ", VtoString(y*a; format))
    println("a / y = ", VtoString(a/y; format))
    println("||b|| = ", StoString(norm(b); format, precision, aligned))
    println("Base vectors:")
    e1 = unitVector(a)
    e2 = unitVector(b)
    println("e1    = ", VtoString(e1; format))
    println("e2    = ", VtoString(e2; format))
    e3 = cross(e1, e2)
    println("e3    = ", VtoString(e3; format), "  where e3 = e1 ?? e2")
    println("e1.e2 = ", StoString(e1*e2; format, precision, aligned), "  angle is obtuse")
    println("e1.e3 = ", StoString(e1*e3; format, precision, aligned), "  angle is right")
    println("e2.e3 = ", StoString(e2*e3; format, precision, aligned), "  angle is right")
    println()
    println("Test ArrayOfPhysicalVectors:")
    println()
    entries = 5
    len = 3
    v??? = newPhysicalVector(len, PASCAL)
    v???1 = newPhysicalScalar(PASCAL)
    set!(v???1, 1)
    v???2 = newPhysicalScalar(PASCAL)
    set!(v???2, 2)
    v???3 = newPhysicalScalar(PASCAL)
    set!(v???3, 3)
    v???[1] = v???1
    v???[2] = v???2
    v???[3] = v???3
    a = newArrayOfPhysicalVectors(entries, v???)
    n = 3
    for i in 2:entries
        v??? = newPhysicalVector(len, PASCAL)
        n += 1
        v???1 = newPhysicalScalar(PASCAL)
        set!(v???1, n)
        n += 1
        v???2 = newPhysicalScalar(PASCAL)
        set!(v???2, n)
        n += 1
        v???3 = newPhysicalScalar(PASCAL)
        set!(v???3, n)
        v???[1] = v???1
        v???[2] = v???2
        v???[3] = v???3
        a[i] = v???
    end
    println("This array of vectors has a length of ", string(a.e), ".")
    for i in 1:entries
        v??? = a[i]
        println("a[", string(i), "] = ", VtoString(v???; format))
    end
    a[3] = newPhysicalVector(len, PASCAL)
    println("resetting the third entry to zeros, one has")
    println("a[3] = ", VtoString(a[3]; format))
    println()
    println("If these answers make sense, then this test passes.")
    return nothing
end

end  # module testPhysicalVectors
