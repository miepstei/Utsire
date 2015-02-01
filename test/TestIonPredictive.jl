using FactCheck

function setup()
    nopen = 2
    k = 5

    tres = 5e-5

    Q = [ [ -3050 50           0      3000   0; ]
          [ 2/3   -(500 + 2/3) 500    0      0; ]
          [ 0     15000        -19000 4000   0; ]
          [ 15    0            50     -2065  2000;]
          [ 0     0            0      10     -10;] ]

    #standard dcp options
    dcProgsOptions = [2 1e-12 1e-12 100 -1e6 0] 

    t = [1e-5 2e-5 3e-5 4e-5 5e-5 6e-5 7e-5 8e-5 9e-5 1e-4]

    data = Dict()
    data["Q"] = Q
    data["tres"] = tres
    data["nopen"] = nopen
    data["k"] = k
    data["dcProgsOptions"] = dcProgsOptions
    data["t"] = t
    data["tcrit"] = 1e-4

    return data
end

function TestIonPredictive()
    #setup code to match DCPROGS
    println("Running tests from CHS 1996...")
    data = setup()
    TestCHSOccupancies(data)
    TestExactSurvivorRecursiveMatrix(data)
    TestExactSurvivorXt(data)
    TestDetWs(data)
    TestDerivDetWs(data)
    TestSurvivorXs(data)
    TestDerivSurvivorXs(data)
    TestApproxSurvivorComponent(data)
    TestFindRoots(data)
    TestOccupancies(data)
    TestIdealAFt(data)
    TestIdealFAt(data)
end

function TestCHSOccupancies(data::Dict)
    occ = zeros(2,)

    #int dcpCHSOccupancies ( double * phi, double * Q, int nopen, int k, double tres, double tcrit, bool initial, double * dcpOptions)
    x = ccall( (:dcpCHSOccupancies,"libdcprogs_wrapper") , Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 , Int64, Float64, Float64, Bool, Ptr{Float64}) , occ , data["Q"], data["nopen"], data["k"], data["tres"], data["tcrit"], true, data["dcProgsOptions"])
   
    facts("Testing CHS occupancies") do
        @fact x => 0
        @fact occ => [0.12 0.88]
    end
    
end


function TestExactSurvivorRecursiveMatrix(data::Dict)
    recursivematrix = zeros(data["nopen"],data["nopen"])

    #( double * recursiveMatrix , double * Q, int nopen, int k, double tres, bool isOpen, int i, int m, int r)

    x = ccall( (:dcpExactSurvivorRecursiveMatrix , "libdcprogs_wrapper") , Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 ,Int64, Float64, Bool, Int64, Int64, Int64) , recursivematrix, data["Q"] , data["nopen"] , data["k"] , data["tres"] , true, 0, 0, 0)

    facts("Testing ExactSurvivor Recursive Matrices") do
        @fact x => 0
        @pending recursivematrix => [2.3 4.5; 5.3 4.1]
    end
end

function TestExactSurvivorXt(data::Dict)
    survivor = zeros(data["nopen"],data["nopen"])
    #( double * survivor, double * Q, int nopen, int k, double tres, double t, bool isOpen)

    x = ccall( (:dcpExactSurvivorXt , "libdcprogs_wrapper") , Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 ,Int64, Float64, Float64, Bool) , survivor, data["Q"] , data["nopen"] , data["k"] , data["tres"] , data["t"][10], true)

    facts("Checking Exact Survivor Function") do
        @fact x => 0
        @pending survivor => [0.0035  0.02; 0.3  0.76] 
    end
end


function TestFindRoots(data::Dict)

    #dcpFindRoots( double * roots, double * Q, int nopen, int k, double tres, bool isOpen, double * dcpOptions )
    roots = zeros(2)
    x = ccall( (:dcpFindRoots , "libdcprogs_wrapper") , Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 ,Int64, Float64, Bool, Ptr{Float64}) ,  roots, data["Q"] , data["nopen"] , data["k"] , data["tres"] , true, data["dcProgsOptions"])
    facts("Checking Open Roots") do
        @fact x => 0
        @pending roots => [0.0035 ; 0.02] 
    end
end

function TestApproxSurvivorComponent(data::Dict)
    #dcpApproxSurvivorComponent( double * survivor, double * root, double * Q, int nopen, int k, double tres, int ith, bool isOpen, double * dcpOptions)
    survivorcomponent = zeros(data["nopen"] , data["nopen"])
    root = zeros(1)
    x = ccall( (:dcpApproxSurvivorComponent , "libdcprogs_wrapper") , Int64 , (Ptr{Float64} , Ptr{Float64} , Ptr{Float64} , Int64 ,Int64, Float64, Int64, Bool, Ptr{Float64}) , survivorcomponent , root, data["Q"] , data["nopen"] , data["k"] , data["tres"] , 0, true, data["dcProgsOptions"])
    facts("Checking Component of Survivor Function AR") do
        @fact x => 0
        @pending root[1] => 0.0035 
        @pending survivorcomponent => [20 10; 14 17]
    end
end

function TestDerivSurvivorXs(data::Dict)
    #dcpDerivSurvivorXs( double * survivor, double * Q, int nopen, int k, double tres, double s, bool isOpen )
    derivsurvivor = zeros(data["nopen"] , data["nopen"])
    x = ccall((:dcpDerivSurvivorXs , "libdcprogs_wrapper"), Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 , Int64 , Float64, Float64, Bool) , derivsurvivor , data["Q"], data["nopen"] , data["k"], data["tres"] , -0.1, true )

    facts("Checking Derivative of SurvivorXs (dAR(s) / ds)") do
        @fact x => 0
        @pending survivor => [20 10; 14 17]
    end
end

function TestSurvivorXs(data::Dict)
    #dcpSurvivorXs ( double * survivor, double * Q, int nopen, int k, double tres, double s, bool isOpen)
    survivor = zeros(data["nopen"] , data["nopen"])
    x = ccall((:dcpSurvivorXs , "libdcprogs_wrapper"), Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 , Int64 , Float64, Float64, Bool) , survivor , data["Q"], data["nopen"] , data["k"], data["tres"] , -0.1, true )

    facts("Checking SurvivorXs (AR(s))") do
        @fact x => 0
        @pending survivor => [20 10; 14 17]
    end

end

function TestDetWs(data::Dict)
    #int dcpDetWs ( double * determinant , double * Q, int nopen, int k, double tres, double s, bool isOpen )

    det = zeros(1)
   #x = ccall((: , "libdcprogs_wrapper"), Int64 , ()  )

    x = ccall((:dcpDetWs , "libdcprogs_wrapper"), Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 , Int64 , Float64, Float64, Bool) , det , data["Q"], data["nopen"] , data["k"], data["tres"] , -0.1, true)

    facts("Checking DetWs (Determinant of W(s))") do
        @fact x => 0
        @pending det[1] => 50
    end
end

function TestDerivDetWs(data::Dict)

# int dcpDerivDetWs ( double * deriv , double * Q, int nopen, int k, double tres, double s, bool isOpen )
    detderiv=zeros(data["nopen"] , data["nopen"] )

    x = ccall((:dcpDerivDetWs , "libdcprogs_wrapper"), Int64 , (Ptr{Float64} , Ptr{Float64} , Int64 , Int64 , Float64 , Float64 , Bool), detderiv , data["Q"], data["nopen"] , data["k"], data["tres"] , -10.0, true )

    facts("Checking the derivative of DetWs") do
        @fact x => 0   
        @pending detderiv => [20 10; 14 17]
    end

end


function TestOccupancies(data::Dict)
   #dcpOccupancies ( double * phi, double * Q, size_t nopen, size_t k, double tau, bool initial, double * dcpOptions )

   phiA = zeros(1, data["nopen"])
   x=ccall((:dcpOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Int64,Int64,Float64,Bool, Ptr{Float64}), phiA, data["Q"] , data["nopen"] , data["k"], data["tres"], true, data["dcProgsOptions"] )
   
   phiF = zeros(data["k"] - data["nopen"])
   x=ccall((:dcpOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Int64,Int64,Float64,Bool, Ptr{Float64}), phiF, data["Q"] , data["nopen"] , data["k"], data["tres"], false, data["dcProgsOptions"] )

   facts("Checking Occupancies") do
       context("Checking phiA") do
           @fact size(phiA) => (1,2)
           @fact phiA => roughly([0.1187288113282500 0.8812711886717499],atol=1e-12)
       end

       context("Checking phiF") do
           @fact size(phiF) => (3,)
           @fact phiF =>  roughly([0.6610020847564215; 0.3153887038404265; 0.0236092114031523],atol=1e-12)
       end
   end
end

function TestIdealAFt(data)
    AFt = zeros(data["nopen"] , data["k"] - data["nopen"])
    x=ccall((:dcpIdealGXYt,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Int64,Int64,Float64,Bool ), AFt, data["Q"] , 2 , 5, data["t"][10], true )

    facts("Checking AF(t[10])") do
        @fact size(AFt) => (2 , 3)
        @fact AFt => roughly([2.099006101003331 2211.3705251584956 0.0;
                        475.5830785962392 0.16792048808026652 0.0],atol=1e-12)
    end

end

function TestIdealFAt(data)
    FAt = zeros(data["k"] - data["nopen"], data["nopen"])
    x=ccall((:dcpIdealGXYt,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Int64,Int64,Float64,Bool ), FAt, data["Q"] , 2 , 5, data["t"][10], false )
    facts("Checking FA(t[10])") do
        @fact size(FAt) => (3 , 2)
        @fact FAt => roughly([2.3528873099014249   2247.824863943168566;
                      12.2101689641311815  29.4110913737678104;
                      0.0135489739973833   0.0201862973212001],atol=1e-12)
    end
end
