using Base.Test

#setup code

nopen = 2
k = 5
Q = [-3050        50  3000      0    0;
          2./3. -1502./3.     0    500    0;
          15       0 -2065     50 2000;
          0     15000  4000 -19000    0;
          0         0    10      0  -10;]

burst_intervals = [0.1 0.2 0.1 0.2 0.15 0.16 0.18 0.05 0.1]
burst_lengths = [3 1 5]
tres = 1e-4
tcrit = 1e-2
useChs = 1
likelihood = Array(Cdouble,2)

#calling library from Julia

println("Testing likelihood example from dcprogs")

x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Float64},Ptr{Float64},Csize_t,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres,tcrit,useChs )

#equivilent from dcprogs in python --113.13965841419255

test_handler(r::Test.Success) = println("Success on $(r.expr) dcprogs likelihood calculation")
test_handler(r::Test.Failure) = error("test failed: $(r.expr)")
test_handler(r::Test.Error)   = rethrow(r)

Test.with_handler(test_handler) do
  @test (abs(likelihood[1] - -113.13965841419255) < 1e-12)
  @test (x==0)
end

println("Testing tres change")

x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Csize_t,Ptr{Csize_t},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, 1e-5,tcrit,useChs )

test_handler(r::Test.Success) = println("Success on $(r.expr) dcprogs tres change")

Test.with_handler(test_handler) do
  @test (abs(likelihood[1] - -336.84228995233161) < 1e-12)
  @test (x==0)
end


println("Testing tcrit change")

x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Csize_t,Ptr{Csize_t},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres,1e-3,useChs )
test_handler(r::Test.Success) = println("Success on $(r.expr) dcprogs tcrit change")

Test.with_handler(test_handler) do
  @test (abs(likelihood[1] - -113.12100468596346) < 1e-12)
  @test (x==0)
end

println("Testing using Chs vector")

x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Csize_t,Ptr{Csize_t},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres, tcrit, 0 )
test_handler(r::Test.Success) = println("Success on $(r.expr) dcprogs chs vector change")


Test.with_handler(test_handler) do
  @test (abs(likelihood[1] - -111.51192635114862) < 1e-12)
  @test (x==0)
end

burst_intervals = [0.1 0.2 0.1 0.2 0.15 0.16 0.18 0.05 0.1 0.10 0.18 0.02 0.1 0.03]
burst_lengths = [3 1 5 5]


println("Testing likelihood example from dcprogs")

x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Csize_t,Ptr{Csize_t},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres,tcrit,useChs )


Test.with_handler(test_handler) do
  @test (abs(likelihood[1] - -127.25641384309243) < 1e-12)
  @test (x==0)
end

