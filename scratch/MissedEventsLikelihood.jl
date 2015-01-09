using MAT


datafile = "/home/michaelepstein/Utsire/data/AchRealData.mat"
file = matopen(datafile)

bursts = read(file,"bursts")
concs = read(file,"bursts")
tcrit = read(file,"tcrit")
tres = read(file,"tres")


close(file)

likelihood = Array(Float64,1)
nopen = 2
k = 5
Q = [-3025        25  3000      0    0;
          2./3. -1502./3.     0    500    0;
          15       0 -2065     50 2000;
          0     15000  4000 -19000    0;
          0         0    10      0  -10;]
burst_number = length(bursts[1]) 
println(string("Bursts = ", burst_number))
burst_lengths = Array(Int32 , burst_number)

for i=1:burst_number
  burst_lengths[i] = length(bursts[1][i])
end

bsts = zeros(Float64,sum(burst_lengths))
counter=1;
for i=1:burst_number
  for j=1:burst_lengths[i] 
    bsts[counter] = bursts[1][i][j]
    counter = counter + 1
  end
end


println(string("number of burst intervals = " ,length(bsts)))
println(string("burst intervals = " ,bsts))
println(string("Length of first burst = ", burst_lengths[1]))
println(string("Length of second burst = ", burst_lengths[2]))
println(string("Length of third burst = ", burst_lengths[3]))
println(string("typeof burst[1] = "),typeof(bursts[1]))


"double* jl_likelihood , double** jl_bursts, size_t burst_number, int* burst_lengths, int sizeAny, double* Q , int nopen, int k, double tau, double tcrit"
#x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Float64},Ptr{Float64},Int32,Ptr{Int32},Int32,Ptr{Float64},Int32,Int32,Float64,Float64),likelihood,bsts,burst_number,burst_lengths,128,Q, nopen, k, tres[1],tcrit[1])
