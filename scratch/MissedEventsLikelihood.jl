using MAT


datafile = "/home/michaelepstein/Utsire/data/AchRealData.mat"
file = matopen(datafile)

bursts = read(file,"bursts")
concs = read(file,"bursts")
tcrit = read(file,"tcrit")
tres = read(file,"tres")
useChs =  read(file,"useChs")
close(file)

likelihood = Array(Float64,3)

nopen = 2

k = 5

Q = [-3025        25  3000      0    0;
          2./3. -1502./3.     0    500    0;
          15       0 -2065     50 2000;
          0     15000  4000 -19000    0;
          0         0    10      0  -10;]

experimentno = length(concs)
println("number of sets = $experimentno")
# number of bursts in each set
burstnos = zeros(Int32,experimentno)
likelihoods = zeros(Float64,experimentno)

#need the indexes for each burst set to delineate bursts
counter=1
for i=1:experimentno
  burst_lengths = zeros(Csize_t,length(bursts[i]))
  for j=1:length(bursts[i])
    burst_lengths[j] = length(bursts[i][j])
  end

  bsts = zeros(sum(burst_lengths))
  counter = 1
  for j = 1:length(bursts[i])
    for m = 1:length(bursts[i][j])
      bsts[counter] = bursts[i][j][m]
      counter += 1
    end  
  end

  println(string("number of burst intervals = " ,length(bsts)))
  println(string("number of bursts = ", sum(burst_lengths)))
  #println(string("burst intervals = " ,bsts))
  
  likelihood = Cdouble[0]
  "extern C int missed_events_likelihood(double* jl_likelihood , double* jl_bursts, size_t interval_length, size_t* burst_lengths, size_t burst_number, double* Q , int nopen, int k, double tau, double tcrit, bool useChs)"
  x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Csize_t,Ptr{Csize_t},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,bsts,length(bsts),burst_lengths,length(burst_lengths),Q, nopen, k, tres[i],tcrit[i],useChs[i])
  likelihoods[i] = likelihood[1]
  println(string("Likelihood for set ", i , " = ", likelihood[1]))
end

println(string("Likelihood for test set ", sum(likelihoods)))
println("")
