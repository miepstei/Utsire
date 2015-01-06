

N=100000 #number of samples
NDATA=30
data = rand(Normal(2,15),NDATA)

sigma = 5 #for sampler

samples = Array(Float64,N,2)
samples[1,1] = -10
samples[1,2] = 50

lognormals = zeros(Float64,N)
proposedlognormals = zeros(Float64,N)
proposals = zeros(Float64,N,2)

model = UnivNormal(samples[1,1],samples[1,2])

for j=1:NDATA
    lognormals[1]+=model.logpdf(data[j])
end

proposedlognormals[1] = lognormals[1]

println("Sampling $N samples")
for i=2:N
    Q = MvNormal(vec(samples[i-1,:]),[sigma 0.0; 0.0 sigma])
    proposal = rand(Q)
    proposals[i,:] = proposal
    model = UnivNormal(proposal[1],proposal[2])
    
    newlognormal = 0
    for j=1:NDATA
        newlognormal+=model.logpdf(data[j])
    end 
    proposedlognormals[i] = newlognormal

    #accept/reject, (symmetric proposal density)
    alpha = min(0 , newlognormal - lognormals[i-1])

    if(alpha > log(rand()))
        samples[i,:] = proposal
        lognormals[i] = newlognormal
    else
        samples[i,:] = samples[i-1,:]
        lognormals[i] = lognormals[i-1]
    end

end
println("Sampling complete")


