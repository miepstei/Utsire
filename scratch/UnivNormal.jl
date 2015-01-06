using Distributions

type UnivNormal
    mu::Float64
    sigma::Float64
    getMu::Function
    getSigma::Function
    pdf::Function
    logpdf::Function

    function UnivNormal(mu,sigma)
        this = new()
        this.mu = mu
        this.sigma = sigma

        this.getMu = function()
            return this.mu
        end

        this.getSigma = function()
            return this.sigma
        end

        this.pdf = function(x)
            norm = Normal(this.mu,this.sigma)
            return pdf(norm,x)
        end

        this.logpdf = function(x)
            norm = Normal(this.mu,this.sigma)
            return logpdf(norm,x)
        end
        return this
    end
end
