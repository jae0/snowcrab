using DataArrays, DataFrames, DataFramesMeta, RCall

R"p = bio.snowcrab::load.environment( year.assessment=2016)
res = biomass.summary.survey.nosa.db(p=p)
region = 'cfasouth' 
biomassindex = res$B[[region]]
catch = res$L[[region]]
yrs = as.numeric(rownames( res$B))
"

R"sb = list(
  O = biomassindex, # observed index of abundance
  Omissing0 = mean(biomassindex, na.rm=TRUE ),
  removals = catch , # removalsches  , assume 20% handling mortality and illegal landings
  removalsmissing0 = mean(catch, na.rm=TRUE ),
  er = 0.2,  # target exploitation rate
  N = length( biomassindex ) , # no years with data
  M = 5, # no years for projections
  MN = length( biomassindex ) + 5,
  ty = which(yrs==2004) ,  # index of the transition year (2004) between spring and fall surveys
  r0= 1,
  K0= 100,
  q0= mean(biomassindex, na.rm=TRUE)/100 ,
  S0= 0.6, # normalised
  cv = 0.5,
  smax =1.25,
  eps = 1e-6
)"

@rget sb

using Distributions, Optim, Mamba
## @everywhere eval(extensions)


samples = 10000
burnin = 250
thin = 1
nchains = 1

sn = Dict{Symbol,Any}()
for i in 1:sb[:N]
  Scur = symbol(string("S", i))
  sn[:($Scur)] = rand( Normal(0.5, 0.1), 1 )
end

inits = [
  Dict{Symbol, Any}(
    :q => sb[:q0],
    :r => sb[:r0],
    :K => sb[:K0],
    :S0 => sb[:S0],
    :S_sd => sb[:cv],
    :O_sd => sb[:cv],
    :S_nodes => sn
  )
  for i in 1:nchains
]

# AR process model needs to be specified outside the model 
S_nodes = Dict{Symbol,Any}()
S_nodes[:S1] = Stochastic( (S0, S_sd) -> LogNormal( log(S0), S_sd ) )
for i in 2:sb[:N]
  global Scur = symbol(string("S", i)),
         Spre = symbol(string("S", i-1))
  @eval begin
    S_nodes[:($Scur)] = Stochastic( 
      ($Spre, r, R, S_sd) -> UnivariateDistribution[
        Sp = $Spre * ( 1.0 + r*(1.0 - $Spre) ) - R[i-1]
        Sp > 0 ? LogNormal( log(Sp), S_sd ) : Flat() 
      ]
    )
  end
end


bd = Model(;
  S_nodes...,
  S = Logical( 1, (S_nodes) -> (
    for i in 1:sb[:N]
      Scur = symbol(string("S", i))
      S[i] = S_nodes[:($Scur)] 
    end )
  ),
  q = Stochastic( (q0, cv) -> LogNormal( log(q0), cv )), 
  r = Stochastic( (r0, cv) -> LogNormal( log(r0), cv )), 
  K = Stochastic( (K0, cv) -> LogNormal( log(K0), cv )), 
  S_sd = Stochastic( (eps) -> Truncated( Cauchy(0, 25), eps, 1 )),
  O_sd = Stochastic( (eps) -> Truncated( Cauchy(0, 25), eps, 1 )),
  R = Logical( 1, (removals,K) -> removals/K, false ),  
  O = Stochastic( 1, (K, q, S, R, O_sd, N) -> UnivariateDistribution[
    begin 
      if i==1            
        Op = K*q*(S[i] - R[i])
      elseif i > 1 && i < ty  
        Op = K*q*(S[i] - R[i-1]) 
      elseif i == ty           
        Op = K*q*(S[i] - (R[i-1] + R[i])/2 ) 
      else         
        Op = K*q*(S[i] - R[i]) 
      end
      Op > 0 ? LogNormal( log(Op), O_sd ) : Flat()
    end
    for i in 1:N] , true ),
  AR = Logical( 1, (r,S) -> ( 1.0 + r*(1.0 - S) ) ), 
  ER = Logical( 1, (R,S) -> R/S ),
  F  = Logical( 1, (ER) -> -log(1-ER) ),
  B  = Logical( 1, (S,K) -> S*K ) ,
  C  = Logical( 1, (R,K) -> R*K ) 
)


draw(bd)

showall(bd)
logpdf(bd)

setsamplers!(bd, [NUTS([:q, :r, :K,  :S_sd, :O_sd])])

sim = mcmc(bd, sb, inits, samples, burnin=burnin, thin=thin, chains=nchains)


gelmandiag(sim1, mpsrf=true, transform=true) |> showall
gewekediag(sim1) |> showall
heideldiag(sim1) |> showall

rafterydiag(sim1) |> showall

describe(sim)
hpd(sim1)
cor(sim1)

## Write to and Read from an External File
write("sim1.jls", sim1)
sim1 = read("sim1.jls", ModelChains)

## Restart the Sampler
sim = mcmc(sim1, 5000)
describe(sim)

## Default summary plot (trace and density plot
p = plot(sim1)

## Write plot to file
draw(p, filename="summaryplot.svg")

## Autocorrelation and running mean plots
p = plot(sim1, [:autocor, :mean], legend=true)
draw(p, nrow=3, ncol=2, filename="autocormeanplot.svg")

## Pairwise contour plots
p = plot(sim1, :contour)
draw(p, nrow=2, ncol=2, filename="contourplot.svg")



















