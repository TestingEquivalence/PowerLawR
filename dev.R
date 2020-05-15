source("PowerLaw.R")
source("tests.R")
source("power.R")
source("read_data.R")
source("simulation.R")
source("size.R")
source("bootstrap_tests.R") 

kmin=20e3
kmax=10e6
scale=10e3
alpha=0.05
nSamples=1000
n=662
beta=2.3
eps=0.08

vpoint=list()
vbeta=rep(0,100)
vdst=rep(0,100)

set.seed(01012020)
for (i in c(1:100)){
  #vpoint[[i]]=uniformRandomStress(kmin=kmin/scale,kmax=kmax/scale,beta,eps)
  vpoint[[i]]=powerLawStress(n,eps,kmin=kmin/scale,kmax=kmax/scale,beta)
  res=nearestPowerLaw(cumsum(vpoint[[i]]),kmin/scale,kmax/scale,1,3)
  vbeta[i]=res$minimum
  vdst[i]=res$objective
}

sums=sapply(vpoint, sum)
View(sums)
mins=sapply(vpoint, min)
mins>=0
maxs=sapply(vpoint, max)
maxs<=1


k=50
vbeta[k]
vdst[k]
pl=powerLawDensity(vbeta[k],kmin/scale,kmax/scale)
df=vpoint[[k]]
diff=cumsum(pl)-cumsum(df)
plot(diff)
res=nearestPowerLaw(cumsum(vpoint[[k]]),kmin/scale,kmax/scale,2,3,tol=0.001)


k=2
p=vpoint[[k]]
v_mineps=c(1:1000)
v_dst=c(1:1000)
v_beta=c(1:1000)

set.seed(01082019)

res=nearestPowerLaw(cumsum(p),kmin/scale,kmax/scale,2,3,tol=0.001)
v=rep(0,length(p))

for (i in c(1:1000)){
  counting=rmultinom(n=1,size=n,prob=p)
  v=v+counting/n
  #res=nearestPowerLaw(cumsum(counting/n),kmin/scale,kmax/scale,2,3,tol=0.001)
  #v_beta[i]=res$minimum
  #v_dst[i]=res$objective
  res=asymptotic_test(alpha,counting,kmin/scale,kmax/scale,scale=1,tol=0.001)
  v_mineps[i]=res[1]
  v_dst[i]=res[2]
  v_beta[i]=res[3]
}

summary(v_mineps)
summary(v_dst)
summary(v_beta)

diff=p-v/1000
plot(diff)
l2(cumsum(p),cumsum(v/1000))
res=nearestPowerLaw(cumsum(v/1000),kmin/scale,kmax/scale,2,3,tol=0.001)


hist(v_mineps)
hist(v_dst)
hist(v_beta)

v=v_mineps<eps
sum(v==TRUE)/1000

citySize=readVector("C:\\Users\\Ostrovski\\Google Drive\\Writing\\PowerLaw\\CitySize\\list_ge.csv")
alpha=0.05
kmin=20e3
kmax=10e6
scale=10e3
nSamples=1000
n=662
eps=0.065
tol=0.001
nSimulation=1000
nDirections=100
minEps=0.0607
maxEps=0.0608

frequency=list2freq(citySize,kmin,kmax,scale)

asymptotic_test(alpha,frequency,kmin,kmax,scale,tol)

set.seed(30062020)
bootstrap_test2(frequency,kmin,kmax,scale,nSimulation,nDirections,eps=.06075,tol)

set.seed(30062020)
res=bootstrap_test2_1(alpha,frequency,kmin,kmax,scale,nSimulation,
                tol,nDirections)  
