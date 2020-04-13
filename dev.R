 
summary(vbeta)
summary(veps)

k=2
pl=powerLawCDF(vbeta[k],kmin,kmax)
pcdf=cumsum(vpoints[[k]])
diff=pl-pcdf
plot(diff)

p=vpoints[[1]]
v_mineps=c(1:1000)
v_dst=c(1:1000)
v_beta=c(1:1000)

set.seed(01082019)

for (i in c(1:1000)){
  counting=rmultinom(n=1,size=n,prob=p)
  res=asymptotic_test(alpha,counting,kmin,kmax,scale=1,tol=0.001)
  v_mineps[i]=res[1]
  v_dst[i]=res[2]
  v_beta[i]=res[3]
}

summary(v_mineps)
summary(v_dst)
summary(v_beta)

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
eps=0.08

frequency=list2freq(citySize,kmin,kmax,scale)
set.seed(30062020)
res=bootstrap_test2(alpha, frequency, kmin, kmax,
                    scale,nSimulation=1000, tol=0.001, nDirections=200)
