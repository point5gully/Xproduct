library(primes)
r=c(294,216,135,98,112,84,245,40)
#c(135,45,64,280,70)
#c(210,144,54,135,4,49)
nr=length(r)
p=prime_factors(r)
b=c(8890560,156800,55566)
#c(3000,3969,640)
#c(6615,15552,420)
nc=length(b)
q=prime_factors(b)

targ<-function(a)sum(abs(apply(a,1,prod)-r))+sum(abs(apply(a,2,prod)-b))

#new objects
pr=c(2,3,5,7)
prm=c(3,2,1,1)
mr=rep(0,nc)
mc=rep(0,nr)
rp=matrix(0,4,nr)
cp=matrix(0,4,nc)
for(i in 1:4){
  for(j in 1:nr)rp[i,j]=sum(p[[j]]==pr[i])
  for(j in 1:nc)cp[i,j]=sum(q[[j]]==pr[i])}

sol=matrix(1,nr,nc)

#zover=function(sol,rp,cp){
  
for(tt in 1:100){

for(f in 1:10){
for(i in 1:nr) #rule 2 saturation 
  for(j in (1:4)[rp[,i]>0]){
      mr=(cp[j,]>0)&(sol[i,]*pr[j]<10)
      if(rp[j,i]>prm[j]*(sum(mr)-1)){
        sol[i,mr==1]=sol[i,mr==1]*pr[j]
        cp[j,mr==1]=cp[j,mr==1]-1
        rp[j,i]=rp[j,i]-sum(mr)}}
for(i in 1:nc)
  for(j in 1:4)
    if(cp[j,i]>0){
      mc=(rp[j,]>0)&(sol[,i]*pr[j]<10)
      if(cp[j,i]>prm[j]*(sum(mc)-1)){
        sol[mc==1,i]=sol[mc==1,i]*pr[j]
        cp[j,i]=cp[j,i]-sum(mc)
        rp[j,mc==1]=rp[j,mc==1]-1}}
}

#rule 3 check if only one allocation is feasible
for(i in (1:nr)[apply(rp,2,sum)>0]){#non-empty row
  alo=as.vector(rep(1:4,rp[,i]))
  mr=(cp[alo[1],]>0)&(sol[i,]*pr[alo[1]]<10)
  clo=NULL
  for(t in 1:sum(mr))clo=rbind(clo,sol[i,]*pr[alo[1]]^((1:nc)==(1:nc)[mr][t]))
  cla=clo
  if(length(alo)>1){
  for(j in 2:length(alo)){
    mr=(cp[alo[j],]>0)&(sol[i,]*pr[alo[j]]<10)
    cla=NULL
    for(t in 1:sum(mr))cla=rbind(cla,t(t(clo)*pr[alo[j]]^((1:nc)==(1:nc)[mr][t])))
    clo=cla}}
  cla=cla[!duplicated(cla),]
  if(sum(!!cla)>nc){
   off=1:(dim(cla)[1])
   for(t in 1:dim(cla)[1]) off[t]=max(cla[t,])<10
   cla=cla[!!off,]
   if(sum(!!off)>1){
    off=1:dim(cla)[1]
    for(t in 1:dim(cla)[1])off[t]=max(apply(rbind(cla[t,],sol[-i,]),2,prod)>b)
    cla=cla[!off,]
   }}
  if(sum(!!cla)==nc){#single choiz for row
    rp[,i]=rep(0,4)
    s=prime_factors(cla/sol[i,])
    sol[i,]=cla
    for(j in 1:nc)
      for(k in s[[j]]){
      l=(1:4)[pr==k]
      cp[l,j]=cp[l,j]-1}
  }else{#single choiz for entry
      s=t(t(cla)/sol[i,])
      for(j in 1:nc)
        if(length(z<-unique(s[,j]))==1){
          w=prime_factors(z)
          for(k in w[[1]]){ 
            l=(1:4)[pr==k]
            cp[l,j]=cp[l,j]-1
            rp[l,i]=rp[l,i]-1}
          sol[i,j]=sol[i,j]*z}
  }}}
 
  return(sol)
#}
