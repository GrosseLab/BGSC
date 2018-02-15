#HERLEITUNG 

### Herleitungs TEST

######### TEST ob log(P) richtig --> JA STIMMT !
get_sum<-function(d){
  s=0
  for(i in 1:length(d$v)){
    s=s+( d$v[i]-d$m )^2  
    a=d$v[i]-d$m
  }
  return(s)
}

da_a0[[3]]=get_sum(da_a0)
da_b0[[3]]=get_sum(da_b0)
da_b1[[3]]=get_sum(da_b1)
da_c0[[3]]=get_sum(da_c0)
da_c1[[3]]=get_sum(da_c1)
da_d0[[3]]=get_sum(da_d0)
da_d1[[3]]=get_sum(da_d1)

lnL_a <- -N * log(sqrt(2*pi)) + N * log(sqrt(tau)) - 0.5*tau*( da_a0[[3]] )
lnL_b <- -N * log(sqrt(2*pi)) + N * log(sqrt(tau)) - 0.5*tau*( da_b0[[3]] ) - 0.5*tau*( da_b1[[3]] )
lnL_c <- -N * log(sqrt(2*pi)) + N * log(sqrt(tau)) - 0.5*tau*( da_c0[[3]] ) - 0.5*tau*( da_c1[[3]] )
lnL_d <- -N * log(sqrt(2*pi)) + N * log(sqrt(tau)) - 0.5*tau*( da_d0[[3]] ) - 0.5*tau*( da_d1[[3]] )
########
############


Ns<-c()
a0<-c(1,2,3,4,5,6);   Ns[1]  <- length(a0);
b0<-c(1,3,5);         Ns[2] <- length(b0);
b1<-c(2,4,6);         Ns[3] <- length(b1);
c0<-c(1,3,5,4);       Ns[4] <- length(c0);
c1<-c(2,6);           Ns[5] <- length(c1);
d0<-c(1,3,5,4,6);     Ns[6] <- length(d0);
d1<-c(2);             Ns[7] <- length(d1);
sub<-list(a0=a0,b0=b0,b1=b1,c0=c0,c1=c1,d0=d0,d1=d1)

da=d[1,]
mu_da<-calc.MUs(da,sub)
da_a0=list(v=da,m=mu_da[1])
da_b0=list(v=da[sub[[2]]],m=mu_da[2])
da_b1=list(v=da[sub[[3]]],m=mu_da[3])
da_c0=list(v=da[sub[[4]]],m=mu_da[4])
da_c1=list(v=da[sub[[5]]],m=mu_da[5])
da_d0=list(v=da[sub[[6]]],m=mu_da[6])
da_d1=list(v=da[sub[[7]]],m=mu_da[7])


N=6 ; Nx=4 ; Ny=2 ;Na=3 ; Nb=3
P.C.0<-((Nx * Ny) / N) * ( (sum(da_c0$v)/Nx) - ((sum(da_c1$v)/Ny) ) )^2
P.B.0<-((Na * Nb) / N) * ( (sum(da_b0$v)/Na) - ((sum(da_b1$v)/Nb) ) )^2
Erg0<-P.C.0 -P.B.0

N=6 ; Nx=4 ; Ny=2 ;Nd=5 ; Ne=1
P.C.0<-((Nx * Ny) / N) * ( (sum(da_c0$v)/Nx) - ((sum(da_c1$v)/Ny) ) )^2
P.D.0<-((Nd * Ne) / N) * ( (sum(da_d0$v)/Nd) - ((sum(da_d1$v)/Ne) ) )^2
Erg0d<-P.C.0 - P.D.0

P.D.0 - P.C.0  

z2=c(da_d0$v,da_d1$v)
## restrucktieren !!! in x und y ( x = - ; y = + )
z=c(da_c0$v,da_c1$v)
# Startpunkt
P.C<-((Nx * Ny) / N) * ( (sum(z[1:Nx])/Nx) - ((sum(z[(Nx+1):N])/Ny) ) )^2
P.B<-((Na * Nb) / N) * ( (sum(z[1:Na])/Na) - ((sum(z[(Na+1):N])/Nb) ) )^2
Erg1<-P.C -P.B  

# Startpunkt nach A, B , C umwandlung
A=sum(z[1:Na])
C=sum(z[(Nx+1):N])
B=sum(z[(Na+1):Nx])       

P.c<-((Nx * Ny) / N) * ( (A+B)/Nx -    C/Ny )^2 
P.b<-((Na * Nb) / N) * (     A/Na - (C+B)/Nb  )^2
Erg2<-P.c -P.b  

## meine Herleitung
v1<-(Na*Ny-Nb*Nx)*(Nb*Ny)
v2<-(Nb*Ny-Na*Nx)*(Na*Ny)
v3<-(Nb*Nx-Na*Ny)*(Na*Nx)
v4<-   2 *( (Ny+Nx)*(Na*Nb*Ny) )
v5<- (-2)*( (Na+Nb)*(Na*Nx*Ny))
Erg3<-( (v1*A^2) + (v2*B^2) + (v3*C^2) + (v4*A*B) + (v5*B*C) ) / (N*Na*Nb*Nx*Ny)    


get.v<-function(Ns){
  Nx=Ns[1] ; Ny=Ns[2] ;Na=Ns[3] ; Nb=Ns[4]
  v1<-(Na*Ny-Nb*Nx)*(Nb*Ny)
  v2<-(Nb*Ny-Na*Nx)*(Na*Ny)
  v3<-(Nb*Nx-Na*Ny)*(Na*Nx)
  v4<-   2 *( (Ny+Nx)*(Na*Nb*Ny) )
  v5<- (-2)*( (Na+Nb)*(Na*Nx*Ny))
  return(c(v1,v2,v3,v4,v5))
}

get.v2<-function(Ns){ # short form of get.v
  x=Ns[1] ; y=Ns[2] ;a=Ns[3] ; b=Ns[4]
  v1<-(a*y-b*x)*(b*y)
  v2<-(b*y-a*x)*(a*y)
  v3<-(b*x-a*y)*(a*x)
  v4<-   2 *( (y+x)*(a*b*y) )
  v5<- (-2)*( (a+b)*(a*x*y))
  return(c(v1,v2,v3,v4,v5))
}

get.v3<-function(Ns){ # after  substitution x=(a+q) ; b=(y+q)
  q=Ns[1]-Ns[3]
  y=Ns[2]
  a=Ns[3]
  x=(a+q)
  b=(y+q)
  
  v1<-   -1 * y     * q * (  a*(y + q) + (y+q)^2 )
  v2<-        y * a *     (  q*(3a- y) + (y+q)^2 - (a+q)^2)     #  y * a *     (  y*(y+q) -a*(a+q) )       # y^3*a + y^2*a*q   - y*a^3 - y*a^2*q    #    
  v3<-            a * q * (  y*(a + q) + (a+q)^2 )  #2 * ( y*a^2*q + a^3*q +a*q^3 +2*a^2*q^2 + y*a*q^2 )
  v4<-    2 * y * a *     (  a*(y + q) + (y+q)^2 )   
  v5<- (-2) * y * a *     (  y*(a + q) + (a+q)^2 )
  return(c(v1,v2,v3,v4,v5))
}
v<-get.v(c(4,2,3,3))
v<-get.v3(c(4,2,3,3))

Erg3<-( (v[1]*A^2) + (v[2]*B^2) + (v[3]*C^2) + (v[4]*A*B) + (v[5]*B*C) ) / (N*Na*Nb*Nx*Ny)    


get.v(c(4,2,5,1))
get.v3(c(4,2,5,1))

v_c_b<-
  v_c_d<-get.v(c(4,2,5,1))

v_b_c<-get.v(c(3,3,4,2))
v_b_d<-get.v(c(3,3,5,1))

v_d_c<-get.v(c(5,1,4,2))
v_d_b<-get.v(c(5,1,3,3))

v_c_b/v_c_b[1]
v_c_d/v_c_d[1]
v_b_c/v_b_c[1]
v_b_d/v_b_d[1]
v_d_c/v_d_c[1]
v_d_b/v_d_b[1]


### BIC -BIC kontrolle ob rechnug stimmt ---> PASST !!!

La<-get.ML.logLs(t(as.matrix(da)),tau)
BICa<-make.IC(La,c(2,3),6)
#c-b
BICa$IC[3]-BICa$IC[2] 
-2*La[3]- (-2*La[2])
Erg0
Erg1
Erg2
Erg3

# c-a 
BICa$IC[3]-BICa$IC[1] 
log(N)+( -2*La[3]- (-2*La[1]))
log(N) - tau * P.C

# c-d 
BICa$IC[3]-BICa$IC[4] 
-2*La[3]- (-2*La[4])
v<-get.v(c(4,2,5,1))
((v[1]*A^2) + (v[2]*B^2) + (v[3]*C^2) + (v[4]*A*B) + (v[5]*B*C) ) / (N*Na*Nb*Nx*Ny)

get.v3(c(4,2,5,1))

### Herleitungs TEST

#  IVO : X, Y Z 
z
Nx<-3; Nz<-1 ; Ny<-2 ;N<-Nx+Nz+Ny

X<-sum(z[1:Nx])
Z<-sum(z[(Nx+1):(Nx+Nz)])
Y<-sum(z[(Nx+Nz+1):N])


P.c<-((Nx+Nz)*Ny)/N * ( (X+Z)/(Nx+Nz) - Y/Ny )^2 
P.b<-((Ny+Nz)*Nx)/N * (  X/Nx         - (Y+Z)/(Ny+Nz) )^2
Erg2<-P.c -P.b


qq<-function(X,Y,Z,x,y,z){
  q<- ( -(-x*(X+Y)+(x+y)*Z)^2/(x*(y+z))+(-Y*(x+z)+y*(X+Z))^2/(y*(x+z)))
  return(q)
}
(1/N) *qq(X,Y,Z,Nx,Ny,Nz)

#  CLAUS : X, Y Z
Na<-3;  Nb<-1 ; Nc<-2  ;N<-Na+Nb+Nc

X<-sum(z[1:Na])
Y<-sum(z[(Na+1):(Na+Nb)])
Z<-sum(z[(Na+Nb+1):(Na+Nb+Nc)])


wolfr<-function(X,Y,Z,a,b,c){
  #((a+q)*y)/((a+q)+y)*((A+B)/(a+q)-C/y)^2 - (a*(y+q))/(a+(y+q))*( A/a - (C+B)/(y+q) )^2
  #((a+b)*c)/((a+b)+c)*((X+Y)/(a+b)-Z/c)^2 - (a*(b+c))/(a+(b+c))*( X/a - (Z+Y)/(b+c) )^2
  
  #qw<-(2*a*A*B+a*B^2-q*A^2)/(a*(a+q)) - (B+C)^2/(q+y) + C^2/y
  qw<-(2*a*X*Y+a*Y^2-b*X^2)/(a*(a+b)) - (Y+Z)^2/(b+c) + Z^2/c
  
  return(qw)
  
}
#wolfr(A,B,C,Na,1,Ny)
wolfr(X,Y,Z,Na,Nb,Nc)


ivo<-function(X,Y,Z,a,b,c){
  ######## !!!!!!!!!!!! ENDRESULTAT !!!!!!!!!!!!
  v1<- -a*b/(a+b)
  v2<- -b*c/(b+c)
  
  return( v1*( X/a - Y/b)^2 - v2* ( Z/c - Y/b)^2 )
  
}

ivo2<-function(X,Y,Z,a,b,c){
  
  N=(a+b+c)
  
  return( (a+b)*c/N * ( (X+Y)/(a+b) - Z/c )^2 - (b+c)*a /N * (X/a - (Y+Z)/(b+c))^2)
  
}

ivo(X,Y,Z,Na,Nb,Nc)
ivo2(X,Y,Z,Na,Nb,Nc)



#############################
N=6
Ns<-c()
a0<-c(1,2,3,4,5,6);   Ns[1]  <- length(a0);
b0<-c(1,3,5);         Ns[2] <- length(b0);
b1<-c(2,4,6);         Ns[3] <- length(b1);
c0<-c(1,3,5,4);       Ns[4] <- length(c0);
c1<-c(2,6);           Ns[5] <- length(c1);
d0<-c(1,3,5,4,6);     Ns[6] <- length(d0);
d1<-c(2);             Ns[7] <- length(d1);
sub<-list(a0=a0,b0=b0,b1=b1,c0=c0,c1=c1,d0=d0,d1=d1)

La<-get.ML.logLs(SF[5:15,],tau)
BICa<-make.IC(La,c(2,3),6)

#c-b
BICa$IC[,3]-BICa$IC[,2] 
-2*La[,3]- (-2*La[,2])
# c-a 
BICa$IC[,3]-BICa$IC[,1] 
log(N)+( -2*La[,3]- (-2*La[,1]))

BICa$IC[,1]-BICa$IC[,3] <0

# c-d 
BICa$IC[,3]-BICa$IC[,4] 
-2*La[,3]- (-2*La[,4])


#  CLAUS : X, Y Z

reSF<-cbind(SF[5:15,sub[[6]]],SF[5:15,sub[[7]]])
colnames(reSF)<-c( "SF767.1","SF767.3","SF767.5","SF767.4","SF767.6","SF767.2")

#b-c
Na<-3;  Nb<-1 ; Nc<-2  ;N<-Na+Nb+Nc
X<-rowSums(reSF[,1:Na])
#Y<-rowSums(reSF[,(Na+1):(Na+Nb)])
Y<-reSF[,(Na+1):(Na+Nb)] # bei Nb<-1
Z<-rowSums(reSF[,(Na+Nb+1):(Na+Nb+Nc)])

ivo<-function(X,Y,Z,a,b,c){
  ######## !!!!!!!!!!!! ENDRESULTAT !!!!!!!!!!!!
  v1<- -a*b/(a+b)
  v2<- -b*c/(b+c)
  
  return( v1*( X/a - Y/b)^2 - v2* ( Z/c - Y/b)^2 )
  
}

#c-d
Na<-4;  Nb<-1 ; Nc<-1  ;N<-Na+Nb+Nc
X<-rowSums(reSF[,1:Na])
#Y<-rowSums(reSF[,(Na+1):(Na+Nb)])
Y<-reSF[,(Na+1):(Na+Nb)] # bei Nb<-1
Z<-reSF[,(Na+Nb+1):(Na+Nb+Nc)]

ivo(X,Y,Z,Na,Nb,Nc)

#TODO
#c-a
Na<-4;  Nb<-2 ; Nc<-0  ;N<-Na+Nb+Nc ## PASST NOCH NICHT !!!! #TODO 
X<-rowMeans(reSF[,1:Na])
Y<-rowMeans(reSF[,(Na+1):(Na+Nb)])

 log(N) - tau*((4*2)/6 *(X*Y)^2) < 0
 (X*Y) > sqrt(log(N) * 1/tau *(6/(4*2)) )
