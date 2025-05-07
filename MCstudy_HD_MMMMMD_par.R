# Computes vertices of C* for simulated data using true function eta
# or non parametrically estimated eta^:
# for sample sizes large enough for eta^ to be monotonic in each of its arguments,
# the inequalities needed are still only those at range(x), range(v0).

# NOTE1: in "vertices" this version only retains unique points

# NOTE2: produces plots for Minkowski average

# NOTE3: should try to answer the following questions:
# 1. is C^MD in C^MMM for each sample?
# 2. are the focal points necessarily in C^MD?  see CORRECTION
# 3. are they necessarily in C^MMM? 

# NOTE4: quantile plots: HD < Q25 und HD > Q75

# NOTE5: dissociate TOL and tolrank

# NOTE6: use "trimmed" extremes: for v0, maxv0min minv0max, where these values are obtained in a preliminary 
#        simulation for the given N and resetting the seed;
#        and for x quantiles atrim/2, 1-atrim/2 .
# In fact I am unable to circumvent the compatibility problem between the bandwidth object defined for the first mc
# and the X.dat for the other mcs
 
rm(list = ls())

mac = 1
# if (mac==0) { 
# #setwd("H:/PART/nutzbar")                    # set working directory
# #setwd("J:/Daten/fla/ManskiTamer/nutzbar")
# #setwd("D:/SauveD/ZEW/HoE_5_6")
# #setwd("E:/Dokumente und Einstellungen/hul/Desktop/MT")
# #setwd("J:/Daten/fla/ManskiTamer/nutzbar/TakingStockMarch2011")
# setwd("C:/Dokumente und Einstellungen/fam-edv5/Desktop/MCstudy_tex/plotsuse/MMDnormtrueeta") }
# if (mac==1) { setwd("/Users/hannesu/Documents/Papers/ManskiTamer/work_october/MMD") }

tstart = date()
cpu0   = proc.time()
library(scatterplot3d)
library(geometry)
# library(chplot)
library(graphics)
library(Matrix)
library(np)
#library(help=np)
set.seed(123456789)  
                                            # set parameters of MC experiment 

NR = 3        # Number of regressors
MC = 1000      # Number of MC iterations
N = 100       # Random sample size
inorm = 0     # 0 --> uniform distr, 1 --> normal distr

inp = 0       # nonparametric estimation of eta ? 1 yes, 0 no (set to 0 when doing MMM)
              # the indented lines concern parameters for nonparametric estimation
	hmax = 10     # if bw$bw[2] > hmax, replace it with rule of thumb bandwidth
      nmult = 10
	bw_new = 0    # if 0, use starting values or set bandwidth
	              # if 1, do not use starting values
	N_ = 200      # previous sample size (used in scaling bw)
	bw_imp = c( 0.1320053, 5.0029208)     # import previous bw as starting values
	bw_imp[1] = bw_imp[1]*(N_/N)^(2/5)  # scale second component down
	bw_imp[2] = bw_imp[2]*(N_/N)^(1/5)  # scale second component down
	bw_c = FALSE   # value for bandwidth.compute 
	reg_t  = "ll"        # ll (local linear) or lc (local constant)
	bw_m   = "cv.ls"     # cv.aic or cv.ls
	cker_t = "gaussian"   # gaussian, epanechnikov, or uniform   
	bw_t   = "fixed"     # fixed, adaptive_nn or generalized_nn (only fixed found to function)
      atrim  = .0   # see NOTE 6 above     

method  = c("MMM","MMD")
mmm = 2       # MMM: mmm=1, MMD: mmm=2
K = 30         # number of powers of abs(x) used as instruments for MMM
MaxVert = 10*(K+10) # maximum number of vertices stored
#TOL = min(N/100,6)
TOL = 6
tol=10^(-TOL)
tolrank = 1e-6
p = cbind(0.25,0.75)

minv0max = 100
maxv0min = -100
if (inp == 1 & inorm == 1) {                 # loop to obtain minv0max and maxv0min for inp == inorm == 1	
	for (i in 1:N) { v = rnorm(N,0,sqrt(2))   # NOTE: R uses sd as secnd param in normal
		v1 = ceiling(v)
		v0 = v1 - 1
            minv0max = min(minv0max,max(v0))
            maxv0min = max(maxv0min,min(v0))
                      }
print("common extremes for v0"); print(c(maxv0min,minv0max))
set.seed(123456789)        }

if (inorm == 0) { dist = "Unif"  }              
if (inorm == 1) { dist = "Norm"  } 
method  = c("MMM","MMD")

#--------------------------------------------
# ColM 
 ColM<-function(z,w) {
#        wz = matrix(NA,dim(w*z),1)
#for i in 1:dim(w*z) {
#        wz[i,1] = w*z 
#        x = cbind(colMeans(w*z),0)
        x = colMeans(w*z)
        return(x)
        }
#--------------------------------------------------------------------------------
 eta<-function(x,v0,i,a,b){   # we use estimated parameters instead of the true ones
    if (i==0) {return(gamm0[1]*v0+gamm0[2]*x+gamm0[3]) }
    if (i==1) {
    v1=v0+1
    v0_ = sd_v*((dnorm(v0/sd_v)-dnorm(v1/sd_v))/(pnorm(v1/sd_v)-pnorm(v0/sd_v)))
   return( gamm[1]*v0_+gamm[2]*x+gamm[3]) }
}
#--------------------------------------------------------------------------------
 vert<-function(allvert) { # may not be the end yet, but should prevent crashes
                           # see program checkvert for explanations
 d = dim(allvert) 
if (d[1]>d[2]) {                
    cv = convhulln(allvert)
    if (ncol(allvert) == 3) {
    vertices = allvert[unique(c(cv[,1],cv[,2],cv[,3])),]} # unique extracts the indices of the vertices
    if (ncol(allvert) == 2) {
    vertices = allvert[unique(c(cv[,1],cv[,2])),]} 
    return(vertices) }
if (d[1]<=d[2]) { return(allvert) }
                          }
#--------------------------------------------------------------------------------
 filtna<-function(allvert) {                 
    i = min(which(is.na(allvert[,1]))-1,dim(allvert)[1])
    vertices = allvert[1:i,]
    return(vertices)
}
#--------------------------------------------------------------------------------
# orders vertices couterclockwise so that vertex 1 has lowest y 
# and in case of a tie the lower x 
ccw <- function(vv) {

nv = nrow(vv)
vert1 = vv[which.min(vv[,2]),]      # find smallest y
vert1=rbind(vert1,c(0,0))
vert1u = vert1[which.min(vert1[1:nrow(vert1)-1,1]),]

a =     matrix(0,nv)             # compute angle of each vertex to v1
for (i in 1:nv) { a[i] = ang(rbind(vert1u,vv[i,])) }

V = vv[order(a),]               # rank vertices according to angle

return(V) }
# ------------------------------------------------------------
# This computes the angle of the line between (x1,y1) and (x2,y2)
# and positive of the X axis. Courtesy of A.Beresteanu & F.Molinari
ang <- function(v) {                 
    x1 = v[1,1]; y1 = v[1,2]
    x2 = v[2,1]; y2 = v[2,2]

px = x2-x1
py = y2-y1

flag = 0

if (py < 0) { py=-py; flag=1 } 

xang=atan2(py,abs(px))

if (px < 0) { xang=pi-xang }

if (flag==1) { xang=2*pi-xang }

return(xang)
}
# --------------------------------------------------------------
# computes Minksum between two polygons A, B
# must be organized according to ccw. Courtesy of A.Beresteanu & F.Molinari
mks <- function(A,B) {
dA = dim(A); dB = dim(B)
if (dA[1]>dA[2]) { Ac = convhulln(A); A=unique(rbind(A[Ac[,1],],A[Ac[,2],])) }

A = ccw(A)
B = ccw(B)
m = dim(A)[1]
n = dim(B)[1]
i = 1
j = 1
M = matrix(0,nrow=1,ncol=2)

tol = 1e-10
aA = ang(rbind(A[m,],A[1,]))
aB = ang(rbind(B[n,],B[1,]))

if (aA > aB) {
    while ( i<m+1 | j<n+1 ) {

        ci = min((i%%(m+1))+1,i)
        ni = min(((i+1)%%(m+1))+1,(i+1))

        cj = min((j%%(n+1))+1,j)
        nj = min(((j+1)%%(n+1))+1,(j+1))


        M = rbind(M,(A[ci,]+B[cj,]))
        aA = ang(rbind(A[ci,],A[ni,]))
        aB = ang(rbind(B[cj,],B[nj,]))

        ii = 0
        jj = 0

        if ( aA <= aB+tol | j==n+1) {
            ii = 1
        }
        if ( aA >= aB-tol | i==m+1) {
            jj = 1
        }

        i = i+ii
        j = j+jj
        i = min(i,m+1)
        j = min(j,n+1)

    }
} 

if (aA <= aB)  {
    while ( i<m+1 | j<n+1 ) {

        ci = min((i%%(m+1))+1,i)
        ni = min(((i+1)%%(m+1))+1,(i+1))

        cj = min((j%%(n+1))+1,j)
        nj = min(((j+1)%%(n+1))+1,(j+1))        

        M = rbind(M,(B[cj,]+A[ci,]))
        aA = ang(rbind(A[ci,],A[ni,]))
        aB = ang(rbind(B[cj,],B[nj,]))

        ii = 0
        jj = 0

        if ( aB <= aA+tol | i==m+1) {
            jj = 1
        }
        if ( aB >= aA-tol | j==n+1) {
            ii = 1
        }

        i = i+ii
        j = j+jj
        i = min(i,m+1)
        j = min(j,n+1)

    }
}
    M = M[2:nrow(M),]
    if (dim(M)[1] > 2) { MCV = convhulln(M,)
                         M = M[unique(c(MCV[,1],MCV[,2])),] }
return(M)
} 
# -----------------------------------------------
# Minimmal distance between the dot P3 and the segment
# between the points P1 and P2. It returns the distance 
# in dist and the closest point on the segment to P3.
# Courtesy of A.Beresteanu & F.Molinari.

dotdist <- function(P1,P2,P3) {
P = c(0,0)
if (P1[1] == P2[1] && P1[2] == P2[2]) {
    P[1] = P1[1]
    P[2] = P1[2]        # The segment is a point
} else {
    u = (P3[1]-P1[1])*(P2[1]-P1[1])+(P3[2]-P1[2])*(P2[2]-P1[2])
      u = u/((P1[1]-P2[1])^2+(P1[2]-P2[2])^2)
      if (u>0 && u<1) {             # The closest point is
        P[1]=P1[1]+u*(P2[1]-P1[1])  # inside the segment
            P[2]=P1[2]+u*(P2[2]-P1[2])  
    } else if (u<0) {
        P[1]=P1[1]
            P[2]=P1[2]
    } else {
        P[1]=P2[1]
            P[2]=P2[2]
    }
}       
dist = sqrt((P[1]-P3[1])^2+(P[2]-P3[2])^2)
return(dist)
}
# -----------------------------------------------
# computes Hausdorff distance between two polygons A, B
# organized according to ccw. Courtesy of A.Beresteanu & F.Molinari

HD <- function(A,B) {

dstarAB = dHD(A,B)
dstarBA = dHD(B,A)

HD = max(dstarAB,dstarBA)
return(HD)
}
# ----------------------------------------------------
# computes the directed Hausdorff distance between two polygons A, B
# organized according to ccw. Courtesy of A.Beresteanu & F.Molinari

dHD <- function(A,B) {

if (dim(t(as.matrix(A)))[1] == 1) { na = 1 } else { na = nrow(A) }
if (dim(t(as.matrix(B)))[1] == 1) { nb = 1 } else { nb = nrow(B) }
                          # starting points of A and B
if (na>1)  { A0 = A[1,] }
if (na==1) { A0 = A  }
if (nb>1)  { B0 = B[1,] }
if (nb==1) { B0 = B  }
#print("na");print(na)
#print("nb");print(nb)
                          # close polygons
A = rbind(A,A0)
#print("A");print(A)
B = rbind(B,B0)
#print("B");print(B)

dstarAB = 0                          # dH(A,B)
for (i in 1:na) {                    # compute d(a,B)
   dist = Inf
   for (j in 1:nb) {
#      print(B[j:(j+1),]);print(A[i,])
    newdist = dotdist(B[j,],B[j+1,],A[i,])
      if (dist > newdist) { dist = newdist }
#      print(" ")
#      print("i");print(i)
#      print("newdist");print(newdist);print("dist");print(dist)
   }
#   print("dstarAB"); print(dstarAB)
   if (dstarAB < dist) { dstarAB= dist }
#   print("dstarAB"); print(dstarAB)
}   
dHD = dstarAB
return(dHD)
}
# -----------------------------------------------

#Working example for 2D Minkowski sum
#A=rbind(c(1,1),c(2,1),c(2,2))
#B=rbind(c(3,1),c(4,1),c(4,2),c(3,2))
#C=mks(A,B)

#if (inp == 1) { tol = 1e-0/sqrt(N) }              # this should still be checked (set it low enough that no crash arises)
focval = matrix(NA,MC,NR)
#if (mmm == 2) { array.mvertices <- array(NA, dim = c(MC, MaxVert, NR))  } # what is this ???
array.mvertices <- array(NA, dim = c(MC, MaxVert, NR))
array.mcvert <- array(NA, dim = c(MC, MaxVert, NR))   # Make array containing final unique vertices
array.mvertices2D <- array(NA, dim = c(MC, MaxVert, 2, 3)) # MC, MaxVert vertices, 2D coordinates, 3 projections
array.r_ <- array(0, dim = c(MC, NR, 2))
mvolume = matrix(NA,MC,4)

# Arrays for Hausdorff distances and their quantiles
array.HD <- array(0, dim = c(MC, 2, 3))
array.HDsorted <- array(0, dim = c(MC, 2, 3))
array.dHD0 <- array(0, dim = c(MC, 2, 3))
array.dHD1 <- array(0, dim = c(MC, 2, 3))

HDquant = matrix(NA,3,ncol(p))

mc = 1
mceff = 0
dvertmin = MaxVert
focvalpb = 0

while (mc < (MC+1)) { # mc loop =======================================================================================
mceff = mceff + 1
print(paste("mc ",mc," mceff ",mceff))

if (inorm == 0) { dist = "Unif" ; x = runif(N,0,5); v = runif(N,-2,3) }              
if (inorm == 1) { dist = "Norm" ; x = rnorm(N,1,2); v = rnorm(N,0,sqrt(2)) }  # NOTE: R uses sd as secnd param in normal

sd_v = sd(v)
sd_v

y = 1 + v - x + rnorm(N,0,1)
v1 = ceiling(v)
v0 = v1 - 1

rx  = range(x)  
rv0 = range(v0)
print("rv0 rx ");print(rv0);print(rx)
if (inp == 1 & inorm == 1) {
rx  = quantile(x,c(atrim/2,1-atrim/2),type=1)                 # use "trimmed" extremes instead of true ones if inp == 1 in the normal case
rv0 = c(maxv0min,minv0max)
print("trimmed rv0 rx ");print(rv0);print(rx)
print(table(v0))
             }
evt = rbind(cbind(rv0[1],rx[1]),cbind(rv0[2],rx[1]),cbind(rv0[1],rx[2]),cbind(rv0[2],rx[2]))
X.eval = data.frame(v0=evt[,1],x=evt[,2])                     # used for np estimation of eta at the four min max combinations
                                                              # but should we use four different banwidths?

model = lm(y~v0+x)                                            # why not more simply use v? It is done below and tried for K=2 N=100 norm
#summary(model)                                               # no change in results, but focalvpb = 34 instead of 30 !
gamm0 = gamm = model$coef[c(2,3,1)]  
gamm[3] = gamm[3] - 0.5*gamm[1]    # focal point CORRECTION 27.6.11. This is correct for uniform.
if (inorm == 1) {pi0 = 0.9600785   # focal point CORRECTION 27.6.11 for normal: see focval_pi0_m0.R 
                 gamm[1] = gamm[1]/pi0
                 }
if (1==2) {
model = lm(y~v+x)                                            # why not more simply use v?
#summary(model)
gamm0 = gamm = model$coef[c(2,3,1)]
          }

if (inp == 1 & mmm == 2) {                                    # estimation for the np model; functions of Li/Racine -> np

# 4. Now cv method adapted to local linear

if ( inorm == 0 )        {                                    # uniform case: here we use a shortcut  
	if (mc == 1) {                                          	  # compute scale factor from first mc
	if (bw_new == 0) {
	bw= npregbw(y~ordered(v0)+x,regtype=reg_t,bwmethod=bw_m,ckertype=cker_t,nmulti=nmult,bandwidth.compute = bw_c,
		bwtype=bw_t,bws=bw_imp,bwscaling=FALSE)
	bw1 = bw                                                      # keep for later use
                 }
	if (bw_new == 1) {
	bw= npregbw(y~ordered(v0)+x,regtype=reg_t,bwmethod=bw_m,ckertype=cker_t,nmulti=nmult,bandwidth.compute = bw_c,bwtype=bw_t)
	bw1 = bw                                                      # keep for later use

	# bw1= npregbw(y~ordered(v0)+x,subset=X.eval[1,],regtype=reg_t,bwmethod=bw_m,ckertype=cker_t,nmulti=5,bandwidth.compute = bw_c,bwtype=bw_t)
	# aborted attempt to compute the bandwidth at one evaluation point ...
	                 }
	if (bw$bw[2] > hmax) { bw$bw[2] = 1.06*sd(x)/N^(1/5)
                       print("ROT"); print(bw$bw) 
	                     }
	print("summary(bw)")
	summary(bw)
	# Sleep for 10 seconds so that we can examine the output...
	# Sys.sleep(10)
	scalex = bw$bw[2]/sd(x)
	             }
	if (mc > 1) { bw$bw[2] = scalex*sd(x) }                       # use scale factor for other mcs
                         }                                     # end uniform case

if ( inorm == 1 )        {                                     # here I don't know (yet) how to use a shortcut
	bw= npregbw(y~ordered(v0)+x,regtype=reg_t,bwmethod=bw_m,ckertype=cker_t,nmulti=nmult,bandwidth.compute = bw_c,
		bwtype=bw_t)
	print("summary(bw)")
	summary(bw)

	if (bw$bw[2] > hmax & mc > 1) { bw$bw = bw1$bw
                       print("PREVIOUS BW USED"); 
	                     }
	bw1 = bw                                                 # keep for later use
                         }                                     # end normal  case

model4 = npreg(newdata=X.eval, bws=bw, gradients=TRUE)         # do np only on 4 evaluation points
# summary(model4)
# npplot(bws=bw,plot.errors.method="asymptotic")               # rotating 3D plot : look at this only once!
# npplot(bws=bw)               					   # rotating 3D plot : look at this only once!
# title(sub=paste("mceff=",mceff))

gradient=model4$grad
print("extreme values of gradients")
print(cbind(min(gradient[,1]),max(gradient[,1]),min(gradient[,2]),max(gradient[,2])))

# actually we should drop mc replications for which we do not have monotonicity!

# now evaluate eta at the four combinations of interest
# eta4 = as.matrix(fitted(npreg(newdata=X.eval, bws=bw))) # useful if estimation at all points
eta4 = as.matrix(predict(model4,newdata=X.eval))          # useful if estimation at 4 points of interest
                }                                              # end nonparametric estimation and prediction

print("focal point");print(gamm)

if (mmm == 2) {                                                # MMD estimator

eta_vec = matrix(c(0), nrow = 4, ncol=1)
mat_inq = matrix(c(1), nrow = 4, ncol=3) # inequalities written as mat_inq * c <= eta_vec:
                                         # the first four inequalities written for v0,
                                         # the four next written for v1
for (i in 1:2) {
for (j in 1:2) { eta_vec[j+2*(i-1)] = eta(rx[i],rv0[j],inorm,i,j) 
                 mat_inq[j+2*(i-1),1:2] = cbind(rv0[j],rx[i])
               }
               }
if (inp == 1) { print("cbind(eta_vec,eta4)") ; print(cbind(eta_vec,eta4)) 
                eta_vec = eta4
              }
eta_vec = rbind(eta_vec,-eta_vec)
#print("eta_vec");#print(eta_vec)

mat_inq = rbind(mat_inq, -mat_inq)
mat_inq[5:8,1] = mat_inq[5:8,1] - 1      # c1 multiplied by v1=v0+1 for those inequalities
#print("mat_inq");#print(mat_inq)

# find vertices : for each combination of 3 equations [triple (i<j<k)]
# find intersection of the three planes
# keep it if it satisfies the inequalities

ri = 8         # number of inequalities
vertices = matrix(0,nrow=1, ncol=6)
vertices[1,1:3] = gamm
check = sum(mat_inq %*% gamm <= eta_vec + tol)
if (check < ri) { print(paste("focal point satisfies only ",check," inequalities out of 8")) 
                  print("inequalities: mat_inq*gamm > eta_vec + tol  ")
                  for (i in 1:ri) { 
                    if (mat_inq[i,] %*% gamm > eta_vec[i]+tol) {
                       print(c(i,mat_inq[i,] %*% gamm, eta_vec[i])) 
                                                           }
                                  }
                }
for (i in 1:(ri-2)) {
  for (j in (i+1):(ri-1)) {
    for (k in (j+1):ri) { # print(mat_inq[c(i,j,k),]);#print(det(mat_inq[c(i,j,k),]))
      if (rankMatrix(mat_inq[c(i,j,k),],tol=tolrank)[1]==3) {
        vertex = solve(mat_inq[c(i,j,k),],eta_vec[c(i,j,k)]) 
        #print(cbind(mat_inq %*% vertex,eta_vec + tol))
        #print(mat_inq %*% vertex <= eta_vec + tol)
        if (sum(mat_inq %*% vertex <= eta_vec + tol) == ri) {
            vertices=rbind(vertices,cbind(t(vertex), matrix(c(i, j, k),1,3))) # had to modify this 21.01.11
        }
      }
}}}

}   # End MMD

if (mmm == 1) {                                                # MMM estimator

# v0 indicators : extremes have at least 5% of observations
a = range(v0)
nind =  a[2]-a[1] + 1
vind = matrix(NA, N, nind)
for (i in a[1]:a[2]) { vind[,i-a[1]+1] = (v0 == i ) }
fv0 = colSums(vind)
#print(rbind((a[1]:a[2]),colSums(vind)))

while (fv0[1] < N*0.05) {
a[1] = a[1]+1
vind[,2] = (v0 <= a[1])
vind = vind[,2:nind]
nind = ncol(vind)
fv0 = colSums(vind)
#print(rbind((a[1]:a[2]),colSums(vind)))
                        }     

while (fv0[nind] < N*0.05) {
a[2] = a[2]-1
vind[,(nind-1)] = (v0 >= a[2])
vind = vind[,1:(nind-1)]
nind = ncol(vind)
fv0 = colSums(vind)
#print(rbind((a[1]:a[2]),colSums(vind)))
                        }   
print(rbind((a[1]:a[2]),colSums(vind)))

# definition of w instruments: v0 indicators + abs(x)^k, k=1, ..., K

if (K == 0) { w = vind } 
if (K > 0)  { w = cbind(vind,abs(x)) }                 # use only abs(x) and indicators
if (K > 1)  { 
for (k in 2:K) { w = cbind(w,abs(x)^k) }                 # add powers of abs(x)
            }
w
#w1[1,]
#w2[1,]
#w3[1,]

eta_vec = eta_vec_np = matrix(c(0), nrow = 4, ncol=1)
mat_inq = matrix(c(1), nrow = 4, ncol=3) # inequalities written as mat_inq * c <= eta_vec:
                                         # the first four inequalities written for v0,
                                         # the four next written for v1

eta_vec = ColM(y,w);
eta_vec = c(eta_vec,-eta_vec);
#print("eta_vec"); print(eta_vec)
rw = dim(w)[2]                                    # number of instruments
ri = 2*rw                                         # number of inequalities
mat_inq = cbind(ColM(v0,w),ColM(x,w),colMeans(w)) # inequalities written as mat_inq * c <= eta_vec:
                                                  # the first rw inequalities written for v0,
mat_inq = rbind(mat_inq, -mat_inq)
mat_inq[(rw+1):(2*rw),1] = -ColM(v1,w)            # the next rw written for v1

# Scale mat_inq and eta_vec
for (i in 1:(dim(mat_inq)[1])) {
eta_vec[i] = eta_vec[i]/sum(abs(mat_inq[i,]))
mat_inq[i,] = mat_inq[i,]/sum(abs(mat_inq[i,]))
if (eta_vec[i]=="NaN") {eta_vec[i]=0;mat_inq[i,]=rep(0,3)}
}

#print("mat_inq");print(mat_inq)

# find vertices : for each combination of 3 equations [triple (i<j<k)]
# find intersection of the three planes
# keep it if it satisfies the inequalities

vertices = matrix(0,nrow=1, ncol=6)
vertices[1,1:3] = gamm
check = sum(mat_inq %*% gamm <= eta_vec + tol)
if (check < ri) { print(paste("focal point satisfies only ",check," inequalities out of ",ri))
                  print("inequalities: mat_inq*gamm > eta_vec + tol  ")
                  for (i in 1:ri) { 
                    if (mat_inq[i,] %*% gamm > eta_vec[i]+tol) {
                       print(paste(i," ",mat_inq[i,] %*% gamm," > ", eta_vec[i])) 
                                                           }
                                  }
                }
for (i in 1:(ri-2)) {
  for (j in (i+1):(ri-1)) {
    for (k in (j+1):ri) { 
    #print(mat_inq[c(i,j,k),]);print(det(mat_inq[c(i,j,k),]))
      if (rankMatrix(mat_inq[c(i,j,k),],tol=tolrank)[1]==3) {
        vertex = solve(mat_inq[c(i,j,k),],eta_vec[c(i,j,k)])
        #print("vertex"); #print(vertex)
        #print(cbind(mat_inq %*% vertex,eta_vec + tol))
        #print(mat_inq %*% vertex <= eta_vec + tol)
        if (sum(mat_inq %*% vertex <= eta_vec + tol) == ri) {
        #print("1");
            vertices=rbind(vertices,cbind(t(vertex), matrix(c(i, j, k),1,3)))
        }
      }
}}}

}   # End MMM

# print("dim(vertices)"); # print(dim(vertices))
# print("vertices"); # print(vertices)          
if (1==2) { # skip plots
plot(vertices[,c(1,2)])
abline(v=gamm[1],h = gamm[2], col = 2)
plot(vertices[,c(1,3)])
abline(v=gamm[1],h = gamm[3], col = 2)
plot(vertices[,c(2,3)])
abline(v=gamm[2],h = gamm[3], col = 2)
scatterplot3d(vertices[,1:NR],highlight.3d=TRUE)
         } # end skip plots
#if (mmm == 1) {save.image(file = paste(method[mmm],"_N",N,"MC",MC,"K",K,dist,"MA_.RData",sep = "")) }    # remove when correcly functioning
#if (mmm == 2) {save.image(file = paste(method[mmm],"_N",N,"MC",MC,dist,"MA_.RData",sep = ""))    }       # remove when correcly functioning

if ((dim(vertices)[1] > 1)|(check==ri)) {   # exclude case where single line in vertices corresponds to focal point
                                            # AND focal point does not satisfy all constraints
if (check<ri) {                             # focal point not in set !? 
vertices = vert(unique(vertices[2:dim(vertices)[1],1:NR]))
focvalpb =  focvalpb + 1
              }
if (check==ri) {                            # focal point is in set 
vertices = vert(unique(vertices[,1:NR])) 
               }
print("vertices");  print(vertices)

# Store mc results
for (i in 1:NR) { array.r_[mc,i,] = range(vertices[,i]) } # bounds in the three dimensions
#print(array.r_)
           numrows = min(MaxVert,dim(vertices)[1])
           array.mvertices[mc,(1:dim(vertices)[1]),] = vertices[,1:NR]
           focval[mc,1:NR]=gamm
mc=mc+1
dvertmin = min(dvertmin,dim(vertices)[1])
                                          } # end if dim(vertices)
} # end mc loop
print("End MC loop")
print(paste("dvertmin ",dvertmin))
print(paste("number of simulations with focal point outside set : ",focvalpb))
print(paste("number of rejected simulations (no vertex found and focal point outside set : ",(mceff-MC)))
if (inp == 1) { print("bw1$bw");print(bw1$bw) }
#---------------------------------------------------- HERE =============================================================
if (1==2) {
if (inp == 0 & mmm == 1) {
load(file = paste(method[mmm],"_N",N,"MC",MC,"K",K,dist,"tol",TOL,"MA.RData",sep = ""))
              }
if (inp == 0 & mmm == 2) {
load(file = paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,"MA.RData",sep = ""))
              }
          }
r_arr = matrix(NA,NR,2)
for (i in 1:NR) {r_arr[i,]  = range(array.r_[,i,])}
print("range"); print(r_arr)
# which(array.r_[,1,1]==0,arr.ind = TRUE)     # useful only in case one of the mc vertices contains NA


plot.new()
# in all cases we produce three plots
if (inorm == 0) {  # begin uniform case

# Define 2D Projections of True Identified Set
true = matrix(c(1.000000, -1.0, 1.5000000,
     1.000000, -1.2, 1.5000000,
     0.800000, -1.0, 1.1000000,
     1.333333, -1.0, 0.8333333,
     1.000000, -0.8, 0.5000000,
     1.000000, -1.0, 0.5000000),nrow=6,ncol=3,byrow = TRUE)

true12 = rbind(c(1.0,-1.2),c(1.0,-0.8),c(1.333333,-1.0),c(0.8,-1.0))
true13 = rbind(c(0.8,1.1),c(1,1.5),c(1.0,0.5),c(1.333333,0.8333333))
true23 = rbind(c(-0.8,0.5),c(-1,0.5),c(-1.0,1.5),c(-1.2,1.5))
true12 = ccw(true12);true13 = ccw(true13);true23 = ccw(true23)
               }   # end uniform case

if (inorm == 1) {   # begin normal case 

true   =  matrix(c(0.9721836077,-1.000000000,0.5206604371,
    0.9721836077,-1.000000000,1.4928440449,
    0.9721836077,-1.000000486,1.0067522411,
    0.9721836077,-0.999999514,1.0067522411,
    1.0163737717,-1.000000000,1.0067522411,
    0.9316759574,-1.000000000,1.0067522411
    ),nrow=6,ncol=3,byrow = TRUE)

true12 = ccw(true[,c(1,2)]);true13 = ccw(true[,c(1,3)]);true23 = ccw(true[,c(2,3)])
               }   # end normal case


par(mfrow = c(2,2),lab=c(1,1,7),las=2)     # get the three plots in one graph (their ordering and labeling still requires some care)
                                           # here we focus on two-dimensional projections
t = 1
for (j in 1:2) { for (k in (j+1):3) {      # begin for (j in 1:2) { for (k in (j+1):3)

for (mc in 1:MC) {
fill = vert(filtna(unique(array.mvertices[mc,,c(j,k)])))
array.mvertices2D[mc,1:dim(fill)[1],,t] = fill
} # end MC-loop


# Obtain Minkowski average
# --------------------------------------------------------
# Minkowski Sum
MA = mks(filtna(unique(array.mvertices2D[1,,,t])),filtna(unique(array.mvertices2D[2,,,t])))
if (MC > 2) {
for (mc in 3:MC) {               # mc loop
  MM = filtna(unique(array.mvertices2D[mc,,,t]))
  MA = mks(MM,unique(MA))  
                }              # end MC loop
             }
# Minkowski Average
MA = MA/MC
print("range(MA[,1], range(MA[,2]")
print(range(MA[,1]))
print(range(MA[,2]))

# --------------------------------------------------------
clab = c("c1","c2","c3")
# if (1==2) {                              # skip mc-plot loop
for (mc in 1:MC) {                         # mc-plot loop
    vertices2D = filtna(unique(array.mvertices2D[mc,,,t]))
    d1 = dim(vertices2D)[1]
    v  = ccw(vertices2D)
    if (d1 > 2) { cv = convhulln(vertices2D,"FA") ; 
                  v = ccw(vertices2D[unique(c(cv$hull[,1],cv$hull[,2])),]) }      
    v = rbind(v,v[1,])
    v1 = v[1:2,]                             # vertices defining segment 1
    
    if (mc==1) {
#     if (t==1) {
        plot(v1[,1],v1[,2],type="l",col="pink",xlim=r_arr[j,],ylim=r_arr[k,],xlab=clab[j],ylab=clab[k])
#               }
#     if (t>1) {
#         bla=r_arr[k,]
#         plot(v1[,1],v1[,2],type="l",col="pink",xlim=r_arr[j,],ylim=c(bla[1]-0.1*bla[1],bla[2]+0.1*bla[2]),xlab=clab[j],ylab=clab[k])
#               }
    if (t==1) {
        if (inp == 0 & mmm == 1) { title(paste(method[mmm],": N=",N," MC=",MC," K=",K," ",dist,sep = "")) }
        if (inp == 0 & mmm == 2) { title(paste(method[mmm],": N=",N," MC=",MC," ",dist,sep = "")) }
        if (inp == 1) { title( main = paste(method[mmm],": N=",N," MC=",MC," ",dist," NP",sep = "") ) }
             }
               }
    if (mc>1) {
    lines(v1[,1],v1[,2],type="l",col="pink")
              }
    for (is in 2:dim(v)[1]) {        # loop over segments
        v1 = v[c(is-1,is),]
        lines(v1[,1],v1[,2], type="l",col="pink") 
        }                                  # end loop over segments
    axis(1,at=round(c(range(array.r_[,j,]),range(MA[,1])), digits=3))
    axis(2,at=round(c(range(array.r_[,k,]),range(MA[,2])), digits=3))
    }                                      # end mc-plot loop
#         }                                # end skip mc-plot loop
# Plot Minkowski Average
    cv = convhulln(MA,"FA") 
    surf_MA = cv$vol
    v1 = MA[cv$hull[1,],]      # vertices defining segment 1
            lines(v1[,1],v1[,2],type="l",xlim=cbind(min(range(MA[,1]),r_arr[1,]),
            max(range(MA[,1]),r_arr[1,])),ylim=cbind(min(range(MA[,2]),r_arr[k,]),max(range(MA[,2]),r_arr[k,])),asp=1)
    for (is in 1:dim(cv$hull)[1]) {        # vertices defining segment i
        v1 = MA[cv$hull[is,],]
        lines(v1[,1],v1[,2], type="l")
                                 }         # end loop over segments
# Plot True Sets
if (t==1) {
    cv = convhulln(true12,"FA") 
    surf_true = cv$vol
    v1 = true12[cv$hull[1,],]      # vertices defining segment 1
            lines(v1[,1],v1[,2],type="l",col="green",xlim=cbind(min(range(true12[,1]),r_arr[1,]),
            max(range(true12[,1]),r_arr[1,])),ylim=cbind(min(range(true12[,2]),r_arr[k,]),max(range(true12[,2]),r_arr[k,])),asp=1)
    for (is in 1:dim(cv$hull)[1]) {        # vertices defining segment i
        v1 = true12[cv$hull[is,],]
        lines(v1[,1],v1[,2], type="l",col="green")
                                 }         # end loop over segments
        title(sub=paste("surfaces ",round(surf_MA,digits=3),
           round(surf_true,digits=3)))
          }
if (t==2) {
    cv = convhulln(true13,"FA") 
    surf_true = cv$vol
    v1 = true13[cv$hull[1,],]      # vertices defining segment 1
            lines(v1[,1],v1[,2],type="l",col="green",xlim=cbind(min(range(true13[,1]),r_arr[1,]),
            max(range(true13[,1]),r_arr[1,])),ylim=cbind(min(range(true13[,2]),r_arr[k,]),max(range(true13[,2]),r_arr[k,])),asp=1)
    for (is in 1:dim(cv$hull)[1]) {        # vertices defining segment i
        v1 = true13[cv$hull[is,],]
        lines(v1[,1],v1[,2], type="l",col="green")
                                 }         # end loop over segments
        title(sub=paste("surfaces ",round(surf_MA,digits=3),
           round(surf_true,digits=3))) 
          }
if (t==3) {
    cv = convhulln(true23,"FA") 
    surf_true = cv$vol
    v1 = true23[cv$hull[1,],]      # vertices defining segment 1
            lines(v1[,1],v1[,2],type="l",col="green",xlim=cbind(min(range(true23[,1]),r_arr[1,]),
            max(range(true23[,1]),r_arr[1,])),ylim=cbind(min(range(true23[,2]),r_arr[k,]),max(range(true23[,2]),r_arr[k,])),asp=1)
    for (is in 1:dim(cv$hull)[1]) {        # vertices defining segment i
        v1 = true23[cv$hull[is,],]
        lines(v1[,1],v1[,2], type="l",col="green")
                                 }         # end loop over segments
   title(sub=paste("surfaces ",round(surf_MA,digits=3),
           round(surf_true,digits=3)))

          }
# -------------------------------------------------
# Compute Hausdorff Distances to true projections, compute quantiles, plot average sets in quantiles or convex hulls
#if (1==2) { # begin skip Hausdorff
for (mc in 1:MC) {
    vertices2D = filtna(unique(array.mvertices2D[mc,,,t]))
    vertices2D = ccw(vertices2D)
    if (t==1) {
    HD0  = HD(true12,vertices2D)
    dHD0 = dHD(true12,vertices2D)
    dHD1 = dHD(vertices2D,true12)
    }
    if (t==2) {
    HD0  = HD(true13,vertices2D)
    dHD0 = dHD(true13,vertices2D)
    dHD1 = dHD(vertices2D,true13)
    }
    if (t==3) {
    HD0  = HD(true23,vertices2D)
    dHD0 = dHD(true23,vertices2D)
    dHD1 = dHD(vertices2D,true23)
    }
    array.HD[mc,,t] = cbind(HD0,mc)
    array.dHD0[mc,,t] = cbind(dHD0,mc)
    array.dHD1[mc,,t] = cbind(dHD1,mc)
    }  # End MC-loop
    array.HDsorted[,,t] = array.HD[order(array.HD[,1,t]),,t]
    HDquant[t,] = quantile(array.HDsorted[,1,t], p)
#} # end skip Hausdorff

#if (1==2) { # skip quantile averages loop
# Generate sets within desired quantiles of HD distances from True, compute and plot Minkowski averages for all quantiles


for (q in 1:ncol(p)) { # begin quantile loop
    if (q == 1) { 
        MCquant = array.mvertices2D[array.HDsorted[which(array.HDsorted[,1,t]<=HDquant[t,q]),2,t],,,t]
    }
    if (q > 1 ) { 
#       MCquant = array.mvertices2D[array.HDsorted[which(array.HDsorted[,1,t]<=HDquant[t,q] & array.HDsorted[,1,t]>HDquant[t,q-1]),2,t],,,t]
        MCquant = array.mvertices2D[array.HDsorted[which(array.HDsorted[,1,t] >HDquant[t,q]),2,t],,,t]

    }
	# Minkowski Sum
    MA = mks(filtna(unique(MCquant[1,,])),filtna(unique(MCquant[2,,])))
    if (dim(MCquant)[1] > 2) {  # MCquant if-loop
        for (mc in 3:dim(MCquant)[1]) {               # mc loop
        MA = mks(filtna(unique(MCquant[mc,,])),filtna(unique(MA)))
                }              # end MC loop    # Minkowski Sum
        }			# end MCquant if-loop
    # Minkowski Average
    MA = MA/dim(MCquant)[1]

    # Plot Minkowski Average
    cv = convhulln(MA,"FA") 
    v1 = MA[cv$hull[1,],]      # vertices defining segment 1
            lines(v1[,1],v1[,2],type="l",lty=2*q+q%%2,xlim=cbind(min(range(MA[,1]),r_arr[1,]),
            max(range(MA[,1]),r_arr[1,])),ylim=cbind(min(range(MA[,2]),r_arr[k,]),max(range(MA[,2]),r_arr[k,])),asp=1)
    for (is in 1:dim(cv$hull)[1]) {        # vertices defining segment i
        v1 = MA[cv$hull[is,],]
        lines(v1[,1],v1[,2], type="l",lty=2*q+q%%2)
                                 }         # end loop over segments
} # end quantile-loop
#} # end skip quantile averages loop
#------------------------------------------------------------------------------------

t = t+1                                    # t = 1 if j=1, k=2, 2 if j=1, k=3, 3 if j=2, k=3                    
}}                                         # end j,k loops
par(mfrow = c(1,1),lab=c(5,5,7))

volumemc=matrix(NA,MC,1)
t = 0
for (mc in 1:MC) {
smallestconvset = filtna(array.mvertices[mc,,])
if (dim(smallestconvset)[1]>2) {
volcv = convhulln(smallestconvset,"FA")   #unique(array....
volumemc[mc] = volcv$vol
t = t+1
}
}
mvolume=sum(filtna(volumemc))/t
print("Mean Volume");print(mvolume)
# print("Volume: Mean block, mean octahedron, mean ratio octahedron/block");print(colMeans(mvolume[,1:3]))
# print("Volume: Standard Deviation");print(sd(marea[,1]))

cker_t
print("means, sd, CI of range of each component and CI of volume")
for (ic in 1:3) {
  for (im in 1:2) {
print(c(ic,mean(array.r_[,ic,im]),sd(array.r_[,ic,im]),quantile(array.r_[,ic,im], cbind(0.05,0.95))))
                  }
                }
print("CI of volume");print(quantile(filtna(volumemc), cbind(0.05,0.95)))

volumemc_cube = matrix(NA,MC,1)
for (mc in 1:MC) { 
smallestconvset = filtna(array.mvertices[mc,,])
if (dim(smallestconvset)[1]>2) {
	volumemc_cube[mc] = prod(array.r_[mc,,2]-array.r_[mc,,1]) 
	}
	}
volovervolcube = filtna(volumemc)/filtna(volumemc_cube)
mvolume_cube=sum(filtna(volumemc_cube))/t
print("Mean of Volume over Cube Volume");print(mean(volovervolcube));print("CI:");print(quantile(volovervolcube,cbind(0.05,0.95)))

print(tstart)
print(date())
tend = date()
print("cpu time in seconds")
proc.time() - cpu0
cpu1 = proc.time()

print(paste("number of simulations with focal point outside set : ",focvalpb))
print(paste("number of rejected simulations (no vertex found and focal point outside set : ",(mceff-MC)))

if (mac == 0) {
if (inp == 1) {
print("bw1$bw");print(bw1$bw)

savePlot(filename = paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,reg_t,bw_t,cker_t,bw_m,".pdf",sep = ""),
         type = c("pdf"),device = dev.cur(),restoreConsole = TRUE) 
save.image(file = paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,reg_t,bw_t,cker_t,bw_m,".RData",sep = ""))
              }

if (inp == 0 & mmm == 1) {
savePlot(filename = paste(method[mmm],"_N",N,"MC",MC,"K",K,dist,"tol",TOL,"MA.pdf",sep = ""),
         type = c("pdf"),device = dev.cur(),restoreConsole = TRUE) 
save.image(file = paste(method[mmm],"_N",N,"MC",MC,"K",K,dist,"tol",TOL,"MA.RData",sep = ""))
              }
if (inp == 0 & mmm == 2) {
savePlot(filename = paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,"MA.pdf",sep = ""),
         type = c("pdf"),device = dev.cur(),restoreConsole = TRUE) 
save.image(file = paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,"MA.RData",sep = ""))
              }
              } #end mac plot loop
if (mac == 1) {
if (inp == 1) {
print("bw1$bw");print(bw1$bw)

# quartz.save(paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,reg_t,bw_t,cker_t,bw_m,".pdf",sep = ""),
#          type = "pdf") 
save.image(file = paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,reg_t,bw_t,cker_t,bw_m,".RData",sep = ""))
              }

if (inp == 0 & mmm == 1) {
# quartz.save(paste(method[mmm],"_N",N,"MC",MC,"K",K,dist,"tol",TOL,"MA.pdf",sep = ""),
#          type = "pdf") 
save.image(file = paste(method[mmm],"_N",N,"MC",MC,"K",K,dist,"tol",TOL,"MA.RData",sep = ""))
              }
if (inp == 0 & mmm == 2) {
# quartz.save(paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,"MA.pdf",sep = ""),
#          type = "pdf") 
save.image(file = paste(method[mmm],"_N",N,"MC",MC,dist,"tol",TOL,"MA.RData",sep = ""))
              }
            } #end mac plot loop

#} # Exclude loop
