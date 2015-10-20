rm(list=ls())
#------------------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------------------
b_h<-0.05			    # Birth rate of hosts
mu_h<-0.045				# Natural death rate of hosts
alpha_h<-0.05			# Per-parasite rate of parasite-induced host mortality

rho<-0.8				# Instantaneous rate of uptake
mu_p<-0.15				# Parasite death
tau_p<-5				# Length of prepatent period
D_p<-exp(-mu_p*tau_p)	# Porportion surviving prepatent period

lambda<-2				# Birth (shedding?) rate of larvae
tau_L<-30				# Development period to infectious stage
mu_L<-0.01					# Death of larvae
D_L<-exp(-mu_L*tau_L)

c.hat<-1				# Migration speed in the absence of parasites
gamma<-0.01				# Per-parasite rate of decline in the migration speed

k<-0.5					# Parasite aggregation

#------------------------------------------------------------------------------------------
# Grid set up
#------------------------------------------------------------------------------------------
x<-c(0:1200)
t<-c(0:1000)
n.x<-length(x)
n.t<-length(t)

H<-matrix(nrow=n.t, ncol=n.x)
P<-matrix(nrow=n.t, ncol=n.x)
L<-matrix(nrow=n.t, ncol=n.x)

# Initial conditions
init<-(max(tau_L, tau_p)+1) #Number of timesteps to initiate
sig<-12; mu<-50
H[1:init,]<-matrix(rep(10000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)
P[1:init,]<-matrix(rep(0.1*H[1,], init),nrow=init, ncol=n.x, byrow=TRUE) #1% prevalence to start
L[1:init,]<-0

#Plot initial conditions
plot(x, H[1,]/max(H, na.rm=TRUE), "l", ylab="", xlab="space", bty="l", yaxt="n", ylim=c(0,1), lwd=3)
lines(x, P[1,]/max(P, na.rm=TRUE), col=2, lwd=2)
lines(x, L[1,], col=3)
	

#------------------------------------------------------------------------------------------
# Simulation
#------------------------------------------------------------------------------------------

for(i in (init+1):n.t){
	C<-round(c.hat)
	H.last<-H[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	P.last<-P[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	
	H[i,]<-H.last+(b_h-mu_h)*H.last-alpha_h*P.last
	
	agg.term<-(P.last/H.last+(P.last^2/H.last^2*(k+1)/k))
	agg.term[is.na(agg.term)]<-0
	
	P[i,]<-P.last
			+rho*D_p*L[(i-1-tau_p),c((n.x+1-C-tau_p):n.x, 1:(n.x-C-tau_p))]*H[(i-1-tau_p),c((n.x+1-C-tau_p):n.x, 1:(n.x-C-tau_p))]
			-(mu_p+mu_h)*P.last
			-alpha_h*H.last*agg.term
	
	L[i,]<-L[i-1,]+lambda*D_L*P[i-1-tau_L,c((n.x+1-C-tau_L):n.x, 1:(n.x-C-tau_L))]-mu_L*L[i-1,]-rho*L[i-1]*H.last
	
}

cat("\n Hosts:", range(H), "\n Parasites", range(P), "\n Larvae", range(L))
#------------------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------------------
par(mfrow=c(3,1))
plot(x, H[1,]/max(H, na.rm=TRUE), "n", ylab="", xlab="space", bty="l", ylim=c(0,1))
for(i in 1:19){
	t<-round(seq(1, n.t, length.out=19))[i]
	lines(x, H[t,]/max(H, na.rm=TRUE), lwd=3, col=paste("#000000", c("05", seq(10,95,5)), sep="")[i])
	}
	
plot(x, P[1,]/max(P, na.rm=TRUE), "n", ylab="", xlab="space", bty="l", ylim=c(0,1))
for(i in 1:19){
	t<-round(seq(1, n.t, length.out=19))[i]
	lines(x, P[t,]/max(P, na.rm=TRUE), lwd=3, col=paste("#FF0000", c("05", seq(10,95,5)), sep="")[i])
	}

plot(x, L[1,]/max(L, na.rm=TRUE), "n", ylab="", xlab="space", bty="l", ylim=c(0,1))
for(i in 1:19){
	t<-round(seq(1, n.t, length.out=19))[i]
	lines(x, L[t,]/max(L, na.rm=TRUE), lwd=3, col=paste("#00FF00", c("05", seq(10,95,5)), sep="")[i])
	}
#------------------------------------------------------------------------------------------
# Animation
#------------------------------------------------------------------------------------------
require(animation)

saveLatex({
	par(mfrow=c(3,1), mar=c(4,4,1,1), oma=c(3,1,2,0))
	ani.options(interval=0.06)
	for(i in seq(6,n.t,20)){
		plot(x, H[i,], "l", ylim=c(0, max(H)), xlab="", ylab="Hosts")
		plot(x, P[i,], "l", col=2, ylim=c(0, max(P)), xlab="", ylab="Parasites")
		plot(x, L[i,], "l", col=3, ylim=c(0, max(L)), xlab="", ylab="Larvae")
		mtext(side=1, outer=TRUE, "Migration route (x)")
		mtext(side=3, outer=TRUE, paste("Time step (t=", i, ")", sep=""))
		}
}, img.name="ParasiteMigration", outdir="~/Documents/Migration Speed/Molnar/Animation", ani.dev = 'pdf', ani.type = 'pdf', ani.height = 6, ani.width = 8, ani.opts='controls,width=8in', pdflatex = '/usr/texbin/pdflatex', caption=c("Abundance of hosts (black), associated parasites (red), and 'free-living' parasite larvae (green) over a migration route from x=0 to x=500.  Hosts (and associated parasites) are migrating at c=1 unit space/unit time. The rest of the model is similar to Molar et al. (2013) but with an equation for host abundance."), overwrite=TRUE, documentclass = paste("\\documentclass{article}", "\\usepackage[landscape,margin=0.3in]{geometry}", sep="\n")
)

#------------------------------------------------------------------------------------------
# Different host states
#------------------------------------------------------------------------------------------

# Additional parameters
b_h<-0.05			    # Birth rate of hosts
mu_h<-0				# Natural death rate of hosts
alpha_h<-0.02			# Per-parasite rate of parasite-induced host mortality

rho<-0.8				# Instantaneous rate of uptake
mu_P<-0.15				# Parasite death
tau_P<-2				# Length of prepatent period
D_P<-exp(-mu_P*tau_P)	# Porportion surviving prepatent period

lambda<-2				# Birth (shedding?) rate of larvae
tau_L<-5				# Development period to infectious stage
mu_L<-0.01					# Death of larvae
D_L<-exp(-mu_L*tau_L)

c.hat<-1				# Migration speed in the absence of parasites
k<-0.5

lambda_0<-0.05
lambda_1<-0.07 # because lambda_1>lambda_0, the population should slow down over time

#------------------------------------------------------------------------------------------

# Setup matrices
x<-c(0:1200)
t<-c(0:1000)
n.x<-length(x)
n.t<-length(t)

H1<-matrix(nrow=n.t, ncol=n.x)
H0<-matrix(nrow=n.t, ncol=n.x)
P1<-matrix(nrow=n.t, ncol=n.x)
P0<-matrix(nrow=n.t, ncol=n.x)
L<-matrix(nrow=n.t, ncol=n.x)

# Initial conditions
init<-(max(tau_L, tau_p)+1) #Number of timesteps to initiate
sig<-12; mu<-50
H1[1:init,]<-matrix(rep(5000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)
H0[1:init,]<-matrix(rep(5000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)

P1[1:init,]<-matrix(rep(0.1*H1[1,], init),nrow=init, ncol=n.x, byrow=TRUE) #1% prevalence to start
P0[1:init,]<-matrix(rep(0.1*H0[1,], init),nrow=init, ncol=n.x, byrow=TRUE) #1% prevalence to start

L[1:init,]<-0

#------------------------------------------------------------------------------------------

# Simulation
for(i in (init+1):n.t){
	C<-round(c.hat)
	H1.last<-H1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	H0.last<-H0[i-1,]
	P1.last<-P1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	P0.last<-P0[i-1,]
	
	H1[i,]<-H1.last+(b_h-mu_h)*H1.last - alpha_h*P1.last - lambda_1*H1.last + lambda_0*H0.last
	H0[i,]<-H0.last+(b_h-mu_h)*H0.last - alpha_h*P0.last - lambda_0*H0.last + lambda_1*H1.last
	
	agg.term1<-(P1.last/H1.last+(P1.last^2/H1.last^2*(k+1)/k))
	agg.term1[is.na(agg.term1)]<-0
	agg.term0<-(P0.last/H0.last+(P0.last^2/H0.last^2*(k+1)/k))
	agg.term0[is.na(agg.term0)]<-0
	
	P1[i,]<-P1.last
			+rho*D_p*L[(i-1-tau_p), c((n.x+1-C-tau_p):n.x, 1:(n.x-C-tau_p))]*H1[(i-1-tau_p),c((n.x+1-C-tau_p):n.x, 1:(n.x-C-tau_p))]
			-(mu_p+mu_h)*P1.last
			-alpha_h*H1.last*agg.term1
			-lambda_1*P1.last + lambda_0*P0.last
			
	P0[i,]<-P0.last
			+rho*D_p*L[(i-1-tau_p),]*H0[(i-1-tau_p),]
			-(mu_p+mu_h)*P0.last
			-alpha_h*H0.last*agg.term0
			-lambda_0*P0.last + lambda_1*P1.last
	
	L[i,]<-L[i-1,]+lambda*D_L*(P1[i-1-tau_L, c((n.x+1-C-tau_L):n.x, 1:(n.x-C-tau_L))]+P0[i-1-tau_L,])-mu_L*L[i-1,]-rho*L[i-1]*(H1.last+H0.last)
	
	H0[i,][which(H0[i,]<1)]<-0
	H1[i,][which(H1[i,]<1)]<-0
	P0[i,][which(P0[i,]<1)]<-0
	P1[i,][which(P1[i,]<1)]<-0
	L[i,][which(L[i,]<1)]<-0
}

cat("\n H1:", range(H1), " H0:", range(H0), "\n P1", range(P1), " P0", range(P0), "\n Larvae", range(L))

#------------------------------------------------------------------------------------------
par(mfrow=c(3,1),mar=c(4,5,2,2))

#for(i in c(100,120,150,180,200)){
plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.2)
lines(x, H1[i,], col=3)
lines(x, H0[i,], col=2)

plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.2)
lines(x, P1[i,], col=3)
lines(x, P0[i,], col=2)

plot(x, L[i,], "l", col=4)
#}