#Constant background source of parasites
b_H<-0#0.05			    # Birth rate of hosts
mu_H<-0.002				# Natural death rate of hosts
alpha<-0.001				# Per-parasite rate of parasite-induced host mortality

lambda_1<-0.06
lambda_0<-0.04

rho<-0.001				# Instantaneous rate of uptake
mu_A<-0.2				# Non-reproductive (A) parasite death
tau_A<-10				# Developmental time of parasites to reproductive stage in host

mu_P<-0.10				# Reproductive (P) Parasite death

beta<-5					# Birth (shedding?) rate of larvae
mu_E<-0.2				# Death of eggs
tau_E<-20				# Development time of eggs

mu_I<-0.3				# Death of infective larvae

c.hat<-1				# Migration speed in the absence of parasites
k<-0.05					# Parasite aggregation


# Setup matrices
x<-c(0:1200)
t<-c(0:1000)
n.x<-length(x)
n.t<-length(t)

#Hosts
H1<-matrix(nrow=n.t, ncol=n.x)
H0<-matrix(nrow=n.t, ncol=n.x)
#Non-reproductive parasites
A1<-matrix(nrow=n.t, ncol=n.x)
A0<-matrix(nrow=n.t, ncol=n.x)
#Reproductive parasites
P1<-matrix(nrow=n.t, ncol=n.x)
P0<-matrix(nrow=n.t, ncol=n.x)
#Eggs (uninfective)
E<-matrix(nrow=n.t, ncol=n.x)
#Infective larvae
I<-matrix(nrow=n.t, ncol=n.x)

# Initial conditions
init<-tau_A+tau_E #Number of timesteps to initiate
sig<-12; mu<-50
H1[1:init,]<-matrix(rep(10000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)
H0[1:init,]<-matrix(rep(10000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)
A1[1:init,]<-0
A0[1:init,]<-0
P1[1:init,]<-0
P0[1:init,]<-0
E[1:init,]<-100
I[1:init,]<-0

# Simulation
for(i in (init+1):n.t){
	C<-round(c.hat)
	x.i<-c((n.x+1-C):n.x, 1:(n.x-C))
	
	H1[i,]<-H1[i-1, x.i] + (b_H-mu_H)*H1[i-1, x.i] - alpha*P1[i-1, x.i] - lambda_1*H1[i-1, x.i] + lambda_0*H0[i-1,]
	
	H0[i,]<-H0[i-1,] + (b_H-mu_H)*H0[i-1,] - alpha*P0[i-1,] - lambda_0*H0[i-1,] + lambda_1*H1[i-1, x.i]
	
	A.agg.term1<-P1[i-1, x.i]*A1[i-1, x.i]/H1[i-1, x.i]
	A.agg.term1[is.na(A.agg.term1)]<-0
	
	A1[i,]<-A1[i-1,x.i] + rho*I[i-1,x.i]*H1[i-1,x.i] - (mu_A + mu_H)*A1[i-1,x.i] - (1/tau_A)*A1[i-1, x.i] - lambda_1*A1[i-1, x.i] + lambda_0*A0[i-1,] - alpha*A.agg.term1
	
	A.agg.term0<-P0[i-1,]*A0[i-1,]/H0[i-1,]
	A.agg.term0[is.na(A.agg.term0)]<-0
	A0[i,]<-A0[i-1,]+rho*I[i-1,]*H0[i-1,] - (mu_A + mu_H)*A0[i-1,] - (1/tau_A)*A0[i-1,] + lambda_1*A1[i-1, x.i] - lambda_0*A0[i-1,] - alpha*A.agg.term0
	
	agg.term1<-(P1[i-1,x.i]/H1[i-1,x.i]+(P1[i-1,x.i]^2/H1[i-1,x.i]^2*(k+1)/k))
	agg.term1[is.na(agg.term1)]<-0
	agg.term0<-(P0[i-1,]/H0[i-1,]+(P0[i-1,]^2/H0[i-1,]^2*(k+1)/k))
	agg.term0[is.na(agg.term0)]<-0
	
	P1[i,]<- P1[i-1,x.i] + (1/tau_A)*A1[i-1, x.i] - (mu_P+mu_H)*P1[i-1, x.i] - alpha*H1[i-1, x.i]*agg.term1 - lambda_1*P1[i-1,x.i] + lambda_0*P0[i-1,]
			
	P0[i,]<- P0[i-1,] + (1/tau_A)*A0[i-1,] - (mu_P+mu_H)*P0[i-1,] - alpha*H0[i-1,]*agg.term0 - lambda_0*P0[i-1,]  + lambda_1*P1[i-1,x.i]
			
	E[i,] <- E[i-1,]+ beta*(P0[i-1,]+P1[i-1,x.i]) - mu_E*E[i-1,] - (1/tau_E)*E[i-1,]
	
	I[i,] <- I[i-1,] + (1/tau_E)*E[i-1,] - mu_I*I[i-1,] - rho*I[i-1,]*(H0[i-1,]+H1[i-1,x.i]) 
	}

cat("\n H1:", range(H1), " H0:", range(H0), "\n A1: ", range(A1), " A0: ", range(A0), "\n P1: ", range(P1), " P0: ", range(P0), "\n Eggs: ", range(E), " I: ", range(I))



for(i in c(12,13,18,20,21,22,26)){
par(mfrow=c(3,2),mar=c(4,5,2,2))
plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.2, main=i)
lines(x, H1[i,], col=3)
lines(x, H0[i,], col=2)

plot(x, A0[i,]+A1[i,], "l", col=grey(0.8), lwd=1.2)
lines(x, A1[i,], col=3)
lines(x, A0[i,], col=2)

plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.2)
lines(x, P1[i,], col=3)
lines(x, P0[i,], col=2)

plot(x, E[i,], "l", col=4)
plot(x, I[i,], "l", col=4)

}
require(animation)

saveLatex({
	ani.options(interval=0.06)
	for(i in seq(init,n.t,20)){
		layout(matrix(c(1,2,3,4,1,2,3,5), nrow=4, ncol=2))
		par(mar=c(4,4,1,1), oma=c(3,1,2,0))
		
		plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Host", ylim=c(0, max(H1+H0)))
		lines(x, H1[i,], col=3)
		lines(x, H0[i,], col=2)
		
		plot(x, A0[i,]+A1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Arrested", ylim=c(0, max(A1+A0)))
		lines(x, A1[i,], col=3)
		lines(x, A0[i,], col=2)

		plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Parasite", ylim=c(0, max(P1+P0)))
		lines(x, P1[i,], col=3)
		lines(x, P0[i,], col=2)
		
		plot(x, E[i,], "l", lwd=1.2, ylab="Eggs", ylim=c(0, max(E)))
		plot(x, I[i,], "l", lwd=1.2, ylab="Infective", ylim=c(0, max(I)))

	}
}, img.name="ParasiteMigration", outdir="~/Documents/Migration Speed/Molnar/Animation/8Eqn", ani.dev = 'pdf', ani.type = 'pdf', ani.height = 8, ani.width = 7, ani.opts='controls', pdflatex = '/usr/texbin/pdflatex', caption=c(""), overwrite=TRUE, documentclass = paste("\\documentclass{article}", "\\usepackage[landscape,margin=0.3in]{geometry}", sep="\n")
)

