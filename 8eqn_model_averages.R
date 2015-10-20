#Constant background source of parasites
b_H<-0.001			    # Birth rate of hosts
mu_H<-0.001				# Natural death rate of hosts
alpha<-0.002				# Per-parasite rate of parasite-induced host mortality

lambda_1<-0.1
lambda_0<-0.08

beta<-0.001				# Instantaneous rate of uptake
mu_A<-0.2				# Non-reproductive (A) parasite death
tau_A<-10				# Developmental time of parasites to reproductive stage in host

mu_P<-0.10				# Reproductive (P) Parasite death

rho<-10					# Birth (shedding?) rate of larvae
mu_E<-0.1				# Death of eggs
tau_E<-20				# Development time of eggs

mu_I<-0.3				# Death of infective larvae

c.hat<-1				# Migration speed in the absence of parasites
k<-0.01					# Parasite aggregation


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
	
	# Hosts
	H1[i,] <- H1[i-1, x.i] + (b_H-mu_H)*H1[i-1, x.i] - alpha*H1[i-1, x.i]*P1[i-1, x.i] - lambda_1*H1[i-1, x.i] + lambda_0*H0[i-1,]
	H0[i,] <- H0[i-1,] + (b_H-mu_H)*H0[i-1,] - alpha*H0[i-1,]*P0[i-1,] - lambda_0*H0[i-1,] + lambda_1*H1[i-1, x.i]
	
	# Arrested (non-reproductive) parasites
	A1[i,] <- A1[i-1,x.i] + beta*I[i-1,x.i] - (mu_A + b_H + 1/tau_A)*A1[i-1,x.i] + lambda_1*(A0[i-1,]-A1[i-1, x.i])
	A0[i,] <- A0[i-1,] + beta*I[i-1,] - (mu_A + b_H + 1/tau_A)*A0[i-1,] + lambda_0*(A1[i-1,x.i]-A0[i-1,])
	
	# Reproductive parasites
	P1[i,] <- P1[i-1,x.i] + (1/tau_A)*A1[i-1, x.i] - (mu_P+b_H)*P1[i-1, x.i] - alpha*(P1[i-1, x.i]+P1[i-1,x.i]^2/k) + lambda_1*(P0[i-1,]-P1[i-1,x.i])		
	P0[i,] <- P0[i-1,] + (1/tau_A)*A0[i-1,] - (mu_P+b_H)*P0[i-1,] - alpha*(P0[i-1,]+P0[i-1,]^2/k) + lambda_0*(P1[i-1,x.i]-P0[i-1,])
			
	# Eggs
	E[i,] <- E[i-1,]+ rho*(H0[i-1,]*P0[i-1,] + H1[i-1,x.i]*P1[i-1,x.i]) - mu_E*E[i-1,] - (1/tau_E)*E[i-1,]
	
	# Infective
	I[i,] <- I[i-1,] + (1/tau_E)*E[i-1,] - mu_I*I[i-1,] - beta*I[i-1,]
	}

cat("\n H1:", range(H1), " H0:", range(H0), "\n A1: ", range(A1), " A0: ", range(A0), "\n P1: ", range(P1), " P0: ", range(P0), "\n Eggs: ", range(E), " I: ", range(I))


# j<-min(which(is.na(A1)==TRUE, arr.ind=TRUE)[,1])
# j<-min(which(A1<0, arr.ind=TRUE)[,1])
# for(i in c(init, init+1, j-50,j-40, j-30, j-20, j-10, j-5, j-3, j-2, j-1, j)){
# par(mfrow=c(3,2),mar=c(4,5,2,2))
# plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.2, main=i)
# lines(x, H1[i,], col=3)
# lines(x, H0[i,], col=2)

# plot(x, A0[i,]+A1[i,], "l", col=grey(0.8), lwd=1.2)
# lines(x, A1[i,], col=3)
# lines(x, A0[i,], col=2)

# plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.2)
# lines(x, P1[i,], col=3)
# lines(x, P0[i,], col=2)

# plot(x, E[i,], "l", col=4)
# plot(x, I[i,], "l", col=4)

# }
require(animation)

saveLatex({
	ani.options(interval=0.06)
	for(i in seq(init,n.t,20)){
		layout(matrix(c(1,2,3,4,1,2,3,5), nrow=4, ncol=2))
		par(mar=c(4,4,1,1), oma=c(3,1,2,0))
		
		plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Hosts", ylim=c(0, max(H1+H0)))
		lines(x, H1[i,], col=3)
		lines(x, H0[i,], col=2)
		
		plot(x, A0[i,]+A1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Arrested parasites per host", ylim=c(0, max(A1+A0)))
		lines(x, A1[i,], col=3)
		lines(x, A0[i,], col=2)

		plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Reproductive parasites per host", ylim=c(0, max(P1+P0)))
		lines(x, P1[i,], col=3)
		lines(x, P0[i,], col=2)
		
		plot(x, E[i,], "l", lwd=1.2, ylab="Eggs", ylim=c(0, max(E)))
		plot(x, I[i,], "l", lwd=1.2, ylab="Infective", ylim=c(0, max(I)))

	}
}, img.name="ParasiteMigration", outdir="~/Documents/Migration Speed/Molnar/Animation/8EqnAvg", ani.dev = 'pdf', ani.type = 'pdf', ani.height = 8, ani.width = 7, ani.opts='controls', pdflatex = '/usr/texbin/pdflatex', caption=c(""), overwrite=TRUE, documentclass = paste("\\documentclass{article}", "\\usepackage[landscape,margin=0.3in]{geometry}", sep="\n")
)

########################################################################################################
