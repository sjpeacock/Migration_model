#Constant background source of parasites
b_h<-0.05			    # Birth rate of hosts
mu_h<-0.045				# Natural death rate of hosts
alpha_h<-0.05			# Per-parasite rate of parasite-induced host mortality

rho<-0.8				# Instantaneous rate of uptake
mu_p<-0.15				# Parasite death

lambda<-2				# Birth (shedding?) rate of larvae
mu_L<-0.01					# Death of larvae

c.hat<-1				# Migration speed in the absence of parasites
k<-0.5					# Parasite aggregation

L<-1

lambda_1<-0.06
lambda_0<-0.04
# Setup matrices
x<-c(0:1200)
t<-c(0:1000)
n.x<-length(x)
n.t<-length(t)

H1<-matrix(nrow=n.t, ncol=n.x)
H0<-matrix(nrow=n.t, ncol=n.x)
P1<-matrix(nrow=n.t, ncol=n.x)
P0<-matrix(nrow=n.t, ncol=n.x)

# Initial conditions
init<-1 #Number of timesteps to initiate
sig<-12; mu<-50
H1[1:init,]<-matrix(rep(5000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)
H0[1:init,]<-matrix(rep(5000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)

P1[1:init,]<-matrix(rep(0.1*H1[1,], init),nrow=init, ncol=n.x, byrow=TRUE) #1% prevalence to start
P0[1:init,]<-matrix(rep(0.1*H0[1,], init),nrow=init, ncol=n.x, byrow=TRUE) #1% prevalence to start


# Simulation
for(i in 2:n.t){
	C<-round(c.hat)
	H1.last<-H1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	H0.last<-H0[i-1,]
	P1.last<-P1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	P0.last<-P0[i-1,]
	
	H1[i,]<-H1.last+(b_h-mu_h)*H1.last - alpha_h*P1.last - lambda_1*H1.last + lambda_0*H0.last
	H0[i,]<-H0.last+(b_h-mu_h)*H0.last - alpha_h*P0.last - lambda_0*H0.last + lambda_1*H1.last
	
	agg.term1<-(P1.last/H1.last+(P1.last^2/H1.last^2*(k+1)/k))
	agg.term1[is.na(agg.term1)|agg.term1==Inf]<-0
	agg.term0<-(P0.last/H0.last+(P0.last^2/H0.last^2*(k+1)/k))
	agg.term0[is.na(agg.term0)|agg.term0==Inf]<-0
	
	P1[i,]<- P1.last + rho*L*H1.last - (mu_p+mu_h)*P1.last - alpha_h*H1.last*agg.term1
			- lambda_1*P1.last + lambda_0*P0.last
			
	P0[i,]<- P0.last + rho*L*H0.last - (mu_p+mu_h)*P0.last - alpha_h*H0.last*agg.term0
			- lambda_0*P0.last + lambda_1*P1.last
	}

cat("\n H1:", range(H1), " H0:", range(H0), "\n P1", range(P1), " P0", range(P0), "\n Larvae", range(L))


par(mfrow=c(2,1),mar=c(4,5,2,2))

for(i in c(1,5,10,100, 150, 200, 500, 800)){
plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.2, ylim=c(0, max(H0+H1)))
lines(x, H1[i,], col=3)
lines(x, H0[i,], col=2)

plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.2, ylim=c(0, max(P0+P1)))
lines(x, P1[i,], col=3)
lines(x, P0[i,], col=2)

}

#####################################################################
 # Lambdas depend on the number of parasites (not really working)
#####################################################################
 
l0<-0.02
l1<-0.02
p<-seq(0,10,0.1)

plot(p, l0*exp(-0.3*p), "l", ylim=c(0,0.1))
lines(p, l1*(1-exp(-0.3*p)), lty=2)

# Setup matrices
x<-c(0:1200)
t<-c(0:1000)
n.x<-length(x)
n.t<-length(t)

H1<-matrix(nrow=n.t, ncol=n.x)
H0<-matrix(nrow=n.t, ncol=n.x)
P1<-matrix(nrow=n.t, ncol=n.x)
P0<-matrix(nrow=n.t, ncol=n.x)

# Initial conditions
init<-1 #Number of timesteps to initiate
sig<-12; mu<-50
H1[1:init,]<-matrix(rep(5000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)
H0[1:init,]<-0
P1[1:init,]<-matrix(rep(0.1*H1[1,], init),nrow=init, ncol=n.x, byrow=TRUE) #1% prevalence to start
P0[1:init,]<-0

# Simulation
for(i in 2:n.t){
	C<-round(c.hat)
	H1.last<-H1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	H0.last<-H0[i-1,]
	P1.last<-P1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	P0.last<-P0[i-1,]
	
	p0<-P0.last/H0.last; p0[is.na(p0)]<-0
	l0<-0.02*exp(-0.3*p0)
	
	p1<-P1.last/H1.last; p1[is.na(p1)]<-0
	l1<-0.02*(1-exp(-0.3*p1))
	
	H1[i,]<-H1.last+(b_h-mu_h)*H1.last - alpha_h*P1.last - l1*H1.last + l0*H0.last
	H0[i,]<-H0.last+(b_h-mu_h)*H0.last - alpha_h*P0.last - l0*H0.last + l1*H1.last
	
	agg.term1<-(P1.last/H1.last+(P1.last^2/H1.last^2*(k+1)/k))
	agg.term1[is.na(agg.term1)|agg.term1==Inf]<-0
	agg.term0<-(P0.last/H0.last+(P0.last^2/H0.last^2*(k+1)/k))
	agg.term0[is.na(agg.term0)|agg.term0==Inf]<-0
	
	P1[i,]<- P1.last + rho*L*H1.last - (mu_p+mu_h)*P1.last - alpha_h*H1.last*agg.term1
			- l1*P1.last + l0*P0.last
			
	P0[i,]<- P0.last + rho*L*H0.last - (mu_p+mu_h)*P0.last - alpha_h*H0.last*agg.term0
			- l0*P0.last + l1*P1.last
	}

cat("\n H1:", range(H1), " H0:", range(H0), "\n P1", range(P1), " P0", range(P0), "\n Larvae", range(L))

par(mfrow=c(2,1),mar=c(4,5,2,2))

for(i in c(1,5,10,15,25,35,50,80,100,150, 500, 800)){
plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.2)#, ylim=c(0, max(H0+H1)))
lines(x, H1[i,], col=3)
lines(x, H0[i,], col=2)
mtext(side=3, i)
plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.2)#, ylim=c(0, max(P0+P1)))
lines(x, P1[i,], col=3)
lines(x, P0[i,], col=2)

}


#####################################################################
 # Lambdas depend on the time
#####################################################################
b_h<-0.05			    # Birth rate of hosts
mu_h<-0.045				# Natural death rate of hosts
alpha_h<-0.05			# Per-parasite rate of parasite-induced host mortality

rho<-0.8				# Instantaneous rate of uptake
mu_p<-0.15				# Parasite death

lambda<-2				# Birth (shedding?) rate of larvae
mu_L<-0.01					# Death of larvae

c.hat<-1				# Migration speed in the absence of parasites
k<-0.5					# Parasite aggregation

# Setup matrices
x<-c(0:1200)
t<-c(0:1600)
n.x<-length(x)
n.t<-length(t)

l<-0.02
L0<-2*l*(1+sin(t/1000*pi)) #prob of starting to move: start small, get big, end small
L1<-2*l*(1+sin((t+1000)/1000*pi)) #prob of stopping, start big, get small, end big
plot(t, L0, "l", ylim=c(0,0.1))
lines(t, L1, lty=2)
legend("topleft", ncol=2, lty=c(1,2), c(expression(paste(lambda[0])), expression(paste(lambda[1]))), bty="n")

H1<-matrix(nrow=n.t, ncol=n.x)
H0<-matrix(nrow=n.t, ncol=n.x)
P1<-matrix(nrow=n.t, ncol=n.x)
P0<-matrix(nrow=n.t, ncol=n.x)

# Initial conditions
init<-1 #Number of timesteps to initiate
sig<-40; mu<-100
H0[1:init,]<-matrix(rep(5000/sqrt(2*pi*(sig^2))*exp(-(x-mu)^2/(2*sig^2)), init), nrow=init, ncol=n.x, byrow=TRUE)
H1[1:init,]<-0
P0[1:init,]<-matrix(rep(0.1*H0[1,], init),nrow=init, ncol=n.x, byrow=TRUE) #1% prevalence to start
P1[1:init,]<-0

# Simulation
for(i in 2:n.t){
	C<-round(c.hat)
	H1.last<-H1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	H0.last<-H0[i-1,]
	P1.last<-P1[i-1, c((n.x+1-C):n.x, 1:(n.x-C))]
	P0.last<-P0[i-1,]
	
	l0<-L0[i]
	l1<-L1[i]
	
	H1[i,]<-H1.last+(b_h-mu_h)*H1.last - alpha_h*P1.last - l1*H1.last + l0*H0.last
	H0[i,]<-H0.last+(b_h-mu_h)*H0.last - alpha_h*P0.last - l0*H0.last + l1*H1.last
	
	agg.term1<-(P1.last/H1.last+(P1.last^2/H1.last^2*(k+1)/k))
	agg.term1[is.na(agg.term1)|agg.term1==Inf]<-0
	agg.term0<-(P0.last/H0.last+(P0.last^2/H0.last^2*(k+1)/k))
	agg.term0[is.na(agg.term0)|agg.term0==Inf]<-0
	
	P1[i,]<- P1.last + rho*L*H1.last - (mu_p+mu_h)*P1.last - alpha_h*H1.last*agg.term1
			- l1*P1.last + l0*P0.last
			
	P0[i,]<- P0.last + rho*L*H0.last - (mu_p+mu_h)*P0.last - alpha_h*H0.last*agg.term0
			- l0*P0.last + l1*P1.last
	}

cat("\n H1:", range(H1), " H0:", range(H0), "\n P1", range(P1), " P0", range(P0), "\n Larvae", range(L))

par(mfrow=c(2,1),mar=c(4,5,2,2))

for(i in c(1,10,80,100,150, 500, 600,800,900,1000,1200,1400,1600)){
plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.2)#, ylim=c(0, max(H0+H1)))
lines(x, H1[i,], col=3)
lines(x, H0[i,], col=2)
mtext(side=3, i)
plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.2)#, ylim=c(0, max(P0+P1)))
lines(x, P1[i,], col=3)
lines(x, P0[i,], col=2)

}

#------------------------------------------------------------------------------------------
require(animation)

saveLatex({
	par(mfrow=c(3,1), mar=c(4,4,1,1), oma=c(3,1,2,0))
	ani.options(interval=0.06)
	for(i in seq(6,n.t,20)){
		plot(t, L0, "l", col=3, ylim=c(0,0.1), ylab="Switching rates", xlab="time")
		lines(t, L1, col=2)
		points(c(i,i), c(L0[i], L1[i]))
		mtext(side=3, paste("time =", i))
		
		plot(x, H0[i,]+H1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Host", ylim=c(0, max(H1+H0)))
		lines(x, H1[i,], col=3)
		lines(x, H0[i,], col=2)
		
		plot(x, P0[i,]+P1[i,], "l", col=grey(0.8), lwd=1.5, ylab="Parasite", ylim=c(0, max(P1+P0)))
		lines(x, P1[i,], col=3)
		lines(x, P0[i,], col=2)
	}
}, img.name="ParasiteMigration", outdir="~/Documents/Migration Speed/Molnar/Animation/SpeedTime", ani.dev = 'pdf', ani.type = 'pdf', ani.height = 6, ani.width = 8, ani.opts='controls,width=8in', pdflatex = '/usr/texbin/pdflatex', caption=c(""), overwrite=TRUE, documentclass = paste("\\documentclass{article}", "\\usepackage[landscape,margin=0.3in]{geometry}", sep="\n")
)
