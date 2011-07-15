library(lsm)
require(Matrix)

# Load Davis' Southern Women dataset
data(davis)

# Fit the model via MCMC with 2 latent sets
K <- 2                    
dims <- list(T=nrow(davis), N=ncol(davis), K=K)

# Set priors (see documentation for more explanation)
priors <- list(theta=c(-2,10),z=c(1,1),epsilon=c(1,dims$N))

# Fit the model via Gibbs sampling
set.seed(123)
niter <- 100
fit <- gibbs(davis,dims,priors,niter=niter)

# Plot the progress of loglikelihood on the training data
pdf("davis-llks.pdf",width=5,height=5)
df <- data.frame(iter=1:niter,llk=fit$llks)
jjplot(llk ~ line() + iter, data=df, theme=jjplot.theme("bw"),xlab="iteration",ylab="training loglikelihood")
dev.off()

# Plot predictive distribution using last 20 iterations of MCMC
subs <- expand.grid(1:dims$T,1:dims$N)
keep <- (niter-20):niter
yhat <- post.predictive(list(fit),subs,keep)
ploty <- disp(t(davis))
plotyhat <- disp(t(yhat))

pdf("davis-y-yhat.pdf",width=8,height=5)
sidebysideplot <- grid.arrange(ploty,plotyhat,ncol=2)
dev.off()

# Plot W and Z \hadamard theta matrices
w <- fit$params$w
thetaz <- fit$params$z * fit$params$theta
rownames(w) <- rownames(davis)
rownames(thetaz) <- colnames(davis)
colnames(w) <- colnames(thetaz) <- paste("Set", 1:fit$dims$K)
plotW <- disp(w)
plotTheta <- disp(thetaz,limits=c(0,max(thetaz)))

pdf("davis-w-z.pdf",width=5,height=5)
grid.arrange(plotW,plotTheta,ncol=2)
dev.off()

# Fit multiple chains
nchains <- 3
fits <- list()
for (i in 1:nchains) 
  fits[[i]] <- gibbs(davis,dims,priors,niter=niter)

# Show posterior predictive estimates averaged over chains
keep <- (niter/2):niter
subs <- expand.grid(1:dims$T,1:dims$N)
yhat <- post.predictive(fits, subs, keep)
plot1 <- disp(t(davis))
plot2 <- disp(t(yhat))
pdf("davis-y-yhat-3chains.pdf",width=8,height=4)
sidebysideplot <- grid.arrange(plot1, plot2, ncol=2)
dev.off()

