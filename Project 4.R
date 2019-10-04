library(data.table)
library(viridis)
library(reshape2)
library(geoR)
library(ggplot2)
library(MASS)
library(cowplot)
library(fields)
library(gtable)
library(ggpubr)
library(gridExtra)
library(akima)
library(party)

n = 100
B = 200
t_j = as.vector(as.matrix(read.table("https://folk.ntnu.no/joeid/emnemodul/traveltimedata.txt")))

set.seed(2606)
i = 1:100

sigma = 0.05
tau = 0.1
eta = 0.1


h = distances::distances(i)[1:n,1:n]
cov_mat = list(sigma^2*(h*eta+1)*exp(-eta*h))
mu = list(0.5 - (1:100)*0.001)
x = list(mvrnorm(B,mu[[1]],cov_mat[[1]]))

t_asim <- function(x, j){
  return(sum(x[1:(50+j)])/cos(atan(40/(50+j))))
}

sigma_t_j = rep(NA, 50)
sigma_x_t_j = matrix(NA, nrow = 50, ncol = n)

t_t = matrix(NA, ncol = 50, nrow = B)
K = matrix(NA, ncol = n, nrow = 50)
epsilon_j = rnorm(50, 0, tau)

rev = F

if(rev){iter = 50:1}else{iter = 1:50}

for (j in iter){
  if (rev){k = 51-j}else{k=j}
  print(k)
  for (b in 1:B){
    t_t[b,j] = t_asim(x[[k]][b,], j)
  }
  t_t[,j] = t_t[,j] + epsilon_j[j]
  sigma_t_j[j] = 1/B*sum((t_t[,j] - mean(t_t[,j]))^2)
  sigma_x_t_j[j,] = rep(0, n)
  for (b in 1:B){
    sigma_x_t_j[j,] = sigma_x_t_j[j,] + (x[[k]][b,] - mean(x[[k]][b,]))*(t_t[b,j]-mean(t_t[,j]))
  }
  sigma_x_t_j[j,] = sigma_x_t_j[j,]/B
  K[j, ] = sigma_x_t_j[j,]/sigma_t_j[j]
  x[[k+1]] = matrix(NA, ncol = 100, nrow = B)
  for (b in 1:B){
    x[[k+1]][b, ] = x[[k]][b, ] + K[j, ]*(t_j[j] - t_t[b, j])
  }
}

cmat = x[[51]]
rownames(cmat) = paste("trial", seq(B), sep="")
colnames(cmat) = paste("time", seq(n), sep="")

dat = as.data.frame(cmat)
dat$trial = rownames(dat)
mdat = melt(dat, id.vars="trial")
mdat$time = as.numeric(gsub("time", "", mdat$variable))


p = ggplot(mdat, aes(x=time, y=value, group=trial, col = trial)) +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  geom_line()+ theme(legend.position = "none")+ xlab("i") + ylab("x")
p


h = distances::distances(i)[1:n,1:n]
cov_mat = list(sigma^2*(h*eta+1)*exp(-eta*h))
mu = list(0.5 - (1:100)*0.001)


for (j in 1:50){
  g_j = c(rep(1,50+j), rep(0,50-j))/cos(atan(40/(50+j)))
  K = cov_mat[[j]]%*%g_j%*%solve(g_j%*%cov_mat[[j]]%*%g_j+tau^2)
  mu[[j+1]] = mu[[j]] + K%*%(t_j[j]-g_j%*%mu[[j]])
  cov_mat[[j+1]] = cov_mat[[j]] - K%*%g_j%*%cov_mat[[j]]
}

mu_est = mu[[51]]
cov_est = cov_mat[[51]]

lb = mu_est - qnorm(0.9)*sqrt(diag(cov_est))
ub = mu_est + qnorm(0.9)*sqrt(diag(cov_est))

df_zero = data.frame(x = seq(1,100), pred = mu_est, lb = lb, ub = ub)

ggplot() +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  geom_line(data = df_zero, aes(x = x, y = pred), size = 1, col = 'red') + 
  geom_line(data = df_zero, aes(x = x, y = lb), size = 1, col = 'blue', linetype = "dashed") + 
  geom_line(data = df_zero, aes(x = x, y = ub), size = 1, col = 'blue', linetype = "dashed") + xlab("i") + ylab("x")
