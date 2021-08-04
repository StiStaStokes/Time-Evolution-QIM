#needed libaries
library(datasets)
library(qsimulatR)
library(tidyverse)
library(latex2exp)
library(plotrix)

#Input Quantum state
prep_state <- function(state, noise){
  wf <- qstate(nbits = length(state))
  for(i in 1:length(state)){
    if(state[i] == 1){ wf <- X(i)*wf }
  }
  wf <- plot.xequal(wf)
  if(!is.na(noise)){ wf@noise = genNoise(nbits=length(state), p=noise, error="X", sigma=0.05)}
  return(wf)
} 
#comparission of the two methods
exakt_trotter_compare <- function(initial_state,nbits,j,m_field,t,t_steps,rept,bc,bit,max_steps){
  steps <- 1:max_steps
  compared_values <- array(dim=c(max_steps,2))
  compared_sd <- array(dim=c(max_steps,1))
  for(i in 1:length(steps)){
    exakt_values <- time_evolution(initial_state,nbits,j,m_field,t,t_steps,rept)$results
    trotter_values <- time_evolution_trotter(initial_state,nbits,J,m_field,t,i,rept,bc)c)$results
    compare_calc<- abs(Trotter_values[2:(i+1),1:4]-exakt_values[2:(i+1),1:4])
	compared_sd[i] <- sd(compare_calc[1:4])
    compared_values[i,1]=mean(compare_calc[1:4],na.rm = TRUE)
    compared_values[i,2]=i
  }
  return(data.frame(Freq.A = compared_values[1:max_steps,1],Freq.B = compared_values[1:max_steps,2], Freq.C = compared_sd[1:max_steps]))
}
  
#----------Variables----------
initial_state <- c(0,0,0,1)
nbits = length(initial_state)
cconst <- 0.8
m_field <- 0.9
t <- 15
t_steps <- 120
rept <- 5000

#----------Different Trotter Steps----------
max_Steps <- 50	#Steps vary from 1 to this
bc <- 1
bit <- 1
t_vals <- trott_step_compare(initial_state,nbits,j,m_field,t,t_steps,rept,bc,bit,max_steps)
#Convert array to data frame for ggplot
d2 <- data.frame(y=double(), x=double(), stringsAsFactors=FALSE)
d2 <- as.data.frame.table(t_vals)
d1 <- reshape(d2, timevar = 'Var2', idvar ="Var1", direction = "wide")
df <- d1[-c(1)]
#For the axis of the log-plot
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
name <- paste('$<n_',i,'>(t)$', sep ='')
#the plot itself
ggplot(df, aes(x = Freq.B, y = Freq.A)) + theme_bw() + theme(axis.title=element_text(size=10,face="bold")) +
  ylab(TeX(name)) +
  xlab(TeX('$\\delta t$')) +
  geom_point(shape=3,size=1,color="#FF4500") + 
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks) + annotation_logticks(sides="b") +
  ggtitle(paste("qbit", bit ))

#----------Different boundary conditions with Trotter----------
bc <- c(1,0,-1) #PBC, OBC, ABC
color <- c("red","blue","black")
mean_values <- array(dim=c(t_steps+1, length(bc)*(nbits+1)))
wf <- list()
#calc. the time evolution with Trotter for all bc in the list
for (i in 1:length(bc)){
  results <- time_evolution_trotter(initial_state,nbits,J,m_field,t,t_steps,rept,bc[i])
  mean_values[1:(t_steps+1), (1+((nbits+1)*(i-1))):((nbits+1)*i)] <- results$results
  wf[i] <- results$WF
}
#print the time evolution
par(mfrow=c(nbits%/%2,nbits - nbits%/%2))
for(i in 1:nbits){
    name <- paste('$<n_',i,'>(t)$', sep ='')
    plot( mean_values[1:(t_steps+1),nbits+1],mean_values[1:(t_steps+1),i],
          xlim= c(0,t+1), ylim=c(0,1.1), axes=FALSE, main = paste("bit", i ), xlab = "t", ylab = TeX(name), pch = 20,
          panel.first=grid(), type="n")
    for(j in 1:length(bc)){
      points(mean_values[1:(t_steps+1),((nbits+1)*j)],mean_values[1:(t_steps+1),i+(5*(j-1))], col=color[j], pch=20)
    }
    axis(1,at = seq(0,t+1,2.5))
    axis(2,at = seq(0,1.1,0.5))
    box()
}

#----------Different J with Trotter----------
J <- c(0.02,0.1,0.4,0.7) 
bit <- 4 #choose bit to be plotted
bc <- 1
mean_values<- array(dim=c(t_steps+1, length(J)*(nbits+1)))
wf <- list()
#calc. and plot
par(mfrow=c(length(J)%/%2,length(J) - length(J)%/%2))	
for (i in 1:length(J)){
  name <- paste('$<n_',bit,'>(t)$', sep ='')
  results <- time_evolution_trotter(initial_state,nbits,J[i],m_field,t,t_steps,rept,bc)
  mean_values[1:(t_steps+1), (1+((nbits+1)*(i-1))):((nbits+1)*i)] <- results$results
  wf[i] <- results$WF
  plot( mean_values[1:(t_steps+1),((nbits+1)*i)],mean_values[1:(t_steps+1),bit+(5*(i-1))],
        xlim= c(0,t+1), ylim=c(0,1.1), axes=FALSE, main = paste("J=", J[i] ), xlab = "t", ylab = TeX(name), pch = 20,
        panel.first=grid(), type="n")
  points(mean_values[1:(t_steps+1),((nbits+1)*i)],mean_values[1:(t_steps+1),bit+(5*(i-1))], col="red", pch=20)
  axis(1,at = seq(0,t+1,2.5))
  axis(2,at = seq(0,1.1,0.5))
  box()
}

#----------Simulation based on the analytical Solution----------
#calc. of the time evoluation
results <- time_evolution_exakt(initial_state,nbits,j,m_field,t,t_steps,rept)
mean_values <- results$results
wf <- results$WF
sd_values <- results$sd
#plot time evolution
par(mfrow=c(nbits%/%2,nbits - nbits%/%2))
for(i in 1:nbits){
  name <- paste('$<n_',i,'>(t)$', sep ='')
  plot( mean_values[1:(t_steps+1),nbits+1],mean_values[1:(t_steps+1),i],
        xlim= c(0,t+1), ylim=c(0,1.1), axes=FALSE, main = paste("bit", i ), xlab = "t", ylab = TeX(name), pch = 20,
        panel.first=grid(), type="n")
  points(mean_values[1:(t_steps+1),nbits+1],mean_values[1:(t_steps+1),i], col="red", pch=20)
  axis(1,at = seq(0,t+1,2.5))
  axis(2,at = seq(0,1.1,0.5))
  box()
}

#----------Comparission of Trotter and analytical method----------
bc <- 1 #Bc for Trotter have to be choosen depending on the amoud of |1> qubits in the input state 
compare_values <- exakt_trotter_compare(initial_state,nbits,j,m_field,t,t_steps,rept,bc,bit,max_steps)
min <- min(compare_values$Freq.A - compare_values$Freq.C)
max <- max(compare_values$Freq.A + compare_values$Freq.C)

#the plot itself
name <- paste("$J t =",cconst*t,"$, $J=", cconst, "$, $t=",t,"$", sep ='')
ggplot(compare_values, aes(x = Freq.B, y = Freq.A)) + theme_bw() + theme(axis.title=element_text(size=10,face="bold")) +
  ylab(TeX('$\\Delta <n>$')) + xlab(TeX('$N_{Trott}')) + scale_x_continuous(breaks=seq(0,Nstep,5)) + 
  scale_y_continuous(breaks = sort(c(seq(0, 1, length.out=5), 0.05)),limits=c(min,max),
                     minor_breaks = sort(c(seq(0,1,length.out = 9)))) +
  geom_vline(xintercept = 4.5,alpha = 0.8) + geom_hline(yintercept = 0.05,alpha = 0.3) +
  geom_point(shape=19,size=1,color="#FF4500") + ggtitle(TeX(name)) + 
  geom_errorbar(aes(ymin=Freq.A-Freq.C, ymax=Freq.A+Freq.C), width=.3,
                position=position_dodge(.9), color="#FF4500")