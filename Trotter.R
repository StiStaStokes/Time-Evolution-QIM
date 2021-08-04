#----------Implementation fo the n-Qubits----------
#----------the Quantum Circuit is explained in my bachelor thesis----------


#The two time evolution operators for nbits
U1 <- function(nbits,J,dt,bc){
  function(wf){
    for (i in 0:(nbits/2-1)){
      wf <- CNOT(c((1+2*i),2*(i+1))) * ( Rx((1+2*i),-2*dt*J) * ( CNOT(c((1+2*i),2*(i+1))) * wf ) )
    }
    if(nbits > 2){
      for (i in 1:(nbits/2-0.5)){
        wf <- CNOT(c(2*i,2*i+1)) * (Rx(2*i,-2*dt*J) * ( CNOT(c(2*i,2*i+1)) * wf ) ) 
      }
    }
    wf <- plot.xequal(wf)
    if(bc != 0){
      wf <- CNOT(c(nbits,1)) * ( Rx((nbits),-bc*2*dt*J) * ( CNOT(c(nbits,1)) * wf ) )
      wf <- plot.xequal(wf)
    }
    return(wf)
  }
}
U2 <- function(nbits,m_field,dt){
  function(wf){
    for(i in 1:nbits){
      wf <- Rz(i,-2*dt*m_field) * wf 
    }
    return(wf)
  }
}
#One Trotter Step
trotter <- function(nbits,J,m_field,t,t_steps,bc){
  function(wf){
      wf <- U1(nbits,J,t/t_steps,bc)(wf)
      wf <- U2(nbits,m_field,t/t_steps)(wf)
    return(wf)
  }
}
#The whole time evoluation with trotter methode 
time_evolution_trotter <- function(initial_state,nbits,J,m_field,t,t_steps,rept,bc){
  mean_values <- array(dim=c(t_steps+1, nbits+1))
  sd_values <- array(dim=c(t_steps+1,nbits+1))
  wf <- prep_state(initial_state, NA)
  for(i in 1:nbits){ 
    mes <- measure(wf, i, rep=rept)$value[1:rept]
    mean_values[1,i] = mean(mes)
	sd_values[1,i]   = std.error(mes)
   }
  mean_values[1,nbits+1] = 0
  sd_values[1,nbits+1] = 0
  time <- seq(0+t/t_steps,t,t/t_steps)
  for(i in 1:t_steps){
    wf <- trotter(nbits,J,m_field,t,t_steps,bc)(wf)
    for(j in 1:nbits){ 
	  mes <- measure(wf, j, rep=rept)$value[1:rept]
	  mean_values[i+1,j] <- mean(mes) 
	  sd_values[i+1,j] <- std.error(mes)
	}
    mean_values[i+1,nbits+1] = (t/t_steps) * i
	sd_values[i+1,nbits+1] = (t/t_steps) * i
  }
  return(list("WF"=wf,"results"=mean_values, "sd"=sd_values))
}

#Compare the solution of the time evolution for different Trotter Steps 
trott_step_compare <- function(initial_state,nbits,j,m_field,t,t_steps,rept,bc,bit,max_steps){
  steps <- 1:max_steps
  t_values <- array(dim=c(length(steps), 2))
  for(i in 1:length(steps)){
    wf <- prep_state(initial_state, NA)
    for(j in 1:steps[i]){
      wf <- trotter(nbits,cconst,m_field,t,steps[i],bc)(wf)
    }
    t_values[i,1] = mean(measure(wf, bit, rep=rept)$value[1:rept])
    t_values[i,2] = (t*cconst)/steps[i]
  }
  return(t_values)
}