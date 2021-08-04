#----------Implementation fo the exakt time evolution for 4 Qubits----------
#----------the Quantum Circuit is explained in my bachelor thesis----------

#The Qskit P gate
P<-  function(bit, angle, sign = +1){
  return(methods::new("sqgate",
                      bit=as.integer(bit),
                      M=array(as.complex(c(1., 0, 0, exp(sign*1i*angle))), dim=c(2,2)),
                      type="P"))
}
#QFFT for 4 qbits 
ft <- function(){
  function(wf){
    wf <- cqgate(bits=c(2,3), gate=Z(3)) * (SWAP(c(2,3)) * wf)
    wf <- plot.xequal(wf)
    
    #1,2
    wf <- P(1,2*pi*0.0/4.0)*wf
    wf <- CNOT(c(1,2))*wf
    wf <- cqgate(bits=c(2,1), gate=H(1)) * wf
    wf <- CNOT(c(1,2))*wf
    wf <- cqgate(bits=c(1,2), gate=Z(1)) * wf
    #3,4
    wf <- P(3,2*pi*0.0/4.0)*wf
    wf <- CNOT(c(3,4))*wf
    wf <- cqgate(bits=c(4,3), gate=H(3)) * wf
    wf <- CNOT(c(3,4))*wf
    wf <- cqgate(bits=c(3,4), gate=Z(3)) * wf
    
    wf <- cqgate(bits=c(2,3), gate=Z(3)) * (SWAP(c(2,3)) * wf)
    wf <- plot.xequal(wf)
    
    #1,2
    wf <- P(1,2*pi*1.0/4.0)*wf
    wf <- CNOT(c(1,2))*wf
    wf <- cqgate(bits=c(2,1), gate=H(1)) * wf
    wf <- CNOT(c(1,2))*wf
    wf <- cqgate(bits=c(1,2), gate=Z(1)) * wf
    #3,4
    wf <- P(3,2*pi*0.0/4.0)*wf
    wf <- CNOT(c(3,4))*wf
    wf <- cqgate(bits=c(4,3), gate=H(3)) * wf
    wf <- CNOT(c(3,4))*wf
    wf <- cqgate(bits=c(3,4), gate=Z(3)) * wf
    return(wf)
  }
}
#Inverse QFFT for 4 qbits
ftt <- function(){
  function(wf){
    #1,2
    wf <- cqgate(bits=c(1,2), gate=Z(2)) * wf
    wf <- CNOT(c(1,2))*wf
    wf <- cqgate(bits=c(2,1), gate=H(1)) * wf
    wf <- CNOT(c(1,2))*wf
    wf <- P(1,2*pi*1.0/4.0, sign = -1)*wf
    #3,4
    wf <- cqgate(bits=c(3,4), gate=Z(4)) * wf
    wf <- CNOT(c(3,4))*wf
    wf <- cqgate(bits=c(4,3), gate=H(3)) * wf
    wf <- CNOT(c(3,4))*wf
    wf <- P(3,2*pi*0.0/4.0, sign = -1)*wf
   
    wf <- cqgate(bits=c(2,3), gate=Z(3)) * (SWAP(c(2,3)) * wf)
    wf <- plot.xequal(wf)
    
    #1,2
    wf <- cqgate(bits=c(1,2), gate=Z(1)) * wf
    wf <- CNOT(c(1,2))*wf
    wf <- cqgate(bits=c(2,1), gate=H(1)) * wf
    wf <- CNOT(c(1,2))*wf
    wf <- P(1,2*pi*0.0/4.0, sign = -1)*wf
    #3,4
    wf <- cqgate(bits=c(3,4), gate=Z(3)) * wf
    wf <- CNOT(c(3,4))*wf
    wf <- cqgate(bits=c(4,3), gate=H(3)) * wf
    wf <- CNOT(c(3,4))*wf
    wf <- P(3,2*pi*0.0/4.0, sign = -1)*wf
    
    wf <- cqgate(bits=c(2,3), gate=Z(3)) * (SWAP(c(2,3)) * wf)
    wf <- plot.xequal(wf)
    
    return(wf)
  }
}
#BT for 4 qubits 
bt <- function(J,m_field){
  function(wf){
    wf <- X(2)*wf
    wf <- CNOT(c(2,1)) * wf
    wf <- cqgate(bits=c(1,2), gate=Rx(2,acos( (m_field - J*cos(2*pi/4))/ sqrt( (m_field - J*cos(2*pi/4))^2 + J^2*sin(2*pi/4)^2 ) )))*wf 
    wf <- CNOT(c(2,1)) * wf
    wf <- X(2)*wf
    wf <- X(4)*wf
    wf <- CNOT(c(4,3)) * wf
    wf <- cqgate(bits=c(3,4), gate=Rx(4,acos( (m_field - J*cos(0*pi/4))/ sqrt( (m_field - J*cos(0*pi/4))^2 + J^2*sin(0*pi/4)^2 ) )))*wf
    wf <- CNOT(c(4,3)) * wf
    wf <- X(4)*wf
    wf <- plot.xequal(wf)
    return(wf)
  }
}
#Inverse BT for 4 qbits
btt <- function(J,m_field){
  function(wf){
    wf <- X(2)*wf
    wf <- CNOT(c(2,1)) * wf
    wf <- cqgate(bits=c(1,2), gate=Rx(2,-acos( (m_field - J*cos(2*pi/4))/ sqrt( (m_field - J*cos(2*pi/4))^2 + J^2*sin(2*pi/4)^2 ) )))*wf 
    wf <- CNOT(c(2,1)) * wf
    wf <- X(2)*wf
    wf <- X(4)*wf
    wf <- CNOT(c(4,3)) * wf
    wf <- cqgate(bits=c(3,4), gate=Rx(4,-acos( (m_field - J*cos(0*pi/4))/ sqrt( (m_field - J*cos(0*pi/4))^2 + J^2*sin(0*pi/4)^2 ) )))*wf
    wf <- CNOT(c(4,3)) * wf
    wf <- X(4)*wf
    wf <- plot.xequal(wf)
    return(wf)
  }
}
#time evolution for 4 qubits
tstep<- function(j,m_field,t){
  function(wf){
    for(i in 1:nbits){
      k = c(2*pi/4,-2*pi/4,0,pi)
      wf <- Rz(i, - 2*sqrt( (m_field - J*cos(k[i]))^2 + J^2*sin(k[i])^2 )*t) * wf
    }
    return(wf)
  }
}
time_evolution_exakt <- function(initial_state,nbits,j,m_field,t,t_steps,rept){
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
      wf <- prep_state(initial_state, NA)
      wf <- ft()(wf)
      wf <- bt(J,m_field)(wf)
      wf <- tstep(J,m_field,t/t_steps*i)(wf)
      wf <- btt(J,m_field)(wf)
      wf <- ftt()(wf)
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
