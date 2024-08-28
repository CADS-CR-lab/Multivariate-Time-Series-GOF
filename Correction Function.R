# Multivariate Time Series, Type 1 error Correction

# The function is built to calculate the correct Type 1 error rates and
# returns one critical value for any combination of (num_of_ts, sample_size, lag)
# Parameters include: number of time series, sample size, and lag
# Number of time series input: 2, 3, 5, 10
# Sample size input: any sample size from 40 to 300
# Lag input: any lag from 2 to (sample_size - 1)
                            
MTS_T1_Correction <- function(num_of_ts, sample_size, lag){ 
  # models for analyzing 2 time series
  if (num_of_ts==2){
    # When sample size is less than or equal to 50
    if (sample_size<=50){
      k2_mod1 <- c(sample_size^1.8,lag^1.3,sample_size^(1.8)*lag^(1.3),sample_size^(2*1.8)*lag^(2*1.3),sample_size^(2*1.8),sample_size^(3*1.8)*lag^(2*1.3),sample_size^(3*1.8)*lag^(3*1.3),lag^(4*1.3),lag^(5*1.3))
      k2_coef1 <- c(-3.708346e-05,2.509301e-03,-1.372297e-06,-8.709501e-11,2.461332e-08,7.446143e-14,-8.681287e-17,1.071161e-09,-2.212193e-12)
      crit_val <- sum(k2_mod1*k2_coef1) + 0.05
      return (crit_val)
    # When sample size is between 51 and 70
    }else if (50<sample_size & sample_size<=70){
      k2_mod2 <- c(sample_size^2,lag^1.1,sample_size^(2)*lag^(1.1),sample_size^(2*2)*lag^(2*1.1),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.1),sample_size^(3*2)*lag^(3*1.1),lag^(4*1.1),lag^(5*1.1))
      k2_coef2 <- c(-8.985535e-06,4.351975e-03,-7.116296e-07,-6.841265e-12,1.503328e-09,1.493836e-15,-3.374209e-18,1.279253e-09,9.205486e-12)
      crit_val <- sum(k2_mod2*k2_coef2) + 0.05
      return (crit_val)
    # When sample size is between 71 and 90
    }else if (70<sample_size & sample_size<=90){
      k2_mod3 <- c(sample_size^2,lag^1.1,sample_size^(2)*lag^(1.1),sample_size^(2*2)*lag^(2*1.1),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.1),sample_size^(3*2)*lag^(3*1.1),lag^(4*1.1),lag^(5*1.1))
      k2_coef3 <- c(-4.630759e-06,3.380715e-03,-3.276978e-07,-1.292280e-12,4.772458e-10,1.702767e-16,-3.122203e-19,3.515224e-10,2.687904e-12)
      crit_val <- sum(k2_mod3*k2_coef3) + 0.05
      return (crit_val)
    # When sample size is between 91 and 120
    }else if (90<sample_size & sample_size<=120){
      k2_mod4 <- c(sample_size^2,lag^2,sample_size^(2)*lag^(2),sample_size^(2*2)*lag^(2*2),sample_size^(2*2),sample_size^(3*2)*lag^(2*2),sample_size^(3*2)*lag^(3*2),lag^(4*2),lag^(5*2))
      k2_coef4 <- c(2.818366e-06,-1.329575e-05,2.709333e-09,-3.982286e-17,-2.048533e-10,1.133300e-21,3.891565e-26,2.761994e-17,-1.399962e-21)
      crit_val <- sum(k2_mod4*k2_coef4) + 0.05
      return (crit_val)
    # When sample size is between 121 and 200
    }else if (120<sample_size & sample_size<=200){
      #c(n^s,m^p,n^(s)*m^(p),n^(2*s)*m^(2*p),n^(2*s),n^(3*s)*m^(2*p),n^(3*s)*m^(3*p),m^(4*p),m^(5*p))
      k2_mod5 <- c(sample_size^0.3,lag^0.9,sample_size^(0.3)*lag^(0.9),sample_size^(2*0.3)*lag^(2*0.9),sample_size^(2*0.3),sample_size^(3*0.3)*lag^(2*0.9),sample_size^(3*0.3)*lag^(3*0.9),lag^(4*0.9),lag^(5*0.9))
      k2_coef5 <- c(2.662781e-03,8.581804e-03,-1.806328e-03,-1.106108e-05,-7.297226e-04,2.905030e-06,-1.190938e-08,7.340158e-09,-3.837861e-12)
      crit_val <- sum(k2_mod5*k2_coef5) + 0.05
      return (crit_val)
    # When sample size is larger than 200
    }else if (sample_size>200){
      k2_mod6 <- c(sample_size^1.5,lag^1.5,sample_size^(1.5)*lag^(1.5),sample_size^(2*1.5)*lag^(2*1.5),sample_size^(2*1.5),sample_size^(3*1.5)*lag^(2*1.5),sample_size^(3*1.5)*lag^(3*1.5),lag^(4*1.5),lag^(5*1.5))
      k2_coef6 <- c(2.818366e-06,-1.329575e-05,2.709333e-09,-3.982286e-17,-2.048533e-10,1.133300e-21,3.891565e-26,2.761994e-17,-1.399962e-21)
      crit_val <- sum(k2_mod6*k2_coef6) + 0.05
      return (crit_val)
    }
  # models for analyzing 3 time series
  }else if(num_of_ts==3){
    # When sample size is less than or equal to 50
    if (sample_size<=50){
      # s=1.3, p=0.7, n<=50
      k3_mod1 <- c(sample_size^1.3,lag^0.7,sample_size^(1.3)*lag^(0.7),sample_size^(2*1.3)*lag^(2*0.7),sample_size^(2*1.3),sample_size^(3*1.3)*lag^(2*0.7),sample_size^(3*1.3)*lag^(3*0.7),lag^(4*0.7),lag^(5*0.7))
      k3_coef1 <- c(-9.189031e-04,6.391232e-02,-5.035098e-04,-1.626540e-07,6.186786e-06,2.127624e-09,-6.148985e-11,-1.233283e-05,9.143425e-07)
      crit_val <- sum(k3_mod1*k3_coef1) + 0.05
      return (crit_val)
    # When sample size is between 51 and 70
    }else if (50<sample_size & sample_size<=70){
      k3_mod2 <- c(sample_size^2,lag^1.1,sample_size^(2)*lag^(1.1),sample_size^(2*2)*lag^(2*1.1),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.1),sample_size^(3*2)*lag^(3*1.1),lag^(4*1.1),lag^(5*1.1))
      k3_coef2 <- c(-1.093557e-05,5.459875e-03,-9.005444e-07,-9.375999e-12,1.864100e-09, 2.020921e-15,-4.579188e-18,2.512191e-09,8.075043e-12)
      crit_val <- sum(k3_mod2*k3_coef2) + 0.05
      return (crit_val)
    # When sample size is between 71 and 90
    }else if (70<sample_size & sample_size<=90){
      k3_mod3 <- c(sample_size^2,lag^1.1,sample_size^(2)*lag^(1.1),sample_size^(2*2)*lag^(2*1.1),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.1),sample_size^(3*2)*lag^(3*1.1),lag^(4*1.1),lag^(5*1.1))
      k3_coef3 <- c(-6.411849e-06,4.572879e-03,-4.727665e-07,-1.792934e-12,6.884537e-10,2.380999e-16,-4.489660e-19,6.629214e-10,3.057809e-12)
      crit_val <- sum(k3_mod3*k3_coef3) + 0.05
      return (crit_val)
    # When sample size is between 91 and 120
    }else if (90<sample_size & sample_size<=120){
      k3_mod4 <- c(sample_size^0.4,lag^2,sample_size^(0.4)*lag^(2),sample_size^(2*0.4)*lag^(2*2),sample_size^(2*0.4),sample_size^(3*0.4)*lag^(2*2),sample_size^(3*0.4)*lag^(3*2),lag^(4*2),lag^(5*2))
      k3_coef4 <- c(1.340743e-02,-1.016702e-04,1.928474e-05,1.322170e-10,-2.001855e-03,-4.137456e-11,1.197990e-15,1.632013e-17,-1.283559e-21)
      crit_val <- sum(k3_mod4*k3_coef4) + 0.05
      return (crit_val)
    # When sample size is between 121 and 200
    }else if (120<sample_size & sample_size<=200){
      #c(n^s,m^p,n^(s)*m^(p),n^(2*s)*m^(2*p),n^(2*s),n^(3*s)*m^(2*p),n^(3*s)*m^(3*p),m^(4*p),m^(5*p))
      k3_mod5 <- c(sample_size^0.3,lag^0.9,sample_size^(0.3)*lag^(0.9),sample_size^(2*0.3)*lag^(2*0.9),sample_size^(2*0.3),sample_size^(3*0.3)*lag^(2*0.9),sample_size^(3*0.3)*lag^(3*0.9),lag^(4*0.9),lag^(5*0.9))
      k3_coef5 <- c(-7.174943e-04,1.332538e-02,-2.923840e-03,-1.701132e-05,5.324415e-05,4.480828e-06,-1.856506e-08,1.252202e-08,-1.341736e-11)
      crit_val <- sum(k3_mod5*k3_coef5) + 0.05
      return (crit_val)
    # When sample size is larger than 200
    }else if (sample_size>200){
      k3_mod6 <- c(sample_size^2,lag^1.7,sample_size^(2)*lag^(1.7),sample_size^(2*2)*lag^(2*1.7),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.7),sample_size^(3*2)*lag^(3*1.7),lag^(4*1.7),lag^(5*1.7))
      k3_coef6 <- c(6.806470e-07,-1.040297e-05,4.166030e-10,-1.187976e-18,-8.639787e-12,7.958527e-24, 1.593097e-30,2.193833e-17,-7.843942e-22)
      crit_val <- sum(k3_mod6*k3_coef6) + 0.05
      return (crit_val)
    }
  # models for analyzing 5 time series
  }else if(num_of_ts==5){
    if (sample_size<=50){
      # When sample size is less than or equal to 50
      k5_mod1 <- c(sample_size^1.5,lag^1.4,sample_size^(1.5)*lag^(1.4),sample_size^(2*1.5)*lag^(2*1.4),sample_size^(2*1.5),sample_size^(3*1.5)*lag^(2*1.4),sample_size^(3*1.5)*lag^(3*1.4),lag^(4*1.4),lag^(5*1.4))
      k5_coef1 <- c(-1.380839e-04,3.109495e-03,-4.981331e-06,-8.967959e-10,2.699682e-07,2.399881e-12,-2.086438e-15,7.256124e-10,-1.360533e-12)
      crit_val <- sum(k5_mod1*k5_coef1) + 0.05
      return (crit_val)
      # When sample size is between 51 and 70
    }else if (50<sample_size & sample_size<=70){
      k5_mod2 <- c(sample_size^2,lag^1.1,sample_size^(2)*lag^(1.1),sample_size^(2*2)*lag^(2*1.1),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.1),sample_size^(3*2)*lag^(3*1.1),lag^(4*1.1),lag^(5*1.1))
      k5_coef2 <- c(-1.591810e-05,8.314120e-03,-1.356343e-06,-1.535674e-11,2.715568e-09,3.241555e-15,-6.979655e-18,5.175389e-09,3.836440e-12)
      crit_val <- sum(k5_mod2*k5_coef2) + 0.05
      return (crit_val)
      # When sample size is between 71 and 90
    }else if (70<sample_size & sample_size<=90){
      k5_mod3 <- c(sample_size^2,lag^1,sample_size^(2)*lag^(1),sample_size^(2*2)*lag^(2*1),sample_size^(2*2),sample_size^(3*2)*lag^(2*1),sample_size^(3*2)*lag^(3*1),lag^(4*1),lag^(5*1))
      k5_coef3 <- c(-1.053974e-05,1.123504e-02,-1.288667e-06,-5.857440e-12,1.188668e-09,8.350555e-16,-2.745446e-18,1.405021e-09,8.048967e-11)
      crit_val <- sum(k5_mod3*k5_coef3) + 0.05
      return (crit_val)
      # When sample size is between 91 and 120
    }else if (90<sample_size & sample_size<=120){
      k5_mod4 <- c(sample_size^0.4,lag^2,sample_size^(0.4)*lag^(2),sample_size^(2*0.4)*lag^(2*2),sample_size^(2*0.4),sample_size^(3*0.4)*lag^(2*2),sample_size^(3*0.4)*lag^(3*2),lag^(4*2),lag^(5*2))
      k5_coef4 <- c(2.258185e-02,-1.793009e-04,3.262378e-05,3.325682e-10,-3.393526e-03,-8.394449e-11,2.326269e-15,1.261262e-17,-1.660552e-21)
      crit_val <- sum(k5_mod4*k5_coef4) + 0.05
      return (crit_val)
      # When sample size is between 91 and 120
    }else if (120<sample_size & sample_size<=200){
      #c(n^s,m^p,n^(s)*m^(p),n^(2*s)*m^(2*p),n^(2*s),n^(3*s)*m^(2*p),n^(3*s)*m^(3*p),m^(4*p),m^(5*p))
      k5_mod5 <- c(sample_size^0.3,lag^0.9,sample_size^(0.3)*lag^(0.9),sample_size^(2*0.3)*lag^(2*0.9),sample_size^(2*0.3),sample_size^(3*0.3)*lag^(2*0.9),sample_size^(3*0.3)*lag^(3*0.9),lag^(4*0.9),lag^(5*0.9))
      k5_coef5 <- c(2.929246e-03,1.602767e-02,-3.472773e-03,-2.137544e-05,-7.864162e-04,5.536952e-06,-2.283467e-08,1.623846e-08,-1.987672e-11)
      crit_val <- sum(k5_mod5*k5_coef5) + 0.05
      return (crit_val)
      # When sample size is larger than 200
    }else if (sample_size>200){
      k5_mod6 <- c(sample_size^2,lag^1.6,sample_size^(2)*lag^(1.6),sample_size^(2*2)*lag^(2*1.6),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.6),sample_size^(3*2)*lag^(3*1.6),lag^(4*1.6),lag^(5*1.6))
      k5_coef6 <- c(8.619332e-07,-2.801229e-05,1.003208e-09,-4.901480e-18,-1.153743e-11,3.303065e-23,-4.621521e-29,2.814372e-16,-1.723768e-20)
      crit_val <- sum(k5_mod6*k5_coef6) + 0.05
      return (crit_val)
    }
  # models for analyzing 10 time series
  }else if(num_of_ts==10){
    if (sample_size<=50){
      # When sample size is less than or equal to 50
      k10_mod1 <- c(sample_size^1.4,lag^0.8,sample_size^(1.4)*lag^(0.8),sample_size^(2*1.4)*lag^(2*0.8),sample_size^(2*1.4),sample_size^(3*1.4)*lag^(2*0.8),sample_size^(3*1.4)*lag^(3*0.8),lag^(4*0.8),lag^(5*0.8))
      k10_coef1 <- c(-1.722358e-03,1.364009e-01,-6.268533e-04,-2.528472e-07,7.190455e-06,1.392121e-09,-1.879418e-11,-2.065772e-06,2.537805e-07)
      crit_val <- sum(k10_mod1*k10_coef1) + 0.05
      return (crit_val)
      # When sample size is between 51 and 70
    }else if (50<sample_size & sample_size<=70){
      k10_mod2 <- c(sample_size^2,lag^0.9,sample_size^(2)*lag^(0.9),sample_size^(2*2)*lag^(2*0.9),sample_size^(2*2),sample_size^(3*2)*lag^(2*0.9),sample_size^(3*2)*lag^(3*0.9),lag^(4*0.9),lag^(5*0.9))
      k10_coef2 <- c(-4.521172e-05,4.374216e-02,-8.761956e-06,-1.156707e-10,8.618181e-09,2.863483e-14,-1.805493e-16,-9.482726e-08,6.187954e-09)
      crit_val <- sum(k10_mod2*k10_coef2) + 0.05
      return (crit_val)
      # When sample size is between 71 and 90
    }else if (70<sample_size & sample_size<=90){
      k10_mod3 <- c(sample_size^2,lag^1,sample_size^(2)*lag^(1),sample_size^(2*2)*lag^(2*1),sample_size^(2*2),sample_size^(3*2)*lag^(2*1),sample_size^(3*2)*lag^(3*1),lag^(4*1),lag^(5*1))
      k10_coef3 <- c(-2.025258e-05,2.128571e-02,-2.325250e-06,-1.193464e-11,2.243568e-09,1.621225e-15,-4.745842e-18,6.133523e-09,1.117398e-10)
      crit_val <- sum(k10_mod3*k10_coef3) + 0.05
      return (crit_val)
      # When sample size is between 91 and 110
    }else if (90<sample_size & sample_size<=120){
      k10_mod4 <- c(sample_size^2,lag^2,sample_size^(2)*lag^(2),sample_size^(2*2)*lag^(2*2),sample_size^(2*2),sample_size^(3*2)*lag^(2*2),sample_size^(3*2)*lag^(3*2),lag^(4*2),lag^(5*2))
      k10_coef4 <- c(6.560008e-06,-2.814827e-05,8.957106e-09,-2.888567e-16,-5.523747e-10,1.195876e-20,3.044985e-25,2.124154e-16,-1.308718e-20)
      crit_val <- sum(k10_mod4*k10_coef4) + 0.05
      return (crit_val)
      # When sample size is between 111 and 200
    }else if (110<sample_size & sample_size<=200){
      k10_mod5 <- c(sample_size^2,lag^1.7,sample_size^(2)*lag^(1.7),sample_size^(2*2)*lag^(2*1.7),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.7),sample_size^(3*2)*lag^(3*1.7),lag^(4*1.7),lag^(5*1.7))
      k10_coef5 <- c(3.954236e-06,-5.704072e-05,4.373655e-09,-5.817395e-17,-1.124866e-10,8.912555e-22,2.764302e-27,8.600325e-16,-6.547621e-20)
      crit_val <- sum(k10_mod5*k10_coef5) + 0.05
      return (crit_val)
      # When sample size is larger than 200
    }else if (sample_size>200){
      k10_mod6 <- c(sample_size^2,lag^1.7,sample_size^(2)*lag^(1.7),sample_size^(2*2)*lag^(2*1.7),sample_size^(2*2),sample_size^(3*2)*lag^(2*1.7),sample_size^(3*2)*lag^(3*1.7),lag^(4*1.7),lag^(5*1.7))
      k10_coef6 <- c(1.589487e-06,-3.609943e-05,1.143797e-09,-2.931991e-18,-2.076415e-11,1.729112e-23,1.783932e-28,6.156433e-17,-2.540740e-21)
      crit_val <- sum(k10_mod6*k10_coef6) + 0.05
      return (crit_val)
    }
  }
}

MTS_T1_Correction(10,80,76)  
