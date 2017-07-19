- generic description of the log-likelihood, and inherent differences of fitting tau, sigma, sigma_hat, etc. with illustrations :  DRW_log_likelihood.ipynb


- simulating long light curves that we then chop at different lengths  : Celerite_Kozlowski_2017_simulate.ipynb 
  
- plotting the results : reading  in the results of simulating one very well sampled, long LC, that is then chopped at various lengths ,,,, Celerite_Kozlowski_2017_plot.ipynb    


- simulating various light curves, and eventually reproducing Chelsea Fig. 15 : Celerite_fit_simulated_data.ipynb  

- fitting with Celerite LCs that Zeljko simulated now in 2017 : the input properties for each of the 1000 light curves  : were  $\tau_{in} = 100 $ days ,  $\sigma_{in} = 0.2  $ mag,   $SF_{\infty} = 0.2*\sqrt{2} = 0.2828$ mag ,  length = $l= 20 \tau$ , dt = 5 days  , random sampling  from a uniform distribution  , yerr = 0.001  mag , 400 points  : Celerite_fit_LC_ZI.ipynb 

- fitting the short simulated DRW  light curves (made by Zeljko in 2013  : there were three flavors,  short, med, long, but I could only find short,  and fit params from Chelsea’s code for all lengths.).  Using Chelsea’s result for short , I fit with Celerite short LCs and plot comparison against Chelsea’s results  : DRW_fit_Chelsea_lc.ipynb  
