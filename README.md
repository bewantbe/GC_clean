GC_clean
========

A collection of GNU Octave (also Matlab compatible) scripts and functions that calculate and investigate Granger Causality(GC). We developed these code for the study of neuroscience (See the last three [*References*](#references)).

The code in this package is capable to compute conditional GC of a thousand variables (Within an hour in usual PC). There are code to compute both time domain and frequency domain GC. Statistical test (p-value, confidence interval) can also be obtained (so far time domain only).

Note that these code are mainly for research purposes. While the correctness is one of our main concern, the robustness for all input cases is not guaranteed (see also *Function overview*).

Function overview
-----------------

There are various ways to calculate GC, some are fast, some are more robust to ill-condition problem (but slow and maybe consume a lots memory).

The calculation of GC usually can be divided into three steps (time domain):

- Compute [autocorrelation](http://en.wikipedia.org/wiki/Autocorrelation) of input time series.
- Use the autocorrelation to solve the [vector autoregression](http://en.wikipedia.org/wiki/Vector_autoregression) problem, and get residual variances.
- Use the residual variances obtained above to calculate GC.

Or, in frequency domain it is:

- Estimate the [spectral density](http://en.wikipedia.org/wiki/Spectral_density_estimation) of input time series.
- Perform the (minimum phase) spectral factorization to get the transfer function and residual variances.
- Use the residual variance obtained above to calculate GC.

Both time domain and frequency domain method can be used to calculate frequency domain GC.

Ill-condition data here is the time series that some variables are highly correlated. The data after a low/high/band pass filter could be Ill-conditioned. Resampling to the pass band can turn it to "good-condition" problem.

### Functions:

See GCcal/readme_GCcal.txt for details. Here list some main functions:

* Calculate GC in time domain (under `GCcal/`).

  - `nGrangerT.m`

        Use Yule-Walker equation for the 2nd step. It's simple and fast (for fewer than hundred variables). Not stable for short data (e.g. <1e4) and ill-condition data.

  - `pos_nGrangerT.m, pos_nGrangerT2.m, pos_nGrangerT_qr.m`

        Solve the finite data point linear least square problem in the 2nd step. Essentially solving the normal equation `"X'*X*A =X'*Y"`. The robustness relys on the '/' operator in Octave (or Matlab) and accumulation round-off error in the 1st step.

        `pos_nGrangerT.m` calculate as the definition. Much more stable than `nGrangerT.m`.

        `pos_nGrangerT2.m` calculate as its definition, but use some trick to do it fast (as fast as `nGrangerT.m`) and use much less memory. Slightly more Round-off error than `pos_nGrangerT.m`.

        `pos_nGrangerT_qr.m` solve the linear least square problem by [QR decomposition](http://en.wikipedia.org/wiki/QR_decomposition) in "X*A=Y" (performed by the Octave operator "/"). Effectively square rooted the [condition number](http://en.wikipedia.org/wiki/Condition_number). This is the most stable (most slow) and accurate method in this package.

  - `nGrangerTfast.m`

        Almost as stable as `pos_nGrangerT2.m`, and mathematically the same as `pos_nGrangerT2.m` for non-singular problem. It is fast (even faster than `nGrangerT.m`) for case of tens and hundreds of variables.

  - `RGrangerTLevinson.m`

        Same method as `nGrangerTfast.m`, but use [Levinson recursion](http://en.wikipedia.org/wiki/Levinson_recursion) to perform the matrix inversion. Much faster than even `nGrangerTfast.m` for the case of hundreds and a thousand variables. Use it as `GC = RGrangerTLevinson( getcovpd(X, m) );`. Stability is worse than `nGrangerT.m`, but still enough for non-ill-condition problem (e.g. cond<1e6).

  - `pairGrangerT.m`

        Calculate pairwise GC. Same stability as `nGrangerT.m`.


* Calculate frequency domain GC (under `GCcal/`).

  - `nGrangerF.m`

        As not stable as `nGrangerT.m`. And very slow for large variables.


* Related functions

  - `GCcal/gc_prob_nonzero.m`

        Get p-value for the corresponding GC. Used for significance test.

  - `GCcal/gc_prob_intv.m`

        Get confidence interval of GC.

  - `GCcal/chooseOrderAuto.m, GCcal/chooseROrderFull.m`

        Get suitable regression order for GC. Based on Akaike information criterion (AIC) or Bayesian information criterion (BIC).
        
        `chooseOrderAuto.m` is a fully automatic routine that use Levinson recursion for "good-condition" problem, and use method like `pos_nGrangerT2.m` for ill-condition problem.
        
        `chooseROrderFull.m` can also return autoregression order.

  - `GCcal_spectrum/mX2S_wnd.m`
  
        Estimate spectral density with window function applied.

  - `GCcal_spectrum/mX2S_nuft.m`

        Estimate spectral density for non-uniformly sampled data.

References<a name="references"></a>
----------

* John Geweke (1982) Measurement of Linear Dependence and Feedback Between Multiple Time Series. Journal of the American Statistical Association, Vol. 77, No. 378, pp. 304-313. [doi:10.2307/2287238](http://www.jstor.org/stable/2287238)

* John Geweke (1984) Measures of Conditional Linear Dependence and Feedback Between Time Series. Journal of the American Statistical Association, Vol. 79, No. 388, pp. 907-915. [doi:10.2307/2288723](http://www.jstor.org/stable/2288723)

* Ding, M., Chen, Y., & Bressler, S. L. (2006). 17 Granger Causality: Basic Theory and Application to Neuroscience. Handbook of time series analysis: recent theoretical developments and applications, 437.

* Zhou D, Xiao Y, Zhang Y, Xu Z, Cai D (2013) Causal and structural connectivity of pulse-coupled nonlinear networks. [Physical Review Letters 111: 054102](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.054102)

* Zhou D, Xiao Y, Zhang Y, Xu Z, Cai D (2014) Granger Causality Network Reconstruction of Conductance-Based Integrate-and-Fire Neuronal Systems. PLoS ONE 9(2): e87636. [doi:10.1371/journal.pone.0087636](http://dx.plos.org/10.1371/journal.pone.0087636)

* Zhou D, Zhang Y, Xiao Y, Cai D (2014) Reliability of the Granger causality inference. New J. Phys. 16 043016. [doi:10.1088/1367-2630/16/4/043016](http://iopscience.iop.org/1367-2630/16/4/043016)

