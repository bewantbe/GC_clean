GC_clean
========

A collection of GNU Octave (also Matlab compatible) scripts and functions that calculate and investigate Granger Causality(GC). We developed these code for the study of neuroscience (See the last three [*References*](#references)).

The code in this package is capable to compute conditional GC of a thousand variables (Within an hour in usual PC). There are code to compute both time domain and frequency domain GC. Statistical test (p-value, confidence interval) can also be obtained (so far time domain only).

Note that these code are mainly for research purposes. While the correctness is one of our main concern, the robustness for all input cases is not guaranteed (see also *Function overview*).

Function overview
-----------------

There are various ways to calculate GC, some are fast, some are more robust to ill-condition problem (but slow and maybe consume a lots memory).

The calculation of GC usually can be divided into three steps (time domain <a name="gc-step"></a>):

1. Compute multivariate [autocorrelation](http://en.wikipedia.org/wiki/Autocorrelation) (also called [correlation function](http://en.wikipedia.org/wiki/Correlation_function)) of input time series.
2. Use the autocorrelation to solve the [vector autoregression](http://en.wikipedia.org/wiki/Vector_autoregression) problem, and get residual variances.
3. Use the residual variances obtained above to calculate GC.

Or, in frequency domain it is:

1. Estimate the [spectral density](http://en.wikipedia.org/wiki/Spectral_density_estimation) of input time series.
2. Perform the (minimum phase) spectral factorization to get the transfer function and residual variances.
3. Use the residual variance obtained above to calculate GC.

Both time domain and frequency domain method can be used to calculate frequency domain GC.

Knowing d

Ill-condition data here is the time series that some variables are highly correlated. The data after a low/high/band pass filter will be ill-conditioned in general. Resampling to the pass band can turn it to "good-condition" problem.

### Functions:

See GCcal/readme_GCcal.txt for details (in Chinese). Here list some main functions:

* Calculate GC in time domain (under `GCcal/`).

  - `nGrangerT.m`

        Use Yule-Walker equation for the 2nd step. It's simple and fast (for fewer than hundred variables). Not stable for short data (e.g. <1e4) and ill-condition data.

  - `pos_nGrangerT.m, pos_nGrangerT2.m, pos_nGrangerT_qr.m, pos_nGrangerT_qrm.m`

        Solve the finite data point linear least square problem in the 2nd step. Essentially solving the normal equation `"X'*X*A =X'*Y"`. The robustness relys on the '/' operator in Octave (or Matlab) and accumulation round-off error in the 1st step.

        `pos_nGrangerT.m` calculate as the definition. More stable than `nGrangerT.m`.

        `pos_nGrangerT2.m` calculate as its definition, but use some trick to do it fast (as fast as `nGrangerT.m`) and use much less memory. Slightly more Round-off error than `pos_nGrangerT.m`.

        `pos_nGrangerT_qr.m` solve the linear least square problem by [QR decomposition](http://en.wikipedia.org/wiki/QR_decomposition) in "X*A=Y" (performed by the operator "/"). Effectively square rooted the [condition number](http://en.wikipedia.org/wiki/Condition_number).

        `pos_nGrangerT_qrm.m` A variation of `pos_nGrangerT_qr.m`, which also put mean value into least square problem (instead of simply subtract it). This is the most stable (most slow) and accurate method in this package.

  - `nGrangerTfast.m`

        Almost as stable as `pos_nGrangerT2.m`, and mathematically equivalent to `pos_nGrangerT2.m`. It is fast (even faster than `nGrangerT.m`) for case of tens and hundreds of variables.

  - `RGrangerTLevinson.m`

        Same method as `nGrangerTfast.m`, but use [Levinson recursion](http://en.wikipedia.org/wiki/Levinson_recursion) to perform the matrix inversion. Much faster than even `nGrangerTfast.m` for the case of hundreds and a thousand variables. Use it as `GC = RGrangerTLevinson( getcovpd(X, m) );`. Stability is worse than `nGrangerT.m`, but still enough for non-ill-condition problem (e.g. cond<1e6).
        Note: to overcome the ill-condition problem, one may whiten the time series fist (see `WhiteningFilter.m`).

  - `pairGrangerT.m`

        Calculate pairwise GC. Same stability as `nGrangerT.m`.


* Calculate frequency domain GC (under `GCcal/`).

  - `nGrangerF.m`

        As not stable as `nGrangerT.m`. And very slow for large variables (not speed optimized, but should be easy to read and understand following [John Geweke (1984)](#references)).


* Related functions

  - `GCcal/gc_prob_nonzero.m`

        Get p-value ("probability" of nonzero) for the corresponding GC. Used for significance test.

  - `GCcal/gc_prob_intv.m`

        Get confidence interval of GC value.

  - `GCcal/chooseOrderAuto.m, GCcal/chooseROrderFull.m`

        Get suitable regression order for GC. Based on Akaike information criterion (AIC) or Bayesian information criterion (BIC).
        
        `chooseOrderAuto.m` is a fully automatic routine that use Levinson recursion for "good-condition" problem, and use method like `pos_nGrangerT2.m` for ill-condition problem.
        
        `chooseROrderFull.m` can also return autoregression order.

  - `GCcal_spectrum/mX2S_wnd.m`
  
        Estimate spectral density with window function applied.

  - `GCcal_spectrum/mX2S_nuft.m`

        Estimate spectral density for non-uniformly sampled data.

* `IFsimu_release_2.1.1.zip, prj_neuron_gc/raster_tuning`

    These are Integrate-and-Fire model neuron simulator. See the readme in `IFsimu_release_2.1.1.zip/exec/Readme.txt` for details. The latest version can be found https://bitbucket.org/bewantbe/ifsimu or alternatively https://github.com/bewantbe/point-neuron-network-simulator.

Speed
-----
(Machine: i5-2400, 8GB)

### Speed of the [1st step](#gc-step)

#### Function `getcovpd.m, getcovzpd.m`.

> O(p^2 * m * L)

len=1e5, sec | p=100 | 200 | 500 | 1000
:-----------:|:-----:|:---:|:---:|:-----:
od=20        | 1.339 | 3.96 | 17.02 | 59.7
od=40        | 2.67  | 7.71 | 32.14 | 118.1

(In many cases, you will need 10 times longer data, which the time above will be 10 times higher)

### Speed of the [2nd and 3rd step](#gc-step)

#### Function `RGrangerT.m` (similarly the 2nd step of `pos_nGrangerT2.m`).

> O(p^4 * m^3)

sec    | p=100 | 200 | 500
:-----:|:-----:|:---:|:---:
od=20  | 24.2  | 243.6  | 6870.3
od=40  | 120.4 | 1358.2 | oom

#### Function `RGrangerTfast.m`.

> O(p^3 * m^3)

sec    | p=100 | 200 | 500
:-----:|:-----:|:---:|:---:
od=20  | 0.721 | 4.06 | 52.1
od=40  | 4.63  | 27.2 | oom

#### Function `RGrangerTLevinson.m`.

> O( p^3 * m^2 * log(m) )

sec    | p=100 | 200 | 500  | 1000
:-----:|:-----:|:---:|:----:|:------:
od=20  | 0.434 | 1.987 | 19.5 | 110.5
od=40  | 1.662 | 8.353 | 76.9 | 433.9

(oom = out of memory)

References<a name="references"></a>
----------

* John Geweke (1982) Measurement of Linear Dependence and Feedback Between Multiple Time Series. Journal of the American Statistical Association, Vol. 77, No. 378, pp. 304-313. [doi:10.2307/2287238](http://www.jstor.org/stable/2287238)

* John Geweke (1984) Measures of Conditional Linear Dependence and Feedback Between Time Series. Journal of the American Statistical Association, Vol. 79, No. 388, pp. 907-915. [doi:10.2307/2288723](http://www.jstor.org/stable/2288723)

* Ding, M., Chen, Y., & Bressler, S. L. (2006). 17 Granger Causality: Basic Theory and Application to Neuroscience. Handbook of time series analysis: recent theoretical developments and applications, 437.

* Zhou D, Xiao Y, Zhang Y, Xu Z, Cai D (2013) Causal and structural connectivity of pulse-coupled nonlinear networks. [Physical Review Letters 111: 054102](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.054102)

* Zhou D, Xiao Y, Zhang Y, Xu Z, Cai D (2014) Granger Causality Network Reconstruction of Conductance-Based Integrate-and-Fire Neuronal Systems. PLoS ONE 9(2): e87636. [doi:10.1371/journal.pone.0087636](http://dx.plos.org/10.1371/journal.pone.0087636)

* Zhou D, Zhang Y, Xiao Y, Cai D (2014) Reliability of the Granger causality inference. New J. Phys. 16 043016. [doi:10.1088/1367-2630/16/4/043016](http://iopscience.iop.org/1367-2630/16/4/043016)

