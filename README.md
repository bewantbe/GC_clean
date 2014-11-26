GC_clean
========

A collection of GNU Octave (also Matlab compatible) scripts and functions that calculate and investigate Granger Causality(GC). We developed these code for the study of neuroscience.

Note that these code are mainly for research purposes. While the correctness is one of our main concern, the robustness of all input cases is not guaranteed (see also *Function overview*).

Function overview
-----------------

There are various ways to calculate GC, some are fast, some are more robust to ill-condition case (but slow and maybe memory hungry).

The calculation of GC usually can be divided into three steps (time domain):

- Compute [autocorrelation](http://en.wikipedia.org/wiki/Autocorrelation) of input time series.
- Use the autocorrelation to solve the [vector autoregression](http://en.wikipedia.org/wiki/Vector_autoregression) problem, and get residual variances.
- Use the residual variances obtained above to calculate GC.

Or, in frequency domain it is:

- Estimate the [spectral density](http://en.wikipedia.org/wiki/Spectral_density_estimation) of input time series.
- Perform the (minimum phase) spectral factorization to get the transfer function and residual variances.
- Use the residual variance obtained above to calculate GC.

Both time domain and frequency domain method can be used to calculate frequency domain GC.

### Functions:

* Calculate GC in time domain.
  - `nGrangerT.m`

        Use Yule-Walker equation for the 2nd step. It's fast (for fewer than hundred variables). bad for short data (<1e4).

  - `pos_nGrangerT.m, pos_nGrangerT2.m, pos_nGrangerT_qr.m`

        Solve the finite data point linear least square problem in the 2nd step. Essentially solving the normal equation "X'*X*A =X'*Y". The robustness relys on the '/' operator in Octave (or Matlab) and accumulation round-off error in the 1st step.

References
----------

* Zhou D, Xiao Y, Zhang Y, Xu Z, Cai D (2013) Causal and structural connectivity of pulse-coupled nonlinear networks. [Physical Review Letters 111: 054102](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.054102)

* Zhou D, Xiao Y, Zhang Y, Xu Z, Cai D (2014) Granger Causality Network Reconstruction of Conductance-Based Integrate-and-Fire Neuronal Systems. PLoS ONE 9(2): e87636. [doi:10.1371/journal.pone.0087636](http://dx.plos.org/10.1371/journal.pone.0087636)

* Zhou D, Zhang Y, Xiao Y, Cai D (2014) Reliability of the Granger causality inference. New J. Phys. 16 043016. [doi:10.1088/1367-2630/16/4/043016](http://iopscience.iop.org/1367-2630/16/4/043016)

