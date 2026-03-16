# Things that should be updated in Spinaker

- spk_rqm.f90: 414: taking the square root of the Frobenius norm is missing (compare line 214)
- spk_rqm.f90: 390: the rq is not needed anymore
- spk_finites.f90: 386: the sqrt assumes the accuracy to be of order 2. However, it is just 1.
