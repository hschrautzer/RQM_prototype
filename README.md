# RQM_prototype

Prototype for Rayleigh Quotient Minimization

## Finite Difference Part

The computation of the action of the Hessian $H\in\mathbbÖÄ$  (2N by 2N) on the thin matrix orthonormal `X` (2N by p) first is
separated into the contributions of the individual columns of $`X=(x_1,...,x_p)`$: `Hx_i`with `xi`beeing of dimension 2N.

How can we compute `Hx_i`? The answer is finite differences:

```
Hx_i ~ [g(m+epsilon*x_i) - g(m)] / epsilon
```

with the gradient `g` of the magnetization `m`
The above does not include the neccesary retractions and transport expressions for the sake of clarity. Note, these are
crucial for correct numerical results. However, the above is a simple finite difference (forward) scheme with 1 stencil
point and a fixed steplength (implemented as `fd_simple_HcolX()` in `magnetization.py`. Therefore, the error will be
large and we can do better.

One possibility (which is used in Spinaker) is the Richardson extropolation.
