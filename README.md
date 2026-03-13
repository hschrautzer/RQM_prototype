# RQM_prototype

Prototype for Rayleigh Quotient Minimization

## Finite Difference Part

Here we describe the computation of the action of the Hessian $H\in\mathbb{R}^{2N\times 2N}$ on the thin orthonormal matrix $X\in\mathbb{R}^{2N\times p}$ , where $p$ indicates the dimension of the subspace of $H$ that the RQM method is trying to estimate.

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
