# RQM_prototype

Prototype for Rayleigh Quotient Minimization

## Finite Difference Part

Here we describe the computation of the action of the Hessian $H\in\mathbb{R}^{2N\times 2N}$ on the thin orthonormal matrix $X\in\mathbb{R}^{2N\times p}$ , where $p$ indicates the dimension of the subspace of $H$ that the RQM method is trying to estimate.

We compute the action on the individual columns of $X$

$$
HX=(Hx_1,\dots, Hx_p)~,
$$

with $x_i\in\mathbb{R}^{2N}$ and $|x_i|=1$. Say $\bm{g}(\bm{m})\in\mathbb{R}^{3N}$ is the energy gradient of a magnetic configuration $\bm{m}=(m_1^x,m_1^y,m_1^z,\dots,m_N^x,m_N^y,m_N^z)\in\mathbb{R}^{2N}$ of $N$ magnetic moments. Note, despite that this gradient is represented in $3N$-dimensional embedding space it is tangent to the magnetic configuration ($\bm{m}\cdot\bm{g}(\bm{m})=0$). It can be projected to the local $2N$-dimensional tangent space using the projection matrix $U\in\mathbb{R}^{3N\times 2N}$ with $\bm{g}_{2N}=U^T\bm{g}\in\mathbb{R}^{2N}$. The gradient is represented in $3N$ embedding space to avoid frequent basis changes in the following.

The idea to compute the action of the Hessian is now expanding the gradient:

$$
\bm{g}(\bm{m}+\epsilon\bm{x}_i)= \bm{g}(\bm{m})+U\epsilon H\bm{x}_i+\mathcal{O}(\epsilon^2)~,
$$

with a small displacement parameter $\epsilon$. Re-arranging yields:

$$
H\bm{x}_i=U^T\frac{\bm{g}(\bm{m}+\epsilon\bm{x}_i)-\bm{g}(\bm{m})}{\epsilon}+\mathcal{O}(\epsilon)~,
$$

which is first-order accurate in $\epsilon$.

**Note**: The above is not correct since we did not care about the curvature of the spherical manifold our magnetic configuration is living on. Denoting the retraction of $\bm{m}$ along the direction $\bm{x}_i$ with the displacement parameter $\epsilon$ by $\bm{m}'=\mathcal{R}_{\epsilon\bm{x}_i}(\bm{m})$ and the parallel transport of a vector $\bm{g}$ along the direction $\bm{x}_i$ with the displacement parameter $\epsilon$ and the anchor $\bm{m}$ by $\mathcal{P}_{\bm{m},\epsilon\bm{x}_i}(\bm{g})$ correctly the above should read:

$$
Hx_i = U^T \frac{\mathcal{P}_{\bm{m}',-\epsilon\bm{x}_i}[\bm{g}(\bm{m}')]-\bm{g}(\bm{m})}{\epsilon}+\mathcal{O}(\epsilon)~.
$$

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
