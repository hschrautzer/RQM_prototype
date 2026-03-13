# RQM_prototype

Prototype for Rayleigh Quotient Minimization

## Finite Difference Part

Here we describe the computation of the action of the Hessian $H\in\mathbb{R}^{2N\times 2N}$ on the thin orthonormal matrix $X\in\mathbb{R}^{2N\times p}$ , where $p$ indicates the dimension of the subspace of $H$ that the RQM method is trying to estimate.

We compute the action on the individual columns of $X$

$$
HX=(Hx_1,\dots, Hx_p)~,
$$

with $x_i\in\mathbb{R}^{2N}$ and $|x_i|=1$. Say $g(m)\in\mathbb{R}^{3N}$ is the energy gradient of a magnetic configuration $m=(m_1^x,m_1^y,m_1^z,\dots,m_N^x,m_N^y,m_N^z)\in\mathbb{R}^{2N}$ of $N$ magnetic moments. Note, despite that this gradient is represented in $3N$-dimensional embedding space it is tangent to the magnetic configuration ($m\cdot g(m)=0$). It can be projected to the local $2N$-dimensional tangent space using the projection matrix $U\in\mathbb{R}^{3N\times 2N}$ with $g_{2N}=U^Tg\in\mathbb{R}^{2N}$. The gradient is represented in $3N$ embedding space to avoid frequent basis changes in the following.

The idea to compute the action of the Hessian is now expanding the gradient:

$$
g(m+\epsilon x_i)= g(m)+U\epsilon H x_i+\mathcal{O}(\epsilon^2)~,
$$

with a small displacement parameter $\epsilon$. Re-arranging yields:

$$
H x_i=U^T\frac{g(m+\epsilon x_i)-g(m)}{\epsilon}+\mathcal{O}(\epsilon)~,
$$

which is first-order accurate in $\epsilon$.

**Note**: The above is not correct since we did not care about the curvature of the spherical manifold our magnetic configuration is living on. Denoting the retraction of $m$ along the direction $x_i$ with the displacement parameter $\epsilon$ by $m'=R_{\epsilon x_i}(m)$ and the parallel transport of a vector $g$ along the direction $x_i$ with the displacement parameter $\epsilon$ and the anchor $m$ by $P_{m,\epsilon x_i}(g)$ correctly the above should read:

$$
Hx_i = U^T \frac{P_{m',-\epsilon x_i}[g(m')]-g(m)}{\epsilon}+\mathcal{O}(\epsilon)~.
$$

The above is implemented in `magnetization.py` as `fd_simple_HcolX()`. However, can we do better than $\epsilon$-accuracy? For that purpose the original implementation in Spinaker uses a Richardson Extrapolation of the finite difference displacement parameter $\epsilon$. For better readability omit the retraction and parallel transport notations in the following.

Consider now two times the same taylor expansion but once for $\epsilon$ and once for $\epsilon/2$:

$$
Hx_i=U^T\frac{g(m+\epsilon x_i)-g(m)}{\epsilon}+C\epsilon+\mathcal{O}(\epsilon^2)=A(\epsilon)+\mathcal{O}(\epsilon^2)\\
Hx_i=U^T\frac{g(m+\frac{\epsilon}{2}x_i)-g(m)}{\epsilon/2}+C\frac{\epsilon}{2}+\mathcal{O}(\epsilon^2)=A(\epsilon/2)+\mathcal{O}(\epsilon^2)~,
$$

where $C$ is some constant yielding the value of the third derivative of the energy and representing the linear error of this scheme and $A(\epsilon)$ represents the approximation including the linear error. This linear error can be eliminated by computing

$$
Hx_i = 2A(\epsilon/2)-A(\epsilon)+\mathcal{O}(\epsilon^2)~.
$$

The benefit is that we can estimate the absolute error of this approximation by

$$
e(\epsilon/2)\approx ||A(\epsilon/2)-A(\epsilon)||~.
$$

Therefore in practice one can define a maximum allowed absolute error that the finite difference calculation may have. Let's say this desired error is $e_R$ and that we started with some $\epsilon$ and $e(\epsilon/2)>e_R$. With $e(\epsilon/2)= C\frac{\epsilon}{2}+\mathcal{O}(\epsilon^2)$ we have $C\approx \frac{2e(\epsilon/2)}{\epsilon}$ the new step to reach the desired error might be:

$$
\epsilon_{new}= 
$$

**Note from Hendrik to Hendrik**: It seems that my error reduction scheme in spinaker thinks that the accuracy of the scheme is p=2 while it is p=1 in reality. Do I want to fix that?

**Note**: This is of course not perfect. If I would have to do this again, I would consider a central-finite difference scheme using Richardson Extrapolation and use of course an estimate of the relative error.



$$
H x_i=U^T\frac{g(m+\frac{\epsilon}{2} x_i)-g(m)}{\frac{\epsilon}{2}}+\mathcal{O}(\epsilon)~,
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
