# Rayleigh Quotient Minimization Prototype

This is the Python prototype for Rayleigh Quotient Minimization implemented in the atomistic spin simulation code [Spinaker](https://journals.aps.org/prb/abstract/10.1103/z673-hhnp) in Fortran that shall be implemented in the micromagnetic code MuMax in Go.

Important reads for the below concepts:

[1]: [Edelman et al. (1994)](https://epubs.siam.org/doi/abs/10.1137/S0895479895290954?casa_token=I4ARBNfZzhwAAAAA:oSEat7Mqh8itgzWGutFX5PzQHDCzOr58Kkz48ZH9zqdKa88pfeKoxiyJ6pYxOLqz-RMkqCJGBOM)

[2]: [Book (free download): Absil (1999): Optimization Algorithms on Matrix Manifolds](https://press.princeton.edu/absil?srsltid=AfmBOorxNB9X7qOZCXIqdUmKyVdfD1EQ6dGdVGSAbjrxLtcPvbSiZR1f)

## Concept

The Hessian of the energy of the system for a particular magnetic configuration $H\in\mathbb{R}^{2N\times 2N}$ is symmetric. Consider the quadratic form:

$$
R(X)=tr(X^THX)~,
$$

which is the generalized Rayleigh Quotient with $X\in St_{2N\times p}$ element of the non-compact Stiefel manifold [2]:

$$
St_{2N\times p}=\{X\in\mathbb{R}^{2N\times p} |X^TX=I_p\}~,
$$

where $I_p$ is the $p\times p$ unit matrix. In other words: The columns of $X$ form an orthonormal $p$-frame in $\mathbb{R}^{2N}$. Take another look at the generalized Rayleigh Quotient for $X=(x_1,...,x_p)$ formed of $p$-columns with $x_i\in\mathbb{R}^{"N}$. Then:

$$
X^THX=\begin{bmatrix}x_1^THx_1 & \dots  \\ \dots & \ddots \end{bmatrix}
$$

and the trace becomes:

$$
R(X)=\sum\limits_{i=1}^p x_i^THx_i
$$

which is the sum over the Rayleigh quotients of the individual orthonormal vectors. Intuitively, one can already guess that the global minimum of $R(X)$ corresponds to the lowest subspace somehow. Let's take a look. Since $H$ is symmetric, we can diagonalize it: $H=D\Lambda D^T$, where $D$ is orthogonal ($D^TD=I_{2N}$) and $\Lambda$ is diagonal with $\Lambda = \text{diag}(\lambda_1,\dots,\lambda_{2N})$. Assume now the eigenvalues are ordered: $\lambda_1\leq \lambda_2\leq \dots\leq\lambda_{2N}$ and substitute into the generalized Rayleigh Quotient:

$$
R(X)=tr(X^TD\Lambda D^T X)
$$

Now define $Y=D^TX$ and since both $D$ and $X$ have orthonormal columns we have $Y^TY=I_p$. Therefore, $Y$ is also a point on the Stiefel manifold (a different one than $X$). Now compute the trace:

$$
R(X)=\sum\limits_{i=1}^{p}(Y^T\Lambda Y)_{ii}=\sum\limits_{i=1}^{p}\sum\limits_{j=1}^{2N}\lambda_j Y_{ji}^2=\sum\limits_{j=1}^{2N}\lambda_{j}\sum\limits_{i=1}^pY_{ji}^2
$$

with defining $||y_i||^2$ as the squared norm of the $i$-th row of $Y$ we have

$$
R(X)=\sum\limits_{i=1}^{2N}\lambda_i ||y_i||^2~.
$$

Finally, with $\sum_{i=1}^{2N}||y_i||^2=tr(Y^TY)=tr(I_p)=p$, we can see that the global minimum of $R(X)$ is exactly given by the $p$ smallest eigenvalues $\lambda_1 ,..., \lambda_p $.

However, consider any orthogonal matrix $Q\in\mathbb{R}^{p\times p}$ then:

$$
R(XQ)=R(X).
$$

So the function "Rayleigh Quotient" is actually defined on the Grassmann manifold, while we optimize it using coordinates on the Stiefel manifold.

### Algorithm

The generalized Rayleigh Quotient $R(X)=tr(X^THX)$ is minimized under the Stiefel constraint $X^TX=I_{p}$. Its extrema correspond to eigenvectors of $H$ and the function *really* depends only on the subspace spanned by $X$ (Grassmann manifold), but the Stiefel manifold provides convenient coordinates.

The **gradient** of the Rayleigh Quotient with respect to $X$ is given by [1]:

$$
\nabla R(X)=2[HX-X(X^THX)].
$$

We implement an [L-BFGS solver](https://www.sciencedirect.com/science/article/pii/S0010465520303696?casa_token=hOH6c63JhsMAAAAA:zHB6lHMSZM5LIwAPapgHN3FITDQXENXsUOB7M9SWJWyrgQpX3VzR08McKNgkI3uK8bZpteL7jw) for the Rayleigh Quotient based on the above gradient and adhering to the below geometrical concepts of the Grassmanian.


## Optimization on Riemannnian matrix manifolds

### Spin Space (Configuration Space)


### Grassmanian

The hierarchy of involved manifolds can be explained nicely by the concept of quotient spaces [1]. Start with the orthogonal group $O(2N)=\{Q\in\mathbb{R}^{2N\times 2N}|Q^TQ=I_{2N}\}$ that consists of all orthonormal based of $\mathbb{R}^{2N}$. Each $X\in St_{2N,p}$ is a set of $p$ orthonormal vectors in $\mathbb{R}^{2N}$. So, every matrix $Q\in O(2N)$ determines a point on the Stiefel manifold but this representation is not unique since any choice can be made for the remaining $2N-p$ columns. Thus, the Stiefel manifold can be represented as a quotient space of the orthogonal group:

$$
St_{2N,p}\cong O(2N)/O(2N-p).
$$

Similarly, the Grassmann manifold is obtained by identifying those matrices in $St_{2N,p}$ whose columns span the same subspace with the equivalency relation $Y\sim Q X$ with $Q\in\mathbb{R}^{p\times p}$ orthogonal.

$$
G_{2N,p}\cong St_{2N,p}/O(p)~.
$$

### Retraction

A projection like retraction on the Grassmannian is simply computing the $X=QR%$ decomposition and setting $X'=Q$. Say we computed the step (e.g. by L-BFGS) $\Delta X$ and we want to advance $X$. Following the geodesic along $\Delta X$ on the Grassmannian with a displacement parameter $\delta $ can be computed using

$$
X(\delta)=XV\cos(\Sigma)V^T + P \sin(\Sigma)V^T
$$

where $\Delta X  = P\Sigma V^T$ is the the compact singular value decomposition of $\Delta X$ with $\Sigma=\text{diag}(\sigma_1,\dots,\sigma_p)$ and things like $\cos(\Sigma)$ means $\text{diag}(\cos(\sigma_1),\dots,\cos(\sigma_p))$.

### Parallel Transport

**Singular Value Decomposition**

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
