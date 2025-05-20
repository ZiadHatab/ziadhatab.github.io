---
title: Jacobian of Eigendecomposition
date: 2025-03-22 12:00:00 +0100
categories: [Tutorial]
tags: [linear-algebra, decomposition, derivative]     # TAG names should always be lowercase
math: true
img_path: ../../../assets/img/posts_img/
image: # ...
---

If you ever get the chance to meet me and bring up matrix decomposition as a topic of discussion, it would be an understatement to say I love matrix decompositions, and the most beautiful one is eigendecomposition. While I've never preferred calculus as much as linear algebra, I find myself constantly amazed when calculus merges with linear algebra in an elegant way, particularly with matrix decompositions.

In this post, I want to discuss computing derivatives of eigendecomposition, specifically the Jacobian. Although this is a well-known topic, I believe sharing my perspective might offer valuable insights. For further reading, Google and Wikipedia provide an excellent start. For even further reading, see references [1]-[3] below.

## The Eigendecomposition
The eigendecomposition expresses a square matrix as a transformation operation on a diagonal matrix. I won't delve deeper into the interpretation of eigendecomposition and spectral analysis of matrices; for interested readers, the [Wikipedia page](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix) offers an excellent introduction with many references.

For a square matrix $\bs{X}\in\mathbb{C}^{N\times N}$, if it's diagonalizable (non-defective), its eigendecomposition can be written as:
\begin{equation}
\bs{X} = \bs{U}\bs{\Lambda}\bs{U}^{-1} = \begin{bmatrix} \vert & \vert & & \vert \\\ \bs{u}_1 & \bs{u}_2 & \cdots & \bs{u}_N \\\ \vert & \vert & & \vert \end{bmatrix}\begin{bmatrix} \lambda_1 & 0 & \cdots & 0 \\\ 0 & \lambda_2 & \cdots & 0 \\\ \vdots & \vdots & \ddots & \vdots \\\ 0 & 0 & \cdots & \lambda_N \end{bmatrix}\begin{bmatrix} \vert & \vert & & \vert \\\ \bs{u}_1 & \bs{u}_2 & \cdots & \bs{u}_N \\\ \vert & \vert & & \vert \end{bmatrix}^{-1}
\label{eq:1}
\end{equation}
where $\bs{U}$ is a square matrix containing the eigenvectors, and $\bs{\Lambda}$ is a diagonal matrix containing the eigenvalues.

Alternatively, Equation \eqref{eq:1} can be expressed by examining individual eigenvectors and their corresponding eigenvalues. This gives us the standard eigenvector/eigenvalue notation, which is valid for non-repeated eigenvalues:
\begin{equation}
\bs{X}\bs{u}_i = \lambda_i\bs{u}_i
\label{eq:2}
\end{equation}

An important consideration is that eigendecomposition exists only if we can find $N$ linearly independent eigenvectors. If not, the matrix is [defective](https://en.wikipedia.org/wiki/Defective_matrix) and can only be "diagonalized" using a [Jordan form](https://en.wikipedia.org/wiki/Jordan_matrix). While defective matrices always have repeated eigenvalues ($\lambda_i=\lambda_j$ for some $i\neq j$), the converse isn't necessarily true. For example, the identity matrix has identical eigenvalues but remains diagonalizable due to its independent eigenvectors. In this post, I will focus on computing the Jacobian of diagonalizable matrices, including those with repeated eigenvalues, while excluding defective matrices from the discussion.

## Uniqueness of Eigendecomposition
With any matrix decomposition, the question of "uniqueness" naturally arises. For eigendecomposition, while the eigenvalues are always unique for a given matrix, the eigenvectors are not. The degree of non-uniqueness in eigenvectors depends primarily on whether there are repeated eigenvalues.

Consider a $4\times4$ matrix with non-repeated eigenvalues. One valid eigendecomposition would be:
\begin{equation}
\bs{X} = \bs{U}\bs{\Lambda}\bs{U}^{-1} = \begin{bmatrix} \vert & \vert & \vert & \vert \\\ \bs{u}_1 & \bs{u}_2 & \bs{u}_3 & \bs{u}_4 \\\ \vert & \vert & \vert & \vert \end{bmatrix}\begin{bmatrix} \lambda_1 & 0 & 0 & 0 \\\ 0 & \lambda_2 & 0 & 0 \\\ 0 & 0 & \lambda_3 & 0 \\\ 0 & 0 & 0 & \lambda_4 \end{bmatrix}\begin{bmatrix} \vert & \vert & \vert & \vert \\\ \bs{u}_1 & \bs{u}_2 & \bs{u}_3 & \bs{u}_4 \\\ \vert & \vert & \vert & \vert \end{bmatrix}^{-1}
\label{eq:3}
\end{equation}

Another valid eigendecomposition using the same eigenvectors can be obtained by scaling with any non-zero diagonal matrix $\bs{D}$:
\begin{equation}
\bs{X} = (\bs{U}\bs{D})\bs{\Lambda}(\bs{U}\bs{D})^{-1} = \bs{U}\underbrace{\bs{D}\bs{\Lambda}\bs{D}^{-1}}_{=\bs{D}\bs{D}^{-1}\bs{\Lambda}}\bs{U}^{-1} = \bs{U}\bs{\Lambda}\bs{U}^{-1}
\label{eq:4}
\end{equation}

For non-repeated eigenvalues, each eigenvector is unique only up to a scalar multiple. Most software libraries normalize eigenvectors to have unit length. However, this normalization still leaves an ambiguous phase factor in the eigenvectors.

The complexity increases significantly when dealing with repeated eigenvalues. While eigenvalues maintain their uniqueness, eigenvectors become more ambiguous. Instead of just a scalar multiplier, the ambiguity involves an $M\times M$ matrix for each distinct eigenvalue with multiplicity $M$ (also known as the [algebraic multiplicity](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors#Algebraic_multiplicity)).

For example, consider a $4\times4$ matrix with one eigenvalue repeated twice:
\begin{equation}
\bs{X} = \bs{U}\bs{\Lambda}\bs{U}^{-1} = \begin{bmatrix} \vert & \vert & \vert & \vert \\\ \bs{u}_1 & \bs{u}_2 & \bs{u}_3 & \bs{u}_4 \\\ \vert & \vert & \vert & \vert \end{bmatrix}\begin{bmatrix} \lambda_1 & 0 & 0 & 0 \\\ 0 & \lambda_1 & 0 & 0 \\\ 0 & 0 & \lambda_3 & 0 \\\ 0 & 0 & 0 & \lambda_4 \end{bmatrix}\begin{bmatrix} \vert & \vert & \vert & \vert \\\ \bs{u}_1 & \bs{u}_2 & \bs{u}_3 & \bs{u}_4 \\\ \vert & \vert & \vert & \vert \end{bmatrix}^{-1}
\label{eq:5}
\end{equation}
where $\lambda_2=\lambda_1$ and the eigenvectors remain linearly independent.

For this case, we can find another decomposition by multiplying the eigenvector matrix by a block diagonal matrix $\bs{D}$ of the following form:
\begin{equation}
\bs{D} = \begin{bmatrix} \bs{G} & \bs{0} & \bs{0} \\\ \bs{0}^T & d_3 & 0 \\\ \bs{0}^T & 0 & d_4 \end{bmatrix}
\label{eq:6}
\end{equation}
where $\bs{G}\in\mathbb{C}^{2\times 2}$ is invertible.

This leads to the same eigendecomposition:
\begin{equation}
\bs{X} = (\bs{U}\bs{D})\bs{\Lambda}(\bs{U}\bs{D})^{-1} = \bs{U}\underbrace{\bs{D}\bs{\Lambda}\bs{D}^{-1}}_{=\bs{D}\bs{D}^{-1}\bs{\Lambda}}\bs{U}^{-1} = \bs{U}\bs{\Lambda}\bs{U}^{-1}
\label{eq:7}
\end{equation}

The non-repeated eigenvalue case is thus a special case where the block diagonal matrix reduces to a standard diagonal matrix. For a more conventional representation, we can express the eigendecomposition using both eigenvectors arranged in a matrix for their common eigenvalue:
\begin{equation}
\bs{X}\begin{bmatrix} \vert & \vert \\\ \bs{u}_1 & \bs{u}_2 \\\ \vert & \vert \end{bmatrix} = \begin{bmatrix} \vert & \vert \\\ \bs{u}_1 & \bs{u}_2 \\\ \vert & \vert \end{bmatrix}\lambda_1
\label{eq:8}
\end{equation}
where multiplying by any invertible $2\times 2$ matrix from the right side maintains the equality.

Understanding this non-uniqueness is essential when computing Jacobians or performing other mathematical operations on eigenvectors. Poorly chosen normalization methods (particularly those involving very small numbers) can result in numerically unstable solutions.

## The Jacobian and Derivative Rules

So, before I start discussing the Jacobian of eigendecomposition, it makes sense that I explain how Jacobians are generally computed and the derivative rules that I will be using along the way.

Starting basically, for simple single-input single-output functions, $f(x): \mathbb{R}\rightarrow\mathbb{C}$, the derivative is straightforward as $\mathrm{d}f/\mathrm{d}x$. Note, I'm intentionally computing derivatives with respect to real-valued parameters (the output can be complex). In general, this notation would work with complex-valued input parameters only if the function is analytic (sometimes called [holomorphic](https://en.wikipedia.org/wiki/Holomorphic_function)). When the analytic condition fails, then you have to split the complex-valued parameters into their real-valued basis and compute over $\mathbb{R}^2$. So, I will stick with input parameters being real-valued, and if you encounter non-analytical complex-valued input parameters in some problems, just split them and everything should work again.

Now, if I expand the dimension of the input parameters to be a vector, while the output is still a scalar, i.e., $f(\bs{x}): \mathbb{R}^N\rightarrow\mathbb{C}$, I can misuse the derivative notation and write:
\begin{equation}
\frac{\mathrm{d}f}{\mathrm{d}\bs{x}} = \begin{bmatrix} \frac{\partial f}{\partial x_1} & \frac{\partial f}{\partial x_2} & \cdots & \frac{\partial f}{\partial x_N} \end{bmatrix}
\label{eq:9}
\end{equation}

You might be looking at Equation \eqref{eq:9} and thinking "you cannot write $\mathrm{d}f/\mathrm{d}\bs{x}$ when $\bs{x}$ is a vector!" - and you would be right! However, arranging the partial derivatives into a vector is actually the correct definition, which is what we call the [Jacobian matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) (for a single-output function, this is the [Gradient](https://en.wikipedia.org/wiki/Gradient)). Therefore, I will use the following notation:
\begin{equation}
\bs{J}_f(\bs{x}) = \begin{bmatrix} \frac{\partial f}{\partial x_1} & \frac{\partial f}{\partial x_2} & \cdots & \frac{\partial f}{\partial x_N} \end{bmatrix}
\label{eq:10}
\end{equation}
where $\bs{J}_f(\bs{x})$ reads "the Jacobian of $f$ with respect to the vector $\bs{x}$".

As a matter of fact, with the Jacobian notation in Equation \eqref{eq:10}, we can easily expand this to account for multiple-output functions, i.e., vector-in, vector-out $\bs{f}(\bs{x}): \mathbb{R}^N\rightarrow\mathbb{C}^M$. Then the Jacobian indeed becomes a matrix as below:
\begin{equation}
\bs{J}_\bs{f}(\bs{x}) = \begin{bmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_N} \\\ \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_N} \\\ \vdots & \vdots & \ddots & \vdots \\\ \frac{\partial f_M}{\partial x_1} & \frac{\partial f_M}{\partial x_2} & \cdots & \frac{\partial f_M}{\partial x_N}
\end{bmatrix}
\label{eq:11}
\end{equation}

Thus, for $M$ functions dependent on $N$ inputs, the Jacobian has dimensions of $M\times N$. For the special case where $M=1$, we have what is known as the [Gradient](https://en.wikipedia.org/wiki/Gradient). Further, when both $M=1$ and $N=1$, this reduces to a [standard derivative](https://en.wikipedia.org/wiki/Derivative). The Jacobian notation elegantly generalizes derivatives to higher dimensions.

In fact, all conventional derivative rules apply to Jacobians. The chain rule is perhaps the most widely used. For example, consider the composition of two vector functions $\bs{f}$ and $\bs{g}$, where $\bs{f}(\bs{g}(\bs{x})): \mathbb{R}^N\rightarrow\mathbb{C}^M$. The Jacobian of this composition can be written as:
\begin{equation}
\bs{J}\_\bs{f}(\bs{x}) = \bs{J}\_\bs{f}(\bs{g})\bs{J}_\bs{g}(\bs{x})
\label{eq:12}
\end{equation}

Furthermore, the product rule holds for Jacobians as well. However, we need to be careful here because we've only defined Jacobians for vector-functions and vector inputs. So far, the only functional product we can use is the inner product. I'll show you how to resolve this using some powerful linear algebra tools, specifically the Kronecker product. For the simple inner-product case (or multiplication with a scaling function), if we have two vector functions $\bs{f}$ and $\bs{g}$, the Jacobian of their product is given by:
\begin{equation}
\bs{J}\_{\bs{f}\cdot\bs{g}}(\bs{x}) =  \bs{g}(\bs{x}) \cdot \bs{J}_{\bs{f}}(\bs{x}) + \bs{f}(\bs{x})\cdot \bs{J}\_{\bs{g}}(\bs{x})
\label{eq:13}
\end{equation}

For more details, I recommend reviewing the Wikipedia pages on [Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant), [chain rule](https://en.wikipedia.org/wiki/Chain_rule), and [product rule](https://en.wikipedia.org/wiki/Product_rule).

To handle matrix-functions and matrix inputs, we need a method to convert them into vectors. For this purpose, I'll introduce a powerful tool that might be unfamiliar to many: the [Kronecker product](https://en.wikipedia.org/wiki/Kronecker_product) combined with the [vectorization operator](https://en.wikipedia.org/wiki/Vectorization_(mathematics)). While I won't explain their element-wise operations here, you can find comprehensive details in their respective Wikipedia pages or references [4]-[5] below.

For Jacobian calculations, we need a systematic way to vectorize matrices and their products so that standard Jacobian properties still apply. The Kronecker product provides an elegant solution, particularly when dealing with the product of three matrices $\bs{X}\bs{Y}\bs{Z}$. Their vectorized form is given by:
\begin{equation}
\vc{\bs{X}\bs{Y}\bs{Z}} = (\bs{Z}^T\otimes\bs{X})\vc{\bs{Y}}
\label{eq:14}
\end{equation}

While Equation \eqref{eq:14} showed how vectorization can be expressed as a matrix-vector product with $\vc{\bs{Y}}$ as the vector, we can also vectorize the other matrices by treating them as the central matrix and treating their neighbors as identity matrices:
\begin{equation}
\vc{\bs{X}\bs{Y}\bs{Z}} = (\bs{Z}^T\otimes\bs{X})\vc{\bs{Y}} = (\bs{I}^T\otimes\bs{X}\bs{Y})\vc{\bs{Z}} = ((\bs{Y}\bs{Z})^T\otimes\bs{I})\vc{\bs{X}}
\label{eq:15}
\end{equation}

> __Bonus fact__: The [Khatri–Rao product](https://en.wikipedia.org/wiki/Khatri%E2%80%93Rao_product) is similar in properties to Kronecker product, but strictly for vectorizing diagonal matrices, e.g.,
\begin{equation}
\vc{\bs{X}\bs{D}\_\bs{y}\bs{Z}} = (\bs{Z}^T\ast\bs{X})\bs{y}
\label{eq:16}
\end{equation}
where $\bs{D}\_\bs{y} = \mathrm{diag}(\bs{y})$ is a diagonal matrix made from the vector $\bs{y}$.
{: .prompt-info }

Using these properties, we can readily expand the Jacobian for matrix functions and their products. Consider two matrix functions $\bs{F}(\bs{X})$ and $\bs{G}(\bs{X})$ that take a matrix input $\bs{X}$. The Jacobian of their product is defined as the Jacobian of their vectorized versions:
\begin{equation}
\bs{J}_{\vc{\bs{F}\bs{G}}}(\vc{\bs{X}}) = (\bs{G}^T(\bs{X})\otimes \bs{I})\bs{J}\_{\vc{\bs{F}}}(\vc{\bs{X}}) + (\bs{I}\otimes \bs{F}(\bs{X}))\bs{J}\_{\vc{\bs{G}}}(\vc{\bs{X}})
\label{eq:17}
\end{equation}

While Equation \eqref{eq:17} may appear complex due to the explicit vectorization notation $\vc{}$, it is actually quite compact. It follows the same pattern as the general vector case shown in Equation \eqref{eq:13}. When computing the partial Jacobian for $\bs{F}$, we treat $\bs{G}$ as constant and use the Kronecker product to vectorize $\bs{F}$. Similarly, when computing the partial Jacobian for $\bs{G}$, we treat $\bs{F}$ as constant and use the Kronecker product to vectorize $\bs{G}$.

For a deeper understanding of matrix derivatives and Jacobians, I recommend the book "Complex-Valued Matrix Derivatives" [6].

## Jacobian of the Eigendecomposition with Non-repeated Eigenvalues

Now we can address the main problem: computing the Jacobian of eigenvalues and eigenvectors. While we could formulate this problem using matrix notation as $\bs{X} = \bs{U}\bs{\Lambda}\bs{U}^{-1}$, it's simpler to start with vector notation. Later, we can combine individual Jacobians to describe all eigenvalues and eigenvectors in a unified form.

Let's begin with the fundamental eigendecomposition equation for a distinct eigenvalue:
\begin{equation}
\bs{X}\bs{u}_i = \bs{u}_i\lambda_i
\label{eq:18}
\end{equation}

Additionally, we need the left eigendecomposition of the same matrix, which we write in transposed form:
\begin{equation}
\bs{X}^T\bs{v}_i = \bs{v}_i\lambda_i \Longrightarrow \bs{v}_i^T\bs{X} = \lambda_i\bs{v}_i^T
\label{eq:19}
\end{equation}
where the left eigenvector matrix $\bs{V}$ relates to the right eigenvector matrix $\bs{U}$ by $\bs{V} = \bs{U}^{-T}$. The need for left eigendecomposition in Equation \eqref{eq:19} will become clear shortly.

Let's start by computing the Jacobian of Equation \eqref{eq:18} and applying the product rule. While I'll omit the independent parameter for brevity, remember that we're taking derivatives with respect to a set of independent parameters:
\begin{equation}
(1\otimes\bs{X})\bs{J}\_{\bs{u}\_i} + (\bs{u}\_i^T\otimes\bs{I})\bs{J}_{\bs{X}} =  (\lambda_i\otimes\bs{I})\bs{J}\_{\bs{u}\_i} + (1\otimes\bs{u}_i)\bs{J}\_{\lambda_i}
\label{eq:20}
\end{equation}

Next, we leverage the left eigendecomposition from Equation \eqref{eq:19}. We multiply both sides of Equation \eqref{eq:20} with $(1\otimes\bs{v}\_i^T)$ and use the Kronecker product property $(\bs{A}\otimes\bs{B})(\bs{C}\otimes\bs{D}) = (\bs{AC})\otimes(\bs{BD})$:
\begin{equation}
(1\otimes\underbrace{\bs{v}\_i^T\bs{X}}\_{=\lambda_i\bs{v}\_i^T})\bs{J}\_{\bs{u}\_i} + (\bs{u}\_i^T\otimes\bs{v}\_i^T)\bs{J}_{\bs{X}} =  (\underbrace{\lambda_i\otimes\bs{v}\_i^T}\_{=1\otimes\lambda_i\bs{v}\_i^T})\bs{J}\_{\bs{u}\_i} + (\underbrace{1\otimes\bs{v}\_i^T\bs{u}_i}\_{=\bs{v}\_i^T\bs{u}_i})\bs{J}\_{\lambda_i}
\label{eq:21}
\end{equation}

In Equation \eqref{eq:21}, I made several simplifications. First, from Equation \eqref{eq:19}, we know that $\bs{v}\_i^T\bs{X} = \lambda_i\bs{v}\_i^T$. Second, since $\lambda_i$ is a scalar, we can commute it in the Kronecker product: $\lambda_i\otimes\bs{v}\_i^T = 1\otimes\lambda_i\bs{v}\_i^T$. Third, the last Kronecker product simplifies to a scalar as both terms are scalars: $1\otimes\bs{v}\_i^T\bs{u}_i = \bs{v}\_i^T\bs{u}_i$. After canceling $(1\otimes\lambda_i\bs{v}\_i^T)\bs{J}\_{\bs{u}\_i}$ from both sides and rearranging to solve for $\bs{J}\_{\lambda_i}$, we get:

\begin{equation}
\boxed{\bs{J}\_{\lambda_i} = \frac{1}{\bs{v}\_i^T\bs{u}\_i}(\bs{u}\_i^T\otimes\bs{v}\_i^T)\bs{J}_\bs{X}}
\label{eq:22}
\end{equation}

Equation \eqref{eq:22} gives us the Jacobian of the $i$th eigenvalue, which depends only on the Jacobian of the matrix itself and the $i$th eigenvector. To get the Jacobian of all eigenvalues, simply stack all individual Jacobians. However, remember that this solution assumes non-repeated eigenvalues. For cases with repeated eigenvalues, see the following section.

Now, let's discuss the Jacobian of the eigenvectors. This derivation is straightforward since we already have the Jacobian of the eigenvalue. By substituting the Jacobian of the eigenvalue into Equation \eqref{eq:20} and collecting common terms, we obtain:
\begin{equation}
(\underbrace{1\otimes\bs{X}}\_{=\bs{X}})\bs{J}\_{\bs{u}\_i} + (\bs{u}\_i^T\otimes\bs{I})\bs{J}\_{\bs{X}} =  (\underbrace{\lambda_i\otimes\bs{I}}\_{=\lambda_i\bs{I}})\bs{J}\_{\bs{u}\_i} + \frac{1}{\bs{v}\_i^T\bs{u}\_i}\underbrace{(1\otimes\bs{u}\_i)(\bs{u}\_i^T\otimes\bs{v}\_i^T)}\_{=(\bs{u}\_i^T\otimes\bs{u}_i\bs{v}\_i^T)}\bs{J}\_\bs{X}
\label{eq:23}
\end{equation}

Rearranging terms with $\bs{J}\_{\bs{u}\_i}$ on the left side yields:
\begin{equation}
\boxed{(\bs{X} - \lambda_i\bs{I})\bs{J}\_{\bs{u}\_i} = \left(\frac{1}{\bs{v}\_i^T\bs{u}\_i}(\bs{u}\_i^T\otimes\bs{u}_i\bs{v}\_i^T) - (\bs{u}\_i^T\otimes\bs{I})\right)\bs{J}\_{\bs{X}}}
\label{eq:24}
\end{equation}

While Equation \eqref{eq:24} isn't the final answer for $\bs{J}\_{\bs{u}\_i}$, this is as far as we can go without additional assumptions. Until now, we never made any assumptions besides dealing with non-defective matrices and having non-repeating eigenvalues. The question of uniqueness becomes relevant at this point.

To solve for $\bs{J}\_{\bs{u}\_i}$, we might consider using the pseudo-inverse of $(\bs{X} - \lambda_i\bs{I})$. However, since we're dealing with non-repeated eigenvalues, $(\bs{X} - \lambda_i\bs{I})$ has rank $N-1$, making it singular. While a pseudo-inverse exists, it omits the nullspace component of the solution. For full-rank matrices, the nullspace contains only the zero vector, but here, $(\bs{X} - \lambda_i\bs{I})$ has one non-zero nullspace vector due to its rank deficiency.

This rank deficiency occurs because eigenvectors for non-repeated eigenvalues are only unique up to scalar multiplication. To resolve this ambiguity, we must determine this scalar. One approach is to norm normalize the eigenvectors, though this leaves a phase ambiguity. Another method is to fix one component of the eigenvector, making its Jacobian zero at that position. The key insight is that we cannot obtain a unique solution for the eigenvector Jacobian without some form of anchoring constraint.

The complete solution for $\bs{J}\_{\bs{u}\_i}$ combines the particular solution using the pseudo-inverse with a scaled nullspace solution:
\begin{equation}
\boxed{ \bs{J}\_{\bs{u}\_i} = (\bs{X} - \lambda_i\bs{I})^{+}\left(\frac{1}{\bs{v}\_i^T\bs{u}\_i}(\bs{u}\_i^T\otimes\bs{u}_i\bs{v}\_i^T) - (\bs{u}\_i^T\otimes\bs{I})\right)\bs{J}\_{\bs{X}} + \alpha (\bs{1}^T\otimes\bs{q}) }
\label{eq:25}
\end{equation}
where $\bs{q}$ is the non-trivial nullspace vector satisfying $(\bs{X} - \lambda_i\bs{I})\bs{q} = \bs{0}$, $\alpha$ is any non-zero scalar, and $\bs{1}^T = [1, 1, \cdots, 1]$ is a vector of ones that repeats the nullspace solution to match the dimensions of $\bs{J}\_{\bs{u}\_i}$. The term $\alpha (\bs{1}^T\otimes\bs{q})$ has the structure:
\begin{equation}
\alpha (\bs{1}^T\otimes\bs{q}) = \alpha \begin{bmatrix} \vert & \vert &  & \vert \\\ \bs{q} & \bs{q} & \cdots & \bs{q} \\\ \vert & \vert &  & \vert \end{bmatrix}
\label{eq:26}
\end{equation}

## The Problem with Repeated Eigenvalues

For cases with repeated eigenvalues, we need a special handling approach. While the Jacobian of eigenvectors remains similar to the previous section (as we still assume non-defective matrices), repeated eigenvalues require different treatment since multiple linearly independent eigenvectors are associated with a common repeated eigenvalue. The derivation starts by considering all associated eigenvectors in the eigendecomposition for the repeated eigenvalue:
\begin{equation}
\bs{X}\underbrace{\begin{bmatrix} \vert & \vert & & \vert \\\ \bs{u}_1 & \bs{u}_2 & \cdots & \bs{u}_M \\\ \vert & \vert & & \vert\end{bmatrix}}\_{\bs{U}_i} = \underbrace{\begin{bmatrix} \vert & \vert & & \vert \\\ \bs{u}_1 & \bs{u}_2 & \cdots & \bs{u}_M \\\ \vert & \vert & & \vert\end{bmatrix}}\_{\bs{U}_i}\lambda_i
\label{eq:27}
\end{equation}
where $M$ is the algebraic multiplicity of the $i$th eigenvalue.

Similarly, we can write the left eigendecomposition as:
\begin{equation}
\bs{X}^T\bs{V}_i = \bs{V}_i\lambda_i \Longrightarrow \bs{V}_i^T\bs{X} = \lambda_i\bs{V}_i^T
\label{eq:28}
\end{equation}
where $\bs{V}_i$ represents a subset of the full left eigenvector matrix corresponding to the $i$th eigenvalue, with the full left eigenvector given by $\bs{V} = \bs{U}^{-T}$.

Following the same approach as in Equation \eqref{eq:20}, we compute the Jacobian of both sides of Equation \eqref{eq:27} and apply the product rule:
\begin{equation}
(\bs{I}\_M\otimes\bs{X})\bs{J}\_{\bs{U}\_i} + (\bs{U}\_i^T\otimes\bs{I}\_N)\bs{J}_{\bs{X}} = (\lambda_i\bs{I}\_M\otimes\bs{I}_N)\bs{J}\_{\bs{U}\_i} + (\bs{I}_M\ast\bs{U}_i)\bs{1}\_M\bs{J}\_{\lambda_i}
\label{eq:29}
\end{equation}
where $\bs{I}_M$ and $\bs{I}_N$ are identity matrices of size $M$ and $N$ respectively, $\bs{1}\_M$ is a vector of ones with length $M$, and $\ast$ denotes the [Khatri–Rao product](https://en.wikipedia.org/wiki/Khatri%E2%80%93Rao_product).

As in the non-repeating case, we multiply both sides from the left with $(\bs{I}\_M\otimes\bs{V}\_i^T)$:
\begin{equation}
(\bs{I}\_M\otimes\underbrace{\bs{V}\_i^T\bs{X}}\_{=\lambda_i\bs{V}\_i^T})\bs{J}\_{\bs{U}\_i} + (\bs{U}\_i^T\otimes\bs{V}\_i^T)\bs{J}_{\bs{X}} = (\underbrace{\lambda_i\bs{I}\_M\otimes\bs{V}\_i^T}\_{=\bs{I}\_M\otimes\lambda_i\bs{V}\_i^T})\bs{J}\_{\bs{U}\_i} + \underbrace{(\bs{I}\_M\otimes\bs{V}\_i^T)(\bs{I}\_M\ast\bs{U}_i)}\_{=(\bs{I}\_M\ast\bs{V}\_i^T\bs{U}_i)}\bs{1}\_M\bs{J}\_{\lambda_i}
\label{eq:30}
\end{equation}

In Equation \eqref{eq:30}, we apply the same simplifications as before, plus the Khatri–Rao and Kronecker product property: $(\bs{A}\otimes\bs{B})(\bs{C}\ast\bs{D}) = (\bs{A}\bs{C})\ast(\bs{B}\bs{D})$.

After canceling the terms associated with $\bs{J}\_{\bs{U}\_i}$ and rearranging the remaining terms with $\bs{J}_{\bs{X}}$ and $\bs{J}\_{\lambda_i}$, the Jacobian of the eigenvalue becomes:
\begin{equation}
\boxed{\bs{J}\_{\lambda_i} = \left((\bs{I}\_M\ast\bs{V}\_i^T\bs{U}\_i)\bs{1}\_M\right)^{+}(\bs{U}\_i^T\otimes\bs{V}\_i^T)\bs{J}\_{\bs{X}}}
\label{eq:31}
\end{equation}

Note that Equation \eqref{eq:31} reduces to Equation \eqref{eq:22} when $M=1$. The pseudo-inverse in Equation \eqref{eq:31} will always exist since it solves for a rank-1 term, and the matrix will always be at least rank-1 (exactly rank-1 when $M=1$).

For the eigenvector Jacobian, since the eigenvectors are linearly independent (remember, non-defective matrix!), we can solve for each eigenvector individually rather than all $M$ simultaneously. We reuse Equation \eqref{eq:20} and substitute the eigenvalue Jacobian from Equation \eqref{eq:31}.

Starting with Equation \eqref{eq:20} and denoting eigenvectors related to repeated eigenvalues with index $k=[1,\ldots,M]$:
\begin{equation}
(\underbrace{1\otimes\bs{X}}\_{=\bs{X}})\bs{J}\_{\bs{u}\_{i,k}} + (\bs{u}\_{i,k}^T\otimes\bs{I})\bs{J}_{\bs{X}} =  (\underbrace{\lambda_i\otimes\bs{I}}\_{=\lambda_i\bs{I}})\bs{J}\_{\bs{u}\_{i,k}} + (\underbrace{1\otimes\bs{u}\_{i,k}}\_{=\bs{u}\_{i,k}})\bs{J}\_{\lambda_i}
\label{eq:32}
\end{equation}

After substituting $\bs{J}\_{\lambda_i}$ from Equation \eqref{eq:31} and simplifying:
\begin{equation}
\boxed{(\bs{X}-\lambda_i\bs{I})\bs{J}\_{\bs{u}\_{i,k}} = \left(\bs{u}\_{i,k}\left((\bs{I}\_M\ast\bs{V}\_i^T\bs{U}\_i)\bs{1}\_M\right)^{+}(\bs{U}\_i^T\otimes\bs{V}\_i^T) - (\bs{u}\_{i,k}^T\otimes\bs{I})\right)\bs{J}_{\bs{X}}}
\label{eq:33}
\end{equation}

As before, solving for $\bs{J}\_{\bs{u}\_{i,k}}$ requires assumptions about the eigenvector. With $M$ degrees of freedom, we need to anchor $M$ values in the eigenvector for a unique solution.

The complete solution for $\bs{J}\_{\bs{u}\_{i,k}}$ is:
\begin{equation}
\boxed{\bs{J}\_{\bs{u}\_{i,k}} = (\bs{X}-\lambda_i\bs{I})^{+}\left(\bs{u}\_{i,k}\left((\bs{I}\_M\ast\bs{V}\_i^T\bs{U}\_i)\bs{1}\_M\right)^{+}(\bs{U}\_i^T\otimes\bs{V}\_i^T) - (\bs{u}\_{i,k}^T\otimes\bs{I})\right)\bs{J}_{\bs{X}} + \bs{1}^T\otimes\bs{Q}\bs{\alpha}}
\label{eq:34}
\end{equation}
where $\bs{Q}$ contains all $M$ non-trivial nullspace solutions of $(\bs{X} - \lambda_i\bs{I})$, i.e., $(\bs{X} - \lambda_i\bs{I})\bs{Q} = \bs{0}$, and $\bs{\alpha}$ contains the arbitrary non-zero scalars needed for a unique solution. The term $\bs{1}^T\otimes\bs{Q}\bs{\alpha}$ has the structure:
\begin{equation}
\bs{1}^T\otimes\bs{Q}\bs{\alpha} = \begin{bmatrix} \vert & \vert &  & \vert \\\ \bs{Q}\bs{\alpha} & \bs{Q}\bs{\alpha} & \cdots & \bs{Q}\bs{\alpha} \\\ \vert & \vert &  & \vert \end{bmatrix}
\label{eq:35}
\end{equation}
where $\bs{Q}$ and $\bs{\alpha}$ are defined as:
\begin{equation}
\bs{Q} = \begin{bmatrix} \vert & \vert &  & \vert \\\ \bs{q}_1 & \bs{q}_2 & \cdots & \bs{q}_M \\\ \vert & \vert &  & \vert \end{bmatrix}; \qquad \bs{\alpha} = \begin{bmatrix} \alpha_1 \\\ \alpha_2 \\\ \vdots \\\ \alpha_M \end{bmatrix}
\label{eq:36}
\end{equation}

## Python Implementation
Below is a sample code implementation for computing the Jacobian of eigendecomposition for non-repeated eigenvalues. Note that the Jacobian of each eigenvector includes a scalar ambiguity. For cases involving repeated eigenvalues, I'll leave the implementation as an exercise for interested readers 😉.

```python
import numpy as np

def jaceig(X, JX=None, return_eig=False):
    """
    Computes the Jacobian of eigenvalues and eigenvectors with respect to matrix perturbations.

    Parameters
    ----------
    X : array_like
        Input matrix for which to compute eigendecomposition Jacobians.
        
    JX : array_like, optional
        Jacobian of the input matrix X. If None, assumes identity matrix of size N*N where N
        is the dimension of X. Default is None.
        
    return_eig : bool, optional
        If True, also returns the eigenvalues and eigenvectors of X.
        Default is False.

    Returns
    -------
    Jlambs : ndarray
        Jacobian matrix of all eigenvalues with respect to matrix perturbations.
        Shape is (N, N*N) where N is dimension of input matrix.
        
    JU : ndarray
        Jacobian matrix of all eigenvectors with respect to matrix perturbations.
        Shape is (N*N, N*N).
        
    lambs : ndarray, optional
        Eigenvalues of input matrix X. Only returned if return_eig=True.
        
    U : ndarray, optional
        Right eigenvectors of input matrix X. Only returned if return_eig=True.
    """
    X = np.atleast_2d(X)
    N = X.shape[0]
    I = np.eye(N)
    JX = JX if JX is not None else np.eye(N*N)

    # compute the eigendecomposition
    lambs, U = np.linalg.eig(X) # right eigenvectors
    V = np.linalg.inv(U).T      # left eigenvectors

    # iterate through each eigenvalue/eigenvector
    Jlambs = []
    JU = []
    for inx in range(N):
        uiT = U[:,inx]  # already transposed 
        viT = V[:,inx]  # already transposed
        # Jacobian of eigenvalue
        norm = viT.dot(uiT)
        Jlambs.append((np.kron(uiT,viT)/norm)@JX)
        # Jacobian of eigenvector
        P = np.linalg.pinv(X - lambs[inx]*I)
        JU.append(P@(np.kron(uiT,np.outer(uiT,viT))/norm - np.kron(uiT,I))@JX)
    
    if return_eig:
        return np.vstack(Jlambs), np.vstack(JU), lambs, U
    else:
        return np.vstack(Jlambs), np.vstack(JU)

if __name__ == '__main__':
    # Example with 4x4 matrix
    X = np.array([[ 4.0,  1.0, -2.0,  2.0],
                  [ 1.0,  3.0,  0.0,  1.0],
                  [-2.0,  0.0,  3.0,  2.0],
                  [ 2.0,  1.0,  2.0,  6.0]])
    
    # compute Jacobians
    Jlambs, JU = jaceig(X)
    
    # Verify results by comparing with finite differences.
    # Be careful with this! If eps is too small, you get rounding error.
    # If eps is too large, you might get different ordering
    # of eigenvalues and eigenvectors. Also, the Jacobian of the eigenvectors 
    # is not unique and will have a sign ambiguity.
    eps = 1e-8
    dX = np.random.randn(*X.shape)
    # compute eigendecomposition at X and X + eps*dX
    lambs1, U1 = np.linalg.eig(X)
    lambs2, U2 = np.linalg.eig(X + eps*dX)
    
    # finite difference solution
    dlambs_fd = (lambs2 - lambs1)/eps
    dU_fd = (U2 - U1).flatten('F')/eps
    
    # solution using the analytical Jacobian
    dlambs = Jlambs@dX.flatten('F')
    dU = JU@dX.flatten('F')

    # compare results
    print('\nFinite difference check:')
    print(f'Max error in eigenvalues = {np.max(np.abs(dlambs - dlambs_fd))}')
    print(f'Max error in eigenvectors = {np.max(np.abs(dU - dU_fd))}')

    # EOF
```

## References

[1] Magnus, Jan R. "On Differentiating Eigenvalues and Eigenvectors." Econometric Theory, vol. 1, no. 2, 1985, pp. 179–91. JSTOR, <http://www.jstor.org/stable/3532409>.

[2] J. H. Wilkinson, The Algebraic Eigenvalue Problem. Oxford: Clarendon Press, 1988. ISBN: 9780198534181.

[3] R. B. Nelson, "Simplified calculation of eigenvector derivatives," AIAA Journal, vol. 14, no. 9, pp. 1201–1205, Sep. 1976, doi: [10.2514/3.7211](https://doi.org/10.2514/3.7211).

[4] J. Brewer, "Kronecker products and matrix calculus in system theory," in IEEE Transactions on Circuits and Systems, vol. 25, no. 9, pp. 772-781, September 1978, doi: [10.1109/TCS.1978.1084534](https://doi.org/10.1109/TCS.1978.1084534).

[5] J. Brewer, "Correction to 'Kronecker Products and Matrix Calculus in System Theory'," in IEEE Transactions on Circuits and Systems, vol. 26, no. 5, pp. 360-360, May 1979, doi: [10.1109/TCS.1979.1084638](https://doi.org/10.1109/TCS.1979.1084638).

[6] A. Hjørungnes, Complex-Valued Matrix Derivatives: With Applications in Signal Processing and Communications. Cambridge: Cambridge University Press, 2011, doi: [10.1017/CBO9780511921490](https://doi.org/10.1017/CBO9780511921490).