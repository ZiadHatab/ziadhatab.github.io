---
title: A Primer on TRL Calibration
date: 2022-12-27 22:00:00 +0100
categories: [Tutorial]
tags: [vna, calibration]     # TAG names should always be lowercase
math: true
img_path: ../../../assets/img/posts_img/
image: TRL_waveguide_GCPW_kit.png
---

If you're looking to get accurate measurements from your device under test (DUT) using a vector network analyzer (VNA), calibrating your VNA is a crucial step. There are many methods out there, but I believe TRL (thru-reflect-line) calibration [1] is one of the best in microwave metrology. It's quite simple to design a TRL kit, and the kit can be manufactured accurately using standard machinery process. 

TRL calibration has some limitations though, such as being frequency limited and requiring the VNA to be described with the error-box model. However, the frequency limitation can be improved with multiline TRL, which I will cover in a future post. Most modern VNAs use a four-receiver architecture and can be accurately modeled by the error-box model in a straightforward manner.

In this post, I'll walk you through the mathematical development of TRL calibration. I also included a section to guide you on how to design your TRL kit yourself. By the end of this post, you'll be able to implement the calibration yourself and design your own TRL kit with confidence. If you're interested, my [github repository](https://github.com/ZiadHatab/trl-calibration) has a python implementation that follows the same notation used here.

Let's get started!

## Error-box modeling of a 2-port VNA

When measuring a 2-port DUT with a VNA, we can describe the overall measurement model in terms of cascading the DUT with unknown network boxes. We refer to such a model as the error-box model (see Fig. 1).

![Illustration of VNA 2-port error-box model](VNA_error_box_illustration.png)
_Fig. 1. Illustration of VNA 2-port error-box model._

The goal of the calibration is to estimate the error-boxes. With the obtained estimates, we can remove their effects by an inverse operation (de-cascading), therefore shifting the measurement plane to the DUT. There are, however, some caveats you need to be aware of when modeling a 2-port VNA with an error-box model:

1. We assume that the ports of the VNA are isolated from each other, i.e., leakage between the ports is negligible. This is generally true in most VNAs.

2. Most VNAs measure S-parameters under the assumption that the source termination at the receiving port is reflectionless. However, source termination is not ideal and we refer to the reflection it cause as _switch terms_. The term _switch_ is quite misleading, as we are dealing with reflection due to termination and from the switches themselves. To know the historical reason behind it, see [2]. In short, if you are dealing with a four receiver architecture, you can measure the forward and reverse switch-terms by connecting a line standard between the ports (e.g., straight connection) and measuring the wave ratios:
\begin{equation}
\Gamma_\mathrm{f} = \left.\frac{a_2}{b_2}\right|\_\text{sourced by port 1}, \quad \Gamma_\mathrm{r} = \left.\frac{a_1}{b_1}\right|_\text{sourced by port 2}\label{eq:1}
\end{equation}
The equations to correct the S-parameters can be found in references [2].

Now that we've cleared up the caveats about the error-box model, we can delve into the math. The best way to describe the error-box model is by the T-parameters (also known as the chain parameters). The VNA is now simplified to a simple cascade of networks, as shown in Fig. 2.

![Block diagram of the error-box model of a 2-port VNA](block_error_box_model.png)
_Fig. 2. Block diagram of the error-box model of a 2-port VNA._

The conversion between the S-parameters and T-parameters is give by

\begin{equation}
\bs{T} = \frac{1}{S_{21}}\begin{bmatrix} S_{12}S_{21}-S_{11}S_{22} & S_{11} \\\ -S_{22} & 1\end{bmatrix}, \qquad
\bs{S} = \frac{1}{T_{22}}\begin{bmatrix}T_{12} & T_{11}T_{22}-T_{12}T_{21} \\\ 1 & -T_{21}\end{bmatrix}\label{eq:2}
\end{equation}

With the help of the T-parameters definition, we can describe the error-box model in Fig. 2 by

\begin{equation}
\bs{M}\_\mathrm{dut} = \underbrace{k_a \begin{bmatrix} a_{11} & a_{12} \\\ a_{21} & 1\end{bmatrix} }\_{\text{left error-box}}\bs{T}\_\mathrm{dut}\underbrace{\begin{bmatrix} b_{11} & b_{12} \\\ b_{21} & 1\end{bmatrix}k_b}_{\text{right error-box}}\label{eq:3}
\end{equation}

The parameters $\\\{k_a, a_{11}, a_{21}, a_{12}, k_b, b_{11}, b_{21}, b_{12}\\\}$ are the unknowns that we want to solve for. You can notice that we don’t really need to solve for $k_a$ and $k_b$ separately, we only need their product, thus we can simplify the model by defining $k=k_ak_b$. Therefore, we only need to solve for 7-unknowns, and this simplifies Eq. \eqref{eq:3} to

\begin{equation}
\bs{M}\_\mathrm{dut} = \underbrace{k_ak_b}\_{k}\underbrace{\begin{bmatrix} a_{11} & a_{12} \\\ a_{21} & 1\end{bmatrix}}\_{\bs{A}}\bs{T}\_\mathrm{dut}\underbrace{\begin{bmatrix} b_{11} & b_{12} \\\ b_{21} & 1\end{bmatrix}}\_{\bs{B}} = k\bs{A}\bs{T}\_\mathrm{dut}\bs{B}\label{eq:4}
\end{equation}

## Solving the error-boxes with TRL standards

Now, that we have the model established, we shift our focused on the details of TRL calibration. As the acronym stands for Thru-Reflect-Line, we are working with three standards. The calibration is done in two steps. The first step we use the the line and the thru standard to form an eigenvalue problem to solve for normalized calibration coefficients. Then, in a second step we use the reflect standard and the thru standard again to undo the normalization.

### Eigenvalue problem

The start by measuring the thru and line standard, described in their T-parameters:

\begin{equation}
\bs{M}\_\mathrm{thru} = k\bs{A}\begin{bmatrix}1 & 0 \\\ 0 & 1\end{bmatrix}\bs{B}, \qquad
\bs{M}_\mathrm{line} = k\bs{A}\begin{bmatrix}e^{-l\gamma} & 0 \\\ 0 & e^{l\gamma}\end{bmatrix}\bs{B}\label{eq:5}
\end{equation}

where $l$ is the length of the line (referenced to the thru) and $\gamma$  is the propagation constant. Note that the thru standard is same as a line standard but with $l=0$. Next, we multiply the inverses of the thru measurement with the line measurement:

\begin{equation}
\bs{M}\_\mathrm{line}\bs{M}^{-1}\_\mathrm{thru} = \bs{A}\begin{bmatrix}e^{-l\gamma} & 0 \\\ 0 & e^{l\gamma}\end{bmatrix}\bs{A}^{-1}\label{eq:6}
\end{equation}

Eq. \eqref{eq:6} is the eigenvalue problem we are looking for, where $e^{-l\gamma}$ and $e^{l\gamma}$ are the eigenvalues and $\bs{A}$ holds the two eigenvectors. We can also write Eq. \eqref{eq:6} in the following way:

\begin{equation}
\bs{M}\_\mathrm{line}\bs{M}^{-1}\_\mathrm{thru}\begin{bmatrix} a_{11} \\\ a_{21}\end{bmatrix} = e^{-l\gamma}\begin{bmatrix} a_{11} \\\ a_{21}\end{bmatrix};\qquad \bs{M}\_\mathrm{line}\bs{M}^{-1}\_\mathrm{thru}\begin{bmatrix} a_{12} \\\ 1\end{bmatrix} = e^{l\gamma}\begin{bmatrix} a_{12}\\\ 1\end{bmatrix}\label{eq:7}
\end{equation}

This means, if apply the [eigendecomposition](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix) of $\bs{M}\_\mathrm{line}\bs{M}^{-1}\_\mathrm{thru}$ we should recover the matrix $\bs{A}$ directly, right?! Well, not completely, we have two issues to deal with:

1. When we apply the eigendecomposition, the order of the eigenvalues is arbitrary. From our point of view, we only get the two eigenvalues as numbers and we don't know which one represents $e^{-l\gamma}$ or $e^{l\gamma}$.

2. You should know from your linear algebra 101, that the eigenvectors are unique up to a scalar multiple. That means, we can only solve the coefficients in $\bs{A}$ as normalized coefficients. Think about it this way, look at Eq. \eqref{eq:7}, if we multiply the eigenvectors with a non-zero scalar, it will cancel out automatically, since the scalar exists also on the other side of the equation. That is why any scalar multiple of the eigenvector is still a valid eigenvector.

To resolve the issue with sorting the eigenvalues, we compare the eigenvalues with an estimate we know. That is, we already have an estimate of $\gamma$ to build $e^{-l\gamma}$ and $e^{l\gamma}$. This estimate doesn't need to be extremely accurate, but should be close enough to the right answer to allow us to distinguish the eigenvalues. Generally, the way this works is that we start with an estimate of $\gamma$ at the first frequency point and solve the eigenvalue problem, and then extract the new $\gamma$ from the eigenvalues. This new $\gamma$ is your new estimate to use in the next frequency point.

Now, the way we solve for the coefficients is by normalizing them such that one element of the vector equals one. The equations in Eq. \eqref{eq:7} become:

\begin{equation}
\bs{M}\_\mathrm{line}\bs{M}^{-1}\_\mathrm{thru}\begin{bmatrix} 1 \\\ a_{21}/a_{11}\end{bmatrix} = e^{-l\gamma}\begin{bmatrix} 1 \\\ a_{21}/a_{11}\end{bmatrix};\qquad \bs{M}\_\mathrm{line}\bs{M}^{-1}\_\mathrm{thru}\begin{bmatrix} a_{12} \\\ 1\end{bmatrix} = e^{l\gamma}\begin{bmatrix} a_{12}\\\ 1 \end{bmatrix} \label{eq:8}
\end{equation}

It's worth noting that $a_{12}$ doesn't require normalization since the eigenvector is already normalized. You may be wondering why I chose to divide the first eigenvector by $a_{11}$ instead of $a_{21}$, which would also be a valid eigenvector. While this is generally true, $a_{21}$ can be zero, whereas $a_{11}$ cannot be zero as transmission is required for a 2-port VNA to function.

Assuming $\bs{v}\_1$ is the computed eigenvector corresponding to $e^{-l\gamma}$ and $\bs{v}\_2$ is the eigenvector corresponding to $e^{l\gamma}$, we can compute $a_{21}/a_{11}$ and $a_{12}$ as follows:

\begin{equation}
a_{21}/a_{11} = \bs{v}\_1[1]/\bs{v}\_1[0]; \qquad a_{12} = \bs{v}\_2[0]/\bs{v}\_2[1]\label{eq:9}
\end{equation}
where the indexing of the vector starts at 0.

Before we continue discussing about how to de-normalize the coefficient $a_{21}/a_{11}$, we first need to solve for the right error-box $\bs{B}$. Similar to Eq. \eqref{eq:6}, we can also write an eigenvalue problem for $\bs{B}$ as follows:

\begin{equation}
(\bs{M}^{-1}\_\mathrm{thru}\bs{M}\_\mathrm{line})^T = \bs{B}^T\begin{bmatrix} e^{-l\gamma} & 0 \\\ 0 & e^{l\gamma}\end{bmatrix}(\bs{B}^T)^{-1}\label{eq:10}
\end{equation}

The transpose $(\cdot)^T$ is used here to push the inverse matrix to the right, otherwise, we would have it on the left side. Now, similar to Eq. \eqref{eq:8}, we can write the eigenvalue problem as:

\begin{equation}
(\bs{M}^{-1}\_\mathrm{thru}\bs{M}\_\mathrm{line})^T\begin{bmatrix} 1 \\\ b_{12}/b_{11}\end{bmatrix} = e^{-l\gamma}\begin{bmatrix} 1 \\\ b_{12}/b_{11}\end{bmatrix};\qquad (\bs{M}^{-1}\_\mathrm{thru}\bs{M}\_\mathrm{line})^T\begin{bmatrix} a_{21} \\\ 1\end{bmatrix} = e^{l\gamma}\begin{bmatrix} a_{21} \\\ 1\end{bmatrix}\label{eq:11}
\end{equation}

where $$\bs{B}^T = \begin{bmatrix}b_{11}&b_{21}\\b_{12}&1\end{bmatrix}$$. 

Now, suppose $\bs{u}\_1$ and $\bs{u}\_2$ are the computed eigenvectors of $(\bs{M}^{-1}\_\mathrm{thru}\bs{M}\_\mathrm{line})^T$ that correspond to the eigenvalues $e^{-l\gamma}$ and $e^{l\gamma}$, respectively. Then $b_{12}/b_{11}$ and $b_{21}$ are solved as:

\begin{equation}
b_{12}/b_{11} = \bs{u}\_1[1]/\bs{u}\_1[0];\qquad b_{21} = \bs{u}\_2[0]/\bs{u}\_2[1]\label{eq:12}
\end{equation}

### Coefficient de-normalization

To finalize the calibration, we need to de-normalize $a_{21}/a_{11}$ and $b_{12}/b_{11}$, and solve for the 7th error term $k$ (I haven't forgotten about it). The first what we can do is to use the Thru standard measurement and express the error-boxes in terms of their normalized coefficients. This looks something like:

\begin{equation}
\bs{M}\_\mathrm{thru} = k\underbrace{\begin{bmatrix} 1 & a_{12} \\\ a_{21}/a_{11} & 1\end{bmatrix}\begin{bmatrix} a_{11} & 0 \\\ 0 & 1\end{bmatrix}}\_{\bs{A}}\underbrace{\begin{bmatrix} b_{11} & 0 \\\ 0 & 1\end{bmatrix}\begin{bmatrix} 1 & b_{12}/b_{11} \\\ b_{21} & 1\end{bmatrix}}_{\bs{B}}\label{eq:13}
\end{equation}

By taking the inverse of the normalized matrices, we get:

\begin{equation}
\begin{bmatrix} 1 & a_{12} \\\ a_{21}/a_{11} & 1\end{bmatrix}^{-1}\bs{M}\_\mathrm{thru}\begin{bmatrix} 1 & b_{12}/b_{11} \\\ b_{21} & 1\end{bmatrix}^{-1} = \begin{bmatrix} ka_{11}b_{11} & 0 \\\ 0 & k\end{bmatrix}\label{eq:14}
\end{equation}

So, from Eq. \eqref{eq:14} we are able to solve for $k$. And with $k$ we can also solve for $a_{11}b_{11}$ as

\begin{equation}
a_{11}b_{11} = ka_{11}b_{11}/k\label{eq:15}
\end{equation}

The next step is to measure the reflect standard. The main property of the reflect standards is that it must be symmetric, that is, both ports must see the same reflect standard. This is illustrated in Fig. 3 below.

![Error-box model of a 2-port VNA when measuring a symmetric reflect standard.](block_error_box_model_symmetric_reflect.png)
_Fig. 3.: Error-box model of a 2-port VNA when measuring a symmetric reflect standard._

In the figure above, $\Gamma$ is the reflection coefficient of the reflect standard. We don’t know its value, but it should be the same at both ports. Generally, it is prefer to use a Short standard as it radiate less than an Open standard. But also any reflective load could also be used. The measured reflection seen at the left error-box is given by

\begin{equation}
\Gamma_a = \frac{a_{12}+a_{11}\Gamma}{1+a_{21}\Gamma}\quad\Longrightarrow\quad \Gamma=\frac{\Gamma_a-a_{12}}{a_{11}-a_{21}\Gamma_a}\label{eq:16}
\end{equation}

A similar equation can be obtained for the right error-box:

\begin{equation}
\Gamma_b = \frac{b_{11}\Gamma-b_{21}}{1-b_{12}\Gamma}\quad\Longrightarrow\quad \Gamma=\frac{\Gamma_b+b_{21}}{b_{11}+b_{12}\Gamma_b}\label{eq:17}
\end{equation}

> You might be wondering how I came with above equations, and why the measured reflection coefficient in the right error-box differs from the left error-box? Well, this comes from the famous formula for the input reflection coefficient measured from left-to-right (like with the left error-box)
\\[
\Gamma_{in} = S_{11} + \frac{S_{12}S_{21}\Gamma}{1-S_{22}\Gamma} = \frac{S_{11}+(S_{12}S_{21}-S_{11}S_{22})\Gamma}{1-S_{22}\Gamma}
\\]
The above equation so ubiquities, you can also find it on [Wikipedia](https://en.wikipedia.org/wiki/Scattering_parameters).
\
Now, you simply substitute the S-parameters with their corresponding T-parameters. Since above equations is defined from left-to-right, you can’t use it as it is for the right error-box (right-to-left). For that, you first flip the S-parameters, i.e., $S_{11} \Leftrightarrow S_{22}$ and $S_{21} \Leftrightarrow S_{12}$, and then substitute with their corresponding T-parameters. That is why, the formula for the left and right reflection coefficient are different.
{: .prompt-info }

Since both reflect standards are the same at both ports, we can define the following ratio:

\begin{equation}
\frac{a_{11}\Gamma}{b_{11}\Gamma} = \frac{a_{11}}{b_{11}} = \frac{\Gamma_a-a_{12}}{1-(a_{21}/a_{11})\Gamma_a}\frac{1+(b_{12}/b_{11})\Gamma_b}{\Gamma_b+b_{21}}\label{eq:18}
\end{equation}

So, now, we have $a_{11}b_{11}$ from the thru measurement and $a_{11}/b_{11}$ from the reflect measurement. Thus, we can solve for $a_{11}$ and $b_{11}$ as follows:

\begin{equation}
a_{11} = \pm\sqrt{\frac{a_{11}}{b_{11}}a_{11}b_{11}}; \qquad b_{11} = \frac{a_{11}b_{11}}{a_{11}}\label{eq:19}
\end{equation}

To resolve the sign ambiguity, we choose the answer that will give a $\Gamma$ closest to the expected value. For example, if the reflect standard is a short, we expect $\Gamma$ to be close to -1. Of course, you want to take any offset into account if your reflect standard is offsetted by some length.

Now, that we have $a_{11}$ and $b_{11}$, the matrices $\bs{A}$ and $\bs{B}$ are obtained as:

\begin{equation}
\bs{A}=\begin{bmatrix} 1 & a_{12} \\\ a_{21}/a_{11} & 1\end{bmatrix}\begin{bmatrix} a_{11} & 0 \\\ 0 & 1\end{bmatrix}; \qquad \bs{B} = \begin{bmatrix} b_{11} & 0 \\\ 0 & 1\end{bmatrix}\begin{bmatrix} 1 & b_{12}/b_{11} \\\ b_{21} & 1\end{bmatrix}\label{eq:20}
\end{equation}

> Everything we derived here is based on the definition of T-parameters as given in Eq. \eqref{eq:2}. Be carful with this, because there are authors that use a different definition for the T-parameters. For example, Michael Steer, in his open access books, he defines the T-parameters as:
\\[
\bs{T} = \frac{1}{S_{21}}\begin{bmatrix} 1 & -S_{22} \\\ S_{11} & S_{12}S_{21}-S_{11}S_{22}\end{bmatrix}; \qquad \bs{S} = \frac{1}{T_{11}}\begin{bmatrix}T_{21} & T_{11}T_{22}-T_{12}T_{21}\\\ 1 & -T_{12}\end{bmatrix}
\\]
As you can see, the order in Michael Steer definition is different, but he is consistent with his ordering in the backward conversion. His equations are correct, but if you use it, you want to modify the equations I presented accordingly.
\
Also, check out Michael Steer’s free books on microwave engineering: <https://repository.lib.ncsu.edu/handle/1840.20/36776>
{: .prompt-warning}

## Applying the calibration

Your DUT can be calibrated in terms of T-parameters using the following equation

\begin{equation}
\bs{T}\_\mathrm{dut} = \frac{1}{k}\bs{A}^{-1}\bs{M}\_\mathrm{dut}\bs{B}^{-1}\label{eq:21}
\end{equation}

Note that you must first correct for switch terms before applying the calibration in T-parameters. Finally, you can convert back to S-parameters by using eq. \eqref{eq:2}.

On the other hand, it is often desired to convert the obtained error terms to the conventional 12-term VNA model, where the switch term is already included. I won't derive the equations here, but you can find the details in Joel Dunsmore's book [3]. The 12 error terms are computed as follows:

* Single port 3 terms forward direction:
\\[
E_\mathrm{DF} = a_{12}; \qquad E_\mathrm{SF} = -a_{21}; \qquad E_\mathrm{RF} = a_{11} - a_{12}a_{21}
\\]

* Single port 3 terms reverse direction:
\\[
E_\mathrm{DR} = -b_{21}; \qquad E_\mathrm{SR} = b_{12}; \qquad E_\mathrm{RR} = b_{11} - b_{12}b_{21}
\\]

* Remaining 3 forward terms:
\\[
E_\mathrm{LF} = E_\mathrm{SR} + \frac{E_\mathrm{RR}\Gamma_\mathrm{f}}{1-E_\mathrm{DR}\Gamma_\mathrm{f}}; \qquad E_\mathrm{TF} = \frac{1}{k(1-E_\mathrm{DR}\Gamma_\mathrm{f})}; \qquad E_\mathrm{XF} = 0
\\]

* Remaining 3 terms in the opposite direction:
\\[
E_\mathrm{LR} = E_\mathrm{SF} + \frac{E_\mathrm{RF}\Gamma_\mathrm{r}}{1-E_\mathrm{DF}\Gamma_\mathrm{r}}; \qquad E_\mathrm{TR} = \frac{kE_\mathrm{RR}E_\mathrm{RF}}{1-E_\mathrm{DF}\Gamma_\mathrm{r}}; \qquad E_\mathrm{XR} = 0
\\]

The equations for applying the 12 error terms to obtain the calibrated measurement can be found in [3].

## Auxiliary operations

Below are some operations that are not strictly needed for the calibration itself, but can be useful to extend the calibration beyond the discussion above.

### Extracting the propagation constant

After you have computed the eigenvalues and have been able to distinguish between the signs of the exponents. You can calculate the propagation constant $\gamma$ the other way around. Suppose the first eigenvalue $\lambda_1$ is associated with $e^{-\gamma l}$ and the second eigenvalue $\lambda_2$ is associated with $e^{+\gamma l}$. Then we define a new ratio as follows:
\begin{equation}
e^{2\gamma l} = \lambda_2/\lambda_1\label{eq:22}
\end{equation}

We extract the exponent using the logarithmic $e$, which gives us
\begin{equation}
\log(\lambda_2/\lambda_1) = 2\gamma^{\mathrm{unwrap}} l+ jn2\pi, \qquad \text{for any } n\in\\{0,\pm1,\pm2,\ldots\\}\label{eq:23}
\end{equation}

The extra phase-wrapping factor you see is the result of taking the logarithm of complex numbers, as is unique only in the principle brach $[-\pi,\pi]$. You can read more about the complex logarithm on [Wikipedia](https://en.wikipedia.org/wiki/Complex_logarithm). For now, we just need to determine the wrap factor to accurately determine $\gamma$. We do this by using the approximation we used to distinguish the eigenvalues. So, the equation to determine $n$ is as follows
\begin{equation}
n = \mathrm{round}\left(
    \frac{\Im{\log(\lambda_2/\lambda_1) - 2\gamma^{\mathrm{est}}l}}{2\pi}
    \right)
    \label{eq:24}
\end{equation}

So the unwrapped $\gamma$ is solved as 
\begin{equation}
\gamma^{\mathrm{unwrap}} = \frac{\log(\lambda_2/\lambda_1) - jn2\pi }{2l}\label{eq:25}
\end{equation}

> If you are wondering which eigenvalues to take, either from the forward direction (left error box $\bs{A}$), or from the backward direction (right error box $\bs{B}$). The answer is that it doesn't matter. Both versions give the same eigenvalues, only the eigenvectors change.
{: .prompt-info }

### Shift reference plane

By default, after calibration, the reference plane is set to the center of the thru standard. Therefore, if you were to apply the calibration, it would be referenced to the center of the thru line. However, one of the advantages of the TRL calibration is that you have access to the propagation constant as derived above. With the propagation constant, you can adjust the error boxes to a new reference plane. So we can update the reference plane with the following equations:

\begin{equation}
\bs{A}\_\mathrm{new} = \bs{A}\begin{bmatrix}e^{-2\gamma d} & 0 \\\\ 0 & 1\end{bmatrix}; \qquad
\bs{B}\_\mathrm{new} = \begin{bmatrix}e^{-2\gamma d} & 0 \\\ 0 & 1\end{bmatrix}\bs{B}; \qquad
k\_\mathrm{new} = e^{2\gamma d}k \label{eq:26}
\end{equation}

where $d$ is the offset as seen from the port. A positive offset shifts the plane away from the ports, while a negative offset shifts the plane toward the ports.

![Illustration of negative plane shift of a TRL calibration done with microstrip line.](illustration_plane_shift.png)
_Fig. 4.: Illustration of negative plane shift of a TRL calibration done with a thru standard as microstrip line._

> I will leave it to you as a homework to figure out how I derived these equations. Hint: start by placing a line of length $2d$ between the error boxes $\bs{A}$ and $\bs{B}$.
{: .prompt-info }

### Renormalize reference impedance

By default, the reference impedance after calibration is set to the characteristic impedance of the line standard. For many types of transmission lines, the characteristic impedance is frequency dependent. The frequency dependency can complicate the interpretation of the calibrated measurements and make it impossible to compare different TRL calibrations of different transmission line types. Therefore, it is often desirable to change the impedance to a different value, most commonly to a frequency-independent value (e.g., 50 ohms).

If the value of the characteristic impedance of the line standard is known, the updated error terms for a new reference impedance can be obtained from the following equations:

\begin{equation}
\bs{A}_\mathrm{new} = \frac{\widehat{\bs{A}}}{\widehat{\bs{A}}[1,1]}, \qquad \text{where, } \widehat{\bs{A}} = \bs{A}\bs{Q}^{nm} \label{eq:27}
\end{equation}

\begin{equation}
\bs{B}_\mathrm{new} = \frac{\widehat{\bs{B}}}{\widehat{\bs{B}}[1,1]}, \qquad \text{where, } \widehat{\bs{B}} = \bs{Q}^{mn}\bs{B} \label{eq:28}
\end{equation}

\begin{equation}
k_\mathrm{new} = k\widehat{\bs{A}}[1,1]\widehat{\bs{B}}[1,1]\label{eq:29}
\end{equation}

The matrix $\bs{Q}^{nm}$ is the impedance transformation matrix [4], which transform the original reference impedance $Z_n$ to the new impedance $Z_m$. The expression for $\bs{Q}^{nm}$ is given as follows (refer to [4] for derivation):  
\begin{equation}
\bs{Q}^{nm} = \left(\bs{Q}^{mn}\right)^{-1} = \frac{1}{\sqrt{1-\Gamma_{nm}^2}}\begin{bmatrix} 1 & \Gamma_{nm} \\\ \Gamma_{nm} & 1\end{bmatrix}, \qquad \text{where, } \Gamma_{nm} = \frac{Z_m-Z_n}{Z_m+Z_n}\label{eq:30}
\end{equation}

To learn more about the origin of the transformation matrix, see section _"3.7 Change of Reference Impedance"_ in [4].

> In the equations \eqref{eq:27} and \eqref{eq:28} I normalized the new transformed error boxes by their last element. This is because I defined $\bs{A}$ and $\bs{B}$ so that their last element is normalized to one (go back to equation \eqref{eq:3}). As a result, the factor $k$ had to be updated accordingly by the factored out elements.
{: .prompt-info }

## Design TRL kit

So this is the most important part of the discussion. Even if you master the math above, if you still make crappy TRL kit, there is no math to fix it. So, the goal of this section is to make you aware of the limitations of TRL kit, which aspects you need to pay attention to, and which of those you need to relax your thinking about.

### Line length from frequency range and vice versa

(on going writing...)

Before we move forward with the discussion below, the term "line length" refer to the length difference between the line and the thru standard, i.e., the length of the line relative to the thru. By definition, the thru standard has a length of zero. Please keep that in mind when ever working with a TRL calibration.

Back to the topic on hand, a very important question you should be asking yourself when designing a TRL is: _What is the frequency range the kit should operate in?_. However, an even more fundamental question you should be asking is _why is there a frequency limit to begin with?_ 

The answer is quite simple. Recall the eigenvalue problem in Eq. \eqref{eq:6}. The eigenvalues depend on $\gamma$, but $\gamma$ itself depend on frequency. This is generally given by the relation:
\begin{equation}
\gamma = \alpha + j\beta = \frac{2\pi f}{c_0}\sqrt{-\epsilon_\mathrm{r,eff}},
\label{eq:31}
\end{equation}
where $\epsilon_\mathrm{r,eff}$ is the relative effective permittivity seen by the wave, which is also a complex-valued number, often written as
\begin{equation}
\epsilon_\mathrm{r,eff} = \epsilon_\mathrm{r,eff}^\prime - j\epsilon_\mathrm{r,eff}^{\prime\prime}
\label{eq:32}
\end{equation}

The real-part of $\gamma$ describe the losses experience by the wave, whereas the imaginary-part describe the phase of the wave. The frequency limitation of TRL comes from the latter term. Therefore, it is easier to explain when we consider lossless case:
\begin{equation}
\gamma_\mathrm{lossless} = j\beta = \frac{2\pi f}{c_0}j\sqrt{\epsilon_\mathrm{r,eff}^\prime}
\label{eq:33}
\end{equation}

Now, because $\gamma$ has a frequency dependency, the complex exponential terms, i.e., the eigenvalues, will rotate with frequency. Since the two eigenvalues have different exponent sign, their rotation will be opposite to one another. There is certain frequencies, where both eigenvalues meet and equal each other. That is when the exponent $\gamma l = \pm j\pi$. Many people call the term $\gamma l$ the electrical length of the transmission line. A visualization of the eigenvalues rotating in the complex plane with respect to frequency is shown in figure below.

![Illustration of the rotation of the eigenvalues with frequency.](trl_eigvals_wide.gif)
_Fig. 5.: Illustration of the rotation of eigenvalues in the complex plane with respect to frequency._

> I created the animation of the eigenvalues with [desmos](https://www.desmos.com/calculator). What do you think will happen to the eigenvalues if we include losses? Fun fact, when you have high losses, the design criteria for the length and frequency get relaxed. Basically, lossless is the worst case scenario to have.
{: .prompt-info }

So, clearly we don't want our TRL kit operate in a frequency range, where it could cross those critical points. Because of the rotation nature of the eigenvalues, we can also visualize the multi-band nature of TRL with sketch below. Where specify a phase region, where we call phase margin. Basically, we stay out from this phase margin region, then we are OK. You can see in the sinusoidal plot who these phase margin region repeat, which also show multiple pass-bands. The choice what deems a good phase margin is up to debate, but the rule of tumbe is to not to use a TRL kit with a phase margin less than $20^{\circ}$.

![Illustration of the rotation of the eigenvalues with frequency.](trl_eigvals_wide.gif)
_Fig. 6.: Illustration of the rotation of eigenvalues in the complex plane with respect to frequency._

Now, to determine the appropriate length given a frequency band, we write Eq. \eqref{eq:33} as an inequality and multiply it with the length to get the electrical length. We bound electrical length to the phase margin, and include the $\pi$ periodicity.

\begin{equation}
\pi n + \pi\frac{\phi}{180}\leq l\frac{2\pi f}{c_0}\sqrt{\epsilon_\mathrm{r,eff}^\prime} \leq \pi n + \left( 1-\frac{\phi}{180}\right)\pi, \qquad n=0,1,2,\ldots
\label{eq:34}
\end{equation}
where $\phi$ is the phase margin in degrees. We can simplify above equation by canceling $\pi$ from all sides. This gives us:
\begin{equation}
n + \frac{\phi}{180}\leq l\frac{2f}{c_0}\sqrt{\epsilon_\mathrm{r,eff}^\prime} \leq n + 1-\frac{\phi}{180}
\label{eq:35}
\end{equation}

Now, we can isolate the frequency term in the middle by multiplying the equation by $c_0/\left(2l\sqrt{\epsilon_\mathrm{r,eff}^\prime}\right)$.
\begin{equation}
\frac{n + \frac{\phi}{180}}{2l\sqrt{\epsilon_\mathrm{r,eff}^\prime}}c_0 \leq f \leq \frac{n + 1-\frac{\phi}{180}}{2l\sqrt{\epsilon_\mathrm{r,eff}^\prime}}c_0
\label{eq:36}
\end{equation}

Therefore, we can define the lower and upper frequency limits given a length $l$, a relative effective permittivity $\epsilon_\mathrm{r,eff}^\prime$, and a phase margin $\phi$ as follows: 
\begin{equation}
f_\mathrm{min} = \frac{n + \frac{\phi}{180}}{2l\sqrt{\epsilon_\mathrm{r,eff}^\prime}}c_0, \qquad  f_\mathrm{max} = \frac{n + 1-\frac{\phi}{180}}{2l\sqrt{\epsilon_\mathrm{r,eff}^\prime}}c_0
\label{eq:37}
\end{equation}
The value of $n$ allow you to move between bands.

Of course, our goal is not only to figure out the frequency range given an already designed TRL kit, but we would also like to design the length given a frequency range specification. Now, observe from Eq. \eqref{eq:37} that the length is common for both frequency limits. Therefore, we can reformulate Eq. \eqref{eq:37} in terms of $l$ for both frequency limits:
\begin{equation}
l = \frac{c_0}{2f_\mathrm{min}\sqrt{\epsilon_\mathrm{r,eff}^\prime}}\left(n + \frac{\phi}{180}\right), \qquad  l = \frac{c_0}{2f_\mathrm{max}\sqrt{\epsilon_\mathrm{r,eff}^\prime}}\left(n + 1 - \frac{\phi}{180}\right)
\label{eq:38}
\end{equation}

Now, since $l$ is the same for both cases, we can equate the equations and try to

### Determining the Physical Length of the Thru Standard

One common question I receive is how to determine the absolute physical length of the thru standard. While we know that calibration is performed with reference to the center of the thru standard, the length of the standard itself can vary. So, how long should it be? Can it be too short or too long?

The answer is not straightforward, but here are some things to consider that will help you decide on the appropriate length:

1. Crosstalk: Ensure that the connectors and ports do not couple with each other. This can happen if the thru standard is too short and any radiation due to imperfections in the connectors/probes couples with each other. This is usually only a problem if you design bad transitions to your transmission line. For PCBs with coaxial connectors, it is advisable to have the total thru length (from the outer edge of the connectors) be three times the length of a connector. This will provide one connector length between the inner edges of the connectors. To be on the safe side, it is best to simulate your structure in EM simulation.

2. Maximum length: This depends on the losses associated with your medium. If you are dealing with a highly lossy material, keep the thru as short as possible without any coupling problems. If you are dealing with a low-loss material, you can make it longer. As a general rule, make it a little longer than the minimum length that causes no problems. If you cannot simulate your structure, you can use an analytical model of transmission (e.g. from [scikit-rf](https://scikit-rf.readthedocs.io/en/latest/api/media/index.html)) and look for the length that will result in 3 dB attenuation. Set the length below this 3 dB length. The exact number depends on manufacturing capabilities and cost.

3. Tips for specific types of transmission lines:
   * Waveguides: Rather than connecting the flanges together, it is recommended that you implement the thru as a line, and move the plane back with the propagation constant obtained from calibration. This helps minimize the effects of imperfections at the connection interface.
   * PCB: The connectors, particularly the transition from the connector to the PCB, can cause problems. If you design your transition without consideration, you are prone to generating surface waves and radiation at the connector interface. One simple thing you can do to reduce these effects is to introduce a taper between your connector and the final transmission line type and smooth out the traces. Avoid sharp corners if possible. While corners may be acceptable at lower frequencies (below 18 GHz), it becomes problematic at higher frequencies.

### What should the reflect standard be?

There are several misconceptions surrounding the implementation of the reflect standard. Some people claim that it must be implemented in the same media as the line and thru standards, while others suggest adding an offset to the reflect or using a single-mode propagation. While these statements are based on some truth, they are often overgeneralized and can create confusion.

The three strict properties that your reflect standard should possess are:

1. It must be reflective (duhhh~).
2. It must be symmetrical from both ports.
3. You must know an estimate of the reflection to resolve sign ambiguity during calibration.

The above conditions are easier said than done, and there are possible problems you may encounter during implementation:

* Replicating the reflect on both ports is not that difficult, but you can make things difficult if you're not careful. For example, an open standard of a poorly designed connector or interface can resonate in the frequency band you are working with. You can still design stable open standards, but if you are unsure, it is better to implement your reflect as a short standard, since shorts tend to radiate and resonate less. If you are still unsure of your design, an EM simulation can help. However, a typical short standard should suffice for most interfaces.

* Adding an offset is common practice. For example, a short standard on a PCB is typically offsetted by some distance from the reference plane to avoid coupling the reflect standard to the probes/connectors, which can cause resonance and make it difficult to replicate the standard on both ports. For most interface types, you probably won't need an offset. If you are doing on-wafer probing, you may want to add the offset to be on the safe side, but even there, it depends on the coupling strength between your probe and the substrate on which the standard is implemented. If you are unsure, adding the offset won't hurt, but it is unnecessary to do so for every TRL kit.

* Some people believe that the reflect standard must contain the same propagation mode as the line and thru standards, otherwise the reflect standard is invalid. This is a misconception that arises from a misunderstanding of the requirements for reflect. However, it could be that a different mode would make the reflect less stable compared to your original mode, which in turn would make it difficult to replicate the standard on both sides. But, this is not always the case. If you have made your line standards as microstrip on a PCB, but you didn't make a reflect standard, you can still do your calibration by leaving the cable open on both ends (assuming it does not resonate or radiate), and it would be a valid reflect standard. You can also use a coaxial short standard from your SOLT box, and it would be a valid reflect standard. Therefore, it is unnecessary to overthink the propagation mode of your reflect standard.

* Regardless of how you make the reflect standard, the most important thing is to know its estimated value to resolve the sign ambiguity. Most people make their reflect standard on the same media as the line and thru standards, making it easy to predict its reflection coefficient. One useful trick to solve this problem is to approximate the reflectance at the first frequency point and use the new estimate extracted from the previous frequency point at each new frequency point. This is particularly helpful when implementing the reflect as open or short standards, which have a reflection coefficient of approximately 1 and -1, respectively, at low frequencies.

## References

[1] G. F. Engen and C. A. Hoer, "Thru-Reflect-Line: An Improved Technique for Calibrating the Dual Six-Port Automatic Network Analyzer," _in IEEE Transactions on Microwave Theory and Techniques_, vol. 27, no. 12, pp. 987-993, Dec. 1979, doi: [10.1109/TMTT.1979.1129778](https://dx.doi.org/10.1109/TMTT.1979.1129778).

[2] R. B. Marks, "Formulations of the Basic Vector Network Analyzer Error Model including Switch-Terms," _50th ARFTG Conference Digest_, 1997, pp. 115-126, doi: [10.1109/ARFTG.1997.327265](https://dx.doi.org/10.1109/ARFTG.1997.327265).

[3] J. Dunsmore, Handbook of Microwave Component Measurements. John Wiley & Sons, Ltd, 2020, doi: [10.1002/9781119477167](https://dx.doi.org/10.1002/9781119477167)

[4] Marks, R. and Williams, D. (1992), A General Waveguide Circuit Theory, Journal of Research (NIST JRES), National Institute of Standards and Technology, Gaithersburg, MD, doi: [10.6028/jres.097.024](https://dx.doi.org/10.6028/jres.097.024)

<!-- EOF -->