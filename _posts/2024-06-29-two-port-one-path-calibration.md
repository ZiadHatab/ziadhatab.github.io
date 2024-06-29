---
title: Two-port/One-path Calibration
date: 2024-06-29 12:00:00 +0100
categories: [Tutorial]
tags: [vna, calibration, s-parameters, SOLT, three-receivers]     # TAG names should always be lowercase
math: true
img_path: ../../../assets/img/posts_img/
image: # TRL_waveguide_GCPW_kit.png
---

If you ever find yourself using a simple vector network analyzer (VNA) that can only take S11 and S21 measurements, the explanation here might be an eye-opener to how you can actually calibrate such a system and what assumptions you can make to get the answer.

I was inspired to write this blog based on the examples provided by the folks at Scikit-rf [1]. Additionally, the Math Reference of METAS VNA Tools briefly explains this method [2]. In summary, it is just half of an SOLT (Short-Open-Load-Thru) calibration. However, you still need fully defined Short, Open, Load (SOL) standards and a fully defined two-port standard (often a Thru standard is used).

Keep in mind, most modern VNAs operate in full-reflectometry mode, and some use half-reflectometry but still use switches to obtain two-port measurements. This tutorial is specific to the case where no switches are used (hence one-path). Because I love linear algebra, I will try to write everything in vector-matrix notations as much as possible. I believe it makes the content more approachable for many people, as opposed to using only algebraic expressions that can sometimes lead to more confusion than clarification. So, let’s get started…

You can find a sample code demonstrating the math presented here on my GitHub: <https://github.com/ZiadHatab/two-port-one-path-calibration>

## The Block Diagram

As with any VNA calibration topic, we begin by examining the block diagram that illustrates the error boxes of the VNA, as shown in the diagram below. Please note that for simplicity, we will disregard any cross-talk terms.

![Error box model of a two-port VNA with one-path source excitation.](two-port_one-path_model.png)
_Error box model of a two-port VNA with one-path source excitation._

Because it is a one-path measurement, we only have one set of measurements of the waves $\hat{a}_1, \hat{b}_1, \hat{b}_2$, which we can write down in combination with the T-parameters description of the networks:

$$
\begin{bmatrix}
\hat{b}_1\\
\hat{a}_1
\end{bmatrix} = \underbrace{k_ak_b}_{k}\boldsymbol{A}\boldsymbol{T}_\mathrm{dut}\boldsymbol{B}\begin{bmatrix}
\hat{a}_2\\
\hat{b}_2
\end{bmatrix}
\label{eq:1}
$$

where $\hat{a}_2$ denotes the incident signal that would have been measured if there was a coupler to measure it at port 2. The error boxes ***A*** and ***B*** are given as follows:

$$
\boldsymbol{A} = \begin{bmatrix}
a_{11} & a_{12}\\
a_{21} & 1
\end{bmatrix}; \qquad \boldsymbol{B} = \begin{bmatrix}
b_{11} & b_{12}\\
b_{21} & 1
\end{bmatrix}
\label{eq:2}
$$

The scalars $k_a$ and $k_b$ are factored from the matrices ***A*** and ***B*** for simplicity as they can be lumped together as $k = k_a k_b$. Lastly, the ***T*** matrix in the middle represents the device under test (DUT) or the standards to be measured.

Currently, the measurement in Equation (1) gives the wave parameters. To eliminate the term $\hat{a}_2$ from the equation, S-parameters are used instead to describe the measurements. The measurement of S11 and S21 is determined by dividing both sides of the equation by $\hat{a}_1$, which gives the following:

$$
\begin{bmatrix}
\hat{b}_1/\hat{a}_1\\
\hat{a}_1/\hat{a}_1
\end{bmatrix} = k\boldsymbol{A}\boldsymbol{T}_\mathrm{dut}\boldsymbol{B}\begin{bmatrix}
\hat{a}_2/\hat{a}_1\\
\hat{b}_2/\hat{a}_1
\end{bmatrix}\quad\Longrightarrow\quad
\begin{bmatrix}
S_{11}\\
1
\end{bmatrix} = k\boldsymbol{A}\boldsymbol{T}_\mathrm{dut}\boldsymbol{B}\begin{bmatrix}
\hat{a}_2/\hat{a}_1\\
S_{21}
\end{bmatrix}
\label{eq:3}
$$

We can further simplify the ratio $\hat{a}\_2/\hat{a}\_1$ by writing it in terms of the termination reflection coefficient $\Gamma_{21}$ as follows:

$$
\frac{\hat{a}_2}{\hat{a}_1} = \frac{\hat{a}_2\hat{b}_2}{\hat{a}_1\hat{b}_2} = S_{21}\Gamma_{21}, \qquad \text{where...}\ \Gamma_{21} = \frac{\hat{a}_2}{\hat{b}_2}
\label{eq:4}
$$

Therefore, the measurements can be written as follows:

$$
\begin{bmatrix}
S_{11}\\
1
\end{bmatrix} = k\boldsymbol{A}\boldsymbol{T}_\mathrm{dut}\boldsymbol{B}\begin{bmatrix}
\Gamma_{21}\\
1
\end{bmatrix}S_{21}
\label{eq:5}
$$

Note that $\Gamma_{21}$ is constant and can be treated as part of the error terms.

## Solving the Error Terms

In total, we are able to solve for 5 terms, which require us to use SOL standards and a thru connection (or any known 2-port device). We start with the SOL standards. These one-port devices allow us to solve for the error terms in the matrix ***A*** which comprises three error terms. The block diagram below illustrates the measurement concept.

![Illustration of one-port measurements.](one-port_model.png)
_Illustration of one-port measurements._

The way we describe the measurement given a connected device is as follows:

$$
S_{11}^{(i)} = \frac{a_{11}\Gamma^{(i)} + a_{12}}{a_{21}\Gamma^{(i)} + 1}
\label{eq:6}
$$

Where $\Gamma^{(i)}$ describes the actual reflection of the connected device (i.e., known value) and $S_{11}^{(i)}$ is the corresponding measurement. The superscript $(\cdot)^{(i)}$ indicate the $i$-th measured device. The relation mentioned above is known as a linear fractional transform, or a Möbius transform. You can check [3] for more information on this. The reason the term $k_a$ does not appear in the above equation is because it actually appears in both the numerator and denominator, which cancels out and therefore not mentioned.

For the solution of the error terms $a_{11}, a_{12}, a_{21}$, they can be solved by solving the linear system of equations below:

$$
\begin{bmatrix} 
\Gamma^{(1)} & 1 & -S_{11}^{(1)}\Gamma^{(1)} \\[.5mm]
\Gamma^{(2)} & 1 & -S_{11}^{(2)}\Gamma^{(2)} \\[.5mm]
\Gamma^{(3)} & 1 & -S_{11}^{(3)}\Gamma^{(3)}
\end{bmatrix}
\begin{bmatrix} a_{11}\\ a_{12}\\ a_{21} \end{bmatrix} = \begin{bmatrix}S_{11}^{(1)}\\S_{11}^{(2)}\\S_{11}^{(3)}\end{bmatrix}
\label{eq:7}
$$

where superscripts 1, 2, and 3 correspond to open, short, and load standards in no particular order. The solution for $a_{11}, a_{12}, a_{21}$ is obtained by taking the inverse of the system matrix.

Now that we have solved for the matrix ***A***, we can move forward and solve for two more error terms with the help of the thru standard (or any known two-port device). Recall the equation describing the measurement of a two-port device under the one-path model:

$$
\begin{bmatrix}
S_{11}\\
1
\end{bmatrix} = k\boldsymbol{A}\boldsymbol{T}_\mathrm{thru}\boldsymbol{B}\begin{bmatrix}
\Gamma_{21}\\
1
\end{bmatrix}S_{21}
\label{eq:8}
$$

Since we already know ***A***, we can multiply its inverse on both sides. We can also multiply the inverse of the thru standard $\boldsymbol{T}_\mathrm{thru}$ on both sides, as it is also assumed to be known. This leaves us with the following:

$$
\boldsymbol{T}_\mathrm{thru}^{-1}\boldsymbol{A}^{-1}\begin{bmatrix}
S_{11}\\
1
\end{bmatrix} = k\boldsymbol{B}\begin{bmatrix}
\Gamma_{21}\\
1
\end{bmatrix}S_{21}
\label{eq:9}
$$

A further simplification, as S21 is not zero since we are measuring a transmissive device (i.e., two-port), we can simplify and divide both sides by S21 and combine the remaining error terms into a single vector of two elements:

$$
\boldsymbol{T}_\mathrm{thru}^{-1}\boldsymbol{A}^{-1}\begin{bmatrix}
S_{11}/S_{21}\\
1/S_{21}
\end{bmatrix} = k\boldsymbol{B}\begin{bmatrix}
\Gamma_{21}\\
1
\end{bmatrix} = \begin{bmatrix}
\alpha_2\\
\beta_2
\end{bmatrix}
\label{eq:10}
$$

Therefore, the terms $\alpha_2$ and $\beta_2$ summarize the remaining two error terms, thus leading us to solve for a total of five error terms. The reason I called them $\alpha$ and $\beta$ will become clear in the next section when talking about applying the calibration. Generally speaking, we cannot actually fully correct for a two-port device with these five error terms. A full correction requires measurements from both sides (i.e., two paths). However, we can make assumptions about the DUT to help simplify the calibration process.

## Apply the Calibration

Before we start with tackling the application of the calibration, I just want to rewrite the measurement model with the usage of the known error terms in a form that we can transform to S-parameter form instead of T-parameters. We start with the basic error box equation we already know:

$$
\begin{bmatrix}
S_{11}\\
1
\end{bmatrix} = k\boldsymbol{A}\boldsymbol{T}_\mathrm{dut}\boldsymbol{B}\begin{bmatrix}
\Gamma_{21}\\
1
\end{bmatrix}S_{21}
\label{eq:11}
$$

From the error box model, we know the matrix $\boldsymbol{A}$ and the two error terms combined from the terms $k, \boldsymbol{B}, \Gamma_{21}$. We can first do some rearrangement by multiplying both sides with the inverse of $\boldsymbol{A}$ and substitute the values for the product $k, \boldsymbol{B}, \Gamma_{21}$ with $\alpha$ and $\beta$ as derived earlier:

$$
\boldsymbol{A}^{-1}\begin{bmatrix}
S_{11}\\
1
\end{bmatrix} = \boldsymbol{T}_\mathrm{dut}\begin{bmatrix}
\alpha_2\\
\beta_2
\end{bmatrix}S_{21}
\label{eq:12}
$$

A final step is to combine the error terms with the S-parameter measurements (S11 and S21) to describe the T-parameters as a linear transformation of a vector, as they were originally defined for waves:

$$
\begin{bmatrix}
\hat{\beta}_1\\
\hat{\alpha}_1
\end{bmatrix} = \boldsymbol{T}_\mathrm{dut}\begin{bmatrix}
\hat{\alpha}_2\\
\hat{\beta}_2
\end{bmatrix}
\label{eq:13}
$$

Where the terms $\alpha$ and $\beta$ on the left-hand side represent the combined S11 measurement and the inverse $\boldsymbol{A}$, and the terms $\alpha$ and $\beta$ on the right-hand side represent the derived error terms combined with the S21 measurement. I chose to call them $\alpha$ and $\beta$ to resemble the a and b waves and how they define T-parameters, as we started with the first equation. However, it's important to note that these terms are not the actual waves in the absolute sense, as any scaling applied to both sides of the equation would still make it valid. Nevertheless, their ratio is unique and can be thought of as waves for the purpose of S-parameter calculations. We can rewrite the above equation in terms of the DUT S-parameters instead of T-parameters as follows:

$$
\begin{bmatrix}
\hat{\beta}_1\\
\hat{\beta}_2
\end{bmatrix} = \boldsymbol{S}_\mathrm{dut}\begin{bmatrix}
\hat{\alpha}_1\\
\hat{\alpha}_2
\end{bmatrix}
\label{eq:14}
$$

Lastly, recall that the terms $\hat{\alpha}$ and $\hat{\beta}$ are defined as follows:

$$
\begin{bmatrix}
\hat{\beta}_1\\
\hat{\alpha}_1
\end{bmatrix} = \boldsymbol{A}^{-1}\begin{bmatrix}
S_{11}\\
1
\end{bmatrix}; \qquad \begin{bmatrix}
\hat{\alpha}_2\\
\hat{\beta}_2
\end{bmatrix}= \begin{bmatrix}
\alpha_2\\
\beta_2
\end{bmatrix}S_{21}
\label{eq:15}
$$

### DUT is a one-port

This case is straightforward. Take the equation we used to describe the measurement of a one-port device and reverse calculate for the actual device, as we already know the error terms $a_{11}, a_{12}, a_{21}$:

$$
S_{11}^{(\mathrm{dut})} = \frac{a_{12}-S_{11}}{S_{11}a_{21}-a_{11}}
\label{eq:16}
$$

where $S_{11}$ is the measured reflection and $S_{11}^{(\mathrm{dut})}$ is the calibrated reflection of the DUT.

An alternative way to see this problem is to go back to the equation (14) where the DUT S-parameters are described with a vector-matrix product:

$$
\begin{bmatrix}
\hat{\beta}_1\\
\hat{\beta}_2
\end{bmatrix} = \underbrace{\begin{bmatrix}S_{11}^{(\mathrm{dut})} & S_{12}^{(\mathrm{dut})}\\
S_{21}^{(\mathrm{dut})} & S_{22}^{(\mathrm{dut})}\end{bmatrix}}_{\boldsymbol{S}_\mathrm{dut}}\begin{bmatrix}
\hat{\alpha}_1\\
\hat{\alpha}_2
\end{bmatrix}
\label{eq:17}
$$

Given that S21 and S12 are zero, and S22 cannot be measured as there is no source at port-2 (but also due to the fact that $\hat{\alpha}_2 = \hat{\beta}_2 = 0$). In the end, we have:

$$
S_{11}^{(\mathrm{dut})} = \frac{\hat{\beta}_1}{\hat{\alpha}_1}
\label{eq:18}
$$

### Partial correction of a two-port DUT

A two-port DUT has four terms, while a one-path measurement delivers only two measurements (S11 and S21). Even with knowledge of the five error terms, we cannot recover all S-parameters of the DUT simultaneously. However, if we make assumptions on two of the S-parameters of the DUT, then we end up with two unknowns that we can solve for with the two measurements.

This method of partial correction of a DUT is often referred to as "Enhanced Response". To be honest, there is nothing enhanced; we are just making assumptions. The best calibration is full two-port when you measure in two paths. The "enhanced" aspect comes from the S11 measurement when the DUT also transmits. In the previous example, we discussed just a simple one-port DUT where S21=S12=0, and the calculation for S11 is straightforward. Now, if there is transmission, the calibrated measurement of S11 and S21 will show dependence on S22 and S12. This can be written explicitly from the matrix equation using $\hat{\alpha}$ and $\hat{\beta}$ notations:

$$
S_{11}^{(\mathrm{dut})} = \frac{\hat{\beta}_1}{\hat{\alpha}_1} - S_{12}^{(\mathrm{dut})}\frac{\hat{\alpha}_2}{\hat{\alpha}_1}; \qquad S_{21}^{(\mathrm{dut})} = \frac{\hat{\beta}_2}{\hat{\alpha}_1} - S_{22}^{(\mathrm{dut})}\frac{\hat{\alpha}_2}{\hat{\alpha}_1}
\label{eq:19}
$$

Therefore, to solve S11 we need to make assumptions about S12; and to solve for S21 we need to make assumptions about S22. I don't think there is a consistent definition of this partial calibration in terms of assumptions. In scikit-rf and METAS notes [1,2], they assume S12=S22=0, which makes sense if you are working with amplifiers.

However, you could also assume S22=0 while assuming reciprocity, thus S12=S21. That is, you first solve for S21 by setting S22=0 and substitute it back in the first equation in place of S12, and then solve for S11. Of course, you could also assume symmetry, where S11=S22 and S12=S21, then you have two equations with just two unknowns (scikit-rf calls this fakeflip [1]). At the end of the day, these are all assumptions, and the best ones that suit your DUT is up to you to decide. What you need to know is that the one-path measurements, after knowing the five error terms, give you at best these two equations mentioned above. More than that, you just assume!!! (or make more measurements, which I will talk about next).

### Full correction of a two-port DUT

This is the rigorous way to perform a two-port calibration with a system that allows you to measure in one path. The trick is to take two measurements of the DUT in two orientations, once forward and the other reversed (flipped). We can write the two equations as follows:

$$
\begin{bmatrix}
\hat{\beta}_1^{(1)}\\
\hat{\beta}_2^{(1)}
\end{bmatrix} = \underbrace{\begin{bmatrix}S_{11}^{(\mathrm{dut})} & S_{12}^{(\mathrm{dut})}\\
S_{21}^{(\mathrm{dut})} & S_{22}^{(\mathrm{dut})}\end{bmatrix}}_{\boldsymbol{S}_\mathrm{dut}}\begin{bmatrix}
\hat{\alpha}_1^{(1)}\\
\hat{\alpha}_2^{(1)}
\end{bmatrix}; \qquad 
\begin{bmatrix}
\hat{\beta}_1^{(2)}\\
\hat{\beta}_2^{(2)}
\end{bmatrix} = \underbrace{\begin{bmatrix}S_{22}^{(\mathrm{dut})} & S_{21}^{(\mathrm{dut})}\\
S_{12}^{(\mathrm{dut})} & S_{11}^{(\mathrm{dut})}\end{bmatrix}}_{\boldsymbol{P}\boldsymbol{S}_\mathrm{dut}\boldsymbol{P}}\begin{bmatrix}
\hat{\alpha}_1^{(2)}\\
\hat{\alpha}_2^{(2)}
\end{bmatrix}
\label{eq:20}
$$

where the superscripts on $\alpha$ and $\beta$ indicate the measurement orientation. By flipping the DUT, its S-parameters are swapped. This action is equivalent to multiplying the original matrix with a permutation matrix from both sides, i.e., matrix ***P***, which is given as follows:

$$
\boldsymbol{P} = \begin{bmatrix} 0 & 1\\
1 & 0\end{bmatrix}; \qquad \boldsymbol{P}^T=\boldsymbol{P}^{-1} = \boldsymbol{P}
\label{eq:21}
$$

To simplify, we can multiply the inverse of $\boldsymbol{P}$ on both sides of the second equation to eliminate one of the permutations:

$$
\boldsymbol{P}\begin{bmatrix}
\hat{\beta}_1^{(2)}\\
\hat{\beta}_2^{(2)}
\end{bmatrix} = \boldsymbol{S}_\mathrm{dut}\boldsymbol{P}\begin{bmatrix}
\hat{\alpha}_1^{(2)}\\
\hat{\alpha}_2^{(2)}
\end{bmatrix}
\label{eq:22}
$$

Because the permutation matrix, when multiplied to a vector, has the effect of reordering the elements of the vector, we can absorb the permutation matrices into the vectors. This simply translates to reordering the elements of the vectors:

$$
\begin{bmatrix}
\hat{\beta}_2^{(2)}\\
\hat{\beta}_1^{(2)}
\end{bmatrix} = \boldsymbol{S}_\mathrm{dut}\begin{bmatrix}
\hat{\alpha}_2^{(2)}\\
\hat{\alpha}_1^{(2)}
\end{bmatrix}
\label{eq:23}
$$

Therefore, now we can combine the first equation with the above equation in a single matrix product as follows:

$$
\begin{bmatrix}
\hat{\beta}_1^{(1)} & \hat{\beta}_2^{(2)}\\
\hat{\beta}_2^{(1)} & \hat{\beta}_1^{(2)}
\end{bmatrix} = \boldsymbol{S}_\mathrm{dut}\begin{bmatrix}
\hat{\alpha}_1^{(1)} & \hat{\alpha}_2^{(2)}\\
\hat{\alpha}_2^{(1)} & \hat{\alpha}_1^{(2)}
\end{bmatrix} \quad \Longrightarrow \quad  \boldsymbol{S}_\mathrm{dut} = \begin{bmatrix}
\hat{\beta}_1^{(1)} & \hat{\beta}_2^{(2)}\\
\hat{\beta}_2^{(1)} & \hat{\beta}_1^{(2)}
\end{bmatrix}\begin{bmatrix}
\hat{\alpha}_1^{(1)} & \hat{\alpha}_2^{(2)}\\
\hat{\alpha}_2^{(1)} & \hat{\alpha}_1^{(2)}
\end{bmatrix}^{-1}
\label{eq:24}
$$

And this is how you perform a full correction of a DUT using a two-port/one-path VNA system.

> **Fun Fact:** If the DUT is symmetric such that S11=S22 and S21=S12, then from just one measurement, you can substitute the second measurement with the same original measurement, and you will get the answer directly (this is the fakeflip 😉).
> 
> 
> $$
> \boldsymbol{S}_\mathrm{dut}^{\mathrm{(sym)}} = \begin{bmatrix}
> \hat{\beta}_1^{(1)} & \hat{\beta}_2^{(1)}\\
> \hat{\beta}_2^{(1)} & \hat{\beta}_1^{(1)}
> \end{bmatrix}\begin{bmatrix}
> \hat{\alpha}_1^{(1)} & \hat{\alpha}_2^{(1)}\\
> \hat{\alpha}_2^{(1)} & \hat{\alpha}_1^{(1)}
> \end{bmatrix}^{-1}\label{eq:25}
> $$
> 

## References

[1] Scikit-rf: <https://scikit-rf.readthedocs.io/en/latest/examples/metrology/TwoPortOnePath%2C%20EnhancedResponse%2C%20and%20FakeFlip.html>

[2] METAS VNA Tools - Math Reference: <https://www.metas.ch/metas/en/home/fabe/hochfrequenz/vna-tools.html>

[3] Z. Hatab, M. E. Gadringer and W. Bösch, "Symmetric-Reciprocal-Match Method for Vector Network Analyzer Calibration," in *IEEE Transactions on Instrumentation and Measurement*, vol. 73, pp. 1-11, 2024, Art no. 1001911, doi: [10.1109/TIM.2024.3350124](https://doi.org/10.1109/TIM.2024.3350124)