---
layout: page
title: HyShock++
description: A solver for the chemistry relaxation equations, in high-temperature post-shock conditions.
img: assets/img/projects/hyshock++/thumbnail.png
importance: 2
category: courses
---

The goal is to study the effect of thermochemical non equilibrium on the relaxation region behind a normal shock wave. The solver was written in python and made extensive use of the Mutation++ library.
Three different models were adopted:

<ul>
    <li>Thermochemical equilibrium, the shock wave is approximated by a step function, the relaxation region has zero thickness.</li>
    <li>Thermodynamic equilibrium, the gas is assumed to be well described by 1 temperature, finite rate chemistry acts on the flow.</li>
    <li>Themochemical non-equilibrium, the internal vibrational degrees of freedom are not in equilibrium with the others. The two temperature approximation is used. </li>
</ul>

#### Equations
s
In the most general case, the solver solves the following steady state equations, which correspond to the 1D euler equations for a reacting mixture. 

$$
    \begin{cases}
        \dfrac{\partial}{\partial x} [\rho u] = 0\\
        \dfrac{\partial}{\partial x} [\rho_i u] = \dot{\omega}_i(\rho_i,\rho e, \rho e^v)\\
        \dfrac{\partial}{\partial x} [\rho u u + P(\rho_i,\rho e, \rho e^v)]= 0\\
        \dfrac{\partial}{\partial x} [\dfrac{1}{2} \rho u^3 + \rho ue + uP(\rho_i,\rho e, \rho e^v)] = 0\\
        \dfrac{\partial}{\partial x} [\rho u e^v] = \Omega_v(\rho_i,\rho e, \rho e^v)\\
        e_i = e_i^{tr}+e_i^{int}
    \end{cases}
$$

#### Solver 

Since the equations have the following form

$$
    \frac{\partial}{\partial x} \textbf{f}(P, Y_i, u, T, T_v) = \textbf{s}(P, Y_i, u, T, T_v)
$$

We can integrate the equations in the numerical domain to get a system of nonlinear equations which can then be solved numerically. Note that this process is very fast and usually requires very few iterations if the solution at the previous cell is used as initial guess. 

#### Code and presentation

The code was tested using different mixtures and post-shock conditions and is available on GitHub

<div class="repo p-2 text-center">
  <a href="https://github.com/omarkahol/hyShockpp" rel="external nofollow noopener" target="_blank">
    <img class="repo-img-light w-100" alt="omarkahol/hyShockpp" src="https://github-readme-stats.vercel.app/api/pin/?username=omarkahol&amp;repo=hyShockpp&amp;theme=default&amp;show_owner=false">
    <img class="repo-img-dark w-100" alt="omarkahol/hyShockpp" src="https://github-readme-stats.vercel.app/api/pin/?username=omarkahol&amp;repo=hyShockpp&amp;theme=dark&amp;show_owner=false">
  </a>
</div>

You can also have a look at the project presentation, which includes all the relevant details and derivations

<embed src="../../assets/pdf/projects/hyshock++/Hypersonic_Flows.pdf" width="100%" height="1000"> 