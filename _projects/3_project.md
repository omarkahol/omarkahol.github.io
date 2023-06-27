---
layout: page
title: BSc Thesis
description: Parallel computing techniques for advection diffusion probles
img: 
importance: 2
category: thesis
---

The goal of this project was to solve the advection diffusion problem using parallele computing techniques

$$
\frac{\partial u}{\partial t} + c \nabla u = k \Delta u
$$

The problem was solved using the FD method and two different parallelization techniques

<ul>
    <li>Shared memory systems using OpenMP</li>
    <li>Distributed memory systems using MPI</li>
</ul>

The codes were implemented in C++ in both 1 and 2 spacial dimensions.

The FE approach, implemented in the deal.ii library, was then used to solve the problem in 2D and 3D using also mesh adaptation techniques. 

{% include figure.html path="assets/img/bscthesis/3d.png" title="solution" class="img-fluid rounded z-depth-1" %}

You can have a look at the thesis, which is in italian.

---
<embed src="../../assets/pdf/projects/bscthesis/tesi.pdf" width="100%" height="1000"> 