---
layout: page
title: CMP++
description: A library for the Bayesian calibration of computer codes using sequential and adaptive approaches.
img: assets/img/projects/cmp/logo.png
importance: 1
category: phd
---

<div class="repo p-2 text-center">
  <a href="https://github.com/omarkahol/cmp" rel="external nofollow noopener" target="_blank">
    <img class="repo-img-light w-100" alt="omarkahol/cmp" src="https://github-readme-stats.vercel.app/api/pin/?username=omarkahol&amp;repo=cmp&amp;theme=default&amp;show_owner=false">
    <img class="repo-img-dark w-100" alt="omarkahol/cmp" src="https://github-readme-stats.vercel.app/api/pin/?username=omarkahol&amp;repo=cmp&amp;theme=dark&amp;show_owner=false">
  </a>
</div>

---
<header class="post-header">
            <h1 class="post-title">Technical documentation <a href="../../assets/pdf/projects/cmp/cmp.pdf" target="_blank" rel="noopener noreferrer" class="float-right"><i class="fas fa-file-pdf"></i></a></h1>
            <p class="post-description">The document containing the technical documentation for the library, install instructions and a test case.</p>
          </header>
---

CMP++ is a c++ library for the calibration of computer codes using Bayesian methods. It supports the use of a full bayesian solution, sequential approaches and modular approaches. In the following sections I will introduce different aspects of the calibration problem and show how the CMP method works.

The following section briefly discusses the Bayesian calibration framework.

$$
\newcommand{\cmp}{\text{CMP}}
\newcommand{\fmp}{\text{FMP}}
\newcommand{\koh}{\text{KOH}}
\newcommand{\map}{\text{MAP}}
\newcommand{\reals}{\mathbb{R}}
\newcommand{\btheta}{\boldsymbol{\theta}}
\newcommand{\bx}{\boldsymbol{x}}
\newcommand{\bt}{\boldsymbol{t}}
\newcommand{\bz}{\boldsymbol{z}}
\newcommand{\bY}{\boldsymbol{\mathcal{Y}}}
\newcommand{\bres}{\boldsymbol{R}_{\btheta}}
\newcommand{\bpsi}{\boldsymbol{\psi}}
\newcommand{\beps}{\boldsymbol{\epsilon}}
\newcommand{\balpha}{\boldsymbol{\alpha}}
\newcommand{\rr}{\bres^T \bres}
\newcommand{\bphi}{\boldsymbol{\phi}_{\boldsymbol{\theta}}}
\newcommand{\gpmean}{\mu_{\boldsymbol{\psi}_z}}
\newcommand{\gpcov}{c_{\boldsymbol{\psi}_z}}
\newcommand{\detf}[1]{\mid #1 \mid}
\newcommand{\muh}{\hat{\mu}}
\newcommand{\data}{\boldsymbol{\mathcal{D}}}
\newcommand{\model}{\mathcal{M}}
\newcommand{\yobs}{\boldsymbol{y}_{\text{obs}}}
\newcommand{\xobs}{\boldsymbol{x}_{\text{obs}}}
\newcommand{\Yobs}{\boldsymbol{\mathcal{Y}}_{\text{obs}}}
\newcommand{\kernel}{K_{\boldsymbol{\psi}}}
\newcommand{\summ}[3]{\displaystyle \sum_{#1=#2}^{#3}}
\newcommand{\dd}{\text{d}}
\newcommand{\intpsi}[1]{\displaystyle\int_{\Psi} \ #1 \ \text{d} \boldsymbol{\psi}}
\newcommand{\intphi}[1]{\displaystyle\int_{\Phi} \ #1 \ \text{d} \boldsymbol{\phi}}
\newcommand{\intt}[1]{\displaystyle \int_{\Theta} \ #1 \ \text{d} \boldsymbol{\theta}}
\newcommand{\likelihood}[1]{\mathcal{L}\left( #1\right)} 
\newcommand{\prior}[1]{\pi \left( #1\right)} 
\newcommand{\p}[1]{\text{p} \left( #1 \right)} 
\newcommand{\pq}[1]{\text{q} \left( #1 \right)} 
\newcommand{\normalpdf}{\mathcal{N}}
\newcommand{\approxp}[2]{\hat{\text{p}}_{#1} \left( #2 \right)}
\newcommand{\argp}[2]{\text{p}_{_{#1}} \left( #2 \right)} 
\newcommand{\psimap}{\boldsymbol{\psi}_{_\text{MAP}} \left( \boldsymbol{\theta} \right)} 
\newcommand{\phimap}{\boldsymbol{\phi}_{_\text{MAP}} \left( \boldsymbol{\theta} \right)} 
\newcommand{\psikoh}{\boldsymbol{\psi}_{_\text{KOH}}}
\newcommand{\Stheta}{S_{\btheta}}
\newcommand{\hessian}[1]{\mathcal{H}_{#1}} 
\DeclareMathOperator*{\argmax}{arg\,max}
\newcommand{\expectation}[2]{\mathbb{E}_{#1} \left[ \ #2 \ \right]}
\newcommand{\trace}[1]{\text{trace} \left( #1 \right)}
$$

## Bayesian calibration framework

Let $$y$$ be the output of a scalar computer model, $$\model$$ which depends on some coordinates, $$\bx \in X \subset \reals^{d_x}$$ and some parameters, $$\btheta \in \Theta \subset \reals^p$$ so that

$$
\begin{equation}
    y = \model(\bx, \btheta) \; .
\end{equation}
$$

The objective of model calibration is, given a set of n experiments $$\data = \left\{ \xobs^i, y^i \right\}_{i = 1 \ldots n}$$ to compute a plausible value (or distribution) of the model parameters.

To do so, one needs to define a proper statistical model that explains the observations. 
As standard in many calibration studies, we model the observations, $$\yobs = \left\{ y^i \right\}_{i = 1 \ldots n}$$ as a set of iid random variables, $$\Yobs$$ that are explained by the following statistical model:

$$
\begin{equation}
    \Yobs = \boldsymbol{\model} + \bz + \beps \;,
    \label{eq:stat_model}
\end{equation}
$$

where we have denoted with $$\beps$$ and $$\bz$$ two random vectors that are representatives of the measurement and model errors, respectively. The random vector $$\boldsymbol{\model}$$ contains the output of the computer model at each observation, $$\boldsymbol{\model} = \left\{\model(\xobs^i, \btheta) \right\}_{i=1 \ldots n}$$.
Measurement noise is considered to follow a normal distribution with zero mean and diagonal covariance, $$\beps|_{\sigma_e} \sim \normalpdf(0, \, \sigma_e^2 \ \mathbb{I}_n)$$. The model discrepancy term is usually modeled with a Gaussian Process $$ \bz |_{\bpsi_z} \sim \text{GP}(\gpmean,\gpcov)$$, where $$\gpmean$$ and $$\gpcov$$ are, respectively, the mean and covariance function of the process. It is suggested to have $$\gpmean(x) = 0$$ for better identifiability. We denote with hyper-parameters the variables $$\bpsi = (\bpsi_z, \sigma_e) \subset \reals^h$$ that describe experimental and model errors. In the literature, one can find different expressions of Eq. \ref{eq:stat_model} that can account also for data inconsistency or multiplicative experimental error. 


According to the frequentist perspective, $$\btheta, \bpsi$$ are deterministic variables and one usually constructs a suitable estimator. In the Bayesian perspective these variables make sense only as possible realizations of random variables, $$\boldsymbol{\Theta}, \boldsymbol{\Psi}$$. Under this framework, one tries to compute the posterior distribution, $$\p{\btheta, \bpsi \mid \data, \model}$$, which is defined as the joint PDF of $$\boldsymbol{\Theta}, \boldsymbol{\Psi}$$ conditioned on the realized experiments and the model choice. This can be done using Bayes' rule:

$$
\begin{equation}
    \p{\bpsi, \btheta | \data, \model} = \dfrac{\likelihood{\data \mid \btheta, \bpsi, \model } \ \prior{\bpsi, \btheta \mid \model}}{\p{\data \mid \model}} \; .
    \label{eq:posterior}
\end{equation}
$$


The term $$\prior{\bpsi, \btheta \mid \model}$$ is called prior and represents the knowledge of the distribution of the parameters before observing the data. The prior can be non-informative or informative and particular care must be taken to define a proper prior function since the results can be strongly affected by it. 


The term $$\likelihood{\data \mid \btheta, \bpsi, \model }$$ is called likelihood and, because of the assumptions made, it is a multivariate normal distribution:

$$
\begin{equation}
    \log \likelihood{\data \mid \btheta, \bpsi, \model } = - \dfrac{n}{2} \log 2 \pi - \dfrac{n}{2} \detf{\kernel + \sigma_e^2 \mathbb{I}_n} -\dfrac{1}{2} \bres^T (\kernel + \sigma_e^2 \mathbb{I}_n)^{-1} \bres \; ,
    \label{eq:ll}
\end{equation}
$$


where we indicated with $$\bres$$ the vector containing the residuals, $$\bres=\yobs - \boldsymbol{\model}$$, and with $$\kernel$$ the covariance matrix obtained from the structure of the model error terms, $$ (\kernel)_{i,j} = \gpcov (\xobs^i, \ \xobs^j)$$.

Finally, the denominator in Eq. \ref{eq:posterior} is usually called model evidence and appears as a normalization constant for the posterior distribution:

$$
\begin{equation}
    \p{\data \mid \model} = \intt{\intpsi{\likelihood{\data \mid \btheta, \bpsi, \model } \ \prior{\bpsi, \btheta \mid \model}}} \;.
    \label{eq:model_evidence}
\end{equation}
$$

Once the posterior PDF is available, it can be used to make predictions at a new coordinate, $$\bx^*$$. The calibrated model prediction is the probability distribution associated with the evaluation of the computer model at the new point, $$\mathcal{Y}_{\text{cal}}^* \sim \p{y_{\text{cal}}^* \mid \data}$$ and can be computed using

$$
\begin{equation}
    \p{y_{\text{cal}}^* \mid \data,\model} = \intt{\p{y_{\text{cal}}^* \mid \btheta, \data,\model} \ \p{\btheta \mid \data,\model}} \; .
    \label{eq:cal_pred}
\end{equation}
$$


Because the computer model is deterministic, $$\p{y_{\text{cal}}^* \mid \btheta, \data, \model} = \delta (y_{\text{cal}}^* - \model(\bx^*, \btheta))$$. In some cases the computer model, if too expensive to evaluate, can be substituted with a cheaper surrogate. The latter is usually a Gaussian Process so the associated uncertainty must be taken into account. In this derivation we used the parameters' posterior,

$$
\begin{equation}
    \p{\btheta \mid \data, \model} = \intpsi{\p{\btheta, \bpsi \mid \data, \model}} \;,
    \label{eq:par_posterior}
\end{equation}
$$


which is the marginal posterior distribution of the model parameters. 


The corrected model prediction regards the probability distribution associated with the prediction of the computer model and model error term at a new point, $$\mathcal{Y}_{\text{corr}}^* \sim \p{y_{\text{corr}}^* \mid \data, \model}$$. Similarly, the distribution can be computed using the rule of marginalization:  

$$
\begin{equation}
    \p{y_{\text{corr}}^* \mid \data,\model} = \intt{\intpsi{\p{y_{\text{corr}}^* \mid \btheta, \bpsi_z, \data,\model} \ \p{\btheta, \bpsi_z \mid \data,\model}}} \; .
    \label{eq:corr_pred}
\end{equation}
$$


In this case, the distribution $$\p{y_{\text{corr}}^* \mid \btheta, \bpsi_z, \data,\model}$$ is normal with mean and variance specified by the Gaussian process predictive equations: 

$$
\begin{equation}
    \left. 
        \begin{aligned}
            \mathcal{Y}_{\text{corr}}^*|_{\btheta,\bpsi_z,\data,\model} \sim& \normalpdf (\mu_\text{pred}, \sigma_\text{pred}^2) \; ,\\
            &\mu_\text{pred} = \model (\bx^*, \btheta) + k_*^T \ \kernel^{-1} \ \bres \;,\\
            &\sigma_\text{pred}^2 = \kernel - k_*^T \ \kernel^{-1} \ k_*^T \;.
        \end{aligned}
    \right.
\end{equation}
$$


Where $$(k_*^T)_i = \gpcov(\bx^*,\xobs^i)$$. 
We remark that, in both cases, the experimental noise random variable can be included if the prediction required concerns the value of a new experiment. 