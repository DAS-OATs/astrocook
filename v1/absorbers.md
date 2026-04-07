---
layout: default
title: Modeling absorbers
parent: Common tasks
nav_order: 3
math: mathjax2
---

# Modeling absorbers
{: .no_toc}

Astrocook includes several tools to model the absorption features observed in quasar spectra. Absorbers are modeled with composite [Voigt profiles](https://en.wikipedia.org/wiki/Voigt_profile), which are fitted to the data using the [Levenberg-Marquardt algorithm](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) for least-squares minimization, as implemented by the [Lmfit](https://lmfit.github.io/lmfit-py/index.html) package.

To model absorbers, the spectrum must be first normalized to [continuum](continuum.md). 🚧

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---

## Voigt-profile modeling

The [Voigt profile](https://en.wikipedia.org/wiki/Voigt_profile) is the convolution of a Gaussian and a Cauchy-Lorentz distribution. In the context of line fitting, it is typically parametrized as a function of three quantities that describe the absorbing medium: its redshift $$z$$, its column density $$N$$, and its Doppler broadening $$b$$ (due to the thermal and/or turbulent motions of the particles). The profile of an absorption line as a function of the wavelength $$\lambda$$ is thus $$e^{-\tau_\lambda}$$, where the opacity $$\tau_\lambda$$ is computed as

$$\tau_\lambda=N\frac{\sqrt{\pi}e^2}{m_ec}\frac{f}{\Delta\nu_b}V(a,u_\lambda),$$

with $$e$$ the electron charge, $$m_e$$ the electron mass, $$c$$ the speed of light, and $$f$$ the oscillator strength of the transition producing the line. The Voigt function $$V(a,u)$$ id defined as

$$V(a,u_\lambda)=\frac{a}{\pi}\int_\infty^\infty\frac{e^{-y^2}}{a^2+(u_\lambda-y)^2}dy,$$

where $$a$$ and $$u_\lambda$$ depend on $$z$$ and $$b$$ as follows:

$$a=\frac{\Gamma}{4\pi b}\frac{\lambda_\mathrm{obs}}{1+z},$$

$$u_\lambda=\frac{c}{b}\left(\frac{\lambda}{\lambda_\mathrm{obs}}-1\right),$$

with $$\Gamma$$ the transition damping constant and $$\lambda_\mathrm{obs}$$ the observed wavelength of the line.

Voigt-profile modeling assumes that an arbitrarily complex absorbers can be analyzed as a superposition of finitely-many discrete components, each one characterized by a single set of parameters $$(z,N,b)$$. This is of course an oversimplification. The decomposition is not always unique, as the parameters may be degenerate across different components and even within a single component ($$N$$ and $$b$$ are degenerate when a line is saturated). It is generally possible to find a composite Voigt profile which is fitting to the data and complies with a chi-squared test or a given information criterion, but its physical interpretation is not always straightforward.

In Astrocook, each component in a Voigt-profile fit is called an *absorption system* (or simply a *system*) and appears as a separate entry in the corresponding [list](structures.md#list-of-absorption-systems). In the list, the column `series` specifies if each system is a single line (e.g. `Ly_a` for Lyman-alpha) or a multiplet of lines (e.g. `Ly`, including all transitions in the Lyman series, or `CIV`, including the two members of the C <span style="font-variant:small-caps;">iv</span> doublet at 154.8204 and 155.0781 nm). When a system is a multiplet, it includes several Voigt profiles centered at different wavelengths, and with different equivalent widths depending on the oscillator strength of the multiplet members. All these profiles are nevertheless described by a single set of parameters $$(z,N,b)$$.

[Here](series.md) you can find a complete list of the ionic transitions used to model absorption features.

## An interactive example 🚧

## Automatize the workflow 🚧
