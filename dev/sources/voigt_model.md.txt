# The Voigt model

This page describes how Astrocook computes the model of an absorption
spectrum and its derivatives. The same model serves the Voigt fit, the
evaluation of alternative models, the profiles drawn in the System
Inspector and the model column of the spectrum. The code is in
`astrocook/fitting/voigt_fitter.py` and `astrocook/fitting/lsf.py`.

## Optical depth of a line

A component has a redshift $z$, a column density $N$ (fitted as
$\log_{10} N$, with $N$ in cm$^{-2}$), a Doppler parameter $b$ and an
optional turbulent term $b_\mathrm{tur}$ (both in km/s). The line width
uses $b_\mathrm{eff} = \sqrt{b^2 + b_\mathrm{tur}^2}$; below, $b$ stands
for $b_\mathrm{eff}$. Each transition of the component, with rest
wavelength $\lambda_0$ (Å), oscillator strength $f$ and damping constant
$\Gamma$ (s$^{-1}$), has

$$
\lambda_c = \lambda_0 (1 + z), \qquad
\Delta_D = \frac{b}{c}\,\lambda_c, \qquad
a = \frac{\Gamma \lambda_0}{4 \pi b}, \qquad
\tau_0 = 1.4974 \times 10^{-15}\,\frac{N f \lambda_0}{b},
$$

and its optical depth at wavelength $\lambda$ is

$$
\tau(\lambda) = \tau_0\, H(x, a), \qquad
x = \frac{\lambda - \lambda_c}{\Delta_D}, \qquad
H(x, a) = \mathrm{Re}\, w(x + i a),
$$

where $w$ is the Faddeeva function (`scipy.special.wofz`). The damping
parameter $a$ uses the rest wavelength, so that it does not depend on
$z$. A line is computed only where its optical depth is above
$10^{-6}$ (`TAU_CUT`). The optical depths of all the lines are summed, and
the transmitted flux is $F = e^{-\tau}$.

## Instrumental profile and pixels

The optical depth is computed on a fine grid, uniform in velocity, that
extends beyond the pixels by five sigma of the instrumental profile. Its
step is the smallest of:

- the narrowest pixel;
- half the sigma of the instrumental profile;
- a third of the $b$ of the narrowest line in the pixels, but never below
  0.5 km/s.

Each pixel $i$, with edges $[v_{\mathrm{lo},i}, v_{\mathrm{hi},i}]$ in
velocity, gets the mean over its own width of the fine flux convolved
with a Gaussian of the resolving power *of that pixel* (the `resol`
column of the spectrum):

$$
M_i = \frac{\int F(v)\, K_i(v)\, dv}{\int K_i(v)\, dv}, \qquad
K_i(v) = \Phi\!\left(\frac{v_{\mathrm{hi},i} - v}{\sigma_i}\right)
       - \Phi\!\left(\frac{v_{\mathrm{lo},i} - v}{\sigma_i}\right),
$$

with $\Phi$ the Gaussian cumulative distribution and
$\sigma_i = c / (R_i\, 2\sqrt{2 \ln 2})$. On the fine grid the integral is
a weighted sum, $M_i = \sum_j W_{ij} F_j$, with weights normalised over
$j$. A component with a resolution of its own is convolved apart, and the
model is the product of the convolved flux of each resolution:
$M_i = \prod_g \sum_j W^{(g)}_{ij} F^{(g)}_j$.

## The fit

The fit minimises $\sum_i r_i^2$ with $r_i = (y_i - M_i) / \sigma_{y,i}$
on the pixels of the fitted lines, with `scipy.optimize.least_squares`
(trust-region reflective, bounded). Linked parameters (for example the
$z$ of the lines of a system, or a $b$ tied to another) are functions of
the free ones, written as expressions.

### The derivatives

The Jacobian $\partial r_i / \partial p_k$ is analytic. Until November
2026 it was computed by finite differences, one model per free
parameter at each iteration, so a fit of $k$ free parameters cost about
$k$ models per iteration.

The derivative of the Faddeeva function is

$$
w'(\zeta) = -2 \zeta\, w(\zeta) + \frac{2 i}{\sqrt{\pi}},
$$

so that, with $w = u + i v$ at $\zeta = x + i a$,

$$
\frac{\partial H}{\partial x} = \mathrm{Re}\, w' = -2 (x u - a v), \qquad
\frac{\partial H}{\partial a} = -\mathrm{Im}\, w' = 2 (x v + a u) - \frac{2}{\sqrt{\pi}}.
$$

Since $\tau_0 \propto N / b$, $a \propto 1 / b$, $x \propto 1 / b$, and
$x = \lambda c / (b \lambda_0 (1 + z)) - c / b$, the derivatives of the
optical depth of a line are

$$
\frac{\partial \tau}{\partial \log_{10} N} = \ln(10)\, \tau, \qquad
\frac{\partial \tau}{\partial z} = -\tau_0\, \frac{\partial H}{\partial x}\,
    \frac{x + c / b}{1 + z}, \qquad
\frac{\partial \tau}{\partial b} = -\frac{\tau_0}{b}
    \left(H + x \frac{\partial H}{\partial x} + a \frac{\partial H}{\partial a}\right),
$$

with $\partial b_\mathrm{eff} / \partial b = b / b_\mathrm{eff}$ and
$\partial b_\mathrm{eff} / \partial b_\mathrm{tur} = b_\mathrm{tur} / b_\mathrm{eff}$
for the two Doppler terms. The convolution is linear, so the derivatives
of the flux are convolved with the same weights as the flux:

$$
\frac{\partial M^{(g)}_i}{\partial p} = -\sum_j W^{(g)}_{ij}\, F^{(g)}_j\,
    \frac{\partial \tau^{(g)}_j}{\partial p},
$$

the groups of resolution are combined by the product rule, and the
linked parameters by the chain rule (the gradient of their expressions,
by central differences). Finally
$\partial r_i / \partial p = -(\partial M_i / \partial p) / \sigma_{y,i}$.

Two approximations, both far below the precision of a fit:

- the cut of each line at $\tau = 10^{-6}$ moves with its parameters; the
  derivative of the cut (samples entering or leaving it) is left out;
- the model takes $b$ no smaller than 0.1 km/s; below it the derivative
  in $b$ is set to zero (the fit keeps $b \geq 1$ km/s).

### How it was checked

Against central differences with small steps, the analytic Jacobian
agrees to $10^{-7}$–$10^{-9}$ of each column (the tests require
$10^{-5}$): at R = 45000 and 8900, without a resolution, with a
resolution per pixel, with two resolutions, with linked parameters and a
free turbulent term, on a damped Ly-$\alpha$ line (log N = 20) and with
fixed components around (`tests/test_voigt_jacobian.py`). The fits end at
the same minimum as with finite differences: on the mock benchmark the
counts are identical and the statistics of the pulls agree within
$4 \times 10^{-4}$. A Jacobian costs about three models whatever the
number of free parameters: a blend of 40 Ly-$\alpha$ lines (120 free
parameters) is fitted in 2.8 s instead of 29 s.
