# Quadrupole Lab — RF/DC & Mathieu Simulator

Hi there :)

This is a personal project dedicated to understanding mass spectrometers, in particular quadrupoles and ion traps.

This website allows physics calculations of an ideal linear quadrupole mass filter. No server, dependencies, build tooling, telemetry, or API keys. Open `index.html` directly or deploy the folder to GitHub Pages.

I will updating and editing this over time. Feel free to have a browse.

**Features**

- RF and DC operating voltages from target m/z, frequency, field radius, and a target **ideal** mass-window width.
- Manual Mathieu `(a, q)` operating point, or RF-only ion-guide mode.
- Mathieu region-I stability diagram calculated from **fourth-order Runge–Kutta** (RK4) Floquet-monodromy integration (not a pre-drawn curve).
- Individual transverse ion trajectories for finite quadrupole length, entrance RF phase, kinetic energy, entrance offset and velocity.
- RF/DC voltage profiles versus m/z, with exportable CSV.

## Quick start / GitHub Pages

1. Create a **public** GitHub repository and upload `index.html`, `app.js`, `quad-math.js`, and `style.css` to the repository root. Upload `README.md` if you want the documentation in GitHub.
2. Open the repository's **Settings → Pages** and choose **Deploy from a branch → main / (root)**.
3. After GitHub Pages is published, use the URL shown in Pages settings, usually `https://USERNAME.github.io/REPOSITORY/`. No server code or Node runtime is needed for visitors.

Run numerical self-tests locally with `node test.js` (Node is only needed for testing, **not** for the website).

## Voltage conventions — important

The code assumes the opposite rod pairs are driven at:

\[
 +[U + V\cos(\Omega t)],\qquad -[U + V\cos(\Omega t)]
\]

where `U` is DC voltage per rod pair relative to the centre potential and `V` is **RF zero-to-peak per rod pair**; \(\Omega=2\pi f\), and the inscribed field radius is \(r_0\). For an ideal hyperbolic quadrupole:

\[
 \Phi(x,y,t)=[U+V\cos(\Omega t)]\frac{x^2-y^2}{r_0^2}.
\]

The dimensionless Mathieu parameters are:

\[
 a=\frac{8 e U}{(m/z)\,u\,r_0^2\Omega^2},\qquad
 q=\frac{4 e V}{(m/z)\,u\,r_0^2\Omega^2}.
\]

Here `(m/z)` is a numeric mass-to-charge ratio in thomsons (`u/e`), `e` is elementary charge and `u` is the atomic mass constant. These equations work for multiply charged ions when `m/z` is represented as `u/e`. A different common convention describes the **differential** rod-pair RF amplitude and DC voltage, with Mathieu prefactors **4 and 2** instead. Do not interchange conventions.

Display definitions:

| Display | Meaning |
|---|---|
| RF amplitude `V` | Zero-to-peak relative to 0 V, on one rod pair |
| DC voltage `U` | DC bias relative to centre potential, on one rod pair |
| Per-rod RF peak-to-peak | `2 V` |
| Differential RF peak-to-peak | `4 V` between the opposite rod pairs |
| Differential DC potential | `2 U` between opposite rod pairs |

## How RK4 numerical stability works

Use dimensionless time \(\tau=\Omega t/2\). In canonical Mathieu form:

\[
 u''+[a-2q\cos(2\tau)]u=0.
\]

Numerically integrate two linearly independent initial solutions `(u, u') = (1, 0)` and `(0, 1)` across the coefficient period `π`. Their final values form the **Floquet monodromy matrix**. For a second-order Mathieu equation, the determinant of that matrix is unity; a bounded, strictly stable solution requires `|trace(M)| < 2`. The transverse x and y directions respectively have `(a, q)` and `(-a, -q)`; changing the sign of `q` is equivalent to shifting RF phase by π, so `trace(-a, q)` may be used for y.

The region-I boundary is calculated by root-finding around the first periodic/antiperiodic characteristic values:

- x-axis boundary: `trace(M(a, q)) = -2`.
- y-axis boundary: `trace(M(-a, q)) = +2`.

The intersection is the apex, expected at approximately **q = 0.705996, a = 0.236994**. The program obtains it numerically, rather than hardcoding the point.

### What target Δ(m/z) means

For constant `a/q` there are two intersections along a line in the first stability region: `qLow`, `qHigh`. Because `q ∝ 1/(m/z)` for fixed rod voltages, the ideal transmitted range for a selected `qCentre` is:

\[
 (m/z)_{low} = (m/z)_{target}\frac{q_{centre}}{q_{high}},\quad
 (m/z)_{high} = (m/z)_{target}\frac{q_{centre}}{q_{low}}.
\]

The target sits at the arithmetic midpoint of the **mass interval** when `qCentre = 2 qLow qHigh / (qLow + qHigh)` (harmonic midpoint of the q limits). A bisection search finds the operating line slope whose ideal mass window equals the requested width. This is a theoretical, infinite-duration stability window, **not** the 10%-height or FWHM of a measured transmission peak.

RF-only mode is an ion guide without a finite upper mass cutoff; it must not report a narrow isolation window.

### Finite-length trajectory

The trajectories use a physically equivalent `+cos`-phase form (see above). User inputs specify a single ion, not an ensemble. The axial drift velocity comes from the supplied **axial energy per charge** `eV/e`, and the number of RF periods sampled is `f·L/vz`. One ion is labelled *transmitted* if its radial displacement remains within the conservative circular acceptance radius `r0`. The true rod aperture differs from a circle, and this is a simplified proxy. A single trajectory cannot predict transmission probability or instrumental peak width.

## Limitations / validation plan

This is an ideal **two-dimensional** Mathieu model plus constant axial speed. It does not solve Laplace/Poisson fields for cylindrical rather than hyperbolic rods; it omits fringe fields, ion–ion Coulomb interactions, RF loading, gas collisions, manufacturing tolerances, space-charge, measured abundance sensitivity, and detector/ion-source transfer functions. Boundary curves are determined by Floquet theory; finite-length trajectories have a much broader acceptance that depends on injection position, angle and phase. The code should not be used to directly specify high-voltage hardware without independent electrical and mechanical design review.

For a more physical model, integrate a measured or SIMION-calculated 3-D electric-field map, launch a distribution of ions with phase/energy/angle spreads, and compute transmission vs m/z; benchmark against known instrument data.

**Checks:** Run `node test.js` for the apex (first region), algebraic RF/DC roundtrip, prescribed passband width, and a known unstable lower-mass ion. Suggested extra validation: grid refinement study (96 vs 192 RK4 steps/period), phase-space convergence and comparing resolution predictions with experimental data.

## References

- Tunstall, J. J., Voo, A. C. C. & Taylor, S. *Computer Simulation of the Mass Filter for a Finite Length Quadrupole*, **Rapid Communications in Mass Spectrometry** 11 (1997), 184–188. [DOI](https://doi.org/10.1002/(SICI)1097-0231(19970131)11:2%3C184::AID-RCM724%3E3.0.CO;2-P); [open PDF](https://www.liverpool.ac.uk/media/livacuk/massspectrometry/pdfs/Computer%2Csimulation%2Cof%2Cthe%2Cmass%2Cfilter%2Cfor%2Ca%2Cfinite%2Clength%2Cquadrupole.pdf).
- GAA Custom Electronics, [quadrupole RF/DC reference calculator](https://www.gaacustom.com/tools.html). Similar feature scope; **this project implements its own equations/solver/UI**.
- Giraud / BYU, [*Simulation and Manufacture of a Quadrupole Mass Filter*](https://marvin.byu.edu/beryllium7/publications/KellenGiraudSeniorThesis2008.pdf).

## Licence

MIT — see `LICENSE`. If publishing, keep the attribution in the licence; replace the project name and add your own name for ownership/attribution as you wish.
