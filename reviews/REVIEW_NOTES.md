I reviewed both `894b0aa` and the newer `1f17671`, using isolated snapshots. **I have not modified the shared checkout.**

The additions are worthwhile. More importantly, I found **two further analytical improvements**: a cancellation that sharpens the log-map decomposition, and an evaluated leading-octupole spectrum extending the existing quadrupole result.

Full derivations and isolated verification scripts are recorded in:
[`REVIEW_NOTES.md`](/tmp/temperature-review-1f17671/REVIEW_NOTES.md).

## 1. Review of the latest additions

### Entire-function transfer formula: retain it

This is a genuine improvement over the older expanded \(\mathrm{Ci}/\mathrm{Si}\) expressions. The separation into two endpoint principal parts is sound, and regularity explains why the subtraction terms disappear after taking the real part.

I reproduced the exact checks for **\(\ell=1,\ldots,6\)**, without numerical quadrature.

### Log-map decomposition: correct and useful

The decomposition
\[
C_\ell^{(4)}
=C_\ell^{L,4}+\mathcal N_\ell+\mathcal M_\ell
\]
correctly separates the cubic log-map response from the universal normalization contribution and the mixed quadratic contribution. In particular,
\[
-\frac{C_\ell^2}{2\pi}
\]
is the correct random-sky-monopole correction.

The supplied test initially exceeded my time limit. Using the **same exact test model**, but truncating the formal \(\epsilon\)-polynomial before each multiplication, all checks for \(\ell=1,\ldots,4\) passed. This suggests a practical improvement to the verification code: avoid expanding orders that will immediately be discarded.

Two wording qualifications would improve the report:

- “All dynamics sits in the log-map spectrum” is too strong: \(\mathcal M_\ell\) also depends on the quadratic dynamics. **All cubic response** is confined to \(C_\ell^{L,4}\).
- The general \(\mathcal N_\ell\) expression is an explicit convergent harmonic series, rather than a finite closed expression.

A compact alternative to its double harmonic sum is
\[
\boxed{
\mathcal N_\ell
=\pi\int_{-1}^{1}P_\ell(u)\,\xi(u)^2\,du
-\frac{C_\ell^2}{2\pi},
\qquad \ell\ge1.
}
\]
This replaces the infinite sum by one angular integral; it does not evaluate that integral.

### Shell locality: useful, but distinguish the two angular integrations

The shell identities and second-order shear formula passed the exact checks. Isolating \(A_2\) as the only genuinely second-order field response entering \(L_2\) is useful.

However, removing the **primordial relative-momentum angle** does not automatically remove the **viewing-angle Legendre projection**. The real-space formulation should explicitly retain that projection unless it has been evaluated.

### Dilation factorization: correct algebra, interpretation needs qualification

The factorization into \(e^{8\epsilon^2m_0}\) and the remaining polynomial is correct. But describing the exponential as merely chart normalization, or treating the dilation-factored ratio as a replacement “physical answer,” goes too far.

For a genuine constant relabelling,
\[
X'=e^cX,\qquad
\zeta'(X')=\zeta(e^{-c}X')-c,
\]
we have
\[
\mu'=\mu-c,\qquad m_j'=e^{-jc}m_j,\qquad m_0'=m_0.
\]
The full spectrum is invariant, but **the stochastic factor \(e^{8\epsilon^2m_0}\) remains**.

Writing the answer using \(\langle\Xi_{\rm phys}^2\rangle\) is a valid factorization. Holding a locally inferred scale fixed, however, is a different conditional specification from holding the primordial chart-shell \(k_0\) fixed.

**For the ensemble specified in the problem, the correction remains \(29/3\), not \(5/3\).** The latter is a useful dilation-factored ratio.

---

## 2. New cancellation: \(\mathcal M_\ell\) starts at \(O(\kappa^8)\)

The latest report gives the power-counting estimate
\[
\mathcal M_\ell=O(\kappa^6).
\]
For these Gaussian initial conditions, that leading contribution actually vanishes.

Write
\[
A_{ij}=(\partial_i\partial_jg)_{\rm TF},
\qquad
B_{ij}=(\partial_i g\,\partial_jg)_{\rm TF}.
\]
At leading gradient order,
\[
\theta(\hat n)=-\tau^2 A_{ij}n^in^j,
\qquad
L_2(\hat n)=\tau^2(4gA_{ij}+2B_{ij})n^in^j.
\]

The Gaussian variables \(g\), \(A\), and \(\nabla g\) are mutually independent. Moreover,
\(\langle g\rangle=0\) and \(\langle B\rangle=0\). Therefore
\[
\left\langle
\theta(\hat n)\theta(\hat n')L_2(\hat m)
\right\rangle_{O(\kappa^6)}=0
\]
for arbitrary directions. This kills both collapsed configurations entering \(\mathcal M_\ell\). Spatial parity eliminates the next, odd-gradient order.

Thus, under the finite-moment gradient expansion,
\[
\boxed{
\mathcal M_\ell=O(\kappa^8),
\qquad
\mathcal N_\ell=O(\kappa^8).
}
\]

**The fourth-order-amplitude temperature and log-temperature spectra agree through total gradient order six—not just order four.**

I checked the cancellation by exact Wick contraction of every independent Cartesian tensor component. This is a genuine cancellation, not a change of notation.

---

## 3. New evaluated result: the leading octupole

The curvature-quadrupole result extends particularly cleanly to the next multipole.

Define
\[
\beta_3=
\frac{(\bar t_o-\bar t_\star)^2(\bar t_o+2\bar t_\star)}{18}.
\]
Then the leading octupole is
\[
\boxed{
\Theta_{3,\mathrm{lead}}(\hat n)
=\beta_3
\big[\nabla_{\langle a}R_{bc\rangle}\big]_o
n^an^bn^c,
}
\]
where curvature, covariant differentiation, and orthonormal components refer to the **seed metric**. Angular brackets denote the fully symmetric trace-free part.

This has an elegant conformally Cartesian expression. With \(f=e^{-\zeta}\),
\[
\boxed{
\nabla_{\langle a}R_{bc\rangle}
=\frac f2
\left[\partial_a\partial_b\partial_c(f^2)\right]_{\rm STF}.
}
\]

It complements the existing quadrupole:
\[
\Theta_{2,\mathrm{lead}}
=\tau^2f\,(\partial_a\partial_bf)_{\rm TF}n^an^b.
\]

### Why the octupole formula holds

In seed Riemann normal coordinates, the leading shear is
\[
\sigma_{ab}=-\frac t3 R_{{\rm TF},ab}.
\]
Expanding curvature along the radial ray gives
\[
\Theta(\hat n)
=\frac13\int_{\bar t_\star}^{\bar t_o}
t\left[
R_{{\rm TF},ab}
+(\bar t_o-t)n^c\nabla_cR_{{\rm TF},ab}
\right]n^an^b\,dt
+O(\kappa^4).
\]
The first term is the quadrupole. The second splits into an octupole and a dipole; the latter is removed by the prescribed final dipole projection. Its time integral is precisely \(\beta_3\).

### Its Gaussian spectrum also resums

For finite moments \(m_0,m_2,m_4,m_6\), the resulting leading-gradient spectrum is
\[
\boxed{
C_{3,\mathrm{lead}}
=
\frac{16\pi}{1225}\beta_3^2\epsilon^2
e^{18\epsilon^2m_0}
\left[
m_6
+28\epsilon^2m_2m_4
+\frac{560}{9}\epsilon^4m_2^3
\right].
}
\]

On the monochromatic shell, with \(v=\epsilon^2\sigma^2\),
\[
\boxed{
C_{3,\mathrm{lead}}
=
\frac{4\pi}{3675}
(x_o-x_\star)^4(x_o+2x_\star)^2
\epsilon^2\sigma^2\,
e^{18v}\left(1+28v+\frac{560}{9}v^2\right).
}
\]

Through the requested fourth amplitude order,
\[
\boxed{
C_{3,\mathrm{lead}}
=\epsilon^2C_3^{(11)}
\left[1+46\epsilon^2\sigma^2\right]
+O(\epsilon^6).
}
\]

The separately checked coefficients satisfy
\[
\frac{C_3^{(22)}}{C_3^{(11)}}=148\sigma^2,
\qquad
\frac{C_3^{(13)}}{C_3^{(11)}}=27\sigma^2,
\]
so \(148/4+27/3=46\).

**Scope:** this evaluates the leading \(O(k_0^6)\) octupole power, not the arbitrary-wavelength spectrum. The exponential resums the leading-gradient expression only.

Checks passed for the tensor identity, Gaussian contractions, angular normalization, linear-transfer limit, and a separate leading-gradient ray calculation in the original coordinates. The full existing ray-generator cross-check timed out; I have not represented it as completed.

---

## 4. A compact leading transfer coefficient for every multipole

For every \(\ell\ge2\),
\[
\boxed{
T_\ell(k)=
-\frac{
k^\ell(\bar t_o-\bar t_\star)^{\ell-1}
[\bar t_o+(\ell-1)\bar t_\star]
}{
3(2\ell+1)!!
}
+O(k^{\ell+2}).
}
\]
This follows directly from the lowest powers of \(s'\) and \(w_\ell\), and supplies a simple boundary check on the more elaborate entire-function expression. I verified it symbolically for \(\ell=2,\ldots,8\).

## Recommendation

The strongest next additions are:

1. **Sharpen \(\mathcal M_\ell\) to \(O(\kappa^8)\).**
2. **Add the curvature-gradient octupole and its evaluated spectrum.**
3. **Qualify the “physical scale” interpretation without altering the correct factorization.**

These improve the actual results. The finite-time factorizations remain useful, but I would not yet claim that the general comparable-momentum fourth-order spectrum has been reduced to the compact evaluated form the problem author may be hoping for.
