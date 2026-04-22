# Helper library for computing F-harmonic forms

The following description can be more detailed in [docs/F-Harmonic-3-forms.pdf]

## Introduction 

Let $`(M, /omega)`$ be a 6-dimensional symplectic Lie group, let $/frak(g)$ its lie algebra, we will compute its Lie algebra cohomology and symplectic cohomology analogues. In general, for a 6-dimensional symplectic manifold, from any 3-form $/phi$ on M, or better any $/phi /in /bigwedge^3 g^*$, one can construct an endomorphism $K(/phi)$ of TM, another 3-form $F(/phi)$ and a scalar function $Q(/phi)$. The endomorphisms $F,K,Q$ are homogeneous polynomials in $/phi$ of degree 2, 3 and 4, respectively.

We say $/phi$ is $F$ harmonic if $d/phi=0$ and $dF(/phi)= 0$. We want to characterize the set of $F$ harmonic 3-forms. Also, we want to consider the following evolution equation of 3-forms $/del_t /phi = d/Lambda_/omega d F(/phi)$, where $/Lambda_/omega$ is the Lefschetz operator of contraction with respect to the symplectic form $/omega$. We want to verify the short-time existence and uniqueness.

## Background

Let \(G\) be a Lie group of real dimension \(n\). All the left-invariant data on \(G\) are determined by their value at the identity, hence they can be viewed as an element of the Lie algebra \(\mathfrak{g}\) or other algebraic constructions built from \(\mathfrak{g}\). After fixing a basis \(\{e_1, e_2, \ldots, e_n\}\) of \(\mathfrak{g}\), the isomorphism class of \(\mathfrak{g}\) is determined by the structure constants \(\{c_{ij}^{k}\}_{1\leq i,j,k\leq n}\), which are defined by \([e_i,e_j] = c_{ij}^k e_k\). In many cases, it is more convenient to work with the dual Lie algebra \(\mathfrak{g}^*\) and its exterior product

\[
\bigwedge^{*}\mathfrak{g}^{*} = \bigoplus_{i = 0}^{n}\bigwedge^{i}\mathfrak{g}^{*}.
\]

The structure of \(\mathfrak{g}\) can be recovered from the exterior derivative \(\mathrm{d}\) defined on \(\mathfrak{g}^{*}\), given by

\[
\mathrm{d}e^{k} = -\frac{1}{2} c_{ij}^{k}e^{i}\wedge e^{j},
\]

where \(\{e^{1},\ldots,e^{n}\}\) and \(c_{ij}^{k}\) are the dual basis and the structure constants associated to the basis \(\{e_{1},e_{2},\ldots,e_{n}\}\). Extending the above defined \(\mathrm{d}\) by linearity and the graded Leibniz rule, one gets a differential \(\mathrm{d}:\bigwedge^{i}\mathfrak{g}^{*}\to \bigwedge^{i + 1}\mathfrak{g}^{*}\) satisfying \(\mathrm{d}^{2} = 0\). The cohomology associated to the cochain complex \((\bigwedge^{*}\mathfrak{g}^{*},\mathrm{d})\) is known as the Lie algebra cohomology of \(\mathfrak{g}\). In many cases, the Lie algebra cohomology of \(\mathfrak{g}\) coincides with the de Rham cohomology of quotients of \(G\) by its cocompact discrete subgroups, including the case when \(G\) is compact or nilpotent.

For the convenience of notation, we shall use a symbol like \((0,0,0,e^{15},0,e^{13})\) to denote the 6-dimensional Lie algebra \(\mathfrak{g}\) defined by a basis \(\{e^{1},\ldots,e^{6}\}\) of \(\mathfrak{g}^{*}\) satisfying

\[
\mathrm{d}e^{1} = \mathrm{d}e^{2} = \mathrm{d}e^{3} = \mathrm{d}e^{5} = 0,\quad \mathrm{d}e^{4} = e^{1}\wedge e^{5},\quad \mathrm{d}e^{6} = e^{1}\wedge e^{3}.
\]

We say \((G,\omega)\) is a symplectic Lie group or \((\mathfrak{g},\omega)\) is a symplectic Lie algebra if \(\omega \in \bigwedge^{2}\mathfrak{g}^{*}\) is a d-closed and non-degenerate 2-form on \(\mathfrak{g}\). For a symplectic Lie group or a symplectic Lie algebra, its dimension must an even number \(2m\). For any \(0\leq k\leq m\), we can define the space of primitive \(k\)-forms as

\[
\mathcal{P}^{k}\mathfrak{g}^{*} = \left\{\alpha \in \bigwedge^{k}\mathfrak{g}^{*}:\omega^{m - k + 1}\wedge \alpha = 0\right\} .
\]

We have the well-known Lefschetz decomposition

\[
\bigwedge^{k}\mathfrak{g}^{*} = \mathcal{P}^{k}\mathfrak{g}^{*}\oplus \left(\omega \wedge \mathcal{P}^{k - 2}\mathfrak{g}^{*}\right)\oplus \left(\omega^{2}\wedge \mathcal{P}^{k - 4}\mathfrak{g}^{*}\right)\oplus \ldots
\]

and the Lefschetz isomorphism

\[
\omega^{m - k}:\bigwedge^{k}\mathfrak{g}^{*}\cong \bigwedge^{2m - k}\mathfrak{g}^{*}.
\]

One can easily show that

\[
\mathrm{d}:\mathcal{P}^{k}\mathfrak{g}^{*}\to \mathcal{P}^{k + 1}\mathfrak{g}^{*}\oplus \left(\omega \wedge \mathcal{P}^{k - 1}\mathfrak{g}^{*}\right),
\]

therefore one decompose \(\mathrm{d} = \partial_{+} + \omega \wedge \partial_{-}\), where \(\partial_{\pm}:\mathcal{P}^{k}\mathfrak{g}^{*}\to \mathcal{P}^{k\pm 1}\mathfrak{g}^{*}\) are first order differential operators mapping primitive forms to primitive forms. The symplectic cohomology \(\mathrm{SH}^{*}\) of \((\mathfrak{g},\omega)\) is defined to be the cohomology of the cochain complex

\[
0\rightarrow \mathcal{P}^{0}\mathfrak{g}^{*}\xrightarrow{\partial_{+}}\mathcal{P}^{1}\mathfrak{g}^{*}\xrightarrow{\partial_{+}}\ldots\xrightarrow{\partial_{+}}\mathcal{P}^{m}\mathfrak{g}^{*}\xrightarrow{\partial_{+}\partial_{-}}\mathcal{P}^{m}\mathfrak{g}^{*}\xrightarrow{\partial_{-}}\mathcal{P}^{m-1}\mathfrak{g}^{*}\xrightarrow{\partial_{-}}\mathcal{P}^{m-1}\mathfrak{s}^{*}\xrightarrow{\partial_{-}}\mathcal{P}^{0}\mathfrak{g}^{*}\to 0.
\]

In this paper, we are concerned with three kinds of cohomologies on \(\mathfrak{g}\). They are the de Rham/Lie algebra cohomology group

\[
\mathrm{H}^{m}(\mathfrak{g};\mathbb{R})\coloneqq \frac{\ker\left(\mathrm{d}:\bigwedge^{m}\mathfrak{g}^{*}\to\bigwedge^{m + 1}\mathfrak{g}^{*}\right)}{\mathrm{im}\left(\mathrm{d}:\bigwedge^{m - 1}\mathfrak{g}^{*}\to\bigwedge^{m}\mathfrak{g}^{*}\right)},
\]

the primitive part of the de Rham/Lie algebra cohomology

\[
\mathrm{PH}^{m}(\mathfrak{g};\mathbb{R})\coloneqq \frac{\mathcal{P}^{m}\mathfrak{g}^{*}\cap\ker\left(\mathrm{d}:\bigwedge^{m}\mathfrak{g}^{*}\to\bigwedge^{m + 1}\mathfrak{g}^{*}\right)}{\mathcal{P}^{m}\mathfrak{g}^{*}\cap\lim \left(\mathrm{d}:\bigwedge^{m - 1}\mathfrak{g}^{*}\to\bigwedge^{m}\mathfrak{g}^{*}\right)},
\]

and the symplectic cohomologies

\[
\mathrm{SH}_{+}^{m}(\mathfrak{g};\mathbb{R})\coloneqq \frac{\ker(\partial_{+}\partial_{-}:\mathcal{P}^{m}\mathfrak{g}^{*}\to\mathcal{P}^{m}\mathfrak{g}^{*})}{\mathrm{im}(\partial_{+}:\mathcal{P}^{m - 1}\mathfrak{g}^{*}\to\mathcal{P}^{m}\mathfrak{g}^{*})}
\]

and

\[
\mathrm{SH}_{-}^{m}(\mathfrak{g};\mathbb{R})\coloneqq \frac{\ker(\partial_{-}:\mathcal{P}^{m}\mathfrak{g}^{*}\to\mathcal{P}^{m - 1}\mathfrak{g}^{*})}{\mathrm{im}(\partial_{+}\partial_{-}:\mathcal{P}^{m}\mathfrak{g}^{*}\to\mathcal{P}^{m}\mathfrak{g}^{*})}.
\]

The Poincare duality says the natural pairing \(\wedge :\mathrm{SH}_{+}^{m}(\mathfrak{g};\mathbb{R})\otimes \mathrm{SH}_{- }^{m}(\mathfrak{g};\mathbb{R})\to \bigwedge^{2m}\mathfrak{g}^{*}\) is nondegenerate if \(\mathfrak{g}\) is unimodular.

Since \(\partial_{+}\partial_{- } = \mathrm{d}\Lambda_{\omega}\mathrm{d}\) and \(\ker \partial_{- } = \ker \mathrm{d}\cap \mathcal{P}^{m}\mathfrak{g}^{*}\), we naturally have the surjective map

\[
\mathrm{SH}_{-}^{m}(\mathfrak{g};\mathbb{R})\to \mathrm{PH}^{m}(\mathfrak{g};\mathbb{R})
\]

and the injective map

\[
\mathrm{PH}^{m}(\mathfrak{g};\mathbb{R})\hookrightarrow \mathrm{H}^{m}(\mathfrak{g};\mathbb{R}).
\]

On the other hand, as \(\mathrm{im}\mathrm{d}\cap \mathcal{P}^{m}\mathfrak{g}^{*}\subset \mathrm{im}\partial_{+}\), we also have a natural map

\[
\mathrm{SH}_{+}^{m}(\mathfrak{g};\mathbb{R})\to \mathrm{PH}^{m}(\mathfrak{g};\mathbb{R}).
/]

## Some Explicit Calculations

Let $V$ be a 6 dimensional real vector space equipped with a symplectic form $/omega$. For any 3-form \(\phi \in \bigwedge^{3}V^{*}\) , we can define \(\mathcal{K}(\phi)\in \operatorname {End}V\otimes \bigwedge^{6}V^{*}\) by

\[
\mathcal{K}(\phi)(v) = -\iota_{v}\phi \wedge \phi \in \bigwedge^{5}V^{*}\cong V\otimes \bigwedge^{6}V^{*}
\]

for any \(v\in V\) . The isomorphism \(\bigwedge^{5}V^{*}\cong V\otimes \bigwedge^{6}V^{*}\) is canonical. If we choose a basis \(\{e_{1},\ldots ,e_{6}\}\) of \(V\) with dual basis \(\{e^{1},\ldots ,e^{6}\}\) of \(V^{*}\) , we have

\[
\mathcal{K}(\phi)(v) = -\sum_{i = 1}^{6}e_{i}\otimes e^{i}\wedge \iota_{v}\phi \wedge \phi \in V\otimes \bigwedge^{6}V^{*}.
\]

In addition, we can define \(\mathcal{F}(\phi)\in \bigwedge^{3}V^{*}\otimes \bigwedge^{6}V^{*}\) by

\[
\mathcal{F}(\phi)(v_{1},v_{2},v_{3}) = -2\phi (\mathcal{K}(\phi)(v_{1}),v_{2},v_{3})
\]

and \(\mathcal{Q}(\phi)\in (\bigwedge^{6}V^{*})^{\otimes 2}\) by

\[
\mathcal{Q}(\phi) = -\phi \wedge \mathcal{F}(\phi).
\]

It is clear from their definitions that \(\mathcal{K}\) \(\mathcal{F}\) , and \(\mathcal{Q}\) are quadratic, cubic, and quartic in \(\phi\) respectively. Moreover, \(\mathcal{K}\) \(\mathcal{F}\) , and \(\mathcal{Q}\) satisfies the following identities:

\[
\begin{array}{r l} 
& {\mathcal{K}(\phi)\circ \mathcal{K}(\phi) = \frac{\mathrm{id}_{V}}{4}\cdot \mathcal{Q}(\phi)\in \mathrm{End}V\otimes (\bigwedge^{6}V^{*})^{\otimes 2},}\\
& {\mathcal{K}(\mathcal{F}(\phi)) = -\mathcal{K}(\phi)\cdot \mathcal{Q}(\phi)\in \mathrm{End}V\otimes (\bigwedge^{6}V^{*})^{\otimes 3},}\\
& {\mathcal{F}(\mathcal{F}(\phi)) = -\phi \cdot \mathcal{Q}^{2}(\phi)\in \bigwedge^{3}V^{*}\otimes (\bigwedge^{6}V^{*})^{\otimes 4}.}
\end{array} \quad (c).
\]

What we presented so far is independent of the symplectic structure. The symplectic form \(\omega\) determines a canonical volume form \(\omega^{3} / 3! \in \bigwedge^{6}V^{*}\) , with which we can view \(\mathcal{K}(\phi)\) , \(\mathcal{F}(\phi)\) and \(\mathcal{Q}(\phi)\) as elements in \(\operatorname {End}V\) , \(\bigwedge^{3}V^{*}\) and \(\mathbb{R}\) respectively.

A notable example is that for \(\phi\) such that \(\mathcal{Q}(\phi)< 0\) , we can define a complex structure \(\mathcal{J}(\phi)\) on \(V\) by

\[
\mathcal{J}(\phi) = \frac{2\mathcal{K}(\phi)}{\sqrt{-\mathcal{Q}(\phi)}}.
\]

The 3- form \(\phi\) is the real part of a complex \((3,0)\) - form \(\phi +i\hat{\phi}\) with respect to \(\mathcal{J}(\phi)\) . The imaginary part \(\hat{\phi}\) is purely determined from \(\phi\) by

\[
\hat{\phi} = \mathcal{J}(\phi)^{*}\phi = \frac{\mathcal{F}(\phi)}{\sqrt{-\mathcal{Q}(\phi)}}.
\]

If we define \(|\phi |^{2}\) by

\[
\phi \wedge \mathcal{F}(\phi) = |\phi |^{4}\frac{\omega^{3}}{3!},
\]

then we have \(\mathcal{K}(\phi) = \frac{|\phi|^{2}}{2}\mathcal{J}(\phi)\) , \(\mathcal{F}(\phi) = |\phi |^{2}\hat{\phi}\) , and \(\mathcal{Q}(\phi) = -|\phi |^{4}\) .

Under the presence of the symplectic form \(\omega\) , we usually only consider \(\phi\) that is primitive with respect to \(\omega\) . In this case, \(\mathcal{F}(\phi)\) is also primitive. If in addition that \(\mathcal{Q}(\phi)< 0\) , we know that \(\omega\) is a \((1,1)\) - form with respect to the complex structure \(\mathcal{J}(\phi)\) and that \(|\phi |^{2}\) is indeed the norm square of \(\phi\) with respect to the metric (not necessarily positive definite) \(g(\cdot ,\cdot) = \omega (\cdot ,\mathcal{J}(\phi)\cdot)\) .

Given a 6- dimensional symplectic vector space \((V,\omega)\) , we can choose a basis \(\{e_{1},\ldots ,e_{6}\}\) of \(V\) and its dual basis \(\{e^{1},\ldots ,e^{6}\}\) of \(V^{*}\) such that the symplectic form \(\omega\) takes the standard form \(\omega = e^{12} + e^{34} + e^{56}\) . Under such a choice, a general 3- form \(\phi \in \bigwedge^{3}V^{*}\) takes the form

\[
\begin{array}{rcl}
{\phi} & = & {A e^{135} + B e^{136} + C e^{145} + D e^{146} + E e^{235} + F e^{236} + G e^{245} + H e^{246}}\\
& {} & {+(I e^{1} + J e^{2})(e^{34} - e^{56}) + (K e^{3} + L e^{4})(e^{12} - e^{56}) + (M e^{5} + N e^{6})(e^{12} - e^{34})}\\
& {} & {+(O e^{1} + P e^{2})(e^{34} + e^{56}) + (Q e^{3} + R e^{4})(e^{12} + e^{56}) + (S e^{5} + T e^{6})(e^{12} + e^{34}),}
\end{array} \quad (2)
\]

where \(A,B,\ldots ,S,T\) are constants. And \(\phi\) is primitive if and only if the last 6 coefficients \(O,P,Q,R,S,T\) are all zero, namely \(\phi\) takes the form

\[
\begin{array}{rcl}
{\phi} & = & {A e^{135} + B e^{136} + C e^{145} + D e^{146} + E e^{235} + F e^{236} + G e^{245} + H e^{246}}\\
& {} & {+(I e^{1} + J e^{2})(e^{34} - e^{56}) + (K e^{3} + L e^{4})(e^{12} - e^{56}) + (M e^{5} + N e^{6})(e^{12} - e^{34}).}
\end{array} \quad (3)
\]

For later use, we shall compute the expression of \(\mathcal{K}(\phi)\) , \(\mathcal{F}(\phi)\) , and \(\mathcal{Q}(\phi)\) for \(\phi\) in (2). The explicit values for $F$ for this choice of base can be found in the file Polynomials.cpp, for example the letter A is sent to A_hat. 

## Examples

We want to study the problem of existence and uniqueness of \(\mathcal{F}\) - harmonic forms, and the long- time behavior of the Type IIA flow on symplectic Lie groups with left invariant data. To be more specific, we would like to address the following two questions.

**Question 1: The existence and uniqueness of \(\mathcal{F}\) - harmonic 3- forms in a fixed cohomology class.**

By definition, a 3- form \(\phi\) is \(\mathcal{F}\) - harmonic if \(\mathrm{d}\phi = 0\) and \(\mathrm{d}\mathcal{F}(\phi) = 0\) . As an analogue of the standard Hodge theory, we would like to address Question 1 in the set- up of a symplectic Lie algebra \((\mathfrak{g},\omega)\) with \(\phi \in \Lambda^3\mathfrak{g}^*\) . Since \(\mathrm{d}\phi = 0\) , we know that \(\phi\) defines a de Rham cohomology class in \(\mathrm{H}^3 (\mathfrak{g};\mathbb{R})\) for a general 3- form, and that \(\phi\) defines a cohomology class in \(\mathrm{PH}^3 (\mathfrak{g};\mathbb{R})\) and \(\mathrm{SH}_{\pm}^{3}(\mathfrak{g};\mathbb{R})\) if \(\phi\) is primitive. We would like to know that, if the cohomology class of \(\phi\) is fixed in the above sense (de Rham for general \(\phi\) , and one of the three kinds of cohomologies for primitive \(\phi\) ), does there exist a unique \(\mathcal{F}\) - harmonic representative.

For a given symplectic Lie algebra \((\mathfrak{g},\omega)\) , to answer Question 1, we proceed along the following steps.

Step 1: Compute \(\mathrm{H}^3 (\mathfrak{g};\mathbb{R})\) , \(\mathrm{PH}^3 (\mathfrak{g};\mathbb{R})\) , and \(\mathrm{SH}_{\pm}^{3}(\mathfrak{g};\mathbb{R})\) explicitly.

Step 2: By assuming \(\phi\) taking the form of (2) or (3), we reduce the \(F\) - harmonic equation to a system of polynomial equations.

Step 3: Determine the existence and uniqueness of solutions to the system of algebraic equations in Step 2 within each cohomology class.

**Question 2: The long- time behavior of the Type IIA flow.**

We would like to investigate the long-time behavior of the Type IIA flow (1) on symplectic Lie groups and their quotients, including whether the flow has long- time existence of finite- time singularities. In its original setup, the initial data is both closed, primitive, and positive. However, our formulation in (1) allows us to extend the initial data to be a general 3- form, not necessarily closed or primitive. In this paper, we will consider such a flow with initial data under three kinds of generalities: (1) a general 3- form; (2) a primitive 3- form; (3) a closed primitive 3- form.

It is worth mentioning that there is little difference between the case of general initial data and that of primitive initial data. The Lefschetz decomposition for 3- form reads

\[
\bigwedge^{3}\mathfrak{g}^{*} = \mathcal{P}^{3}\mathfrak{g}^{*}\oplus \omega \wedge \mathfrak{g}^{*}.
\]

For any 1- form \(\alpha \in \mathfrak{g}^{*}\) , we have that

\[
\mathrm{d}\Lambda_{\omega}\mathrm{d}(\omega \wedge \alpha) = \mathrm{d}\Lambda_{\omega}(\omega \wedge \mathrm{d}\alpha) = \mathrm{d}(\mathrm{d}\alpha +\omega \wedge \Lambda_{\omega}\mathrm{d}\alpha) = 0
\]

since \(\Lambda_{\omega}\mathrm{d}\alpha\) is a constant. This computation implies that \(\mathrm{d}\Lambda_{\omega}\mathrm{d}\mathcal{F}(\phi)\) has only primitive components, hence the non- primitive components of \(\phi\) are stationary under the Type IIA flow (1) with invariant data.

To address Question 2, we first reduce the Type IIA flow to an ODE system. Then, we find its stationary points. Finally, we analyze the long- time behavior of this ODE system.

## What this project was able to do

This project was able to compute some calculations to answer of the previous two questions with the classification of 6-dimensional symplectic algebras found in the paper Symplectic or contact structures on Lie groups. [docs/symplectic-lie-groups.pdf]. There one can find 23 different symplectic lie algebras. In the file [res/input.txt] one can find the description of these lie algebras using the canonical \(\omega = e^{12} + e^{34} + e^{56}\).
