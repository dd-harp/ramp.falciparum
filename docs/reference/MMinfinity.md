# MMinfinity

The queuing model \\M/M/\infty\\ tracks the MoI in a cohort of humans as
it ages. It assumes a time- and age-dependent hazard rate for infection,
called the force of infection (FoI, \\h\_\tau(a)\\). Infections do not
affect each other, and each one clears independently at the rate \\r\\.

Let \\\zeta_i\\ the fraction of the population with MoI = i, then
\$\$\frac{d\zeta_0}{da}= -h\_\tau(a) \zeta_0 + r \zeta_1\$\$ and for
\\i\geq 1\\ \$\$\frac{d\zeta_i}{da}= h\_\tau(a) \left( \zeta\_{i-1} -
\zeta_i \right) - ri \zeta_i + r(i+1)\zeta\_{i+1}\$\$
