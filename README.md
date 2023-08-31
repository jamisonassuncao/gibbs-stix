# Gibbs free energy

The Gibbs free energy $\mathcal{G}$ of a multiphase assemblage is

$$
\mathcal{G}(P,T,\vec{n})=\sum_{i}^{species}n_i\mu_i(P,T,\vec{n})
$$

where $P$, $T$, $n_i$, $\mu_i$ are, respectively, pressure, temperature, amount and chemical potential, and $\vec{n}$ is the vector containing all of the $n_i$.

The chemical potential is assumed to be

$$
\mu_i=
\mathcal{G}_i-RT\sum_{k}^{sites}
\Big[S_{ik}\ln{N_k}-\sum_{j}^{c}s_{ijk}\ln{N_{jk}}\Big] \\
-\sum_{\beta> \alpha}^{species}W_{i\alpha\beta}
(\delta_{i\alpha}-\phi_{\alpha})
(\delta_{i\beta}-\phi_{\beta})
\\=
\mathcal{G}_i-RT\sum_{k}^{sites}
\Big[
\Big(\sum_{j}^{c}s_{ijk}\Big)\ln{\Big(\sum_{j}^{c}N_{jk}\Big)}-\\
\sum_{j}^{c}\Big(s_{ijk}\ln{\Big(\sum_{i}^{species}s_{ijk}n_i}\Big)\Big)\Big] \\
-\sum_{\beta>\alpha}^{species}W_{i\alpha\beta}
(\delta_{i\alpha}-\phi_{\alpha})
(\delta_{i\beta}-\phi_{\beta})
$$

where $N_{jk}$, $N_k$ and $S_{ik}$ are, respectively, the number of atoms of component $j$ on site $k$, the total number of atoms on site $k$, and the sum over the stoichiometric coefficients of species $i$ at site $k$, $s_{ijk}$ is the stoichiometric coefficient of component $j$ on site $k$ in species $i$, and $c$ is the number of components, $\mathcal{G}$ is the Gibbs free energy in cation-ordered pure form, $R$ is the gas constant, and $T$ is the temperature.




