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
\Big[S_{ik}\ln{N_k}-\sum_{j}^{comp.}s_{ijk}\ln{N_{jk}}\Big]

$$

$$
-\sum_{\beta> \alpha}^{species}W_{i\alpha\beta}
(\delta_{i\alpha}-\phi_{\alpha})
(\delta_{i\beta}-\phi_{\beta})
$$

where

$$
N_{jk}=\sum_{i}^{species}s_{ijk}n_i
$$

$$
N_k=\sum_{j}^{comp.}N_{jk}
$$

$$
S_{ik}=\sum_{j}^{comp.}s_{ijk}
$$

are, respectively, the number of atoms of component $j$ on site $k$, the total number of atoms on site $k$, and the sum over the stoichiometric coefficients of species $i$ at site $k$, $s_{ijk}$ is the stoichiometric coefficient of component $j$ on site $k$ in species $i$, and $comp.$ is the number of components, $\mathcal{G}$ is the Gibbs free energy in cation-ordered pure form, $R$ is the gas constant, and $T$ is the temperature.

$$
fo=Mg_2SiO_4=2MgO+1SiO_2\\
fa=Fe_2SiO_4=2FeO+1SiO_2
$$

$$
S_{ik}=S_{fo,k} = s_{fo,SiO_2,k}+s_{fo,MgO,k}
$$

$$
N_{jk}=N_{SiO_2,k}=s_{fo,SiO_2,k}\times n_{fo}+s_{fa,SiO_2,k}\times n_{fa}
$$

$$
N_k=s_{fo,SiO_2,k}\times n_{fo}+s_{fa,SiO_2,k}\times n_{fa}\\
+s_{fo,MgO,k}\times n_{fo}+s_{fa,MgO,k}\times n_{fa}\\
+s_{fo,FeO,k}\times n_{fo}+s_{fa,FeO,k}\times n_{fa}
$$

$$
\Sigma (s_{ijk} \ln{N_jk})=\\s_{fo,SiO_2,k}\times \ln{(s_{fo,SiO_2,k}\times n_{fo}+s_{fa,SiO_2,k}\times n_{fa})}\\
+s_{fo,MgO,k}\times \ln{(s_{fo,MgO,k}\times n_{fo}+s_{fa,MgO,k}\times n_{fa})}\\
+s_{fo,FeO,k}\times \ln{(s_{fo,FeO,k}\times n_{fo}+s_{fa,FeO,k}\times n_{fa})}
$$


