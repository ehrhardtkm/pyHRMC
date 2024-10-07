Theory 
===

## Monte Carlo and loss function

The H/RMC algorithm employs the Metropolis Markov chain Monte Carlo method to stochastically sample atomistic configurations of a structure. The simulations converge by minimizing a loss function, which follows the general form:

$$
\chi^{2}_{total} = \sum_{i}^{n} \sum_{j} \frac{(A_{i}(j)_{target} - A_{i}(j)_{calc})^{2}} {W_{i}} .
$$




Here, the total loss $\chi^{2}_{total}$ is summed over the *n* constraints that guide the simulation. Each term in the loss function represents a squared difference between the calculated value *$A_{calc}$* and the target value *$A_{target}$*, scaled by a weight term *$W_{i}$*. The number of terms in this function depends on the constraints that are defined by the user. Typically, in RMC simulations, this loss function includes an error term that fits the calculated pair distribution function (PDF) to experimental PDF data 

$$
\chi^{2}_{error} = \sum_{i}^{n_{bin}}\frac{(G(r_{i})_{target} - G(r_{i})_{calc})^{2}} {\sigma_{}} .
$$

where the squared difference between the experimental/target PDF $G(r)_{target}$ and the calculated PDF $G(r)_{calc}$ of the Monte Carlo structure is scaled with *$\sigma$*, the experimental error that acts as the weight term, and is summed over $n_{bin}$, the number of bins of $G(r)_{target}$. 

In pyHRMC, the PDF is always assumed to be the electron pair distribution function (ePDF). This distinction is critical because the scattering cross-sections used to calculate the theoretical electron scattering—and, consequently, the PDF—depend on the ePDF of the simulated structure.

While RMC is effective at matching a structure's calculated ePDF to experimental data, it cannot produce a unique solution. In other words, RMC may arrive at one of many different atomistic configurations can produce approximately the same ePDF. . This ambiguity can be somewhat mitigated by choosing simulation cells large enough to yield a statistically robust distribution of coordination environments. However, under-constrained RMC simulations may still produce unrealistic configurations.

## Hybrid Reverse Monte Carlo (HRMC)
To address this limitation, HRMC was introduced by Opletal et al. in 2002 (see ref 1). In HRMC, an additional energy term is included in the total loss function, leading to the following form:

$$
\chi^{2}_{total} = \chi^{2}_{error} + \frac { \Delta E} {k_{B}T}
$$

Similar to the error term, the energy term is scaled by a temperature variable (T). This addition results in a Boltzmann-type probability distribution, where the algorithm may accept a step even when the energy loss term increases.

**BUG**
The program will always accept the proposed Monte Carlo step if the new loss value $\chi^{2}_{new}$ is less than old loss value $\chi^{2}_{old}$. Otherwise, the probability of accepting a move is:

$$
P = e^{\chi^{2}_{old} - \chi^{2}_{new}}
$$

## Quench schemes
As noted above, the weight term in each loss function affects the Boltzmann-type probability distribution of accepting a "worse" step. In traditional RMC, where the PDF error term is the sole contributor to the loss function $\chi^{2}_{total}$, this weight $\sigma$ is typically kept static throughout the simulation. However, initializing the simulation with greater weight terms and gradually reducing them can help optimize convergence. This approach allows the structure to explore a broader structural space early on while still refining toward a low loss value as the simulation progresses.

This following quenching scheme was introduced by Opletal et al. and has been implemented in pyHRMC:

$$
T_{n+1} = T_{n} * q_{T}^{n} 
$$

$$
\sigma_{n+1} = \sigma_{n} * q_{\sigma}^{n} ,
$$

where n is the number of steps in millions, and q is a quenching parameter between 0-1. Because H/RMC is performed over millions of steps and gradual quenching is necessary for proper convergence, q is typically very close to 1.

When selecting values for these quenching parameters in pyHRMC, note the distinction between a [step and an iteration](https://ehrhardtkm.github.io/pyHRMC/runfile_setup/). Quenching will occur after each *iteration*, so greater parallelization will result in a slightly coarser quenching scheme.

***

For further details and examples in the literature, we suggest the following references:

1. [Opletal, G., Petersen, T., O’Malley, B., Snook, I., McCulloch, D. G., Marks, N. A., & Yarovsky, I. Hybrid approach for generating realistic amorphous carbon structure using metropolis and reverse Monte Carlo. Molecular Simulation, 2002, 28(10–11), 927–938. https://doi.org/10.1080/089270204000002584](https://www.tandfonline.com/doi/abs/10.1080/089270204000002584?casa_token=Il9_dIbtCoIAAAAA:FSzPAgSfbQ2FLkG7_C5VcxNbZhqkeCQJLu-GlYR4s9L090xkQ5KbigcF8idYFZ56FoHBulHD11TD)
2. [G. Opletal, T.C. Petersen, B. O'Malley, I.K. Snook, D.G. McCulloch, I. Yarovsky. HRMC: Hybrid Reverse Monte Carlo method with silicon and carbon potentials. Comput. Phys. Commun., 2008, 178 (10), 777-787. https://doi.org/10.1016/j.cpc.2007.12.007.](https://www.sciencedirect.com/science/article/pii/S0010465508000039)
3. [Amir H. Farmahini, Suresh K. Bhatia. Hybrid Reverse Monte Carlo simulation of amorphous carbon: Distinguishing between competing structures obtained using different modeling protocols. Carbon, 2015, 83, 53-70. https://doi.org/10.1016/j.carbon.2014.11.013.](https://www.sciencedirect.com/science/article/pii/S0008622314010859)
 
