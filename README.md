# Objective

The diversity of mutational spectra among organisms has only recently been linked to long-term evolutionary change. In particular, recent work has shown that a reversal from a biased mutational spectrum enhances access to beneficial mutations. We set out to demonstrate that the mutational spectrum will naturally reverse over evolutionary time. To do this, we construct an adaptive walk simulation that evolves the transition bias concurrently with the population itself. Our results demonstrate that a drive for reversals in the mutational spectrum exists, but does not persist indefinetly. Instead the evolutionary benefit for mutational reversals appears to lessen as a population increases its fitness.

# Summary of Simulation

The sequence we evolved was randomly initialized and had length $N = 100$. This simulation assumed a genetically uniform population where selection is stronger than mutation. Since our population is assumed to be uniform, the sequence we evolve is a representation of the population, and we use the terms `organism' and `population' interchangeably. Our simulation has one degree of epistasis ($K = 1$), where each locus is epistatically coupled to $1/N \times 100 \%$ of the loci (since we allow for the possibility of a locus being its own epistatic partner). 
