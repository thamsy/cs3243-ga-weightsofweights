import java.util.List;
import java.util.Vector;

import org.jgap.Configuration;
import org.jgap.Gene;
import org.jgap.IChromosome;
import org.jgap.ICompositeGene;
import org.jgap.IGeneticOperatorConstraint;
import org.jgap.InvalidConfigurationException;
import org.jgap.Population;
import org.jgap.RandomGenerator;
import org.jgap.impl.CrossoverOperator;

public class UniformCrossoverOperator extends CrossoverOperator {

	public UniformCrossoverOperator(Configuration a_configuration, int a_desiredCrossoverRate)
			throws InvalidConfigurationException {
		super(a_configuration, a_desiredCrossoverRate);
	}

	@Override
	public void operate(final Population a_population, final List a_candidateChromosomes) {
		// Work out the number of crossovers that should be performed.
		// -----------------------------------------------------------
		int size = Math.min(getConfiguration().getPopulationSize(), a_population.size());
		int numCrossovers = 0;
		if (getCrossOverRate() >= 0) {
			numCrossovers = size / getCrossOverRate() * a_population.getChromosome(0).getGenes().length;
		} else {
			numCrossovers = (int) (size * getCrossOverRatePercent());
		}
		RandomGenerator generator = getConfiguration().getRandomGenerator();
		IGeneticOperatorConstraint constraint = getConfiguration().getJGAPFactory().getGeneticOperatorConstraint();
		// For each crossover, grab two random chromosomes, pick a random
		// locus (gene location), and then swap that gene and all genes
		// to the "right" (those with greater loci) of that gene between
		// the two chromosomes.
		// --------------------------------------------------------------
		int index1, index2;
		//
		for (int i = 0; i < numCrossovers; i++) {
			index1 = generator.nextInt(size);
			index2 = generator.nextInt(size);
			IChromosome chrom1 = a_population.getChromosome(index1);
			IChromosome chrom2 = a_population.getChromosome(index2);
			// Verify that crossover is allowed.
			// ---------------------------------
			if (!isXoverNewAge() && chrom1.getAge() < 1 && chrom2.getAge() < 1) {
				// Crossing over two newly created chromosomes is not seen as
				// helpful
				// here.
				// ------------------------------------------------------------------
				continue;
			}
			if (constraint != null) {
				List v = new Vector();
				v.add(chrom1);
				v.add(chrom2);
				if (!constraint.isValid(a_population, v, this)) {
					// Constraint forbids crossing over.
					// ---------------------------------
					continue;
				}
			}
			// Clone the chromosomes.
			// ----------------------
			IChromosome firstMate = (IChromosome) chrom1.clone();
			IChromosome secondMate = (IChromosome) chrom2.clone();
			// In case monitoring is active, support it.
			// -----------------------------------------
			if (m_monitorActive) {
				firstMate.setUniqueIDTemplate(chrom1.getUniqueID(), 1);
				firstMate.setUniqueIDTemplate(chrom2.getUniqueID(), 2);
				secondMate.setUniqueIDTemplate(chrom1.getUniqueID(), 1);
				secondMate.setUniqueIDTemplate(chrom2.getUniqueID(), 2);
			}
			// Cross over the chromosomes.
			// ---------------------------
			doCrossover(firstMate, secondMate, a_candidateChromosomes, generator);
		}
	}

	protected void doCrossover(IChromosome firstMate, IChromosome secondMate, List a_candidateChromosomes,
			RandomGenerator generator) {
		Gene[] firstGenes = firstMate.getGenes();
		Gene[] secondGenes = secondMate.getGenes();
		int locus = generator.nextInt(firstGenes.length);
		// Swap the genes.
		// ---------------
		Gene gene1;
		Gene gene2;
		Object firstAllele;

		// Make a distinction for ICompositeGene for the first gene.
		// ---------------------------------------------------------
		int index = 0;
		if (firstGenes[locus] instanceof ICompositeGene) {
			// Randomly determine gene to be considered.
			// -----------------------------------------
			index = generator.nextInt(firstGenes[locus].size());
			gene1 = ((ICompositeGene) firstGenes[locus]).geneAt(index);
		} else {
			gene1 = firstGenes[locus];
		}
		// Make a distinction for the second gene if CompositeGene.
		// --------------------------------------------------------
		if (secondGenes[locus] instanceof ICompositeGene) {
			gene2 = ((ICompositeGene) secondGenes[locus]).geneAt(index);
		} else {
			gene2 = secondGenes[locus];
		}
		if (m_monitorActive) {
			gene1.setUniqueIDTemplate(gene2.getUniqueID(), 1);
			gene2.setUniqueIDTemplate(gene1.getUniqueID(), 1);
		}
		firstAllele = gene1.getAllele();
		gene1.setAllele(gene2.getAllele());
		gene2.setAllele(firstAllele);

		// Add the modified chromosomes to the candidate pool so that
		// they'll be considered for natural selection during the next
		// phase of evolution.
		// -----------------------------------------------------------
		a_candidateChromosomes.add(firstMate);
		a_candidateChromosomes.add(secondMate);
	}

}
