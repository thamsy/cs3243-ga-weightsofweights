import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.jgap.Chromosome;
import org.jgap.Configuration;
import org.jgap.DefaultFitnessEvaluator;
import org.jgap.FitnessFunction;
import org.jgap.Gene;
import org.jgap.Genotype;
import org.jgap.IChromosome;
import org.jgap.InvalidConfigurationException;
import org.jgap.UnsupportedRepresentationException;
import org.jgap.audit.EvolutionMonitor;
import org.jgap.data.DataTreeBuilder;
import org.jgap.data.IDataCreators;
import org.jgap.event.EventManager;
import org.jgap.impl.DoubleGene;
import org.jgap.impl.GaussianRandomGenerator;
import org.jgap.impl.MutationOperator;
import org.jgap.impl.StandardPostSelector;
import org.jgap.xml.XMLDocumentBuilder;
import org.jgap.xml.XMLManager;
import org.w3c.dom.Document;

public class HeuristicWeightsLearning {

	// The total number of times we'll let the population evolve.
	private static final int MAX_ALLOWED_EVOLUTIONS = 45000;
	private static final int POPULATION_SIZE = 1000;
	public static EvolutionMonitor m_monitor;

	private static void doLearning(int popSize, int numEvo, boolean a_doMonitor) throws Exception {

		// ----------- Configuration ----------

		Configuration conf = new Configuration();
		// conf.setPreservFittestIndividual(true);
		conf.setKeepPopulationSizeConstant(false);
		conf.setPopulationSize(popSize);

		conf.setRandomGenerator(new GaussianRandomGenerator());
		conf.setEventManager(new EventManager());
		conf.setFitnessEvaluator(new DefaultFitnessEvaluator());
		conf.addNaturalSelector(new StandardPostSelector(conf), true);
		conf.addGeneticOperator(new MutationOperator(conf, 10));
		conf.addGeneticOperator(new UniformCrossoverOperator(conf, 1));

		// Set Fitness Function
		FitnessFunction fitFunc = new HeuristicWeightsFitnessFuction();
		conf.setFitnessFunction(fitFunc);

		// Turn on monitoring/auditing of evolution progress.
		if (a_doMonitor) {
			m_monitor = new EvolutionMonitor();
			conf.setMonitor(m_monitor);
		}

		// Set Genes
		Gene[] sampleGenes = new Gene[9];
		sampleGenes[0] = new DoubleGene(conf, -1d, 1d); // MinHeight for H1
		sampleGenes[1] = new DoubleGene(conf, -1d, 1d); // MaxHeight for H1
		sampleGenes[2] = new DoubleGene(conf, -1d, 1d); // Constant for H1
		sampleGenes[3] = new DoubleGene(conf, -1d, 1d); // MinHeight for H2
		sampleGenes[4] = new DoubleGene(conf, -1d, 1d); // MaxHeight for H2
		sampleGenes[5] = new DoubleGene(conf, -1d, 1d); // Constant for H2
		sampleGenes[6] = new DoubleGene(conf, -1d, 1d); // MinHeight for H3
		sampleGenes[7] = new DoubleGene(conf, -1d, 1d); // MaxHeight for H3
		sampleGenes[8] = new DoubleGene(conf, -1d, 1d); // Constant for H3
		IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
		conf.setSampleChromosome(sampleChromosome);
		Genotype population;

		try {
			Document doc = XMLManager.readFile(new File("GA_mutate10_co1.xml"));
			population = XMLManager.getGenotypeFromDocument(conf, doc);
		} catch (UnsupportedRepresentationException uex) {
			// JGAP codebase might have changed between two consecutive runs.
			// --------------------------------------------------------------
			population = Genotype.randomInitialGenotype(conf);
		} catch (FileNotFoundException fex) {
			population = Genotype.randomInitialGenotype(conf);
		}

		// ------ Start -------

		long startTime = System.currentTimeMillis();
		for (int i = 0; i < numEvo; i++) {
			if (m_monitor != null) {
				population.evolve(m_monitor);
			} else {
				population.evolve();
			}
			System.out.printf("Evolution no. %d, Current Fittest Fitness Value is: %f%n", i + 1,
					population.getFittestChromosome().getFitnessValue());

			// Build Document
			DataTreeBuilder builder = DataTreeBuilder.getInstance();
			IDataCreators doc2 = builder.representGenotypeAsDocument(population);
			// create XML document from generated tree
			XMLDocumentBuilder docbuilder = new XMLDocumentBuilder();
			Document xmlDoc = (Document) docbuilder.buildDocument(doc2);
			XMLManager.writeFile(xmlDoc, new File("GA_mutate10_co1.xml"));

		}
		long endTime = System.currentTimeMillis();
		System.out.println("Total evolution time: " + (endTime - startTime) + " ms");

		// Print results
		IChromosome bestSolutionSoFar = population.getFittestChromosome();
		System.out.println("The best solution has a fitness value of " + bestSolutionSoFar.getFitnessValue());
		System.out.println("It contains the following: ");
		System.out.printf("\tW1 is %f%n", (double) bestSolutionSoFar.getGene(0).getAllele());
		System.out.printf("\tW2 is %f%n", (double) bestSolutionSoFar.getGene(1).getAllele());
		System.out.printf("\tW3 is %f%n", (double) bestSolutionSoFar.getGene(2).getAllele());
		System.out.printf("\tW4 is %f%n", (double) bestSolutionSoFar.getGene(3).getAllele());
		System.out.printf("\tW5 is %f%n", (double) bestSolutionSoFar.getGene(4).getAllele());
		System.out.printf("\tW6 is %f%n", (double) bestSolutionSoFar.getGene(5).getAllele());
		System.out.printf("\tW7 is %f%n", (double) bestSolutionSoFar.getGene(6).getAllele());
		System.out.printf("\tW8 is %f%n", (double) bestSolutionSoFar.getGene(7).getAllele());
		System.out.printf("\tW9 is %f%n", (double) bestSolutionSoFar.getGene(8).getAllele());

		// Save in file
		try {
			PrintWriter writer = new PrintWriter("results.txt", "UTF-8");
			writer.println("The best solution has a fitness value of " + bestSolutionSoFar.getFitnessValue());
			writer.println("It contains the following: ");
			writer.printf("\tW1 is %f%n", (double) bestSolutionSoFar.getGene(0).getAllele());
			writer.printf("\tW2 is %f%n", (double) bestSolutionSoFar.getGene(1).getAllele());
			writer.printf("\tW3 is %f%n", (double) bestSolutionSoFar.getGene(2).getAllele());
			writer.printf("\tW4 is %f%n", (double) bestSolutionSoFar.getGene(3).getAllele());
			writer.printf("\tW5 is %f%n", (double) bestSolutionSoFar.getGene(4).getAllele());
			writer.printf("\tW6 is %f%n", (double) bestSolutionSoFar.getGene(5).getAllele());
			writer.printf("\tW7 is %f%n", (double) bestSolutionSoFar.getGene(6).getAllele());
			writer.printf("\tW8 is %f%n", (double) bestSolutionSoFar.getGene(7).getAllele());
			writer.printf("\tW9 is %f%n", (double) bestSolutionSoFar.getGene(8).getAllele());
			writer.close();
		} catch (IOException e) {
			// do something
		}
	}

	public static void main(String args[]) {
		try {
			int popSize = Integer.parseInt(args[0]);
			int numEvo = Integer.parseInt(args[1]);
			doLearning(popSize, numEvo, false);
		} catch (NumberFormatException e) {
			System.out.println("Arguments must be a valid integer value");
			System.exit(1);
		} catch (InvalidConfigurationException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
