import org.jgap.FitnessFunction;
import org.jgap.Gene;
import org.jgap.IChromosome;

public class HeuristicWeightsFitnessFuction extends FitnessFunction {

	public static final int MAX_PIECES = 200000;
	public static final int NUM_GAMES = 10;

	// Returns value of fitness function
	@Override
	public double evaluate(IChromosome a_subject) {

		// Get Weights
		Gene[] genes = a_subject.getGenes();
		double[] weights = new double[9];
		for (int i = 0; i < genes.length; i++) {
			weights[i] = (double) genes[i].getAllele();
		}

		double totalValue = 0.0;

		for (int i = 0; i < NUM_GAMES; i++) {
			// Initialize Game
			State s = new State();
			// new TFrame(s);
			double numPiecesPlayed = 0.0d;
			while (!s.hasLost() && numPiecesPlayed < MAX_PIECES) {
				s.makeMove(findNextMove(s, weights));
				// s.draw();
				// s.drawNext(0, 0);
				numPiecesPlayed++;
				// try {
				// Thread.sleep(10);
				// } catch (InterruptedException e) {
				// e.printStackTrace();
				// }
			}
			totalValue += (double) s.getRowsCleared() + numPiecesPlayed / 1000d;
		}

		return totalValue;
	}

	// Determine which move to make
	private int[] findNextMove(State s, double[] weights) {
		boolean[][] newField = new boolean[State.ROWS][State.COLS];
		int[] newTop = new int[State.COLS];
		int bestRot = 0;
		int bestPos = 0;
		double bestValue = -Double.MAX_VALUE;

		int nextPiece = s.getNextPiece();
		int[][] legalMoves = s.legalMoves();
		for (int i = 0; i < legalMoves.length; i++) {
			int rot = legalMoves[i][State.ORIENT];
			int pos = legalMoves[i][State.SLOT];
			int rowsCleared = performMove(s, newField, newTop, nextPiece, rot, pos);
			int holes = 0;
			int bumpiness = 0;
			int maxHeight = 0;
			int minHeight = State.ROWS;

			for (int c = 0; c < State.COLS; ++c) {
				boolean blocked = false;
				// total height
				maxHeight = Math.max(maxHeight, s.getTop()[c]);
				minHeight = Math.min(minHeight, s.getTop()[c]);
				// sum of difference of consecutive heights
				if (c > 0)
					bumpiness += (newTop[c] - newTop[c - 1]) * (newTop[c] - newTop[c - 1]);
				for (int r = State.ROWS - 1; r >= 0; --r) {
					if (newField[r][c]) {
						blocked = true;
					} else if (!newField[r][c] && blocked) {
						// number of holes
						holes += 1;
					}
				}
			}

			double value = calculateValueOfField(weights, maxHeight, minHeight, rowsCleared, holes, bumpiness);
			if (value > bestValue) {
				bestValue = value;
				bestRot = rot;
				bestPos = pos;
			}

		}

		return new int[] { bestRot, bestPos };

	}

	private int performMove(State s, boolean[][] newField, int[] newTop, int piece, int rot, int pos) {
		// Perform Deep Copy
		for (int i = 0; i < State.ROWS; i++) {
			for (int j = 0; j < State.COLS; j++) {
				newField[i][j] = s.getField()[i][j] != 0;
			}
		}
		for (int j = 0; j < State.COLS; j++) {
			newTop[j] = s.getTop()[j];
		}

		// height if the first column makes contact
		int height = newTop[pos] - State.getpBottom()[piece][rot][0];
		// for each column beyond the first in the piece
		for (int c = 0; c < State.getpWidth()[piece][rot]; c++) {
			height = Math.max(height, newTop[pos + c] - State.getpBottom()[piece][rot][c]);
		}

		// check if game ended
		if (height + State.getpHeight()[piece][rot] >= State.ROWS) {
			// cout << pieceIndex << " " << rotationIndex << " " << leftPosition
			// << endl;
			// cout << height << " " << pHeight[pieceIndex][rotationIndex] << "
			// " << top[leftPosition] << endl;
			// cout << " You lost :( " << height +
			// pHeight[pieceIndex][rotationIndex] << endl;
			return -10;
		}

		// for each column in the piece - fill in the appropriate blocks
		for (int i = 0; i < State.getpWidth()[piece][rot]; i++) {

			// from bottom to top of brick
			for (int h = height + State.getpBottom()[piece][rot][i]; h < height + State.getpTop()[piece][rot][i]; h++) {
				newField[h][i + pos] = true;
			}
		}

		// adjust top
		for (int c = 0; c < State.getpWidth()[piece][rot]; c++) {
			newTop[pos + c] = height + State.getpTop()[piece][rot][c];
		}

		int rowsCleared = 0;

		// check for full rows - starting at the top
		for (int r = height + State.getpHeight()[piece][rot] - 1; r >= height; r--) {
			// check all columns in the row
			boolean full = true;
			for (int c = 0; c < State.COLS; c++) {
				if (newField[r][c] == false) {
					full = false;
					break;
				}
			}
			// if the row was full - remove it and slide above stuff down
			if (full) {
				rowsCleared++;
				// for each column
				for (int c = 0; c < State.COLS; c++) {
					// slide down all bricks
					for (int i = r; i < newTop[c]; i++) {
						newField[i][c] = newField[i + 1][c];
					}
					// lower the top
					newTop[c]--;
					while (newTop[c] >= 1 && newField[newTop[c] - 1][c] == false)
						newTop[c]--;
				}
			}
		}
		return rowsCleared;
	}

	private double calculateValueOfField(double[] weights, int maxHeight, int minHeight, int rowsCleared, int holes,
			int bumpiness) {
		return (weights[0] * maxHeight + weights[1] * (maxHeight - minHeight) + weights[2]) * (double) rowsCleared
				* Math.abs(rowsCleared) / 5
				+ (weights[3] * maxHeight + weights[4] * (maxHeight - minHeight) + weights[5]) * (double) holes / 10
				+ (weights[6] * maxHeight + weights[7] * (maxHeight - minHeight) + weights[8]) * (double) bumpiness
						/ 1000;
	}

	// ----- Main Heuristics ------
	// // H1
	// private int getCompleteLines(State state) {
	//
	// }
	//
	// // H2
	// private int getHoles(State state) {
	//
	// }
	//
	// // H3
	// private int getBumpiness(State state) {
	//
	// }

}
