# Tetris AI

Genetic Algorithm, learning weights of game state, used for weights of features.

## Setting up

1. Use Eclipse to import project
2. Add java files from CS3230 java.zip to src
3. Add genetic algorithm jar to project. (JGAP 3.6.3 from https://sourceforge.net/projects/jgap/files/jgap/)

## Running

Run the HeuristicWeightsLearning class with two arguments: (1) int of population size, (2) int of number of evolutions
`java HeuristicWeightsLearning 1000 50`


## TODO
* make tournament selection or some selection criteria
* make function-ed features
* make feature based GA