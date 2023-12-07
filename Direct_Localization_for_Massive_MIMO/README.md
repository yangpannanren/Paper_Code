# Codes for paper 《Direct Localization for Massive MIMO》

To generate all figures run file `run.sh`.

## Folders

- The folder `Algorithms` contains the 5 positioning algorithms tested in the paper.
- The folder `Routines` contains some minor functions necessary for running the algorithms and executing the figures.

## Some comments on the code files
1. Fig1.m generates a random scenario and estimates the source position. So every time that is run it will generate a different result different than the Fig. 1 in the paper.
2. These are mostly Monte Carlo simulations and as such they take quite some time to run. To increase the speed but also the noise, lower the number of 'MonteCarloRuns'.