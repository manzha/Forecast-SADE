# Forecast-SADE
This algorithm are for financial forecasting. It uses a Simulated Annealing algorithm and Differential Evolution.
It still has lot of work to do and can me improved in different ways, especially sadepar which is a parallel version.
It has hardcoded the population size, and the initial and final temperature, as it was manually tuned.

Compile:
gcc sade.c -o sade.o -lm

For sadepar you must have openMP.
