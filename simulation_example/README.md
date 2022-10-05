# Example simulation

It is often useful to know what the expected topology weights are under a given species history.
The python script shows how you can simulate genealogies using [msprime](https://tskit.dev/msprime/docs/stable/intro.html) and the directly compute the weights.
Although visualisation is possible witin python, here we instead export the weights and use the R script to plot them.

The simulation model is defined using [demes](https://popsim-consortium.github.io/demes-docs/main/introduction.html) and [demesdraw](https://grahamgower.github.io/demesdraw/latest/quickstart.html)

You will need to have the above python packages installed, along with the standard dependencies for Twisst.
You will also need to have the `twisst.py` and `plot_twisst.R` scripts in the same directory as these are imported for the analysis.
