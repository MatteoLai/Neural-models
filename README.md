# Neural-models
This repo contains two separate models:

- "Two_neurons_interaction" simulates the interaction between two neurons, modeling each neuron with the **Integrate & fire** model, which models the cell membrane potential as a lossy integrator circuit (see Hodgkin-Huxley model) that generates a single spike when the membrane potential is above a threshold. To consider the refractory period, the threshold is modeled with a first-order equation.

- "EEGwave_Simulation" generate synthetic signals able to simulate some EEG rhythms, implementing the **Jansen-Rit** model, which considers the interaction between a population of pyramidal neurons and populations of inhibitory and excitatory interneurons.
