# Negative Probability: Forward-Backward Algorithm

Here is the project's website for our paper titled "Negative Probabilities in Gene Regulatory Networks" 

The paper can be accessed at [insert link to the paper].


In this project, we introduce a salient algorithm: the Forward-Backward (Sinkhorn-like) algorithm, which seeks the transition matrix with the allowance of negative transition probability. Additionally, in order to verify the results/minimizer, we also provide the gradient descent method.

To implement the two methods independently, see the code "forward_backward.m" for the Forward-Backward method and "gradient_descent.m" for gradient descent. We provide four examples for each of the methods. Specifically, the examples include three/four/then-node networks and a random large-scale network with a user-defined number of nodes (100 in the default setting, its connectivity is parameterized by $x$). The convergence comparison can also be found in "convergence.m".


i) The data for the large-scale example can be found in the "example_data" folder to reproduce the experiments conducted in the paper.

ii) All the functions used in the paper are available in the "function" folder.

Enjoy!

![alt text](https://github.com/dytroshut/negative-probability-forward-backward/blob/main/gene_network.png)
