# CS4106 Genetic Algorithm
### [Adam Aherne](https://github.com/underwaterjesus "Adam's GitHub") - [Alan Finnin](https://github.com/alanfinnin "Alan's GitHub") - [Daniel Dalton](https://github.com/ddalton98 "Daniel's GitHub") - [William Cummins](https://github.com/Willc200 "William's GitHub")
CS4106 Genetic Algorithm project. Use a genetic algorithm to minimize graph edge crossings.<br/>
## Project Aim
Placing all nodes of a graph on the perimeter of a circle, order the nodes to have the lowest number of edge crossings.
## Description
**For more a more in-depth description see [CS4106_Project_20.pdf](https://github.com/underwaterjesus/CS4106_Genetic_Algorithm/blob/master/CS4106_Project_20.pdf
"CS4106 Project Specification") in the main git repository.<br/>
Final deliverable named [is12159603.java](https://github.com/underwaterjesus/CS4106_Genetic_Algorithm/blob/master/is12159603.java "is12159603.java"). This is to match project specifications.**<br/>
Given an input file, input.txt, representing an edge list, we first create an adjacency matrix.
The user is then asked to input the population size(P), number of generations(G), crossover rate(Cr), and mutation rate(Mu).
These must all be validated and checked against constraints in the project specification.<br/>
The user must also select between two fitness functions to use.<br/>
The program cycles through the genetic process for the given number of generations,
outputting the best performer to the console each generation. The graph will also be visually presented to the user in a JFrame window.<br/>
A member of the population is an array representing the order that the nodes are to be placed on the circle.
## Fitness Functions
The user will be able to pick from two different fitness functions. One is given in
the project sepcification, the other has been chosen from relevant academic literature.<br/><br/>
The fitness function described in the specification is the combined length of all edges in the graph.<br/>
The lower the total, the higher the fitness of the individual. While the function was described, no implementation details were given.
Writing this algorithm was part of the project.<br/><br/>
The second fitness function we have chosen to use is an implementation of the **AngGA** algorithm.
More details can be found [here](https://ulir.ul.ie/bitstream/handle/10344/5395/Eaton_2016_GA.pdf?sequence=1
"A GA-Inspired Approach to the Reduction of EdgeCrossings in Force-Directed Layouts").
This paper is by Farshad Toosi(farshad.toosi@cit.ie), Nikola Nikolov(nikola.nikolov@ul.ie) and Malachy Eaton(malachy.eaton@ul.ie). All contact details correct as of 12th May 2020.<br/>
The lower the total of this function, the higher the fitness of the individual. To implement this algorithm it was necessary to calculate the
angles formed by all connected edges, the length of these edges, and the betweeness centrality of each node.
## Graph Visualization
All graphs are visualized using an algorithm given in the project specification.
## Results
Current testing has produced a best result of one(1) edge crossing by the fitness function specified in the project, and a best result of one(1) edge crossing by the AngGA algorithm.
