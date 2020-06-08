## Fast and Parallelizable Algorithms for Submodular Maximization

This is the implementation of all algorithms evaluated in the paper,
for the cardinality-constrained maximum cut application (weighted or unweighted)
and revenue maximization application.

### Dependencies 
GNU g++, GNU make utility

### Binaries
. `maxcut`, the main program to run all algorithms for cardinality-constrained maximum cut application
. `revmax`, the main program to run all algorithms for revenue maximization application
. `preproc`, for preprocessing edge lists to custom binary format

### Format of input graph

The programs take as input a custom binary format; the program
`preproc` is provided to convert an undirected edge list format
to the binary format. The first line of the edge list must be of the form
```
<Number of Nodes> <1 if weighted, 0 otherwise> 0
```
Each successive line describes an edge
```
<From id> <To id> <Weight (only if weighted)>
```
The node ids must be nonnegative integers; they are not restricted to lie in any range.
The `Weight` must be an unsigned integer.

If `graph.el` is an edge list in the above format, then
```
preproc graph.el graph.bin
```
generates the `graph.bin` file in the correct binary format for input to the other programs.

### Parameters
```
Options: 
-G <graph filename in binary format>
-k <cardinality constraint>
-A [run AdaptiveSimpleThreshold]
-M [run AdaptiveNonmonotoneMax]
-L [run AdaptiveThresholdGreedy]
-F [run FastInterlaceGreedy]
-T [run IteratedGreedy]
-Q [run FastRandomGreedy]
-B [run Blits]
-o <outputFileName>
-N <repetitions>
-e <epsilon (default 0.1)>
-d <delta (default 0.1)>
-v [verbose]>
-f [fast mode (for algs. using ThresholdSample)]>
-r [report round information]>
```
### Example
```
bash gen-table.bash
```
This command should generate a table similar to the one below,
for cardinality-constrained maxcut
on Barabasi-Albert random graph with n=968 and k=200.
```
           Algorithm            Objective              Queries               Rounds
      IteratedGreedy                 1992               309840                  800
 FastInterlaceGreedy                 1847               133589               139106
    FastRandomGreedy                 1635               173702                  200
   AdaNonmonotoneMax                 1701               476073                   23
     AdaSimpleThresh                 1694               499716                   30
     AdaThreshGreedy                 1978                96104                  157
```