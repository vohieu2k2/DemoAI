## Nearly Linear-Time, Deterministic Algorithm for Maximizing (Non-Monotone) Submodular Functions Under Cardinality Constraint

This is the implementation of all algorithms evaluated in the paper,
for the cardinality-constrained maximum cut application (weighted or unweighted).

### Dependencies 
GNU g++, GNU make utility

### Binaries
. `maxcut`, the main program to run all algorithms
. `preproc`, for preprocessing edge lists to custom binary format
. `er`, to generate Erdos-Renyi edges lists
. `ba`, to generate Barabasi-Albert scale-free edge lists

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

### Parameters of `maxcut`
```
-G <graph filename in binary format>
-k <cardinality constraint>
-I [run InterlaceGreedy]
-F [run FastInterlaceGreedy]
-T [run Gupta et al. (2010)]
-R [run RandomizedGreedy]
-Q [run FastRandomizedGreedy]
-B [run Blits]
-l [turn off stealing behavior of (Fast)InterlaceGreedy]
-o <outputFileName>
-N <repetitions>
-e <epsilon (default 0.1)>
```
### Example
```
./er 1000 0.5 er1k.el #creates ER edge list with n = 1000, p=0.5
./preproc er1k.el er1k.bin #creates binary input file
./maxcut -G er1k.bin -k 50 -F #Runs FastInterlaceGreedy for cardinality-constrained maxcut problem
```