# Using Graphgen
Generate {±1}-weighted complete graphs using rudy by G. Rinaldi.

A useful graph generator is provided by G. Rinaldi. 
It makes many types of graphs with machine independent random number generator.
This script downloads it and makes 100 benchmarking graph instances of order N which are uniformly ±1-weighted complete graphs.

You need C compiler. I used gitbash to run this shell command.
You can generate graphs with `./graphgen.sh <graph_order>`, e.g. `./graphgen.sh 2000`. 
