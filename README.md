# Collapse graph 
### Required arguments

* ```-i``` - path to infile in  .fasta format
* ```-k``` - kmer size
* ```-v``` - you can choose 'full' or 'nick' view of graph
* ```-c``` - collapse graph or not
* ```-outgraph``` - path to outfile with graph in .dot format
* ```-outfasta``` - path to outfile with assembly in .fasta format (only when you want to collapse graph)

### Running

Open terminal and run:
```
hw5.py -i reads.fasta -k 55 -v nick -c -outgraph my_graph.dot -outfasta my_assembly.fasta
```
or
```
hw5.py -i reads.fasta -k 55 -v full -outgraph my_graph.dot 
```


## Author

* **Adel Gazizova** - *BI student* - [flowdel](https://github.com/flowdel)
