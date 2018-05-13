from Bio import SeqIO
from graphviz import Digraph
import argparse

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2
        self.coverage = 0
        self.n = 2
    
    def calc_coverage(self,c1,c2):
        self.coverage = c1+c2



class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.k = k
        self.graph = Digraph()
        
    def add_read(self,read):
        read_lng = len(read)
        if read_lng < self.k:
            return
        
        kmer = read[:k]
        
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)

        for next_kmer_indx in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer[-1])
            self.vertices[next_kmer].in_edges[kmer] = new_edge
            self.vertices[kmer].out_edges[next_kmer] = new_edge
            

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
    
    def graphviz(self, view, fgraph):

        for vertex in self.vertices:
            if view == 'full':
                self.graph.node(vertex, label=vertex)
            elif view == 'nick':
                self.graph.node(vertex, label=str(self.vertices[vertex].coverage))
            for next_vertex in self.vertices[vertex].out_edges:
                if view == 'full':
                    self.graph.edge(vertex, next_vertex, label = self.vertices[vertex].out_edges[next_vertex].seq)
                elif view == 'nick':
                    self.graph.edge(vertex, next_vertex, label = str(self.vertices[vertex].out_edges[next_vertex].coverage)+','+
                                    str(self.vertices[vertex].out_edges[next_vertex].n))
                    
        self.graph.render(filename=fgraph, view=True)
    
    def collapse_graph(self):
        
        needed_vertices = []
        for vertex in self.vertices: 
            if len(self.vertices[vertex].in_edges) == 1 and len(self.vertices[vertex].out_edges) == 1:
                needed_vertices.append(self.vertices[vertex].seq)
    
        for i in range(len(needed_vertices)):
            if needed_vertices[i] in self.vertices and len(self.vertices) > 2:
                in_vertex = [j for j in self.vertices[needed_vertices[i]].in_edges.keys()][0]
                out_vertex = [m for m in self.vertices[needed_vertices[i]].out_edges.keys()][0]
                new_edge = Edge(self.vertices[needed_vertices[i]].in_edges[in_vertex].seq, self.vertices[needed_vertices[i]].out_edges[out_vertex].seq[-1])
                self.vertices[in_vertex].out_edges[out_vertex] = new_edge
                self.vertices[out_vertex].in_edges[in_vertex] = new_edge
                # edge length
                self.vertices[in_vertex].out_edges[out_vertex].n = self.vertices[out_vertex].in_edges[needed_vertices[i]].n + self.vertices[in_vertex].out_edges[needed_vertices[i]].n - 1  
                self.vertices[out_vertex].in_edges[in_vertex].n = self.vertices[out_vertex].in_edges[needed_vertices[i]].n + self.vertices[in_vertex].out_edges[needed_vertices[i]].n - 1
                # edge coverage
                self.vertices[in_vertex].out_edges[out_vertex].coverage = self.vertices[in_vertex].out_edges[needed_vertices[i]].coverage + self.vertices[out_vertex].coverage
                self.vertices[out_vertex].in_edges[in_vertex].coverage = self.vertices[in_vertex].out_edges[needed_vertices[i]].coverage + self.vertices[out_vertex].coverage
                
                del self.vertices[out_vertex].in_edges[needed_vertices[i]]
                del self.vertices[in_vertex].out_edges[needed_vertices[i]]
                del self.vertices[needed_vertices[i]]
            else:
                continue
    
    def calc_final_edge_coverage(self):
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex].coverage = self.vertices[current_vertex].out_edges[next_vertex].coverage/self.vertices[current_vertex].out_edges[next_vertex].n
        
    
    def write_fasta(self, fasta, k):
        
        a = 0
        with open(fasta, 'w') as f:
            for vertex in self.vertices:
                for next_vertex in self.vertices[vertex].out_edges:
                    if len(self.vertices[vertex].out_edges[next_vertex].seq) > k+1:
                        a += 1
                        f.write('>')
                        f.write(str(a))
                        f.write('\n')
                        f.write(self.vertices[vertex].out_edges[next_vertex].seq)
                        f.write('\n')
        

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Graph')
    
    parser.add_argument('-i', '--infile', help='fasta file', type=str, required=True)
    parser.add_argument('-outgraph', '--outfile_graph', help='outfile with graph in .dot format', type=str, required=True)
    parser.add_argument('-outfasta', '--outfile_fasta', help='outfile with assembly in .fasta format', type=str, required=True)
    parser.add_argument('-k', '--kmer_size', help='kmer size', type=int, default=15)
    parser.add_argument('-v', '--view', help='choose full or nick view', type=str, required=True)
    parser.add_argument('-c', '--collapse', help='choose full or nick view', action='store_true', default=True)
    
    
    args = parser.parse_args()
    infile = args.infile
    fgraph = args.outfile_graph
    k = args.kmer_size
    view = args.view
    collapse = args.collapse
    fasta = args.outfile_fasta
    my_graph = Graph(k)
    
    with open(infile, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            read = str(record.seq)
            my_graph.add_read(read)
            my_graph.add_read(str(record.reverse_complement().seq))

    my_graph.calc_init_edge_coverage()
    if collapse:
        my_graph.collapse_graph()
    my_graph.calc_final_edge_coverage()
    my_graph.graphviz(view, fgraph)
    my_graph.write_fasta(fasta, k)
