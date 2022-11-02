# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
from matplotlib import pyplot as plt
import Bio
from Bio import SeqIO


#the function takes a fasta file and extracts the genome as a string
def extractFasta( path_to_file ):
    sequence = ''
    with open(path_to_file, mode='r') as handle:
        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        for record in SeqIO.parse(handle, 'fasta'):
            # Extract individual parts of the FASTA record
            identifier = record.id
            description = record.description
            sequence = record.seq

            # Example: adapt to extract features you are interested in
            print('----------------------------------------------------------')
            print('Processing the record {}:'.format(identifier))
            print('Its description is: \n{}'.format(description))
            amount_of_nucleotides = len(sequence)
            print('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides))

    sequence = str(sequence)
    print(type(sequence))
    #print(sequence)
    return sequence[0:49]

#the function cuts the sequence to kmers of length k and step
#returns a list with all the sequences
def kmers(k,step,sequence):
    # string_name[start:end:step]
    list1 = []
    for index in range(0,len(sequence)+1):
        #print(sequence[index:index+k:step],'\n')
        list1.append(sequence[index:index+k:step])
    return list1

# This function creates a tree
# Each node is a sequence from the kmersList
def resetTree(tree,kmersList):
    tree.create_node("Idle","idle") # root node
    tree.create_node(kmersList[0], "k1", parent="idle")
    for index in range(2,len(kmersList)):
        tree.create_node(kmersList[index-1], str("k"+str(index)), parent=str("k"+str((index-1))))
    #tree.show()

#def addToTree(tree,kmersList):
    # BFS algorithm

#first refernce of BFS
def bfs(graph, root):

    visited, queue = set(), collections.deque([root])
    visited.add(root)

    while queue:

        # Dequeue a vertex from queue
        vertex = queue.popleft()
        print(str(vertex) + " ", end="")

        # If not visited, mark it as visited, and
        # enqueue it
        for neighbour in graph[vertex]:
            if neighbour not in visited:
                visited.add(neighbour)
                queue.append(neighbour)


def tree_to_graph(self, root, parent, graph):
    graph[root.val] = []

    if (root.left != None):
        graph[root.val].append(root.left.val)
        self.tree_to_graph(root.left, root, graph)

    if (root.right != None):
        graph[root.val].append(root.right.val)
        self.tree_to_graph(root.right, root, graph)

    if (parent != None):
        graph[root.val].append(parent.val)

    return graph

def bfs(self, graph, target, K, level):
    queue = [[target.val, 0]]
    K_list = []
    visited = []

    while (queue):
        node, level = queue.pop(0)
        visited.append(node)

        if (level == K):
            K_list.append(node)
            # print("K_list", K_list)

        elif (level < K):

            if (graph[node] != []):

                for x in graph[node]:
                    if (x not in visited):
                        queue.append([x, level + 1])
                        # print("queue", queue)

    return K_list






if __name__ == '__main__':
    # input_file =
    # output_file =
    # fasta_sequenses = SeqIO.parse(open(input_file),'fasta')
    # with open(output_file) as out_file:
    #     for fasta in fasta_sequences:
    #         name, sequence = fasta.id, str(fasta.seq)
    #         new_sequence = some_function(sequence)
    #         write_fasta(out_file)

     # Install the biopython library (if not already installed)

    # Import parts of Biopython
    from Bio import SeqIO
    #import anytree
    from treelib import *
    import collections

    # File path to your FASTA file
    path_to_file = 'C:/Users/User/PycharmProjects/PhylogeneticTreeBuild/fasta/EPI_ISL_402124.fasta'  # <--- substitute by your local path
    sequence = extractFasta(path_to_file)
    kmers_list = kmers(3,1,sequence)
    print(kmers_list)
    print(type(kmers_list[7]))

    tree = Tree()
    resetTree(tree, kmers_list)
    tree.show()
    #addToTree(tree,)
    #print(bfs(tree,tree.root))
    print(tree.expand_tree())

