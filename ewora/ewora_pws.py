import numpy as np
import numpy as np
import itertools as it
import scipy.special.binom as binom
import re


class Node:
    """
    This Class describes the Node of genes
    """
    def __init__(self, genes, degree, ng_total=None):

        self._genes = genes
        self._degree = degree
        self._ng_node = len(genes)
        self._ng_total = ng_total

    @property
    def ng_node(self):
        """ Number of member genes"""
        return self._ng_node

    @property
    def degree(self):
        """ Node degree """
        return self._degree

    @property
    def genes(self):
        """ Member genes """
        return self._genes

    @property
    def ng_total(self):
        """ Total number of gene in the analysis """
        return self._ng_total

    @ng_total.setter
    def ng_total(self, ng_total):
        """ Set the total number of gene in the analysis """
        self._ng_total = ng_total

    def get_prob(self, ng_diff):
        """
        Return the probability to be
        :param ng_diff:
        :return:
        """
        return 1 - binom(self.ng_total - self.ng_node, ng_diff) / \
                   binom(self.ng_total, ng_diff)

    def is_affected(self, g_diff):
        intersection = [gene for gene in self.genes if gene in g_diff]
        return True


class Pathway:
    """
    This class describes a pathway with nodes
    """
    def __init__(self, name):
        self._name = name
        # create an empty list of nodes
        self._nodes = []


    @property
    def name(self):
        """ Total number of gene in the analysis """
        return self._name

    @property
    def n_nodes(self):
        return len(self._nodes)

    def add_node(self, genes, degree, ng_total=None):
        self._nodes.append(Node(genes, degree, ng_total))

    def set_ng_total(self, ng_total):
        """
        This function set the total number of nodes
        """
        [node.ng_total(ng_total) for node in self._nodes]

    @property
    def degrees(self):
        """
        This function returns the list of degrees of nodes
        """
        return [node.degree for node in self._nodes]

    def get_probs(self, ng_diff):
        """
        Returns probabilities of nodes to be afected fy the list of DiffEexpr genes
        :param ng_diff: number of differentially expressed genes
        """
        return [node.get_prob(ng_diff) for node in self._nodes]

    def get_affected(self, g_diff):
        """
        Return True-False list whether weach nod is affected fy the list if DE genes
        :param g_diff:
        :return:
        """
        return [node.is_affected(g_diff) for node in self._nodes]



def define_pathway(file_with_pathways = 'data/kegg_formatted.txt'):
    """
    This function is a parser of file with pathways
    Sorry, this function is not flexible to original format
    :param file_with_pathways:
    :return:
    """
    # Does the file exists?

    reading_states = ['Pathway', 'Genes', 'Nodes']
    with open(file_with_pathways, "r") as f:
        # f = open(file_with_pathways, "r")
        # lines = f.readlines()
        # line = lines[3]

        node_cnt = n_nodes = 0
        reading_state = reading_states[0]
        pathways = []
        genes = []
        for line in f:
            if len(line) == 0:
                # if line if empty - get next line
                continue



            # Parser itself
            if reading_state == reading_states[0]:
                # read Pathway name
                name = re.findall(r'"([\w\s-]+)"', line)
                pathways.append(Pathway(name))
                reading_state = reading_states[1]
            elif reading_state == reading_states[1]:
                # read Gene names
                tmp = line.split('\t')
                #tmp.remove('','\n') !!!!!

                genes.append(tmp[2:])
                reading_state = reading_states[2]
            else:
                if n_nodes == node_cnt == 0:
                    n_nodes = int(re.findall(r'\t([\d]+)', line)[0])
                else:
                    # Read nodes while their count reach the number of nodes
                    tmp = line.split('\t')
                    #tmp.remove('', '\n') !!!!!
                    degree = int(tmp[1])
                    genes = tmp[2:]
                    pathway = pathways[-1]
                    # agg node to last pathway
                    pathway.add_node(genes, degree)
                    node_cnt += 1
                    if n_nodes == node_cnt:
                        # stop reading pathway and be ready to read next
                        node_cnt = 0
                        n_nodes = 0
                        reading_state = reading_states[0]

        # Get unique genes
        genes = set(genes)
        ng_total = len(genes)
        [pathway.set_ng_total(ng_total) for pathway in pathways]



























