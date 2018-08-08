"""
Coalescent simuation of a pathogen within a host tree, given a constant
effective population size within hosts and a transmission bottleneck to
one lineage.
Assume one pathogen lineage sampled per host.
"""

# TODO:  transmission bottleneck should only affect one of two hosts (H1, H2)
# TODO:  Pick one host at random, and force coalescence of all its pathogen lineages
# TODO:  Then take that ancestor and carry it over to ancestor of H1 and H2 (H1,2)
# TODO:  Carry over ALL lineages from the other host into H1,2

from ete3 import Tree, TreeNode
import argparse
import random
from itertools import combinations
from sys import float_info

class NestedCoalescent():
    def __init__(self, hosttree, c_rate, m_rate, p_cospec, verbose=False):
        self.verbose = verbose
        self.hosttree = hosttree

        self.c_rate = c_rate  # coalescent rate
        self.m_rate = m_rate  # migration (host-switch) rate
        self.p_cospec = p_cospec  # probability of co-speciation

        # pathogen objects
        self.extant_p = []  # active pathogen lineages
        self.not_yet_sampled_p = []  # pathogen lineages higher in the tree
        self.not_extant_p = []

        self.extant_h = []  # sampled hosts
        self.not_yet_sampled_h = []  # host tips higher in the tree
        self.not_extant_h = []  # host nodes that have coalesced

        # simulation time frame is defined by inter-node intervals in host tree
        self.host_nodes = []

        self.locations = {}  # {host: [p1, p2, ...]}
        self.choices = {}  # {(p1, p2}: host}

        # process host tree
        self.get_node_heights()


    def __str__(self):
        """ Display contents of extant lists for debugging """
        output = 'Pathogens:\n'
        for x in self.extant_p:
            output += '  n %s  h %f  d %f  host %s  extant\n' % (x.name, x.height, x.dist, x.host.name)
        for x in self.not_extant_p:
            output += '  n %s  h %f  d %f  host %s  not extant\n' % (x.name, x.height, x.dist, x.host.name)
        for x in self.not_yet_sampled_p:
            output += '  n %s  h %f  d %f  host %s  not yet sampled\n' % (x.name, x.height, x.dist, x.host.name)

        output += 'Hosts:\n'
        for x in self.extant_h:
            output += '  %s %f %f  extant\n' % (x.name, x.height, x.dist)
        for x in self.not_extant_h:
            output += '  %s %f %f  not extant\n' % (x.name, x.height, x.dist)
        for x in self.not_yet_sampled_h:
            output += '  %s %f %f  not yet sampled\n' % (x.name, x.height, x.dist)

        return output

    def get_nonterminals(self):
        """
        Return an iterator over all internal nodes of a Tree object
        :param tree:
        :return:
        """
        for node in self.hosttree.traverse():
            if not node.is_leaf():
                yield node


    def get_node_heights(self):
        """
        Calculate node heights for all nodes of the tree.
        Annotate nodes with heights in-place.
        :param tree:  ete3 Tree object
        """

        # reset lists
        self.extant_h = []
        self.not_extant_h = []
        self.not_yet_sampled_h = []

        # the total height of the tree is the maximum distance from root to any tip
        root = self.hosttree.get_tree_root()
        apex = self.hosttree.get_farthest_leaf()[0]
        total_height = self.hosttree.get_distance(root, apex)

        # iterate through all nodes in the tree
        self.host_nodes = []
        for node in self.hosttree.traverse():
            # the node's depth is its distance from the root
            depth = self.hosttree.get_distance(root, node)
            # the node's height is the difference between its depth and the total height
            height = total_height-depth
            node.add_feature('height', height)  # modify TreeNode in place            

            # we do not allow zero branch length
            # node.dist = node.dist + (np.finfo(float).eps)
            if (node.dist == 0):
                greaterThan = False
                while greaterThan == False:
                    node.dist += float_info.epsilon  # add a small amount
                    if (node.height+node.dist<=node.height):
                        greaterThan = False
                    else:
                        greaterThan = True
                node.height += node.dist

            if node.is_leaf():
                # keep track of leaf nodes
                if node.height == 0:
                    self.extant_h.append(node)
                    continue  # do not append to host_nodes list
                else:
                    self.not_yet_sampled_h.append(node)

            # store the host node and its height
            self.host_nodes.append((node.height, node))

        self.host_nodes.sort()  # order by ascending node height


    def initialize_pathogen_tree(self):
        """
        Initialize one pathogen lineage per host tip
        dist records height that pathogen lineage was started
        TODO: relax this assumption - needs some way to input
        """
        # reset containers
        self.extant_p = []  # pathogen lineages that have not coalesced
        self.not_yet_sampled_p = [] # pathogen lineages higher in the tree
        for i, host_tip in enumerate(self.hosttree.get_leaves()):
            pnode = TreeNode(name=host_tip.name+'_P', dist=0)
            pnode.add_features(height=host_tip.height, host=host_tip)
            if host_tip.height == 0:
                self.extant_p.append(pnode)
            else:
                self.not_yet_sampled_p.append(pnode)

        # FIXME: may be clearer to explicitly call each function at a higher level,
        # FIXME: instead of having a chain of function calls like we have now (AP)
        #self.get_pairs()  # also calls get_locations()


    def choose2(self, n):
        """ Shorthand for n choose 2 - probably overkill """
        return n * (n-1) / 2.


    def get_locations(self):
        """
        Collect host locations of all extant pathogen lineages into dictionary of
          host1: [pathogen1, pathogen2, ...]  key-value pairs
        """
        self.locations = {}  # reset dictionary
        for node in self.extant_p:
            if node.host not in self.locations:
                self.locations.update({node.host: []})
            self.locations[node.host].append(node)


    def get_pairs(self):
        """
        Extract all pairs of pathogen lineages that may coalesce.
        """
        self.get_locations()
        self.choices = {}
        for host, pathogens in self.locations.iteritems():
            if len(pathogens) > 1:
                for pair in combinations(pathogens, 2):
                    self.choices.update({pair: host})  # pairs of pathogens in same host


    def coalesce_paths(self, child_paths, t0):
        """
        Create a new TreeNode and assign a given list of child nodes and its host node.
        :param child_paths:  A list of TreeNodes in the pathogen tree.
        :param t0:  Time of pathogen coalescence as height
        :return:  A tuple containing:
            1. TreeNode object for the new pathogen lineage.
            2. updated extant list
        """
        assert len(child_paths)==2, 'Can only coalesce 2 pathogen lineages at a time'
        p1, p2 = child_paths

        assert p1 in self.extant_p and p2 in self.extant_p, 'Both pathogen lineages must be extant'
        assert p1.host == p2.host, 'Can only coalesce pathogen lineages in the same host'
        host = p1.host

        assert p1.height < t0 and p2.height < t0, \
            'Pathogen lineage heights %f %f cannot exceed coalescent event %f' % (p1.height, p2.height, t0)

        # create new pathogen lineage
        new_path = TreeNode(name='_'.join([x.name for x in child_paths]), dist=0)
        new_path.add_features(host=host, height=t0)

        # cast child_paths as a List because ete3.Tree.children requires it
        new_path.children = list(child_paths)
        self.extant_p.append(new_path)

        # coalesced pathogen lineages are no longer extant
        for node in child_paths:
            node.up = new_path
            node.dist = t0 - node.height  # when node was created, we stored the height
            self.extant_p.remove(node)
            self.not_extant_p.append(node)

        return new_path


    def coalesce_hosts(self, host_node):
        """
        With coalescence of host lineages, transfer pathogen lineages carried by
        the respective hosts to the ancestor.

        :param host_node:  TreeNode object representing ancestral host
        :param extant:  List of extant pathogen lineages
        :param prob:  Probability of cospeciation.
        :return:  A TreeNode object for the new pathogen lineage
        """
        c_hosts = host_node.children
        assert len(c_hosts) == 2, "Error: function assumes binary tree"
        h1, h2 = c_hosts

        # label the ancestral host node
        host_node.name = h1.name + '_' + h2.name

        # extract affected pathogen lineages in each descendant (H1, H2) of the host node
        p1 = filter(lambda x: x.host == h1, self.extant_p)
        p2 = filter(lambda x: x.host == h2, self.extant_p)

        if len(p1)>0 and len(p2) > 0 and random.uniform(0,1) < self.p_cospec:
            # cospeciation - pathogen lineages carried in H1 and H2 coalesce in host node
            # TODO: What if there are multiple pathogen lineages in H1 and/or H2?
            # Possibilities:  (1) select one random pair of pathogen lineages to coalesce (only 1 cospeciation)
            #  (2) every pair of pathogen lineages in H1/H2 has probability of cospeciation
            #      This makes it possible for 3 or more path. lineages to coalesce at once
            # Current implementation (below) assumes (1).

            pick1 = random.sample(p1, 1)[0]  # returns a list
            pick2 = random.sample(p2, 1)[0]
            pick1.host = host_node  # relocate these pathogen lineages to ancestral host
            pick2.host = host_node
            to_coalesce = [pick1, pick2]
            self.coalesce_paths(to_coalesce, t0=host_node.height)

        # carry over all lineages to the ancestral host
        for node in p1+p2:
            node.host = host_node

        # update host lists
        self.extant_h.remove(h1)
        self.not_extant_h.append(h1)
        self.extant_h.remove(h2)
        self.not_extant_h.append(h2)
        self.extant_h.append(host_node)


    def switch_hosts(self, t0, seed=None):
        """
        Select an extant pathogen lineage at random and reassign its host
        :return:
        """
        assert len(self.extant_h) > 1, "Error: attempted to switch between one host"
        if seed:
            random.seed(seed)
        pick_p = random.choice(self.extant_p)  # select an extant pathogen lineage at random
        pick_h = pick_p.host
        while pick_h == pick_p.host:
            pick_h = random.choice(self.extant_h)

        # add a node of degree size 2 to annotate host switch event in tree
        pick_p.dist = t0 - pick_p.height
        next_p = Tree(name=pick_p.name+'_m%s-%sm' % (pick_p.host.name, pick_h.name), dist=0)
        next_p.add_features(host=pick_h, height=t0)
        pick_p.up = next_p
        next_p.children = [pick_p]

        self.extant_p.remove(pick_p)
        self.extant_p.append(next_p)
        self.not_extant_p.append(pick_p)

    def coalesce_within_root(self, host_node):
        """
        Coalesce all pathogen lineages in the root branch of the host tree until
        there is only one left.
        """
        height = host_node.height
        while len(self.extant_p) > 1 and len(self.choices) >= 1:
            if self.verbose:
                print self
            self.get_pairs()
            if len(self.choices) == 0:
                #
                return
            height += random.expovariate(len(self.choices)*self.c_rate)
            cpaths = random.choice(self.choices.keys())
            self.coalesce_paths(cpaths, t0=height)
        if self.verbose:
            print self


    def simulate(self, seed=None):
        """
        Main function.  Perform a coalescent simulation of lineages within
        the host tree provided as the first argument.  This implementation assumes
        that all host and pathogen tips have been sampled at the same time.
        :return:  an ete3 tree representing a pathogen tree
        """

        # reset lists
        self.extant_p = []
        self.not_extant_p = []
        self.not_yet_sampled_p = []
        self.extant_h = []
        self.not_extant_h = []
        self.not_yet_sampled_h = []

        if seed:
            random.seed(seed)  # for debugging and unit testing

        # populates host_nodes list
        self.get_node_heights()  # extant_h initialized with zero-height host tips
        next_height, next_node = self.host_nodes.pop(0)

        # populate host tips with pathogen lineages
        self.initialize_pathogen_tree()

        this_height = 0  # initialize at the end of the tree (most recent tip)
        while True:
            if self.verbose:
                # report on state
                print self

            # update pathogen locations (hosts) and record pairs
            self.get_pairs()

            # total rate of pathogen coalescence or host-switch events within this interval
            lambd_mig = len(self.extant_p) * self.m_rate if len(self.extant_h) > 1 else 0.
            lambd_tot = len(self.choices) * self.c_rate + lambd_mig

            # draw waiting time
            wait = random.expovariate(lambd_tot) if lambd_tot > 0 else None
            if wait is None or wait > (next_height - this_height):
                # waiting time exceeds next host node height
                if next_node.is_leaf():
                    if self.verbose:
                        print ('sampled host tip', next_node.name)
                    self.extant_h.append(next_node)
                    self.not_yet_sampled_h.remove(next_node)

                    # activate pathogen lineages in this host tip
                    nodes_to_move = filter(lambda x: x.host==next_node, self.not_yet_sampled_p)
                    self.extant_p.extend(nodes_to_move)
                    for node in nodes_to_move:
                        self.not_yet_sampled_p.remove(node)
                else:
                    if self.verbose:
                        print 'coalesce hosts at', next_node
                    self.coalesce_hosts(next_node)

                this_height = next_height  # move to top of current interval
                if len(self.host_nodes) == 0:
                    break
                next_height, next_node = self.host_nodes.pop(0)  # retrieve next host node
                continue

            # ELSE there is either a host-switch or pathogen coalescence event
            if random.uniform(0, 1) < lambd_mig/lambd_tot:
                if self.verbose:
                    print 'switch hosts'
                self.switch_hosts(t0=this_height+wait)  # host-switch event
            else:
                # coalescence of pathogen lineages within a host
                if self.verbose:
                    print 'coalesce within hosts'
                c_paths = random.choice(self.choices.keys())  # randomly select a pair that can coalesce
                self.coalesce_paths(c_paths, t0=this_height+wait)

            this_height += wait  # move up the tree

        # coalesce remaining pathogen lineages in the last host
        self.get_pairs()
        self.coalesce_within_root(next_node)

        return self.extant_p[0].get_tree_root()  # should correspond to root of pathogen tree



def main():
    parser = argparse.ArgumentParser(description='Coalescent simulation of pathogen tree within host tree.')
    parser.add_argument('tree', type=argparse.FileType('rU'),
                        help='<input> File containing Newick string representation of host tree')
    parser.add_argument('L', type=float, help='<input> Coalescence rate of pathogen within host')
    parser.add_argument('M', type=float, help='<input> Migration (host-switch) rate between hosts')
    parser.add_argument('P', type=float, help='<input> Probability of co-speciation')
    parser.add_argument('outfile', type=argparse.FileType('w'), help='<output> File to write Newick output including internal nodes.')
    parser.add_argument('outfile2', type=argparse.FileType('w'), help='<output> File to write Newick output excluding internal nodes.')
    parser.add_argument('--verbose', action='store_true', help='<option> verbose output')
    parser.add_argument('-seed', type=int, help='<optional> set random seed', default=None)
    parser.add_argument('-n', type=int, help='Number of replicate simulations, defaults to 1', default=1)

    args = parser.parse_args()

    hosttree = Tree(args.tree.name)
    nc = NestedCoalescent(hosttree, args.L, args.M, args.P, verbose=args.verbose)

    for rep in range(args.n):
        simtree = nc.simulate(seed=args.seed)
        args.outfile.write(simtree.write(format=1)+"\n")
        args.outfile2.write(simtree.write(format=5)+"\n")

if __name__ == '__main__':
    main()
