'''
    This file contains implementation of N Tree.
    We utilised N-Tree to reconstruct the phyogenetic tree and obtain distances
    between nodes as this functionality was not available with Biopython package.
    Used for multichildren tree reconstruction and distance between nodes calculation.
'''


class Node:
    def __init__(self, parent, line):
        self.children = []
        self.parent = parent
        self.line = line
        if line and "name" in line:
            self.name = self.get_name(line)
        else:
            self.name = ""

    def get_name(self, line):
        line = line[line.index("name"):]
        line = line.split("=")
        return line[1][1:-2]


'''
    Builds the tree given directory with XML file or directly passing the XML string.
'''


class NTree:
    def __init__(self, dir=""):
        self.root = None
        self.all_nodes = {}
        self.dir = dir

    def build_tree(self, lines=""):
        if self.dir != "":
            sample_data_file = open(self.dir, 'rw+')
            lines = sample_data_file.readlines()
        lines = lines.split("\n")
        root_line = lines[1]
        root_node = Node(None, root_line)
        self.all_nodes[root_line] = root_node
        previous_node = root_node
        current_offset = 1

        for line in lines[2:]:
            new_offset = line.count("    ")
            step = new_offset - current_offset
            current_offset = new_offset
            if step > 0:
                parent = previous_node
                new_node = Node(parent, line)
                previous_node = new_node
                parent.children.append(new_node)
            elif step == 0:
                new_node = Node(parent, line)
                previous_node = new_node
                parent.children.append(new_node)
            elif step < 0:
                i = 0
                while i != abs(step):
                    parent = parent.parent
                    i += 1
                new_node = Node(parent, line)
                previous_node = new_node
                parent.children.append(new_node)

        self.root = root_node

    '''
        Initialises the traversal of full tree to find all nodes
    '''

    def get_all_nodes(self):
        self._find_all_nodes(self.root)
        return self.all_nodes

    '''
        Does the actual search of all nodes
    '''

    def _find_all_nodes(self, root):
        if root:
            self.all_nodes[root.name] = root
            for child in root.children:
                self._find_all_nodes(child)
        else:
            return None

    '''
        Gets paths to the root
    '''

    def get_path_to_root(self, node):
        path = []
        while node.parent != self.root:
            path.append(node)
            node = node.parent
        path.append(node)
        path.append(self.root)
        return path

    '''
        Finds first common ancestor of two nodes given their path to the root
    '''

    def find_first_ca(self, path_a, path_b):
        path_a = path_a[::-1]
        path_b = path_b[::-1]

        for i in xrange(min(len(path_a), len(path_b))):
            if path_a[i] != path_b[i]:
                return path_a[i - 1]

    '''
        Numbers of steps required to get to teh common ancestor
    '''

    def steps_to_ca(self, ca, path):
        i = 0
        for el in path:
            if el != ca:
                i += 1
            else:
                return i

    '''
        Finds distance between two nodes
    '''

    def get_distance_between_nodes(self, node_a, node_b):
        path_a = self.get_path_to_root(node_a)
        path_b = self.get_path_to_root(node_b)

        common_ancestor = self.find_first_ca(path_a, path_b)
        a_steps = self.steps_to_ca(common_ancestor, path_a)
        b_steps = self.steps_to_ca(common_ancestor, path_b)

        return a_steps + b_steps
