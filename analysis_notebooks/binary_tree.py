'''
    This file contains implementation of Binary Tree.
    We utilised B-Tree to reconstruct the phyogenetic tree and obtain distances
    between nodes as this functionality was not available with Biopython package.
'''


class Node:
    def __init__(self, parent, line):
        self.left = None
        self.right = None
        self.parent = parent
        self.line = line
        if line and "name" in line:
            self.name = self.get_name(line)
        else:
            self.name = ""

    '''
        Extracts node name out of XML line.
    '''

    def get_name(self, line):
        line = line[line.index("name"):]
        line = line.split("=")
        return line[1][1:-2]


class BTree:
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

        # XML processing of the tree
        for line in lines[2:]:
            new_offset = line.count("    ")
            step = new_offset - current_offset
            current_offset = new_offset
            if step > 0:
                parent = previous_node
                left_child = Node(parent, line)
                previous_node = left_child
                parent.left = left_child
            elif step == 0:
                right_child = Node(parent, line)
                previous_node = right_child
                parent.right = right_child
            elif step < 0:
                i = 0
                while i != abs(step):
                    parent = parent.parent
                    i += 1
                right_child = Node(parent, line)
                previous_node = right_child
                parent.right = right_child

        self.root = root_node

    '''
        Traverses the tree and prints all nodes.
    '''

    def _print(self, root):
        if root:
            if root.left:
                self._print(root.left)
            else:
                print root.line
            if root.right:
                self._print(root.right)
            else:
                print root.line

    '''
        Initiates tree traversal and gets all nodes
    '''

    def get_all_nodes(self):
        self._find_all_nodes(self.root)
        return self.all_nodes

    '''
        Actual tree traversal
    '''

    def _find_all_nodes(self, root):

        if root:
            self.all_nodes[root.name] = root
            if root.left:
                self._find_all_nodes(root.left)
            if root.right:
                self._find_all_nodes(root.right)

    '''
        Get the nodes which have to be visited to get to teh root
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
        Finds the first common ancestor for two nodes given their path to the root
    '''

    def find_first_ca(self, path_a, path_b):
        path_a = path_a[::-1]
        path_b = path_b[::-1]

        for i in xrange(min(len(path_a), len(path_b))):
            if path_a[i] != path_b[i]:
                return path_a[i - 1]

    '''
        Get number of steps to the common ancestor.
    '''

    def steps_to_ca(self, ca, path):
        i = 0
        for el in path:
            if el != ca:
                i += 1
            else:
                return i

    '''
        Get the actual distance between two nodes
    '''

    def get_distance_between_nodes(self, node_a, node_b):
        path_a = self.get_path_to_root(node_a)
        path_b = self.get_path_to_root(node_b)

        common_ancestor = self.find_first_ca(path_a, path_b)
        a_steps = self.steps_to_ca(common_ancestor, path_a)
        b_steps = self.steps_to_ca(common_ancestor, path_b)

        return a_steps + b_steps
