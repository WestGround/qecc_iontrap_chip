from rustworkx import PyDiGraph
from Node import Node

'''
    Logical quantum circuit
'''
class Circuit:
    def __init__(self, num_qubits):
        self.num_qubits = num_qubits
        self.dag = PyDiGraph()
        self.qubit_list = []
        self.last_gates = []

        for qubit in range(num_qubits):
            node_index = self.dag.add_node(Node('qubit', [qubit], []))
            self.qubit_list.append(node_index)
            self.last_gates.append(node_index)

    def append_gate(self, name, qubits, parameter=[], duration=1, clifford=True):
        node = Node(name, qubits, parameter, duration, clifford)
        self.append_node(node)

    def append_node(self, node):
        node_index = self.dag.add_node(node)
        for qubit in node.qubits:
            last_gate = self.last_gates[qubit]
            self.dag.add_edge(last_gate, node_index, None)
            self.last_gates[qubit] = node_index

    def get_next_gate(self, qubit_index: int, node_index: int) -> int:
        next_candidate = []
        for next_index in self.dag.successor_indices(node_index):
            if qubit_index in self.dag[next_index].qubits:
                next_candidate.append(next_index)

        if len(next_candidate)==1:
            return next_candidate[0]
        elif len(next_candidate)==2:
            assert(self.dag[node_index].name=='cx' or self.dag[node_index].name=='ms')
            gate1 = self.dag[next_candidate[0]]
            gate2 = self.dag[next_candidate[1]]
            ### case 1) not ambiguous
            if qubit_index not in gate1.qubits:
                assert(qubit_index in gate2.qubits)
                return next_candidate[1]
            elif qubit_index not in gate2.qubits:
                assert(qubit_index in gate1.qubits)
                return next_candidate[0]
            ### case 2) ambiguous with different gate
            elif next_candidate[0]!=next_candidate[1]:
                curr_node = self.dag[node_index]
                if curr_node.qubits.index(qubit_index)==0:
                    pair_index = curr_node.qubits[1]
                else:
                    pair_index = curr_node.qubits[0]
                if qubit_index in gate1.qubits and pair_index in gate1.qubits:
                    return next_candidate[1]
                else:
                    assert(qubit_index in gate2.qubits and pair_index in gate2.qubits)
                    return next_candidate[0]
            ### case 3) ambiguous with same gate
            else:
                return next_candidate[0]
        else:
            assert(len(next_candidate)==0)
            return None

    def get_prev_gate(self, qubit_index: int, node_index: int) -> int:
        prev_candidate = []
        for prev_index in self.dag.predecessor_indices(node_index):
            if qubit_index in self.dag[prev_index].qubits:
                prev_candidate.append(prev_index)

        if len(prev_candidate)==1:
            return prev_candidate[0]
        elif len(prev_candidate)==2:
            gate1 = self.dag[prev_candidate[0]]
            gate2 = self.dag[prev_candidate[1]]
            ### case 1) 2Q gate preceded by 1Q gate and 2Q gate
            if prev_candidate[0]!=prev_candidate[1]:
                curr_node = self.dag[node_index]
                pair_index = curr_node.qubits[(curr_node.qubits.index(qubit_index)+1)%2]

                if qubit_index in gate1.qubits and \
                  pair_index in gate1.qubits:
                    return prev_candidate[1]
                else:
                    assert(qubit_index in gate2.qubits and pair_index in gate2.qubits and pair_index not in gate1.qubits)
                    return prev_candidate[0]
            ### case 2) 2Q gate preceded by another 2Q gate that is applied to the same two qubits
            else:
                return prev_candidate[0]
        else:
            assert(len(prev_candidate)==0)
            return None
        
    def print_circuit(self, abstract: bool):
        for qubit_index in range(self.num_qubits):
            qubit_line = "q%d"%(qubit_index)
            curr_index = self.get_next_gate(qubit_index, qubit_index)
            while curr_index!=None:
                curr_node = self.dag[curr_index]
                if abstract:
                    if curr_node.name in ['x', 'y', 'z', 'h', 's']:
                        qubit_line += " %s"%curr_node.name
                    elif curr_node.name=='sdg':
                        qubit_line += " d"
                    elif curr_node.name=='rz':
                        qubit_line += " r"
                    elif curr_node.name=='cx':
                        qubit_line += " c"
                    else:
                        raise Exception(f"Invalid gate name %s"%curr_node.name)
                else:
                    if curr_node.name=='gpi':
                        qubit_line += " g"
                    elif curr_node.name=='gpi2':
                        qubit_line += " p"
                    elif curr_node.name=='ms':
                        qubit_line += " m"
                    elif curr_node.name=='rz' or curr_node.name=='vz':
                        qubit_line += " r"
                    else:
                        raise Exception(f"Invalid gate name %s"%curr_node.name)
                curr_index = self.get_next_gate(qubit_index, curr_index)
            #qubit_line += '\n'
            print(qubit_line)

    def circuit_to_txt(self, filepath: str, qec: bool, approx_factor: int) -> int:
        output_file = open(filepath, 'w+')
        output_file.write("%d %d"%(self.num_qubits, approx_factor))
        gate_cnt = 0

        qubit_stack = []
        pending_gate = []
        for qubit_index in range(self.num_qubits):
            pending_gate.append(self.get_next_gate(qubit_index, qubit_index))

        for qubit_index in range(self.num_qubits):
            curr_qubit_index = qubit_index
            curr_index = pending_gate[curr_qubit_index]

            while curr_index!=None:
                curr_node = self.dag[curr_index]
                if curr_node.name in ['gpi', 'gpi2']:
                    output_file.write("\n%s %d"%(curr_node.name, curr_qubit_index))
                    gate_cnt += 1
                    next_index = self.get_next_gate(curr_qubit_index, curr_index)
                    pending_gate[curr_qubit_index] = next_index
                    curr_index = next_index
                elif curr_node.name=='rz' and qec:
                    assert(isinstance(curr_node.parameter[0], int))
                    output_file.write("\n%s %d %d"%(curr_node.name, curr_qubit_index, curr_node.parameter[0]))
                    gate_cnt += 1
                    next_index = self.get_next_gate(curr_qubit_index, curr_index)
                    pending_gate[curr_qubit_index] = next_index
                    curr_index = next_index
                elif curr_node.name=='rz' and not qec:
                    assert(isinstance(curr_node.parameter[0], float))
                    output_file.write("\n%s %d"%(curr_node.name, curr_qubit_index))
                    gate_cnt += 1
                    next_index = self.get_next_gate(curr_qubit_index, curr_index)
                    assert(next_index==None)
                    pending_gate[curr_qubit_index] = next_index
                    curr_index = next_index
                elif curr_node.name=='ms':
                    pair_qubit_index = curr_node.qubits[(curr_node.qubits.index(curr_qubit_index)+1)%2]
                    if qubit_stack and qubit_stack[-1]==(pair_qubit_index, curr_index):
                        output_file.write("\n%s %d %d"%(curr_node.name, curr_node.qubits[0], curr_node.qubits[1]))
                        qubit_stack.pop()
                        gate_cnt += 1
                        pending_gate[curr_qubit_index] = self.get_next_gate(curr_qubit_index, curr_index)
                        pending_gate[pair_qubit_index] = self.get_next_gate(pair_qubit_index, curr_index)
                    else:
                        qubit_stack.append((curr_qubit_index, curr_index))
                    curr_qubit_index = pair_qubit_index
                    curr_index = pending_gate[pair_qubit_index]
                else:
                    raise Exception(f'Invalid gate instruction %s %s'%(curr_node.name, curr_node.qubits))
                
        for qubit_index in range(self.num_qubits):
            assert(pending_gate[qubit_index]==None)

        output_file.close()

        return gate_cnt