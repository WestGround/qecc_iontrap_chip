import math

from rustworkx.visualization import graphviz_draw

from Circuit import Circuit
from Node import Node

class AbsCircuitOptimizer:
    def __init__(self, circuit: Circuit):
        self.circuit = circuit

        self.h_rule = []
        self.cancel_rule = []
        self.merge_rule = []
        self.commute_rule = []
        self.initialize_h_ruleset()
        self.initialize_cancelling_ruleset()
        self.initialize_merging_ruleset()
        #self.initialize_commuting_ruleset()

    def abs_opt(self) -> Circuit:
        optimized = True

        cnt = 0
        while optimized:
            cnt += 1
            optimized = False
            if cnt%3==0:
                print("Start %dth H reduction"%cnt)
            optimized |= self.H_reduction()
            if cnt%3==0:
                print("Start %dth gate cancellation"%cnt)
            optimized |= self.gate_cancellation()

        return self.circuit

    def H_reduction(self) -> bool:
        result = False
        for qubit_index in range(self.circuit.num_qubits):
            curr_index = self.circuit.get_next_gate(qubit_index, qubit_index)
            while curr_index!=None:
                next_index = self.circuit.get_next_gate(qubit_index, curr_index)
                ### Pattern Match
                result |= self.update_h_rule(qubit_index, curr_index)

                curr_index = next_index
            self.reset_h_pointer()

        return result

    def initialize_h_ruleset(self):
        self.h_rule.append([['h', 's', 'h'], ['sdg', 'h', 'sdg'], 0])                    ### rule 0
        self.h_rule.append([['h', 'sdg', 'h'], ['s', 'h', 's'], 0])                     ### rule 1
        self.h_rule.append([['h', 's', 'cx', 'sdg', 'h'], ['s', 'cx', 'sdg'], 0])        ### rule 2
        self.h_rule.append([['h', 'sdg', 'cx', 's', 'h'], ['sdg', 'cx', 's'], 0])        ### rule 3
        self.h_rule.append([['h', 'cx', 'h', 'h', 'h'], ['cx'], 0])                    ### rule 4

    def reset_h_pointer(self):
        for rule_index in range(len(self.h_rule)):
            self.h_rule[rule_index][-1] = 0

    def update_h_rule(self, qubit_index: int, node_index: int) -> bool:
        opt_flag = False
        curr_node = self.circuit.dag[node_index]
        for rule_index in range(len(self.h_rule)):
            before_opt = self.h_rule[rule_index][0]
            after_opt_rev = self.h_rule[rule_index][1]
            pointer = self.h_rule[rule_index][2]

            if curr_node.name==before_opt[pointer]:
                ### Additional Check for rule 2, 3
                if (rule_index==2 or rule_index==3) and pointer==2:
                    assert(curr_node.name=='cx')
                    if curr_node.qubits.index(qubit_index)==0:
                        pointer = 0
                    else:
                        pointer += 1
                else:
                    pointer += 1
            else:
                pointer = 0

            ### Additional Check for rule 4
            if rule_index==4 and pointer==3:
                tqgate_index = self.circuit.get_prev_gate(qubit_index, node_index)
                tqgate = self.circuit.dag[tqgate_index]
                assert(tqgate.name=='cx')
                pair_index = tqgate.qubits[(tqgate.qubits.index(qubit_index)+1)%2]
                pair_prev_index = self.circuit.get_prev_gate(pair_index, tqgate_index)
                pair_prev_node = self.circuit.dag[pair_prev_index]
                pair_next_index = self.circuit.get_next_gate(pair_index, tqgate_index)
                if pair_prev_node.name=='h' and pair_next_index!=None and self.circuit.dag[pair_next_index].name=='h':
                    pointer += 2
                else:
                    pointer = 0

            if pointer==len(before_opt) and (rule_index==0 or rule_index==1):
                opt_flag = True
                curr_index = node_index
                for gate_type in after_opt_rev:
                    curr_node = self.circuit.dag[curr_index]
                    curr_node.name = gate_type
                    curr_index = self.circuit.get_prev_gate(qubit_index, curr_index)
            elif pointer==len(before_opt) and (rule_index==2 or rule_index==3):
                opt_flag = True
                ### Remove the current H gate
                prev_index = self.circuit.get_prev_gate(qubit_index, node_index)
                assert(self.circuit.dag[node_index].name=='h')
                self.circuit.dag.remove_edge(prev_index, node_index)
                next_index = self.circuit.get_next_gate(qubit_index, node_index)
                if next_index!=None:
                    self.circuit.dag.remove_edge(node_index, next_index)
                    self.circuit.dag.add_edge(prev_index, next_index, None)
                self.circuit.dag.remove_node(node_index)
                ### Convert S, sdg to sdg, S respectively
                curr_index = prev_index
                for gate_type in after_opt_rev:
                    curr_node = self.circuit.dag[curr_index]
                    if gate_type=='cx':
                        assert(curr_node.name=='cx')
                    curr_node.name = gate_type
                    curr_index = self.circuit.get_prev_gate(qubit_index, curr_index)
                ### Remove the H gate
                assert(self.circuit.dag[curr_index].name=='h')
                next_index = self.circuit.get_next_gate(qubit_index, curr_index)
                self.circuit.dag.remove_edge(curr_index, next_index)
                prev_index = self.circuit.get_prev_gate(qubit_index, curr_index)
                if prev_index!=None:
                    self.circuit.dag.remove_edge(prev_index, curr_index)
                    self.circuit.dag.add_edge(prev_index, next_index, None)
                self.circuit.dag.remove_node(curr_index)
            elif pointer==len(before_opt) and rule_index==4:
                assert(pointer==5)
                opt_flag = True
                cnot_index = self.circuit.get_prev_gate(qubit_index, node_index)
                cnot_gate = self.circuit.dag[cnot_index]
                assert(cnot_gate.name=='cx')
                pair_index = cnot_gate.qubits[(cnot_gate.qubits.index(qubit_index)+1)%2]

                h2_index = node_index
                h1_index = self.circuit.get_prev_gate(qubit_index, cnot_index)

                h3_index = self.circuit.get_prev_gate(pair_index, cnot_index)
                h4_index = self.circuit.get_next_gate(pair_index, cnot_index)

                assert(self.circuit.dag[h1_index].name=='h' and self.circuit.dag[h2_index].name=='h' 
                  and self.circuit.dag[h3_index].name=='h' and self.circuit.dag[h4_index].name=='h')

                self.circuit.dag.remove_edge(h1_index, cnot_index)
                h1_prev_index = self.circuit.get_prev_gate(qubit_index, h1_index)
                self.circuit.dag.add_edge(h1_prev_index, cnot_index, None)
                
                self.circuit.dag.remove_edge(cnot_index, h2_index)
                h2_next_index = self.circuit.get_next_gate(qubit_index, h2_index)
                if h2_next_index!=None:
                    self.circuit.dag.remove_edge(h2_index, h2_next_index)
                    self.circuit.dag.add_edge(cnot_index, h2_next_index, None)

                self.circuit.dag.remove_edge(h3_index, cnot_index)
                h3_prev_index = self.circuit.get_prev_gate(pair_index, h3_index)
                self.circuit.dag.add_edge(h3_prev_index, cnot_index, None)

                self.circuit.dag.remove_edge(cnot_index, h4_index)
                if self.circuit.dag.successor_indices(h4_index):
                    h4_next_index = self.circuit.dag.successor_indices(h4_index)[0]
                    self.circuit.dag.remove_edge(h4_index, h4_next_index)
                    self.circuit.dag.add_edge(cnot_index, h4_next_index, None)

                self.circuit.dag.remove_node(h1_index)
                self.circuit.dag.remove_node(h2_index)
                self.circuit.dag.remove_node(h3_index)
                self.circuit.dag.remove_node(h4_index)

                self.circuit.dag[cnot_index].qubits[0], self.circuit.dag[cnot_index].qubits[1] = \
                  self.circuit.dag[cnot_index].qubits[1], self.circuit.dag[cnot_index].qubits[0]

            if opt_flag:
                self.reset_h_pointer()
                break
            else:
                self.h_rule[rule_index][2] = pointer

        return opt_flag

    def gate_cancellation(self) -> bool:
        for qubit_index in range(self.circuit.num_qubits):
            opt_flag = False

            prev_index = qubit_index
            curr_index = self.circuit.get_next_gate(qubit_index, prev_index)
            if not curr_index:
                continue

            while curr_index!=None:
                prev_node = self.circuit.dag[prev_index]
                curr_node = self.circuit.dag[curr_index]
                ### Single-qubit gate Cancellation
                if (prev_node.name, curr_node.name) in self.cancel_rule:
                    opt_flag = True
                    self.circuit.dag.remove_edge(prev_index, curr_index)

                    assert(len(self.circuit.dag.predecessor_indices(prev_index))==1)
                    pprev_index = self.circuit.get_prev_gate(qubit_index, prev_index)
                    self.circuit.dag.remove_edge(pprev_index, prev_index)
                    self.circuit.dag.remove_node(prev_index)
                    next_index = self.circuit.get_next_gate(qubit_index, curr_index)
                    if next_index:
                        self.circuit.dag.remove_edge(curr_index, next_index)
                        self.circuit.dag.add_edge(pprev_index, next_index, None)
                        self.circuit.dag.remove_node(curr_index)
                        prev_index = next_index
                        curr_index = self.circuit.get_next_gate(qubit_index, prev_index)
                        continue
                    else:
                        self.circuit.dag.remove_node(curr_index)
                        break
                ### Two-qubit gate Cancellation - qubit_index : control
                elif prev_node.name=='cx':
                    control_index = prev_node.qubits[0]
                    target_index = prev_node.qubits[1]
                    ### Z rotation can be commutated in control qubit
                    control_next_index = self.circuit.get_next_gate(control_index, prev_index)
                    control_adjacent = True
                    while control_next_index!=None:
                        control_next_node = self.circuit.dag[control_next_index]
                        if control_next_node.name in ['s', 'sdg', 'z', 'rz']:
                            control_adjacent = False
                            control_next_index = self.circuit.get_next_gate(control_index, control_next_index)
                        else:
                            break
                    ### X rotation can be commutated in target qubit
                    target_next_index = self.circuit.get_next_gate(target_index, prev_index)
                    target_adjacent = True
                    while target_next_index!=None:
                        target_next_node = self.circuit.dag[target_next_index]
                        if target_next_node.name=='x':
                            target_adjacent = False
                            target_next_index = self.circuit.get_next_gate(target_index, target_next_index)
                        else:
                            break

                    tqgate_opt_flag = False
                    if control_next_index==target_next_index and control_next_index!=None:
                        cancel_index = control_next_index
                        cancel_node = self.circuit.dag[cancel_index]
                        assert(cancel_node.name=='cx')
                        if prev_node.qubits[0]==cancel_node.qubits[0] and \
                          prev_node.qubits[1]==cancel_node.qubits[1]:
                            tqgate_opt_flag = True

                    if not tqgate_opt_flag:
                        prev_index = curr_index
                        curr_index = self.circuit.get_next_gate(qubit_index, prev_index)
                        continue    # Do not consider merge for CX gates

                    control_pprev_index = self.circuit.get_prev_gate(control_index, prev_index)
                    control_nprev_index = self.circuit.get_next_gate(control_index, prev_index)
                    control_pcurr_index = self.circuit.get_prev_gate(control_index, cancel_index)
                    control_ncurr_index = self.circuit.get_next_gate(control_index, cancel_index)

                    target_pprev_index = self.circuit.get_prev_gate(target_index, prev_index)
                    target_nprev_index = self.circuit.get_next_gate(target_index, prev_index)
                    target_pcurr_index = self.circuit.get_prev_gate(target_index, cancel_index)
                    target_ncurr_index = self.circuit.get_next_gate(target_index, cancel_index)

                    if control_adjacent:
                        assert(control_nprev_index==cancel_index and control_pcurr_index==prev_index)
                        self.circuit.dag.remove_edge(control_pprev_index, prev_index)
                        self.circuit.dag.remove_edge(prev_index, cancel_index)
                        if control_ncurr_index:
                            self.circuit.dag.remove_edge(cancel_index, control_ncurr_index)
                            self.circuit.dag.add_edge(control_pprev_index, control_ncurr_index, None)
                    else:
                        assert(control_nprev_index!=cancel_index and control_pcurr_index!=prev_index)
                        self.circuit.dag.remove_edge(control_pprev_index, prev_index)
                        self.circuit.dag.remove_edge(prev_index, control_nprev_index)
                        self.circuit.dag.remove_edge(control_pcurr_index, cancel_index)
                        self.circuit.dag.add_edge(control_pprev_index, control_nprev_index, None)
                        if control_ncurr_index:
                            self.circuit.dag.remove_edge(cancel_index, control_ncurr_index)
                            self.circuit.dag.add_edge(control_pcurr_index, control_ncurr_index, None)

                    if target_adjacent:
                        assert(target_nprev_index==cancel_index and target_pcurr_index==prev_index)
                        self.circuit.dag.remove_edge(target_pprev_index, prev_index)
                        self.circuit.dag.remove_edge(prev_index, cancel_index)
                        if target_ncurr_index:
                            self.circuit.dag.remove_edge(cancel_index, target_ncurr_index)
                            self.circuit.dag.add_edge(target_pprev_index, target_ncurr_index, None)
                    else:
                        assert(target_nprev_index!=cancel_index and target_pcurr_index!=prev_index)
                        self.circuit.dag.remove_edge(target_pprev_index, prev_index)
                        self.circuit.dag.remove_edge(prev_index, target_nprev_index)
                        self.circuit.dag.remove_edge(target_pcurr_index, cancel_index)
                        self.circuit.dag.add_edge(target_pprev_index, target_nprev_index, None)
                        if target_ncurr_index:
                            self.circuit.dag.remove_edge(cancel_index, target_ncurr_index)
                            self.circuit.dag.add_edge(target_pcurr_index, target_ncurr_index, None)

                    self.circuit.dag.remove_node(prev_index)
                    self.circuit.dag.remove_node(cancel_index)

                    if qubit_index==control_index and control_adjacent:
                        prev_index = control_ncurr_index
                    elif qubit_index==control_index and not control_adjacent:
                        prev_index = control_nprev_index
                        assert(prev_index!=None)
                    elif qubit_index==target_index and target_adjacent:
                        prev_index = target_ncurr_index
                    else:
                        assert(qubit_index==target_index and not target_adjacent)
                        prev_index = target_nprev_index
                        assert(prev_index!=None)

                    if prev_index!=None:
                        curr_index = self.circuit.get_next_gate(qubit_index, prev_index)
                    else:
                        curr_index = None
                    continue

                ### Merge Rule
                merge_opt_flag = False
                for pattern, merged in self.merge_rule:
                    if prev_node.name!=pattern[0] or curr_node.name!=pattern[1]:
                        continue

                    merge_opt_flag = True
                    if merged!='rz' or (prev_node.parameter[0]+curr_node.parameter[0])<math.pi/2:
                        prev_node.name = merged
                        if merged=='rz':
                            prev_node.parameter[0] = prev_node.parameter[0] + curr_node.parameter[0]
                        self.circuit.dag.remove_edge(prev_index, curr_index)
                        next_index = self.circuit.get_next_gate(qubit_index, curr_index)
                        if next_index:
                            self.circuit.dag.remove_edge(curr_index, next_index)
                            self.circuit.dag.add_edge(prev_index, next_index, None)
                        self.circuit.dag.remove_node(curr_index)
                        curr_index = next_index
                    else:
                        assert(prev_node.parameter[0]+curr_node.parameter[0] >= math.pi/2)
                        assert(prev_node.parameter[0]+curr_node.parameter[0] < math.pi)
                        new_param = prev_node.parameter[0] + curr_node.parameter[0] - math.pi/2
                        prev_node.name = 's'
                        prev_node.parameter = []

                        if new_param>0:
                            curr_node.name = 'rz'
                            curr_node.parameter = [new_param]
                            
                            prev_index = curr_index
                            curr_index = self.circuit.get_next_gate(qubit_index, prev_index)
                        else:
                            self.circuit.dag.remove_edge(prev_index, curr_index)
                            next_index = self.circuit.get_next_gate(qubit_index, curr_index)
                            if next_index:
                                self.circuit.dag.remove_edge(curr_index, next_index)
                                self.circuit.dag.add_edge(prev_index, next_index, None)
                            self.circuit.dag.remove_node(curr_index)
                            curr_index = next_index
                    break

                opt_flag |= merge_opt_flag
                if not merge_opt_flag:
                    prev_index = curr_index
                    curr_index = self.circuit.get_next_gate(qubit_index, prev_index)

        return opt_flag

    def initialize_merging_ruleset(self):
        self.merge_rule.append((('rz', 'rz'), 'rz'))
        self.merge_rule.append((('s', 's'), 'z'))
        self.merge_rule.append((('sdg', 'sdg'), 'z'))

    def initialize_cancelling_ruleset(self):
        self.cancel_rule.append(('x', 'x'))
        self.cancel_rule.append(('y', 'y'))
        self.cancel_rule.append(('z', 'z'))
        self.cancel_rule.append(('h', 'h'))
        self.cancel_rule.append(('s', 'sdg'))
        self.cancel_rule.append(('sdg', 's'))


class NativeConverter:
    def __init__(self):
        self.conversion_rule = {}
        self.initialize_conversion_rule()

    def initialize_conversion_rule(self):
        self.conversion_rule['x'] = [('gpi', 0)]
        self.conversion_rule['y'] = [('gpi', math.pi/2)]
        self.conversion_rule['z'] = [('vz', math.pi)]
        self.conversion_rule['s'] = [('vz', math.pi/2)]
        self.conversion_rule['sdg'] = [('vz', 3*math.pi/2)]
        self.conversion_rule['h'] = [('gpi2', 3*math.pi/2), ('vz', math.pi)]

    def convert_to_native(self, circuit: Circuit) -> Circuit:
        """ Make a new circuit consisting of native gates, from the given (abstract) circuit """
        native_circuit = Circuit(circuit.num_qubits)

        pending_gate = []
        for qubit_index in range(circuit.num_qubits):
            pending_gate.append(circuit.get_next_gate(qubit_index, qubit_index))
        stack = []

        for q in range(circuit.num_qubits):
            if pending_gate[q]!=None:
                stack.append((q, pending_gate[q]))
            else:
                continue

            while stack:
                qubit_index, curr_node_index = stack.pop()
                assert(pending_gate[qubit_index]==curr_node_index)
                curr_node = circuit.dag[curr_node_index]
                if curr_node.name=='cx':
                    if stack and stack[-1][1]==curr_node_index:
                        pair_qubit_index, pair_node_index = stack.pop()
                        assert(pair_qubit_index==curr_node.qubits[(curr_node.qubits.index(qubit_index)+1)%2])
                        control_qubit_index = curr_node.qubits[0]
                        target_qubit_index = curr_node.qubits[1]
                        native_circuit.append_gate('gpi2', [control_qubit_index], [math.pi/2])
                        native_circuit.append_gate('ms', [control_qubit_index, target_qubit_index], [0, 0, math.pi/2], duration=5)
                        native_circuit.append_gate('gpi2', [control_qubit_index], [math.pi])
                        native_circuit.append_gate('gpi2', [target_qubit_index], [math.pi])
                        native_circuit.append_gate('gpi2', [control_qubit_index], [3*math.pi/2])

                        next_node_index = circuit.get_next_gate(qubit_index, curr_node_index)
                        pending_gate[qubit_index] = next_node_index

                        pair_next_index = circuit.get_next_gate(pair_qubit_index, pair_node_index)
                        pending_gate[pair_qubit_index] = pair_next_index
                        if pair_next_index!=None:
                            stack.append((pair_qubit_index, pair_next_index))
                    else:
                        pair_qubit_index = curr_node.qubits[(curr_node.qubits.index(qubit_index)+1)%2]
                        stack.append((qubit_index, curr_node_index))
                        stack.append((pair_qubit_index, pending_gate[pair_qubit_index]))
                    continue
                elif curr_node.name=='rz':
                    param = curr_node.parameter[0]
                    assert(param>0 and param<math.pi/2)
                    native_circuit.append_gate('rz', [qubit_index], [param])
                else:
                    assert(curr_node.name in self.conversion_rule.keys())
                    for native_gate_name, param in self.conversion_rule[curr_node.name]:
                        native_circuit.append_gate(native_gate_name, [qubit_index], [param])

                next_index = circuit.get_next_gate(qubit_index, curr_node_index)
                pending_gate[qubit_index] = next_index
                if next_index!=None:
                    stack.append((qubit_index, next_index))

        return native_circuit

class NativeCircuitOptimizer:
    def __init__(self, circuit: Circuit):
        self.circuit = circuit

    def native_opt(self, qec: bool):
        opt_flag = True
        opt_cnt = 0
        while opt_flag:
            opt_flag = False
            opt_flag |= self.vz_elimination(qec)
            opt_flag |= self.gpi_cancellation(qec)
            opt_cnt += 1

            assert((opt_flag and opt_cnt==1) or (not opt_flag and opt_cnt==2))

        return self.circuit

    def gpi_cancellation(self, qec: bool) -> bool:
        opt_flag = False

        for qubit_index in range(self.circuit.num_qubits):
            accumulated_phase = 0.0
            prev_index = qubit_index
            curr_index = self.circuit.get_next_gate(qubit_index, prev_index)
            if curr_index!=None:
                next_index = self.circuit.get_next_gate(qubit_index, curr_index)
            else:
                next_index = None

            while curr_index!=None:
                curr_node = self.circuit.dag[curr_index]
                if next_index!=None:
                    next_node = self.circuit.dag[next_index]
                else:
                    next_node = None

                if curr_node.name=='gpi' and next_node!=None and next_node.name=='gpi':
                    opt_flag = True
                    new_param = 2*(next_node.parameter[0] - curr_node.parameter[0])
                    accumulated_phase += new_param
                    while accumulated_phase<0.0:
                        accumulated_phase += (2*math.pi)
                    while accumulated_phase>=2*math.pi:
                        accumulated_phase -= (2*math.pi)

                    if qec:
                        assert(accumulated_phase in [n*math.pi for n in range(0, 2)])

                    self.circuit.dag.remove_edge(prev_index, curr_index)
                    self.circuit.dag.remove_edge(curr_index, next_index)
                    nnext_index = self.circuit.get_next_gate(qubit_index, next_index)
                    if nnext_index!=None:
                        self.circuit.dag.remove_edge(next_index, nnext_index)
                        self.circuit.dag.add_edge(prev_index, nnext_index, None)
                    self.circuit.dag.remove_node(curr_index)
                    self.circuit.dag.remove_node(next_index)

                    curr_index = nnext_index
                    if curr_index!=None:
                        next_index = self.circuit.get_next_gate(qubit_index, curr_index)
                    else:
                        next_index = None

                    continue
                elif (curr_node.name=='gpi' or curr_node.name=='gpi2'):
                    new_param = curr_node.parameter[0] - accumulated_phase
                    while new_param<0:
                        new_param += (2*math.pi)
                    while new_param>=2*math.pi:
                        new_param -= (2*math.pi)
                    if qec:
                        assert(new_param in [n*math.pi/2 for n in range(0, 4)])

                    curr_node.parameter[0] = new_param
                elif curr_node.name=='ms':
                    param_index = curr_node.qubits.index(qubit_index)
                    assert(param_index==0 or param_index==1)
                    curr_node.parameter[param_index] -= accumulated_phase
                elif curr_node.name=='rz':
                    ### Do nothing, just pass by commutation
                    ### For no qec circuit, rz can appear only at the end of rz
                    if qec and next_index!=None:
                        assert(curr_node.parameter[0]>0 and curr_node.parameter[0]<math.pi/2)
                    elif qec and next_index==None:
                        assert(self.circuit.get_next_gate(qubit_index, curr_index)==None)
                    else:
                        assert(not qec and next_index==None and curr_node.parameter[0]>0 and curr_node.parameter[0]<2*math.pi)
                else:
                    ### vz cannot appear in gpi cancellation, if vz_elimination is previously executed
                    raise Exception(f'Invalid gate %s'%curr_node.name)

                prev_index = curr_index
                curr_index = next_index
                if curr_index!=None:
                    next_index = self.circuit.get_next_gate(qubit_index, curr_index)

            assert(accumulated_phase>=0.0 and accumulated_phase<2*math.pi)
            if qec:
                assert(accumulated_phase in [(math.pi/2)*n for n in range(0, 4)])

            last_node = self.circuit.dag[prev_index]
            assert(self.circuit.get_next_gate(qubit_index, prev_index)==None)
            if accumulated_phase!=0.0 and last_node.name=='rz':
                new_param = last_node.parameter[0] + accumulated_phase
                while new_param<0:
                    new_param += (2*math.pi)
                while new_param>=(2*math.pi):
                    new_param -= (2*math.pi)
                if new_param!=0.0:
                    last_node.parameter[0] = new_param
                else:
                    pprev_index = self.circuit.get_prev_gate(qubit_index, prev_index)
                    self.circuit.dag.remove_edge(pprev_index, prev_index)
                    self.circuit.dag.remove_node(prev_index)
            elif accumulated_phase!=0.0:
                new_node = Node('rz', [qubit_index], [accumulated_phase])
                new_index = self.circuit.dag.add_node(new_node)
                assert(self.circuit.get_next_gate(qubit_index, prev_index)==None)
                self.circuit.dag.add_edge(prev_index, new_index, None)

        return opt_flag

    def vz_elimination(self, qec: bool) -> bool:
        opt_flag = False

        for qubit_index in range(self.circuit.num_qubits):
            accumulated_phase = 0.0
            prev_index = qubit_index
            curr_index = self.circuit.get_next_gate(qubit_index, prev_index)
            if curr_index!=None:
                next_index = self.circuit.get_next_gate(qubit_index, curr_index)
            else:
                next_index = None

            while curr_index!=None:
                curr_node = self.circuit.dag[curr_index]
                if curr_node.name=='gpi' or curr_node.name=='gpi2':
                    new_param = curr_node.parameter[0] - accumulated_phase
                    while new_param<0:
                        new_param += (2*math.pi)
                    while new_param>=2*math.pi:
                        new_param -= (2*math.pi)
                    if qec:
                        assert(new_param in [(math.pi/2)*n for n in range(0, 4)])

                    curr_node.parameter[0] = new_param
                    prev_index = curr_index
                elif curr_node.name=='vz':
                    opt_flag = True
                    assert(curr_node.parameter[0] in [(math.pi/2)*n for n in range(1, 4)])
                    ### VZ Gate is deleted regardless of qec or not
                    accumulated_phase += curr_node.parameter[0]
                    while accumulated_phase<0:
                        accumulated_phase += (2*math.pi)
                    while accumulated_phase>=2*math.pi:
                        accumulated_phase -= (2*math.pi)
                    if qec:
                        assert(accumulated_phase in [(math.pi/2)*n for n in range(0, 4)])

                    self.circuit.dag.remove_edge(prev_index, curr_index)
                    if next_index!=None:
                        self.circuit.dag.remove_edge(curr_index, next_index)
                        self.circuit.dag.add_edge(prev_index, next_index, None)
                    self.circuit.dag.remove_node(curr_index)
                elif curr_node.name=='ms':
                    if qec:
                        assert(accumulated_phase in [n*math.pi/2 for n in range(0, 4)])
                    param_index = curr_node.qubits.index(qubit_index)
                    assert(param_index==0 or param_index==1)
                    curr_node.parameter[param_index] -= accumulated_phase
                    prev_index = curr_index
                elif curr_node.name=='rz':
                    assert(curr_node.parameter[0]>0 and curr_node.parameter[0]<2*math.pi)

                    if not qec and next_index!=None:
                        opt_flag = True
                        ### Also treated as vz
                        accumulated_phase += curr_node.parameter[0]
                        while accumulated_phase<0:
                            accumulated_phase += (2*math.pi)
                        while accumulated_phase>=2*math.pi:
                            accumulated_phase -= (2*math.pi)

                        self.circuit.dag.remove_edge(prev_index, curr_index)
                        if next_index!=None:
                            self.circuit.dag.remove_edge(curr_index, next_index)
                            self.circuit.dag.add_edge(prev_index, next_index, None)
                        self.circuit.dag.remove_node(curr_index)
                    else:
                        ### Moving VZ commutes with current rz, passing without any modification
                        prev_index = curr_index
                else:
                    raise Exception(f"Invalid gate %s"%curr_node.name)

                ### prev_index should be updated in the if-else statements above
                curr_index = next_index
                if curr_index!=None:
                    next_index = self.circuit.get_next_gate(qubit_index, curr_index)

            assert(accumulated_phase>=0.0 and accumulated_phase<2*math.pi)
            if qec:
                assert(accumulated_phase in [(math.pi/2)*n for n in range(0, 4)])
            
            last_node = self.circuit.dag[prev_index]
            if accumulated_phase!=0.0 and last_node.name=='rz':
                new_param = last_node.parameter[0] + accumulated_phase
                while new_param<0:
                    new_param += (2*math.pi)
                while new_param>=(2*math.pi):
                    new_param -= (2*math.pi)
                if new_param!=0.0:
                    last_node.parameter[0] = new_param
                else:
                    pprev_index = self.circuit.get_prev_gate(qubit_index, prev_index)
                    self.circuit.dag.remove_edge(pprev_index, prev_index)
                    self.circuit.dag.remove_node(prev_index)
            elif accumulated_phase!=0.0:
                new_node = Node('rz', [qubit_index], [accumulated_phase])
                new_index = self.circuit.dag.add_node(new_node)
                assert(self.circuit.get_next_gate(qubit_index, prev_index)==None)
                self.circuit.dag.add_edge(prev_index, new_index, None)

        return opt_flag

def rz_approximation(native_circ: Circuit, threshold: int) -> (Circuit, int):
    gate_cnt = 0
    pending_gate = []
    for qubit_index in range(native_circ.num_qubits):
        pending_gate.append(native_circ.get_next_gate(qubit_index, qubit_index))

    qubit_stack = []
    for qubit_index in range(native_circ.num_qubits):
        curr_qubit_index = qubit_index
        curr_index = pending_gate[curr_qubit_index]

        while curr_index!=None:
            curr_node = native_circ.dag[curr_index]
            next_index = native_circ.get_next_gate(curr_qubit_index, curr_index)
            if curr_node.name in ['gpi', 'gpi2']:
                gate_cnt += 1
                pending_gate[curr_qubit_index] = next_index
                curr_index = next_index
            elif curr_node.name=='rz' and next_index!=None:
                assert(isinstance(curr_node.parameter[0], float))
                assert(curr_node.parameter[0]>0 and curr_node.parameter[0]<(math.pi/2))

                prev_index = native_circ.get_prev_gate(curr_qubit_index, curr_index)                
                native_circ.dag.remove_edge(prev_index, curr_index)
                native_circ.dag.remove_edge(curr_index, next_index)

                ### Discretize angle of Rz
                gate_applied = [False for _ in range(threshold+1)]
                applied_cnt = 0
                for n in range(2, threshold+1):
                    divisor = 2**n
                    if curr_node.parameter[0]>=(math.pi/divisor):
                        curr_node.parameter[0] -= (math.pi/divisor)
                        gate_applied[n] = True
                        applied_cnt += 1

                if applied_cnt>int((threshold-1)/2):
                    for n in range(2, threshold+1):
                        gate_applied[n] = not gate_applied[n]

                prev_rz_index = prev_index
                for n in range(2, threshold+1):
                    if gate_applied[n]:
                        new_rz_index = native_circ.dag.add_node(Node('rz', [curr_qubit_index], [n-1], clifford=False))
                        native_circ.dag.add_edge(prev_rz_index, new_rz_index, None)
                        prev_rz_index = new_rz_index
                        gate_cnt += 1
                native_circ.dag.add_edge(prev_rz_index, next_index, None)
                native_circ.dag.remove_node(curr_index)

                pending_gate[curr_qubit_index] = next_index
                curr_index = next_index
            elif curr_node.name=='rz':
                assert(next_index==None)
                assert(isinstance(curr_node.parameter[0], float))
                assert(curr_node.parameter[0]>0 and curr_node.parameter[0]<(2*math.pi))

                prev_index = native_circ.get_prev_gate(curr_qubit_index, curr_index)
                native_circ.dag.remove_edge(prev_index, curr_index)

                if curr_node.parameter[0]>=math.pi:
                    curr_node.parameter[0] -= math.pi
                if curr_node.parameter[0]>=(math.pi/2):
                    curr_node.parameter[0] -= (math.pi/2)

                ### Discretize angle of Rz
                gate_applied = [False for _ in range(threshold+1)]
                applied_cnt = 0
                for n in range(2, threshold+1):
                    divisor = 2**n
                    if curr_node.parameter[0]>=(math.pi/divisor):
                        curr_node.parameter[0] -= (math.pi/divisor)
                        gate_applied[n] = True
                        applied_cnt += 1

                if applied_cnt>int((threshold-1)/2):
                    for n in range(2, threshold+1):
                        gate_applied[n] = not gate_applied[n]

                prev_rz_index = prev_index
                for n in range(2, threshold+1):
                    if gate_applied[n]:
                        new_rz_index = native_circ.dag.add_node(Node('rz', [curr_qubit_index], [n-1], clifford=False))
                        native_circ.dag.add_edge(prev_rz_index, new_rz_index, None)
                        prev_rz_index = new_rz_index
                        gate_cnt += 1
                native_circ.dag.remove_node(curr_index)

                pending_gate[curr_qubit_index] = next_index
                curr_index = next_index
            elif curr_node.name=='ms':
                pair_qubit_index = curr_node.qubits[(curr_node.qubits.index(curr_qubit_index)+1)%2]
                if qubit_stack and qubit_stack[-1]==(pair_qubit_index, curr_index):
                    qubit_stack.pop()
                    gate_cnt += 1
                    pending_gate[curr_qubit_index] = native_circ.get_next_gate(curr_qubit_index, curr_index)
                    pending_gate[pair_qubit_index] = native_circ.get_next_gate(pair_qubit_index, curr_index)
                else:
                    qubit_stack.append((curr_qubit_index, curr_index))
                curr_qubit_index = pair_qubit_index
                curr_index = pending_gate[pair_qubit_index]
            else:
                raise Exception(f'Invalid gate instruction %s %s'%(curr_node.name, curr_node.qubits))
            
    for qubit_index in range(native_circ.num_qubits):
        assert(pending_gate[qubit_index]==None)

    return native_circ, gate_cnt