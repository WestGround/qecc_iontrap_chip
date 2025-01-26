import os, math, datetime
from copy import deepcopy

try:
    from qiskit import QuantumCircuit, transpile
except:
    print("Qiskit version should be lower than 1.0")
try:
    from qiskit.circuit.library import UnitaryGate
except:
    print("Qiskit version should be over 1.0")

from qiskit.quantum_info.operators import Operator
from qiskit.quantum_info import Statevector

from benchmark import *
from Circuit import Circuit
from CircuitOpt import AbsCircuitOptimizer, NativeConverter, NativeCircuitOptimizer, rz_approximation

_filepath = os.path.abspath(__file__)
_dirname = os.path.dirname(_filepath)

def main(algorithm: str, N: int, approx_factor: int):
    ### 1. Generate Qiskit circuit from benchmark
    print("%s (%s %d(af %d)) Start generating Qiskit"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    qiskit_circ: QuantumCircuit
    if algorithm.lower()=="ae":
        qiskit_circ = amplitude_estimation(N, approx_factor)
    elif algorithm.lower()=="bv":
        s_int = 1
        for _ in range(0, N-2, 2):
            s_int = 4*s_int + 1
        qiskit_circ = bernstein_vazirani(N, s_int)
    elif algorithm.lower()=="grover":
        s_int = 1
        for _ in range(0, N-2, 2):
            s_int = 4*s_int + 1
        num_iter = int(math.pi*math.sqrt(math.pow(2, N))/4)
        qiskit_circ = grover_search(N, s_int, num_iter, approx_factor)
    elif algorithm.lower()=="hs":
        qiskit_circ = hamiltonian(N)
    elif algorithm.lower()=="mc":
        mu = 0.5
        qiskit_circ = monte_carlo(N, approx_factor, mu)
    elif algorithm.lower()=="pe":
        theta = 0.5
        qiskit_circ = phase_estimation(N, theta)
    elif algorithm.lower()=="qft":
        s_int = 1
        for _ in range(0, N-2, 2):
            s_int = 4*s_int + 1
        qiskit_circ = qft(N, s_int)
    elif algorithm.lower()=="vqe":
        qiskit_circ = vqe(N)
    print("%s (%s %d(af %d)) Finish generating Qiskit"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))

    print("%s (%s %d(af %d)) Start transpiling"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    qiskit_circ = transpile(qiskit_circ, basis_gates=['cx', 'x', 'y', 'z', 'rz','h'])
    print("%s (%s %d(af %d)) Finish transpiling"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))

    ### 2. Convert Qiskit to Abstract DAG
    print("%s (%s %d(af %d)) Start Qiskit to DAG"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    abs_circ = qiskit_to_circuit(qiskit_circ)
    print("%s (%s %d(af %d) Finish Qiskit to DAG"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))

    ### 3. Abstract Optimization
    print("%s (%s %d(af %d)) Start abstract optimization"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    abs_opter = AbsCircuitOptimizer(abs_circ)
    abs_circ = abs_opter.abs_opt()
    print("%s (%s %d(af %d)) Finish abstract optimization"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))

    ### 4. Convert to Native
    print("%s (%s %d(af %d)) Start native conversion"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    native_circ = NativeConverter().convert_to_native(abs_circ)
    print("%s (%s %d(af %d)) Finish native conversion"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))

    ### 5. Native Optimization
    print("%s (%s %d(af %d)) Start native optimization with QEC"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    native_circ_qec = deepcopy(native_circ)
    native_opter_qec = NativeCircuitOptimizer(native_circ_qec)
    native_circ_qec = native_opter_qec.native_opt(qec=True)
    native_circ_qec, gate_cnt_qec1 = rz_approximation(native_circ_qec, threshold=13)
    print("%s (%s %d(af %d)) Finish native optimization with QEC"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))

    print("%s (%s %d(af %d)) Start native optimization without QEC"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    native_circ_noqec = deepcopy(native_circ)
    native_opter_noqec = NativeCircuitOptimizer(native_circ_noqec)
    native_circ_noqec = native_opter_noqec.native_opt(qec=False)
    print("%s (%s %d(af %d)) Finish native optimization without QEC"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))

    ### 5-1. Intermediate step for state comparison
    ### Only works for circuits with small number of qubits
    ### To pass the following test, rz_approximation in line 87 should be deleted.
    ### Deleting line 87 disables circuit_to_txt() in line 118, so use it for validity check only.
    """print("%s (%s %d) Start reconstructing Qiskit"%(datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N))
    recon_qiskit_qec = native_qiskit_reconstructor(native_circ_qec)
    recon_qiskit_noqec = native_qiskit_reconstructor(native_circ_noqec)
    print("%s (%s %d) Finish reconstructing Qiskit"%(datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N))

    print("%s (%s %d) Start comparison"%(datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N))
    orig_state = Statevector.from_instruction(qiskit_circ)
    recon_state_qec = Statevector.from_instruction(recon_qiskit_qec)
    recon_state_noqec = Statevector.from_instruction(recon_qiskit_noqec)
    assert(orig_state.equiv(recon_state_qec))
    assert(orig_state.equiv(recon_state_noqec))
    print("%s (%s %d) Finish comparison"%(datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N))"""

    ### 6. Write to file
    print("%s (%s %d(af %d)) Start writing to file"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    gate_cnt_qec2 = native_circ_qec.circuit_to_txt(\
        os.path.join(_dirname, "circuit", "%s(%d)_af(%d)_qec.txt"%(algorithm, N, approx_factor)), True, approx_factor)
    gate_cnt_noqec = native_circ_noqec.circuit_to_txt(\
        os.path.join(_dirname, "circuit", "%s(%d)_af(%d)_noqec.txt"%(algorithm, N, approx_factor)), False, approx_factor)

    assert(gate_cnt_qec1==gate_cnt_qec2)
    print("%s (%s %d(af %d)) Finish writing to file"%\
          (datetime.datetime.today().strftime("[%H:%M:%S]"), algorithm, N, approx_factor))
    
def native_qiskit_reconstructor(native_circ: Circuit) -> QuantumCircuit:
    new_qiskit_qc = QuantumCircuit(native_circ.num_qubits)
    pending_gate = []
    for i in range(native_circ.num_qubits):
        pending_gate.append(native_circ.get_next_gate(i, i))
    stack = []

    for q in range(native_circ.num_qubits):
        if pending_gate[q]!=None:
            stack.append((q, pending_gate[q]))
        else:
            continue

        while stack:
            qubit_index, curr_index = stack.pop()
            assert(pending_gate[qubit_index]==curr_index)
            curr_node = native_circ.dag[curr_index]

            if curr_node.name=='ms':
                if stack and stack[-1][1]==curr_index:
                    pair_index, curr_index = stack.pop()
                    assert(pair_index==curr_node.qubits[(curr_node.qubits.index(qubit_index)+1)%2])
                    control_index = curr_node.qubits[0]
                    target_index = curr_node.qubits[1]

                    phi0 = curr_node.parameter[0]
                    phi1 = curr_node.parameter[1]
                    T = curr_node.parameter[2]
                    assert(T==math.pi/2)

                    S = phi0 + phi1
                    D = phi0 - phi1
                    T /= 2

                    matrix = [[math.cos(T), 0, 0, complex((-1)*math.sin(S)*math.sin(T), (-1)*(math.sin(T)*math.cos(S)))],
                      [0, math.cos(T), complex(math.sin(T)*math.sin(D), (-1)*math.sin(T)*math.cos(D)), 0],
                      [0, complex((-1)*math.sin(T)*math.sin(D), (-1)*math.sin(T)*math.cos(D)), math.cos(T), 0],
                      [complex(math.sin(S)*math.sin(T), (-1)*(math.sin(T)*math.cos(S))), 0, 0, math.cos(T)]]
                    #ms_gate = UnitaryGate(matrix)
                    ms_gate = Operator(matrix)
                    new_qiskit_qc.append(ms_gate, [control_index, target_index])

                    next_index = native_circ.get_next_gate(qubit_index, curr_index)
                    pending_gate[qubit_index] = next_index

                    pair_next_index = native_circ.get_next_gate(pair_index, curr_index)
                    pending_gate[pair_index] = pair_next_index
                    if pair_next_index!=None:
                        stack.append((pair_index, pair_next_index))
                else:
                    pair_index = curr_node.qubits[(curr_node.qubits.index(qubit_index)+1)%2]
                    stack.append((qubit_index, curr_index))
                    stack.append((pair_index, pending_gate[pair_index]))
                continue

            if curr_node.name=='gpi':
                param = curr_node.parameter[0]
                new_qiskit_qc.u(math.pi, param, math.pi-param, qubit_index)
            elif curr_node.name=='gpi2':
                param = curr_node.parameter[0]
                new_qiskit_qc.u(math.pi/2, param-math.pi/2, math.pi/2-param, qubit_index)
            elif curr_node.name=='vz' or curr_node.name=='rz':
                param = curr_node.parameter[0]
                new_qiskit_qc.rz(param, qubit_index)
            else:
                raise Exception(f'Invalid gate type: %s'%curr_node.name)
            
            next_index = native_circ.get_next_gate(qubit_index, curr_index)
            pending_gate[qubit_index] = next_index
            if next_index!=None:
                stack.append((qubit_index, next_index))

    return new_qiskit_qc

def qiskit_to_circuit(qiskit_qc: QuantumCircuit) -> Circuit:
    qubit_offset = {}
    for qubit in qiskit_qc.qubits:
        qubit_offset[qubit] = len(qubit_offset)
    circuit = Circuit(len(qubit_offset))

    for entry in qiskit_qc.data:
        name = entry[0].name
        qubit = entry[1][0]
        offset = qubit_offset[qubit]

        if name.lower() in ['x', 'y', 'z', 'h']:
            assert(len(entry[1])==1)
            circuit.append_gate(name.lower(), [offset], [])
        elif name.lower()=='rz':
            assert(len(entry[1])==1)
            assert(len(entry[0].params)==1)
            param = float(entry[0].params[0])
            param = param - int(param/(2*math.pi))*(2*math.pi)
            while param<0.0:
                param += (2*math.pi)
            while param>=2*math.pi:
                param -= (2*math.pi)
            if param>=3*math.pi/2:
                circuit.append_gate('sdg', [offset], [])
                param -= 3*math.pi/2
            elif param>=math.pi:
                circuit.append_gate('z', [offset], [])
                param -= math.pi
            elif param>=math.pi/2:
                circuit.append_gate('s', [offset], [])
                param -= math.pi/2
            assert(param>=0.0 and param<math.pi)
            if param!=0.0:
                circuit.append_gate('rz', [offset], [param])
        elif name.lower() in ['cx' or 'cnot']:
            assert(len(entry[1])==2)
            qubit2 = entry[1][1]
            offset2 = qubit_offset[qubit2]
            circuit.append_gate('cx', [offset, offset2], [])
        elif name=="measure" or name=="barrier":
            continue
        else:
            print("invalid gate : %s q%d"%(name, offset))
            exit(-1)

    return circuit

if __name__=="__main__":    
    native_set = {
                  "vqe": [(8, 1)],
                  "qft": [(8, 1)],
                  "ae": [(8, 1)],
                  "bv": [(8, 1)],
                  "pe": [(8, 1)],
                  "hs": [(8, 1)],
                  "grover": [(8, 1)],
                  "mc": [(8, 1)],
                }

    for algorithm, num_qubit_range in native_set.items():
        for N, af in num_qubit_range:
            main(algorithm, N, af)