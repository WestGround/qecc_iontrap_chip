"""
Variational Quantum Eigensolver Benchmark Program - Qiskit
"""

import json
import os
import sys
import time

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, transpile
try:
    from qiskit.opflow import PauliTrotterEvolution, Suzuki
    from qiskit.opflow.primitive_ops import PauliSumOp
except:
    print("Need Qiskit version less than 1.0")

_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path[1:1] = [_dirname]

# Benchmark Name
benchmark_name = "VQE Simulation"

verbose= False

# saved circuits for display
QC_ = None
Hf_ = None
CO_ = None

################### Circuit Definition #######################################

# Construct a Qiskit circuit for VQE Energy evaluation with UCCSD ansatz
# param: n_spin_orbs - The number of spin orbitals.
# return: return a Qiskit circuit for this VQE ansatz
def vqe(n_spin_orbs, circuit_id=0, method=1):

    # number of alpha spin orbitals
    norb_a = int(n_spin_orbs / 2)

    # construct the Hamiltonian
    qubit_op = ReadHamiltonian(n_spin_orbs)

    # allocate qubits
    num_qubits = n_spin_orbs
    na = int(num_qubits/4)
    nb = int(num_qubits/4)

    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr, name=f"vqe-ansatz({method})-{num_qubits}-{circuit_id}")

    # initialize the HF state
    Hf = HartreeFock(num_qubits, na, nb)
    qc.append(Hf, qr)

    # form the list of single and double excitations 
    excitationList = []
    for occ_a in range(na):
        for vir_a in range(na, norb_a):
            excitationList.append((occ_a, vir_a))

    for occ_b in range(norb_a, norb_a+nb):
        for vir_b in range(norb_a+nb, n_spin_orbs):
            excitationList.append((occ_b, vir_b))

    for occ_a in range(na):
        for vir_a in range(na, norb_a):
            for occ_b in range(norb_a, norb_a+nb):
                for vir_b in range(norb_a+nb, n_spin_orbs):
                    excitationList.append((occ_a, vir_a, occ_b, vir_b))

    # get cluster operators in Paulis
    pauli_list = readPauliExcitation(n_spin_orbs, circuit_id)

    # loop over the Pauli operators
    for index, PauliOp in enumerate(pauli_list):
        # get circuit for exp(-iP)
        cluster_qc = ClusterOperatorCircuit(PauliOp, excitationList[index])

        # add to ansatz
        qc.append(cluster_qc, [i for i in range(cluster_qc.num_qubits)])

    # method 1, only compute the last term in the Hamiltonian
    if method == 1:
        # last term in Hamiltonian
        qc_with_mea, is_diag = ExpectationCircuit(qc, qubit_op[1], num_qubits)

        # return the circuit
        return qc_with_mea

    # now we need to add the measurement parts to the circuit
    # circuit list 
    qc_list = []
    diag = []
    off_diag = []
    global normalization
    normalization = 0.0

    # add the first non-identity term
    identity_qc = qc.copy()
    identity_qc.measure_all()
    qc_list.append(identity_qc) # add to circuit list
    diag.append(qubit_op[1])
    normalization += abs(qubit_op[1].coeffs[0]) # add to normalization factor
    diag_coeff = abs(qubit_op[1].coeffs[0]) # add to coefficients of diagonal terms

    # loop over rest of terms 
    for index, p in enumerate(qubit_op[2:]):
        
        # get the circuit with expectation measurements
        qc_with_mea, is_diag = ExpectationCircuit(qc, p, num_qubits)

        # accumulate normalization 
        normalization += abs(p.coeffs[0])

        # add to circuit list if non-diagonal
        if not is_diag:
            qc_list.append(qc_with_mea)
        else:
            diag_coeff += abs(p.coeffs[0])

        # diagonal term
        if is_diag:
            diag.append(p)
        # off-diagonal term
        else:
            off_diag.append(p)

    # modify the name of diagonal circuit
    qc_list[0].name = qubit_op[1].primitive.to_list()[0][0] + " " + str(np.real(diag_coeff))
    normalization /= len(qc_list)
    return qc_list

# Function that constructs the circuit for a given cluster operator
def ClusterOperatorCircuit(pauli_op, excitationIndex):
    
    # compute exp(-iP)
    exp_ip = pauli_op.exp_i()

    # Trotter approximation
    qc_op = PauliTrotterEvolution(trotter_mode=Suzuki(order=1, reps=1)).convert(exp_ip)

    # convert to circuit
    qc = qc_op.to_circuit(); qc.name = f'Cluster Op {excitationIndex}'
    global CO_
    if CO_ == None or qc.num_qubits <= 4:
        if qc.num_qubits < 7: CO_ = qc

    # return this circuit
    #qc.draw(output='mpl', filename='clusterop')
    return qc


# Function that adds expectation measurements to the raw circuits
def ExpectationCircuit(qc, pauli, nqubit, method=2):

    # copy the unrotated circuit
    raw_qc = qc.copy()

    # whether this term is diagonal
    is_diag = True

    # primitive Pauli string
    PauliString = pauli.primitive.to_list()[0][0]

    # coefficient
    coeff = pauli.coeffs[0]

    # basis rotation
    for i, p in enumerate(PauliString):
    
        target_qubit = nqubit - i - 1 
        if (p == "X"):
            is_diag = False
            raw_qc.h(target_qubit)
        elif (p == "Y"):
            raw_qc.sdg(target_qubit)
            raw_qc.h(target_qubit)
            is_diag = False

    # perform measurements - do not measure
    # raw_qc.measure_all()

    # name of this circuit
    raw_qc.name = PauliString + " " + str(np.real(coeff))

    # save circuit
    global QC_
    if QC_ == None or nqubit <= 4:
        if nqubit < 7: QC_ = raw_qc

    #raw_qc.draw(output='mpl', filename="exp_circ")
    return raw_qc, is_diag

# Function that implements the Hartree-Fock state 
def HartreeFock(norb, na, nb):

    # initialize the quantum circuit
    qc = QuantumCircuit(norb, name="Hf")
    
    # alpha electrons
    for ia in range(na):
        qc.x(ia)

    # beta electrons
    for ib in range(nb):
        qc.x(ib+int(norb/2))

    # Save smaller circuit
    global Hf_
    if Hf_ == None or norb <= 4:
        if norb < 7: Hf_ = qc

    # return the circuit
    #qc.draw(output='mpl', filename='hf')
    return qc

################ Helper Functions

# Function that converts a list of single and double excitation operators to Pauli operators
def readPauliExcitation(norb, circuit_id=0):

    # load pre-computed data
    filename = os.path.join(os.path.dirname(__file__), f'vqe/ansatz_{norb}.txt')
    with open(filename) as f:
        data = f.read()
    ansatz_dict = json.loads(data)

    # initialize Pauli list
    pauli_list = []

    # current coefficients 
    cur_coeff = 1e5

    # current Pauli list 
    cur_list = []

    # loop over excitations
    for ext in ansatz_dict:

        if cur_coeff > 1e4:
            cur_coeff = ansatz_dict[ext]
            cur_list = [(ext, ansatz_dict[ext])]
        elif abs(abs(ansatz_dict[ext]) - abs(cur_coeff)) > 1e-4:
            pauli_list.append(PauliSumOp.from_list(cur_list))
            cur_coeff = ansatz_dict[ext]
            cur_list = [(ext, ansatz_dict[ext])]
        else:
            cur_list.append((ext, ansatz_dict[ext]))
        
    # add the last term
    pauli_list.append(PauliSumOp.from_list(cur_list))

    # return Pauli list
    return pauli_list

# Get the Hamiltonian by reading in pre-computed file
def ReadHamiltonian(nqubit):

    # load pre-computed data
    filename = os.path.join(os.path.dirname(__file__), f'vqe/hamiltonian_{nqubit}.txt')
    with open(filename) as f:
        data = f.read()
    ham_dict = json.loads(data)

    # pauli list 
    pauli_list = []
    for p in ham_dict:
        pauli_list.append( (p, ham_dict[p]) )

    # build Hamiltonian
    ham = PauliSumOp.from_list(pauli_list)

    # return Hamiltonian
    return ham

# if main, execute methods     
if __name__ == "__main__":
    n = 16
    qc = transpile(vqe(n), basis_gates=['x', 'y', 'z', 'h', 'rz', 'cx'])
    qc.draw(output='mpl', filename="../vqe_%d"%n)