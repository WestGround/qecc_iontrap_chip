"""
Amplitude Estimation Benchmark Program via Phase Estimation - Qiskit
"""

import copy
import sys
import time
import os

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path[1:1] = [_dirname]
from QFT import inv_qft_gate
from MonteCarlo import mcx

# Benchmark Name
benchmark_name = "Amplitude Estimation"

np.random.seed(0)

verbose = False

# saved subcircuits circuits for printing
A_ = None
Q_ = None
cQ_ = None
QC_ = None
QFTI_ = None

############### Circuit Definition

def amplitude_estimation(num_qubits: int, approx_factor: int, \
                         s_int=None, psi_zero=None, psi_one=None):
    num_state_qubits = 1
    num_counting_qubits = num_qubits - num_state_qubits - 1

    if s_int==None:
        s_int = 2**(num_counting_qubits-1)
    theta = s_int * np.pi / (2**num_counting_qubits)
    precision = int(num_counting_qubits / (np.log2(10))) + 2
    a = round(np.sin(theta)**2, precision)

    if num_state_qubits>=3:
        qr_state = QuantumRegister(num_state_qubits+1)
        qr_counting = QuantumRegister(num_counting_qubits)
        qr_mcx_ancilla = QuantumRegister(num_state_qubits-2)
        qc = QuantumCircuit(qr_counting, qr_state, qr_mcx_ancilla, name=f"qae-{num_qubits}-{a}")
    else:
        qr_state = QuantumRegister(num_state_qubits+1)
        qr_counting = QuantumRegister(num_counting_qubits)
        #cr = ClassicalRegister(num_counting_qubits)
        qc = QuantumCircuit(qr_counting, qr_state, name=f"qae-{num_qubits}-{a}")

    # create the Amplitude Generator circuit
    A = A_gen(num_state_qubits, a, psi_zero, psi_one)

    # create the Quantum Operator circuit and a controlled version of it
    cQ, Q = Ctrl_Q(num_state_qubits, A)
    
    # save small example subcircuits for visualization
    global A_, Q_, cQ_, QFTI_
    if (cQ_ and Q_) == None or num_state_qubits <= 6:
        if num_state_qubits < 9: cQ_ = cQ; Q_ = Q; A_ = A
    if QFTI_ == None or num_qubits <= 5:
        if num_qubits < 9: QFTI_ = inv_qft_gate(num_counting_qubits)

    # Prepare state from A, and counting qubits with H transform 
    qc.append(A, [qr_state[i] for i in range(num_state_qubits+1)])
    for i in range(num_counting_qubits):
        qc.h(qr_counting[i])
    
    repeat = 1
    if num_state_qubits>=3:
        for j in reversed(range(num_counting_qubits)):
            iter = int(repeat/approx_factor)
            for _ in range(iter):
                qc.append(cQ, [qr_counting[j]] + [qr_state[l] for l in range(num_state_qubits+1)] +\
                          [qr_mcx_ancilla[m] for m in range(num_state_qubits-2)])
            repeat *= 2
    else:
        for j in reversed(range(num_counting_qubits)):
            iter = int(repeat/approx_factor)
            for _ in range(iter):
                qc.append(cQ, [qr_counting[j]] + [qr_state[l] for l in range(num_state_qubits+1)])
            repeat *= 2
    
    qc.barrier()

    # inverse quantum Fourier transform only on counting qubits
    qc.append(inv_qft_gate(num_counting_qubits), qr_counting)

    qc.barrier()

    # measure counting qubits - do not measure
    # qc.measure([qr_counting[m] for m in range(num_counting_qubits)], list(range(num_counting_qubits)))

    # save smaller circuit example for display
    global QC_
    if QC_ == None or num_qubits <= 5:
        if num_qubits < 9: QC_ = qc
    return qc
    
# Construct A operator that takes |0>_{n+1} to sqrt(1-a) |psi_0>|0> + sqrt(a) |psi_1>|1>
def A_gen(num_state_qubits, a, psi_zero=None, psi_one=None):

    if psi_zero==None:
        psi_zero = '0'*num_state_qubits
    if psi_one==None:
        psi_one = '1'*num_state_qubits
        
    theta = 2 * np.arcsin(np.sqrt(a))
    # Let the objective be qubit index n; state is on qubits 0 through n-1
    qc_A = QuantumCircuit(num_state_qubits+1, name=f"A")
    
    # takes state to |0>_{n} (sqrt(1-a) |0> + sqrt(a) |1>)
    qc_A.ry(theta, num_state_qubits)
    
    # takes state to sqrt(1-a) |psi_0>|0> + sqrt(a) |0>_{n}|1>
    qc_A.x(num_state_qubits)
    for i in range(num_state_qubits):
        if psi_zero[i]=='1':
            qc_A.cx(num_state_qubits,i)
    qc_A.x(num_state_qubits)
    
    # takes state to sqrt(1-a) |psi_0>|0> + sqrt(a) |psi_1>|1>
    for i in range(num_state_qubits):
        if psi_one[i]=='1':
            qc_A.cx(num_state_qubits,i)
    
    return qc_A

# Construct the grover-like operator and a controlled version of it
def Ctrl_Q(num_state_qubits, A_circ):

    # index n is the objective qubit, and indexes 0 through n-1 are state qubits
    # Basically, q0 ~ q(N-1) are control and q(N) is target in MCX
    # N=1 => cx(0, 1) / N=2 => Toffoli(0, 1, 2) / N>=3 => MCX with aux qubits
    if num_state_qubits>=3:
        qr_main = QuantumRegister(num_state_qubits+1)
        qr_mcx_ancilla = QuantumRegister(num_state_qubits-2)
        qc = QuantumCircuit(qr_main, qr_mcx_ancilla, name=f"Q")
    else:
        qc = QuantumCircuit(num_state_qubits+1, name=f"Q")
    
    temp_A = copy.copy(A_circ)
    A_gate = temp_A.to_gate()
    A_gate_inv = temp_A.inverse().to_gate()
    
    ### Each cycle in Q applies in order: -S_chi, A_circ_inverse, S_0, A_circ 
    # -S_chi
    qc.x(num_state_qubits)
    qc.z(num_state_qubits)
    qc.x(num_state_qubits)
        
    # A_circ_inverse
    qc.append(A_gate_inv, [i for i in range(num_state_qubits+1)])
        
    # S_0
    for i in range(num_state_qubits+1):
        qc.x(i)
    qc.h(num_state_qubits)
    
    if num_state_qubits>=3:
        mcx(qc, [x for x in range(num_state_qubits)], num_state_qubits, \
          [(num_state_qubits+1)+x for x in range(num_state_qubits-2)])
    else:
        qc.mcx([x for x in range(num_state_qubits)], num_state_qubits)

    qc.h(num_state_qubits)
    for i in range(num_state_qubits+1):
        qc.x(i)
        
    # A_circ
    qc.append(A_gate, [i for i in range(num_state_qubits+1)])
    
    # Create a gate out of the Q operator
    qc.to_gate(label='Q')
    
    # and also a controlled version of it
    Ctrl_Q_ = qc.control(1)
    
    # and return both
    return Ctrl_Q_, qc

# if main, execute method
if __name__ == '__main__':
    amplitude_estimation()
