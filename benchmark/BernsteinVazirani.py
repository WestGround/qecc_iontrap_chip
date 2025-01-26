"""
Bernstein-Vazirani Benchmark Program - Qiskit
"""

import sys
import time
import os

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path[1:1] = [_dirname]

# Benchmark Name
benchmark_name = "Bernstein-Vazirani"

np.random.seed(0)

verbose = False

# Variable for number of resets to perform after mid circuit measurements
num_resets = 1

# saved circuits for display
QC_ = None
Uf_ = None

############### Circuit Definition

def create_oracle(num_qubits, input_size, secret_int):
    # Initialize first n qubits and single ancilla qubit
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr, name=f"Uf")

    # perform CX for each qubit that matches a bit in secret string
    s = ('{0:0' + str(input_size) + 'b}').format(secret_int)
    for i_qubit in range(input_size):
        if s[input_size - 1 - i_qubit] == '1':
            qc.cx(qr[i_qubit], qr[input_size])

    return qc

def bernstein_vazirani(num_qubits, secret_int, method = 1):
    
    # size of input is one less than available qubits
    input_size = num_qubits - 1

    if method == 1:
        # allocate qubits
        qr = QuantumRegister(num_qubits);
        qc = QuantumCircuit(qr, name=f"bv({method})-{num_qubits}-{secret_int}")

        # put ancilla in |1> state
        qc.x(qr[input_size])

        # start with Hadamard on all qubits, including ancilla
        for i_qubit in range(num_qubits):
             qc.h(qr[i_qubit])

        qc.barrier()

        #generate Uf oracle
        Uf = create_oracle(num_qubits, input_size, secret_int)
        qc.append(Uf,qr)

        qc.barrier()

        # start with Hadamard on all qubits, including ancilla
        for i_qubit in range(num_qubits):
             qc.h(qr[i_qubit])

        # uncompute ancilla qubit, not necessary for algorithm
        qc.x(qr[input_size])

        qc.barrier()

        # measure all data qubits
        #for i in range(input_size):
        #    qc.measure(i, i)

        global Uf_
        if Uf_ == None or num_qubits <= 6:
            if num_qubits < 9: Uf_ = Uf

    elif method == 2:
        # allocate qubits
        qr = QuantumRegister(2); qc = QuantumCircuit(qr, name="main")

        # put ancilla in |-> state
        qc.x(qr[1])
        qc.h(qr[1])

        qc.barrier()

        # perform CX for each qubit that matches a bit in secret string
        s = ('{0:0' + str(input_size) + 'b}').format(secret_int)
        for i in range(input_size):
            if s[input_size - 1 - i] == '1':
                qc.h(qr[0])
                qc.cx(qr[0], qr[1])
                qc.h(qr[0])
            #qc.measure(qr[0], cr[i])

            # Perform num_resets reset operations
            qc.reset([0]*num_resets)

    # save smaller circuit example for display
    global QC_
    if QC_ == None or num_qubits <= 6:
        if num_qubits < 9: QC_ = qc

    # return a handle on the circuit
    #qc.draw(output='mpl', filename="bv(%d)"%num_qubits)
    return qc

# if main, execute method
if __name__ == '__main__':
    bernstein_vazirani(1024, 0)
   
