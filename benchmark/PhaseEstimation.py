"""
Phase Estimation Benchmark Program - Qiskit
"""

import sys
import os
import math, random
from decimal import Decimal, getcontext

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path[1:1] = [_dirname]
from QFT import inv_qft_gate

# Benchmark Name
benchmark_name = "Phase Estimation"

np.random.seed(0)
getcontext().prec = 30000

verbose = False

# saved subcircuits circuits for printing
QC_ = None
QFTI_ = None
U_ = None

############### Circuit Definition

def phase_estimation(num_qubits, theta):
    
    qr = QuantumRegister(num_qubits)
    
    num_counting_qubits = num_qubits - 1 # only 1 state qubit
    
    #cr = ClassicalRegister(num_counting_qubits)
    qc = QuantumCircuit(qr, name=f"qpe-{num_qubits}-{theta}")

    # initialize counting qubits in superposition
    for i in range(num_counting_qubits):
        qc.h(qr[i])

    # change to |1> in state qubit, so phase will be applied by cphase gate
    qc.x(num_counting_qubits)

    qc.barrier()

    repeat = 1
    for j in reversed(range(num_counting_qubits)):
        # controlled operation: adds phase exp(i*2*pi*theta*repeat) to the state |1>
        #                       does nothing to state |0>
        cp, _ = CPhase(2*np.pi*theta, repeat)
        qc.append(cp, [j, num_counting_qubits])
        repeat *= 2

    #Define global U operator as the phase operator
    _, U = CPhase(2*np.pi*theta, 1)

    qc.barrier()
    
    # inverse quantum Fourier transform only on counting qubits
    qc.append(inv_qft_gate(num_counting_qubits), qr[:num_counting_qubits])
    
    qc.barrier()
    
    # measure counting qubits - do not measure
    # qc.measure([qr[m] for m in range(num_counting_qubits)], list(range(num_counting_qubits)))

    # save smaller circuit example for display
    global QC_, U_, QFTI_
    if QC_ == None or num_qubits <= 5:
        if num_qubits < 9: QC_ = qc
    if U_ == None or num_qubits <= 5:
        if num_qubits < 9: U_ = U
    if QFTI_ == None or num_qubits <= 5:
        if num_qubits < 9: QFTI_ = inv_qft_gate(num_counting_qubits)
    return qc

#Construct the phase gates and include matching gate representation as readme circuit
def CPhase(angle, exponent):
    qc = QuantumCircuit(1, name=f"U^{exponent}")
    phase_angle = float(Decimal(exponent)) * angle
    if phase_angle==float('inf') or phase_angle==float('-inf'):
        phase_angle = 2*(math.pi)*random.random()
    else:
        phase_angle = phase_angle - int(phase_angle/(2*math.pi))*(2*math.pi)
    qc.p(phase_angle, 0)
    phase_gate = qc.to_gate().control(1)

    return phase_gate, qc

# if main, execute method
if __name__ == '__main__':
    theta = 2*math.pi*random.random()
    phase_estimation(10, theta)
