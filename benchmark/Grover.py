"""
Grover's Search Benchmark Program - Qiskit
"""

import sys, os
import math

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

from MonteCarlo import mcx

_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path[1:1] = [_dirname]

# Benchmark Name
benchmark_name = "Grover's Search"

np.random.seed(0)

verbose = False

# saved circuits for display
QC_ = None
grover_oracle = None
diffusion_operator = None
 
# for validating the implementation of an mcx shim  
_use_mcx_shim = False 
called = False
called2 = False

############### Circuit Definition

def grover_search(num_qubits, marked_item, n_iterations, approx_factor):

    # allocate qubits
    #cr = ClassicalRegister(num_qubits);
    if num_qubits>=4:
        qr = QuantumRegister(num_qubits+(num_qubits-3))
    else:
        qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr, name=f"grovers-{num_qubits}-{marked_item}")

    # Start with Hadamard on all qubits
    for i_qubit in range(num_qubits):
        qc.h(qr[i_qubit])

    # loop over the estimated number of iterations
    iter = int(n_iterations/approx_factor)
    for _ in range(iter):

        qc.barrier()
    
        # add the grover oracle
        qc.append(add_grover_oracle(num_qubits, marked_item).to_instruction(), qr)
        
        # add the diffusion operator
        qc.append(add_diffusion_operator(num_qubits).to_instruction(), qr)

    qc.barrier()
        
    # measure all qubits - do not measure
    # qc.measure(qr, cr)

    # save smaller circuit example for display
    global QC_    
    if QC_ == None or num_qubits <= 5:
        if num_qubits < 9: QC_ = qc
        
    # return a handle on the circuit
    return qc

############## Grover Oracle

def add_grover_oracle(num_qubits, marked_item):
    global grover_oracle
    
    marked_item_bits = format(marked_item, f"0{num_qubits}b")[::-1]

    qr = QuantumRegister(num_qubits)
    if num_qubits>=4:
        qr_mcx_ancilla = QuantumRegister(num_qubits-3)
        qc = QuantumCircuit(qr, qr_mcx_ancilla, name="oracle")
    else:
        qc = QuantumCircuit(qr, name="oracle")

    for (q, bit) in enumerate(marked_item_bits):
        if not int(bit):
            qc.x(q)

    qc.h(num_qubits - 1)
    
    if num_qubits>=4:
        mcx(qc, [x for x in range(num_qubits - 1)], num_qubits - 1, [(num_qubits)+x for x in range(num_qubits-3)])
    else:
        qc.mcx([x for x in range(num_qubits - 1)], num_qubits - 1)
        
    qc.h(num_qubits - 1)

    qc.barrier()

    for (q, bit) in enumerate(marked_item_bits):
        if not int(bit):
            qc.x(q)

    if grover_oracle == None or num_qubits <= 5:
        if num_qubits < 9: grover_oracle = qc

    return qc

############## Grover Diffusion Operator

def add_diffusion_operator(num_qubits):
    global diffusion_operator

    qr = QuantumRegister(num_qubits)
    if num_qubits>=4:
        qr_mcx_ancilla = QuantumRegister(num_qubits-3)
        qc = QuantumCircuit(qr, qr_mcx_ancilla, name="diffuser")
    else:
        qc = QuantumCircuit(qr, name="diffuser")

    for i_qubit in range(num_qubits):
        qc.h(qr[i_qubit])
    for i_qubit in range(num_qubits):
        qc.x(qr[i_qubit])
    qc.h(num_qubits - 1)

    if num_qubits>=4:
        mcx(qc, [x for x in range(num_qubits - 1)], num_qubits - 1, [(num_qubits)+x for x in range(num_qubits-3)])
    else:
        qc.mcx([x for x in range(num_qubits - 1)], num_qubits - 1)
        
    qc.h(num_qubits - 1)

    qc.barrier()

    for i_qubit in range(num_qubits):
        qc.x(qr[i_qubit])
    for i_qubit in range(num_qubits):
        qc.h(qr[i_qubit])
        
    if diffusion_operator == None or num_qubits <= 5:
        if num_qubits < 9: diffusion_operator = qc

    return qc

############### MCX shim

# single cx / cu1 unit for mcx implementation
def add_cx_unit(qc: QuantumCircuit, cxcu1_unit, controls, target):
    num_controls = len(controls)
    i_qubit = cxcu1_unit[1]
    j_qubit = cxcu1_unit[0]
    theta = cxcu1_unit[2]
    
    if j_qubit != None:
        qc.cx(controls[j_qubit], controls[i_qubit]) 
    #qc.cu1(theta, controls[i_qubit], target)
    qc.cu(0, 0, theta, 0, controls[i_qubit], target)

    i_qubit = i_qubit - 1
    if j_qubit == None:
        j_qubit = i_qubit + 1
    else:
        j_qubit = j_qubit - 1
        
    if theta < 0:
        theta = -theta
    
    new_units = []
    if i_qubit >= 0:
        new_units += [ [ j_qubit, i_qubit, -theta ] ]
        new_units += [ [ num_controls - 1, i_qubit, theta ] ]
        
    return new_units

# mcx recursion loop 
def add_cxcu1_units(qc, cxcu1_units, controls, target):
    new_units = []
    for cxcu1_unit in cxcu1_units:
        new_units += add_cx_unit(qc, cxcu1_unit, controls, target)
    cxcu1_units.clear()
    return new_units

# mcx gate implementation: brute force and inefficent
# start with a single CU1 on last control and target
# and recursively expand for each additional control
def add_mcx(qc, controls, target):
    num_controls = len(controls)
    theta = np.pi / 2**num_controls
    qc.h(target)
    cxcu1_units = [ [ None, num_controls - 1, theta] ]
    while len(cxcu1_units) > 0:
        cxcu1_units += add_cxcu1_units(qc, cxcu1_units, controls, target)
    qc.h(target)

# if main, execute method
if __name__ == '__main__':
    N = 10
    iter = math.pi/4*math.sqrt(2**N)
    grover_search(N, 0, iter)
