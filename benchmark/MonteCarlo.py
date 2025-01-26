"""
Monte Carlo Sampling Benchmark Program via Amplitude Estimation- Qiskit
"""

import copy
import functools
import sys
import time
import os

import numpy as np
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.polynomial import polyfit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.library.standard_gates.ry import RYGate

_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path[1:1] = [_dirname]
import mc_utils as mc_utils
from QFT import inv_qft_gate

# Benchmark Name
benchmark_name = "Monte Carlo Sampling"

np.random.seed(0)

# default function is f(x) = x^2
f_of_X = functools.partial(mc_utils.power_f, power=2)

# default distribution is gaussian distribution
p_distribution = mc_utils.gaussian_dist

verbose = False

# saved circuits and subcircuits for display
A_ = None
Q_ = None
cQ_ = None
QC_ = None
R_ = None
F_ = None
QFTI_ = None

############### Circuit Definition

def monte_carlo(num_qubits: int, approx_factor:int, mu, \
                epsilon=0.05, degree=2, method=2):
    num_state_qubits = 1
    num_counting_qubits = num_qubits - num_state_qubits - 1

    target_dist = p_distribution(num_state_qubits, mu)
    f = functools.partial(f_of_X, num_state_qubits=num_state_qubits)

    A_qr = QuantumRegister(num_state_qubits+1)
    A = QuantumCircuit(A_qr, name=f"A")

    # initialize R and F circuits
    R_qr = QuantumRegister(num_state_qubits+1)
    F_qr = QuantumRegister(num_state_qubits+1)
    R = QuantumCircuit(R_qr, name=f"R")
    F = QuantumCircuit(F_qr, name=f"F")
    
    # method 1 takes in the abitrary function f and arbitrary dist
    if method == 1:
        state_prep(R, R_qr, target_dist, num_state_qubits)
        f_on_objective(F, F_qr, f, epsilon=epsilon, degree=degree)
    # method 2 chooses to have lower circuit depth by choosing specific f and dist
    elif method == 2:
        uniform_prep(R, R_qr, num_state_qubits)
        square_on_objective(F, F_qr)

    # append R and F circuits to A
    A.append(R.to_gate(), A_qr)
    A.append(F.to_gate(), A_qr)

    # run AE subroutine given our A composed of R and F
    qc = AE_Subroutine(num_state_qubits, num_counting_qubits, approx_factor, A, method)

    # save smaller circuit example for display
    global QC_, R_, F_
    if QC_ == None or num_qubits <= 5:
        if num_qubits < 9: QC_ = qc
    if (R_ and F_) == None or num_state_qubits <= 3:
        if num_state_qubits < 5: R_ = R; F_ = F

    return qc    
                    
def f_on_objective(qc, qr, f, epsilon=0.05, degree=2):
    """
    Assume last qubit is the objective. Function f is evaluated on first n-1 qubits
    """
    num_state_qubits = qc.num_qubits - 1
    c_star = (2*epsilon)**(1/(degree+1))
    
    f_ = functools.partial(f, num_state_qubits=num_state_qubits)
    zeta_ = functools.partial(mc_utils.zeta_from_f, func=f_, epsilon=epsilon, degree=degree, c=c_star)
    
    x_eval = np.linspace(0.0, 2**(num_state_qubits) - 1, num= degree+1)
    poly = Polynomial(polyfit(x_eval, zeta_(x_eval), degree))
    
    b_exp = mc_utils.binary_expansion(num_state_qubits, poly)
    
    for controls in b_exp.keys():
        theta = 2*b_exp[controls]
        controls = list(controls)
        if len(controls)==0:
            qc.ry(-theta, qr[num_state_qubits])
        else:
            # define a MCRY gate:
            # this does the same thing as qc.mcry, but is clearer in the circuit printing
            MCRY = RYGate(-theta).control(len(controls))
            qc.append(MCRY, [*(qr[i] for i in controls), qr[num_state_qubits]])

def square_on_objective(qc, qr):
    """
    Assume last qubit is the objective.
    Shifted square wave function: if x is even, f(x) = 0; if x i s odd, f(x) = 1
    """
    num_state_qubits = qc.num_qubits - 1
    for control in range(num_state_qubits):
        qc.cx(control, num_state_qubits)

def state_prep(qc, qr, target_dist, num_state_qubits):
    """
    Use controlled Ry gates to construct the superposition Sum \sqrt{p_i} |i>
    """
    r_probs = mc_utils.region_probs(target_dist, num_state_qubits)
    regions = r_probs.keys()
    r_norm = {}
    
    for r in regions:
        num_controls = len(r) - 1
        super_key = r[:num_controls]

        if super_key=='':
            r_norm[super_key] = 1
        elif super_key == '1':
            r_norm[super_key] = r_probs[super_key]
            r_norm['0'] = 1-r_probs[super_key]
        else:
            try:
                r_norm[super_key] = r_probs[super_key]
                
            except KeyError:
                r_norm[super_key] = r_norm[super_key[:num_controls-1]] - r_probs[super_key[:num_controls-1] + '1']
        
        
        norm = r_norm[super_key]
        p = 0
        if norm != 0:
            p = r_probs[r] / norm
        theta = 2*np.arcsin(np.sqrt(p))
        
        if r == '1':
            qc.ry(-theta, num_state_qubits-1)
        else:
            controls = [qr[num_state_qubits-1 - i] for i in range(num_controls)]
            
            # define a MCRY gate:
            # this does the same thing as qc.mcry, but is clearer in the circuit printing
            MCRY = RYGate(-theta).control(num_controls, ctrl_state=r[:-1])
            qc.append(MCRY, [*controls, qr[num_state_qubits-1 - num_controls]])

def uniform_prep(qc, qr, num_state_qubits):
    """
    Generates a uniform distribution over all states
    """
    for i in range(num_state_qubits):
        qc.h(i)
            
def AE_Subroutine(num_state_qubits, num_counting_qubits, approx_factor, A_circuit, method):

    num_qubits = num_state_qubits + num_counting_qubits
    
    if num_state_qubits>=3:
        qr_state = QuantumRegister(num_state_qubits+1)
        qr_counting = QuantumRegister(num_counting_qubits)
        qr_mcx_ancilla = QuantumRegister(num_state_qubits-2)
        qc = QuantumCircuit(qr_state, qr_counting, qr_mcx_ancilla, name=f"qmc({method})-{num_qubits}-{0}")
    else:
        qr_state = QuantumRegister(num_state_qubits+1)
        qr_counting = QuantumRegister(num_counting_qubits)
        #cr = ClassicalRegister(num_counting_qubits)
        qc = QuantumCircuit(qr_state, qr_counting, name=f"qmc({method})-{num_qubits}-{0}")

    A = A_circuit
    cQ, Q = Ctrl_Q(num_state_qubits, A)

    # save small example subcircuits for visualization
    global A_, Q_, cQ_, QFTI_
    if (cQ_ and Q_) == None or num_state_qubits <= 6:
        if num_state_qubits < 9: cQ_ = cQ; Q_ = Q
    if A_ == None or num_state_qubits <= 3:
        if num_state_qubits < 5: A_ = A
    if QFTI_ == None or num_counting_qubits <= 3:
        if num_counting_qubits < 4: QFTI_ = inv_qft_gate(num_counting_qubits)

    # Prepare state from A, and counting qubits with H transform 
    qc.append(A, qr_state)
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
    
    #qc.measure([qr_counting[m] for m in range(num_counting_qubits)], list(range(num_counting_qubits)))
    
    return qc
            
            
###############################
   
# Construct the grover-like operator and a controlled version of it
def Ctrl_Q(num_state_qubits, A_circ):

    # index n is the objective qubit, and indexes 0 through n-1 are state qubits
    if num_state_qubits>=3:
        qc = QuantumCircuit((num_state_qubits+1)+(num_state_qubits-2), name=f"Q")
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

    return Ctrl_Q_, qc

def mcx(qc: QuantumCircuit, control_qubits: list, target_qubit: int, ancilla: list):
    assert(len(control_qubits)>2)
    assert(len(control_qubits)==(len(ancilla)+2))
    num_control = len(control_qubits)

    qc.ccx(control_qubits[0], control_qubits[1], ancilla[0])

    for i in range(num_control-3):
        qc.ccx(control_qubits[2+i], ancilla[i], ancilla[i+1])

    qc.ccx(control_qubits[num_control-1], ancilla[num_control-3], target_qubit)

    for i in reversed(range(num_control-3)):
        qc.ccx(control_qubits[2+i], ancilla[i], ancilla[i+1])

    qc.ccx(control_qubits[0], control_qubits[1], ancilla[0])

# if main, execute method
if __name__ == '__main__':
    monte_carlo(10, mu)
