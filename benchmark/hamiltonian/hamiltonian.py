import os, math
import numpy as np

from qiskit import QuantumCircuit, QuantumRegister

np.random.seed(0)

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

def xxyyzz_gate(tau):
    alpha = tau
    beta = tau
    gamma = tau
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr)
    qc.rz(math.pi/2, qr[1])
    qc.cx(qr[1], qr[0])
    qc.rz(math.pi*gamma - math.pi/2, qr[0])
    qc.ry(math.pi/2 - math.pi*alpha, qr[1])
    qc.cx(qr[0], qr[1])
    qc.ry(math.pi*beta - math.pi/2, qr[1])
    qc.cx(qr[1], qr[0])
    qc.rz((-1)*math.pi/2, qr[0])
    return qc

def hamiltonian(n, K=3, t=1e-2, w = 10):
    '''
    :param n: Number of qubits
    :param K: The Trotterization order
    :param t: duration of simulation
    '''
    tau = t / K

    qc = QuantumCircuit(n, n)
    h_x = np.random.random(n) * 2 - 1
    h_z = np.random.random(n) * 2 - 1

    for qubit in range(0, n, 2):
        qc.x(qubit)
    for trotter_step in range(K):
        for i in range(n):
            angle = 2*tau*w*h_x[i]
            if angle<0:
                angle += (2*math.pi)*(abs(int(angle/(2*math.pi)))+1)
            elif angle>=(2*math.pi):
                angle -= (2*math.pi)*(abs(int(angle/(2*math.pi))))
            assert(angle>=0 and angle<2*math.pi)
            qc.rx(angle, i)
        for i in range(n):
            angle = 2*tau*w*h_z[i]
            if angle<0:
                angle += (2*math.pi)*(abs(int(angle/(2*math.pi)))+1)
            elif angle>=(2*math.pi):
                angle -= (2*math.pi)*(abs(int(angle/(2*math.pi))))
            assert(angle>=0 and angle<2*math.pi)
            qc.rz(angle, i)
        for j in range(2):
            for i in range(j, n - 1, 2):
                i0 = i
                i1 = (i+1)%n
                # original - qc.append(xxyyzz_gate(tau).to_gate(), [i, (i + 1) % n])
                qc.rz(math.pi/2, i1)
                qc.cx(i1, i0)
                angle = math.pi*tau - math.pi/2
                if angle<0:
                    angle += (2*math.pi)*(abs(int(angle/(2*math.pi)))+1)
                elif angle>=(2*math.pi):
                    angle -= (2*math.pi)*(abs(int(angle/(2*math.pi))))
                assert(angle>=0 and angle<2*math.pi)
                qc.rz(angle, i0)
                qc.ry((-1)*angle, i1)
                qc.cx(i0, i1)
                qc.ry(angle, i1)
                qc.cx(i1, i0)
                qc.rz((-1)*math.pi/2, i0)
    return qc

if __name__ == '__main__':
    qc = hamiltonian(100)
