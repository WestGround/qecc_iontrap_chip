import dataclasses

"""
1. gate error probability simulation
 - variable : gate error probability
 - constant : trap size, num_qubit, approximation factor

2. qubit number simulation
 - variable : number of qubits
 - constant : gate error probability, trap size, approximation factor

3. trap size simulation
 - variable : trap size
 - constant : gate error probability, number of qubits, approximation factor
"""

@dataclasses.dataclass
class Experiment:
    algorithm: str
    num_qubits: [(int, int)]
    trap_size: int

sched_exp = [
    Experiment("qft", [(8, 1)], 6)
]

trap_size_exp = [
    Experiment('qft', [(8, 1)], 2),
    Experiment('bv', [(8, 1)], 20),
    Experiment('hs', [(8, 1)], 2),
    Experiment('pe', [(8, 1)], 2),

    Experiment('mc', [(8, 1)], 2),
    Experiment('ae', [(8, 1)], 2),
    Experiment('vqe', [(8, 1)], 2),
    Experiment('grover', [(8, 1)], 2)
]

error_rate_exp = [
    Experiment('qft', [(8, 1)], 6),
    Experiment('bv', [(8, 1)], 43),
    Experiment('hs', [(8, 1)], 2),
    Experiment('pe', [(8, 1)], 4),

    Experiment('mc', [(8, 1)], 3),
    Experiment('ae', [(8, 1)], 2),
    Experiment('vqe', [(8, 1)], 12), 
    Experiment('grover', [(8, 1)], 2)
]

num_qubit_exp = [
    Experiment("qft", [(8, 1)], 6),
    Experiment("bv", [(8, 1)], 43),
    Experiment("hs", [(8, 1)], 2),
    Experiment("pe", [(8, 1)], 4),

    Experiment("ae", [(8, 1)], 2),
    Experiment("mc", [(8, 1)], 3),
    Experiment("grover", [(8, 1)], 2),
    Experiment("vqe", [(8, 1)], 12),
]