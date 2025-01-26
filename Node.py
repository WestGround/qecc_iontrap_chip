import math
from dataclasses import dataclass
from dataclasses import field

@dataclass
class Node:
    name: str
    qubits : [int]
    parameter: [float] = field(default_factory=list)
    duration: [int] = 1
    clifford: bool = True

    def __str__(self):
        applied_qubit = ""
        for qubit_index in self.qubits:
            applied_qubit += "Q%d "%qubit_index

        param_str = ""
        for param in self.parameter:
            if isinstance(param, float):
                param_str += f"{param:.4f}({param/math.pi:.4f}π) "
            else:
                assert(isinstance(param, int))
                param_str += f"{param} "
        param_str = param_str[:-1]

        return f"{self.name} {param_str} {applied_qubit}"
