use std::collections::HashMap;
use std::fmt::Write;
use crate::circuit::NodeType::H;

#[derive(Clone)]
pub struct Circuit {
    pub qubits: Vec<Vec<u8>>,
    pub cx: HashMap<(usize, usize), (usize, usize)>,
    pub parameters: Vec<f64>,
    pub depths: HashMap<(usize, usize), usize>,
    pub approx_factor: usize,
}

impl Circuit {
    pub fn new(num_qubits: usize) -> Self {
        let qubits = vec![vec![]; num_qubits];
        let cx = HashMap::new();
        let parameters = vec![];
        let depths = HashMap::new();
        let approx_factor = 1;
        Self { qubits, cx, parameters, depths, approx_factor }
    }
    pub fn from_str(input: &str) -> Self {
        let mut lines = input.lines();
        let mut first_line = lines.next().unwrap();
        let mut parts = first_line.split_whitespace();
        let num_qubits: usize = parts.next().unwrap().parse().unwrap();
        let approx_factor: usize = parts.next().unwrap().parse().unwrap();
        let mut qubits = vec![vec![]; num_qubits];
        let mut cx = HashMap::new();
        let mut parameters = vec![];
        for line in lines {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let node_type = match parts[0] {
                "x" => NodeType::X,
                "y" => NodeType::Y,
                "z" => NodeType::Z,
                "h" => NodeType::H,
                "s" => NodeType::S,
                "sdg" => NodeType::SDG,
                "rz" => NodeType::RZ,
                "cx" => NodeType::CX,
                name => panic!("Unknown gate: {}", name)
            };
            let node = match node_type {
                NodeType::CX =>
                    Node { node_type, target: [parts[1].parse().unwrap(), parts[2].parse().unwrap()], parameter: 0.0 },
                NodeType::RZ =>
                    Node { node_type, target: [parts[1].parse().unwrap(), 0], parameter: parts[2].parse().unwrap() },
                _ =>
                    Node { node_type, target: [parts[1].parse().unwrap(), 0], parameter: 0.0 }
            };
            match node.node_type {
                NodeType::CX => {
                    cx.insert((node.target[0], qubits[node.target[0]].len()), (node.target[1], qubits[node.target[1]].len()));
                    cx.insert((node.target[1], qubits[node.target[1]].len()), (node.target[0], qubits[node.target[0]].len()));
                    qubits[node.target[0]].push(b'c');
                    qubits[node.target[1]].push(b'e');
                }
                NodeType::RZ => {
                    qubits[node.target[0]].push(b'r');
                    parameters.push(node.parameter);
                }
                _ => qubits[node.target[0]].push(node.node_type.into_byte())
            }
        }
        Self { qubits, cx, parameters, depths: HashMap::new(), approx_factor }
    }

    pub fn from_str_native(input: &str) -> Self {
        let mut lines = input.lines();
        let mut first_line = lines.next().unwrap();
        let mut parts = first_line.split_whitespace();
        let num_qubits: usize = parts.next().unwrap().parse().unwrap();
        let approx_factor: usize = parts.next().unwrap().parse().unwrap();
        let mut qubits = vec![vec![]; num_qubits];
        let mut cx = HashMap::new();
        let mut depths = HashMap::new();
        for line in lines {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let target: usize = parts[1].parse().unwrap();
            match parts[0] {
                "gpi" => qubits[target].push(b'g'),
                "gpi2" => qubits[target].push(b'p'),
                "rz" => {
                    match parts.get(2) {
                        Some(depth) => {
                            let depth = depth.parse().unwrap();
                            depths.insert((target, qubits[target].len()), depth);
                            qubits[target].push(b'r');
                        }
                        None => { // NO-QEC
                            qubits[target].push(b'r');
                        }
                    }
                }
                "ms" => {
                    let target2: usize = parts[2].parse().unwrap();
                    cx.insert((target, qubits[target].len()), (target2, qubits[target2].len()));
                    cx.insert((target2, qubits[target2].len()), (target, qubits[target].len()));
                    qubits[target].push(b'm');
                    qubits[target2].push(b'm');
                }
                name => panic!("Unknown gate: {}", name)
            }
        }
        Self { qubits, cx, parameters: vec![], depths, approx_factor }
    }

    pub fn into_str(&self) -> String {
        let mut ret = format!("{}\n", self.qubits.len());
        let mut pos = vec![0; self.qubits.len()];
        let mut stack: Vec<_> = self.qubits.iter().map(|q| q.len()).enumerate().collect();
        while !stack.is_empty() {
            let &(qubit_index, node_limit) = stack.last().unwrap();
            if pos[qubit_index] == node_limit {
                stack.pop();
                continue;
            } else {
                let node_index = pos[qubit_index];
                let node = self.qubits[qubit_index][node_index];
                match node {
                    b'g' => {
                        write!(ret, "gpi {}\n", qubit_index).unwrap();
                        pos[qubit_index] += 1;
                    }
                    b'p' => {
                        write!(ret, "gpi2 {}\n", qubit_index).unwrap();
                        pos[qubit_index] += 1;
                    }
                    b'm' => {
                        let partner = self.cx[&(qubit_index, node_index)];
                        if pos[partner.0] == partner.1 {
                            write!(ret, "ms {} {}\n", qubit_index, partner.0).unwrap();
                            pos[qubit_index] += 1;
                            pos[partner.0] += 1;
                        } else {
                            stack.push((partner.0, partner.1));
                            continue;
                        }
                    }
                    b'r' => {
                        write!(ret, "rz {} {}\n", qubit_index, self.depths[&(qubit_index, node_index)]).unwrap();
                        pos[qubit_index] += 1;
                    }
                    _ => unreachable!()
                }
            }
        }
        ret
    }
}
#[derive(PartialEq, Copy, Clone, Debug)]
enum NodeType {
    X,
    Y,
    Z,
    H,
    S,
    SDG,
    RZ,
    CX,
}

impl NodeType {
    fn into_byte(self) -> u8 {
        match self {
            NodeType::X => b'x',
            NodeType::Y => b'y',
            NodeType::Z => b'z',
            NodeType::H => b'h',
            NodeType::S => b's',
            NodeType::SDG => b'd',
            NodeType::RZ => b'r',
            NodeType::CX => b'c',
        }
    }
}
#[derive(PartialEq, Copy, Clone, Debug)]
struct Node {
    node_type: NodeType,
    target: [usize; 2],
    parameter: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_empty_circuit() {
        Circuit::from_str("5\n");
    }

    #[test]
    fn test_single_qubit() {
        let circuit = Circuit::from_str("5\nh 0");
        assert_eq!(circuit.qubits[0][0], b'h');
        let circuit = Circuit::from_str("5\nh 1");
        assert_eq!(circuit.qubits[1][0], b'h');
    }

    #[test]
    fn test_two_qubit() {
        let circuit = Circuit::from_str("5\ncx 0 1");
        assert_eq!(circuit.qubits[0][0], b'c');
        assert_eq!(circuit.qubits[1][0], b'e');
        assert_eq!(circuit.cx[&(0, 0)], (1, 0));
        assert_eq!(circuit.cx[&(1, 0)], (0, 0));
    }

    #[test]
    fn test_rz() {
        let circuit = Circuit::from_str("5\nrz 0 2");
        assert_eq!(circuit.qubits[0][0], b'r');
        assert_eq!(circuit.parameters[0], 2.0);
    }
}