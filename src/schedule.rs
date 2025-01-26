use std::collections::{HashMap, VecDeque};
use rand::Rng;
use crate::circuit::Circuit;
use crate::{SQ_TIME, TQ_TIME, SWAP_TIME, SHUTTLE_TIME, MEASURE_TIME};

type Physical = usize;
type Logical = usize;

struct CircuitState {
    offset: usize,
    shuttle_right: bool,
}


fn apply_non_clifford(qubit: &mut Vec<String>, timestamp: &mut Vec<usize>, depth: usize, sector_time: &mut usize) {
    let mut rng = rand::thread_rng();
    for _ in 0..depth {
        qubit.push("p".to_string());
        timestamp.push(*sector_time);
        *sector_time += SQ_TIME;
        qubit.push("m".to_string());
        timestamp.push(*sector_time);
        *sector_time += TQ_TIME;
        qubit.push("p".to_string());
        timestamp.push(*sector_time);
        *sector_time += SQ_TIME;
        qubit.push("p".to_string());
        timestamp.push(*sector_time);
        *sector_time += SQ_TIME;
        qubit.push("M".to_string());
        timestamp.push(*sector_time);
        *sector_time += MEASURE_TIME;
        if rng.gen_bool(0.5) { // 50% success
            break;
        }
    }
}

pub fn mark_scheduler(circuit: &Circuit, sector_size: usize, empty_sector: usize, qec_time: usize) -> (Schedule, usize) {
    let num_qubits = circuit.qubits.len();
    let mut sched = Schedule { qubits: vec![vec![]; num_qubits], timestamp: vec![vec![]; num_qubits], cx: HashMap::default() };
    let (mut physical_to_logical, mut logical_to_physical) = initialize_qubit_pos(num_qubits, circuit);
    let mut state = CircuitState { offset: 0, shuttle_right: true };
    let mut max_shift = 0;
    if num_qubits%(sector_size+1)==0{
        max_shift = empty_sector * sector_size + (empty_sector-1);
    } else {
        max_shift = (sector_size + 1) - num_qubits % (sector_size + 1) + empty_sector * (sector_size + 1) - 1;
    }
    let sector_num = (num_qubits + max_shift) / (sector_size + 1) * 2 + 1;
    let mut start_time = 0;
    let mut current_nodes = vec![0; num_qubits];

    let mut frontier = VecDeque::new();
    for i in 0..num_qubits {
        if circuit.qubits[i].is_empty() {
            // Empty circuit
            continue;
        }
        if circuit.qubits[i][0] == b'm' {
            // MS gate
            let partner = circuit.cx[&(i, 0)];
            if partner.0 > i && partner.1 == 0 {
                frontier.push_back((i, 0));
            }
        } else {
            frontier.push_back((i, 0))
        }
    }

    while !frontier.is_empty() { // Schedule gates
        let mut new_frontier = VecDeque::new();
        let mut sector_time = vec![start_time; sector_num];
        while !frontier.is_empty() {
            let (qubit_index, node_index) = frontier.pop_front().unwrap();
            let node = circuit.qubits[qubit_index][node_index];
            let current_sector = compute_current_sector(logical_to_physical[qubit_index] + state.offset, sector_size);
            match node {
                b'g' | b'p' => { // GPI, GPI2
                    if current_sector % 2 == 0 { // Horizontal sector
                        sched.qubits[qubit_index].push((node as char).to_string());
                        sched.timestamp[qubit_index].push(sector_time[current_sector]);
                        sector_time[current_sector] += SQ_TIME;
                        advance(&circuit, qubit_index, &mut current_nodes, &mut frontier);
                    } else {
                        new_frontier.push_back((qubit_index, node_index));
                    }
                }
                b'm' => { // MS
                    let partner = circuit.cx[&(qubit_index, node_index)];
                    let partner_sector = compute_current_sector(logical_to_physical[partner.0] + state.offset, sector_size);
                    if current_sector == partner_sector {
                        //sched.qubits[qubit_index].push(format!("m{}", partner.0));
                        //sched.qubits[partner.0].push(format!("m{}", qubit_index));
                        sched.qubits[qubit_index].push(format!("m"));
                        sched.qubits[partner.0].push(format!("m"));
                        sched.cx.insert((qubit_index, sched.qubits[qubit_index].len() - 1),
                                        (partner.0, sched.qubits[partner.0].len() - 1));
                        sched.timestamp[qubit_index].push(sector_time[current_sector]);
                        sched.timestamp[partner.0].push(sector_time[current_sector]);
                        sector_time[current_sector] += TQ_TIME;
                        advance(&circuit, qubit_index, &mut current_nodes, &mut frontier);
                        advance(&circuit, partner.0, &mut current_nodes, &mut frontier);
                    } else {
                        new_frontier.push_back((qubit_index, node_index));
                    }
                }
                b'r' => { // RZ (non Clifford)
                    if current_sector % 2 == 1 {
                        apply_non_clifford(&mut sched.qubits[qubit_index], &mut sched.timestamp[qubit_index], circuit.depths[&(qubit_index, node_index)], &mut sector_time[current_sector]);
                        advance(&circuit, qubit_index, &mut current_nodes, &mut frontier);
                    } else {
                        new_frontier.push_back((qubit_index, node_index));
                    }
                }
                _ => panic!("Unknown node: {}", node)
            }
        }

        frontier = new_frontier;

        let mut move_qubit = vec![false; num_qubits];
        let mut stop_qubit = vec![false; num_qubits];

        for i in 0..frontier.len() { // Mark move and stop
            let (qubit_index, node_index) = frontier[i];
            let node = circuit.qubits[qubit_index][node_index];
            match node {
                b'g' | b'p' => {} // Clifford inside vertical
                b'm' => {
                    let partner = circuit.cx[&(qubit_index, node_index)];
                    let physical = logical_to_physical[qubit_index];
                    let partner_physical = logical_to_physical[partner.0];
                    if physical < partner_physical {
                        if state.shuttle_right {
                            move_qubit[logical_to_physical[qubit_index]] = true;
                            stop_qubit[logical_to_physical[partner.0]] = true;
                        } else {
                            move_qubit[logical_to_physical[partner.0]] = true;
                            stop_qubit[logical_to_physical[qubit_index]] = true;
                        }
                    } else {
                        if state.shuttle_right {
                            move_qubit[logical_to_physical[partner.0]] = true;
                            stop_qubit[logical_to_physical[qubit_index]] = true;
                        } else {
                            move_qubit[logical_to_physical[qubit_index]] = true;
                            stop_qubit[logical_to_physical[partner.0]] = true;
                        }
                    }
                }
                b'r' => move_qubit[logical_to_physical[qubit_index]] = true, // Non Clifford
                _ => unreachable!()
            }
        }

        let (mut first_phys, mut last_phys, mut is_horizontal) = compute_phys(sector_size, state.offset);

        if !state.shuttle_right && state.offset % (sector_size + 1) != 0 && is_horizontal { // If shuttle left and first sector is not full, skip first sector
            first_phys = last_phys;
            last_phys += 1;
            is_horizontal = false;
        }

        while last_phys <= num_qubits {
            if !is_horizontal {
                first_phys = last_phys;
                last_phys += sector_size;
                is_horizontal = true;
                continue;
            }
            let default_shuttle = if state.shuttle_right { last_phys - 1 } else { first_phys };
            let shuttle_physical = choose_shuttle(first_phys, last_phys, &move_qubit, &stop_qubit, state.shuttle_right);
            if shuttle_physical != default_shuttle {
                let shuttle_logical = physical_to_logical[shuttle_physical];
                let default_logical = physical_to_logical[default_shuttle];

                logical_to_physical[shuttle_logical] = default_shuttle;
                logical_to_physical[default_logical] = shuttle_physical;
                physical_to_logical[default_shuttle] = shuttle_logical;
                physical_to_logical[shuttle_physical] = default_logical;

                let current_sector = compute_current_sector(default_shuttle, sector_size);
                //sched.qubits[shuttle_logical].push(format!("s{}", default_logical));
                //sched.qubits[default_logical].push(format!("s{}", shuttle_logical));
                sched.qubits[shuttle_logical].push(format!("s"));
                sched.qubits[default_logical].push(format!("s"));
                sched.timestamp[shuttle_logical].push(sector_time[current_sector]);
                sched.timestamp[default_logical].push(sector_time[current_sector]);
                sector_time[current_sector] += SWAP_TIME;
            }
            first_phys = last_phys;
            last_phys += 1;
            is_horizontal = false;
        }

        if !state.shuttle_right && first_phys < num_qubits && is_horizontal {
            let default_shuttle = first_phys;
            let shuttle_physical = choose_shuttle(first_phys, num_qubits, &move_qubit, &stop_qubit, state.shuttle_right);
            if shuttle_physical != default_shuttle {
                let shuttle_logical = physical_to_logical[shuttle_physical];
                let default_logical = physical_to_logical[default_shuttle];

                logical_to_physical[shuttle_logical] = default_shuttle;
                logical_to_physical[default_logical] = shuttle_physical;
                physical_to_logical[default_shuttle] = shuttle_logical;
                physical_to_logical[shuttle_physical] = default_logical;

                let current_sector = compute_current_sector(default_shuttle, sector_size);
                //sched.qubits[shuttle_logical].push(format!("s{}", default_logical));
                //sched.qubits[default_logical].push(format!("s{}", shuttle_logical));
                sched.qubits[shuttle_logical].push(format!("s"));
                sched.qubits[default_logical].push(format!("s"));
                sched.timestamp[shuttle_logical].push(sector_time[current_sector]);
                sched.timestamp[default_logical].push(sector_time[current_sector]);
                sector_time[current_sector] += SWAP_TIME;
            }
        }

        start_time = sector_time.into_iter().max().unwrap();

        // QEC
        let (mut first_phys, mut last_phys, mut is_horizontal) = compute_phys(sector_size, state.offset);
        if is_horizontal {
            first_phys = last_phys;
            last_phys += 1;
            is_horizontal = false;
        }
        while last_phys <= num_qubits {
            sched.qubits[physical_to_logical[first_phys]].push("q".to_string());
            sched.timestamp[physical_to_logical[first_phys]].push(start_time);
            first_phys = last_phys + sector_size;
            last_phys += sector_size + 1;
        }
        start_time += qec_time;

        // Shuttle
        for qubit in 0..num_qubits {
            sched.qubits[qubit].push("h".to_string());
            sched.timestamp[qubit].push(start_time);
        }
        start_time += SHUTTLE_TIME;

        if state.shuttle_right {
            state.offset += 1;
        } else {
            state.offset -= 1;
        }
        if state.offset == 0 || state.offset == max_shift {
            state.shuttle_right = !state.shuttle_right;
        }

        if start_time > 987654321000 {
            println!("{} TIME OUT", sector_size);
            break;
        }
    }
    for qubit in 0..num_qubits { // Push qec for error generator
        sched.qubits[qubit].push("q".to_string());
        sched.timestamp[qubit].push(start_time);
    }

    (sched, start_time)
}

/// Scheduler without QEC (input: gpi, gpi2, ms)
pub fn pmark_scheduler(circuit: &Circuit, sector_size: usize, empty_sector: usize) -> (crate::schedule::Schedule, usize) {
    let num_qubits = circuit.qubits.len();
    let mut sched = Schedule { qubits: vec![vec![]; num_qubits], timestamp: vec![vec![]; num_qubits], cx: HashMap::default() };
    let (mut physical_to_logical, mut logical_to_physical) = initialize_qubit_pos(num_qubits, circuit);
    let mut state = CircuitState { offset: 0, shuttle_right: true };
    let mut max_shift = 0;
    if num_qubits%sector_size==0{
        max_shift = empty_sector * sector_size;
    } else {
        max_shift = sector_size - num_qubits % sector_size + empty_sector * sector_size;
    }
    let sector_num = (max_shift + num_qubits) / sector_size;
    let mut start_time = 0;
    let mut current_nodes = vec![0; num_qubits];
    let mut frontier = VecDeque::new();
    for i in 0..num_qubits {
        if circuit.qubits[i].is_empty() {
            // Empty circuit
            continue;
        }
        if circuit.qubits[i][0] == b'm' {
            // MS gate
            let partner = circuit.cx[&(i, 0)];
            if partner.0 > i && partner.1 == 0 {
                frontier.push_back((i, 0));
            }
        } else {
            frontier.push_back((i, 0))
        }
    }

    while !frontier.is_empty() { // Schedule gates
        let mut new_frontier = VecDeque::new();
        let mut sector_time = vec![start_time; sector_num];
        while !frontier.is_empty() {
            let (qubit_index, node_index) = frontier.pop_front().unwrap();
            let node = circuit.qubits[qubit_index][node_index];
            let current_sector = (logical_to_physical[qubit_index] + state.offset) / sector_size;
            match node {
                b'g' | b'p' | b'r' => { // GPI, GPI2, rz
                    sched.qubits[qubit_index].push((node as char).to_string());
                    sched.timestamp[qubit_index].push(sector_time[current_sector]);
                    sector_time[current_sector] += SQ_TIME;
                    advance(&circuit, qubit_index, &mut current_nodes, &mut frontier);
                }
                b'm' => { // MS
                    let partner = circuit.cx[&(qubit_index, node_index)];
                    let partner_sector = (logical_to_physical[partner.0] + state.offset) / sector_size;
                    if current_sector == partner_sector {
                        //#sched.qubits[qubit_index].push(format!("m{}", partner.0));
                        //sched.qubits[partner.0].push(format!("m{}", qubit_index));
                        sched.qubits[qubit_index].push(format!("m"));
                        sched.qubits[partner.0].push(format!("m"));
                        sched.cx.insert((qubit_index, sched.qubits[qubit_index].len() - 1),
                                        (partner.0, sched.qubits[partner.0].len() - 1));
                        sched.timestamp[qubit_index].push(sector_time[current_sector]);
                        sched.timestamp[partner.0].push(sector_time[current_sector]);
                        sector_time[current_sector] += TQ_TIME;
                        advance(&circuit, qubit_index, &mut current_nodes, &mut frontier);
                        advance(&circuit, partner.0, &mut current_nodes, &mut frontier);
                    } else {
                        new_frontier.push_back((qubit_index, node_index));
                    }
                }
                _ => panic!("Unknown node: {}", node as char)
            }
        }

        frontier = new_frontier;

        let mut move_qubit = vec![false; num_qubits];
        let mut stop_qubit = vec![false; num_qubits];

        for i in 0..frontier.len() { // Mark move and stop
            let (qubit_index, node_index) = frontier[i];
            let node = circuit.qubits[qubit_index][node_index];
            match node {
                //b'g' | b'p' => {} // No-QEC (executable)
                b'm' => {
                    let partner = circuit.cx[&(qubit_index, node_index)];
                    let physical = logical_to_physical[qubit_index];
                    let partner_physical = logical_to_physical[partner.0];
                    if physical < partner_physical {
                        if state.shuttle_right {
                            move_qubit[logical_to_physical[qubit_index]] = true;
                            stop_qubit[logical_to_physical[partner.0]] = true;
                        } else {
                            move_qubit[logical_to_physical[partner.0]] = true;
                            stop_qubit[logical_to_physical[qubit_index]] = true;
                        }
                    } else {
                        assert!(physical > partner_physical);
                        if state.shuttle_right {
                            move_qubit[logical_to_physical[partner.0]] = true;
                            stop_qubit[logical_to_physical[qubit_index]] = true;
                        } else {
                            move_qubit[logical_to_physical[qubit_index]] = true;
                            stop_qubit[logical_to_physical[partner.0]] = true;
                        }
                    }
                }
                _ => unreachable!()
            }
        }

        let mut first_phys = 0;
        let mut last_phys = sector_size - (state.offset % sector_size);

        if !state.shuttle_right && state.offset % sector_size != 0 { // If shuttle left and first sector is not full, skip first sector
            first_phys = last_phys;
            last_phys += sector_size;
        }

        while last_phys <= num_qubits {
            let default_shuttle = if state.shuttle_right { last_phys - 1 } else { first_phys };
            let shuttle_physical = choose_shuttle(first_phys, last_phys, &move_qubit, &stop_qubit, state.shuttle_right);
            if shuttle_physical != default_shuttle {
                let shuttle_logical = physical_to_logical[shuttle_physical];
                let default_logical = physical_to_logical[default_shuttle];

                logical_to_physical[shuttle_logical] = default_shuttle;
                logical_to_physical[default_logical] = shuttle_physical;
                physical_to_logical[default_shuttle] = shuttle_logical;
                physical_to_logical[shuttle_physical] = default_logical;

                let current_sector = first_phys / sector_size;
                //sched.qubits[shuttle_logical].push(format!("s{}", default_logical));
                //sched.qubits[default_logical].push(format!("s{}", shuttle_logical));
                sched.qubits[shuttle_logical].push(format!("s"));
                sched.qubits[default_logical].push(format!("s"));
                sched.timestamp[shuttle_logical].push(sector_time[current_sector]);
                sched.timestamp[default_logical].push(sector_time[current_sector]);
                sector_time[current_sector] += SWAP_TIME;
            }
            first_phys = last_phys;
            last_phys += sector_size;
        }

        if !state.shuttle_right && first_phys < num_qubits {
            let default_shuttle = first_phys;
            let shuttle_physical = choose_shuttle(first_phys, num_qubits, &move_qubit, &stop_qubit, state.shuttle_right);
            if shuttle_physical != default_shuttle {
                let shuttle_logical = physical_to_logical[shuttle_physical];
                let default_logical = physical_to_logical[default_shuttle];

                logical_to_physical[shuttle_logical] = default_shuttle;
                logical_to_physical[default_logical] = shuttle_physical;
                physical_to_logical[default_shuttle] = shuttle_logical;
                physical_to_logical[shuttle_physical] = default_logical;

                let current_sector = first_phys / sector_size;
                //sched.qubits[shuttle_logical].push(format!("s{}", default_logical));
                //sched.qubits[default_logical].push(format!("s{}", shuttle_logical));
                sched.qubits[shuttle_logical].push(format!("s"));
                sched.qubits[default_logical].push(format!("s"));
                sched.timestamp[shuttle_logical].push(sector_time[current_sector]);
                sched.timestamp[default_logical].push(sector_time[current_sector]);
                sector_time[current_sector] += SWAP_TIME;
            }
        }

        start_time = sector_time.into_iter().max().unwrap();

        // Shuttle
        for qubit in 0..num_qubits {
            sched.qubits[qubit].push("h".to_string());
            sched.timestamp[qubit].push(start_time);
        }

        start_time += SHUTTLE_TIME;
        if state.shuttle_right {
            state.offset += 1;
        } else {
            state.offset -= 1;
        }
        if state.offset == 0 || state.offset == max_shift {
            state.shuttle_right = !state.shuttle_right;
        }

        if start_time > 987654321000 {
            println!("PMARK {} TIME OUT", sector_size);
            break;
        }
    }

    for qubit in 0..num_qubits { // Push qec for error generator
        sched.qubits[qubit].push("q".to_string());
        sched.timestamp[qubit].push(start_time);
    }
    (sched, start_time)
}

fn compute_phys(sector_size: usize, offset: usize) -> (usize, usize, bool) {
    let first_phys = 0;
    let offset_mod = offset % (sector_size + 1);
    let last_phys;
    let is_horizontal;
    if offset_mod != sector_size { // Horizontal sector
        last_phys = sector_size - offset_mod;
        is_horizontal = true;
    } else { // Vertical sector
        last_phys = 1;
        is_horizontal = false;
    }
    (first_phys, last_phys, is_horizontal)
}
fn is_executable(circuit: &Circuit, index: (usize, usize), current_nodes: &Vec<usize>) -> bool {
    let partner = circuit.cx[&index];
    current_nodes[partner.0] == partner.1
}

fn advance(circuit: &Circuit, qubit_index: usize, current_nodes: &mut Vec<usize>, frontier: &mut VecDeque<(usize, usize)>) {
    current_nodes[qubit_index] += 1;
    let node_index = current_nodes[qubit_index];
    if node_index >= circuit.qubits[qubit_index].len() {
        return;
    }
    match circuit.qubits[qubit_index][node_index] {
        b'g' | b'p' | b'r' => frontier.push_back((qubit_index, node_index)),
        b'm' => if is_executable(&circuit, (qubit_index, node_index), &current_nodes) {
            frontier.push_back((qubit_index, node_index));
        }
        _ => panic!("Unknown node: {}", circuit.qubits[qubit_index][node_index])
    }
}

fn choose_shuttle(first: Physical, last: Physical, move_qubit: &[bool], stop_qubit: &[bool], shuttle_right: bool) -> Physical {
    if shuttle_right {
        for i in (first..last).rev() {
            if move_qubit[i] {
                return i;
            }
        }
        for i in (first..last).rev() {
            if !stop_qubit[i] {
                return i;
            }
        }
        return last - 1;
    } else {
        for i in first..last {
            if move_qubit[i] {
                return i;
            }
        }
        for i in first..last {
            if !stop_qubit[i] {
                return i;
            }
        }
        return first;
    }
}

fn compute_current_sector(physical: Physical, sector_size: usize) -> usize {
    let mut sector_num = physical / (sector_size + 1) * 2;
    if physical % (sector_size + 1) == sector_size {
        sector_num += 1;
    }
    sector_num
}

pub fn ntcf_scheduler(circuit: &Circuit, sector_size: usize, empty_sector: usize, qec_time: usize) -> (Schedule, usize) {
    let num_qubits = circuit.qubits.len();
    let mut sched = Schedule { qubits: vec![vec![]; num_qubits], timestamp: vec![vec![]; num_qubits], cx: HashMap::default() };
    let (mut physical_to_logical, mut logical_to_physical) = initialize_qubit_pos(num_qubits, circuit);
    let mut state = CircuitState { offset: 0, shuttle_right: true };
    let mut max_shift = 0;
    if num_qubits%(sector_size+1)==0{
        max_shift = empty_sector * sector_size + (empty_sector-1);
    } else {
        max_shift = (sector_size + 1) - num_qubits % (sector_size + 1) + empty_sector * (sector_size + 1) - 1;
    }
    let mut start_time = 0;
    let mut current_nodes = vec![0; num_qubits];

    let mut swap_count = 0;
    let mut for_r = 0;
    let mut for_close = 0;
    let mut for_far = 0;

    while current_nodes.iter().enumerate().any(|(qubit, pos)| (*pos as usize) < circuit.qubits[qubit].len()) { // Until all qubit is at the end of the circuit
        let mut time_limit = start_time + qec_time; // All gate have to finish before time_limit (max of vertical execution)

        // Iterate vertical sectors
        let (mut first_phys, mut last_phys, mut is_horizontal) = compute_phys(sector_size, state.offset);
        if is_horizontal {
            first_phys = last_phys;
            last_phys += 1;
            is_horizontal = false;
        }
        while first_phys < num_qubits { // first_phys: qubit in a vertical sector
            let qubit_index = physical_to_logical[first_phys];
            let node_index = current_nodes[qubit_index];
            let mut sector_time = start_time;
            if node_index < circuit.qubits[qubit_index].len() {
                match circuit.qubits[qubit_index][node_index] {
                    b'r' => {
                        apply_non_clifford(&mut sched.qubits[qubit_index], &mut sched.timestamp[qubit_index], circuit.depths[&(qubit_index, node_index)], &mut sector_time);
                        current_nodes[qubit_index] += 1;
                        //println!("R q{}", qubit_index);
                    }
                    b'g' | b'p' => { // Don't execute (policy)
                    }
                    b'm' => {} // Can't execute
                    _ => unreachable!()
                }
            }
            // Apply QEC
            sched.qubits[physical_to_logical[first_phys]].push("q".to_string());
            sched.timestamp[physical_to_logical[first_phys]].push(sector_time);
            sector_time += qec_time;

            time_limit = std::cmp::max(time_limit, sector_time);
            first_phys += sector_size + 1;
        }

        let horizontal_limit = std::cmp::min(time_limit, start_time + 2 * qec_time);
        // Iterate horizontal sectors
        let (mut first_phys, mut last_phys, mut is_horizontal) = compute_phys(sector_size, state.offset);
        if !is_horizontal {
            first_phys = last_phys;
            last_phys += sector_size;
            is_horizontal = true;
        }
        while first_phys < num_qubits {
            let mut sector_time = start_time;
            let qubits: Vec<usize> = if state.shuttle_right { // If shuttle_right consider right to left
                (first_phys..std::cmp::min(last_phys, num_qubits)).rev().collect()
            } else {
                (first_phys..std::cmp::min(last_phys, num_qubits)).collect()
            };
            while sector_time < horizontal_limit { // Execute as many gate as possible before horizontal_limit
                let mut changed = false;
                for &physical_qubit in &qubits {
                    let logical_qubit = physical_to_logical[physical_qubit];
                    loop {
                        if current_nodes[logical_qubit] == circuit.qubits[logical_qubit].len() { // This qubit is finished
                            break;
                        }
                        let node_index = current_nodes[logical_qubit];
                        let node = circuit.qubits[logical_qubit][node_index];
                        match node {
                            b'g' | b'p' => {
                                if sector_time + SQ_TIME <= horizontal_limit {
                                    sched.qubits[logical_qubit].push((node as char).to_string());
                                    sched.timestamp[logical_qubit].push(sector_time);
                                    sector_time += SQ_TIME;
                                    current_nodes[logical_qubit] += 1;
                                    changed = true;
                                    //println!("{} q{}", node as char, logical_qubit);
                                } else {
                                    break;
                                }
                            }
                            b'm' => {
                                if sector_time + TQ_TIME > horizontal_limit {
                                    break;
                                }
                                let partner = circuit.cx[&(logical_qubit, node_index)];
                                if current_nodes[partner.0] != partner.1 { // Partner has unfinished gates before ms
                                    break;
                                }
                                let current_sector = compute_current_sector(physical_qubit + state.offset, sector_size);
                                let partner_sector = compute_current_sector(logical_to_physical[partner.0] + state.offset, sector_size);
                                if current_sector != partner_sector {
                                    break;
                                }
                                //sched.qubits[logical_qubit].push(format!("m{}", partner.0));
                                //sched.qubits[partner.0].push(format!("m{}", logical_qubit));
                                sched.qubits[logical_qubit].push(format!("m"));
                                sched.qubits[partner.0].push(format!("m"));
                                sched.cx.insert((logical_qubit, sched.qubits[logical_qubit].len() - 1),
                                                (partner.0, sched.qubits[partner.0].len() - 1));
                                sched.timestamp[logical_qubit].push(sector_time);
                                sched.timestamp[partner.0].push(sector_time);
                                sector_time += TQ_TIME;
                                current_nodes[logical_qubit] += 1;
                                current_nodes[partner.0] += 1;
                                changed = true;
                                //println!("{} q{} q{}", node as char, logical_qubit, partner.0);
                            }
                            b'r' => break,
                            _ => unreachable!()
                        }
                    }
                }
                if !changed {
                    break;
                }
            }

            if sector_time + SWAP_TIME <= time_limit { // Insert swap if possible
                let default_shuttle = qubits[0];
                let (shuttle_physical, flag) = choose_ntcf(&qubits, &physical_to_logical, &logical_to_physical, &current_nodes, &circuit, sector_size);
                if default_shuttle != shuttle_physical {
                    let default_logical = physical_to_logical[default_shuttle];
                    let shuttle_logical = physical_to_logical[shuttle_physical];

                    physical_to_logical[default_shuttle] = shuttle_logical;
                    physical_to_logical[shuttle_physical] = default_logical;
                    logical_to_physical[shuttle_logical] = default_shuttle;
                    logical_to_physical[default_logical] = shuttle_physical;

                    //sched.qubits[shuttle_logical].push(format!("s{}", default_logical));
                    //sched.qubits[default_logical].push(format!("s{}", shuttle_logical));
                    sched.qubits[shuttle_logical].push(format!("s"));
                    sched.qubits[default_logical].push(format!("s"));
                    sched.timestamp[default_logical].push(sector_time);
                    sched.timestamp[shuttle_logical].push(sector_time);
                    swap_count += 1;
                    assert!(flag>0);
                    let mut m = "";
                    if flag==1 {
                        m = "R";
                        for_r += 1;
                    } else if flag==2 {
                        m = "close";
                        for_close += 1;
                    } else {
                        assert!(flag==3);
                        m = "far";
                        for_far += 1;
                    }
                    //println!("{}th SWAP for {} q{} q{}", swap_count, m, shuttle_logical, default_logical);
                }
            }
            first_phys = last_phys + 1;
            last_phys += sector_size + 1;
        }


        start_time = time_limit;
        // Shuttle
        for qubit in 0..num_qubits {
            sched.qubits[qubit].push("h".to_string());
            sched.timestamp[qubit].push(start_time);
        }
        start_time += SHUTTLE_TIME;

        if state.shuttle_right {
            state.offset += 1;
        } else {
            state.offset -= 1;
        }
        // Shuttle direction change
        if state.offset == 0 || state.offset == max_shift {
            state.shuttle_right = !state.shuttle_right;
        }

        if 5<3 {
        let mut config_line = format!("");
        let mut pointer = 0;
        for i in (0..state.offset) {
            if pointer%(sector_size+1)==0 || pointer%(sector_size+1)==sector_size {
                config_line += &format!("| ");
            }
            pointer += 1;
            config_line += &format!("O ");
        }
        for i in (0..num_qubits) {
            if pointer%(sector_size+1)==0 || pointer%(sector_size+1)==sector_size {
                config_line += &format!("| ");
            }
            pointer += 1;
            config_line += &format!("q{} ", physical_to_logical[i]);
        }
        for i in (0..(2*sector_size*empty_sector)) {
            if pointer%(sector_size+1)==0 || pointer%(sector_size+1)==sector_size {
                config_line += &format!("| ");
            }
            pointer += 1;
            config_line += &format!("O ");
        }
        println!("{}", config_line);
        if start_time > 987654321000 {
            println!("NTCF {} TIME OUT", sector_size);
            break;
        }
    }
    }

    for qubit in 0..num_qubits { // Push qec for error generator
        sched.qubits[qubit].push("q".to_string());
        sched.timestamp[qubit].push(start_time);
    }

    (sched, start_time)
}

fn choose_ntcf(qubits: &Vec<usize>, physical_to_logical: &Vec<usize>, logical_to_physical: &Vec<usize>, current_nodes: &Vec<usize>, circuit: &Circuit, sector_size: usize) -> (usize, usize) {
    if qubits.len() == 1 { // Vertical sector
        return (qubits[0], 0);
    }
    for &qubit in qubits { // Choose T
        let logical_qubit = physical_to_logical[qubit];
        let node_index = current_nodes[logical_qubit];
        if node_index == circuit.qubits[logical_qubit].len() { // Finished qubit
            continue;
        }
        let node = circuit.qubits[logical_qubit][node_index];
        match node {
            b'r' => {
                //println!("For R q{}", logical_qubit);
                return (qubit, 1);}, // T gate
            b'g' | b'p' => {}
            b'm' => {}
            _ => unreachable!()
        }
    }
    let mut candidate = None;
    let shuttle_right = qubits[0] > qubits[1];
    for &qubit in qubits { // Choose close
        let logical_qubit = physical_to_logical[qubit];
        let node_index = current_nodes[logical_qubit];
        if node_index == circuit.qubits[logical_qubit].len() { // Finished qubit
            continue;
        }
        let node = circuit.qubits[logical_qubit][node_index];
        match node {
            b'm' => {
                let partner = circuit.cx[&(logical_qubit, node_index)];
                let physical_partner = logical_to_physical[partner.0];
                if current_nodes[partner.0] != partner.1 { // Partner qubit has preceding gates
                    continue;
                }
                if qubit.abs_diff(physical_partner) < sector_size { // Close gate
                    if physical_partner > qubit && shuttle_right || physical_partner < qubit && !shuttle_right {
                        //println!("For Close q{} q{}", logical_qubit, partner.0);
                        return (qubit, 2);
                    }
                } else if candidate.is_none() && shuttle_right && physical_partner >= qubit + sector_size { // Far gate
                    candidate = Some(qubit);
                    //println!("For Far q{} q{}", logical_qubit, partner.0);
                } else if candidate.is_none() && !shuttle_right && qubit >= physical_partner + sector_size { // Far gate
                    candidate = Some(qubit);
                    //println!("For Far q{} q{}", logical_qubit, partner.0);
                }
            }
            b'p' | b'g' => {}
            _ => unreachable!()
        }
    }
    if let Some(candidate) = candidate {
        (candidate, 3)
    } else {
        (qubits[0], 0)
    }
}

fn initialize_qubit_pos(num_qubits: usize, circuit: &Circuit) -> (Vec<Logical>, Vec<Physical>) {
    let mut mapping = vec![0; num_qubits];
    for i in 0..num_qubits {
        mapping[i] = i;
    }
    (mapping.clone(), mapping)
}

pub struct Schedule {
    pub qubits: Vec<Vec<String>>,
    pub timestamp: Vec<Vec<usize>>,
    pub cx: HashMap<(usize, usize), (usize, usize)>,
}