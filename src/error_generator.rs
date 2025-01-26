use std::collections::HashMap;
use std::mem::swap;
use crate::{C17_QECTIME, C7_CODE, C7_QECTIME};
use crate::{SQ_ERR, TQ_ERR, SWAP_ERR, SHUTTLE_ERR, MEASURE_ERR, QEC_ERR};
use crate::{SHUTTLE_TIME, SQ_TIME, SWAP_TIME, TQ_TIME, MEASURE_TIME};
use crate::schedule::Schedule;

/// Compute success probability of given schedule
pub fn generate_error_depol(sched: &Schedule, code: (usize, usize, usize), two_qubit_error: f64) -> f64 {
    let (n, _, d) = code;
    let mut success_prob = 1.0;
    let mut binom = Binom::new();

    let one_qubit_error = two_qubit_error * SQ_ERR;
    let swap_error = two_qubit_error * SWAP_ERR;
    let shuttle_error = two_qubit_error * SHUTTLE_ERR;
    let measure_error = two_qubit_error * MEASURE_ERR;
    let qec_error = two_qubit_error * QEC_ERR;

    for qubit_index in 0..sched.qubits.len() {
        let qubit = &sched.qubits[qubit_index];

        let mut time_last_qec = if n == 7 { -1 * C7_QECTIME as i32 } else { -1 * C17_QECTIME as i32 };
        let mut error_prob = 0.0;

        for gate_index in 0..qubit.len() {
            let gate = &qubit[gate_index];
            match gate.as_str() {
                "g" | "p" | "r" => { // GPI, GPI2, rz
                    error_prob = error_prob + one_qubit_error - 4.0*one_qubit_error*error_prob/3.0;
                }
                "m" => { // MS
                    error_prob = error_prob + two_qubit_error*4.0/5.0 - 16.0*two_qubit_error*error_prob/15.0;
                }
                "s" => { // Swap
                    error_prob = error_prob + swap_error*4.0/5.0 - 16.0*swap_error*error_prob/15.0;
                }
                "h" => { // Shuttle
                    error_prob = error_prob + shuttle_error - 4.0*shuttle_error*error_prob/3.0;
                }
                "M" => {
                    error_prob = error_prob + measure_error - 4.0*measure_error*error_prob/3.0;
                }
                "q" => { // QEC
                    let mut correctable_prob = 0.0;
                    if error_prob > 1.0 {
                        panic!("Error: error_prob is less than 0.0, which is invalid.");
                    }

                    let not_flip_prob = 1.0-(2.0*error_prob/3.0);
                    for num_errors in 0..=d / 2 {
                        correctable_prob += binom.compute(n, num_errors) * not_flip_prob.powi((n - num_errors) as i32) * (1.0 - not_flip_prob).powi(num_errors as i32);
                    }
                    success_prob *= (correctable_prob*correctable_prob);

                    error_prob = qec_error;
                }
                _ => unreachable!()
            }
            if gate_index==(qubit.len()-1) {
                assert_eq!(*gate, "q".to_string(), "The variable 'gate' is not exactly 'q'");
            }
        }
    }
    success_prob
}

pub fn overhead_cnt(sched: &Schedule, code: (usize, usize, usize), two_qubit_error: f64, title: String) -> (usize, usize, usize, f64, f64, f64, f64, f64, f64) {
    let (n, _, d) = code;
    let mut success_prob = 1.0;
    let mut binom = Binom::new();

    let one_qubit_error = two_qubit_error * SQ_ERR;
    let swap_error = two_qubit_error * SWAP_ERR;
    let shuttle_error = two_qubit_error * SHUTTLE_ERR;
    let measure_error = two_qubit_error * MEASURE_ERR;
    let qec_error = two_qubit_error * QEC_ERR;

    let mut swap_count = 0;
    let mut shuttle_count = 0;
    let mut qec_count = 0;

    let mut interqec_count = 0;
    let mut interqec_count_mean = 0.0;
    let mut interqec_count_std = 0.0;
    let mut interqec_counts = Vec::new();

    let mut interqec_swapcount = 0;
    let mut interqec_swapcount_mean = 0.0;
    let mut interqec_swapcount_std = 0.0;
    let mut interqec_swapcounts = Vec::new();

    let mut interqec_shuttlecount = 0;
    let mut interqec_shuttlecount_mean = 0.0;
    let mut interqec_shuttlecount_std = 0.0;
    let mut interqec_shuttlecounts = Vec::new();

    println!("Start of {}", title);
    for qubit_index in 0..sched.qubits.len() {
        assert_eq!(interqec_count, 0, "count err");
        assert_eq!(interqec_shuttlecount, 0, "shuttle count err");
        assert_eq!(interqec_swapcount, 0, "swap count err");
        let qubit = &sched.qubits[qubit_index];

        let mut time_last_qec = if n == 7 { -1 * C7_QECTIME as i32 } else { -1 * C17_QECTIME as i32 };
        let mut error_prob = 0.0;

        for gate_index in 0..qubit.len() {
            let gate = &qubit[gate_index];
            match gate.as_str() {
                "g" | "p" | "r" => { // GPI, GPI2, rz
                    error_prob = error_prob + one_qubit_error - 4.0*one_qubit_error*error_prob/3.0;
                    interqec_count += 1;
                }
                "m" => { // MS
                    error_prob = error_prob + two_qubit_error*4.0/5.0 - 16.0*two_qubit_error*error_prob/15.0;
                    interqec_count += 1;
                }
                "s" => { // Swap
                    error_prob = error_prob + swap_error*4.0/5.0 - 16.0*swap_error*error_prob/15.0;
                    swap_count += 1;
                    interqec_count += 1;
                    interqec_swapcount += 1;
                }
                "h" => { // Shuttle
                    error_prob = error_prob + shuttle_error - 4.0*shuttle_error*error_prob/3.0;
                    shuttle_count += 1;
                    interqec_count += 1;
                    interqec_shuttlecount += 1;
                }
                "M" => {
                    error_prob = error_prob + measure_error - 4.0*measure_error*error_prob/3.0;
                    interqec_count += 1;
                }
                "q" => {
                    let mut correctable_prob = 0.0;
                    if error_prob > 1.0 {
                        panic!("Error: error_prob is less than 0.0, which is invalid.");
                    }

                    let not_flip_prob = 1.0-(2.0*error_prob/3.0);
                    for num_errors in 0..=d / 2 {
                        correctable_prob += binom.compute(n, num_errors) * not_flip_prob.powi((n - num_errors) as i32) * (1.0 - not_flip_prob).powi(num_errors as i32);
                    }
                    success_prob *= (correctable_prob*correctable_prob);

                    error_prob = qec_error;

                    if 5<3{
                    if (correctable_prob * 10_000.0).round() / 10_000.0 != 1.000 {
                        println!("> {} {} {}", interqec_count, interqec_shuttlecount, interqec_swapcount);
                    } else if (correctable_prob * 100_000.0).round() / 100_000.0 != 1.000 {
                        println!(">>> {} {} {}", interqec_count, interqec_shuttlecount, interqec_swapcount);
                    }}

                    interqec_counts.push(interqec_count as f64);
                    interqec_swapcounts.push(interqec_swapcount as f64);
                    interqec_shuttlecounts.push(interqec_shuttlecount as f64);
                    qec_count += 1;
                    interqec_count = 0;
                    interqec_swapcount = 0;
                    interqec_shuttlecount = 0;
                }
                _ => unreachable!(),
            }
        }
    }
    // Calculate mean
    assert_eq!(interqec_counts.len(), qec_count, "final count err");
    assert_eq!(interqec_shuttlecounts.len(), qec_count, "finial shuttle count err");
    assert_eq!(interqec_swapcounts.len(), qec_count, "final swap count err");
    let interqec_count_mean = interqec_counts.iter().sum::<f64>() / interqec_counts.len() as f64;
    let interqec_shuttlecount_mean = interqec_shuttlecounts.iter().sum::<f64>() / interqec_shuttlecounts.len() as f64;
    let interqec_swapcount_mean = interqec_swapcounts.iter().sum::<f64>() / interqec_swapcounts.len() as f64;

    // Calculate standard deviation
    let variance = interqec_counts
        .iter()
        .map(|&count| (count - interqec_count_mean).powi(2))
        .sum::<f64>() / interqec_counts.len() as f64;
    let interqec_count_std = variance.sqrt();

    let variance = interqec_shuttlecounts
        .iter()
        .map(|&count| (count - interqec_shuttlecount_mean).powi(2))
        .sum::<f64>() / interqec_shuttlecounts.len() as f64;
    let interqec_shuttlecount_std = variance.sqrt();

    let variance = interqec_swapcounts
        .iter()
        .map(|&count| (count - interqec_swapcount_mean).powi(2))
        .sum::<f64>() / interqec_swapcounts.len() as f64;
    let interqec_swapcount_std = variance.sqrt();

    (
        shuttle_count,
        swap_count,
        qec_count,
        interqec_count_mean,
        interqec_count_std,
        interqec_shuttlecount_mean,
        interqec_shuttlecount_std,
        interqec_swapcount_mean,
        interqec_swapcount_std,
    )
}

struct Binom {
    memo: HashMap<(usize, usize), f64>,
}
impl Binom {
    fn new() -> Self {
        Self { memo: HashMap::default() }
    }
    fn compute(&mut self, n: usize, k: usize) -> f64 {
        if k == 0 || k == n {
            1.0
        } else if let Some(ret) = self.memo.get(&(n, k)) {
            *ret
        } else {
            let mut ret = self.compute(n - 1, k - 1);
            ret += self.compute(n - 1, k);
            self.memo.insert((n, k), ret);
            ret
        }
    }
}
