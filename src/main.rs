mod circuit;
mod schedule;
mod error_generator;

use std::io;
use std::io::Read;
use crate::circuit::Circuit;
use std::env;
use qecc::*;
use qecc::schedule::ntcf_scheduler;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

fn main() {
    let mut buf = vec![];
    io::stdin().lock().read_to_end(&mut buf).unwrap();

    let mut circuit = Circuit::from_str_native(std::str::from_utf8(&buf).unwrap());


    let empty_sector = 3;
    let mut default_sector_size = 1;
    let default_error_rate = 1e-4;

    let mut num_iter = 0;

    let indiv_experiment = env::var("INDIV_EXP").is_ok();
    let error_rate_experiment = env::var("RATE_EXP").is_ok();
    let sector_size_experiment = env::var("SIZE_EXP").is_ok();
    let sched_experiment = env::var("SCHED_EXP").is_ok();

    let qec_used = env::var("QEC").is_ok();
    let qec_not_used = env::var("NO_QEC").is_ok();
    // Error rate experiment
    if error_rate_experiment && qec_used {
        default_sector_size = env::var("TRAP_SIZE")
            .expect("environment variable TRAP_SIZE should be set")
            .parse::<usize>()
            .expect("Fail to parse TRAP_SIZE into integer");
        num_iter = env::var("NUM_ITER")
            .expect("environment variable NUM_ITER should be set")
            .parse::<usize>()
            .expect("Fail to parse NUM_ITER into integer");
        error_rate_exp_qec(&circuit, default_error_rate, empty_sector, default_sector_size, num_iter);
    } else if error_rate_experiment && qec_not_used {
        default_sector_size = env::var("TRAP_SIZE")
            .expect("environment variable TRAP_SIZE should be set")
            .parse::<usize>()
            .expect("Fail to parse TRAP_SIZE into integer");
        num_iter = env::var("NUM_ITER")
            .expect("environment variable NUM_ITER should be set")
            .parse::<usize>()
            .expect("Fail to parse NUM_ITER into integer");
        error_rate_exp_noqec(&circuit, default_error_rate, empty_sector, default_sector_size, num_iter);
    }

    // Sector size experiment
    else if sector_size_experiment && qec_used {
        num_iter = env::var("NUM_ITER")
            .expect("environment variable NUM_ITER should be set")
            .parse::<usize>()
            .expect("Fail to parse NUM_ITER into integer");
        sector_size_exp_qec_thread(&circuit, default_error_rate, empty_sector, num_iter);
    } else if sector_size_experiment && qec_not_used {
        num_iter = env::var("NUM_ITER")
            .expect("environment variable NUM_ITER should be set")
            .parse::<usize>()
            .expect("Fail to parse NUM_ITER into integer");
        sector_size_exp_noqec_thread(&circuit, default_error_rate, empty_sector, num_iter);
    }

    // Individual experiment (for num_qubit experiment)
    else if indiv_experiment && qec_used {
        default_sector_size = env::var("TRAP_SIZE")
            .expect("environment variable TRAP_SIZE should be set")
            .parse::<usize>()
            .expect("Fail to parse TRAP_SIZE into integer");
        num_iter = env::var("NUM_ITER")
            .expect("environment variable NUM_ITER should be set")
            .parse::<usize>()
            .expect("Fail to parse NUM_ITER into integer");
        indiv_exp_qec(&circuit, default_error_rate, default_sector_size, empty_sector, num_iter);
    } else if indiv_experiment && qec_not_used {
        default_sector_size = env::var("TRAP_SIZE")
            .expect("environment variable TRAP_SIZE should be set")
            .parse::<usize>()
            .expect("Fail to parse TRAP_SIZE into integer");
        num_iter = env::var("NUM_ITER")
            .expect("environment variable NUM_ITER should be set")
            .parse::<usize>()
            .expect("Fail to parse NUM_ITER into integer");
        indiv_exp_noqec(&circuit, default_error_rate, default_sector_size, empty_sector, num_iter);
    }

    // Sched experiment
    else if sched_experiment {
        default_sector_size = env::var("TRAP_SIZE")
            .expect("environment variable TRAP_SIZE should be set")
            .parse::<usize>()
            .expect("Fail to parse TRAP_SIZE into integer");
        let sched = schedule::ntcf_scheduler(&circuit, default_sector_size, empty_sector, C17_QECTIME);
        //let sched = schedule::pmark_scheduler(&no_qec_circuit, default_sector_size, empty_sector);
        for (index, qubit) in sched.0.qubits.iter().enumerate() {
            print!("q{} ", index);
            for gate in qubit {
                print!("{} ", String::from(gate));
            }
            println!();
        }
    }

    else {
        unreachable!();
    }
}

fn error_rate_exp_qec(circuit: &Circuit, error_rate: f64, empty_sector: usize, default_sector_size: usize, num_iter: usize) {
    let num_tick = 26;
    let error_rates: Vec<_> = (0..=num_tick).map(|i| i as f64 * 5e-6).collect(); // (2e-5, 4e-5, ..., 5e-4)
    let mut mark_rate = vec![vec![Vec::new(); num_iter]; 3];
    let mut mark_runtime = vec![vec![Vec::new(); num_iter]; 3];
    let mut ntcf_rate = vec![vec![Vec::new(); num_iter]; 3];
    let mut ntcf_runtime = vec![vec![Vec::new(); num_iter]; 3];

    for code_index in 0..3 {
        for iter in 0..num_iter {
            let (ntcf_sched, runtime) = schedule::ntcf_scheduler(&circuit, default_sector_size, empty_sector, CODE_TIMES[code_index]);
            ntcf_runtime[code_index][iter].push(runtime*circuit.approx_factor);
            let (mark_sched, runtime) = schedule::mark_scheduler(&circuit, default_sector_size, empty_sector, CODE_TIMES[code_index]);
            mark_runtime[code_index][iter].push(runtime*circuit.approx_factor);
            for &error_rate in error_rates.iter() {
                let mark_rate_val = error_generator::generate_error_depol(&mark_sched, CODES[code_index], error_rate);
                let ntcf_rate_val = error_generator::generate_error_depol(&ntcf_sched, CODES[code_index], error_rate);
                mark_rate[code_index][iter].push(mark_rate_val.powi(circuit.approx_factor as i32));
                ntcf_rate[code_index][iter].push(ntcf_rate_val.powi(circuit.approx_factor as i32));
            }
        }
    }

    // Print data
    println!("{}", num_iter);
    let mut x_axis_line = format!("error_rate");
    for i in 0..=num_tick {
        x_axis_line += &format!("; {}", error_rates[i]);
    }
    println!("{}", x_axis_line);

    for code_index in 0..3 {
        for iter in 0..num_iter {
            let mut mark_line = format!("Mark {:?}_{:?}", CODES[code_index], iter);
            for i in 0..=num_tick {
                mark_line += &format!("; {}", mark_rate[code_index][iter][i]);
            }
            // mark_line += &format!("; runtime {:?}", mark_runtime[code_index]);
            println!("{}", mark_line);
        }
    }

    for code_index in 0..3 {
        for iter in 0..num_iter {
            let mut ntcf_line = format!("NTCF {:?}_{:?}", CODES[code_index], iter);
            for i in 0..=num_tick {
                ntcf_line += &format!("; {}", ntcf_rate[code_index][iter][i]);
            }
            // ntcf_line += &format!("; runtime {:?}", ntcf_runtime[code_index]);
            println!("{}", ntcf_line);
        }
    }
}

fn error_rate_exp_noqec(no_qec_circuit: &Circuit, error_rate: f64, empty_sector: usize, default_sector_size: usize, num_iter: usize) {
    let num_tick = 26;
    let error_rates: Vec<_> = (0..=num_tick).map(|i| i as f64 * 5e-6).collect(); // (2e-5, 4e-5, ..., 5e-4)
    let mut no_qec_rate = vec![vec![]; num_iter];
    let mut no_qec_runtime = vec![vec![]; num_iter];
    let mut long_rate = vec![vec![]; num_iter];
    let mut long_runtime = vec![vec![]; num_iter];

    for iter in 0..num_iter {
        let (no_qec_sched, runtime) = schedule::pmark_scheduler(&no_qec_circuit, default_sector_size, empty_sector);
        no_qec_runtime[iter].push(runtime*no_qec_circuit.approx_factor);
        for &error_rate in error_rates.iter() {
            let no_qec_rate_val = error_generator::generate_error_depol(&no_qec_sched, (1, 1, 1), error_rate);
            no_qec_rate[iter].push(no_qec_rate_val.powi(no_qec_circuit.approx_factor as i32));
        }

        let (long_sched, runtime) = schedule::pmark_scheduler(&no_qec_circuit, 100, empty_sector);
        long_runtime[iter].push(runtime*no_qec_circuit.approx_factor);
        for &error_rate in error_rates.iter() {
            let long_rate_val = error_generator::generate_error_depol(&long_sched, (1, 1, 1), error_rate);
            long_rate[iter].push(long_rate_val.powi(no_qec_circuit.approx_factor as i32));
        }
    }

    // Print data
    println!("{}", num_iter);
    let mut x_axis_line = format!("error_rate");
    for i in 0..=num_tick {
        x_axis_line += &format!("; {}", error_rates[i]);
    }
    println!("{}", x_axis_line);

    for iter in 0..num_iter {
        let mut no_qec_line = format!("No QEC (short chain)_{:?}", iter);
        for i in 0..=num_tick {
            no_qec_line += &format!("; {}", no_qec_rate[iter][i]);
        }
        //no_qec_line += &format!("; {:?}", no_qec_runtime);
        println!("{}", no_qec_line);

        let mut no_qec_long_line = format!("No QEC (long chain)_{:?}", iter);
        for i in 0..=num_tick {
            no_qec_long_line += &format!("; {}", long_rate[iter][i]);
        }
        //no_qec_long_line += &format!("; {:?}", long_runtime);
        println!("{}", no_qec_long_line);
    }
}

fn sector_size_exp_qec_thread(circuit: &Circuit, error_rate: f64, empty_sector: usize, num_iter: usize) {
    let sector_sizes: Vec<_> = (2..5).collect();
    
    // Wrapping data structures in Arc<Mutex<>> for shared mutable access
    let mark_rate = Arc::new(Mutex::new(vec![vec![Vec::new(); num_iter]; 3]));
    let mark_runtime = Arc::new(Mutex::new(vec![vec![Vec::new(); num_iter]; 3]));
    let ntcf_rate = Arc::new(Mutex::new(vec![vec![Vec::new(); num_iter]; 3]));
    let ntcf_runtime = Arc::new(Mutex::new(vec![vec![Vec::new(); num_iter]; 3]));

    // Outer loop parallelized over iter
    (0..num_iter).into_par_iter().for_each(|iter| {
        for code_index in 0..3 {
            for &sector_size in sector_sizes.iter() {
                let (ntcf_sched, runtime) = schedule::ntcf_scheduler(&circuit, sector_size, empty_sector, CODE_TIMES[code_index]);
                {
                    let mut ntcf_runtime_lock = ntcf_runtime.lock().unwrap();
                    ntcf_runtime_lock[code_index][iter].push(runtime * circuit.approx_factor);
                }

                let (mark_sched, runtime) = schedule::mark_scheduler(&circuit, sector_size, empty_sector, CODE_TIMES[code_index]);
                {
                    let mut mark_runtime_lock = mark_runtime.lock().unwrap();
                    mark_runtime_lock[code_index][iter].push(runtime * circuit.approx_factor);
                }

                let mark_rate_val = error_generator::generate_error_depol(&mark_sched, CODES[code_index], error_rate);
                let ntcf_rate_val = error_generator::generate_error_depol(&ntcf_sched, CODES[code_index], error_rate);
                {
                    let mut mark_rate_lock = mark_rate.lock().unwrap();
                    mark_rate_lock[code_index][iter].push(mark_rate_val.powi(circuit.approx_factor as i32));

                    let mut ntcf_rate_lock = ntcf_rate.lock().unwrap();
                    ntcf_rate_lock[code_index][iter].push(ntcf_rate_val.powi(circuit.approx_factor as i32));
                }
            }
        }
    });

    // Print data (same as before)
    println!("{}", num_iter);
    let mut x_axis_line = format!("sector_size");
    for &sector_size in sector_sizes.iter() {
        x_axis_line += &format!("; {}", sector_size);
    }
    println!("{}", x_axis_line);

    // Printing the results in the same manner
    for code_index in 0..3 {
        for iter in 0..num_iter {
            let mark_rate_lock = mark_rate.lock().unwrap();
            let mut mark_rate_line = format!("Mark {:?}_{:?}", CODES[code_index], iter);
            for i in 0..mark_rate_lock[code_index][iter].len() {
                mark_rate_line += &format!("; {}", mark_rate_lock[code_index][iter][i]);
            }
            println!("{}", mark_rate_line);
        }
    }

    for code_index in 0..3 {
        for iter in 0..num_iter {
            let ntcf_rate_lock = ntcf_rate.lock().unwrap();
            let mut ntcf_rate_line = format!("NTCF {:?}_{:?}", CODES[code_index], iter);
            for i in 0..ntcf_rate_lock[code_index][iter].len() {
                ntcf_rate_line += &format!("; {}", ntcf_rate_lock[code_index][iter][i]);
            }
            println!("{}", ntcf_rate_line);
        }
    }

    println!("");

    // Same print block for runtime values
    for code_index in 0..3 {
        for iter in 0..num_iter {
            let mark_runtime_lock = mark_runtime.lock().unwrap();
            let mut mark_runtime_line = format!("Mark {:?}_{:?}", CODES[code_index], iter);
            for i in 0..mark_runtime_lock[code_index][iter].len() {
                mark_runtime_line += &format!("; {}", mark_runtime_lock[code_index][iter][i]);
            }
            println!("{}", mark_runtime_line);
        }
    }

    for code_index in 0..3 {
        for iter in 0..num_iter {
            let ntcf_runtime_lock = ntcf_runtime.lock().unwrap();
            let mut ntcf_runtime_line = format!("NTCF {:?}_{:?}", CODES[code_index], iter);
            for i in 0..ntcf_runtime_lock[code_index][iter].len() {
                ntcf_runtime_line += &format!("; {}", ntcf_runtime_lock[code_index][iter][i]);
            }
            println!("{}", ntcf_runtime_line);
        }
    }
}

fn sector_size_exp_noqec_thread(no_qec_circuit: &Circuit, error_rate: f64, empty_sector: usize, num_iter: usize) {
    let sector_sizes: Vec<_> = (2..5).collect();

    // Wrapping data structures in Arc<Mutex<>> for shared mutable access
    let no_qec_rate = Arc::new(Mutex::new(vec![Vec::new(); num_iter]));
    let no_qec_runtime = Arc::new(Mutex::new(vec![Vec::new(); num_iter]));

    // Outer loop parallelized over iter
    (0..num_iter).into_par_iter().for_each(|iter| {
        for &sector_size in sector_sizes.iter() {
            // Calculate pmark_scheduler
            let (no_qec_sched, runtime) = schedule::pmark_scheduler(&no_qec_circuit, sector_size, empty_sector);

            // Lock and push runtime data for no_qec_runtime
            {
                let mut no_qec_runtime_lock = no_qec_runtime.lock().unwrap();
                no_qec_runtime_lock[iter].push(runtime * no_qec_circuit.approx_factor);
            }

            // Error generation for no_qec
            let no_qec_rate_val = error_generator::generate_error_depol(&no_qec_sched, (1, 1, 1), error_rate);

            // Lock and push rate data for no_qec_rate
            {
                let mut no_qec_rate_lock = no_qec_rate.lock().unwrap();
                no_qec_rate_lock[iter].push(no_qec_rate_val.powi(no_qec_circuit.approx_factor as i32));
            }
        }
    });

    // Print data (same as before)
    println!("{}", num_iter);
    let mut x_axis_line = format!("sector_size");
    for &sector_size in sector_sizes.iter() {
        x_axis_line += &format!("; {}", sector_size);
    }
    println!("{}", x_axis_line);

    // Print no_qec_rate data
    for iter in 0..num_iter {
        let no_qec_rate_lock = no_qec_rate.lock().unwrap();
        let mut noqec_rate_line = format!("No QEC (short chain)_{:?}", iter);
        for i in 0..no_qec_rate_lock[iter].len() {
            noqec_rate_line += &format!("; {:?}", no_qec_rate_lock[iter][i]);
        }
        println!("{}", noqec_rate_line);
    }

    // Print no_qec_runtime data
    for iter in 0..num_iter {
        let no_qec_runtime_lock = no_qec_runtime.lock().unwrap();
        let mut noqec_runtime_line = format!("No QEC (short chain)_{:?}", iter);
        for i in 0..no_qec_runtime_lock[iter].len() {
            noqec_runtime_line += &format!("; {:?}", no_qec_runtime_lock[iter][i]);
        }
        println!("{}", noqec_runtime_line);
    }
}

fn indiv_exp_qec(circuit: &Circuit, error_rate: f64, sector_size: usize, empty_sector: usize, num_iter: usize) {
    println!("{}", num_iter);
    for code_index in 0..3 {
        let mut ntcf_rate_line = format!("NTCF {:?}", CODES[code_index]);
        let mut ntcf_time_line = format!("NTCF {:?}", CODES[code_index]);
        let mut mark_rate_line = format!("Mark {:?}", CODES[code_index]);
        let mut mark_time_line = format!("Mark {:?}", CODES[code_index]);
        for iter in 0..num_iter {
            let (ntcf_sched, runtime) = schedule::ntcf_scheduler(&circuit, sector_size, empty_sector, CODE_TIMES[code_index]);
            let rate = error_generator::generate_error_depol(&ntcf_sched, CODES[code_index], error_rate);
            ntcf_rate_line += &format!("; {}", rate.powi(circuit.approx_factor as i32));
            ntcf_time_line += &format!("; {}", runtime*(circuit.approx_factor));

            let (mark_sched, runtime) = schedule::mark_scheduler(&circuit, sector_size, empty_sector, CODE_TIMES[code_index]);
            let rate = error_generator::generate_error_depol(&mark_sched, CODES[code_index], error_rate);
            mark_rate_line += &format!("; {}", rate.powi(circuit.approx_factor as i32));
            mark_time_line += &format!("; {}", runtime*(circuit.approx_factor));
        }
        println!("{}\n{}", mark_rate_line, mark_time_line);
        println!("{}\n{}", ntcf_rate_line, ntcf_time_line);
    }
}

fn indiv_exp_noqec(no_qec_circuit: &Circuit, error_rate: f64, sector_size: usize, empty_sector: usize, num_iter: usize) {
    println!("{}", num_iter);
    let mut short_rate_line = format!("No QEC (short chain)");
    let mut short_time_line = format!("No QEC (short chain)");
    let mut long_rate_line = format!("No QEC (long chain)");
    let mut long_time_line = format!("No QEC (long chain)");
    for iter in 0..num_iter {
        let (no_qec_sched, runtime) = schedule::pmark_scheduler(&no_qec_circuit, sector_size, empty_sector);
        let rate = error_generator::generate_error_depol(&no_qec_sched, (1, 1, 1), error_rate);
        short_rate_line += &format!("; {}", rate.powi(no_qec_circuit.approx_factor as i32));
        short_time_line += &format!("; {}", runtime*(no_qec_circuit.approx_factor));

        let (long_sched, runtime) = schedule::pmark_scheduler(&no_qec_circuit, 100, empty_sector);
        let rate = error_generator::generate_error_depol(&long_sched, (1, 1, 1), error_rate);
        long_rate_line += &format!("; {}", rate.powi(no_qec_circuit.approx_factor as i32));
        long_time_line += &format!("; {}", runtime*(no_qec_circuit.approx_factor));
    }
    println!("{}\n{}", short_rate_line, short_time_line);
    println!("{}\n{}", long_rate_line, long_time_line);
}

fn default_exp(circuit: &Circuit, no_qec_circuit: &Circuit, error_rate: f64, sector_size: usize, empty_sector: usize) {
    println!("No-QEC");
    let (no_qec_sched, runtime) = schedule::pmark_scheduler(&no_qec_circuit, sector_size, empty_sector);
    let rate = error_generator::generate_error_depol(&no_qec_sched, (1, 1, 1), error_rate);
    println!("{}", rate.powi(circuit.approx_factor as i32));
    println!("{}", runtime*circuit.approx_factor);

    let (long_sched, runtime) = schedule::pmark_scheduler(&no_qec_circuit, 100, empty_sector);
    let rate = error_generator::generate_error_depol(&long_sched, (1, 1, 1), error_rate);
    println!("No-QEC(long chain)");
    println!("{}", rate.powi(circuit.approx_factor as i32));
    println!("{}", runtime*circuit.approx_factor);
    for code_index in 0..3 {
        let (ntcf_sched, runtime) = schedule::ntcf_scheduler(&circuit, sector_size, empty_sector, CODE_TIMES[code_index]);
        let rate = error_generator::generate_error_depol(&ntcf_sched, CODES[code_index], error_rate);
        println!("NTCF {:?}", CODES[code_index]);
        println!("{}", rate.powi(circuit.approx_factor as i32));
        println!("{}", runtime*circuit.approx_factor);
        let (mark_sched, runtime) = schedule::mark_scheduler(&circuit, sector_size, empty_sector, CODE_TIMES[code_index]);
        let rate = error_generator::generate_error_depol(&mark_sched, CODES[code_index], error_rate);
        // println!("{}", std::str::from_utf8(&mark_sched.qubits[0]).unwrap());
        println!("Mark {:?}", CODES[code_index]);
        println!("{}", rate.powi(circuit.approx_factor as i32));
        println!("{}", runtime*circuit.approx_factor);
    }
}
