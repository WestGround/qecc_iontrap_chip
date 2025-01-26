# Code for [Title]

This repository contains a Rust implementation of a noisy quantum circuit simulator. It is designed to support circuit optimization, scheduling, and simulation.

---

## Features

### 1. Compilation
- **Gate Reduction**: Optimizes circuits by focusing on reducing `Rz` gates.
  - Implementation in the `optimize` module.
- **Hardware Native Conversion**: Converts circuits into hardware-native gates.
  - Implementation in the `compile` module.

### 2. Scheduling
- **Schedule Generation**: Produces execution schedules based on an optimized circuit and hardware specifications (sector size, quantum error correction used, etc.).
  - Implementation in the `schedule` module.

### 3. Circuit Simulation
- **Success Probability Calculation**: Computes the success probability of a given schedule. Note that results may vary across runs due to noise simulation.
  - Implementation in the `error_generator` module.

---

## Usage

```bash
($FLAGS) cargo run --release - < ($INPUT) > ($OUTPUT)
```
- **Flags**
  - RATE_EXP: Schedule with varying error rate.
  - SIZE_EXP: Schedule with varying sector size.
  - QEC: Enable quantum error correction.
  - NO_QEC: Disable quantum error correction.
  
Additional flags may be supported depending on specific needs.
