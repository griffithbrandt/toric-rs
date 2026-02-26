/// toric_rs - High-performance QEC decoder benchmark suite
///
/// Usage:
///   toric_rs simulate  -d 5 -p 0.06
///   toric_rs benchmark -d 3,5,7,9 -t 5000
///   toric_rs threshold -d 3,5,7 --pmin 0.02 --pmax 0.14 --pstep 0.01
///
/// Author: Griffith Brandt

use std::env;
use std::time::Instant;
use toric_rs::*;
use rand_pcg::Pcg64;
use rand::SeedableRng;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_usage();
        return;
    }

    match args[1].as_str() {
        "simulate" => cmd_simulate(&args[2..]),
        "benchmark" => cmd_benchmark(&args[2..]),
        "threshold" => cmd_threshold(&args[2..]),
        "help" | "--help" | "-h" => print_usage(),
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            print_usage();
        }
    }
}

fn print_usage() {
    println!("toric_rs - High-performance toric code QEC decoder\n");
    println!("COMMANDS:");
    println!("  simulate    Run a single decode cycle");
    println!("  benchmark   Compare Union-Find vs Greedy decoder speed");
    println!("  threshold   Run full threshold sweep and generate CSV\n");
    println!("SIMULATE OPTIONS:");
    println!("  -d <int>    Code distance (default: 5)");
    println!("  -p <float>  Physical error rate (default: 0.06)\n");
    println!("BENCHMARK OPTIONS:");
    println!("  -d <ints>   Comma-separated distances (default: 3,5,7,9)");
    println!("  -p <float>  Error rate for speed test (default: 0.05)");
    println!("  -t <int>    Trials per config (default: 5000)\n");
    println!("THRESHOLD OPTIONS:");
    println!("  -d <ints>   Comma-separated distances (default: 3,5,7)");
    println!("  --pmin      Min error rate (default: 0.02)");
    println!("  --pmax      Max error rate (default: 0.14)");
    println!("  --pstep     Error rate step (default: 0.01)");
    println!("  -t <int>    Trials per point (default: 5000)");
    println!("  -o <file>   Output CSV path (default: threshold_data.csv)\n");
    println!("EXAMPLES:");
    println!("  toric_rs simulate -d 7 -p 0.04");
    println!("  toric_rs benchmark -d 3,5,7,9,11 -t 10000");
    println!("  toric_rs threshold -d 3,5,7,9 -t 10000 -o results.csv");
}

// -------------------------------------------------------------------------
// Arg parsing helpers
// -------------------------------------------------------------------------

fn get_arg(args: &[String], flag: &str, default: &str) -> String {
    for i in 0..args.len() {
        if args[i] == flag && i + 1 < args.len() {
            return args[i + 1].clone();
        }
    }
    default.to_string()
}

fn parse_distances(s: &str) -> Vec<usize> {
    s.split(',').filter_map(|x| x.trim().parse().ok()).collect()
}

// -------------------------------------------------------------------------
// Commands
// -------------------------------------------------------------------------

fn cmd_simulate(args: &[String]) {
    let d: usize = get_arg(args, "-d", "5").parse().unwrap_or(5);
    let p: f64 = get_arg(args, "-p", "0.06").parse().unwrap_or(0.06);
    let seed: u64 = get_arg(args, "--seed", "42").parse().unwrap_or(42);

    let code = ToricCode::new(d);
    let uf = UnionFindDecoder::new(&code);
    let mut rng = Pcg64::seed_from_u64(seed);

    let z_errors = bitflip_noise(code.n_qubits, p, &mut rng);
    let syndrome = code.x_syndrome(&z_errors);

    let n_errors: usize = z_errors.iter().map(|&x| x as usize).sum();
    let n_defects: usize = syndrome.iter().map(|&x| x as usize).sum();

    let start = Instant::now();
    let correction = uf.decode(&syndrome);
    let decode_time = start.elapsed();

    let n_corr: usize = correction.iter().map(|&x| x as usize).sum();
    let residual: Vec<u8> = z_errors.iter()
        .zip(correction.iter())
        .map(|(&e, &c)| e ^ c)
        .collect();
    let n_res: usize = residual.iter().map(|&x| x as usize).sum();
    let logical = code.has_logical_error(&residual);

    println!();
    println!("========================================");
    println!("  TORIC-RS SINGLE DECODE CYCLE");
    println!("========================================");
    println!();
    println!("  Code:            Toric d={}", d);
    println!("  Data qubits:     {}", code.n_qubits);
    println!("  Stabilizers:     {}", code.n_stabs);
    println!("  Error rate:      {:.3}", p);
    println!();
    println!("  Errors injected: {} / {} qubits", n_errors, code.n_qubits);
    println!("  Syndrome weight: {} / {} stabilizers", n_defects, code.n_stabs);
    println!("  Correction wt:   {} qubits", n_corr);
    println!("  Residual weight: {}", n_res);
    println!("  Decode time:     {:.1} us", decode_time.as_nanos() as f64 / 1000.0);
    println!();
    if logical {
        println!("  RESULT: LOGICAL ERROR (decode failed)");
    } else {
        println!("  RESULT: SUCCESSFULLY CORRECTED");
    }
    println!();
}

fn cmd_benchmark(args: &[String]) {
    let distances = parse_distances(&get_arg(args, "-d", "3,5,7,9"));
    let p: f64 = get_arg(args, "-p", "0.05").parse().unwrap_or(0.05);
    let trials: usize = get_arg(args, "-t", "5000").parse().unwrap_or(5000);

    println!();
    println!("========================================");
    println!("  TORIC-RS DECODER BENCHMARK");
    println!("========================================");
    println!();
    println!("  Distances:  {:?}", distances);
    println!("  Error rate: {:.3}", p);
    println!("  Trials:     {}", trials);
    println!();

    // Header
    println!(
        "  {:>5}  {:>8}  {:>14}  {:>10}  {:>14}  {:>10}",
        "d", "n_qubits", "UF p_L", "UF time", "Greedy p_L", "Greedy time"
    );
    println!("  {}", "-".repeat(70));

    for &d in &distances {
        let uf_result = run_benchmark(d, p, trials, "union-find", 42);
        let gr_result = run_benchmark(d, p, trials, "greedy", 42);

        let uf_time = if uf_result.mean_decode_ns > 1_000_000.0 {
            format!("{:.2} ms", uf_result.mean_decode_ns / 1_000_000.0)
        } else {
            format!("{:.1} us", uf_result.mean_decode_ns / 1_000.0)
        };
        let gr_time = if gr_result.mean_decode_ns > 1_000_000.0 {
            format!("{:.2} ms", gr_result.mean_decode_ns / 1_000_000.0)
        } else {
            format!("{:.1} us", gr_result.mean_decode_ns / 1_000.0)
        };

        let code = ToricCode::new(d);
        println!(
            "  {:>5}  {:>8}  {:>14.5}  {:>10}  {:>14.5}  {:>10}",
            d, code.n_qubits,
            uf_result.logical_error_rate, uf_time,
            gr_result.logical_error_rate, gr_time,
        );
    }

    println!();

    // Speedup summary
    println!("  SPEEDUP SUMMARY (Union-Find vs Greedy):");
    println!("  {}", "-".repeat(45));
    for &d in &distances {
        let uf = run_benchmark(d, p, trials, "union-find", 42);
        let gr = run_benchmark(d, p, trials, "greedy", 42);
        let speedup = gr.mean_decode_ns / uf.mean_decode_ns.max(1.0);
        let accuracy_diff = gr.logical_error_rate - uf.logical_error_rate;
        println!(
            "  d={:>2}: UF is {:.1}x faster, {:.3} better p_L",
            d, speedup, accuracy_diff.abs()
        );
    }
    println!();
}

fn cmd_threshold(args: &[String]) {
    let distances = parse_distances(&get_arg(args, "-d", "3,5,7"));
    let p_min: f64 = get_arg(args, "--pmin", "0.02").parse().unwrap_or(0.02);
    let p_max: f64 = get_arg(args, "--pmax", "0.14").parse().unwrap_or(0.14);
    let p_step: f64 = get_arg(args, "--pstep", "0.01").parse().unwrap_or(0.01);
    let trials: usize = get_arg(args, "-t", "5000").parse().unwrap_or(5000);
    let output = get_arg(args, "-o", "threshold_data.csv");

    // Build error rate list
    let mut error_rates = Vec::new();
    let mut p = p_min;
    while p <= p_max + 1e-9 {
        error_rates.push((p * 10000.0).round() / 10000.0);
        p += p_step;
    }

    let total = distances.len() * error_rates.len();

    println!();
    println!("========================================");
    println!("  TORIC-RS THRESHOLD SWEEP");
    println!("========================================");
    println!();
    println!("  Distances:   {:?}", distances);
    println!("  Error rates: {} points ({:.3} to {:.3})", error_rates.len(), p_min, p_max);
    println!("  Trials:      {} per point", trials);
    println!("  Total runs:  {}", total);
    println!("  Decoder:     Union-Find");
    println!();

    let mut csv_lines = vec!["distance,error_rate,logical_error_rate,n_trials,n_logical_errors,mean_decode_ns".to_string()];
    let mut count = 0;

    let sweep_start = Instant::now();

    for &d in &distances {
        for &er in &error_rates {
            count += 1;
            let seed = 42 + d as u64 * 1000 + (er * 10000.0) as u64;
            let result = run_benchmark(d, er, trials, "union-find", seed);

            let marker = if result.logical_error_rate < er { "+" } else { "!" };

            println!(
                "  [{:>3}/{}] d={} p={:.3} -> p_L={:.5} ({}/{}) [{:.0} ns/decode] {}",
                count, total, d, er,
                result.logical_error_rate,
                result.n_logical_errors, result.n_trials,
                result.mean_decode_ns,
                marker,
            );

            csv_lines.push(format!(
                "{},{:.4},{:.6},{},{},{:.1}",
                d, er, result.logical_error_rate,
                result.n_trials, result.n_logical_errors,
                result.mean_decode_ns,
            ));
        }
    }

    let total_time = sweep_start.elapsed();

    // Write CSV
    std::fs::write(&output, csv_lines.join("\n")).expect("Failed to write CSV");

    println!();
    println!("  Sweep completed in {:.1}s", total_time.as_secs_f64());
    println!("  Results saved to {}", output);

    // Print summary table
    println!();
    println!("  THRESHOLD SUMMARY:");
    println!("  {}", "-".repeat(55));
    println!(
        "  {:>5}  {:>10}  {:>10}  {:>10}  {:>10}",
        "d", "p=0.03", "p=0.06", "p=0.09", "p=0.12"
    );
    println!("  {}", "-".repeat(55));
    for &d in &distances {
        let mut row = format!("  {:>5}", d);
        for target_p in [0.03, 0.06, 0.09, 0.12] {
            let seed = 42 + d as u64 * 1000 + (target_p * 10000.0) as u64;
            // Find closest result
            let closest = error_rates.iter()
                .min_by(|a, b| {
                    ((**a - target_p).abs()).partial_cmp(&((**b - target_p).abs())).unwrap()
                })
                .unwrap();
            if (*closest - target_p).abs() < 0.005 {
                let r = run_benchmark(d, *closest, trials, "union-find", seed);
                row.push_str(&format!("  {:>10.5}", r.logical_error_rate));
            } else {
                row.push_str("       N/A");
            }
        }
        println!("{}", row);
    }
    println!();
}
