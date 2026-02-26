/// toric_rs - High-performance toric code QEC simulation
///
/// Implements the toric code, noise models, Union-Find decoder, and
/// Monte Carlo benchmarking for quantum error correction research.
///
/// Author: Griffith Brandt

use rand::Rng;
use rand_pcg::Pcg64;
use rand::SeedableRng;
use std::time::Instant;

// =========================================================================
// Toric Code
// =========================================================================

/// Distance-d toric code on a d x d lattice with periodic boundaries.
///
/// Qubits live on edges: 2*d^2 total (d^2 horizontal + d^2 vertical).
/// X-stabilizers at vertices, Z-stabilizers at faces.
pub struct ToricCode {
    pub d: usize,
    pub n_qubits: usize,
    pub n_stabs: usize,
    /// Each X-stabilizer touches 4 qubit indices
    x_stab_qubits: Vec<[usize; 4]>,
    /// For each qubit, which X-stabilizers it belongs to (at most 2)
    #[allow(dead_code)]
    qubit_x_stabs: Vec<[usize; 2]>,
}

impl ToricCode {
    pub fn new(d: usize) -> Self {
        assert!(d >= 2, "Code distance must be >= 2");
        let n_qubits = 2 * d * d;
        let n_stabs = d * d;

        let mut x_stab_qubits = Vec::with_capacity(n_stabs);
        for i in 0..d {
            for j in 0..d {
                // Star operator at vertex (i,j) touches 4 edges:
                // right h-edge, left h-edge, down v-edge, up v-edge
                x_stab_qubits.push([
                    h_index(d, i, j),
                    h_index(d, i, (j + d - 1) % d),
                    v_index(d, i, j),
                    v_index(d, (i + d - 1) % d, j),
                ]);
            }
        }

        // Build reverse map: for each qubit, which 2 X-stabs it touches
        let mut qubit_x_stabs = vec![[0usize; 2]; n_qubits];
        let mut qubit_stab_count = vec![0usize; n_qubits];
        for (s, qubits) in x_stab_qubits.iter().enumerate() {
            for &q in qubits {
                let idx = qubit_stab_count[q];
                assert!(idx < 2, "Qubit touches more than 2 X-stabs");
                qubit_x_stabs[q][idx] = s;
                qubit_stab_count[q] += 1;
            }
        }

        ToricCode { d, n_qubits, n_stabs, x_stab_qubits, qubit_x_stabs }
    }

    /// Compute X-stabilizer syndrome from Z errors.
    /// syndrome[s] = 1 if stabilizer s is triggered.
    pub fn x_syndrome(&self, z_errors: &[u8]) -> Vec<u8> {
        let mut syndrome = vec![0u8; self.n_stabs];
        for (s, qubits) in self.x_stab_qubits.iter().enumerate() {
            let mut parity = 0u8;
            for &q in qubits {
                parity ^= z_errors[q];
            }
            syndrome[s] = parity;
        }
        syndrome
    }

    /// Check if residual Z errors form a non-trivial homology cycle.
    /// Tests winding numbers on two independent transverse cuts.
    pub fn has_logical_error(&self, residual: &[u8]) -> bool {
        let d = self.d;
        // Horizontal winding: h-edges crossing vertical cut at col 0
        let mut h_wind = 0u8;
        for i in 0..d {
            h_wind ^= residual[h_index(d, i, 0)];
        }
        // Vertical winding: v-edges crossing horizontal cut at row 0
        let mut v_wind = 0u8;
        for j in 0..d {
            v_wind ^= residual[v_index(d, 0, j)];
        }
        h_wind != 0 || v_wind != 0
    }

    /// Get (row, col) coordinates of a stabilizer.
    pub fn stab_coords(&self, s: usize) -> (usize, usize) {
        (s / self.d, s % self.d)
    }

    /// Manhattan distance between two vertices on the torus.
    pub fn toric_dist(&self, a: usize, b: usize) -> usize {
        let (ai, aj) = self.stab_coords(a);
        let (bi, bj) = self.stab_coords(b);
        let d = self.d;
        let di = (ai as isize - bi as isize).unsigned_abs();
        let dj = (aj as isize - bj as isize).unsigned_abs();
        di.min(d - di) + dj.min(d - dj)
    }
}

#[inline]
fn h_index(d: usize, i: usize, j: usize) -> usize {
    (i % d) * d + (j % d)
}

#[inline]
fn v_index(d: usize, i: usize, j: usize) -> usize {
    d * d + (i % d) * d + (j % d)
}

// =========================================================================
// Noise Models
// =========================================================================

/// Apply independent bit-flip (Z) errors with probability p.
pub fn bitflip_noise(n: usize, p: f64, rng: &mut Pcg64) -> Vec<u8> {
    (0..n).map(|_| if rng.gen::<f64>() < p { 1 } else { 0 }).collect()
}

// =========================================================================
// Union-Find Decoder
// =========================================================================

/// Weighted Union-Find data structure with path compression.
struct UF {
    parent: Vec<usize>,
    rank: Vec<usize>,
    size: Vec<usize>,
}

impl UF {
    fn new(n: usize) -> Self {
        UF {
            parent: (0..n).collect(),
            rank: vec![0; n],
            size: vec![1; n],
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]]; // path halving
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) -> usize {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb { return ra; }
        // Union by rank
        let (big, small) = if self.rank[ra] >= self.rank[rb] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.parent[small] = big;
        self.size[big] += self.size[small];
        if self.rank[big] == self.rank[small] {
            self.rank[big] += 1;
        }
        big
    }
}

/// Union-Find decoder for toric codes (Delfosse-Nickerson style).
///
/// Algorithm:
/// 1. Each syndrome defect starts as its own cluster
/// 2. Grow all odd-parity clusters by one half-edge simultaneously
/// 3. When two clusters merge, update parity (XOR of defect counts)
/// 4. Stop when all clusters have even parity
/// 5. Within each cluster, pair defects and correct along shortest paths
///
/// The clustering phase runs in nearly O(n * alpha(n)) time.
pub struct UnionFindDecoder<'a> {
    code: &'a ToricCode,
}

impl<'a> UnionFindDecoder<'a> {
    pub fn new(code: &'a ToricCode) -> Self {
        UnionFindDecoder { code }
    }

    pub fn decode(&self, syndrome: &[u8]) -> Vec<u8> {
        let code = self.code;
        let d = code.d;
        let n_stabs = code.n_stabs;
        let mut correction = vec![0u8; code.n_qubits];

        // Collect defect locations
        let defects: Vec<usize> = (0..n_stabs)
            .filter(|&s| syndrome[s] != 0)
            .collect();

        if defects.is_empty() {
            return correction;
        }

        // Phase 1: Union-Find clustering
        // Grow odd-parity clusters until all are even
        let mut uf = UF::new(n_stabs);
        let mut cluster_parity = vec![0u8; n_stabs];
        for &s in &defects {
            cluster_parity[s] = 1;
        }

        // Build adjacency: stabilizer neighbors on the torus
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_stabs];
        for i in 0..d {
            for j in 0..d {
                let s = i * d + j;
                let right = i * d + (j + 1) % d;
                let down = ((i + 1) % d) * d + j;
                adj[s].push(right);
                adj[right].push(s);
                adj[s].push(down);
                adj[down].push(s);
            }
        }
        // Deduplicate
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }

        let mut growth_radius = vec![0usize; n_stabs];

        for radius in 1..=d {
            // Check if any cluster is still odd
            let mut has_odd = false;
            for &s in &defects {
                let root = uf.find(s);
                if cluster_parity[root] & 1 != 0 {
                    has_odd = true;
                    break;
                }
            }
            if !has_odd { break; }

            // For each edge between adjacent stabs, try to merge
            for i in 0..d {
                for j in 0..d {
                    let s = i * d + j;
                    let neighbors = [(i * d + (j + 1) % d), (((i + 1) % d) * d + j)];

                    for &nb in &neighbors {
                        let ra = uf.find(s);
                        let rb = uf.find(nb);
                        if ra == rb { continue; }

                        let a_odd = cluster_parity[ra] & 1 != 0;
                        let b_odd = cluster_parity[rb] & 1 != 0;

                        // Grow odd clusters
                        let a_reach = if a_odd { radius } else { growth_radius[ra] };
                        let b_reach = if b_odd { radius } else { growth_radius[rb] };

                        if a_reach + b_reach >= radius {
                            let pa = cluster_parity[ra];
                            let pb = cluster_parity[rb];
                            let new_root = uf.union(ra, rb);
                            cluster_parity[new_root] = pa ^ pb;
                            growth_radius[new_root] =
                                growth_radius[ra].max(growth_radius[rb]).max(radius);
                        }
                    }
                }
            }

            // Update growth radius for odd clusters
            for &s in &defects {
                let root = uf.find(s);
                if cluster_parity[root] & 1 != 0 {
                    growth_radius[root] = growth_radius[root].max(radius);
                }
            }
        }

        // Phase 2: Within each even-parity cluster, greedily pair defects
        // and correct along shortest paths
        let mut cluster_defects: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for &s in &defects {
            let root = uf.find(s);
            cluster_defects.entry(root).or_default().push(s);
        }

        for (_root, mut cdefs) in cluster_defects {
            // Greedily pair closest defects
            while cdefs.len() >= 2 {
                let mut best_dist = usize::MAX;
                let mut best_i = 0;
                let mut best_j = 1;
                for i in 0..cdefs.len() {
                    for j in (i + 1)..cdefs.len() {
                        let dist = code.toric_dist(cdefs[i], cdefs[j]);
                        if dist < best_dist {
                            best_dist = dist;
                            best_i = i;
                            best_j = j;
                        }
                    }
                }

                let a = cdefs[best_i];
                let b = cdefs[best_j];
                correct_path(code, &mut correction, a, b);

                cdefs.remove(best_j);
                cdefs.remove(best_i);
            }
        }

        correction
    }
}

/// Flip qubits along the shortest path between two stabilizers on the torus.
fn correct_path(code: &ToricCode, correction: &mut [u8], a: usize, b: usize) {
    let d = code.d;
    let (mut ci, cj_start) = code.stab_coords(a);
    let (bi, bj) = code.stab_coords(b);

    // Vertical moves
    let di = bi as isize - ci as isize;
    let step_i: isize = if di.unsigned_abs() > d / 2 {
        if di > 0 { -1 } else { 1 }
    } else if di > 0 { 1 } else if di < 0 { -1 } else { 0 };
    let vert_dist = if di.unsigned_abs() > d / 2 {
        d - di.unsigned_abs()
    } else {
        di.unsigned_abs()
    };

    for _ in 0..vert_dist {
        if step_i > 0 {
            correction[v_index(d, ci, cj_start)] ^= 1;
            ci = (ci + 1) % d;
        } else {
            ci = (ci + d - 1) % d;
            correction[v_index(d, ci, cj_start)] ^= 1;
        }
    }

    // Horizontal moves
    let mut cj = cj_start;
    let dj = bj as isize - cj as isize;
    let step_j: isize = if dj.unsigned_abs() > d / 2 {
        if dj > 0 { -1 } else { 1 }
    } else if dj > 0 { 1 } else if dj < 0 { -1 } else { 0 };
    let horiz_dist = if dj.unsigned_abs() > d / 2 {
        d - dj.unsigned_abs()
    } else {
        dj.unsigned_abs()
    };

    for _ in 0..horiz_dist {
        if step_j > 0 {
            correction[h_index(d, bi, cj)] ^= 1;
            cj = (cj + 1) % d;
        } else {
            cj = (cj + d - 1) % d;
            correction[h_index(d, bi, cj)] ^= 1;
        }
    }
}

// =========================================================================
// Greedy Decoder (baseline)
// =========================================================================

/// Simple nearest-neighbor greedy decoder for performance comparison.
pub struct GreedyDecoder<'a> {
    code: &'a ToricCode,
}

impl<'a> GreedyDecoder<'a> {
    pub fn new(code: &'a ToricCode) -> Self {
        GreedyDecoder { code }
    }

    pub fn decode(&self, syndrome: &[u8]) -> Vec<u8> {
        let code = self.code;
        let mut correction = vec![0u8; code.n_qubits];

        let mut defects: Vec<usize> = (0..code.n_stabs)
            .filter(|&s| syndrome[s] != 0)
            .collect();

        while defects.len() >= 2 {
            let mut best_dist = usize::MAX;
            let mut best_i = 0;
            let mut best_j = 1;
            for i in 0..defects.len() {
                for j in (i + 1)..defects.len() {
                    let dist = code.toric_dist(defects[i], defects[j]);
                    if dist < best_dist {
                        best_dist = dist;
                        best_i = i;
                        best_j = j;
                    }
                }
            }

            let a = defects[best_i];
            let b = defects[best_j];
            correct_path(code, &mut correction, a, b);

            defects.remove(best_j);
            defects.remove(best_i);
        }

        correction
    }
}

// =========================================================================
// Benchmarking
// =========================================================================

/// Results from a Monte Carlo QEC benchmark.
#[derive(Clone)]
pub struct BenchmarkResult {
    pub distance: usize,
    pub error_rate: f64,
    pub n_trials: usize,
    pub n_logical_errors: usize,
    pub logical_error_rate: f64,
    pub decoder_name: String,
    pub mean_decode_ns: f64,
}

/// Run a Monte Carlo benchmark for a given decoder.
pub fn run_benchmark(
    distance: usize,
    error_rate: f64,
    n_trials: usize,
    decoder: &str,
    seed: u64,
) -> BenchmarkResult {
    let code = ToricCode::new(distance);
    let mut rng = Pcg64::seed_from_u64(seed);

    let uf_dec = UnionFindDecoder::new(&code);
    let greedy_dec = GreedyDecoder::new(&code);

    let mut logical_errors = 0usize;
    let mut total_decode_ns = 0u128;

    for _ in 0..n_trials {
        let z_errors = bitflip_noise(code.n_qubits, error_rate, &mut rng);
        let syndrome = code.x_syndrome(&z_errors);

        let start = Instant::now();
        let correction = match decoder {
            "union-find" => uf_dec.decode(&syndrome),
            "greedy" => greedy_dec.decode(&syndrome),
            _ => panic!("Unknown decoder: {}", decoder),
        };
        total_decode_ns += start.elapsed().as_nanos();

        // Check for logical error
        let residual: Vec<u8> = z_errors.iter()
            .zip(correction.iter())
            .map(|(&e, &c)| e ^ c)
            .collect();

        if code.has_logical_error(&residual) {
            logical_errors += 1;
        }
    }

    BenchmarkResult {
        distance,
        error_rate,
        n_trials,
        n_logical_errors: logical_errors,
        logical_error_rate: logical_errors as f64 / n_trials as f64,
        decoder_name: decoder.to_string(),
        mean_decode_ns: total_decode_ns as f64 / n_trials as f64,
    }
}
