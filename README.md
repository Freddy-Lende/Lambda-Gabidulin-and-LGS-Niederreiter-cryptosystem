# Lambda-Gabidulin-and-LGS-Niederreiter-cryptosystem
This script computes the dimensions of the left and right stabilizer algebras of subcodes of λ-Gabidulin codes, either by fixing or varying the parent λ-Gabidulin code. It further evaluates the complexity of relevant cryptanalytic attacks, estimates key and ciphertext sizes, and implements the LGS-Niederreiter cryptosystem.

## Overview

This repository provides a performance **SageMath framework** for the study of **λ-Gabidulin codes**, their subcodes, and their security in the **rank-metric cryptography setting**.

It focuses on characterizing the algebraic structure of these codes via **stabilizer algebras**, and on evaluating their resistance against state-of-the-art cryptanalysis.

The project is designed to support the design of **post-quantum cryptographic primitives**, in particular the **LGS-Niederreiter cryptosystem**, by providing precise security estimates and structural analysis tools.

---

## Main Features

### 1. Algebraic Structural Analysis
- **Stabilizer computation**: efficient algorithms for determining the left and right stabilizer algebras of subcodes.
- **Matrix system construction**: tools to build \( M_S^d \) and \( M_S^g \) systems for analyzing code extensions.

---

### 2. Multi-Field Support (\( q \geq 2 \))

This framework overcomes standard SageMath limitations (NTL/Givaro) and supports arbitrary extensions \( \mathbb{F}_{q^m} / \mathbb{F}_q \).

- **Relative mapping engine**: formal vector-space morphisms (`to_V`) ensure compatibility across finite field backends.
- **Robust vectorization**: supports \( q = 2, 8, 16, \dots \) for both key generation and decryption workflows.

---

### 3. Cryptanalytic Suite

The framework provides a comprehensive evaluation of **rank-metric attack complexities**, including:

- **Hybrid-Kernel attacks**: complexity models based on guessing + kernel search strategies.
- **Support Minors (SM) attacks**: implementation of models from Bardet et al. (2020) and Bros et al. (2022).
- **New Hybrid attacks**: implementation of models from  H. Beeloo-Sauerbier Couv ´et al. (202).
  
---

### 4. LGS-Niederreiter Cryptosystem

- Implementation of **Niederreiter-type encryption/decryption scheme**.
- Automatic estimation of parameters for for all proposal parameters.
- Integrated evaluation of key sizes and ciphertext expansion.

---

## Project Structure

src/ # Core algebraic and coding theory algorithms
stabilizers/ # Left and right stabilizer computation
codes/ # λ-Gabidulin code construction
crypto/ # LGS-Niederreiter cryptosystem
attacks/ # Cryptanalytic models and complexity evaluation
utils/ # Finite field tools and vector space mappings
examples/ # Experimental cases and benchmarks

---

## Research Context

This implementation supports experimental and theoretical research in:

- Rank-metric code-based cryptography
- Structural invariants of λ-Gabidulin codes
- Algebraic attacks and their complexity analysis
- Post-quantum cryptographic scheme design

---

## Requirements

- SageMath ≥ 9.x
- Python 3.x (embedded in Sage)
---

### Usage Example

- left and right stabilizer algebras of subcodes
run_self_tests()
params = suggested_test_parameters(max_m=12, q_values=(2,8,16))
params

results, summary = run_experiment(
    q=2,
    delta=1,
    m=5,
    n=5,
    k=1,
    kp=4, #k_prime
    N=500,
    seed=None,
    verbose=True,
    experiment_mode="fixed_parent_random_subcodes"
)

- LGS-Niederreiter cryptosystem
p = 2
s_exp = 1
q = p**s_exp
Fq_global = GF(q, names=('z',)); (z,) = Fq_global._first_ngens(1)

n = 34
m = n 
k = 24
k_prime = 800 
delta = 1

- Automatic estimation of Attack Complexity, key and ciphertext sizes for all proposal parameters.


## Authors
Freddy Lendé Metouké, Hervé Talé Kalachi, Hermann Tchatchiem Kamche, Ousmane
Ndiaye, Sélestin Ndjeya

Research implementation for academic work on rank-metric cryptography and λ-Gabidulin codes.

---

## License

MIT License (recommended for academic and open research use)



