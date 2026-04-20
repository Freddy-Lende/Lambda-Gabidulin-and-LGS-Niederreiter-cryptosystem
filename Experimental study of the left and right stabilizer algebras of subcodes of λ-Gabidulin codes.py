from sage.all import *
from collections import Counter
import random as pyrandom
import time

# ============================================================
# lambda-Gabidulin experimental script
# Fixed basis B, deterministic mu, random support g
# Documented version with q-values support: 2, 8, 16 by default
# ============================================================
#
# Important design choices
# ------------------------
# 1) The expansion basis B is fixed: we use the default F_q-basis returned
#    by Sage for K = GF(q^m).
#
# 2) The vector mu = lambda^{-1} is chosen canonically with rank_q(mu)=delta.
#
# 3) The support g is random, with rank_q(g)=n.
#
# 4) The public q-ary subcode is generated as
#        G_e = P * G_vec
#    with P of full row rank over F_q.
#
# 5) We compute only the left and right stabilizer dimensions, in order
#    to keep the script readable and pedagogical.
#
# 6) Nothing is executed automatically when the file is loaded with load(...).
#    After loading the file, the user explicitly runs:
#        run_self_tests()
#        suggested_test_parameters(...)
#        run_experiment(...)
#
# Two experimental modes
# ----------------------
#   experiment_mode = "new_parent_each_time"
#       A fresh parent object is built at each iteration.
#
#   experiment_mode = "fixed_parent_random_subcodes"
#       One parent object is built once, then only the public subcode changes.
#
# Supported q values
# ------------------
# This script is intended in particular for q in {2, 8, 16} as in the paper.
# More generally, any prime power q can be used.
# ============================================================


def seed_all(seed):
    r"""
    Seed both Sage and Python random generators.

    This is useful to make experiments reproducible.

    INPUT:
        - ``seed`` -- integer or ``None``

    OUTPUT:
        - the integer seed, or ``None``

    EXAMPLE:
        sage: seed_all(123)
        123
    """
    if seed is None:
        return None
    s = int(seed)
    set_random_seed(s)
    pyrandom.seed(s)
    return s


def check_q_prime_power(q):
    r"""
    Check that q is a prime power.

    INPUT:
        - ``q`` -- integer

    OUTPUT:
        - ``True`` or ``False``

    EXAMPLE:
        sage: check_q_prime_power(2)
        True
        sage: check_q_prime_power(8)
        True
        sage: check_q_prime_power(16)
        True
        sage: check_q_prime_power(12)
        False
    """
    q = ZZ(q)
    return q > 1 and is_prime_power(q)


def setup_fields(q, m):
    r"""
    Build the base field and the extension field.

    INPUT:
        - ``q`` -- prime power
        - ``m`` -- extension degree

    OUTPUT:
        - ``Fq``     -- the field `GF(q)`
        - ``K``      -- an explicit degree-`m` extension of `Fq`
        - ``B``      -- the polynomial basis `(1,z,...,z^{m-1})`
        - ``to_B``   -- set to ``None`` in this implementation
        - ``from_B`` -- set to ``None`` in this implementation

    IMPORTANT:
        For q = 8 or q = 16, relying on ``K.vector_space(Fq, map=True)``
        can be fragile depending on how the extension is represented.
        To avoid that issue, we use the polynomial basis directly.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(2, 4)
        sage: (Fq.cardinality(), K.cardinality(), len(B))
        (2, 16, 4)

        sage: Fq, K, B, to_B, from_B = setup_fields(8, 2)
        sage: (Fq.cardinality(), K.cardinality(), len(B))
        (8, 64, 2)

        sage: Fq, K, B, to_B, from_B = setup_fields(16, 2)
        sage: (Fq.cardinality(), K.cardinality(), len(B))
        (16, 256, 2)
    """
    q = ZZ(q)
    m = ZZ(m)

    if not check_q_prime_power(q):
        raise ValueError("q must be a prime power.")
    if m < 1:
        raise ValueError("m must be >= 1.")

    Fq = GF(q, name='a')

    if m == 1:
        K = Fq
        z = K.gen() if hasattr(K, "gen") else K.one()
        B = [K.one()]
    else:
        R = PolynomialRing(Fq, 'X')
        f = R.irreducible_element(m)
        K = Fq.extension(f, names='z')
        z = K.gen()
        B = [z**i for i in range(m)]

    # In this script we do not need explicit linear maps to/from coordinates.
    to_B = None
    from_B = None
    return Fq, K, B, to_B, from_B


def coords_in_B(x, Fq, K, to_B):
    r"""
    Coordinates of an element of `K` in the fixed polynomial basis
    `B = (1,z,...,z^{m-1})`.

    INPUT:
        - ``x``   -- element coercible to `K`
        - ``Fq``  -- base field
        - ``K``   -- extension field
        - ``to_B`` -- ignored in this implementation

    OUTPUT:
        - vector in `F_q^m`

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(8, 2)
        sage: len(coords_in_B(K(1), Fq, K, to_B))
        2
    """
    x = K(x)
    m = 1 if K.degree() is None else int(K.degree())

    if m == 1:
        return vector(Fq, [Fq(x)])

    # Depending on the concrete Sage representation of finite field elements,
    # one may have .lift() but not .polynomial().
    if hasattr(x, "lift"):
        rep = x.lift()
        coeffs = list(rep.list()) if hasattr(rep, "list") else list(rep)
    elif hasattr(x, "polynomial"):
        coeffs = list(x.polynomial().list())
    elif hasattr(x, "list"):
        coeffs = list(x.list())
    else:
        raise TypeError("Unable to extract coordinates in the polynomial basis for this field element representation.")

    coeffs += [Fq(0)] * (m - len(coeffs))
    return vector(Fq, coeffs[:m])


def phi_vec_B(c, Fq, K, to_B):
    r"""
    Vector expansion `phi_B^{vec}(c)` of a word `c in K^n`.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(8, 2)
        sage: c = [K(1), K.gen()]
        sage: len(phi_vec_B(c, Fq, K, to_B))
        4
    """
    data = []
    for x in c:
        data.extend(list(coords_in_B(x, Fq, K, to_B)))
    return vector(Fq, data)


def Fold(v, Fq, m, n):
    r"""
    Fold a vector of length `mn` into an `m x n` matrix by column blocks.

    EXAMPLE:
        sage: F2 = GF(2)
        sage: v = vector(F2, [1,0, 0,1, 1,1])
        sage: Fold(v, F2, 2, 3)
        [1 0 1]
        [0 1 1]
    """
    M = matrix(Fq, m, n)
    for j in range(n):
        block = v[j*m:(j+1)*m]
        for i in range(m):
            M[i, j] = block[i]
    return M


def Unfold(M):
    r"""
    Reverse of :func:`Fold`.

    EXAMPLE:
        sage: F2 = GF(2)
        sage: M = matrix(F2, [[1,0,1],[0,1,1]])
        sage: list(Unfold(M))
        [1, 0, 0, 1, 1, 1]
    """
    Fq = M.base_ring()
    data = []
    for j in range(M.ncols()):
        for i in range(M.nrows()):
            data.append(Fq(M[i, j]))
    return vector(Fq, data)


def rank_weight(word, Fq, K, to_B):
    r"""
    Rank weight of a word over `K` with respect to the base field `F_q`.

    INPUT:
        - ``word`` -- iterable of elements of `K`
        - ``Fq``   -- base field
        - ``K``    -- extension field
        - ``to_B`` -- ignored in this implementation

    OUTPUT:
        - integer rank over `F_q`

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(8, 2)
        sage: rank_weight([K(1), K.gen()], Fq, K, to_B) in [1,2]
        True
    """
    A = matrix(Fq, [list(coords_in_B(x, Fq, K, to_B)) for x in word])
    return A.rank()


def random_full_rank_matrix(Fq, r, c):
    r"""
    Return a random `r x c` matrix of full row rank over `F_q`.

    EXAMPLE:
        sage: F8 = GF(8)
        sage: P = random_full_rank_matrix(F8, 2, 5)
        sage: P.rank()
        2
    """
    while True:
        P = random_matrix(Fq, r, c)
        if P.rank() == r:
            return P


def element_from_coords(B, col):
    r"""
    Build an element of `K` from its coordinate column in basis `B`.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(16, 2)
        sage: x = element_from_coords(B, vector(Fq, [1,0]))
        sage: x == B[0]
        True
    """
    K = B[0].parent()
    x = K.zero()
    for i in range(len(B)):
        x += B[i] * col[i]
    return x


def random_matrix_rank_r_no_zero_columns(Fq, rows, cols, r):
    r"""
    Return a random matrix of prescribed rank and with no zero column.

    EXAMPLE:
        sage: F8 = GF(8)
        sage: M = random_matrix_rank_r_no_zero_columns(F8, 3, 3, 2)
        sage: M.rank()
        2
        sage: all(M.column(j) != 0 for j in range(M.ncols()))
        True
    """
    if not (1 <= r <= min(rows, cols)):
        raise ValueError("Invalid rank r.")
    while True:
        M = random_matrix(Fq, rows, cols)
        if M.rank() == r and all(M.column(j) != 0 for j in range(cols)):
            return M


def random_support_vector_g(Fq, K, B, n, seed=None):
    r"""
    Generate a random support vector `g = (g_1,...,g_n)` of rank `n` over `F_q`.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(8, 3)
        sage: g = random_support_vector_g(Fq, K, B, 3, seed=7)
        sage: len(g)
        3
        sage: rank_weight(g, Fq, K, to_B)
        3
    """
    seed_all(seed)
    m = len(B)
    if not (1 <= n <= m):
        raise ValueError("Need 1 <= n <= m.")
    M = random_matrix_rank_r_no_zero_columns(Fq, m, n, n)
    return vector(K, [element_from_coords(B, M.column(j)) for j in range(n)])


def canonical_mu_of_rank_delta(Fq, K, B, delta, n):
    r"""
    Build a canonical vector `mu = lambda^{-1}` with `rank_q(mu)=delta`.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(16, 2)
        sage: mu = canonical_mu_of_rank_delta(Fq, K, B, 1, 2)
        sage: len(mu)
        2
        sage: rank_weight(mu, Fq, K, to_B)
        1
    """
    m = len(B)
    if not (1 <= delta <= min(m, n)):
        raise ValueError("Need 1 <= delta <= min(m,n).")
    values = list(B[:delta])
    s = K.zero()
    for i in range(delta):
        s += B[i]
    while len(values) < n:
        values.append(s)
    return vector(K, values[:n])


def lambda_from_mu(mu):
    r"""
    Recover `lambda` from `mu = lambda^{-1}` coordinatewise.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(8, 2)
        sage: mu = canonical_mu_of_rank_delta(Fq, K, B, 1, 2)
        sage: lam = lambda_from_mu(mu)
        sage: all(lam[i] * mu[i] == 1 for i in range(len(mu)))
        True
    """
    K = mu.base_ring()
    return vector(K, [K(x)**(-1) for x in mu])


def G_lambda_generator(K, q, g, lam, k):
    r"""
    Generator matrix of `G_lambda(g,k)`.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(8, 2)
        sage: g = random_support_vector_g(Fq, K, B, 2, seed=1)
        sage: mu = canonical_mu_of_rank_delta(Fq, K, B, 1, 2)
        sage: lam = lambda_from_mu(mu)
        sage: G = G_lambda_generator(K, 8, g, lam, 1)
        sage: (G.nrows(), G.ncols())
        (1, 2)
    """
    n = len(g)
    return matrix(K, k, n, lambda i, j: lam[j] * (g[j] ** (q**i)))


def G_vec_from_G_lambda(G_lam, B, Fq, K, to_B):
    r"""
    Build the q-ary expansion `G_vec` of the parent code.

    EXAMPLE:
        sage: Fq, K, B, to_B, from_B = setup_fields(8, 2)
        sage: g = random_support_vector_g(Fq, K, B, 2, seed=1)
        sage: mu = canonical_mu_of_rank_delta(Fq, K, B, 1, 2)
        sage: lam = lambda_from_mu(mu)
        sage: G_lam = G_lambda_generator(K, 8, g, lam, 1)
        sage: G_vec = G_vec_from_G_lambda(G_lam, B, Fq, K, to_B)
        sage: G_vec.nrows()
        2
    """
    rows = []
    for row in G_lam.rows():
        row_list = list(row)
        for b in B:
            rows.append(phi_vec_B([b * x for x in row_list], Fq, K, to_B))
    return matrix(Fq, rows).row_space().basis_matrix()


def t_pub(n, k, delta):
    r"""
    Public decoding parameter from the paper:

        t_pub = floor((n-k)/(2*delta)).

    EXAMPLE:
        sage: t_pub(6, 2, 1)
        2
    """
    return (n - k) // (2 * delta)


def check_parameters(q, delta, m, n, k, kp, require_non_multiple=True):
    r"""
    Check the basic parameter constraints used in the experiments.

    EXAMPLE:
        sage: check_parameters(8, 1, 2, 2, 1, 1)
        True
        sage: check_parameters(12, 1, 2, 2, 1, 1)
        False
    """
    if not check_q_prime_power(q):
        return False
    if not (1 <= n <= m):
        return False
    if not (1 <= k <= n):
        return False
    if not (1 <= kp <= k*m):
        return False
    if t_pub(n, k, delta) < 1:
        return False
    if require_non_multiple and (kp % m == 0):
        return False
    return True


def suggested_test_parameters(max_m=12,
                              only_n_equals_m=True,
                              target_tpub_values=(2, 3),
                              delta_values=(1,),
                              q_values=(2, 8, 16)):
    r"""
    Suggest moderate test parameters for exact stabilizer computations.

    EXAMPLE:
        sage: params = suggested_test_parameters(4, q_values=(2,8,16))
        sage: all(p["q"] in [2,8,16] for p in params)
        True
        sage: all(p["t_pub"] in [2,3] for p in params)
        True
    """
    params = []
    for q in q_values:
        if not check_q_prime_power(q):
            continue
        for m in range(2, max_m + 1):
            n_values = [m] if only_n_equals_m else list(range(2, m + 1))
            for n in n_values:
                for delta in delta_values:
                    for k in range(1, n + 1):
                        tp = t_pub(n, k, delta)
                        if tp not in target_tpub_values:
                            continue
                        for s in [1, 2]:
                            kp = k * m - s
                            if check_parameters(q, delta, m, n, k, kp):
                                params.append({
                                    "q": q,
                                    "delta": delta,
                                    "m": m,
                                    "n": n,
                                    "k": k,
                                    "kp": kp,
                                    "t_pub": tp
                                })

    seen = set()
    filtered = []
    for p in params:
        key = (p["q"], p["delta"], p["m"], p["n"], p["k"], p["kp"], p["t_pub"])
        if key not in seen:
            seen.add(key)
            filtered.append(p)
    return filtered


def print_suggested_parameters(params):
    r"""
    Pretty-print a list of suggested parameter dictionaries.

    EXAMPLE:
        sage: params = suggested_test_parameters(3, q_values=(2,), target_tpub_values=(1,))
        sage: print_suggested_parameters(params)  # doctest: +ELLIPSIS
        q delta m n k kp t_pub
        ...
    """
    print("q delta m n k kp t_pub")
    for p in params:
        print(p["q"], p["delta"], p["m"], p["n"], p["k"], p["kp"], p["t_pub"])


def dim_right_stabilizer(G_e, H_e, Fq, m, n):
    r"""
    Exact dimension of the right stabilizer of the matrix image of `G_e`.

    EXAMPLE:
        sage: F8 = GF(8)
        sage: G_e = matrix(F8, [[1,0,0,1]])
        sage: H_e = LinearCode(G_e).parity_check_matrix()
        sage: dim_right_stabilizer(G_e, H_e, F8, 2, 2) >= 1
        True
    """
    I_n = identity_matrix(Fq, n)
    rows = []
    for h in H_e.rows():
        for g in G_e.rows():
            M = Fold(g, Fq, m, n)
            rows.append(h * (I_n.tensor_product(M)))
    if len(rows) == 0:
        return n*n
    S = matrix(Fq, rows)
    return S.ncols() - S.rank()


def dim_left_stabilizer(G_e, H_e, Fq, m, n):
    r"""
    Exact dimension of the left stabilizer of the matrix image of `G_e`.

    EXAMPLE:
        sage: F16 = GF(16)
        sage: G_e = matrix(F16, [[1,0,0,1]])
        sage: H_e = LinearCode(G_e).parity_check_matrix()
        sage: dim_left_stabilizer(G_e, H_e, F16, 2, 2) >= 1
        True
    """
    I_m = identity_matrix(Fq, m)
    rows = []
    for h in H_e.rows():
        for g in G_e.rows():
            M = Fold(g, Fq, m, n)
            rows.append(h * (M.transpose().tensor_product(I_m)))
    if len(rows) == 0:
        return m*m
    S = matrix(Fq, rows)
    return S.ncols() - S.rank()


def build_parent_data(q, delta, m, n, k, seed=None):
    r"""
    Build one parent lambda-Gabidulin object.

    EXAMPLE:
        sage: parent = build_parent_data(8, 1, 2, 2, 1, seed=123)
        sage: parent["G_vec"].nrows()
        2
        sage: parent["delta_observed"]
        1
    """
    seed_all(seed)

    if not check_q_prime_power(q):
        raise ValueError("q must be a prime power.")
    if not (1 <= n <= m and 1 <= k <= n):
        raise ValueError("Need 1 <= k <= n <= m.")

    Fq, K, B, to_B, from_B = setup_fields(q, m)

    g   = random_support_vector_g(Fq, K, B, n, seed=None if seed is None else int(seed) + 101)
    mu  = canonical_mu_of_rank_delta(Fq, K, B, delta, n)
    lam = lambda_from_mu(mu)

    G_lam = G_lambda_generator(K, q, g, lam, k)
    G_vec = G_vec_from_G_lambda(G_lam, B, Fq, K, to_B)

    if G_vec.nrows() != k * m:
        raise RuntimeError("Unexpected rank for G_vec.")

    return {
        "q": q,
        "delta": delta,
        "m": m,
        "n": n,
        "k": k,
        "Fq": Fq,
        "K": K,
        "B": B,
        "to_B": to_B,
        "from_B": from_B,
        "g": g,
        "mu": mu,
        "lam": lam,
        "G_lam": G_lam,
        "G_vec": G_vec,
        "delta_observed": int(rank_weight(mu, Fq, K, to_B)),
        "rank_g_observed": int(rank_weight(g, Fq, K, to_B)),
    }


def one_random_subcode_from_parent(parent_data, kp, seed=None):
    r"""
    Generate one random public q-ary subcode from a fixed parent object.

    EXAMPLE:
        sage: parent = build_parent_data(16, 1, 2, 2, 1, seed=123)
        sage: out = one_random_subcode_from_parent(parent, 1, seed=456)
        sage: "dim_right_stab" in out and "dim_left_stab" in out
        True
    """
    seed_all(seed)

    q = parent_data["q"]
    delta = parent_data["delta"]
    m = parent_data["m"]
    n = parent_data["n"]
    k = parent_data["k"]
    Fq = parent_data["Fq"]
    G_vec = parent_data["G_vec"]

    if not (1 <= kp <= k*m):
        raise ValueError("Need 1 <= k' <= k*m.")

    P   = random_full_rank_matrix(Fq, kp, k*m)
    G_e = (P * G_vec).row_space().basis_matrix()

    if G_e.nrows() != kp:
        raise RuntimeError("Unexpected rank for G_e.")

    H_e = LinearCode(G_e).parity_check_matrix()

    dim_RS = dim_right_stabilizer(G_e, H_e, Fq, m, n)
    dim_LS = dim_left_stabilizer(G_e, H_e, Fq, m, n)

    return {
        "q": q,
        "delta": delta,
        "m": m,
        "n": n,
        "k": k,
        "kp": kp,
        "t_pub": t_pub(n, k, delta),
        "delta_observed": parent_data["delta_observed"],
        "rank_g_observed": parent_data["rank_g_observed"],
        "dim_right_stab": int(dim_RS),
        "dim_left_stab": int(dim_LS),
        "is_right_stab_trivial": (dim_RS == 1),
        "is_left_stab_trivial": (dim_LS == 1),
    }


def one_sample(q, delta, m, n, k, kp, seed=None):
    r"""
    Generate one full sample in the ``new_parent_each_time`` protocol.

    EXAMPLE:
        sage: out = one_sample(8, 1, 2, 2, 1, 1, seed=123)
        sage: out["delta_observed"]
        1
    """
    parent_data = build_parent_data(q=q, delta=delta, m=m, n=n, k=k, seed=seed)
    return one_random_subcode_from_parent(parent_data, kp=kp, seed=None if seed is None else int(seed) + 303)


def summarize_results(results):
    r"""
    Summarize a list of experimental results.

    EXAMPLE:
        sage: summary = summarize_results([{"dim_right_stab":1,"dim_left_stab":1,"is_right_stab_trivial":True,"is_left_stab_trivial":True,"delta_observed":1,"rank_g_observed":2}, {"dim_right_stab":2,"dim_left_stab":1,"is_right_stab_trivial":False,"is_left_stab_trivial":True,"delta_observed":1,"rank_g_observed":2}])
        sage: summary["N"]
        2
    """
    if len(results) == 0:
        return {"N": 0}

    pairs = [(r["dim_right_stab"], r["dim_left_stab"]) for r in results]
    stats = Counter(pairs)

    summary = {
        "N": len(results),
        "dimension_pairs": stats,
        "Pr_right_trivial": RR(sum(int(r["is_right_stab_trivial"]) for r in results) / len(results)),
        "Pr_left_trivial":  RR(sum(int(r["is_left_stab_trivial"])  for r in results) / len(results)),
        "mean_dim_right_stab": RR(sum(r["dim_right_stab"] for r in results) / len(results)),
        "mean_dim_left_stab":  RR(sum(r["dim_left_stab"]  for r in results) / len(results)),
        "delta_observed_values": sorted(set(r["delta_observed"] for r in results)),
        "rank_g_observed_values": sorted(set(r["rank_g_observed"] for r in results)),
    }

    if "experiment_mode" in results[0]:
        summary["experiment_mode"] = results[0]["experiment_mode"]

    return summary


def print_summary(summary):
    r"""
    Pretty-print a summary dictionary.

    EXAMPLE:
        sage: print_summary({"N":0})
        ========================================================================
        SUMMARY
        ========================================================================
        N = 0
    """
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    if "experiment_mode" in summary:
        print("experiment_mode       =", summary["experiment_mode"])
    print("N =", summary["N"])
    if summary["N"] == 0:
        return
    print("delta_observed_values  =", summary["delta_observed_values"])
    print("rank_g_observed_values =", summary["rank_g_observed_values"])
    print("Pr_right_trivial       =", summary["Pr_right_trivial"])
    print("Pr_left_trivial        =", summary["Pr_left_trivial"])
    print("mean_dim_right_stab    =", summary["mean_dim_right_stab"])
    print("mean_dim_left_stab     =", summary["mean_dim_left_stab"])
    print("-" * 72)
    print("{:<24} | {:<10}".format("(dim_RS, dim_LS)", "frequency"))
    print("-" * 72)
    for dims, freq in sorted(summary["dimension_pairs"].items()):
        print("{:<24} | {:<10}".format(str(dims), freq))


def run_experiment(q, delta, m, n, k, kp, N=20, seed=None, verbose=True,
                   experiment_mode="new_parent_each_time"):
    r"""
    Run `N` samples for fixed parameters.

    EXAMPLE:
        sage: results, summary = run_experiment(8, 1, 2, 2, 1, 1, N=2, seed=11, verbose=False, experiment_mode="new_parent_each_time")
        sage: len(results)
        2
        sage: results2, summary2 = run_experiment(16, 1, 2, 2, 1, 1, N=2, seed=11, verbose=False, experiment_mode="fixed_parent_random_subcodes")
        sage: len(results2)
        2
    """
    if not check_parameters(q, delta, m, n, k, kp):
        raise ValueError("The tuple (q, delta, m, n, k, kp) does not satisfy the chosen constraints.")

    if experiment_mode not in ["new_parent_each_time", "fixed_parent_random_subcodes"]:
        raise ValueError("Unknown experiment_mode.")

    t0 = time.time()
    results = []

    if experiment_mode == "fixed_parent_random_subcodes":
        parent_seed = None if seed is None else int(seed) + 1
        parent_data = build_parent_data(q=q, delta=delta, m=m, n=n, k=k, seed=parent_seed)
    else:
        parent_data = None

    for s in range(N):
        sample_seed = None if seed is None else int(seed) + s + 1

        if experiment_mode == "new_parent_each_time":
            out = one_sample(q=q, delta=delta, m=m, n=n, k=k, kp=kp, seed=sample_seed)
        else:
            out = one_random_subcode_from_parent(parent_data=parent_data, kp=kp, seed=None if seed is None else int(seed) + 1000 + s)

        out["sample"] = s + 1
        out["experiment_mode"] = experiment_mode
        results.append(out)

        if verbose:
            print("[{}/{}] done".format(s + 1, N))

    elapsed = time.time() - t0
    summary = summarize_results(results)
    summary["elapsed_seconds"] = RR(elapsed)

    if experiment_mode == "fixed_parent_random_subcodes":
        summary["fixed_parent_delta_observed"] = parent_data["delta_observed"]
        summary["fixed_parent_rank_g_observed"] = parent_data["rank_g_observed"]

    if verbose:
        print_summary(summary)
        print("elapsed_seconds =", summary["elapsed_seconds"])

    return results, summary


def run_self_tests(verbose=True):
    r"""
    Run lightweight self-tests on the main functions.

    The tests intentionally include q = 2, q = 8 and q = 16.

    EXAMPLE:
        sage: run_self_tests(verbose=False)
    """
    # q = 2
    F2, K2, B2, to_B2, from_B2 = setup_fields(2, 3)
    assert len(B2) == 3
    assert len(coords_in_B(K2(1), F2, K2, to_B2)) == 3

    # q = 8
    F8, K8, B8, to_B8, from_B8 = setup_fields(8, 2)
    assert len(B8) == 2
    mu8 = canonical_mu_of_rank_delta(F8, K8, B8, 1, 2)
    assert rank_weight(mu8, F8, K8, to_B8) == 1
    g8 = random_support_vector_g(F8, K8, B8, 2, seed=7)
    assert rank_weight(g8, F8, K8, to_B8) == 2

    # q = 16
    F16, K16, B16, to_B16, from_B16 = setup_fields(16, 2)
    assert len(B16) == 2
    mu16 = canonical_mu_of_rank_delta(F16, K16, B16, 1, 2)
    assert rank_weight(mu16, F16, K16, to_B16) == 1
    g16 = random_support_vector_g(F16, K16, B16, 2, seed=9)
    assert rank_weight(g16, F16, K16, to_B16) == 2

    # Fold / Unfold
    v = vector(F2, [1,0,0,1,1,1])
    M = Fold(v, F2, 2, 3)
    assert list(Unfold(M)) == list(v)

    # Small stabilizer sanity checks
    G_e8 = matrix(F8, [[1,0,0,1]])
    H_e8 = LinearCode(G_e8).parity_check_matrix()
    assert dim_right_stabilizer(G_e8, H_e8, F8, 2, 2) >= 1
    assert dim_left_stabilizer(G_e8, H_e8, F8, 2, 2) >= 1

    # Parent / subcode / experiment
    parent8 = build_parent_data(8, 1, 2, 2, 1, seed=123)
    assert parent8["delta_observed"] == 1
    assert parent8["rank_g_observed"] == 2

    out8 = one_sample(8, 1, 3, 3, 1, 1, seed=123)
    assert "dim_right_stab" in out8 and "dim_left_stab" in out8

    results16, summary16 = run_experiment(16, 1, 3, 3, 1, 1, N=2, seed=11, verbose=False, experiment_mode="fixed_parent_random_subcodes")
    assert len(results16) == 2
    assert summary16["N"] == 2

    if verbose:
        print("All self-tests passed.")


def usage_examples():
    r"""
    Print a few ready-to-use commands.

    EXAMPLE:
        sage: usage_examples()  # doctest: +ELLIPSIS
        ...
    """
    print("Load the script, then run for example:\n")
    print('run_self_tests()')
    print()
    print('params = suggested_test_parameters(max_m=6, q_values=(2,8,16))')
    print('print_suggested_parameters(params)')
    print()
    print('results, summary = run_experiment(q=8, delta=1, m=3, n=3, k=1, kp=1, N=10, seed=20260401, verbose=True, experiment_mode="new_parent_each_time")')
    print()
    print('results, summary = run_experiment(q=16, delta=1, m=3, n=3, k=1, kp=1, N=10, seed=20260401, verbose=True, experiment_mode="fixed_parent_random_subcodes")')
    
    
run_self_tests()
params = suggested_test_parameters(max_m=12, q_values=(2,8,16))
params

results, summary = run_experiment(
    q=2,
    delta=1,
    m=5,
    n=5,
    k=1,
    kp=4,
    N=500,
    seed=None,
    verbose=True,
    experiment_mode="fixed_parent_random_subcodes"
)
