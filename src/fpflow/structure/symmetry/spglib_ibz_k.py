#region modules
from ase import Atoms 
import numpy as np 
#endregion

#region variables
#endregion

#region functions
def ir_kpoints_with_qshift(atoms, kgrid, qshift, symprec=1e-5, eps=1e-10):
    """
    BerkeleyGW kgrid.x logic (NO time-reversal), using spglib only for rotations.

    Parameters
    ----------
    atoms : ase.Atoms
        Periodic structure (PBC must be True in all three directions).
    kgrid : (int, int, int)
        Uniform mesh (Γ-centered; dk=(0,0,0)).
    qshift : (float, float, float)
        q shift in fractional (crystal) coords (same units as the k list).
    symprec : float
        spglib symmetry tolerance (Å).
    eps : float
        Numerical tolerance for equality.

    Returns
    -------
    kpoints : (N, 3) float ndarray
        Fractional coords in [0,1) of the irreducible (k+q) list (BGW ordering).
    weights : (N,) int ndarray
        Corresponding multiplicities.
    """
    # ---------- basic checks ----------
    if not np.all(atoms.pbc):
        raise ValueError("Atoms must be periodic in all three directions.")
    nk = np.asarray(kgrid, int)
    if nk.shape != (3,) or np.any(nk <= 0):
        raise ValueError("kgrid must be three positive integers.")
    q = np.asarray(qshift, float)

    # ---------- helpers mimicking BGW behavior ----------
    def wrap01(x):
        x = np.asarray(x, float)
        return (x % 1.0 + 1.0) % 1.0

    def k_canon(x):
        # BGW calls k_range() on every candidate before compares and storage.
        # Using [0,1) canonicalization is sufficient to make direct comparisons reliable.
        return wrap01(x)

    def allclose(a, b):
        return np.all(np.abs(np.asarray(a) - np.asarray(b)) < eps, axis=-1)

    def lexi_sort_rows(K):
        # Stable lexicographic sort by (kx, ky, kz) with tolerance.
        # Build sortable keys; mergesort keeps stability.
        scale = 10**max(0, int(-np.floor(np.log10(max(eps, 1e-15)))))
        keys = np.round(K * scale).astype(np.int64)
        return np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))

    # ---------- symmetry rotations from spglib (integer matrices in crystal basis) ----------
    try:
        import spglib as spg
    except Exception as e:
        raise ImportError("This function requires the 'spglib' Python package.") from e

    lattice = np.array(atoms.cell.array, float)
    scaled_pos = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    sym = spg.get_symmetry((lattice, scaled_pos, numbers), symprec=symprec)
    R = np.asarray(sym["rotations"], int)  # (nr,3,3)
    nr = R.shape[0]

    # ---------- full Γ-centered mesh ----------
    nkf = nk[0] * nk[1] * nk[2]
    kmesh = np.empty((nkf, 3), float)
    m = 0
    for i1 in range(nk[0]):
        for i2 in range(nk[1]):
            for i3 in range(nk[2]):
                kmesh[m] = k_canon([i1 / nk[0], i2 / nk[1], i3 / nk[2]])
                m += 1

    # ---------- FOLD: reduce full mesh by all rotations (NO TRS) ----------
    kweight = np.ones(nkf, int)
    kfold   = np.zeros(nkf, int)  # points with weight==0 point to their representative

    for ik1 in range(nkf):
        k1 = kmesh[ik1]
        matched = False
        for ir in range(nr):
            kR = k_canon(R[ir] @ k1)
            for ik2 in range(ik1):
                if allclose(kR, kmesh[ik2]):
                    # fold ik1 into the chain of ik2
                    kweight[ik1] = 0
                    rep = ik2
                    while kweight[rep] == 0:
                        rep = kfold[rep]
                    kweight[rep] += 1
                    kfold[ik1] = rep
                    matched = True
                    break
            if matched:
                break

    reps_idx = [i for i in range(nkf) if kweight[i] > 0]
    K_rep = kmesh[reps_idx].copy()  # already canonical
    W_rep = kweight[reps_idx].copy()
    order = lexi_sort_rows(K_rep)
    K_rep = K_rep[order]
    W_rep = W_rep[order]

    # ---------- if q == 0: done ----------
    if np.all(np.abs(q) < eps):
        return K_rep, W_rep

    # ---------- DQUNFOLD: unique orbit of K_rep under all rotations (NO TRS) ----------
    K_unf = []
    for k0 in K_rep:
        for ir in range(nr):
            kR = k_canon(R[ir] @ k0)
            # keep unique (direct compare because all canonical)
            if not any(allclose(kR, p) for p in K_unf):
                K_unf.append(kR)
    K_unf = np.array(K_unf, float)
    W_unf = np.ones(len(K_unf), int)

    # ---------- DQSUBGRP: rotations that leave q invariant (STRICT equality, no wrapping) ----------
    # For tiny q, only ops that keep the chosen crystal axis fixed survive.
    keep = []
    for ir in range(nr):
        Rq = R[ir] @ q
        if np.all(np.abs(Rq - q) < eps):
            keep.append(ir)
    Rq = R[keep] if keep else np.zeros((0, 3, 3), int)
    # Use integer inverse; for these rotation matrices R^{-1} = R^T.
    Rq_inv = np.transpose(Rq, (0, 2, 1))

    # ---------- DQFOLD: fold K_unf by q-subgroup (NO TRS), accumulate weights ----------
    reps_q = []
    wts_q  = []
    for i in range(len(K_unf)):
        matched = False
        if len(reps_q) > 0 and len(Rq_inv) > 0:
            for Rinv in Rq_inv:
                ktest = k_canon(Rinv @ K_unf[i])
                for j, r in enumerate(reps_q):
                    if allclose(ktest, r):
                        wts_q[j] += W_unf[i]
                        matched = True
                        break
                if matched:
                    break
        if not matched:
            reps_q.append(K_unf[i])
            wts_q.append(W_unf[i])

    K_q = np.array(reps_q, float) if reps_q else np.zeros((0, 3))
    W_q = np.array(wts_q,  int) if wts_q  else np.zeros((0,), int)

    # ---------- DQSORT, then add q ONCE, wrap to [0,1) ----------
    order = lexi_sort_rows(K_q)
    K_q = K_q[order]
    W_q = W_q[order]

    K_out = wrap01(K_q + q)
    W_out = W_q

    print(K_out)
    print(W_out)

    return K_out, W_out

#endregion

#region classes
#endregion