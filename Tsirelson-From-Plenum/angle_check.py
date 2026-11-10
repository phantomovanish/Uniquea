# The CHSH cap as a dial set by local non-commutativity (check filed November 10, 2026).
# For +/-1 observables (A_i^2 = B_j^2 = I) with opposite wings commuting:
#   S^2 = 4I - [A0,A1] (x) [B0,B1]          (Landau, Phys. Lett. A 120 (1987) 54)
#   max over states |<S>| = ||S|| = sqrt(4 + ||[A0,A1]|| ||[B0,B1]||) = 2 sqrt(1 + sin(a) sin(b)),
#   with sin(a) = ||[A0,A1]||/2 (Jordan's lemma: the pair splits into 2x2 blocks at angles a_k).
# Commuting (classical) readouts: a = 0, cap 2.  Anticommuting: a = pi/2, cap 2 sqrt 2 (Tsirelson, 1980).
import numpy as np
rng = np.random.default_rng(20261110)
def rand_pm1(d):
    """Random +/-1 observable: U diag(+/-1) U^dagger."""
    X = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    Q, _ = np.linalg.qr(X)
    s = rng.choice([-1, 1], size=d)
    if abs(s.sum()) == d: s[0] = -s[0]
    return Q @ np.diag(s) @ Q.conj().T
def check(d, trials=200):
    worst = 0.0
    for _ in range(trials):
        A0, A1, B0, B1 = (rand_pm1(d) for _ in range(4))
        S = np.kron(A0, B0 + B1) + np.kron(A1, B0 - B1)
        smax = np.max(np.abs(np.linalg.eigvalsh(S)))
        nA = np.linalg.norm(A0@A1 - A1@A0, 2); nB = np.linalg.norm(B0@B1 - B1@B0, 2)
        worst = max(worst, abs(smax - np.sqrt(4 + nA*nB)))
    return worst
for d in (2, 3, 4, 6):
    print(f"dimension {d} per wing: max |S_max - sqrt(4 + ||[A0,A1]|| ||[B0,B1]||)| over 200 random sets = {check(d):.1e}")
# The dial on qubits: A1 at angle a from A0 (Bloch-sphere angle), B likewise.
Z = np.diag([1, -1]).astype(complex); X = np.array([[0, 1], [1, 0]], complex)
obs = lambda t: np.cos(t)*Z + np.sin(t)*X
for a in (0, np.pi/6, np.pi/4, np.pi/3, np.pi/2):
    A0, A1, B0, B1 = obs(0), obs(a), obs(0), obs(a)
    S = np.kron(A0, B0 + B1) + np.kron(A1, B0 - B1)
    print(f"a = b = {a/np.pi:.3f} pi: S_max = {np.max(np.abs(np.linalg.eigvalsh(S))):.6f}, 2 sqrt(1 + sin^2 a) = {2*np.sqrt(1+np.sin(a)**2):.6f}")
