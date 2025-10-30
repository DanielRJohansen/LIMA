import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# -----------------------------
# Physical constants
# -----------------------------
e = 1.602176634e-19
epsilon0 = 8.8541878128e-12
q_O = -0.834 * e
kappa = 3.0 / 1.0e-9  # 3 nm^-1 = 3e9 m^-1

# -----------------------------
# Range of distances (0–1 nm)
# -----------------------------
r_nm = np.linspace(0.01, 1.0, 400)
r_m = r_nm * 1e-9
x = kappa * r_m  # dimensionless κr

# -----------------------------
# Exact erfc scalar
# -----------------------------
scalar_exact = erfc(x) + (2.0 / np.sqrt(np.pi)) * x * np.exp(-x**2)

# -----------------------------
# Smooth blended approximation
# -----------------------------
def g_approx(x):
    INV_SQRT_PI = 1.0 / np.sqrt(np.pi)
    out = np.zeros_like(x)

    # Blend region
    small = x <= 1.2
    large = ~small

    # Small-x Taylor
    xs = x[small]
    x2s = xs**2; x3 = xs*x2s; x5 = x3*x2s; x7 = x5*x2s
    out[small] = 1.0 + 2.0*INV_SQRT_PI*((-2/3)*x3 + (2/5)*x5 - (1/7)*x7)

    # Large-x asymptotic
    xl = x[large]
    invx = 1.0 / xl
    invx3 = invx**3
    out[large] = np.exp(-xl**2) * (2*xl + invx - 0.5*invx3) * INV_SQRT_PI

    return out

scalar_approx = g_approx(x)

# -----------------------------
# Forces
# -----------------------------
F_coul = (1/(4*np.pi*epsilon0)) * abs(q_O*q_O) / (r_m**2)
F_exact = F_coul * scalar_exact
F_approx = F_coul * scalar_approx

# -----------------------------
# Relative error
# -----------------------------
rel_err = np.abs((scalar_approx - scalar_exact) / scalar_exact)

# -----------------------------
# Plotting
# -----------------------------
fig, axes = plt.subplots(3, 1, figsize=(7, 10))

# (1) erfc scalar
axes[0].plot(x, scalar_exact, label="Exact", color="tab:blue")
axes[0].plot(x, scalar_approx, '--', label="Approximation", color="tab:orange")
axes[0].set_xlabel("x = κr")
axes[0].set_ylabel("g(x)")
axes[0].set_title("Ewald erfc-scalar vs Approximation")
axes[0].legend()
axes[0].grid(True)

# (2) resulting force
axes[1].plot(r_nm, F_exact, label="Exact screened force", color="tab:red")
axes[1].plot(r_nm, F_approx, '--', label="Approximation", color="tab:purple")
axes[1].plot(r_nm, F_coul, ':', label="Bare Coulomb", color="tab:gray")
axes[1].set_xlabel("Distance r [nm]")
axes[1].set_ylabel("Force [N]")
axes[1].set_yscale("log")
axes[1].set_title("Resulting Coulomb Force (Oxygen–Oxygen)")
axes[1].legend()
axes[1].grid(True)

# (3) error of approximation
axes[2].plot(x, rel_err, color="tab:green")
axes[2].set_xlabel("x = κr")
axes[2].set_ylabel("Relative error |Δg| / g_exact")
axes[2].set_yscale("log")
axes[2].set_title("Approximation Error vs Exact g(x)")
axes[2].grid(True, which="both")

plt.tight_layout()
plt.show()
