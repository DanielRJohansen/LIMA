import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

# --- Physical constants ---
e = 1.602176634e-19          # Elementary charge [C]
epsilon0 = 8.8541878128e-12  # Vacuum permittivity [C^2/(N·m^2)]
q_O = -0.834 * e             # Oxygen partial charge in TIP3P [C]
kappa = 3 / 1.0e-9         # Ewald kappa ≈ 3.5 nm^-1  -> 3.5 / nm = 3.5e9 / m

# --- Distance range (nm) ---
r_nm = np.linspace(0.1, 1.4, 400)
r_m = r_nm * 1e-9

# --- Coulomb force magnitude (bare) ---
# F = (1/(4πϵ0)) * |q1*q2| / r^2
F_coul = (1/(4*np.pi*epsilon0)) * abs(q_O*q_O) / (r_m**2)

# --- Ewald erfc scalar ---
x = kappa * r_m
scalar = erfc(x) + (2.0/np.sqrt(np.pi)) * x * np.exp(-x**2)

# --- Screened force magnitude ---
F_ewald = F_coul * scalar


# --- Print values every 0.1 nm ---
print(f"{'r (nm)':>8} | {'erfc scalar':>12} | {'F_coul [N]':>12} | {'F_ewald [N]':>12}")
print("-" * 55)
for r_target in np.arange(0.1, 1.51, 0.1):
    idx = (np.abs(r_nm - r_target)).argmin()
    print(f"{r_nm[idx]:8.2f} | {scalar[idx]:12.6e} | {F_coul[idx]:12.6e} | {F_ewald[idx]:12.6e}")


# --- Plot ---
fig, ax1 = plt.subplots(figsize=(7,5))

ax1.set_xlabel("Distance r [nm]")
ax1.set_ylabel("Force magnitude [N]", color='tab:red')
ax1.plot(r_nm, F_coul, '--', color='tab:orange', label='Bare Coulomb')
ax1.plot(r_nm, F_ewald, '-', color='tab:red', label='Ewald-screened force')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.set_yscale('log')

ax2 = ax1.twinx()
ax2.set_ylabel("Ewald erfc-scalar", color='tab:blue')
ax2.plot(r_nm, scalar, color='tab:blue', label='erfc-scalar')
ax2.tick_params(axis='y', labelcolor='tab:blue')

fig.tight_layout()
fig.suptitle("Coulomb Force & Ewald Screening (Oxygen–Oxygen)", y=1.03, fontsize=12)
ax1.legend(loc="lower left")
ax2.legend(loc="upper right")



plt.show()
