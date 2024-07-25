import numpy as np
import matplotlib.pyplot as plt



def ShowLJ():
    # Constants
    epsilon = 1.20  # Depth of the potential well
    sigma = 0.2  # Finite distance at which the inter-particle potential is zero

    # Lennard-Jones potential function
    def lennard_jones(r, epsilon, sigma):
        return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

    # Range of r values (distance)
    r = np.linspace(0.01, 0.3, 500)

    # Calculate potential
    V = lennard_jones(r, epsilon, sigma)

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(r, V, label=f'$\epsilon$={epsilon}, $\sigma$={sigma}')
    plt.xlabel('Distance (nm)')
    plt.ylabel('Potential Energy')
    plt.title('Lennard-Jones Potential')
    plt.ylim(-epsilon, 2 * epsilon)  # Setting y-limit for better visualization
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.grid(True)
    plt.show()