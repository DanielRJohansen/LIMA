import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import Chebyshev
from scipy.special import erfc


def fit_polynomial(func, domain, degree, samples=2000):
    """
    func:    function f(x) returning y
    domain:  (a, b)
    degree:  polynomial degree
    samples: number of sample points for fitting
    """
    a, b = domain
    x = np.linspace(a, b, samples)
    y = func(x)

    cheb_fit = Chebyshev.fit(x, y, deg=degree, domain=[a, b])

    # Convert Chebyshev polynomial to power basis polynomial
    power_poly = cheb_fit.convert(kind=np.polynomial.Polynomial)

    # Evaluate approximation & error on sample grid
    y_approx = power_poly(x)
    absErr = np.abs(y_approx - y)
    relErr = np.abs((y_approx - y) / y)

    return {
        "degree": degree,
        "x": x,
        "y_true": y,
        "y_approx": y_approx,
        "absError": absErr,
        "relError": relErr,
        "coeffs": power_poly.coef,
    }


def plot_approximations(all_results, title="Function Approximation"):
    """
    all_results: list of lists
       [
         [results for domain 1],
         [results for domain 2],
         ...
       ]
    """

    n = len(all_results)
    fig, axes = plt.subplots(n, 3, figsize=(20, 5 * n))

    # Normalize shape when n = 1
    if n == 1:
        axes = np.array([axes])

    for row, results in enumerate(all_results):
        ax1, ax2, ax3 = axes[row]

        # True curve
        ax1.plot(results[0]["x"], results[0]["y_true"], label="True", linewidth=2)

        # Approximations
        for res in results:
            ax1.plot(res["x"], res["y_approx"], '--', label=f"Deg {res['degree']}")

        ax1.set_title(f"{title} — Domain {row+1}")
        ax1.set_xlabel("x")
        ax1.set_ylabel("f(x)")
        ax1.legend()

        # Error
        for res in results:
            ax2.plot(res["x"], res["relError"], label=f"Error deg {res['degree']}")

        ax2.set_title(f"Relative Error — Domain {row+1}")
        ax2.set_yscale("log")
        ax2.set_xlabel("x")
        ax2.set_ylabel("|f - p|")
        ax2.legend()

        for res in results:
            ax3.plot(res["x"], res["absError"], label=f"Error deg {res['degree']}")

        ax3.set_title(f"Absolute Error — Domain {row+1}")
        ax3.set_yscale("log")
        ax3.set_xlabel("x")
        ax3.set_ylabel("|f - p|")
        ax3.legend()

    plt.tight_layout()
    plt.show()


def ErfcScalarFunc(rSq):
    kappa = 2.5
    r = np.sqrt(rSq)
    erfc_term = erfc(kappa * r)
    exp_term = np.exp(-(kappa**2) * r * r)
    return erfc_term + (2*kappa / np.sqrt(np.pi)) * r * exp_term

def InvLenCubeFunc(distSq):
    invLen = 1. / np.sqrt(distSq)
    return invLen * invLen * invLen

def CoulumbForce(distSq):
    modifiedCoulumbConstant = 14.924181
    return InvLenCubeFunc(distSq) * ErfcScalarFunc(distSq) * modifiedCoulumbConstant

def CoulumbPotential(distSq):
    kappa = 2.5
    r = np.sqrt(distSq)
    modifiedCoulumbConstant = 14.924181
    return 1./r * erfc(r*kappa) * modifiedCoulumbConstant * 0.5

def DoApproximation(func, domain, degrees):
    results = []
    for deg in degrees:
        res = fit_polynomial(func, domain, degree=deg, samples=2000)
        results.append(res)
        
        coeffs = res["coeffs"]
        n = len(coeffs)

        # Print coefficients nicely
        print(f"\n=== Domain {domain[0]}-{domain[1]} Degree {deg} Coefficients ===")
        print("Max relative error:", np.max(res["relError"]))
        print(f"static constexpr std::array<float, {n}> coeffs_domain_{domain[0]}_{domain[1]}" + "{")
        #Print cuda code:
        for i, c in enumerate(coeffs):
            print(f"{c:.10f},")
            #print(f"constexpr float a{i} = {c:.10f};")
        print("};")
        # Build the polynomial using Horner form with fmaf
        # Example: fmaf(x, fmaf(x, fmaf(x, a3, a2), a1), a0)
        #expr = f"a{n-1}"
        #for i in range(n-2, -1, -1):
        #    expr = f"fmaf(distSq, {expr}, a{i})"

        #print("const float invLenCubedTimesErfcScalarApprox =")
        #print(f"\t{expr};\n")



        #for i, c in enumerate(res["coeffs"]):
        #    print(f"a{i} = {c:.10f}")
        

    return results





def ApproxErfcScalar():
    cutoff = 1.2
    domain = (0.0, cutoff * cutoff)
    degrees = [6, 7, 8, 10]
    
    DoApproximation(ErfcScalarFunc, domain, degrees)

def ApproxInvlenCube():
    cutoff = 1.2
    domain = (0.02, cutoff*cutoff)
    degrees = [6, 12]

    DoApproximation(InvLenCubeFunc, domain, degrees)

def ApproxCoulumbForce():
    #cutoff = 1.2    
    #domainNear = (0.01, 0.1)
    #domainMedium = (0.1, 0.4)
    #domainFar = (0.4, 1.5)

    domains = [
        (0.1, 0.51),
        (0.49, 1.5)
    ]
    degrees = [6, 8, 10, 20]
    
    res = []
    for domain in domains:
        res.append(DoApproximation(CoulumbForce, domain, degrees))

    plot_approximations(res, title="Approximation")

def ApproxCoulumbPotential():
    domains = [
        (0.1, 0.51),
        (0.49, 1.5)
    ]
    degrees = [7, 9, 10]

    res = []
    for domain in domains:
        res.append(DoApproximation(CoulumbPotential, domain, degrees))

    plot_approximations(res, title="Approximation")

def ApproxSqrt():
    
    domain =(0.01,1.5)
    degrees = [6, 8]
    res = DoApproximation(np.sqrt, domain, degrees)
    plot_approximations([res], title="Approximation")
    
# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    #ApproxErfcScalar()
    #ApproxInvlenCube()
    #ApproxSqrt()
    ApproxCoulumbForce()
    #ApproxCoulumbPotential()
