import math
import matplotlib.pyplot as plt

# =========================
# Input data
# =========================
data = [
    {"name": "This work",       "value": 10.99, "error": 1.31},
    {"name": "OMEG 1992",       "value": 7.50,  "error": 1.00},
    {"name": "OMEG 1998",       "value": 10.00, "error": 2.20},
    {"name": "CLAS 2016",       "value": 21.30, "error": 4.40},
    {"name": "VES 1995 + PDG",  "value": 18.60, "error": 6.90},
    {"name": "PDG fit",         "value": 8.55,  "error": 1.44},
]

reference = data[0]


# =========================
# Utility functions
# =========================
def sigma_difference(x1, s1, x2, s2):
    """
    Difference between two measurements in units of combined sigma:
        n_sigma = |x1 - x2| / sqrt(s1^2 + s2^2)
    """
    return abs(x1 - x2) / math.sqrt(s1**2 + s2**2)


def chi2_pair(x1, s1, x2, s2):
    """
    Chi-square for compatibility of two independent Gaussian measurements:
        chi2 = (x1 - x2)^2 / (s1^2 + s2^2)
    with 1 degree of freedom.
    """
    return (x1 - x2)**2 / (s1**2 + s2**2)


def p_value_1dof(chi2):
    """
    p-value for chi2 with 1 degree of freedom.
    For 1 dof:
        p = 1 - erf(sqrt(chi2 / 2))
    """
    return 1.0 - math.erf(math.sqrt(chi2 / 2.0))


def intervals_overlap(x1, s1, x2, s2, nsigma=1.0):
    """
    Check overlap of nsigma intervals.
    """
    a1, b1 = x1 - nsigma * s1, x1 + nsigma * s1
    a2, b2 = x2 - nsigma * s2, x2 + nsigma * s2
    return max(a1, a2) <= min(b1, b2)


def compatibility_label(nsigma):
    if nsigma < 1:
        return "excellent"
    elif nsigma < 2:
        return "good"
    elif nsigma < 3:
        return "marginal"
    else:
        return "tension"


# =========================
# Print comparison table
# =========================
print("\nComparison with reference measurement")
print(f"Reference: {reference['name']} : R = {reference['value']:.2f} ± {reference['error']:.2f}\n")

header = (
    f"{'Measurement':<18}"
    f"{'R':>8}"
    f"{'err':>8}"
    f"{'Δ/σ':>10}"
    f"{'chi2':>10}"
    f"{'p-value':>12}"
    f"{'1σ overlap':>14}"
    f"{'2σ overlap':>14}"
    f"{'status':>12}"
)
print(header)
print("-" * len(header))

for d in data[1:]:
    nsig = sigma_difference(reference["value"], reference["error"], d["value"], d["error"])
    chi2 = chi2_pair(reference["value"], reference["error"], d["value"], d["error"])
    pval = p_value_1dof(chi2)
    ov1 = intervals_overlap(reference["value"], reference["error"], d["value"], d["error"], nsigma=1)
    ov2 = intervals_overlap(reference["value"], reference["error"], d["value"], d["error"], nsigma=2)

    print(
        f"{d['name']:<18}"
        f"{d['value']:>8.2f}"
        f"{d['error']:>8.2f}"
        f"{nsig:>10.2f}"
        f"{chi2:>10.2f}"
        f"{pval:>12.4f}"
        f"{str(ov1):>14}"
        f"{str(ov2):>14}"
        f"{compatibility_label(nsig):>12}"
    )

# =========================
# Weighted average of previous measurements
# =========================
previous = data[1:]
weights = [1.0 / (d["error"] ** 2) for d in previous]
weighted_mean = sum(w * d["value"] for w, d in zip(weights, previous)) / sum(weights)
weighted_err = math.sqrt(1.0 / sum(weights))

nsig_avg = sigma_difference(reference["value"], reference["error"], weighted_mean, weighted_err)
chi2_avg = chi2_pair(reference["value"], reference["error"], weighted_mean, weighted_err)
pval_avg = p_value_1dof(chi2_avg)

print("\nWeighted average of previous measurements:")
print(f"R_prev = {weighted_mean:.2f} ± {weighted_err:.2f}")
print(f"Comparison with reference: Δ/σ = {nsig_avg:.2f}, chi2 = {chi2_avg:.2f}, p-value = {pval_avg:.4f}")

# =========================
# Plot
# =========================
labels = [d["name"] for d in data]
values = [d["value"] for d in data]
errors = [d["error"] for d in data]
y = list(range(len(data)))

plt.figure(figsize=(10, 5.5))
plt.errorbar(values, y, xerr=errors, fmt='o', capsize=4, color='navy', ecolor='black')

# Highlight reference measurement
plt.axvspan(reference["value"] - reference["error"],
            reference["value"] + reference["error"],
            color='red', alpha=0.15, label='This work ±1σ')

plt.axvline(reference["value"], color='red', linestyle='--', linewidth=1.5, label='This work central value')

plt.yticks(y, labels)
plt.xlabel("R")
plt.title("Comparison of R measurements")
plt.grid(axis='x', linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
