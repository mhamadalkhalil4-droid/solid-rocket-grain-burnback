"""Internal ballistics of a solid rocket motor with tubular grains: burning-surface evolution,
chamber pressure and thrust history, and the neutrality of the thrust profile.

Grain geometries (all cylindrical, outer radius Ro, length L, inhibited where stated):
  * internal-burning tube           — burns radially outward only
  * internal + end burning          — radial + both ends (regressive ends, progressive bore)
  * tapered-end tube (angle θ)      — bore with conical ends; the taper shifts the balance

Pressure from the steady mass balance  m_dot_gen = m_dot_nozzle:
    rho_p * Ab * r = Pc * At / c*,   with r = a * Pc^n   (Saint-Robert law)
  =>  Pc = (rho_p * a * c* * Ab / At) ** (1 / (1 - n))
Thrust  F = CF * Pc * At.

Follow-up in Python of my 4th-year project (HIAST, 2024), which did this in MATLAB.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# propellant (APCP-like) and motor
rho_p, a, n, cstar, CF = 1750.0, 3.3e-5, 0.35, 1550.0, 1.55   # r = a·Pc^n with Pc in Pa, r in m/s
Ro, L, At = 0.10, 0.80, np.pi * 0.025**2                      # grain outer radius, length, throat area [m, m, m²]
Ri0 = 0.03                                                    # initial bore radius

def burn_surface(geom, w, theta=0.0):
    """Burning surface Ab(w) for regressed web w [m]. Returns 0 once burnt out."""
    Ri = Ri0 + w
    if Ri >= Ro: return 0.0
    if geom == "internal":
        return 2 * np.pi * Ri * L
    if geom == "internal_end":
        Leff = L - 2 * w
        return 2 * np.pi * Ri * Leff + 2 * np.pi * (Ro**2 - Ri**2) if Leff > 0 else 0.0
    if geom == "tapered":                                     # conical ends of half-angle theta, burning normal to the surface
        Leff = L - 2 * w
        cone_len = (Ro - Ri) / np.tan(theta) if theta > 0 else 0.0
        cone_len = min(cone_len, Leff / 2)
        slant = cone_len / np.cos(theta)
        r_top = Ri + cone_len * np.tan(theta)
        A_cone = 2 * np.pi * (Ri + r_top) / 2 * slant
        return 2 * np.pi * Ri * max(Leff - 2 * cone_len, 0) + 2 * A_cone
    raise ValueError

def simulate(geom, theta=0.0, dt=1e-3):
    w, t, hist = 0.0, 0.0, []
    while True:
        Ab = burn_surface(geom, w, theta)
        if Ab <= 0: break
        Pc = (rho_p * a * cstar * Ab / At) ** (1 / (1 - n))
        r = a * Pc**n
        hist.append((t, w, Ab, Pc, CF * Pc * At))
        w += r * dt; t += dt
    return np.array(hist)

cases = {"internal-burning tube": simulate("internal"),
         "internal + end burning": simulate("internal_end"),
         "tapered ends θ = 15°": simulate("tapered", np.deg2rad(15)),
         "tapered ends θ = 30°": simulate("tapered", np.deg2rad(30))}

print(f"{'grain':26s} {'burn time':>9s} {'Pc max':>8s} {'Pc mean':>8s} {'neutrality N':>13s} {'total impulse':>14s}")
for name, h in cases.items():
    t, w, Ab, Pc, F = h.T
    N = (Ab.max() - Ab.min()) / Ab.mean()              # neutrality coefficient: 0 = perfectly neutral
    I = np.trapezoid(F, t)
    print(f"{name:26s} {t[-1]:8.2f} s {Pc.max()/1e6:7.2f} MPa {Pc.mean()/1e6:7.2f} MPa {N:13.3f} {I/1e3:11.1f} kN·s")

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
for name, h in cases.items():
    t, w, Ab, Pc, F = h.T
    ax[0].plot(w * 1e3, Ab, label=name); ax[1].plot(t, Pc / 1e6, label=name); ax[2].plot(t, F / 1e3, label=name)
ax[0].set_xlabel("web burnt [mm]"); ax[0].set_ylabel("burning surface [m²]"); ax[0].set_title("Burning surface vs web")
ax[1].set_xlabel("t [s]"); ax[1].set_ylabel("chamber pressure [MPa]"); ax[1].set_title("Chamber pressure"); ax[1].legend(fontsize=8)
ax[2].set_xlabel("t [s]"); ax[2].set_ylabel("thrust [kN]"); ax[2].set_title("Thrust")
for a_ in ax: a_.grid(alpha=.3)
plt.tight_layout(); plt.savefig("figures/burnback.png", dpi=130)

# sweep: neutrality vs L/D and taper angle
LDs = np.linspace(2, 8, 13); thetas = np.deg2rad(np.linspace(5, 45, 9)); Nmap = np.zeros((len(thetas), len(LDs)))
L0 = L
for i, th in enumerate(thetas):
    for j, ld in enumerate(LDs):
        L = ld * 2 * Ro; h = simulate("tapered", th, dt=2e-3); Ab = h[:, 2]; Nmap[i, j] = (Ab.max() - Ab.min()) / Ab.mean()
L = L0
plt.figure(figsize=(6.5, 4.5)); cs = plt.contourf(LDs, np.rad2deg(thetas), Nmap, levels=15, cmap="viridis"); plt.colorbar(cs, label="neutrality coefficient N")
plt.xlabel("L/D"); plt.ylabel("taper half-angle θ [°]"); plt.title("Neutrality of tapered-end grains"); plt.tight_layout(); plt.savefig("figures/neutrality_map.png", dpi=130)
