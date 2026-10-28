# Solid rocket motor: tubular grain burnback and thrust neutrality

Python follow-up of my 4th-year project (HIAST, 2024, originally in MATLAB). Simulates the regression
of tubular propellant grains, the resulting burning surface, chamber pressure (steady mass balance with
Saint-Robert's law `r = a·Pc^n`) and thrust, and quantifies how *neutral* the thrust profile is.

Neutrality coefficient used: `N = (Ab,max − Ab,min) / Ab,mean` — 0 is a perfectly flat thrust curve.

## Result (Ro = 100 mm, bore 30 mm, L = 800 mm, APCP-like propellant)

| Grain | Burn time | Pc max | N (lower = flatter) | Total impulse |
|---|---|---|---|---|
| Internal-burning tube | 12.6 s | 5.1 MPa | 1.14 (strongly progressive) | 96.1 kN·s |
| Internal + end burning | 12.2 s | 3.8 MPa | 0.65 | 96.2 kN·s |
| **Tapered ends, θ = 15°** | 11.9 s | 3.8 MPa | **0.43** | 98.8 kN·s |
| Tapered ends, θ = 30° | 12.3 s | 3.8 MPa | 0.62 | 94.1 kN·s |

The end-burning surfaces regress (shrink) while the bore grows, so combining them flattens the curve;
a modest taper does it best. The map below sweeps L/D and taper angle — the same study as the 2024
project, now reproducible in one command.

![burnback](figures/burnback.png)
![neutrality map](figures/neutrality_map.png)

## Run

```bash
pip install numpy matplotlib
python burnback.py
```

## Limits and next steps

Zero-dimensional (no axial pressure drop, no erosive burning), constant `c*` and `CF`, geometric
burnback only. Next: erosive burning (Lenoir–Robillard), star and finocyl grains via level-set
burnback, and a transient chamber-filling model for ignition.
