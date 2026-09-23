import math
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

# Schematic geometry for readability (not to scale)
A = np.array([-3.0, 0.0])
B = np.array([ 3.0, 0.0])
C = np.array([ 0.0, 8.0])

fig, ax = plt.subplots(figsize=(7, 7))

# Triangle
ax.plot([A[0], B[0]], [A[1], B[1]], linewidth=2)
ax.plot([A[0], C[0]], [A[1], C[1]], linewidth=2)
ax.plot([B[0], C[0]], [B[1], C[1]], linewidth=2)

# Points
ax.scatter([A[0], B[0], C[0]], [A[1], B[1], C[1]], s=50)

# Labels
ax.text(A[0]-0.35, A[1]-0.45, "A", fontsize=15)
ax.text(B[0]+0.15, B[1]-0.45, "B", fontsize=15)
ax.text(C[0]+0.15, C[1]+0.2, "C", fontsize=15)

def ang(v):
    return math.degrees(math.atan2(v[1], v[0]))

# Angle arcs
arc_C = Arc(C, width=1.8, height=1.8, angle=0,
            theta1=ang(A-C), theta2=ang(B-C), linewidth=1.5)
ax.add_patch(arc_C)
ax.text(C[0]-0.55, C[1]-0.9, r"$\theta_{ACB}$", fontsize=14)

theta1_A = ang(B-A)
theta2_A = ang(C-A)
arc_A = Arc(A, width=1.8, height=1.8, angle=0,
            theta1=min(theta1_A, theta2_A), theta2=max(theta1_A, theta2_A),
            linewidth=1.5)
ax.add_patch(arc_A)
ax.text(A[0]+0.7, A[1]+0.55, r"$\theta_{CAB}$", fontsize=14)

theta1_B = ang(C-B)
theta2_B = ang(A-B)
arc_B = Arc(B, width=1.8, height=1.8, angle=0,
            theta1=min(theta1_B, theta2_B), theta2=max(theta1_B, theta2_B),
            linewidth=1.5)
ax.add_patch(arc_B)
ax.text(B[0]-2.2, B[1]+0.55, r"$\theta_{CBA}$", fontsize=14)

# Descriptive text
ax.text(0, -0.95, "A, B: two stars", ha="center", fontsize=13)
ax.text(0, 8.95, "C: distant observer", ha="center", fontsize=13)
ax.text(0, 5.0, r"$\theta_{ACB} \ll \theta_{CAB},\ \theta_{CBA}$",
        ha="center", fontsize=15)
ax.text(0, 4.15, r"$AC \approx BC \gg AB$", ha="center", fontsize=13)
ax.text(0, 3.4, "schematic (not to scale)", ha="center", fontsize=12)

ax.set_xlim(-4.6, 4.6)
ax.set_ylim(-1.5, 10.0)
ax.set_aspect("equal", adjustable="box")
ax.axis("off")

fig.savefig("observer_geometry_theta_acb_schematic_ja.svg", format="svg", bbox_inches="tight")
fig.savefig("observer_geometry_theta_acb_schematic_ja.png", format="png", dpi=200, bbox_inches="tight")
plt.show()
