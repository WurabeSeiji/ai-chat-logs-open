import numpy as np
import matplotlib.pyplot as plt

def generate_baumkuchen_2d():
    n = 10      # Width of the space
    depth = 5   # Thickness of the Baumkuchen layer
    R_out = 50.0 # Curvature radius at the outer crust
    R_in = R_out - depth
    
    # Create Flat Grid Lines (Mesh spacing 1)
    x_flat = np.linspace(-n/2, n/2, n+1)
    y_flat = np.linspace(R_in, R_out, depth+1)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # --- SUBPLOT 1: Flat Index Space ---
    ax1 = axes[0]
    ax1.set_title("Phase 1: Flat Index Space (Square Mesh)")
    
    # Draw horizontal 'thickness' lines
    for y in y_flat:
        ax1.plot([x_flat[0], x_flat[-1]], [y, y], color='blue', alpha=0.5, linewidth=1.5)
    # Draw vertical 'spatial' lines
    for x in x_flat:
        ax1.plot([x, x], [y_flat[0], y_flat[-1]], color='blue', alpha=0.5, linewidth=1.5)
        
    ax1.set_aspect('equal', adjustable='box')
    ax1.set_xlabel("Spatial Index $\\xi^1$")
    ax1.set_ylabel("Orthogonal Depth $\\xi^r$")
    ax1.grid(True, linestyle=':', alpha=0.4)
    # Match the Y scale slightly to visually align with the second plot
    ax1.set_ylim(R_in - 1, R_out + 1)
    ax1.set_xlim(-n/2 - 1, n/2 + 1)
    
    # --- SUBPLOT 2: Curved Physical Space (Baumkuchen) ---
    ax2 = axes[1]
    ax2.set_title(f"Phase 2: Curved Physical Space (Baumkuchen Slice R={int(R_out)})")
    
    # Convert perfectly to polar Baumkuchen
    # At the outer crust R_out, arc length strictly equals x_flat
    # This means the angular width of each cell is exactly d(theta) = 1 / R_out
    theta_vals = x_flat / R_out
    
    # Draw concentric arcs (horizontal lines warped by curvature)
    theta_dense = np.linspace(theta_vals[0], theta_vals[-1], 200)
    for r in y_flat:
        px = r * np.sin(theta_dense)
        pz = r * np.cos(theta_dense)
        # Highlight the "crust" layers and inner layers
        lw = 2.5 if (r == R_out or r == R_in) else 1.0
        c = 'saddlebrown' if r == R_out else 'chocolate'
        ax2.plot(px, pz, color=c, linewidth=lw, alpha=0.8)
        
    # Draw radial thickness lines (vertical lines pointing to origin)
    # These lines show that deeper space gets horizontally squeezed!
    for th in theta_vals:
        px = [R_in * np.sin(th), R_out * np.sin(th)]
        pz = [R_in * np.cos(th), R_out * np.cos(th)]
        lw = 2.5 if (th == theta_vals[0] or th == theta_vals[-1]) else 1.0
        ax2.plot(px, pz, color='chocolate', linewidth=lw, alpha=0.8)
        
    # Set TRUE equal aspect ratio to prevent fake "exaggerated distortion"
    # This guarantees the curve visually represents a true R=50 circle!
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_xlabel("Physical Space X")
    ax2.set_ylabel("Orthogonal Curvature Axis W")
    ax2.grid(True, linestyle=':', alpha=0.4)
    # Dynamically fit the limits gracefully around the gentle curve
    ax2.set_ylim(R_in - 1, R_out + 1)
    ax2.set_xlim(-n/2 - 1, n/2 + 1)
    
    plt.suptitle("Figure 0: Conceptual Cross-Section of Geometric Spatial Distortion\nA Flat Square Mesh (Spacing 1) bent into a 'Baumkuchen' layer with Orthogonal Curvature R=50", fontsize=14, y=0.95)
    plt.savefig('baumkuchen_geometry.png', dpi=150, bbox_inches='tight')
    print("Saved baumkuchen_geometry.png")

if __name__ == '__main__':
    generate_baumkuchen_2d()
