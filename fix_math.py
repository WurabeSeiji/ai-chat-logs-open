import glob

def fix_file(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()

    # The mapping function replaces X_phys[:, 0] = xi[:, 0] * g
    # With X_phys[:, 0] = (R / c) * np.arcsin(v / c)
    # Be careful to get exact replacement
    # I will look for X_phys[:, 0] = xi[:, 0] * g
    
    old_code = "X_phys[:, 0] = xi[:, 0] * g"
    new_code = "X_phys[:, 0] = (R / c) * np.arcsin(v / c)  # True time integral"
    
    if old_code in content:
        content = content.replace(old_code, new_code)
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Fixed {filename}")

fix_file('evaluate_tensor.py')
fix_file('generate_fig2.py')
