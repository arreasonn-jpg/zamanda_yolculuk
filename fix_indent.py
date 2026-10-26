with open('src/zt005/energy_budget.py', 'r') as f:
    content = f.read()

content = content.replace('        E_flat = E.reshape(-1, 3)', '    E_flat = E.reshape(-1, 3)')
content = content.replace('        B_flat = B.reshape(-1, 3)', '    B_flat = B.reshape(-1, 3)')
content = content.replace('        T_dict = em_stress_energy(E_flat, B_flat)', '    T_dict = em_stress_energy(E_flat, B_flat)')
content = content.replace('        T = T_dict["T"]  # (Nx*Ny*Nz, 4, 4)', '    T = T_dict["T"]  # (Nx*Ny*Nz, 4, 4)')
content = content.replace('        # Reshape back to grid', '    # Reshape back to grid')
content = content.replace('        T_grid = T.reshape(n, n, n, 4, 4)', '    T_grid = T.reshape(n, n, n, 4, 4)')
content = content.replace('        # Add time dimension for stored energy computation', '    # Add time dimension for stored energy computation')
content = content.replace('        T_grid = T_grid[np.newaxis, ...]  # (1, Nx, Ny, Nz, 4, 4)', '    T_grid = T_grid[np.newaxis, ...]  # (1, Nx, Ny, Nz, 4, 4)')

with open('src/zt005/energy_budget.py', 'w') as f:
    f.write(content)
print('Fixed')