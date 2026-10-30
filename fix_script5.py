with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix the project_phi_band_covariant call - it returns a tuple
old_projection = """        # Projection: covariant phi + momentum correction in band
        # Compute Hamiltonian constraint residual for fine grid
        H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h=dx_f)
        state_f_phi = project_phi_band_covariant(state_f.phi, state_f.gtilde, dx_f, H_f, 
                                                  n_layers=n_layers, mask_threshold=mask_threshold)
        state_f = BSSNState(
            phi=state_f_phi, gtilde=state_f.gtilde,
            K=state_f.K, Atilde=state_f.Atilde,
            Gamma=state_f.Gamma, alpha=state_f.alpha,
            beta=state_f.beta
        )
        # momentum_band_correct needs At, gt, h
        At_new, _, _ = momentum_band_correct(state_f.Atilde, state_f.gtilde, dx_f, 
                                              n_layers=n_layers, passes=1)
        state_f = BSSNState(
            phi=state_f.phi, gtilde=state_f.gtilde,
            K=state_f.K, Atilde=At_new,
            Gamma=state_f.Gamma, alpha=state_f.alpha,
            beta=state_f.beta
        )"""

new_projection = """        # Projection: covariant phi + momentum correction in band
        # Compute Hamiltonian constraint residual for fine grid
        H_f = hamiltonian_constraint_grid(state_f.phi, state_f.gtilde, state_f.K, state_f.Atilde, h=dx_f)
        state_f_phi, _ = project_phi_band_covariant(state_f.phi, state_f.gtilde, dx_f, H_f, 
                                                  n_layers=n_layers, mask_threshold=mask_threshold)
        state_f = BSSNState(
            phi=state_f_phi, gtilde=state_f.gtilde,
            K=state_f.K, Atilde=state_f.Atilde,
            Gamma=state_f.Gamma, alpha=state_f.alpha,
            beta=state_f.beta
        )
        # momentum_band_correct needs At, gt, h
        At_new, _, _ = momentum_band_correct(state_f.Atilde, state_f.gtilde, dx_f, 
                                              n_layers=n_layers, passes=1)
        state_f = BSSNState(
            phi=state_f.phi, gtilde=state_f.gtilde,
            K=state_f.K, Atilde=At_new,
            Gamma=state_f.Gamma, alpha=state_f.alpha,
            beta=state_f.beta
        )"""

content = content.replace(old_projection, new_projection)

with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'w', encoding='utf-8') as f:
    f.write(content)
print('Fixed projection return value handling')