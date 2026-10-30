with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix the remaining apply_sommerfeld_boundary call
content = content.replace(
    'apply_sommerfeld_boundary(state_f)',
    '# apply_sommerfeld_boundary(state_f)  # Already handled by sommerfeld_all above'
)

with open('scripts/test_amr_bl_moving_patch_Lc160.py', 'w', encoding='utf-8') as f:
    f.write(content)
print('Fixed remaining Sommerfeld boundary call')