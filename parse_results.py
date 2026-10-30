import json
data = json.load(open('results/validated/amr_bl_moving/amr_bl_moving_proj.json'))
print("Keys:", data.keys())
print("Length:", len(data['steps']))
for i in range(len(data['steps'])):
    if i < 30 or i > len(data['steps']) - 10:
        print(f'Step {i}: H_c={data["H_c_max"][i]}, H_f={data["H_f_max"][i]}, COM={data["com_x"][i]}, sep={data["separation"][i]}')