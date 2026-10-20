import subprocess, pathlib, json, re
root = pathlib.Path(".")
mods = ["amr_interp","amr_bh_track","refinement","refinement_tensor",
        "bssn_amr_consistent","matter_scalar","matter_complex","matter_ptum",
        "matter_proca","matter_brans_dicke","sc_condensate","exotic_matter",
        "phase_transition_common","ptum_field","retarded_gr","dual_green_solvers",
        "green_convergence","green_solvers","solver_"]
dirs = ["src","tests","scripts",".github"]
pat = lambda m: re.compile(r"(?:from\s+\.+\s*import|from\s+\.+" + m + r"[\w\s\.]*import|from\s+zt005[\.\w]*" + m + r"\s*import|from\s+\.+" + m + r"|zt005"+chr(46)+m+r"|import\s+zt005"+chr(46)+m+r")")
out = {}
for m in mods:
    users = []
    for d in dirs:
        for p in (root / d).rglob("*.py"):
            t = p.read_text(errors="ignore")
            if pat(m).search(t):
                users.append(str(p))
    out[m] = users
print(json.dumps({k: [len(v), *[(u) for u in v]] for k, v in out.items()}, indent=1))
