import json
from pathlib import Path
from .model import ActiveCylinder,Backpack,Drive
from .scaling_study import source_convergence,current_scaling
OUT=Path('results/exploratory/zt006_2_results.json')
def main():
    cyl,b,drive=ActiveCylinder(),Backpack(),Drive(); res={'stage':'ZT-006.2','source_convergence':{},'current_scaling':{}}
    print('='*88); print('ZAMANDA YOLCULUK — ZT-006.2'); print('SOURCE CONVERGENCE + CURRENT SCALING'); print('='*88)
    for name in ['G1','G2','G3','G4']:
        print(f'[{name}] source convergence'); rows=source_convergence(name,cyl,b,drive); res['source_convergence'][name]=rows
        for r in rows: print(f"  N={r['source_points']:4d} h00={r['h00']:.6e} |h0i|={r['h0i_norm']:.6e} rel_change={r['relative_change_h00']}")
        print(f'[{name}] current scaling'); rows,s00,s0=current_scaling(name,cyl,b,drive); res['current_scaling'][name]={'rows':rows,'log_slope_h00_vs_I':s00,'log_slope_h0i_vs_I':s0}
        print(f'  slope h00 vs I = {s00:.5f}'); print(f'  slope |h0i| vs I = {s0:.5f}')
    OUT.write_text(json.dumps(res,indent=2)); print(f'Results written to {OUT}'); print('SCIENTIFIC STATUS: scaling/convergence study; no Earth source and no CTC conclusion.')
if __name__=='__main__': main()
