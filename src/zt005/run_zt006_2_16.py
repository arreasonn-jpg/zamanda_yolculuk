import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .batch_em_source import compute_T_batch
from .gauss_cyl_quadrature import gauss_cylindrical_nodes
from .cyl_quadrature import EINSTEIN_GREEN_FACTOR
from .cyl_quadrature import quadrature_integral_many
from .spatial_gr import inverse_trace_reverse
from .metric_validation import metric_invariants,timelike_check,causality_margin
OUT=Path('results/exploratory/zt006_2_16_results.json')

def obs():
    return np.array([[0.,0.,0.],[.25,0,0],[.5,0,0],[0,.25,0],[0,0,.5],[.75,0,0]])

def metric(name,ob,cyl,b,d,grid):
    pts,w=gauss_cylindrical_nodes(cyl.radius_m,cyl.height_m,*grid)
    T=compute_T_batch(name,pts,b,d,h=.005)
    hb=quadrature_integral_many(ob,pts,T,w,softening_m=.0025)
    h=np.stack([inverse_trace_reverse(x) for x in hb])
    return h,len(pts)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--mode',choices=['quick','full'],default='quick'); a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive(); ob=obs()
    grids=[(8,16,8),(12,24,12)] if a.mode=='quick' else [(12,24,12),(16,32,16),(20,40,20)]
    print('='*108); print('ZAMANDA YOLCULUK — ZT-006.2.16'); print('SOURCE-MODEL GATE + LINEARIZED METRIC VALIDATION'); print('='*108)
    print(f'Einstein Green factor: {EINSTEIN_GREEN_FACTOR:.6e}')
    results={}
    for name in ['G1','G2','G3','G4']:
        rows=[]; print(f'\n[{name}]')
        for grid in grids:
            h,N=metric(name,ob,cyl,b,d,grid); g=h[0]; inv=metric_invariants(g); tc=timelike_check(g,[1,0,0,0])
            row={'grid':list(grid),'source_points':N,'center_h00':float(g[0,0]),'center_h0i_norm':float(np.linalg.norm(g[0,1:])),'metric_invariants':inv,'rest_frame_timelike_test':tc,'causality_margin':causality_margin(g)}
            rows.append(row)
            print(f"  grid={grid} N={N} h00={row['center_h00']:.6e} |h0i|={row['center_h0i_norm']:.6e} det={inv['det']:.6e} timelike={tc['timelike']}")
        results[name]=rows
    OUT.write_text(json.dumps({'stage':'ZT-006.2.16','mode':a.mode,'observers':ob.tolist(),'results':results},indent=2),encoding='utf-8')
    print(f'\nResults written to {OUT}')
    print('SCIENTIFIC STATUS: linearized metric validation only; no CTC conclusion.')
if __name__=='__main__': main()
