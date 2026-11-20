#!/usr/bin/env python3
"""Paired synthetic frequency sensitivity; other parameters stay fixed."""
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
import model


def main():
    results = []
    for multiplier in (.9, 1.1):
        p=replace(model.P, omega=model.P.omega*multiplier)
        for option in ('A','B'):
            print(f'Frequency multiplier {multiplier}, option {option}',flush=True)
            t,z=model.run(option,p)
            last=z[-p.samples_per_cycle-1:]
            H,R=model.ledger(z,option,p)
            results.append(dict(frequency_multiplier=multiplier,omega=p.omega,option=option,
                final_cycle_theta_min=float(np.min(last[:,0])),final_cycle_theta_max=float(np.max(last[:,0])),
                final_cycle_rms_theta=float(np.sqrt(np.mean(last[:,0]**2))),
                final_cycle_fraction_in_transition=float(np.mean((abs(last[:,0])>p.transition_lo)&(abs(last[:,0])<p.transition_hi))),
                final_cycle_closure_norm=float(np.linalg.norm(last[-1,:16]-last[0,:16])),
                cumulative_net_controller_work=float(np.sum(z[-1,model.POS]-z[-1,model.NEG])),
                cumulative_passive_bath=float(z[-1,model.LOSS]),
                max_abs_ledger_residual=float(np.max(abs(R)))))
    outfile=Path(__file__).resolve().parent/'output'/'sensitivity.json'
    outfile.parent.mkdir(exist_ok=True)
    outfile.write_text(json.dumps(results,indent=2)+'\n')
    print('Wrote',outfile,flush=True)

if __name__ == '__main__':
    main()
