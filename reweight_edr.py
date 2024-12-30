#
# Written by Satoshi Ono
#
import mdtraj as md
import pandas as pd
import numpy as np
import panedr

np.random.seed(seed=123456)

edr = panedr.edr_to_df('md_all.edr',verbose=False)
traj = md.load('md_fit_all.xtc',top='md_Ref.gro')
#traj = md.load('md_noSolpbc.xtc',top='md_Ref.gro')

frames=len(edr.Potential)

P_E_T = pd.read_table('P_E_T_300.dat',
                      sep=r'\s+',header=None,names=['ene','Prob'])
P_E_T['pdf'] = np.exp(P_E_T['Prob'])
binsize = P_E_T.ene[1] - P_E_T.ene[0]

enemin = P_E_T.ene.min()

iene = ((edr.Potential - enemin)/binsize).astype(int)
iene = np.clip(iene,0, len(P_E_T) - 1)
rand = np.random.rand(frames)
mask = rand <= P_E_T.pdf.values[iene]

rew_edr=edr[mask]
rew_traj=traj[mask]

rew_edr.to_csv('reweight.csv',index=False)
rew_traj.save_xtc("300K.xtc")
