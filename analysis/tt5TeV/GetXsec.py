import numpy as np
import sys, os

nomXsc = 66.8
lumi_rel = 0.015
args = sys.argv[1:]
if len(args) != 3:
  print('Usage: python GetXsec.py mu syst stat')
  exit()

# Get values
val, syst, stat = args
val = float(val)
syst = float(syst)
stat = float(stat)
syst_rel = syst/val
stat_rel = stat/val

# Remove lumi from syst
syst_rel = np.sqrt(syst_rel * syst_rel - lumi_rel*lumi_rel)

xsec = nomXsc*val
xsec_syst = xsec*syst_rel
xsec_stat = xsec*stat_rel
xsec_lumi = xsec*lumi_rel

tot_rel = np.sqrt(syst_rel*syst_rel + stat_rel*stat_rel + lumi_rel*lumi_rel)
tot = xsec*tot_rel

print(f'xsec =   {xsec:.2f} +/- {tot:.2f} ({tot_rel*100:.2f} %) (tot)')

print(f'     +/- {xsec_syst:.2f} ({syst_rel*100:.2f} %) (syst)')
print(f'     +/- {xsec_stat:.2f} ({stat_rel*100:.2f} %) (stat)')
print(f'     +/- {xsec_lumi:.2f} ({lumi_rel*100:.2f} %) (lumi)')

print(f'     +/- {tot:.2f} ({tot_rel*100:.2f} %) (total)')
