import ROOT as r
import sys, os
from copy import deepcopy

doSum = True # Summative norm approach

def histToGraph(hist):
    n   = hist.GetNbinsX()
    ret = r.TGraphErrors()
    for b in xrange(1,n+1):
        if hist.GetBinContent(b) == 0 and hist.GetBinError(b) == 0: continue
        ret.Set(ret.GetN() + 1)
        ret.SetPoint(ret.GetN() - 1,      hist.GetXaxis().GetBinCenter(b),      hist.GetBinContent(b))
        ret.SetPointError(ret.GetN() - 1, 0.5 * hist.GetXaxis().GetBinWidth(b), hist.GetBinError(b))
    ret.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    ret.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    return ret


def applySmoothing(path, sd, ch, lv, fixOneSides = False, defaultZeroVal = 1e-5, verbose = False):
    if verbose: print(f"> Applying preset smoothing for file {path}.")
    f = r.TFile(path, 'READ')
    outhistos = []
    hnames = [k.GetName() for k in f.GetListOfKeys()]
    for name in hnames:
        if   "_" not in name or "data_obs" in name:
            outhistos.append(deepcopy(f.Get(name).Clone(name)))
        elif "Up" not in name: continue
        else:
            processed = False
            iU = name.split("_")[-1].replace("Up", "")
            pr = name.split("_")[0]
            if sd.get(iU, None):
                if sd[iU].get(ch, None):
                    if sd[iU][ch].get(lv, None):
                        if sd[iU][ch][lv].get(pr, None):
                            ty  = sd[iU][ch][lv][pr][0]
                            sym = sd[iU][ch][lv][pr][1]
                            processed = True
                            if   ty.lower() == "norm":
    #                            print(name, pr, iU)
                                outhistos.append(deepcopy(f.Get(name).Clone(name)))                       # -2: up
                                outhistos.append(deepcopy(f.Get(name.replace("Up", "Down")).Clone(name.replace("Up", "Down")))) # -1: down
                                tmpup  = deepcopy(f.Get(name).Rebin(f.Get(name).GetNbinsX()).Clone("tmpup"))
                                tmpdn  = deepcopy(f.Get(name.replace("Up", "Down")).Rebin(f.Get(name.replace("Up", "Down")).GetNbinsX()).Clone("tmpdn"))
                                tmpnom = deepcopy(f.Get(pr).Rebin(f.Get(pr).GetNbinsX()).Clone("tmpnom"))
#                                print(tmpup.GetNbinsX(), tmpdn.GetNbinsX(), tmpnom.GetNbinsX())
                                print(tmpup.GetBinContent(1), tmpdn.GetBinContent(1), tmpnom.GetBinContent(1))
                                tmpup.Divide(tmpnom); tmpdn.Divide(tmpnom)
                                
                                isOneSided = False
                                if (tmpup.GetBinContent(1) - 1) * (tmpdn.GetBinContent(1) - 1) > 0:
                                    if verbose and fixOneSides:
                                        print("\t- WARNING: {ch}-{lv} norm. smoothing for unc. {iU} in process {pr} is one-sided. It will be fixed and symmetrised.")
                                    elif verbose:
                                        print("\t- WARNING: {ch}-{lv} norm. smoothing for unc. {iU} in process {pr} is one-sided. It won't be fixed.")
                                    isOneSided = True
                                print(tmpup.GetBinContent(1), tmpdn.GetBinContent(1))
                                if sym or (isOneSided and fixOneSides):   # Symm
                                    # Get the final ratio
                                    goodratioup = tmpup.GetBinContent(1)
                                    if goodratioup < 1: goodratioup = 1/goodratioup
                                    goodratiodn = tmpdn.GetBinContent(1)
                                    if goodratiodn < 1: goodratiodn = 1/goodratiodn
                                    
                                    finalratio = (goodratioup + goodratiodn)/2
                                    
                                    for iB in range(1, tmpnom.GetNbinsX() + 1):
                                        outhistos[-1].SetBinContent(iB, f.Get(pr).GetBinContent(iB) / finalratio if not doSum else f.Get(pr).GetBinContent(iB) * (2 - finalratio) if f.Get(pr).GetBinContent(iB) * (2 - finalratio) > defaultZeroVal else defaultZeroVal)
                                        outhistos[-1].SetBinError(iB, 0.)
                                        outhistos[-2].SetBinContent(iB, f.Get(pr).GetBinContent(iB) * finalratio)
                                        outhistos[-2].SetBinError(iB, 0.)
                                
                                else:                   # No symm
                                    for iB in range(1, tmpnom.GetNbinsX() + 1):
                                        outhistos[-1].SetBinContent(iB, tmpdn.GetBinContent(1) * f.Get(pr).GetBinContent(iB))
                                        outhistos[-1].SetBinError(iB, 0.)
                                        outhistos[-2].SetBinContent(iB, tmpup.GetBinContent(1) * f.Get(pr).GetBinContent(iB))
                                        outhistos[-2].SetBinError(iB, 0.)
                                
                            elif ty.lower() == "fit":
                                outhistos.append(deepcopy(f.Get(name).Clone()))                       # -2: up
                                outhistos.append(deepcopy(f.Get(name.replace("Up", "Down")).Clone())) # -1: down
                                tmpup  = deepcopy(f.Get(name).Clone())
                                tmpdn  = deepcopy(f.Get(name.replace("Up", "Down")).Clone())
                                tmpnom = deepcopy(f.Get(pr).Clone())
                                
                                tmpup.Divide(tmpnom); tmpdn.Divide(tmpnom)
                                
                                if len(sd[iU][ch][lv][pr]) < 3:
                                    raise RuntimeError("FATAL: {ch}-{lv} smoothing for unc. {iU} in process {pr} has no order specified for the polynom.")
                                order = sd[iU][ch][lv][pr][2]
                                
                                graphu = histToGraph(f.Get(name))
                                nu     = graphu.GetN()
                                xminu  = graphu.GetX()[0]   - graphu.GetErrorXlow(0);
                                xmaxu  = graphu.GetX()[n-1] + graphu.GetErrorXhigh(n-1)
                                polyu  = r.TF1("poly", "pol%d" % order, xminu, xmaxu)
                                resu   = graphu.Fit(polyu, "QN0S EX0")

                                graphd = histToGraph(f.Get(name).replace("Up", "Down"))
                                nd     = graphd.GetN()
                                xmind  = graphd.GetX()[0]   - graphd.GetErrorXlow(0);
                                xmaxd  = graphd.GetX()[n-1] + graphd.GetErrorXhigh(n-1)
                                polyd  = r.TF1("poly", "pol%d" % order, xmind, xmaxd)
                                resd   = graphd.Fit(polyd, "QN0S EX0")

                                for iB in range(1, outhistos[-1].GetNbinsX() + 1):
                                    outhistos[-2].SetBinContent(iB, polyu[0].Eval(outhistos[-2].GetXaxis().GetBinCenter(iB)) * tmpnom.GetBinContent(iB))
                                    outhistos[-1].SetBinContent(iB, polyd[1].Eval(outhistos[-1].GetXaxis().GetBinCenter(iB)) * tmpnom.GetBinContent(iB))
                                    
                                    if sym:
                                        isNormal = None
                                        if   outhistos[-2].GetBinContent(iB) > tmpnom.GetBinContent(iB) and outhistos[-1].GetBinContent(iB) < tmpnom.GetBinContent(iB):
                                            isNormal = True
                                        elif outhistos[-2].GetBinContent(iB) < tmpnom.GetBinContent(iB) and outhistos[-1].GetBinContent(iB) > tmpnom.GetBinContent(iB):
                                            isNormal = False
                                        else:
                                            if outhistos[-2].GetBinContent(iB) > outhistos[-1].GetBinContent(iB):
                                                isNormal = True
                                            else:
                                                isNormal = False
                                        thevar = (abs(outhistos[-2].GetBinContent(iB) - tmpnom.GetBinContent(iB)) + abs(outhistos[-1].GetBinContent(iB) - tmpnom.GetBinContent(iB)))/2.
                                        outhistos[-2].SetBinContent(iB, tmpnom.GetBinContent(iB) + thevar if isNormal else tmpnom.GetBinContent(iB) - thevar if thevar < tmpnom.GetBinContent(iB) else defaultZeroVal)
                                        outhistos[-1].SetBinContent(iB, tmpnom.GetBinContent(iB) + thevar if not isNormal else tmpnom.GetBinContent(iB) - thevar if thevar < tmpnom.GetBinContent(iB) else defaultZeroVal)
                                
                                
                            else:
                                raise RuntimeError(f"FATAL: {ch}-{lv} smoothing for unc. {iU} is of type {ty} that is not supported.")

            if not processed:
                outhistos.append(deepcopy(f.Get(name).Clone(name)))
                outhistos.append(deepcopy(f.Get(name.replace("Up", "Down")).Clone(name.replace("Up", "Down")))) # -1: down
    
    
    f.Close(); del f
    if verbose: print(f"\t- Smoothing finished!")
    os.system("mv " + path + " " + path.replace(".root", "_presmooth.root"))
    of = r.TFile(path, 'RECREATE')
    for iH in outhistos: iH.Write()
    of.Close(); del of
    return
