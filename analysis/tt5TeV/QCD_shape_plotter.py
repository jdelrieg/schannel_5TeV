import ROOT
import matplotlib.pyplot as plt
import mplhep as hep

# Set MPLHEP style
hep.style.use("CMS")

def plot_histograms(nominal, variations, canvas_title="Histograms", save_to=None,channel='canal'):
    """
    Plots given TH1D histograms using mplhep, highlighting the nominal histogram
    and showing pairs of up/down variations with distinct alpha transparency.

    Parameters:
        nominal (ROOT.TH1D): The nominal histogram (e.g., QCD).
        variations (list of tuples): List of 4 pairs of up and down variation histograms.
                                     Each pair is a tuple: (up_hist, down_hist).
        canvas_title (str, optional): Title of the matplotlib canvas.
        save_to (str, optional): If specified, saves the plot to this file.

    Returns:
        None
    """
    # Create a figure
    plt.figure(figsize=(20, 14))
    plt.title(canvas_title)

    # Plot the nominal histogram
    if isinstance(nominal, ROOT.TH1D):
        n_bins = nominal.GetNbinsX()
        bin_edges = [nominal.GetBinLowEdge(i) for i in range(1, n_bins + 2)]
        contents = [nominal.GetBinContent(i) for i in range(1, n_bins + 1)]
        plt.step(bin_edges[:-1], contents, where="post", label=nominal.GetName(), linewidth=3, color="black")

    # Plot the variation histograms
    colors = ["blue", "green", "red", "orange"]  # Colors for each pair
    for idx, (up_hist, down_hist) in enumerate(variations):
        if idx>0: continue
        if isinstance(up_hist, ROOT.TH1D):
            n_bins = up_hist.GetNbinsX()
            bin_edges = [up_hist.GetBinLowEdge(i) for i in range(1, n_bins + 2)]
            contents = [up_hist.GetBinContent(i) for i in range(1, n_bins + 1)]
            plt.step(bin_edges[:-1], contents, where="post", label=f"{up_hist.GetName()[8:10]} (Up)", linewidth=1.5,
                     color="green", alpha=1)
        
        if isinstance(down_hist, ROOT.TH1D):
            n_bins = down_hist.GetNbinsX()
            bin_edges = [down_hist.GetBinLowEdge(i) for i in range(1, n_bins + 2)]
            contents = [down_hist.GetBinContent(i) for i in range(1, n_bins + 1)]
            plt.step(bin_edges[:-1], contents, where="post", label=f"{down_hist.GetName()[8:10]} (Down)",linewidth=1.5,
                     color="red", alpha=0.7)
                          #colors[idx]
    # Add legend
    plt.legend()
    plt.text(0.05, 0.95, channel, color="blue",
         horizontalalignment='left', verticalalignment='top',
         transform=plt.gca().transAxes)

    # Add labels
    plt.xlabel('MVA Score') if canal[-4:]=='2j1b' else plt.xlabel(r'$|\eta_{\mathrm{u},0}|$')
    plt.ylabel("Events")
    hep.cms.label(llabel='Internal', data=False, lumi=0.302,year=2017,com=5.02)

    # Show or save the plot
    if save_to:
        plt.savefig(save_to)
        print(f"Plot saved to {save_to}")
    else:
        plt.show()

# Example usage:
if __name__ == "__main__":
    canal="e_minus_2j1b"
    file = ROOT.TFile.Open(f"/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb_0QCDlepton/temp_cards/including_shapes/MVAscore_relaxed_b10_{canal}.root")  # Replace with your ROOT file path
    if canal[-4:]=='2j1b':
        file = ROOT.TFile.Open(f"/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb_btagEff/temp_cards/MVA_2j1b_absu0eta_3jnb/MVAscore_relaxed_b10_{canal}.root")   # Replace with your ROOT file path
    else:
        file = ROOT.TFile.Open(f"/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb_btagEff/temp_cards/MVA_2j1b_absu0eta_3jnb/absu0eta_{canal}.root")
        
    nominal = file.Get("QCD")
    variations = [
        (file.Get("QCD_QCD_hhUp"), file.Get("QCD_QCD_hhDown")),
        (file.Get("QCD_QCD_lhUp"), file.Get("QCD_QCD_lhDown")),
        (file.Get("QCD_QCD_hlUp"), file.Get("QCD_QCD_hlDown")),
        (file.Get("QCD_QCD_llUp"), file.Get("QCD_QCD_llDown"))
    ]

    variations = [
        (file.Get("QCD_QCD_shapeUp"), file.Get("QCD_QCD_shapeDown")),
        (file.Get("QCD_QCD_shapeUp"), file.Get("QCD_QCD_shapeDown")),
        (file.Get("QCD_QCD_shapeUp"), file.Get("QCD_QCD_shapeDown")),
        (file.Get("QCD_QCD_shapeUp"), file.Get("QCD_QCD_shapeDown"))
    ]
    plot_histograms(nominal, variations, canvas_title="QCD Variations", save_to=f"split_charge_goodJECs_mistag_comb_btagEff/temp_cards/MVA_2j1b_absu0eta_3jnb/{canal}.png",channel=canal)
    file.Close()
