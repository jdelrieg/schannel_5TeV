import ROOT
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

def plot_root_histograms(file_path, output_dir="plots_fitseguro"):
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Open the ROOT file
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print("Error: Could not open file.")
        return
    
    # Define histogram categories
    mc_processes = ["QCD", "DY", "WJetsH", "WJetsL", "tW", "tt", "tbarchan", "tchan"]
    data_process = "data_obs"
    ignored_histograms = ["TotalBkg", "TotalProcs", "TotalSig"]
    
    # Get subfolder names
    subfolders = [key.GetName() for key in root_file.GetListOfKeys()]
    subfolders = [name for name in subfolders if name.startswith("ch")]  # Filter relevant folders
    
    for folder in subfolders:
        print(f"Processing folder: {folder}")
        dir = root_file.Get(folder)
        if not dir:
            print(f"Warning: Could not access {folder}")
            continue
        
        # Define plotting parameters
        fig = plt.figure(figsize=(8, 8))
        gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax = fig.add_subplot(gs[0])  # Upper pad
        ratio_ax = fig.add_subplot(gs[1], sharex=ax)  # Lower pad
        
        mc_histograms = []
        mc_labels = []
        mc_colors = ['blue', 'red', 'green', 'orange', 'purple', 'cyan', 'pink', 'pink']  # Define some colors
        mc_colors = ['#717581', '#94A4A2', '#3F90DA', '#92DADD', '#A96B59', '#BD1F01', '#FFA90E', '#FFA90E']
        hist_data = None
        bin_edges = None
        mc_total_uncertainty = None
        
        for i, hist_name in enumerate(mc_processes + [data_process]):
            hist = dir.Get(hist_name)
            if not hist:
                print(f"Warning: {hist_name} not found in {folder}")
                continue
            
            # Convert TH1D to numpy arrays
            bins = hist.GetNbinsX()
            bin_edges = np.array([hist.GetBinLowEdge(b+1) for b in range(bins+1)])
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute bin centers
            y = np.array([hist.GetBinContent(b+1) for b in range(bins)])
            y_err = np.array([hist.GetBinError(b+1) for b in range(bins)])
            
            if hist_name in mc_processes:
                mc_histograms.append((y, mc_colors[i % len(mc_colors)], hist_name))
                if mc_total_uncertainty is None:
                    mc_total_uncertainty = y_err ** 2
                else:
                    mc_total_uncertainty += y_err ** 2
            elif hist_name == data_process:
                hist_data = (bin_centers, y, y_err)
        
        # Plot MC histograms as stacked bars
        bottom = np.zeros_like(mc_histograms[0][0]) if mc_histograms else None
        for y, color, label in mc_histograms:
            ax.bar(bin_centers, y, width=np.diff(bin_edges), bottom=bottom, color=color, label=label, alpha=0.7)
            bottom += y
        
        # Plot MC uncertainty as horizontal patches
        if mc_total_uncertainty is not None:
            mc_total_uncertainty = np.sqrt(mc_total_uncertainty)
            for i in range(len(bin_centers)):
                ax.bar(bin_centers[i], 2 * mc_total_uncertainty[i], width=np.diff(bin_edges)[i],bottom=bottom[i] - mc_total_uncertainty[i], color='none', edgecolor='gray', hatch="\/\/", label='_nolegend_' if i > 0 else'Unc.',linewidth=0)
        
        # Plot data points centered in bins
        if hist_data:
            x, y, y_err = hist_data
            ax.errorbar(x, y, yerr=y_err, fmt='ko', label='Data', markersize=5)
        
        # Calculate and plot the ratio
        if hist_data and bottom is not None:
            ratio = y / bottom
            ratio_err = y_err / bottom
            ratio_ax.errorbar(bin_centers, ratio, yerr=ratio_err, fmt='ko', markersize=5)
            ratio_ax.axhline(1, color='black', linestyle='--', linewidth=1)  # Reference line at ratio = 1
            
    # Add uncertainty bands to the ratio plot
           
        ratio_uncertainty_upper = (bottom + mc_total_uncertainty) / bottom
        ratio_uncertainty_lower = (bottom - mc_total_uncertainty) / bottom
        ratio_uncertainty = ratio_uncertainty_upper - ratio_uncertainty_lower


    # Use fill_between to plot the uncertainty bands
        for i in range(len(bin_centers)):
            ratio_ax.bar(
                bin_centers[i],  # Bin center
                ratio_uncertainty[i],  # Height of the uncertainty band
                width=np.diff(bin_edges)[i],  # Bin width
                bottom=ratio_uncertainty_lower[i],  # Start at the lower uncertainty
                color='none', edgecolor='gray', hatch="\/\/", linewidth=0,
                label='_nolegend_' if i > 0 else 'MC Uncertainty'  # Add label only once
                )

            
        ratio_ax.set_ylabel("Data/Pred.",fontsize=16)
        ratio_ax.set_ylim(0.5, 1.5)
        
        # Finalize the plot
        #ax.set_title(folder)
        folder_title_map = {"ch1_prefit": r"$e^{-}$2j1b","ch2_prefit": r"$e^{-}$3j1b","ch3_prefit": r"$e^{-}$3j2b","ch4_prefit": r"$e^{+}$2j1b","ch5_prefit": r"$e^{+}$3j1b",
         "ch6_prefit": r"$e^{+}$3j2b","ch7_prefit": r"$\mu^{-}$2j1b","ch8_prefit": r"$\mu^{-}$3j1b","ch9_prefit": r"$\mu^{-}$3j2b","ch10_prefit": r"$\mu^{+}$2j1b","ch11_prefit": r"$\mu^{+}$3j1b",
         "ch12_prefit": r"$\mu^{+}$3j2b","ch1_postfit": r"$e^{-}$2j1b","ch2_postfit": r"$e^{-}$3j1b","ch3_postfit": r"$e^{-}$3j2b","ch4_postfit": r"$e^{+}$2j1b","ch5_postfit": r"$e^{+}$3j1b",
         "ch6_postfit": r"$e^{+}$3j2b","ch7_postfit": r"$\mu^{-}$2j1b","ch8_postfit": r"$\mu^{-}$3j1b","ch9_postfit": r"$\mu^{-}$3j2b","ch10_postfit": r"$\mu^{+}$2j1b","ch11_postfit": r"$\mu^{+}$3j1b",
         "ch12_postfit": r"$\mu^{+}$3j2b"}        
        
        
        ax.set_ylabel("Events",fontsize=16)

        # Reverse the legend order
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1][:3]+handles[::-1][4:], labels[::-1][:3]+labels[::-1][4:], loc='upper right', fontsize=14,frameon=False)

        # Remove x-ticks from the upper pad
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        ax.set_xlim(0,8)
        
        if folder[:4]=='ch1_':ax.set_ylim(0,30)
        if folder[:4]=='ch4_':ax.set_ylim(0,50)
        if folder[:4]=='ch7_':ax.set_ylim(0,60)
        if folder[:4]=='ch10':ax.set_ylim(0,70)
        
        ratio_ax.set_xlim(0,8)
        
        ratio_ax.set_xlabel("MVA score",fontsize=16) if folder[:4] in ['ch1_','ch4_','ch7_','ch10'] else ratio_ax.set_xlabel(r"$|\eta_{u_{0}}|$",fontsize=16)
        ratio_ax.set_xticks(ticks=[0,1,2,3,4,5,6,7,8],labels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],fontsize=16)  if folder[:4] in ['ch1_','ch4_','ch7_','ch10'] else ratio_ax.set_xticks(ticks=[0,2,4,6,8],labels=[0,1,2,3,4],fontsize=16)
        
        ratio_ax.tick_params(axis='y',which='both',labelsize=16) 
        ax.tick_params(axis='y',which='both',labelsize=16)   
        
        CMS  = plt.text(0., 1., r"$\bf{CMS}$ ", fontsize=23,fontfamily='TeX Gyre Heros', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes) #Preliminary taba aqui 25
        Preliminary =plt.text(0.15, 1., r"$\it{Preliminary}$", fontsize=18,fontfamily='TeX Gyre Heros', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)         #19
        lumi = plt.text(1., 1., r"%1.0f $\mathrm{pb}^{-1}$ (%s)" % (302, "5.02 TeV"),fontsize=18, horizontalalignment='right', verticalalignment='bottom',transform=ax.transAxes)
        lab = plt.text(0.03, .98, folder_title_map.get(folder,folder), fontsize=16, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)      
        
        # Save the plot
        output_path = os.path.join(output_dir, f"{folder}.png")
        plt.savefig(output_path)
        plt.close()
        print(f"Saved plot: {output_path}")
    
    root_file.Close()

# Example usage
plot_root_histograms("/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/ht_low/temp_cards/carpetasalida/e_minuse_plusm_minusm_plus/Confit_seguro.root")
