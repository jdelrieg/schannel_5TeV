import ROOT
import re
import os

# Path to the directory containing files
path = '/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb/temp_cards/QCD_withshapes_quitando_separando/'

# Function to update histogram names in a ROOT file
def update_root_histograms(file_path, old_names, new_names):
    file = ROOT.TFile.Open(file_path, "UPDATE")
    if not file:
        print(f"Error: Cannot open ROOT file {file_path}")
        return
    
    for old_name, new_name in zip(old_names, new_names):
        hist = file.Get(old_name)
        if hist:
            hist.SetName(new_name)
            hist.Write("", ROOT.TObject.kOverwrite)
        else:
            print(f"Warning: Histogram {old_name} not found in {file_path}.")
    
    file.Close()

# Function to modify lines in a text file
def modify_text_file(file_path, keywords, modified_line, new_line):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if any(keyword in line for keyword in keywords):
                for keyword in keywords:
                    if keyword in line:
                        match = re.match(rf"({keyword})\s+shape\s+(\S+)\s+(\S+)(.*)", line)
                        if match:
                            k, x, y, rest_of_line = match.groups()
                            rest_of_line = rest_of_line.strip()
                            modified = modified_line(k, x, '-', rest_of_line)
                            new = new_line(k, '-', y, rest_of_line)
                            file.write(modified + '\n')
                            file.write(new + '\n')
                            break  # Stop processing further for this line
            else:
                file.write(line)

# List of histogram name changes
old_names = [
    "QCD_QCD_llUp", "QCD_QCD_llDown","QCD_QCD_lhUp", "QCD_QCD_lhDown",
    "QCD_QCD_hlUp", "QCD_QCD_hlDown","QCD_QCD_hhUp", "QCD_QCD_hhDown"
    
]
cat="_e3j2b"
new_names = [
    f"tt_RelJER{cat}Up", f"tchan_RelJER{cat}Up", f"tW_RelJER{cat}Up", f"WJets_RelJER{cat}Up",f"DY_RelJER{cat}Up",
    f"tt_RelJER{cat}Down", f"tchan_RelJER{cat}Down", f"tW_RelJER{cat}Down", f"WJets_RelJER{cat}Down",f"DY_RelJER{cat}Down"
]

# Define keywords and modifications
keywords = ["ISR", "hdamp"]
modified_line = lambda k, x, y, rest: f"{k}_top      shape        {x:<12} {y:<12} {rest}"
new_line = lambda k, x, y, rest: f"{k}_tbar     shape        {x:<12} {y:<12} {rest}"

# Process all files in the directory
for file_name in os.listdir(path):
    full_file_path = os.path.join(path, file_name)
    
    if file_name.endswith('.root'):
        if file_name[21]=="e":cat="_e"
        elif (file_name[21]=="m") & (file_name[23]=="p"):cat="_"+file_name[21:32]
        else: cat="_"+file_name[21:33]
        new_names = [ f"QCD_QCD_ll{cat}Up", f"QCD_QCD_ll{cat}Down",f"QCD_QCD_lh{cat}Up", f"QCD_QCD_lh{cat}Down",f"QCD_QCD_hl{cat}Up", f"QCD_QCD_hl{cat}Down",f"QCD_QCD_hh{cat}Up", f"QCD_QCD_hh{cat}Down"]
        print(f"Processing ROOT file: {full_file_path}")
        update_root_histograms(full_file_path, old_names, new_names)
    
    elif file_name.endswith('.txt'):
        print(f"Processing text file: {full_file_path}")
        #modify_text_file(full_file_path, keywords, modified_line, new_line)
