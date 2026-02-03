import ROOT
import re
import os

# Path to the directory containing files
path = '/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/ht_low/temp_cards/carpetasalida/'

####															FIRST PART: DECORRELATE ISR AND HDAMP IN TCHAN AND TBARCHAN										#####################################################

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
    "tchan_ISRUp", "tchan_ISRDown", "tbarchan_ISRUp", "tbarchan_ISRDown",
    "tchan_hdampUp", "tchan_hdampDown", "tbarchan_hdampUp", "tbarchan_hdampDown"
]

new_names = [
    "tchan_ISR_topUp", "tchan_ISR_topDown", "tbarchan_ISR_tbarUp", "tbarchan_ISR_tbarDown",
    "tchan_hdamp_topUp", "tchan_hdamp_topDown", "tbarchan_hdamp_tbarUp", "tbarchan_hdamp_tbarDown"
]

# Define keywords and modifications
keywords = ["ISR", "hdamp"]
modified_line = lambda k, x, y, rest: f"{k}_top      shape        {x:<12} {y:<12} {rest}"
new_line = lambda k, x, y, rest: f"{k}_tbar     shape        {x:<12} {y:<12} {rest}"


# Process all files in the directory
for file_name in os.listdir(path):
    
    full_file_path = os.path.join(path, file_name)
    
    
    if file_name.endswith('.root'):
        print(f"Processing ROOT file: {full_file_path}")
        update_root_histograms(full_file_path, old_names, new_names)
    
    elif file_name.endswith('.txt'):
        print(f"Processing text file: {full_file_path}")
        modify_text_file(full_file_path, keywords, modified_line, new_line)


####													END OF	FIRST PART: DECORRELATE ISR AND HDAMP IN TCHAN AND TBARCHAN										#####################################################

#exit()

####															SECOND PART: DECORRELATE QCD NORM AND SHAPE IN QCD PROCESSES										#####################################################

def modify_text_file_QCD(file_path, keywords, suffix):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    with open(file_path, 'w') as file:
        for line in lines:
            # Ensure we don't modify structured table rows
            if line.strip().startswith("process") or line.strip().startswith("bin") or line.strip().startswith("rate"):
                file.write(line)
                continue

            # Check if any keyword appears at the start of a valid line
            if any(re.match(rf"^{keyword}(\s+|\S)", line) for keyword in keywords):
                for keyword in keywords:
                    if keyword in line:
                        match = re.match(rf"({keyword}\S*)\s+(\S+)\s+(.*)", line)
                        if match:
                            k, label, rest_of_line = match.groups()
                            modified = f"{k}_{suffix}    {label}    {rest_of_line}"
                            file.write(modified + "\n")
                            break  # Stop processing further for this line
            else:
                file.write(line)  # Write unmodified lines


keywords_QCD = ["QCD", "QCD_shape"]  #For txt files
old_names_QCD=["QCD_QCD_shapeUp","QCD_QCD_shapeDown"]
new_names_QCD=["QCD_QCD_shape_{placeholder}Up","QCD_QCD_shape_{placeholder}Down"] #For .root files


for file_name in os.listdir(path):    
    full_file_path = os.path.join(path, file_name)
    
    suffix="e_"
    if "_m_" in full_file_path: suffix= "m_"
    if "_2j1b" in full_file_path: suffix=suffix+"2j1b"
    if "_3j1b" in full_file_path: suffix=suffix+"3j1b"
    if "_3j2b" in full_file_path: suffix=suffix+"3j2b"     
    
    
    if file_name.endswith('.root'):
        print(f"Processing ROOT file: {full_file_path}")
        formatted_names = [name.format(placeholder=suffix) for name in new_names_QCD]
        update_root_histograms(full_file_path, old_names_QCD, formatted_names)
    
    elif file_name.endswith('.txt'):
        print(f"Processing text file: {full_file_path}")
        modify_text_file_QCD(full_file_path, keywords_QCD, suffix)

