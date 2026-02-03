import ROOT
import re
path='/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs/temp_cards/MVA_splittchan_PDFtry/'

#####Extra if want to decorrelate hdamp, ISR
file = ROOT.TFile.Open(path+"MVAscore_relaxed_b10_e_plus_2j1b.root", "UPDATE")

old_names=["tchan_ISRUp","tchan_ISRDown","tbarchan_ISRUp","tbarchan_ISRDown","tchan_hdampUp","tchan_hdampDown","tbarchan_hdampUp","tbarchan_hdampDown"]

new_names = ["tchan_ISR_topUp","tchan_ISR_topDown","tbarchan_ISR_tbarUp","tbarchan_ISR_tbarDown","tchan_hdamp_topUp","tchan_hdamp_topDown","tbarchan_hdamp_tbarUp","tbarchan_hdamp_tbarDown"]

new_name=0
for old_name in old_names:
# Retrieve the histogram
    hist = file.Get(old_name)
    if hist:
    # Change the name of the histogram
        hist.SetName(new_names[new_name])

    # Write the modified histogram back to the file
        hist.Write("", ROOT.TObject.kOverwrite)

    # Optionally, remove the old key if necessary
        #file.Delete(f"{old_name};1")
    new_name=new_name+1

# Close the file
file.Close()



# Define keywords and modifications
keywords = ["ISR", "hdamp"]
modified_line = lambda k, x, y, rest: f"{k}_top      shape        {x:<12} {y:<12} {rest}"
new_line = lambda k, x, y, rest: f"{k}_tbar     shape        {x:<12} {y:<12} {rest}"

# Open the file and read all lines
file_path = path + "dat_MVAscore_relaxed_b10_e_plus_2j1b.txt"
with open(file_path, 'r') as file:
    lines = file.readlines()

# Open the file in write mode and modify content
with open(file_path, 'w') as file:
    for line in lines:
        # Check if the line contains any of the keywords
        if any(keyword in line for keyword in keywords):
            for keyword in keywords:
                if keyword in line:
                    # Use regex to match the keyword, shape, first value, second value, and the rest of the line
                    match = re.match(rf"({keyword})\s+shape\s+(\S+)\s+(\S+)(.*)", line)
                    if match:
                        k, x, y, rest_of_line = match.groups()
                        
                        # Ensure proper spacing for formatting
                        rest_of_line = rest_of_line.strip()
                        modified = modified_line(k, x, '-', rest_of_line)
                        new = new_line(k, '-', y, rest_of_line)
                        
                        # Write the modified line and the new line
                        file.write(modified + '\n')
                        file.write(new + '\n')
                        break  # Once the keyword is processed, break out of the loop
        else:
            file.write(line)
