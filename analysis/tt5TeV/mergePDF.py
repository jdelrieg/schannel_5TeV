import os, sys
head = "\documentclass{article}\n\\usepackage{geometry}\n\\usepackage{graphicx}\n\geometry{a4paper, landscape, margin=0in}\n\n\\begin{document}\n"
end = "\n\end{document}\n"

def GetPath(var, lev, chan):
  return f"/nfs/fanae/user/juanr/www/public/tt5TeV/ljets/5jul/1l/{chan}/{lev}/{var}_{chan}_{lev}.pdf"

def AddFigures(var, lev, chan='l'):
  if len(var) != 6: 
    print("WARNING: wrong number of variables")
    return

  text = "\\begin{figure}[h]\n\centering\n"
  for i,v in enumerate(var):
    text += "\includegraphics[width=0.325\\textwidth]{%s}"%GetPath(v, lev, chan)
    text += "\\\\ \n" if i == 2 else "\n"
  text += "\end{figure}\n\n"
  return text


def ProducePDF(var, levels, channels):
  out = head
  for c in channels:
    for l in levels:
      out += AddFigures(var, l, c)
  out += end
  with open('temp.tex', 'w') as f:
    f.write(out)
  os.system('pdflatex temp.tex')
  os.system('mv temp.pdf plots.pdf')
  os.system('rm temp.aux temp.tex temp.log')
  print('>> Plots saved to plots.pdf')

var = ["counts", "medianDRjj", "minDRjj", "mjj", "mt", "sumallpt"]
levels = ['2j1b', '3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
channels = ['l']
ProducePDF(var, levels, channels)
