processes = ['TTTo2L2Nu', 'TTJets_DiLept', 'TTToSemiLeptonic', 'TT', 'WW', 'WWTo2L2Nu', 'WZ', 'WZTo2L2Q', 'WZTo3LNu', 'WZTo1L3Nu', 'ZZ', 'ZZTo2L2Nu', 'ZZTo2L2Q', 'ZZTo2Q2Nu', 'ZZTo4L', 'tbarW_noFullHad', 'tW_noFullHad', 'ST_tW_antitop_5f_inclusiveDecays', 'ST_tW_top_5f_inclusiveDecays', 'DYJetsToLL_M_10to50_MLM', 'DYJetsToLL_M_50_MLM', 'WJetsToLNu']


xsec13 = {
  'DY10to50' : 18610.0,
  'DY50' : 6025.2,
  'WJets' : 61526.7,
  'tt' : 831.76,
  'tW' : 35.85,
  'WW' : 115,
  'WWTo2L2Nu' : 12.178,
  'WZ' : 47.13,
  'WZTo2L2Q' : 5.595,
  'WZTo1L3Nu' : 4.42965,
  'WZTo3LNu' : 4.42965,
  'ZZ' : 16.523,
  'ZZTo2L2Nu' : 0.564,
  'ZZTo2L2Q' : 3.28,
  'ZZTo2Q2Nu' : 4.04,
  'ZZTo4L' : 1.256,
}


def AddXsec(text, process, value):
  nspaces = 35
  while len(process) < nspaces: process = process + ' '
  text = text + '\n' + process + ' : %1.4f'%value 
  return text

text = ''

# W+Jets
#text = AddXsec(text, 'WJetsToLNu', ( / 52850.) * xsec13['WJetsToLNu'])

# DY
text = AddXsec(text, 'DYJetsToLL_M_50_MLM',  (5729. / 5334.) * xsec13['DY50'] ) 

#text = AddXsec(text, 'DYJetsToLL_M_10to50_MLM', ( / 15890.) * xsec13['DYJetsToLL_M_10to50_MLM'])


## TT 
tt13p6 = 922.45
toDilep = 88.28769753/831.76
toSemilep = 365.3994209/831.76
text = AddXsec(text, 'TT', tt13p6)
text = AddXsec(text, 'TTTo2L2Nu', tt13p6*toDilep)
text = AddXsec(text, 'TTJets_DiLept', tt13p6*toDilep)
text = AddXsec(text, 'TTToSemiLeptonic', tt13p6*toSemilep)

# tW
toNoFullHad = 19.4674104/35.85
tW13 = xsec13['tW']
tW14 = 84.4/2
tW13p6 = (tW14 - tW13)/10*6. + tW13
text = AddXsec(text, 'tW inclusive', tW13p6)
text = AddXsec(text, 'tbarW_noFullHad', tW13p6*toNoFullHad)
text = AddXsec(text, 'tW_noFullHad', tW13p6*toNoFullHad)
text = AddXsec(text, 'ST_tW_antitop_5f_inclusiveDecays', tW13p6)
text = AddXsec(text, 'ST_tW_top_5f_inclusiveDecays', tW13p6)




# WW
WW13 = 122.470; WW13p5 = 128.673; WW14 = 134.969
WW13p6 = WW13p5 + (WW14-WW13p5)/5
text = AddXsec(text, 'WW', WW13p5)
text = AddXsec(text, 'WWTo2L2Nu', WW13p5*( xsec13['WWTo2L2Nu']/xsec13['WW']) )

# WZ
WZ13 = 50.776; WZ13p5 = 53.727; WZ14 = 56.480
WZ13p6 = WZ13p5 + (WZ14-WZ13p5)/5
text = AddXsec(text, 'WZ', WZ13p6)
text = AddXsec(text, 'WZTo2L2Q', WZ13p6 * ( xsec13['WZTo2L2Q']/xsec13['WZ'] ) )
text = AddXsec(text, 'WZTo3LNu', WZ13p6 * ( xsec13['WZTo3LNu']/xsec13['WZ'] ) )
text = AddXsec(text, 'WZTo1L3Nu', WZ13p6 * ( xsec13['WZTo1L3Nu']/xsec13['WZ'] ) )

# ZZ
ZZ13 = 15.326; ZZ13p5 = 16.276; ZZ14 = 16.843
ZZ13p6 = ZZ13p5 + (ZZ14-ZZ13p5)/5
text = AddXsec(text, 'ZZ', ZZ13p6)
text = AddXsec(text, 'ZZTo2L2Nu',  ZZ13p6 * ( xsec13['ZZTo2L2Nu']/xsec13['ZZ'] ) )
text = AddXsec(text, 'ZZTo2L2Q',   ZZ13p6 * ( xsec13['ZZTo2L2Q' ]/xsec13['ZZ'] ) )
text = AddXsec(text, 'ZZTo2Q2Nu',  ZZ13p6 * ( xsec13['ZZTo2Q2Nu']/xsec13['ZZ'] ) )
text = AddXsec(text, 'ZZTo4L',     ZZ13p6 * ( xsec13['ZZTo4L'   ]/xsec13['ZZ'] ) )


print(text)









