import os, sys, argparse

POIs   = ["r_tch"]#,"r_tchbar"] 

individualList = ['Lumi','btagSF','trigSF','mc_stat','FSR','prefire','JES','tW','WJets','elecSF','muonSF','QCD','hdamp',
                  'MET_UnclusteredEnergy','Scales','UE','DY','ISR',"JER", "PDF", "VV", "tchan", "Nonprompt",]

globalList = ['Lumi', 'systematics']#

combineargs = "--rMin 0.5 --rMax 2.0 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=5000000 --robustFit 1"# --setRobustFitStrategy 2"
                                                                                                                                                                                                                                                                              #    
#basecommand = 'combineTool.py -M MultiDimFit --floatOtherPOIs=1 -m 125 --split-points 1 --expectSignal=1 --saveInactivePOI 1 {algosettings} {combineargs}   {parallel} {queue} {asimov} {extra}'
#El de arriba es el original
basecommand = 'combineTool.py -M MultiDimFit  -m 125 --split-points 1 --setParameters r_tch=1  {algosettings} {combineargs}   {parallel} {queue} {asimov} {extra}'
#El de arriba el que se hace al fit conjunto                                           #o bar segun tenga que ser
#basecommand = 'combineTool.py -M MultiDimFit  -m 125 --split-points 1 --setParameters r_tchbar=1,r_tch=1 --redefineSignalPOI r_tchbar,r_tch  --floatOtherPOIs=1 {algosettings} {combineargs}   {parallel} {queue} {asimov} {extra}'
#El de arriba el que se  hace a cada uno de los fits separados


removelist  = []#"prop_binch13_bin0", "JES", "Nonprompt", "VV", "PDF1", "PDF2", "PDF3", "PDF4", "PDF5", "PDF6", "PDF7", "PDF8", "PDF9", "PDF10",
               #"PDF11", "PDF12", "PDF13", "PDF14", "PDF15", "PDF16", "PDF17", "PDF18", "PDF19", "PDF20", "PDF21", "PDF22", "PDF23", "PDF24", 
               #"PDF25", "PDF26", "PDF27", "PDF28", "PDF29", "PDF30", "PDF31", "PDF32",]

systsGroup = {
 'systematics': [ 
	
    #e_minus m_minus (antitop)
	#"AbsMPF", "AbsScale", "AbsStat", "DY", "ECAL", "FSR", "Flavor", "Frag", "HCAL", "ISR", "JER", "L3Res", "MC", "MET_UnclusteredEnergy", "PDF1", "PDF10", "PDF100", "PDF101", "PDF102", "PDF11", "PDF12", "PDF13", "PDF14", "PDF15", "PDF16", "PDF17", "PDF18", "PDF19", "PDF2", "PDF20", "PDF21", "PDF22", "PDF23", "PDF24", "PDF25", "PDF26", "PDF27", "PDF28", "PDF29", "PDF3", "PDF30", "PDF31", "PDF32", "PDF33", "PDF34", "PDF35", "PDF36", "PDF37", "PDF38", "PDF39", "PDF4", "PDF40", "PDF41", "PDF42", "PDF43", "PDF44", "PDF45", "PDF46", "PDF47", "PDF48", "PDF49", "PDF5", "PDF50", "PDF51", "PDF52", "PDF53", "PDF54", "PDF55", "PDF56", "PDF57", "PDF58", "PDF59", "PDF6", "PDF60", "PDF61", "PDF62", "PDF63", "PDF64", "PDF65", "PDF66", "PDF67", "PDF68", "PDF69", "PDF7", "PDF70", "PDF71", "PDF72", "PDF73", "PDF74", "PDF75", "PDF76", "PDF77", "PDF78", "PDF79", "PDF8", "PDF80", "PDF81", "PDF82", "PDF83", "PDF84", "PDF85", "PDF86", "PDF87", "PDF88", "PDF89", "PDF9", "PDF90", "PDF91", "PDF92", "PDF93", "PDF94", "PDF95", "PDF96", "PDF97", "PDF98", "PDF99", "QCD", "RelBal", "RelJER", "RelPt", "RelStat", "Scales", "UE", "WJetsH", "WJetsL", "btagSF", "elecSF", "hdamp", "muonSF", "prefire", "prop_binch1_bin0", "prop_binch1_bin1", "prop_binch1_bin2", "prop_binch1_bin3", "prop_binch1_bin4", "prop_binch1_bin5", "prop_binch1_bin6", "prop_binch1_bin7", "prop_binch2_bin0", "prop_binch2_bin1", "prop_binch2_bin2", "prop_binch2_bin3", "prop_binch2_bin4", "prop_binch2_bin5", "prop_binch2_bin6", "prop_binch2_bin7", "prop_binch3_bin0", "prop_binch3_bin1", "prop_binch3_bin2", "prop_binch3_bin3_QCD", "prop_binch3_bin3_WJetsH", "prop_binch3_bin3_tW", "prop_binch3_bin3_tchan", "prop_binch3_bin3_tt", "prop_binch3_bin4_QCD", "prop_binch3_bin4_WJetsH", "prop_binch3_bin4_tW", "prop_binch3_bin4_tchan", "prop_binch3_bin4_tt", "prop_binch3_bin5_QCD", "prop_binch3_bin5_WJetsH", "prop_binch3_bin5_tW", "prop_binch3_bin5_tchan", "prop_binch3_bin5_tt", "prop_binch3_bin6", "prop_binch3_bin7", "prop_binch4_bin0", "prop_binch4_bin1", "prop_binch4_bin2", "prop_binch4_bin3", "prop_binch4_bin4", "prop_binch4_bin5", "prop_binch4_bin6", "prop_binch4_bin7", "prop_binch5_bin0", "prop_binch5_bin1", "prop_binch5_bin2", "prop_binch5_bin3", "prop_binch5_bin4", "prop_binch5_bin5", "prop_binch5_bin6", "prop_binch5_bin7", "prop_binch6_bin0", "prop_binch6_bin1", "prop_binch6_bin2", "prop_binch6_bin3", "prop_binch6_bin4", "prop_binch6_bin5", "prop_binch6_bin6", "prop_binch6_bin7", "tW", "trigSF", "tt"
	"AbsMPF", "AbsScale", "AbsStat", "DY", "ECAL", "FSR", "Flavor", "Frag", "HCAL", "ISR_tbar", "ISR_top", "JER", "L3Res", "MC", "MET_UnclusteredEnergy", "PDF_tbar1", "PDF_tbar10", "PDF_tbar100", "PDF_tbar101", "PDF_tbar102", "PDF_tbar11", "PDF_tbar12", "PDF_tbar13", "PDF_tbar14", "PDF_tbar15", "PDF_tbar16", "PDF_tbar17", "PDF_tbar18", "PDF_tbar19", "PDF_tbar2", "PDF_tbar20", "PDF_tbar21", "PDF_tbar22", "PDF_tbar23", "PDF_tbar24", "PDF_tbar25", "PDF_tbar26", "PDF_tbar27", "PDF_tbar28", "PDF_tbar29", "PDF_tbar3", "PDF_tbar30", "PDF_tbar31", "PDF_tbar32", "PDF_tbar33", "PDF_tbar34", "PDF_tbar35", "PDF_tbar36", "PDF_tbar37", "PDF_tbar38", "PDF_tbar39", "PDF_tbar4", "PDF_tbar40", "PDF_tbar41", "PDF_tbar42", "PDF_tbar43", "PDF_tbar44", "PDF_tbar45", "PDF_tbar46", "PDF_tbar47", "PDF_tbar48", "PDF_tbar49", "PDF_tbar5", "PDF_tbar50", "PDF_tbar51", "PDF_tbar52", "PDF_tbar53", "PDF_tbar54", "PDF_tbar55", "PDF_tbar56", "PDF_tbar57", "PDF_tbar58", "PDF_tbar59", "PDF_tbar6", "PDF_tbar60", "PDF_tbar61", "PDF_tbar62", "PDF_tbar63", "PDF_tbar64", "PDF_tbar65", "PDF_tbar66", "PDF_tbar67", "PDF_tbar68", "PDF_tbar69", "PDF_tbar7", "PDF_tbar70", "PDF_tbar71", "PDF_tbar72", "PDF_tbar73", "PDF_tbar74", "PDF_tbar75", "PDF_tbar76", "PDF_tbar77", "PDF_tbar78", "PDF_tbar79", "PDF_tbar8", "PDF_tbar80", "PDF_tbar81", "PDF_tbar82", "PDF_tbar83", "PDF_tbar84", "PDF_tbar85", "PDF_tbar86", "PDF_tbar87", "PDF_tbar88", "PDF_tbar89", "PDF_tbar9", "PDF_tbar90", "PDF_tbar91", "PDF_tbar92", "PDF_tbar93", "PDF_tbar94", "PDF_tbar95", "PDF_tbar96", "PDF_tbar97", "PDF_tbar98", "PDF_tbar99", "PDF_top1", "PDF_top10", "PDF_top100", "PDF_top101", "PDF_top102", "PDF_top11", "PDF_top12", "PDF_top13", "PDF_top14", "PDF_top15", "PDF_top16", "PDF_top17", "PDF_top18", "PDF_top19", "PDF_top2", "PDF_top20", "PDF_top21", "PDF_top22", "PDF_top23", "PDF_top24", "PDF_top25", "PDF_top26", "PDF_top27", "PDF_top28", "PDF_top29", "PDF_top3", "PDF_top30", "PDF_top31", "PDF_top32", "PDF_top33", "PDF_top34", "PDF_top35", "PDF_top36", "PDF_top37", "PDF_top38", "PDF_top39", "PDF_top4", "PDF_top40", "PDF_top41", "PDF_top42", "PDF_top43", "PDF_top44", "PDF_top45", "PDF_top46", "PDF_top47", "PDF_top48", "PDF_top49", "PDF_top5", "PDF_top50", "PDF_top51", "PDF_top52", "PDF_top53", "PDF_top54", "PDF_top55", "PDF_top56", "PDF_top57", "PDF_top58", "PDF_top59", "PDF_top6", "PDF_top60", "PDF_top61", "PDF_top62", "PDF_top63", "PDF_top64", "PDF_top65", "PDF_top66", "PDF_top67", "PDF_top68", "PDF_top69", "PDF_top7", "PDF_top70", "PDF_top71", "PDF_top72", "PDF_top73", "PDF_top74", "PDF_top75", "PDF_top76", "PDF_top77", "PDF_top78", "PDF_top79", "PDF_top8", "PDF_top80", "PDF_top81", "PDF_top82", "PDF_top83", "PDF_top84", "PDF_top85", "PDF_top86", "PDF_top87", "PDF_top88", "PDF_top89", "PDF_top9", "PDF_top90", "PDF_top91", "PDF_top92", "PDF_top93", "PDF_top94", "PDF_top95", "PDF_top96", "PDF_top97", "PDF_top98", "PDF_top99", "QCD_e_2j1b","QCD_e_3j1b","QCD_e_3j2b","QCD_m_2j1b","QCD_m_3j1b","QCD_m_3j2b", "RelBal", "RelJER", "RelPt", "RelStat", "Scales_tbar", "Scales_top", "UE", "WJetsH", "WJetsL", "btagSFbc","btagSFlight", "elecSF", "hdamp_tbar", "hdamp_top", "muonSF", "prefire", "prop_binch10_bin0", "prop_binch10_bin1", "prop_binch10_bin2", "prop_binch10_bin3", "prop_binch10_bin4", "prop_binch10_bin5", "prop_binch10_bin6", "prop_binch10_bin7", "prop_binch11_bin0", "prop_binch11_bin1", "prop_binch11_bin2", "prop_binch11_bin3", "prop_binch11_bin4", "prop_binch11_bin5", "prop_binch11_bin6", "prop_binch11_bin7", "prop_binch12_bin0", "prop_binch12_bin1", "prop_binch12_bin2", "prop_binch12_bin3", "prop_binch12_bin4", "prop_binch12_bin5", "prop_binch12_bin6", "prop_binch12_bin7", "prop_binch1_bin0", "prop_binch1_bin1", "prop_binch1_bin2", "prop_binch1_bin3", "prop_binch1_bin4", "prop_binch1_bin5", "prop_binch1_bin6", "prop_binch1_bin7", "prop_binch2_bin0", "prop_binch2_bin1", "prop_binch2_bin2", "prop_binch2_bin3", "prop_binch2_bin4", "prop_binch2_bin5", "prop_binch2_bin6", "prop_binch2_bin7", "prop_binch3_bin0", "prop_binch3_bin1", "prop_binch3_bin2", "prop_binch3_bin3_QCD", "prop_binch3_bin3_WJetsH", "prop_binch3_bin3_tW", "prop_binch3_bin3_tbarchan", "prop_binch3_bin3_tchan", "prop_binch3_bin3_tt", "prop_binch3_bin4_QCD", "prop_binch3_bin4_WJetsH", "prop_binch3_bin4_tW", "prop_binch3_bin4_tbarchan", "prop_binch3_bin4_tchan", "prop_binch3_bin4_tt", "prop_binch3_bin5_QCD", "prop_binch3_bin5_WJetsH", "prop_binch3_bin5_tW", "prop_binch3_bin5_tbarchan", "prop_binch3_bin5_tchan", "prop_binch3_bin5_tt", "prop_binch3_bin6", "prop_binch3_bin7", "prop_binch4_bin0", "prop_binch4_bin1", "prop_binch4_bin2", "prop_binch4_bin3", "prop_binch4_bin4", "prop_binch4_bin5", "prop_binch4_bin6", "prop_binch4_bin7", "prop_binch5_bin0", "prop_binch5_bin1", "prop_binch5_bin2", "prop_binch5_bin3", "prop_binch5_bin4", "prop_binch5_bin5", "prop_binch5_bin6", "prop_binch5_bin7", "prop_binch6_bin0", "prop_binch6_bin1", "prop_binch6_bin2", "prop_binch6_bin3", "prop_binch6_bin4", "prop_binch6_bin5", "prop_binch6_bin6", "prop_binch6_bin7_QCD", "prop_binch6_bin7_WJetsH", "prop_binch6_bin7_tW", "prop_binch6_bin7_tbarchan", "prop_binch6_bin7_tchan", "prop_binch6_bin7_tt", "prop_binch7_bin0", "prop_binch7_bin1", "prop_binch7_bin2", "prop_binch7_bin3", "prop_binch7_bin4", "prop_binch7_bin5", "prop_binch7_bin6", "prop_binch7_bin7", "prop_binch8_bin0", "prop_binch8_bin1", "prop_binch8_bin2", "prop_binch8_bin3", "prop_binch8_bin4", "prop_binch8_bin5", "prop_binch8_bin6", "prop_binch8_bin7", "prop_binch9_bin0", "prop_binch9_bin1", "prop_binch9_bin2", "prop_binch9_bin3", "prop_binch9_bin4", "prop_binch9_bin5", "prop_binch9_bin6", "prop_binch9_bin7", "tW", "trigSF", "tt", "QCD_shape_e_2j1b", "QCD_shape_e_3j1b", "QCD_shape_e_3j2b", "QCD_shape_m_2j1b", "QCD_shape_m_3j1b", "QCD_shape_m_3j2b"
	
	


    
    #e_plus m_plus (top)
	#"AbsMPF", "AbsScale", "AbsStat", "DY", "ECAL", "FSR", "Flavor", "Frag", "HCAL", "ISR", "JER", "L3Res", "MC", "MET_UnclusteredEnergy", "PDF1", "PDF10", "PDF100", "PDF101", "PDF102", "PDF11", "PDF12", "PDF13", "PDF14", "PDF15", "PDF16", "PDF17", "PDF18", "PDF19", "PDF2", "PDF20", "PDF21", "PDF22", "PDF23", "PDF24", "PDF25", "PDF26", "PDF27", "PDF28", "PDF29", "PDF3", "PDF30", "PDF31", "PDF32", "PDF33", "PDF34", "PDF35", "PDF36", "PDF37", "PDF38", "PDF39", "PDF4", "PDF40", "PDF41", "PDF42", "PDF43", "PDF44", "PDF45", "PDF46", "PDF47", "PDF48", "PDF49", "PDF5", "PDF50", "PDF51", "PDF52", "PDF53", "PDF54", "PDF55", "PDF56", "PDF57", "PDF58", "PDF59", "PDF6", "PDF60", "PDF61", "PDF62", "PDF63", "PDF64", "PDF65", "PDF66", "PDF67", "PDF68", "PDF69", "PDF7", "PDF70", "PDF71", "PDF72", "PDF73", "PDF74", "PDF75", "PDF76", "PDF77", "PDF78", "PDF79", "PDF8", "PDF80", "PDF81", "PDF82", "PDF83", "PDF84", "PDF85", "PDF86", "PDF87", "PDF88", "PDF89", "PDF9", "PDF90", "PDF91", "PDF92", "PDF93", "PDF94", "PDF95", "PDF96", "PDF97", "PDF98", "PDF99", "QCD", "RelBal", "RelJER", "RelPt", "RelStat", "Scales", "UE", "WJetsH", "WJetsL", "btagSF", "elecSF", "hdamp", "muonSF", "prefire", "prop_binch1_bin0", "prop_binch1_bin1", "prop_binch1_bin2", "prop_binch1_bin3", "prop_binch1_bin4", "prop_binch1_bin5", "prop_binch1_bin6", "prop_binch1_bin7", "prop_binch2_bin0", "prop_binch2_bin1", "prop_binch2_bin2", "prop_binch2_bin3", "prop_binch2_bin4", "prop_binch2_bin5", "prop_binch2_bin6", "prop_binch2_bin7", "prop_binch3_bin0", "prop_binch3_bin1", "prop_binch3_bin2", "prop_binch3_bin3", "prop_binch3_bin4", "prop_binch3_bin5", "prop_binch3_bin6", "prop_binch3_bin7_QCD", "prop_binch3_bin7_WJetsH", "prop_binch3_bin7_tW", "prop_binch3_bin7_tchan", "prop_binch3_bin7_tt", "prop_binch4_bin0", "prop_binch4_bin1", "prop_binch4_bin2", "prop_binch4_bin3", "prop_binch4_bin4", "prop_binch4_bin5", "prop_binch4_bin6", "prop_binch4_bin7", "prop_binch5_bin0", "prop_binch5_bin1", "prop_binch5_bin2", "prop_binch5_bin3", "prop_binch5_bin4", "prop_binch5_bin5", "prop_binch5_bin6", "prop_binch5_bin7", "prop_binch6_bin0", "prop_binch6_bin1", "prop_binch6_bin2", "prop_binch6_bin3", "prop_binch6_bin4", "prop_binch6_bin5", "prop_binch6_bin6", "prop_binch6_bin7", "tW", "trigSF", "tt"
	
	#Both
	#"AbsMPF", "AbsScale", "AbsStat", "DY", "ECAL", "FSR", "Flavor", "Frag", "HCAL", "ISR", "JER", "L3Res", "MC", "MET_UnclusteredEnergy", "PDF1", "PDF10", "PDF100", "PDF101", "PDF102", "PDF11", "PDF12", "PDF13", "PDF14", "PDF15", "PDF16", "PDF17", "PDF18", "PDF19", "PDF2", "PDF20", "PDF21", "PDF22", "PDF23", "PDF24", "PDF25", "PDF26", "PDF27", "PDF28", "PDF29", "PDF3", "PDF30", "PDF31", "PDF32", "PDF33", "PDF34", "PDF35", "PDF36", "PDF37", "PDF38", "PDF39", "PDF4", "PDF40", "PDF41", "PDF42", "PDF43", "PDF44", "PDF45", "PDF46", "PDF47", "PDF48", "PDF49", "PDF5", "PDF50", "PDF51", "PDF52", "PDF53", "PDF54", "PDF55", "PDF56", "PDF57", "PDF58", "PDF59", "PDF6", "PDF60", "PDF61", "PDF62", "PDF63", "PDF64", "PDF65", "PDF66", "PDF67", "PDF68", "PDF69", "PDF7", "PDF70", "PDF71", "PDF72", "PDF73", "PDF74", "PDF75", "PDF76", "PDF77", "PDF78", "PDF79", "PDF8", "PDF80", "PDF81", "PDF82", "PDF83", "PDF84", "PDF85", "PDF86", "PDF87", "PDF88", "PDF89", "PDF9", "PDF90", "PDF91", "PDF92", "PDF93", "PDF94", "PDF95", "PDF96", "PDF97", "PDF98", "PDF99", "QCD", "RelBal", "RelJER", "RelPt", "RelStat", "Scales", "UE", "WJetsH", "WJetsL", "btagSF", "elecSF", "hdamp", "muonSF", "prefire", "prop_binch10_bin0", "prop_binch10_bin1", "prop_binch10_bin2", "prop_binch10_bin3", "prop_binch10_bin4", "prop_binch10_bin5", "prop_binch10_bin6", "prop_binch10_bin7", "prop_binch11_bin0", "prop_binch11_bin1", "prop_binch11_bin2", "prop_binch11_bin3", "prop_binch11_bin4", "prop_binch11_bin5", "prop_binch11_bin6", "prop_binch11_bin7", "prop_binch12_bin0", "prop_binch12_bin1", "prop_binch12_bin2", "prop_binch12_bin3", "prop_binch12_bin4", "prop_binch12_bin5", "prop_binch12_bin6", "prop_binch12_bin7", "prop_binch1_bin0", "prop_binch1_bin1", "prop_binch1_bin2", "prop_binch1_bin3", "prop_binch1_bin4", "prop_binch1_bin5", "prop_binch1_bin6", "prop_binch1_bin7", "prop_binch2_bin0", "prop_binch2_bin1", "prop_binch2_bin2", "prop_binch2_bin3", "prop_binch2_bin4", "prop_binch2_bin5", "prop_binch2_bin6", "prop_binch2_bin7", "prop_binch3_bin0", "prop_binch3_bin1", "prop_binch3_bin2", "prop_binch3_bin3_QCD", "prop_binch3_bin3_WJetsH", "prop_binch3_bin3_tW", "prop_binch3_bin3_tchan", "prop_binch3_bin3_tt", "prop_binch3_bin4_QCD", "prop_binch3_bin4_WJetsH", "prop_binch3_bin4_tW", "prop_binch3_bin4_tchan", "prop_binch3_bin4_tt", "prop_binch3_bin5_QCD", "prop_binch3_bin5_WJetsH", "prop_binch3_bin5_tW", "prop_binch3_bin5_tchan", "prop_binch3_bin5_tt", "prop_binch3_bin6", "prop_binch3_bin7", "prop_binch4_bin0", "prop_binch4_bin1", "prop_binch4_bin2", "prop_binch4_bin3", "prop_binch4_bin4", "prop_binch4_bin5", "prop_binch4_bin6", "prop_binch4_bin7", "prop_binch5_bin0", "prop_binch5_bin1", "prop_binch5_bin2", "prop_binch5_bin3", "prop_binch5_bin4", "prop_binch5_bin5", "prop_binch5_bin6", "prop_binch5_bin7", "prop_binch6_bin0", "prop_binch6_bin1", "prop_binch6_bin2", "prop_binch6_bin3", "prop_binch6_bin4", "prop_binch6_bin5", "prop_binch6_bin6", "prop_binch6_bin7_QCD", "prop_binch6_bin7_WJetsH", "prop_binch6_bin7_tW", "prop_binch6_bin7_tchan", "prop_binch6_bin7_tt", "prop_binch7_bin0", "prop_binch7_bin1", "prop_binch7_bin2", "prop_binch7_bin3", "prop_binch7_bin4", "prop_binch7_bin5", "prop_binch7_bin6", "prop_binch7_bin7", "prop_binch8_bin0", "prop_binch8_bin1", "prop_binch8_bin2", "prop_binch8_bin3", "prop_binch8_bin4", "prop_binch8_bin5", "prop_binch8_bin6", "prop_binch8_bin7", "prop_binch9_bin0", "prop_binch9_bin1", "prop_binch9_bin2", "prop_binch9_bin3", "prop_binch9_bin4", "prop_binch9_bin5", "prop_binch9_bin6", "prop_binch9_bin7", "tW", "trigSF", "tt"






	],

 'Lumi': [
	"Lumi",
   ],
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inpath',    '-i',  metavar = 'inpath',     dest = "inpath",   required = False, default = "./temp/cards/")
    parser.add_argument('--channel',   '-ch', metavar = 'channel',   dest = "channel",   required = False, default = "e")
    parser.add_argument('--region',    '-r', metavar = 'region',     dest = "region",    required = False, default = "1j1t")
    parser.add_argument('--doDilep',   '-d', action  = "store_true", dest = "dodlp",     required = False, default = False)
    parser.add_argument('--pretend',   '-p',  action  = "store_true", dest = "pretend",  required = False, default = False)
    parser.add_argument('--nPoints',   '-P',  metavar = 'npoints',    dest = "npoints",  required = False, default = 100, type = int)
#    parser.add_argument('--queue',     '-q',  action  = "store_true", dest = "queue",    required = False, default = False)
    parser.add_argument('--ncores',    '-j',  metavar = 'ncores',     dest = "ncores",   required = False, default = None)
    parser.add_argument('--grUnc',     '-g',  action  = "store_true", dest = "grUnc",    required = False, default = False)
    parser.add_argument('--lumiSyst',  '-lS', action  = "store_true", dest = "lumisyst",    required = False, default = False)
    parser.add_argument('--useData',   '-O', action  = "store_true",  dest = "usedata", required = False, default = False)
    parser.add_argument('--step',   '-s',  dest = "step", required = False, default = 0, type = float)

    args     = parser.parse_args()
    pretend  = args.pretend
    inpath   = args.inpath
    npoints  = args.npoints
#    doBatch  = args.queue
    ncores   = args.ncores
    grUnc    = args.grUnc
    channel  = args.channel
    region   = args.region
    doDilep  = args.dodlp
    doLumiSyst = args.lumisyst
    doAsimov = not args.usedata
    step = args.step
    
    groupList= individualList
    
    if doLumiSyst:
        groupList = globalList
        groupList.remove("Lumi")
        systsGroup["systematics"].append("Lumi")
    elif grUnc:
        groupList = globalList
    
    if not doDilep:
        for iG in groupList:
            for iU in removelist:
                if iU in systsGroup[iG]:
                    systsGroup[iG].remove(iU)

    thecard = None
    # First, find card
    if "," not in region and "," not in channel and not doDilep:
        thecard = inpath + "/" + ("absu0eta" if iR in ["3j1b","3j2b"] else "MVAscore_relaxed_b10") + "_" + channel + "_" + region + ".root"
    else:
        thecard = inpath + "/" + channel.replace(",", "") + "/combcard_" + region.replace(",", "") + "dilep" * doDilep + ".root"
    
    if not os.path.isfile(thecard):
        raise RuntimeError("FATAL: card " + thecard + " not found!")
    
    outdir = inpath + "/" + channel.replace(",", "") + "/unctable_" + region.replace(",", "") + "dilep" * doDilep + ("/globallumisyst" if doLumiSyst else "/global" if grUnc else "/individual")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    print "\n> Beginning fitting..."
    for poi in POIs:
		if step == 0:
			print "\t- POI", poi
			cumulative      = [x for x in ["r_tch"] if x != poi]#,"r_tchbar"
			print "\t- Beginning base fit..."
			nomcomm = basecommand.format(algosettings = "--algo grid --points " + str(npoints),
										 queue        = "",# if not doBatch else "--job-mode slurm --task-name nominal_" + poi,
										 parallel     = "" if not ncores  else "--parallel " + str(ncores),
										 extra        = '-n nominal_{p} {card} -P {p} '.format(p = poi, card = thecard, c = ",".join(cumulative)),
										 combineargs  = combineargs,
										 asimov       = "-t -1" if doAsimov else "",)
			nomcomm = "cd " + outdir + "; " + nomcomm + "; cd -"
			print "Command:", nomcomm
			if not pretend: os.system(nomcomm)


			print "\t- Beginning workspace saving..."
			gridcomm = basecommand.format(algosettings = "--algo none",
										  queue        = "",
										  parallel     = "",
										  extra        = '-n bestfit_{p} --saveWorkspace {card} -P {p}'.format(p = poi, card = thecard),
										  combineargs  = combineargs,
										  asimov       = "-t -1" if doAsimov else "",)
			gridcomm = "cd " + outdir + "; " + gridcomm + "; cd -"
			print "Command:", gridcomm
			if not pretend: os.system(gridcomm)


			for group in groupList:
				print "\t- Fitting group", group
				cumulative += systsGroup[group]
				thecomm = basecommand.format(algosettings = "--algo grid --points " + str(npoints),
											 queue        = "",# if not doBatch else "--job-mode slurm --task-name " + group + "_" + poi,
											 parallel     = "" if not ncores  else "--parallel " + str(ncores),
											 extra        = '-P {p} -n {g}_{p} higgsCombinebestfit_{p}.MultiDimFit.mH125.root --snapshotName MultiDimFit --freezeParameters {c}'.format(p = poi, g = group, c = ",".join(cumulative)),#--snapshotName MultiDimFit
											 combineargs  = combineargs,
											 asimov       = "-t -1" if doAsimov else "",)
				thecomm = "cd " + outdir + "; " + thecomm + "; cd -"
				print "Command:", thecomm
				if not pretend: os.system(thecomm)
		if step == 1:
			# Merging --setParameters r=1 
			print "\n> Beginning merging..."
			fileList        = []
			haddcomm = 'hadd higgsCombinenominal_{p}.MultiDimFit.mH125.root higgsCombinenominal_{p}.POINTS.*.MultiDimFit.mH125.root'.format(p = poi)
			haddcomm = "cd " + outdir + "; " + haddcomm + "; cd -"
			print "Command:", haddcomm
			if not pretend: os.system(haddcomm)

			for gr in groupList:
				tmpcomm = 'hadd higgsCombine{gp}.MultiDimFit.mH125.root higgsCombine{gp}.POINTS.*.MultiDimFit.mH125.root'.format(gp = gr + '_' + poi)
				fileList.append( "'higgsCombine{gp}.MultiDimFit.mH125.root:Freeze += {g}:{i}'".format(gp = gr + '_' + poi,
																									  g  = gr,
																									  i  = groupList.index(gr)))
				tmpcomm = "cd " + outdir + "; " + tmpcomm + "; cd -"
				print "Command:", tmpcomm
				if not pretend: os.system(tmpcomm)

			plotcomm = "plot1DScan.py higgsCombinenominal_{p}.MultiDimFit.mH125.root --others {l1} --breakdown {l2},stat --POI {p} ".format(
				p  = poi,
				l1 = ' '.join(fileList),
				l2 = ','.join(groupList))
			
			plotcomm = "cd " + outdir + "; " + plotcomm + "; cd -"
			print "Command:", plotcomm
			if not pretend: os.system(plotcomm)

			thecomm = 'cd {od}; mv scan.pdf ../scan_{p}.pdf; mv scan.png ../scan_{p}.png; cd -'.format(p = ("globallumisyst_" if doLumiSyst else "global_" if grUnc else "individual_") + poi,
																									   od = outdir)
			print 'Command:', thecomm
			if not pretend: os.system(thecomm)
