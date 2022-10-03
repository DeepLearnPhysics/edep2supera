#!/usr/bin/env python3
import sys
from ROOT import TChain
import supera
from edep2supera import edep2supera
import ROOT

# bbox_cfg = {
#     "BBoxSize": "[ 740,  320,  530]",
#     "VoxelSize": "[0.4,   0.4,  0.4]",
#     "BBoxBottom": "[-370, -160,  400]",

#     # BBoxSize:   [ 740,  320,  530]   # cm
#     # BBoxTop:    [ 370,  160,  930]   # cm
#     # BBoxBottom: [-370, -160,  400]   # cm
#     # # the 4mm voxel pitch is from conversation w/ Dan D.
#     # VoxelSize:  [0.4,   0.4,  0.4]     # cm

# }
# label_cfg = {
#     "LogLevel": "VERBOSE",

#     "UseSimEnergyDeposit":       "True",  # currently unused but required parameter
#     "UseSimEnergyDepositPoints": "False", # ditto

# }

# edep2supera_cfg = {
#     "LogLevel": "VERBOSE",
#     "ActiveDetectors": "[TPCActive_shape]",
# }

edep2supera_cfg = ROOT.std.map("std::string,std::string")()
edep2supera_cfg["ActiveDetectors"]="[TPCActive_shape]"
edep2supera_cfg["LogLevel"]="VERBOSE"

label_cfg = ROOT.std.map("std::string,std::string")()
label_cfg["UseSimEnergyDeposit"]="True"
label_cfg["UseSimEnergyDepositPoints"]="False"
label_cfg["LogLevel"]="VERBOSE"


bbox_cfg = ROOT.std.map("std::string,std::string")()
bbox_cfg["BBoxSize"]="[]"#"[ 740,  320,  530]"
bbox_cfg["VoxelSize"]="[0.4,   0.4,  0.4]"
bbox_cfg["BBoxBottom"]="[]"#"[-370, -160,  400]"



if len(sys.argv) < 2:
   print ('Usage: python',sys.argv[0],'[LARCV_FILE1 LARCV_FILE2 ...]')
   sys.exit(1)


# s_driver = supera.Driver()
# s_driver.ConfigureBBoxAlgorithm("BBoxInteraction", bbox_cfg)
# s_driver.ConfigureLabelAlgorithm("LArTPCMLReco3D", label_cfg)

e2s_driver = edep2supera.SuperaDriver()
e2s_driver.Configure("edep2supera",edep2supera_cfg)
e2s_driver.ConfigureBBoxAlgorithm("BBoxInteraction", bbox_cfg)
e2s_driver.ConfigureLabelAlgorithm("LArTPCMLReco3D", label_cfg)


#ch = TChain('EDepSimEvents')
#print('Adding input:', sys.argv[1])
#ch.AddFile(sys.argv[1])
#print("Chain has", ch.GetEntries(), "entries")
#sys.stdout.flush()

detect = ROOT.TFile.Open(sys.argv[1], "READ") 
events = detect.Get("EDepSimEvents")

for ev in range(events):
   #print 
   #print("considering event:", entry)
   #sys.stdout.flush()
   #bytes = ch.GetEntry(entry)
   #if bytes < 1: break
   #ev = ch.Event
   EventInput=e2s_driver.ReadEvent(ev)
   e2s_driver.GenerateImageMeta(EventInput)
   print("meta created")
   e2s_driver.VoxelizeEvent(ev, e2s_driver.Meta(), EventInput)
   print("voxelized")
   e2s_driver.GenerateLabel(EventInput)




