#!/bin/bash

#~~~~~~~~~~~~for PbPb~~~~~~~~~~~~~
root -l -b -q 'Tree_Analyzer.C("inputfile_PbPbMC_Tracks_Skims.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/quick_Check_Files", "PbPb", 1)'
#root -l -b -q 'Tree_Analyzer.C("inputfile_PbPbData_Tracks_Skims.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/quick_Check_Files", "PbPb", 0)'
#root -l -b -q 'Tree_Analyzer.C("inputfile_PbPbMC_Tracks.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/quick_Check_Files", "PbPb", 1)'
#root -l -b -q 'Tree_Analyzer.C("inputfile_PbPbMC_Tracks_Jussi.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/quick_Check_Files", "PbPb", 1)'
#root -l -b -q 'Tree_Analyzer.C("inputfile_PbPbData_Tracks.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/quick_Check_Files", "PbPb", 0)'

#~~~~~~~~~~~~for pp~~~~~~~~~~~~~
#root -l -b -q 'Tree_Analyzer.C("inputfile_ppMC_Jussi.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/quick_Check_Files", "pp", 1)'
#root -l -b -q 'Tree_Analyzer.C("inputfile_ppData_Jussi.txt", 0, "/eos/user/r/rpradhan/Research/ATUIC/Check_RootFile/Jussi/MC/rootfiles/quick_Check_Files", "pp", 0)'
