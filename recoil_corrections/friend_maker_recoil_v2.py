#from sys import base_prefix
import ROOT
import argparse
import os
from tqdm import tqdm
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="produce friend of input ntuple with corrected met")   
    parser.add_argument('-W', '--wsdir', default='/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/final_run323778/',
                        help='directory of RooFitWorkSpace')
    parser.add_argument('-C', '--corr', default='Zmm', help='Correction Samples: Zmm, Zee, ...?')
    parser.add_argument('-P', '--process', default='DY', help='Process: DY, TTToSemi, TTTo2L2Nu, W0J, W1J or W2J')
    parser.add_argument('-F', '--finalstate', default='mm', help='Final state: mm, ee, mmet or emet')
    parser.add_argument('-I', '--inpath', default='/ceph/moh/CROWN_samples/EarlyRun3_V08/ntuples/2018/', help='path to input samples')
    parser.add_argument('-O', '--outpath', default='/ceph/jdriesch/CROWN_samples/EarlyRun3_V08/friends/2018/', help='path to output')
    parser.add_argument('--overwrite', default=False)
    args = parser.parse_args()
    return args

def load_fit_results(pathName, nBins):
    # loop over parallel and perpendicular component of hadronic recoil and load fit results

    workspace = []

    for i in ["1","2"]:
        file = ROOT.TFile(pathName+"pdfsU" + i + ".root")
        workspace.append(file.Get("pdfsU"+i))

        for j in range(nBins):
            pdf = workspace[-1].pdf("sig_" + str(j))
            myX = workspace[-1].var("u_" + str(j))
            cdf = pdf.createCdf(myX)
            workspace[-1].Import(cdf, ROOT.RooFit.RecycleConflictNodes() ,ROOT.RooFit.Silence())    

    # TODO: implement statistical uncertainty calculation?
    return workspace


def invertCdf(uP, mc_cdf, mc_x, data_cdf, data_x):
    # calculate p-value of mc distribution
    if (uP < mc_x.getMin()) or (uP > mc_x.getMax()): 
        return uP
    mc_x.setVal(uP)
    data_x.setVal(uP)
    data_cdf.getVal()
    # calculate inverse value of mc p-value
    pVal = data_cdf.findRoot(data_x, data_x.getMin(), data_x.getMax(), mc_cdf.getVal())
    data_x.setVal(pVal)
    return pVal


def prep(inpath):    
    # load mc to correct
    rdf = ROOT.RDataFrame("ntuple", inpath)
    
    # check how many leptons in final state and adjust total leptonic momentum accordingly
    finalstate = inpath.split("/")[-2]

    if finalstate == "ee" or finalstate == "mm":
        rdf = rdf.Define("lepPx", "pt_1*cos(phi_1) + pt_2*cos(phi_2)")
        rdf = rdf.Define("lepPy", "pt_1*sin(phi_1) + pt_2*sin(phi_2)")
        rdf = rdf.Define("lepPt", "sqrt(lepPx*lepPx + lepPy*lepPy)")
        rdf = rdf.Define("lepPhi", "atan2(lepPy, lepPx)")

    else:
        rdf = rdf.Define("lepPt", "pt_1")
        rdf = rdf.Define("lepPhi", "phi_1")

    rdf = rdf.Define("uPx", "met_uncorrected*cos(metphi_uncorrected) + lepPt*cos(lepPhi)") # sign is wrong but consistently -> met is correct
    rdf = rdf.Define("uPy", "met_uncorrected*sin(metphi_uncorrected) + lepPt*sin(lepPhi)")
    #rdf = rdf.Define("uPt", "sqrt(uPx*uPx + uPy*uPy)")
    #rdf = rdf.Define("uPhi", "atan2(uPy, uPx)") # angle is wrong by pi as sign is wrong
    rdf = rdf.Define("uP1", "- (uPx*cos(genbosonphi) + uPy*sin(genbosonphi))") # maybe use reco-phi instead
    rdf = rdf.Define("uP2", "uPx*sin(genbosonphi) - uPy*cos(genbosonphi)")
    
    return rdf


def calculateMET(rdf):
    rdf = rdf.Define("metPxcorr", "-lepPt*cos(lepPhi) - uP1corr*cos(genbosonphi) + uP2corr*sin(genbosonphi)") # TODO check if name already used
    rdf = rdf.Define("metPycorr", "-lepPt*sin(lepPhi) - uP1corr*sin(genbosonphi) - uP2corr*cos(genbosonphi)")

    rdf = rdf.Define("metPtcorr", "sqrt(metPxcorr*metPxcorr + metPycorr*metPycorr)")
    rdf = rdf.Define("metPhicorr", "atan2(metPycorr, metPxcorr)")

    return rdf


def main():
    args = parse_args()

    ws_dict = {"Zmm": {"data": "Zmm_data_triple/", "mc": "Zmm_sim_triple/"}, }

    zPtBinEdges = [0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000]

    samplePathMC = args.wsdir + ws_dict[args.corr]['mc']
    samplePathData = args.wsdir + ws_dict[args.corr]['data']
    samplePathToCorr = args.wsdir + ws_dict[args.corr]['mc']

    basepathin = args.inpath
    basepathout = args.outpath

    procs = {"DY": "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
               "TTToSemi": "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X",
               "TTTo2L2Nu": "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X",
               "W0J": "WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
               "W1J": "WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
               "W2J": "WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X"}
    
    nBins = len(zPtBinEdges)-1
    dataWorkspace = load_fit_results(samplePathData, nBins) 
    mcWorkspace = load_fit_results(samplePathMC, nBins)
    toCorrWorkspace = load_fit_results(samplePathToCorr, nBins)

    finalstate = args.finalstate
    proc = args.process
    inpath = basepathin + procs[proc] + "/" + finalstate + "/" + procs[proc]
    outdir = basepathout + procs[proc] + "/" + finalstate + "/"


    #for i in range(100): # loop over root files
    for i in [100,101]: # loop over root files
        infile = inpath + "_"+str(i)+".root"
        outfile= outdir + infile.split("/")[-1]
        if not os.path.exists(infile):
            print("Error: requested root file does not exist.")
            break
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            print("Created new directory: ", outdir)
            
        rdf = prep(infile)

        columns = ["uP1", "uP2", "lepPt", "lepPhi", "genbosonpt", "genbosonphi"]
        fromrdf = rdf.AsNumpy(columns = columns)
        tordf = {"uP1corr": fromrdf["uP1"], "uP2corr": fromrdf["uP2"], "lepPt": fromrdf["lepPt"], "lepPhi": fromrdf["lepPhi"], "genbosonphi": fromrdf["genbosonphi"]}
        corr = {}

        for k in range(nBins):
            corr[str(zPtBinEdges[k])+" to "+str(zPtBinEdges[k+1])] = []

        for j in tqdm(range(len(fromrdf["genbosonpt"]))):
            # Get pT bin of genboson
            kBin = nBins-1
            for k in range(nBins):
                if fromrdf["genbosonpt"][j] < zPtBinEdges[k]: # TODO check if right order
                    kBin = k-1
                    break
            sBin = str(kBin)

            # MET corrections for parallel and perpendicular recoil component
            for k in [0,1]:
                data_pdf   = dataWorkspace[k].pdf('sig_'+sBin)
                mc_pdf     = mcWorkspace[k].pdf('sig_'+sBin)
                toCorr_pdf = toCorrWorkspace[k].pdf('sig_'+sBin)

                data_cdf   = dataWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                mc_cdf     = mcWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf = toCorrWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')

                data_xuP   = dataWorkspace[k].var('u_'+sBin)
                mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                toCorr_xuP = toCorrWorkspace[k].var('u_'+sBin)

                mc_uPValZlike = invertCdf(fromrdf["uP"+str(k+1)][j], toCorr_cdf, toCorr_xuP, mc_cdf, mc_xuP) # TODO check order of matrix
                data_uPValZlike = invertCdf(mc_uPValZlike, mc_cdf, mc_xuP, data_cdf, data_xuP)

                tordf["uP"+str(k+1)+"corr"][j] += data_uPValZlike - mc_uPValZlike

                """
                if k:
                    corr[str(zPtBinEdges[kBin])+" to "+str(zPtBinEdges[kBin+1])].append(data_uPValZlike - mc_uPValZlike)
                """
        """ 
        for k in range(nBins):
            mean = np.mean(corr[str(zPtBinEdges[k])+" to "+str(zPtBinEdges[k+1])])
            std = np.std(corr[str(zPtBinEdges[k])+" to "+str(zPtBinEdges[k+1])])
            print("recoil correction for bin "+str(zPtBinEdges[k])+" to "+str(zPtBinEdges[k+1])+": ", mean, "+/-", std)
        """
        rdf_tosave = ROOT.RDF.MakeNumpyDataFrame(tordf)
        rdf_tosave = calculateMET(rdf_tosave)     

        if not os.path.exists(outfile) or args.overwrite:
            rdf_tosave.Snapshot("ntuple", outfile, ["uP1corr", "uP2corr", "metPtcorr", "metPhicorr"])

        if i%10==0: 
            print("done {} of {}".format(i, 100))


    print("great success")


if __name__=='__main__':
    main()