import ROOT
import math
import glob
import os
import argparse


def invariant_mass(p4s):
    total = ROOT.TLorentzVector()
    for p4 in p4s:
        total += p4
    return total.M()

def get_lepton_p4(pt, eta, phi, mass):
    p4 = ROOT.TLorentzVector()
    p4.SetPtEtaPhiM(pt, eta, phi, mass)
    return p4

def get_jet_p4(pt, eta, phi, mass):
    p4 = ROOT.TLorentzVector()
    p4.SetPtEtaPhiM(pt, eta, phi, mass)
    return p4

def deltaPhi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi: dphi -= 2*math.pi
    while dphi < -math.pi: dphi += 2*math.pi
    return abs(dphi)

def w_transverse_mass(pt_lep, phi_lep, met, met_phi):
    dphi = deltaPhi(phi_lep, met_phi)
    return math.sqrt(2 * pt_lep * met * (1 - math.cos(dphi)))

def mT2_higgs(pt_lep1, phi_lep1, pt_lep2, phi_lep2, met, met_phi):
    lep_p4 = get_lepton_p4(third_lepton['pt'], third_lepton['eta'], third_lepton['phi'], third_lepton['mass'])
    jet1_p4 = get_jet_p4(jets[i1]['pt'], jets[i1]['eta'], jets[i1]['phi'], jets[i1]['mass'])
    jet2_p4 = get_jet_p4(jets[i2]['pt'], jets[i2]['eta'], jets[i2]['phi'], jets[i2]['mass'])
    vis_p4 = lep_p4 + jet1_p4 + jet2_p4
    m_vis = vis_p4.M()
    pt_vis = vis_p4.Pt()
    px_vis = vis_p4.Px()
    py_vis = vis_p4.Py()
    px_miss = event.MET_pt * math.cos(event.MET_phi)
    py_miss = event.MET_pt * math.sin(event.MET_phi)
    et_vis = math.sqrt(m_vis**2 + pt_vis**2)
    et_miss = event.MET_pt
    mt2 = (et_vis + et_miss)**2 - ((px_vis + px_miss)**2 + (py_vis + py_miss)**2)
    mt_higgs = math.sqrt(mt2) if mt2 > 0 else 0.
    return mt_higgs


parser = argparse.ArgumentParser(description="Process NanoAOD ROOT files.")
parser.add_argument("--i", type=str, required=True, help="Path to the directory containing NanoAOD ROOT files.")
parser.add_argument("--f", type=str, required=True, help="Path to the plotting directory.")
args = parser.parse_args()

input_dir = args.i
save_dir = args.f
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
if not os.path.isdir(input_dir):
    raise FileNotFoundError(f"The directory '{input_dir}' does not exist.")

# Create output file and histograms
out = ROOT.TFile(str(save_dir)+"/reco_WZ_mass.root", "RECREATE")
h_leps = ROOT.TH1F("h_leps", "Reconstructed number of leptons; n_{l}; Events", 20, 0, 20)
h_jets = ROOT.TH1F("h_jets", "Reconstructed number of jets; n_{j}; Events", 20, 0, 20)
h_mZ = ROOT.TH1F("h_mZ", "Reconstructed Z mass; m_{ll} [GeV]; Events", 60, 60, 120)
h_mW_had = ROOT.TH1F("h_mW_had", "Reconstructed hadronic W mass; m_{jj} [GeV]; Events", 50, 0, 500)
h_mW_lep_T = ROOT.TH1F("h_mW_lep_T", "Reconstructed leptonic W_T mass; m_{LNu} [GeV]; Events", 50, 0, 200)
h_mH_T = ROOT.TH1F("h_mH_T", "Reconstructed H_T mass; m_{H} [GeV]; Events", 60, 0, 600)

n_Z_candidates = 0
n_W_lep_candidates = 0
n_W_had_candidates = 0
n_H_candidates = 0

root_files = glob.glob(os.path.join(input_dir, "*.root"))
for root_file_iter, root_file_name in enumerate(root_files):
    if root_file_iter % 10 == 0:
        print(f"Processing file {root_file_iter + 1}/{len(root_files)}: {root_file_name}")
    f = ROOT.TFile.Open(root_file_name)
    tree = f.Get("Events")
    for event in tree:
        # --- Lepton selection (as in your table) ---
        leptons = []
        for i in range(event.nElectron):
            # if event.Electron_pt[i] > 25 and abs(event.Electron_eta[i]) < 2.5:
            leptons.append({'pt': event.Electron_pt[i], 'eta': event.Electron_eta[i], 'phi': event.Electron_phi[i], 'mass': 0.000511, 'charge': event.Electron_charge[i], 'pdgId': 11})
        for i in range(event.nMuon):
            # if event.Muon_pt[i] > 15 and abs(event.Muon_eta[i]) < 2.4:
            leptons.append({'pt': event.Muon_pt[i], 'eta': event.Muon_eta[i], 'phi': event.Muon_phi[i], 'mass': 0.105, 'charge': event.Muon_charge[i], 'pdgId': 13})
        leptons = sorted(leptons, key=lambda x: -x['pt'])
        n_leptons = sorted([lep for lep in leptons if lep['pt'] > 10], key=lambda x: -x['pt'])
        h_leps.Fill(len(n_leptons))
        if len(leptons) < 3: continue
        if leptons[0]['pt'] < 25 or leptons[1]['pt'] < 20 or leptons[2]['pt'] < 15: continue
        if len(leptons) > 3 and leptons[3]['pt'] > 10: continue

        # Min(mll) > 12 for all lepton pairs
        pass_mll = True
        for i in range(len(leptons)):
            for j in range(i+1, len(leptons)):
                l1 = get_lepton_p4(leptons[i]['pt'], leptons[i]['eta'], leptons[i]['phi'], leptons[i]['mass'])
                l2 = get_lepton_p4(leptons[j]['pt'], leptons[j]['eta'], leptons[j]['phi'], leptons[j]['mass'])
                if (l1 + l2).M() < 12:
                    pass_mll = False
        if not pass_mll: continue
        if abs(sum([lep['charge'] for lep in leptons])) != 1: continue

        # --- Z candidate: OSSF pair with |mll - mZ| < 25 ---
        z_mass = 91.1876
        zcands = []
        for i in range(len(leptons)):
            for j in range(i+1, len(leptons)):
                if leptons[i]['pdgId'] != leptons[j]['pdgId']: continue
                if leptons[i]['charge'] * leptons[j]['charge'] > 0: continue
                l1 = get_lepton_p4(leptons[i]['pt'], leptons[i]['eta'], leptons[i]['phi'], leptons[i]['mass'])
                l2 = get_lepton_p4(leptons[j]['pt'], leptons[j]['eta'], leptons[j]['phi'], leptons[j]['mass'])
                mll = (l1 + l2).M()
                if abs(mll - z_mass) < 25:
                    zcands.append((i, j, mll))
        if len(zcands) == 0: continue
        zcand = min(zcands, key=lambda x: abs(x[2] - z_mass))
        z_leptons = [leptons[zcand[0]], leptons[zcand[1]]]
        h_mZ.Fill(zcand[2])
        n_Z_candidates += 1

        # # --- b-jet veto ---
        # has_bjet = False
        # for i in range(event.nJet):
        #     if event.Jet_pt[i] > 20 and event.Jet_btagDeepB[i] > 0.4184:
        #         has_bjet = True
        # if has_bjet: continue

        # --- Zγ veto: |m3l - mZ| > 20 GeV ---
        third_lepton = [lep for k, lep in enumerate(leptons) if k not in [zcand[0], zcand[1]]][0]
        l3_p4 = get_lepton_p4(third_lepton['pt'], third_lepton['eta'], third_lepton['phi'], third_lepton['mass'])
        z1_p4 = get_lepton_p4(z_leptons[0]['pt'], z_leptons[0]['eta'], z_leptons[0]['phi'], z_leptons[0]['mass'])
        z2_p4 = get_lepton_p4(z_leptons[1]['pt'], z_leptons[1]['eta'], z_leptons[1]['phi'], z_leptons[1]['mass'])
        m3l = (z1_p4 + z2_p4 + l3_p4).M()
        if abs(m3l - z_mass) < 20: continue

        # --- Jet selection ---
        jets = []
        for i in range(event.nJet):
            if event.Jet_pt[i] > 30 and abs(event.Jet_eta[i]) < 4.7:
                jets.append({'pt': event.Jet_pt[i], 'eta': event.Jet_eta[i], 'phi': event.Jet_phi[i], 'mass': event.Jet_mass[i]})
        h_jets.Fill(len(jets))

        # --- Signal region selection ---
        signal_region = False
        if len(jets) == 1:
            # Δφ(l + MET, j) < π/2
            l_met_px = third_lepton['pt']*math.cos(third_lepton['phi']) + event.MET_pt*math.cos(event.MET_phi)
            l_met_py = third_lepton['pt']*math.sin(third_lepton['phi']) + event.MET_pt*math.sin(event.MET_phi)
            l_met_phi = math.atan2(l_met_py, l_met_px)
            dphi = deltaPhi(l_met_phi, jets[0]['phi'])
            if dphi < math.pi/2:
                signal_region = True
        elif len(jets) >= 2:
            # Δφ(l + MET, jj) < π/2
            px_jj = jets[0]['pt']*math.cos(jets[0]['phi']) + jets[1]['pt']*math.cos(jets[1]['phi'])
            py_jj = jets[0]['pt']*math.sin(jets[0]['phi']) + jets[1]['pt']*math.sin(jets[1]['phi'])
            phi_jj = math.atan2(py_jj, px_jj)
            l_met_px = third_lepton['pt']*math.cos(third_lepton['phi']) + event.MET_pt*math.cos(event.MET_phi)
            l_met_py = third_lepton['pt']*math.sin(third_lepton['phi']) + event.MET_pt*math.sin(event.MET_phi)
            l_met_phi = math.atan2(l_met_py, l_met_px)
            dphi = deltaPhi(l_met_phi, phi_jj)
            if dphi < math.pi/2:
                signal_region = True

        if not signal_region:
            continue 

        # --- W (leptonic) reconstruction ---
        lep_p4 = get_lepton_p4(third_lepton['pt'], third_lepton['eta'], third_lepton['phi'], third_lepton['mass'])
        nu_p4 = ROOT.TLorentzVector()
        nu_p4.SetPtEtaPhiM(event.MET_pt, 0, event.MET_phi, 0)
        w_lep_mass = (lep_p4 + nu_p4).M()
        w_lep_mt = w_transverse_mass(third_lepton['pt'], third_lepton['phi'], event.MET_pt, event.MET_phi)
        h_mW_lep_T.Fill(w_lep_mt)
        n_W_lep_candidates += 1

        # --- W (hadronic) reconstruction ---
        best_wjj_mass = None
        best_wjj_pair = None
        min_dm = 999.

        had_w_p4 = None

        if len(jets) == 1:
            # 1-jet signal region: the jet is the hadronic W
            had_w_p4 = get_jet_p4(jets[0]['pt'], jets[0]['eta'], jets[0]['phi'], jets[0]['mass'])
            h_mW_had.Fill(had_w_p4.M())
            n_W_had_candidates += 1

        elif len(jets) >= 2:
            # 2-jet signal region: leading two jets are the hadronic W
            had_w_p4 = get_jet_p4(jets[0]['pt'], jets[0]['eta'], jets[0]['phi'], jets[0]['mass']) + \
                    get_jet_p4(jets[1]['pt'], jets[1]['eta'], jets[1]['phi'], jets[1]['mass'])
            h_mW_had.Fill(had_w_p4.M())
            n_W_had_candidates += 1

        if had_w_p4 is not None:
            # Visible system: lepton + hadronic W
            vis_p4 = lep_p4 + had_w_p4
            m_vis = vis_p4.M()
            pt_vis = vis_p4.Pt()
            px_vis = vis_p4.Px()
            py_vis = vis_p4.Py()
            px_miss = event.MET_pt * math.cos(event.MET_phi)
            py_miss = event.MET_pt * math.sin(event.MET_phi)
            et_vis = math.sqrt(m_vis**2 + pt_vis**2)
            et_miss = event.MET_pt
            mt2 = (et_vis + et_miss)**2 - ((px_vis + px_miss)**2 + (py_vis + py_miss)**2)
            mt_higgs = math.sqrt(mt2) if mt2 > 0 else 0.
            h_mH_T.Fill(mt_higgs)
            n_H_candidates += 1
    f.Close()

out.Write()
out.Close()

print(f"Number of Z candidates: {n_Z_candidates}")
print(f"Number of leptonic W candidates: {n_W_lep_candidates}")
print(f"Number of hadronic W candidates: {n_W_had_candidates}")
print(f"Number of H candidates: {n_H_candidates}")
print("Histograms saved to reco_WZ_mass.root")
