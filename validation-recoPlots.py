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

# Create output file and histograms
out = ROOT.TFile("reco_WZ_mass.root", "RECREATE")
h_mZ = ROOT.TH1F("h_mZ", "Reconstructed Z mass; m_{ll} [GeV]; Events", 60, 60, 120)
h_mW_lep = ROOT.TH1F("h_mW_lep", "Reconstructed leptonic W mass; m_{l#nu} [GeV]; Events", 60, 0, 200)
h_leps = ROOT.TH1F("h_leps", "Reconstructed number of leptons; n_{l}; Events", 60, 0, 200)
h_jets = ROOT.TH1F("h_jets", "Reconstructed number of jetss; n_{j}; Events", 60, 0, 200)
h_mW_had = ROOT.TH1F("h_mW_had", "Reconstructed hadronic W mass; m_{jj} [GeV]; Events", 60, 0, 200)

n_W_leptonic = 0
n_W_hadronic = 0
n_Z_leptonic = 0
n_Z_hadronic = 0

parser = argparse.ArgumentParser(description="Process NanoAOD ROOT files.")
parser.add_argument("--i", type=str, required=True, help="Path to the directory containing NanoAOD ROOT files.")
args = parser.parse_args()

input_dir = args.i
if not os.path.isdir(input_dir):
    raise FileNotFoundError(f"The directory '{input_dir}' does not exist.")

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
        n_Z_leptonic += 1

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
            continue  # Only fill W mass histograms for signal region

        # --- W (leptonic) reconstruction ---
        lep_p4 = get_lepton_p4(third_lepton['pt'], third_lepton['eta'], third_lepton['phi'], third_lepton['mass'])
        nu_p4 = ROOT.TLorentzVector()
        nu_p4.SetPtEtaPhiM(event.MET_pt, 0, event.MET_phi, 0)
        w_lep_mass = (lep_p4 + nu_p4).M()
        h_mW_lep.Fill(w_lep_mass)
        n_W_leptonic += 1

        # --- W (hadronic) reconstruction ---
        if len(jets) >= 2:
            best_wjj = None
            best_dm = 999
            for i in range(len(jets)):
                for j in range(i+1, len(jets)):
                    j1 = get_jet_p4(jets[i]['pt'], jets[i]['eta'], jets[i]['phi'], jets[i]['mass'])
                    j2 = get_jet_p4(jets[j]['pt'], jets[j]['eta'], jets[j]['phi'], jets[j]['mass'])
                    m_jj = (j1 + j2).M()
                    if abs(m_jj - 80.4) < best_dm:
                        best_dm = abs(m_jj - 80.4)
                        best_wjj = m_jj
            if best_wjj:
                h_mW_had.Fill(best_wjj)
                n_W_hadronic += 1
    f.Close()

h_mW_combined = h_mW_lep.Clone("h_mW_combined")
h_mW_combined.SetTitle("Combined W mass (leptonic+hadronic); m_{W} [GeV]; Events")
h_mW_combined.Add(h_mW_had)

out.WriteTObject(h_mW_combined)
out.Write()
out.Close()

print(f"Number of leptonically decaying Z's: {n_Z_leptonic}")
print(f"Number of leptonically decaying W's: {n_W_leptonic}")
print(f"Number of hadronically decaying W's: {n_W_hadronic}")
print("Histograms saved to reco_WZ_mass.root")
