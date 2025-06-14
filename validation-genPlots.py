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

def solid_angle(v1, v2):
    return v1.Angle(v2.Vect())


h_wmass = ROOT.TH1F("h_wmass", "Gen W mass;Mass [GeV];Events", 100, 0, 200)
h_zmass = ROOT.TH1F("h_zmass", "Gen Z mass;Mass [GeV];Events", 100, 0, 200)
h_angle_WW = ROOT.TH1F("h_angle_WW", "Solid angle between Ws;Angle [rad];Events", 64, 0, math.pi)
h_angle_ZH = ROOT.TH1F("h_angle_ZH", "Solid angle between Z and H;Angle [rad];Events", 64, 0, math.pi)

n_H = 0
n_W = 0
n_WfromH = 0
n_W_leptonic = 0
n_W_hadronic = 0
n_Z = 0
n_Z_leptonic = 0
n_Z_hadronic = 0

FROM_HARD_PROCESS = 8
IS_FIRST_COPY = 12

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
    file = ROOT.TFile.Open(root_file_name)
    tree = file.Get("Events")
    for event in tree:
        nGenPart = event.nGenPart
        pdgId = list(event.GenPart_pdgId)
        statusFlags = list(event.GenPart_statusFlags)
        pt = list(event.GenPart_pt)
        eta = list(event.GenPart_eta)
        phi = list(event.GenPart_phi)
        mass = list(event.GenPart_mass)
        motherIdx = list(event.GenPart_genPartIdxMother)

        w_candidates = []
        z_candidates = []
        h_candidates = []
        all_particles = []

        for i in range(nGenPart):
            flag = statusFlags[i]
            abs_pdg = abs(pdgId[i])
            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i])
            all_particles.append((i, pdgId[i], motherIdx[i], p4))

            if not (flag & (1 << FROM_HARD_PROCESS)):
                continue
            if not (flag & (1 << IS_FIRST_COPY)):
                continue

            if abs_pdg == 25:
                h_candidates.append((i, p4))
                n_H += 1
            elif abs_pdg == 24:
                w_candidates.append((i, p4, motherIdx[i], pdgId[i]))
                n_W += 1
            elif abs_pdg == 23:
                z_candidates.append((i, p4, motherIdx[i]))
                n_Z += 1

        # ----------- Higgs -> WW -----------
        for i, w1, m1, pdg1 in w_candidates:
            for j, w2, m2, pdg2 in w_candidates:
                if i >= j: continue
                if m1 == m2 and m1 >= 0 and abs(pdgId[m1]) == 25:
                    # Fill mass and angle
                    h_wmass.Fill(invariant_mass([w1, w2]))
                    n_WfromH += 1
                    h_angle_WW.Fill(solid_angle(w1, w2))

                    # For each W, check decay mode
                    for w_idx, w_p4, w_mother, w_pdg in [(i, w1, m1, pdg1), (j, w2, m2, pdg2)]:
                        # Find daughters of this W
                        daughters = [p for p in all_particles if p[2] == w_idx]
                        # Leptonic: has e, mu, tau daughter
                        if any(abs(d[1]) in [11,13,15] for d in daughters):
                            n_W_leptonic += 1
                        # Hadronic: has quark daughter
                        elif any(abs(d[1]) in range(1,7) for d in daughters):
                            n_W_hadronic += 1
                    break  # Only one unique pair per event

        # ----------- Z -> ll or hadrons -----------
        for i, z, m1 in z_candidates:
            # Find all daughters of this Z
            daughters = [p for p in all_particles if p[2] == i]
            # Leptonic: at least one lepton daughter
            if any(abs(d[1]) in [11,13,15] for d in daughters):
                n_Z_leptonic += 1
            # Hadronic: at least one quark daughter
            elif any(abs(d[1]) in range(1,7) for d in daughters):
                n_Z_hadronic += 1

            # Try to find the Higgs in the same event for Z-H angle
            if h_candidates:
                h_p4 = h_candidates[0][1]
                h_angle_ZH.Fill(solid_angle(z, h_p4))

            # For mass plot, look for lepton pairs from this Z
            z_leptons = [d[3] for d in daughters if abs(d[1]) in [11,13,15]]
            if len(z_leptons) == 2:
                h_zmass.Fill(invariant_mass(z_leptons))
    file.Close()
out = ROOT.TFile("gen_WZ_mass.root", "RECREATE")
h_wmass.Write()
h_zmass.Write()
h_angle_WW.Write()
h_angle_ZH.Write()
out.Close()

print(f'Total Higgs candidates: {n_H}, total W candidates: {n_W}')
print(f'Total W candidates from H: {n_WfromH}\n Num of leptonic W decays: {n_W_leptonic}\n Num of hadronic W decays: {n_W_hadronic}')
print(f'Total Z candidates: {n_Z}\n Num of leptonic Z decays: {n_Z_leptonic}\n Num of hadronic Z decays: {n_Z_hadronic}')
print("Validation plots saved to validation_plots.root")
