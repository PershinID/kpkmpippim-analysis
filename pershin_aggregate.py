from math import sqrt

import urllib.parse
import urllib.request
import csv
import codecs
import json

from copy import deepcopy

import os

from ROOT import TFile, TGraph, TGraphErrors, TDirectory, TDirectoryFile, TSpline3

def fill_data_info_energy(data_info):
    url_base = "https://cmd.inp.nsk.su/~compton/gogs/compton/tables/raw/online/"
    for version in data_info:
        for year in data_info[version]["years"]:
            elabels = dict()
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabels[ data_info[version]["years"][year]["elabels"][elabel]["compton-energy-point"] ] = elabel
            
            url_year  = urllib.parse.urljoin(url_base, data_info[version]["years"][year]["energy-path"])
            with urllib.request.urlopen(url_year) as file_year:
                ## energy_point,first_run,last_run,mean_energy,mean_energy_stat_err,mean_energy_sys_err,mean_spread,mean_spread_stat_err
                csv_year  = csv.DictReader(codecs.iterdecode(file_year, "utf-8"))
                for row in csv_year:
                    if row['energy_point'] in elabels:
                        data_info[version]["years"][year]["elabels"][ elabels[ row['energy_point'] ] ].update({
                            "mean-energy": float(row['mean_energy']),
                            "mean-energy-stat-error": float(row['mean_energy_stat_err']),
                            "mean-energy-sys-error": float(row['mean_energy_sys_err']),
                            "mean-energy-spread": float(row['mean_spread']),
                            "mean-energy-spread-stat-error": float(row['mean_spread_stat_err']),
                            "scan-first-run": int(row['first_run']),
                            "scan-last-run": int(row['last_run']),
                        })


def fill_data_info_lum(data_info):
    for version in data_info:
        seasons = dict()
        for year in data_info[version]["years"]:
            seasons[ data_info[version]["years"][year]["lum-season-name"] ] = year
        
        path_lumin = f"/storeA/ryzhenenkov/lum_{version}.dat"
        with open(path_lumin, "r") as file_lumin:
            file_lumin.__next__()
            season_name = ''
            elabels = dict()
            for line in file_lumin:
                if not season_name:
                    ## finding "season" word
                    if "season" not in line: ## season id = 1 (HIGH 2011)
                        continue
                    
                    ## finding season name
                    season_name = line[ line.find('(') + 1 : line.rfind(')')] ## (HIGH 2011)
                    if season_name not in seasons:
                        season_name = ''
                        continue
                    if season_name in seasons and "zero version" in line: ## season id = 15 (PHI2024) zero version
                        season_name = ''
                        continue

                    ## finding energy points of this season
                    for elabel in data_info[version]["years"][ seasons[season_name] ]["elabels"]:
                        elabels[
                            data_info[version]["years"][ seasons[season_name] ]["elabels"][elabel]["lum-energy-point"]
                        ] = elabel
                    continue
                
                ## scanning until line isn't empty
                if line.isspace():
                    season_name = ''
                    elabels = dict()
                    continue
                row = line.split()
                if row[0] in elabels:
                    data_info[version]["years"][ seasons[season_name] ]["elabels"][
                        elabels[ row[0] ]
                    ]["off-lum"] = float(row[1])
                    data_info[version]["years"][ seasons[season_name] ]["elabels"][
                        elabels[ row[0] ]
                    ]["off-lum-stat-err"] = float(row[2])
                    data_info[version]["years"][ seasons[season_name] ]["elabels"][
                        elabels[ row[0] ]
                    ]["scan-tr-ph-path"] = row[3]
                    

final_cut_number = -1

## just for v9 version
def fill_data_info_sim_cuts(data_info):
    prelim_dir = "/store11/idpershin/kpkmpippim/prelim_cuts_new"
    final_dir  = "/store11/idpershin/kpkmpippim/pershin_cut_wo_nph_wo_phth"
    
    ## finding sample files
    prelim_hists_paths = filter(
        lambda path: path.find('sim') != -1,
        filter(
            lambda path: path.find('hists') != -1,
            os.listdir(prelim_dir)
        )
    )
    final_hists_paths = filter(
        lambda path: path.find('sim') != -1,
        filter(
            lambda path: path.find('hists') != -1,
            os.listdir(final_dir)
        )
    )
    if not prelim_hists_paths or not final_hists_paths: ## no hists files present
        return

    ## determining a list of used cuts
    prelim_sample_file = TFile.Open(os.path.normpath(prelim_dir + '/' + next(prelim_hists_paths)), 'read')
    final_sample_file = TFile.Open(os.path.normpath(final_dir + '/' + next(final_hists_paths)), 'read')
    prelim_cuts_dirs_names = list(map(
        lambda obj: obj.GetName(),
        filter(
            lambda obj: isinstance(obj, TDirectory),
            map(
                lambda key: prelim_sample_file.Get( key.GetName() ),
                prelim_sample_file.GetListOfKeys()
            )
        )
    ))
    prelim_cuts_dirs_names.sort()
    prelim_cuts_dirs_names.insert(0, '0_no-cut') ## prelim_cuts_dirs_names: ['0_no-cut', '0', '1_nt', '2_tnhit', ...]
    prelim_sample_file.Close()
    
    final_cuts_dirs_names = list(map(
        lambda obj: obj.GetName(),
        filter(
            lambda obj: isinstance(obj, TDirectory),
            map(
                lambda key: final_sample_file.Get( key.GetName() ),
                final_sample_file.GetListOfKeys() ## list of TKeys
            ) ## list of TObjects
        ) ## list of TDirectorys
    ))
    final_cuts_dirs_names.sort()
    final_cuts_dirs_names.pop(0)
    final_cuts_dirs_names.insert(0, '0_no-cut') ## final_cuts_dirs_names:  ['0_no-cut', '1_KpKmPipPimLklhd', ...]
    final_sample_file.Close()

    ## Creating a cuts statistics package
    for version in data_info:
        for year in data_info[version]:
            for elabel in data_info[version][year]["elabels"]:
                ## simulation selection feedback
                hists_file = TFile.Open(prelim_dir + f'/sim_prelim_hists{year}_tr_ph_fc_e{data_info[version][year]["elabels"][elabel]["scan-energy-point"]}_{version}.root', 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                prelim_entries_numbers = graph_entries_selected.GetY()
                hists_file.Close()
                prelim_cuts_statistics = list(zip(prelim_cuts_dirs_names, prelim_entries_numbers)) ## [("0_no-cut", 100000), ("0", 100000), ...]
                
                hists_file = TFile.Open(final_dir + f'/sim_hists{year}_tr_ph_fc_e{data_info[version][year]["elabels"][elabel]["scan-energy-point"]}_{version}.root', 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                final_entries_numbers = graph_entries_selected.GetY()
                hists_file.Close()
                final_cuts_statistics = list(zip(final_cuts_dirs_names, final_entries_numbers)) ## [("0_no-cut", 1000), (0_KpKmPipPimLklhd, 200), ...]
                final_cuts_statistics.pop(0)

                cuts_statistics = prelim_cuts_statistics + final_cuts_statistics
                data_info[version][year]["elabels"][elabel].update({
                    "sim-entries-selected": dict(map(
                        lambda i_name_entries: (f"{i_name_entries[0]}_{i_name_entries[1][0]}", int(i_name_entries[1][1])),
                        enumerate(filter(
                            lambda name_entries: name_entries[0] != '',
                            map(
                                lambda name_entries: (name_entries[0].lstrip('_0123456789'), name_entries[1]),
                                cuts_statistics ## [(), (), ..., (), (), ...]
                            ) ## [()]
                        ))
                    ))
                })
                data_info[version][year]["elabels"][elabel].update({
                    "selection-efficiency": cuts_statistics[final_cut_number][1] / cuts_statistics[0][1],
                    "selection-efficiency-stat-err": sqrt(cuts_statistics[final_cut_number][1]) / cuts_statistics[0][1],
                })


def fill_data_info_multihad_cuts(data_info):
    prelim_dir = "/store11/idpershin/kpkmpippim/prelim_cuts_new"
    final_dir  = "/store11/idpershin/kpkmpippim/pershin_cut_wo_nph_wo_phth"
    
    ## finding sample files
    prelim_hists_paths = filter(
        lambda path: path.find('multihad') != -1,
        filter(
            lambda path: path.find('hists') != -1,
            os.listdir(prelim_dir)
        )
    )
    final_hists_paths = filter(
        lambda path: path.find('multihad') != -1,
        filter(
            lambda path: path.find('hists') != -1,
            os.listdir(final_dir)
        )
    )
    if not prelim_hists_paths or not final_hists_paths: ## no hists files present
        return

    ## determining a list of used cuts
    prelim_sample_file = TFile.Open(os.path.normpath(prelim_dir + '/' + next(prelim_hists_paths)), 'read')
    final_sample_file = TFile.Open(os.path.normpath(final_dir + '/' + next(final_hists_paths)), 'read')
    prelim_cuts_dirs_names = list(map(
        lambda obj: obj.GetName(),
        filter(
            lambda obj: isinstance(obj, TDirectory),
            map(
                lambda key: prelim_sample_file.Get( key.GetName() ),
                prelim_sample_file.GetListOfKeys()
            )
        )
    ))
    prelim_cuts_dirs_names.sort()
    prelim_cuts_dirs_names.insert(0, '0_no-cut') ## prelim_cuts_dirs_names: ['0_no-cut', '0', '1_nt', '2_tnhit', ...]
    prelim_sample_file.Close()

    final_cuts_dirs_names = list(map(
        lambda obj: obj.GetName(),
        filter(
            lambda obj: isinstance(obj, TDirectory),
            map(
                lambda key: final_sample_file.Get( key.GetName() ),
                final_sample_file.GetListOfKeys() ## list of TKeys
            ) ## list of TObjects
        ) ## list of TDirectorys
    ))
    final_cuts_dirs_names.sort()
    final_cuts_dirs_names.pop(0)
    final_cuts_dirs_names.insert(0, '0_no-cut') ## final_cuts_dirs_names:  ['0_no-cut', '1_KpKmPipPimLklhd', ...]
    final_sample_file.Close()

    ## Creating a cuts statistics package
    for version in data_info:
        for year in data_info[version]:
            for elabel in data_info[version][year]["elabels"]:
                ## simulation selection feedback
                hists_file = TFile.Open(prelim_dir + f'/multihad_prelim_hists{year}_tr_ph_fc_e{data_info[version][year]["elabels"][elabel]["scan-energy-point"]}_{version}.root', 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                prelim_entries_numbers = graph_entries_selected.GetY()
                hists_file.Close()
                prelim_cuts_statistics = list(zip(prelim_cuts_dirs_names, prelim_entries_numbers)) ## [("0_no-cut", 100000), ("0", 100000), ...]
                prelim_cuts_statistics.pop(0)
                prelim_cuts_statistics[0] = ('0_no-cut', prelim_cuts_statistics[0][1])

                hists_file = TFile.Open(final_dir + f'/multihad_hists{year}_tr_ph_fc_e{data_info[version][year]["elabels"][elabel]["scan-energy-point"]}_{version}.root', 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                final_entries_numbers = graph_entries_selected.GetY()
                hists_file.Close()
                final_cuts_statistics = list(zip(final_cuts_dirs_names, final_entries_numbers)) ## [("0_no-cut", 1000), (0_KpKmPipPimLklhd, 200), ...]
                final_cuts_statistics.pop(0)

                cuts_statistics = prelim_cuts_statistics + final_cuts_statistics
                data_info[version][year]["elabels"][elabel].update({
                    "multihad-background-entries-selected": dict(map(
                        lambda i_name_entries: (f"{i_name_entries[0]}_{i_name_entries[1][0]}", int(i_name_entries[1][1])),
                        enumerate(filter(
                            lambda name_entries: name_entries[0] != '',
                            map(
                                lambda name_entries: (name_entries[0].lstrip('_0123456789'), name_entries[1]),
                                cuts_statistics ## [(), (), ..., (), (), ...]
                            ) ## [()]
                        ))
                    ))
                })
                data_info[version][year]["elabels"][elabel].update({
                    "simulation-background-share": cuts_statistics[final_cut_number][1] / cuts_statistics[0][1]
                })
    
            
def fill_data_info_scan_cuts(data_info):
    prelim_dir = "/store11/idpershin/kpkmpippim/prelim_cuts_new"
    final_dir  = "/store11/idpershin/kpkmpippim/pershin_cut_wo_nph_wo_phth"
    
    ## finding sample files
    prelim_hists_paths = filter(
        lambda path: path.find('sim') == -1 and path.find('multihad') == -1,
        filter(
            lambda path: path.find('hists') != -1,
            os.listdir(prelim_dir)
        )
    )
    final_hists_paths = filter(
        lambda path: path.find('sim') == -1 and path.find('multihad') == -1,
        filter(
            lambda path: path.find('hists') != -1,
            os.listdir(final_dir)
        )
    )
    if not prelim_hists_paths or not final_hists_paths: ## no hists files present
        return

    ## determining a list of used cuts
    prelim_sample_file = TFile.Open(os.path.normpath(prelim_dir + '/' + next(prelim_hists_paths)), 'read')
    final_sample_file = TFile.Open(os.path.normpath(final_dir + '/' + next(final_hists_paths)), 'read')
    prelim_cuts_dirs_names = list(map(
        lambda obj: obj.GetName(),
        filter(
            lambda obj: isinstance(obj, TDirectory),
            map(
                lambda key: prelim_sample_file.Get( key.GetName() ),
                prelim_sample_file.GetListOfKeys()
            )
        )
    ))
    prelim_cuts_dirs_names.sort()
    prelim_cuts_dirs_names.insert(0, '0_no-cut') ## prelim_cuts_dirs_names: ['0_no-cut', '0', '1_nt', '2_tnhit', ...]
    prelim_sample_file.Close()
    
    final_cuts_dirs_names = list(map(
        lambda obj: obj.GetName(),
        filter(
            lambda obj: isinstance(obj, TDirectory),
            map(
                lambda key: final_sample_file.Get( key.GetName() ),
                final_sample_file.GetListOfKeys() ## list of TKeys
            ) ## list of TObjects
        ) ## list of TDirectorys
    ))
    final_cuts_dirs_names.sort()
    final_cuts_dirs_names.pop(0)
    final_cuts_dirs_names.insert(0, '0_no-cut') ## final_cuts_dirs_names:  ['0_no-cut', '1_KpKmPipPimLklhd', ...]
    final_sample_file.Close()

    ## Creating a cuts statistics package
    for version in data_info:
        for year in data_info[version]:
            for elabel in data_info[version][year]["elabels"]:
                ## simulation selection feedback
                hists_path = os.path.normpath(prelim_dir + f'/prelim_hists{year}_tr_ph_fc_e{data_info[version][year]["elabels"][elabel]["scan-energy-point"]}_{version}.root')
                if not os.path.exists(hists_path):
                    data_info[version][year]["elabels"][elabel]["has-results"] = False
                    continue
                hists_file = TFile.Open(hists_path, 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                prelim_entries_numbers = graph_entries_selected.GetY()
                hists_file.Close()
                prelim_cuts_statistics = list(zip(prelim_cuts_dirs_names, prelim_entries_numbers)) ## [("0_no-cut", 100000), ("0", 100000), ...]

                hists_path = os.path.normpath(final_dir + f'/hists{year}_tr_ph_fc_e{data_info[version][year]["elabels"][elabel]["scan-energy-point"]}_{version}.root')
                if not os.path.exists(hists_path):
                    data_info[version][year]["elabels"][elabel]["has-results"] = False
                    continue
                hists_file = TFile.Open(hists_path, 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                final_entries_numbers = graph_entries_selected.GetY()
                hists_file.Close()
                final_cuts_statistics = list(zip(final_cuts_dirs_names, final_entries_numbers)) ## [("0_no-cut", 1000), (0_KpKmPipPimLklhd, 200), ...]
                final_cuts_statistics.pop(0)

                cuts_statistics = prelim_cuts_statistics + final_cuts_statistics
                data_info[version][year]["elabels"][elabel].update({
                    "scan-entries-selected": dict(map(
                        lambda i_name_entries: (f"{i_name_entries[0]}_{i_name_entries[1][0]}", int(i_name_entries[1][1])),
                        enumerate(filter(
                            lambda name_entries: name_entries[0] != '',
                            map(
                                lambda name_entries: (name_entries[0].lstrip('_0123456789'), name_entries[1]),
                                cuts_statistics ## [(), (), ..., (), (), ...]
                            ) ## [()]
                        ))
                    ))
                })
                data_info[version][year]["elabels"][elabel].update({
                    "n-events": int(cuts_statistics[final_cut_number][1])
                })
                data_info[version][year]["elabels"][elabel]["has-results"] = True


def fill_data_info_efficiency(data_info, final_dir):
    simulation_initial_events_number = 100000
    
    ## Creating a cuts statistics package
    for version in data_info:
        for year in data_info[version]["years"]:
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabel_data = data_info[version]["years"][year]["elabels"][elabel]
                energy_point = elabel_data["scan-energy-point"]
                
                hists_path = final_dir + f'/sim_hists{year}_tr_ph_fc_e{energy_point}_{version}.root'
                if not os.path.exists(hists_path):
                    elabel_data["selection-efficiency"] = None
                    elabel_data["selection-efficiency-stat-err"] = None
                    continue
                hists_file = TFile.Open(hists_path, 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                simulation_final_events_number = int(graph_entries_selected.GetY()[-1])
                hists_file.Close()

                elabel_data["selection-efficiency"] = simulation_final_events_number / simulation_initial_events_number
                elabel_data.pop("selection-efficiency-err", None)
                elabel_data["selection-efficiency-stat-err"] = sqrt(simulation_final_events_number) / simulation_initial_events_number


def fill_data_info_events_number(data_info, final_dir):
    ## Creating a cuts statistics package
    for version in data_info:
        for year in data_info[version]["years"]:
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabel_data = data_info[version]["years"][year]["elabels"][elabel]
                energy_point = elabel_data["scan-energy-point"]
                
                hists_path = final_dir + f'/hists{year}_tr_ph_fc_e{energy_point}_{version}.root'
                if not os.path.exists(hists_path):
                    elabel_data["n-events"] = None
                    continue
                hists_file = TFile.Open(hists_path, 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                events_number = int(graph_entries_selected.GetY()[-1])
                hists_file.Close()

                elabel_data["n-events"] = events_number


def fill_data_info_background_share(data_info, final_dir):
    multihad_initial_events_number = 500000
    
    ## Creating a cuts statistics package
    for version in data_info:
        for year in data_info[version]["years"]:
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabel_data = data_info[version]["years"][year]["elabels"][elabel]
                energy_point = elabel_data["scan-energy-point"]
                
                hists_path = final_dir + f'/multihad_hists{year}_tr_ph_fc_e{energy_point}_{version}.root'
                if not os.path.exists(hists_path):
                    elabel_data["simulation-background-share"] = None
                    continue
                hists_file = TFile.Open(hists_path, 'read')
                graph_entries_selected = hists_file.Get('g_entries_selected')
                multihad_final_events_number = int(graph_entries_selected.GetY()[-1])
                hists_file.Close()

                elabel_data["simulation-background-share"] = multihad_final_events_number / multihad_initial_events_number

                
def calculate_visible_cross_section(data_info):
    for version in data_info:
        for year in data_info[version]["years"]:
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabel_data = data_info[version]["years"][year]["elabels"][elabel]
                
                if elabel_data["n-events"] == None or elabel_data["selection-efficiency"] == None: continue
                n_events = elabel_data["n-events"]
                eff = elabel_data["selection-efficiency"]
                eff_error = elabel_data["selection-efficiency-stat-err"]
                lumin = elabel_data["off-lum"]
                lumin_error = elabel_data["off-lum-stat-err"]
                
                vcs = n_events / (eff * lumin)
                vcs_error = sqrt(n_events / ((eff * lumin) ** 2) + (vcs * (eff_error / eff)) ** 2 + (vcs * (lumin_error / lumin)) ** 2)
                elabel_data["visible-cross-section"] = vcs
                elabel_data["visible-cross-section-stat-err"] = vcs_error


def calculate_born_cross_section(data_info):
    track_detection_efficiency_correction = 0.97355
    path_radcorr = "/spoolA/idpershin/analysis/kpkmpippim/radiative_correction.root"
    file_radcorr = TFile.Open(path_radcorr, 'read')
    spl_radcorr = TSpline3("spl_radcorr", file_radcorr.Get("g_radcorr"))
    
    for version in data_info:
        for year in data_info[version]["years"]:
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabel_data = data_info[version]["years"][year]["elabels"][elabel]
                
                if elabel_data["visible-cross-section"] == None: continue
                vcs = elabel_data["visible-cross-section"]
                vcs_error = elabel_data["visible-cross-section-stat-err"]
                elabel_data["born-cross-section"] = vcs / spl_radcorr.Eval(
                    elabel_data["mean-energy"] * 2.e-3
                ) / track_detection_efficiency_correction
                elabel_data["born-cross-section-stat-err"] = vcs_error / spl_radcorr.Eval(
                    elabel_data["mean-energy"] * 2.e-3
                ) / track_detection_efficiency_correction

    file_radcorr.Close()
                

def draw_born_cross_section(data_info):
    bcs_path = "/spoolA/idpershin/analysis/kpkmpippim/born_cross_section_pershin.root"

    bcs_lists = {}
    for version in data_info:
        for year in data_info[version]["years"]:
            bcs_lists[year] = []
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabel_data = data_info[version]["years"][year]["elabels"][elabel]
                
                if elabel_data["born-cross-section"] == None: continue
                bcs_lists[year].append((
                    elabel_data["mean-energy"],
                    elabel_data["mean-energy-spread"],
                    elabel_data["born-cross-section"],
                    elabel_data["born-cross-section-stat-err"],
                ))
            bcs_lists[year].sort(key = lambda item: item[0])
    
    file_bcs = TFile.Open(bcs_path, 'recreate')
    graphs_bcs = []

    for version in data_info:
        for year in data_info[version]["years"]:
            graphs_bcs.append(TGraphErrors())
            graphs_bcs[-1].SetName(f"g_bcs_{year}")
            graphs_bcs[-1].SetTitle(f"g_bcs_{year}")
            for i, (energy, energy_spread, bcs, bcs_err) in enumerate(bcs_lists[year]):
                graphs_bcs[-1].SetPoint(i, energy * 2.e-3, bcs)
                graphs_bcs[-1].SetPointError(i, energy_spread * 2.e-3, bcs_err)

    for graph in graphs_bcs:
        graph.Write()
    file_bcs.Save()
    file_bcs.Close()
    
    
def draw_efficiency(data_info):
    eff_path = "/spoolA/idpershin/analysis/kpkmpippim/efficiency_pershin.root"

    eff_list = []
    for version in data_info:
        for year in data_info[version]["years"]:
            for elabel in data_info[version]["years"][year]["elabels"]:
                elabel_data = data_info[version]["years"][year]["elabels"][elabel]
                
                eff_list.append((
                    elabel_data["mean-energy"],
                    elabel_data["mean-energy-spread"],
                    elabel_data["selection-efficiency"],
                    elabel_data["selection-efficiency-stat-err"],
                ))
    eff_list.sort(key = lambda item: item[0])

    file_eff = TFile.Open(eff_path, 'recreate')
    g_eff = TGraphErrors()
    g_eff.SetName("g_eff")
    g_eff.SetTitle("Efficiency")
    g_eff.GetXaxis().SetTitle("E_{c.m.}, Mev")
    g_eff.GetYaxis().SetTitle("Efficiency, %")

    for i, (energy, energy_spread, eff, eff_err) in enumerate(eff_list):
        g_eff.SetPoint(i, energy * 2.e-3, eff * 100)
        g_eff.SetPointError(i, energy_spread * 2.e-3, eff_err * 100)
    
    g_eff.Write()
    file_eff.Save()
    file_eff.Close()


def draw_rad_corr(data_other):
    path_radcorr = "/spoolA/idpershin/analysis/kpkmpippim/radiative_correction.root"
    file_radcorr = TFile.Open(path_radcorr, 'recreate')

    g_radcorr = TGraph()
    g_radcorr.SetName("g_radcorr")
    g_radcorr.SetTitle("Radiative corrections")
    for i, elabel in enumerate(data_other["cmd-3"]):
        g_radcorr.SetPoint(i,
            data_other["cmd-3"][elabel]["mean-energy"],
            data_other["cmd-3"][elabel]["radiative-correction"],
        )
    g_radcorr.Write()
    file_radcorr.Save()
    file_radcorr.Close()
    

def draw_cross_section_other_exps(data_other):
    vcs_path = "/spoolA/idpershin/analysis/kpkmpippim/cross_section_data_others.root"
    
    file_vcs = TFile.Open(vcs_path, 'recreate')
    
    g_babar = TGraphErrors()
    g_babar.SetName("g_babar")
    g_babar.SetTitle("BaBar cross section")
    for i, elabel in enumerate(data_other["babar"]):
        g_babar.SetPoint(i,
            data_other["babar"][elabel]["mean-energy"],
            data_other["babar"][elabel]["cross-section"],
        )
        g_babar.SetPointError(i,
            data_other["babar"][elabel]["mean-energy-spread"],
            data_other["babar"][elabel]["cross-section-stat-error"],
        )
    g_babar.Write()

    g_cmd = TGraphErrors()
    g_cmd.SetName("g_cmd")
    g_cmd.SetTitle("CMD-3 cross section (2016)")
    for i, elabel in enumerate(data_other["cmd-3"]):
        g_cmd.SetPoint(i,
            data_other["cmd-3"][elabel]["mean-energy"],
            data_other["cmd-3"][elabel]["cross-section"],
        )
        g_cmd.SetPointError(i,
            data_other["cmd-3"][elabel]["mean-energy-spread"],
            data_other["cmd-3"][elabel]["cross-section-stat-error"],
        )
    
    g_cmd.Write()
    file_vcs.Save()
    file_vcs.Close()
    
    
def sum_number_of_events(data_info):
    n_full_by_year = []
    for version in data_info:
        for year in data_info[version]["years"]:
            n_full_by_year.append(0)
            for elabel in data_info[version]["years"][year]["elabels"]:
                n_full_by_year[-1] += data_info[version]["years"][year]["elabels"][elabel].get("n-events", 0)
            data_info[version]["years"][year]["n-events"] = n_full_by_year[-1]
        data_info[version]["n-events"] = sum(n_full_by_year)
    
    
if __name__ == "__main__":
    path_info = "/spoolA/idpershin/analysis/kpkmpippim/data_info_pershin.json"
    final_dir  = "/store11/idpershin/kpkmpippim/pershin_cut_w_TotalP-DeltaE-2"
#    path_other = "cross_section_data_others.json"
    
    with open(path_info, 'r') as file_info:
        json_info = json.load(file_info)
#    with open(path_other, 'r') as file_other:
#        json_other = json.load(file_other)

#    fill_data_info_sim_cuts(json_info)
#    fill_data_info_multihad_cuts(json_info)
#    fill_data_info_scan_cuts(json_info)
#    draw_efficiency(json_info)
#    draw_rad_corr(json_other)
#    draw_visible_cross_section(json_info)
#    draw_cross_section_other_exps(json_other)
    fill_data_info_efficiency(json_info, final_dir)
    fill_data_info_events_number(json_info, final_dir)
    fill_data_info_background_share(json_info, final_dir)
    calculate_visible_cross_section(json_info)
    calculate_born_cross_section(json_info)
    draw_born_cross_section(json_info)
    sum_number_of_events(json_info)
    
    with open(path_info, 'w') as file_info:
        json.dump(json_info, file_info, indent=4)
