from datetime import date

import json

from argparse import ArgumentParser

import os
from concurrent.futures import ProcessPoolExecutor, wait

from Analyses.DynamicsAnalysis import DynamicsAnalysis

def process_single(version, year, energy_point, is_sim):
    input_dir = "/store11/idpershin/kpkmpippim/final_cuts_new"
    hists_dir = "/store11/idpershin/kpkmpippim/dynamics/pershin"
    if not os.path.exists(input_dir): raise OSError(f"Input directory does not exist: {input_dir}")
    if not os.path.exists(hists_dir): raise OSError(f"Histograms directory does not exist: {hists_dir}")

    prefix = ''
    if is_sim:
        prefix = 'sim_'

    input_path = f"{input_dir}/{prefix}final_cut{year}_tr_ph_fc_e{energy_point}_{version}.root"
    hists_path = f"{hists_dir}/{prefix}dynamics{year}_tr_ph_fc_e{energy_point}_{version}.root"
    if not os.path.exists(input_path): raise OSError(f"Input path does not exist: {input_path}")

    dynamics_analysis = DynamicsAnalysis(
        input_path, hists_path
    )
    dynamics_analysis.addHistogram('h_KpPimInvarMass_KmPipInvarMass')
    dynamics_analysis.addHistogram('h_KpKmInvarMass')
    dynamics_analysis.addHistogram('h_PipPimInvarMass')
    dynamics_analysis.addHistogram('h_PipPimInvarMass_KpKmInvarMass')
    dynamics_analysis.addHistogram('h_KpPimAngle_KmPipAngle')
    dynamics_analysis.addHistogram('h_KpKmAngle')
    dynamics_analysis.addHistogram('h_PipPimAngle')
    dynamics_analysis.loop()
    dynamics_analysis.dumpToFile()
    dynamics_analysis.close()


def process_all(is_sim_included = True):
    path_info = "/spoolA/idpershin/analysis/kpkmpippim/data_info_cmd3.json"
    with open(path_info, 'r') as file_info:
        json_info = json.load(file_info)
        
    with ProcessPoolExecutor(max_workers = 8) as executor:
        futures = []
        for version in json_info:
            for year in json_info[version]:
                for elabel in json_info[version][year]["elabels"]:
                    energy_point = json_info[version][year]["elabels"][elabel]["scan-energy-point"]

                    futures.append(
                        executor.submit(
                            process_single,
                            version,
                            year,
                            energy_point,
                            is_sim = False
                        )
                    )

                    if is_sim_included:
                        futures.append(
                            executor.submit(
                                process_single,
                                version,
                                year,
                                energy_point,
                                is_sim = True
                            )
                        )
        wait(futures)

    
if __name__ == '__main__':
    ##Parsing input arguments
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest = 'mode')
    
    parser_single = subparsers.add_parser('single', help = "Process one energy point")
    parser_single.add_argument('--version', default = 'v9', choices = ['v9'], help = 'Version of CMD-3 data tree')
    parser_single.add_argument('--year', choices = ['2019', '2020', '2021', '2022', '2023'], required = True)
    parser_single.add_argument('--energy', required = True)
    parser_single.add_argument('--is-sim', action = 'store_true')

    parser_all = subparsers.add_parser('all', help = "Process all available energy points")
    parser_all.add_argument('--is-sim-included', action = 'store_true')
    args = parser.parse_args()

    if args.mode == 'single':
        process_single(args.version, args.year, args.energy, args.is_sim)

    if args.mode == 'all':
        process_all(args.is_sim_included)
