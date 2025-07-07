from datetime import date

import json

from argparse import ArgumentParser

import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, wait

from Analyses.KinfitAnalysis import KinfitAnalysis

def process_single(version, year, energy_point, is_sim, is_multihad, kf_type):
    input_dir = "/store11/idpershin/kpkmpippim/prelim_cuts_new"
    output_dir = "/store11/idpershin/kpkmpippim/kinfit"
    if not os.path.exists(input_dir): raise OSError(f"Input directory does not exist: {input_dir}")
    if not os.path.exists(output_dir): raise OSError(f"Output directory does not exist: {output_dir}")

    prefix = ''
    if is_sim:
        prefix = 'sim_'
    elif is_multihad:
        prefix = 'multihad_'

    input_path = f"{input_dir}/{prefix}prelim_cut{year}_tr_ph_fc_e{energy_point}_{version}.root"
    kinfit_2k2pi_path = f"{output_dir}/{prefix}kinfit_2k2pi_{kf_type.lower()}{year}_tr_ph_fc_e{energy_point}_{version}.root"
    kinfit_4pi_path = f"{output_dir}/{prefix}kinfit_4pi{year}_tr_ph_fc_e{energy_point}_{version}.root"
    hists_path = f"{output_dir}/{prefix}hists_{kf_type.lower()}{year}_tr_ph_fc_e{energy_point}_{version}.root"
    log_path = f"{output_dir}/{prefix}{kf_type.lower()}_tr_ph_fc_y{year}_e{energy_point}_{version}_%s.log" % date.today().isoformat()
    if not os.path.exists(input_path): raise OSError(f"Input path does not exist: {input_path}")

    result = subprocess.run(["kpkmpippim-kinfit-2chk2chpi-exe", "--input-path", input_path, "--output-path", kinfit_2k2pi_path, "--kinfit-type", kf_type])
    if result.returncode != 0: raise RuntimeError(f"kpkmpippim-kinfit-2chk2chpi-exe failed with code {result.returncode}")
    result = subprocess.run(["kpkmpippim-kinfit-4chpi-exe", "--input-path", input_path, "--output-path", kinfit_4pi_path])
    if result.returncode != 0: raise RuntimeError(f"kpkmpippim-kinfit-4chpi-exe failed with code {result.returncode}")
    
    analysis = KinfitAnalysis(
        input_path,
        kinfit_2k2pi_path,
        kinfit_4pi_path,
        hists_path,
        log_path = log_path
    )
    analysis.addHistogram('h_KpKmPipPimLklhd')
    analysis.addHistogram('h_TotalP_DeltaE')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
    analysis.addHistogram('h_PiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiMissMass2')
    analysis.addHistogram('h_PiPiPiPiMissMass2')
    analysis.addHistogram('h_KKPiPiMissMass2')
    analysis.addHistogram('h_KKMissMass')
    analysis.addHistogram('h_PiPiMissMass')
    analysis.addHistogram('h_KpKmPipPimKinfitChi2')
    analysis.addHistogram('h_PipPimPipPimKinfitChi2')
    analysis.addHistogram('h_TotalP_DeltaE_KF')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi_KF')
    analysis.addHistogram('h_PiPiPiMissMass2_KF')
    analysis.addHistogram('h_KPiPiMissMass2_KF')
    analysis.addHistogram('h_PiPiPiPiMissMass2_KF')
    analysis.addHistogram('h_KKPiPiMissMass2_KF')
    analysis.addHistogram('h_KKMissMass_KF')
    analysis.addHistogram('h_PiPiMissMass_KF')
    if is_multihad: analysis.addHistogram('h_finalstate_id')
    analysis.loop()
    
    analysis.addCut('KpKmPipPimKinfitChi2')
    analysis.addHistogram('h_KpKmPipPimLklhd')
    analysis.addHistogram('h_TotalP_DeltaE')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
    analysis.addHistogram('h_PiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiMissMass2')
    analysis.addHistogram('h_PiPiPiPiMissMass2')
    analysis.addHistogram('h_KKPiPiMissMass2')
    analysis.addHistogram('h_KKMissMass')
    analysis.addHistogram('h_PiPiMissMass')
    analysis.addHistogram('h_KpKmPipPimKinfitChi2')
    analysis.addHistogram('h_PipPimPipPimKinfitChi2')
    analysis.addHistogram('h_TotalP_DeltaE_KF')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi_KF')
    analysis.addHistogram('h_PiPiPiMissMass2_KF')
    analysis.addHistogram('h_KPiPiMissMass2_KF')
    analysis.addHistogram('h_PiPiPiPiMissMass2_KF')
    analysis.addHistogram('h_KKPiPiMissMass2_KF')
    analysis.addHistogram('h_KKMissMass_KF')
    analysis.addHistogram('h_PiPiMissMass_KF')
    if is_multihad: analysis.addHistogram('h_finalstate_id')
    analysis.loop()

    analysis.addCut('PipPimPipPimKinfitChi2')
    analysis.addHistogram('h_KpKmPipPimLklhd')
    analysis.addHistogram('h_TotalP_DeltaE')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
    analysis.addHistogram('h_PiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiMissMass2')
    analysis.addHistogram('h_PiPiPiPiMissMass2')
    analysis.addHistogram('h_KKPiPiMissMass2')
    analysis.addHistogram('h_KKMissMass')
    analysis.addHistogram('h_PiPiMissMass')
    analysis.addHistogram('h_KpKmPipPimKinfitChi2')
    analysis.addHistogram('h_PipPimPipPimKinfitChi2')
    analysis.addHistogram('h_TotalP_DeltaE_KF')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi_KF')
    analysis.addHistogram('h_PiPiPiMissMass2_KF')
    analysis.addHistogram('h_KPiPiMissMass2_KF')
    analysis.addHistogram('h_PiPiPiPiMissMass2_KF')
    analysis.addHistogram('h_KKPiPiMissMass2_KF')
    analysis.addHistogram('h_KKMissMass_KF')
    analysis.addHistogram('h_PiPiMissMass_KF')
    if is_multihad: analysis.addHistogram('h_finalstate_id')
    analysis.loop()
    
    analysis.dumpToFile()
    analysis.close()


def process_all(is_sim_included = True, is_multihad_included = True, kf_type = "DCLklhd"):
    path_info = "/spoolA/idpershin/analysis/kpkmpippim/data_info_cmd3.json"
    with open(path_info, 'r') as file_info:
        json_info = json.load(file_info)

    with ProcessPoolExecutor(max_workers = 10) as executor:
        futures = []
        for version in json_info:
            for year in json_info[version]["years"]:
                for elabel in json_info[version]["years"][year]["elabels"]:
                    energy_point = json_info[version]["years"][year]["elabels"][elabel]["scan-energy-point"]

                    futures.append(
                        executor.submit(
                            process_single,
                            version,
                            year,
                            energy_point,
                            is_sim = False, is_multihad = False,
                            kf_type = kf_type,
                        )
                    )

                    if is_sim_included:
                        futures.append(
                            executor.submit(
                                process_single,
                                version,
                                year,
                                energy_point,
                                is_sim = True, is_multihad = False,
                                kf_type = kf_type,
                            )
                        )
   
                    if is_multihad_included:
                        futures.append(
                            executor.submit(
                                process_single,
                                version,
                                year,
                                energy_point,
                                is_sim = False, is_multihad = True,
                                kf_type = kf_type,
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
    parser_single.add_argument('--is-multihad', action = 'store_true')
    parser_single.add_argument('--kinfit-type', choices = ("Permut", "DCLklhd"), required = True)

    parser_all = subparsers.add_parser('all', help = "Process all available energy points")
    parser_all.add_argument('--is-sim-included', action = 'store_true')
    parser_all.add_argument('--is-multihad-included', action = 'store_true')
    parser_all.add_argument('--kinfit-type', choices = ("Permut", "DCLklhd"), required = True)
    args = parser.parse_args()

    if args.mode == 'single':
        process_single(args.version, args.year, args.energy, args.is_sim, args.is_multihad, args.kinfit_type)

    if args.mode == 'all':
        process_all(args.is_sim_included, args.is_multihad_included, args.kinfit_type)
