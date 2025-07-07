from datetime import date

import json

from argparse import ArgumentParser

import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, wait

from Analyses.IntermediateAnalysis import IntermediateAnalysis

def process_single(version, year, energy_point, is_sim, is_multihad):
    input_dir = "/home/idpershin/analysis/kpkmpippim/prelim_cuts_new"
    output_dir = "/home/idpershin/analysis/kpkmpippim/pershin_cut_w_TotalP-DeltaE-2"

    prefix = ''
    if is_sim:
        prefix = 'sim_'
    elif is_multihad:
        prefix = 'multihad_'

    if is_multihad:
        final_states = [12, 15, 29, 32] # 12 is for K+K-pi+pi-, 15 for K+KSpi-pi0, 29 for K+K-eta, 32 for K+K-omega

        input_dir = "/home/idpershin/analysis/kpkmpippim/prelim_cuts_new"
            
        for final_state in final_states:
            input_path = f"{input_dir}/{prefix}prelim_cut{year}_tr_ph_fc_e{energy_point}_{version}.root"
            hists_path = f"{input_dir}/{prefix}fs{final_state}_hists{year}_tr_ph_fc_e{energy_point}_{version}.root"
            output_path = f"{input_dir}/{prefix}fs{final_state}_cut{year}_tr_ph_fc_e{energy_point}_{version}.root"

            analysis = IntermediateAnalysis(
                input_path,
                analysis_path = hists_path,
                output_path = output_path,
            )
            analysis.addCutNew(
                f'finalstate-id_{final_state}',
                lambda: analysis.Variables["finalstate_id"].Content != final_state
            )
            analysis.addHistogram('h_tptot_tdedx')
            analysis.addHistogram('h_phen')
            analysis.addHistogram('h_phth')
            analysis.addHistogram('h_KpKmPipPimLklhd')
            analysis.addHistogram('h_TotalP_DeltaE')
            analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
            analysis.addHistogram('h_PiPiPiMissMass2')
            analysis.addHistogram('h_KPiPiMissMass2')
            analysis.addHistogram('h_PiPiPiPiMissMass2')
            analysis.addHistogram('h_KPiPiPiMissMass2')
            analysis.addHistogram('h_KKPiPiMissMass2')
            analysis.addHistogram('h_KKMissMass')
            analysis.addHistogram('h_PiPiMissMass')
            analysis.loop()
            analysis.dumpToFile()
            analysis.close()

    input_path = f"{input_dir}/{prefix}prelim_cut{year}_tr_ph_fc_e{energy_point}_{version}.root"
    hists_path = f"{output_dir}/{prefix}hists{year}_tr_ph_fc_e{energy_point}_{version}.root"
    output_path = f"{output_dir}/{prefix}cut{year}_tr_ph_fc_e{energy_point}_{version}.root"
    if not os.path.exists(input_path): raise OSError(f"Input path does not exist: {input_path}")
            
    analysis = IntermediateAnalysis(
        input_path,
        analysis_path = hists_path,
        output_path = output_path,
    )

    analysis.addHistogram('h_tptot_tdedx')
    analysis.addHistogram('h_phen')
    analysis.addHistogram('h_phth')
    analysis.addHistogram('h_KpKmPipPimLklhd')
    analysis.addHistogram('h_TotalP_DeltaE')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
    analysis.addHistogram('h_PiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiMissMass2')
    analysis.addHistogram('h_PiPiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiPiMissMass2')
    analysis.addHistogram('h_KKPiPiMissMass2')
    analysis.addHistogram('h_KKMissMass')
    analysis.addHistogram('h_PiPiMissMass')
    if is_multihad: analysis.addHistogram('h_finalstate_id')
    analysis.loop()
            
    analysis.addCut('KpKmPipPimLklhd')
    analysis.addHistogram('h_tptot_tdedx')
    analysis.addHistogram('h_phen')
    analysis.addHistogram('h_phth')
    analysis.addHistogram('h_KpKmPipPimLklhd')
    analysis.addHistogram('h_TotalP_DeltaE')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
    analysis.addHistogram('h_PiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiMissMass2')
    analysis.addHistogram('h_PiPiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiPiMissMass2')
    analysis.addHistogram('h_KKPiPiMissMass2')
    analysis.addHistogram('h_KKMissMass')
    analysis.addHistogram('h_PiPiMissMass')
    if is_multihad: analysis.addHistogram('h_finalstate_id')
    analysis.loop()

    analysis.addCutNew(
        'TotalP-DeltaE-2',
        lambda: analysis.Variables["DeltaE"].Content > -200. - 1.8 * analysis.Variables["TotalP"].Content
    )
    analysis.addHistogram('h_tptot_tdedx')
    analysis.addHistogram('h_phen')
    analysis.addHistogram('h_phth')
    analysis.addHistogram('h_KpKmPipPimLklhd')
    analysis.addHistogram('h_TotalP_DeltaE')
    analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
    analysis.addHistogram('h_PiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiMissMass2')
    analysis.addHistogram('h_PiPiPiPiMissMass2')
    analysis.addHistogram('h_KPiPiPiMissMass2')
    analysis.addHistogram('h_KKPiPiMissMass2')
    analysis.addHistogram('h_KKMissMass')
    analysis.addHistogram('h_PiPiMissMass')
    if is_multihad: analysis.addHistogram('h_finalstate_id')
    analysis.loop()
    
    analysis.dumpToFile()
    analysis.close()

    if is_multihad:
        final_states = [12, 15, 29, 32] # 12 is for K+K-pi+pi-, 15 for K+KSpi-pi0, 29 for K+K-eta, 32 for K+K-omega
        
        input_dir = "/home/idpershin/analysis/kpkmpippim/pershin_cut_w_TotalP-DeltaE-2"
            
        for final_state in final_states:
            input_path = f"{input_dir}/{prefix}cut{year}_tr_ph_fc_e{energy_point}_{version}.root"
            hists_path = f"{output_dir}/{prefix}fs{final_state}_hists{year}_tr_ph_fc_e{energy_point}_{version}.root"
            output_path = f"{output_dir}/{prefix}fs{final_state}_cut{year}_tr_ph_fc_e{energy_point}_{version}.root"

            analysis = IntermediateAnalysis(
                input_path,
                analysis_path = hists_path,
                output_path = output_path,
            )
            analysis.addCutNew(
                f'finalstate-id_{final_state}',
                lambda: analysis.Variables["finalstate_id"].Content != final_state
            )
            analysis.addHistogram('h_tptot_tdedx')
            analysis.addHistogram('h_phen')
            analysis.addHistogram('h_phth')
            analysis.addHistogram('h_KpKmPipPimLklhd')
            analysis.addHistogram('h_TotalP_DeltaE')
            analysis.addHistogram('h_TotalP_DeltaEKKPiPi')
            analysis.addHistogram('h_PiPiPiMissMass2')
            analysis.addHistogram('h_KPiPiMissMass2')
            analysis.addHistogram('h_PiPiPiPiMissMass2')
            analysis.addHistogram('h_KPiPiPiMissMass2')
            analysis.addHistogram('h_KKPiPiMissMass2')
            analysis.addHistogram('h_KKMissMass')
            analysis.addHistogram('h_PiPiMissMass')
            analysis.loop()
            analysis.dumpToFile()
            analysis.close()


def process_all(is_sim_included = True, is_multihad_included = True):
    path_info = "/home/idpershin/analysis/kpkmpippim/data_info_cmd3.json"
    with open(path_info, 'r') as file_info:
        json_info = json.load(file_info)

    with ProcessPoolExecutor(max_workers = 10) as executor:
        futures = []
        for version in json_info:
            for year in json_info[version]["years"]:
                for elabel in json_info[version]["years"][year]["elabels"]:
                    elabel_data = json_info[version]["years"][year]["elabels"][elabel]
                    energy_point = elabel_data["scan-energy-point"]

                    futures.append(
                        executor.submit(
                            process_single,
                            version,
                            year,
                            energy_point,
                            is_sim = False, is_multihad = False
                        )
                    )

                    if is_sim_included:
                        futures.append(
                            executor.submit(
                                process_single,
                                version,
                                year,
                                energy_point,
                                is_sim = True, is_multihad = False
                            )
                        )
   
                    if is_multihad_included:
                        futures.append(
                            executor.submit(
                                process_single,
                                version,
                                year,
                                energy_point,
                                is_sim = False, is_multihad = True
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

    parser_all = subparsers.add_parser('all', help = "Process all available energy points")
    parser_all.add_argument('--is-sim-included', action = 'store_true')
    parser_all.add_argument('--is-multihad-included', action = 'store_true')
    args = parser.parse_args()

    if args.mode == 'single':
        process_single(args.version, args.year, args.energy, args.is_sim, args.is_multihad)

    if args.mode == 'all':
        process_all(args.is_sim_included, args.is_multihad_included)
