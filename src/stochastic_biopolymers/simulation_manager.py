'''
TO DO:
- Add docstrings to all methods
- Find a way to read dw and ddw from simulation_run
'''

import os
import pandas as pd
import chaospy as cp
from datetime import datetime
from stochastic_biopolymers.sampler import Sampler
from stochastic_biopolymers.filament import Filament
# from stochastic_biopolymers.network import Network
# from stochastic_biopolymers.abaqus import AbaqusSimulation

class SimulationManager:
    def __init__(self, base_mat_props: dict, study_props_info: dict,
                 scale: str = None, def_info: dict = None,
                 n_samples: int = 1, joint_dist: cp.J = None, sampling_method: str = 'random',
                 stratify: bool = True, train_ratio: float = 0.8):
        self.scale = scale
        self.n_samples = n_samples
        self.timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        self.base_dir = os.getcwd()
        self.sampling_method = sampling_method
        self.stratify = stratify
        self.train_ratio = train_ratio
        if self.scale == "element":
            self.abaqus_dir = os.path.join(self.base_dir, 'test_in_abaqus')
            self.sim_dir = os.path.join(self.abaqus_dir, 'data', self.timestamp)
            # self.simulation = AbaqusSimulation(job_name='cube_anl', 
            #                                 material_file='material_param.inp', 
            #                                 base_dir=self.base_dir, 
            #                                 abaqus_dir=self.abaqus_dir, 
            #                                 sim_dir=self.sim_dir)
        elif self.scale == "GP":
            self.sim_dir = os.path.join(self.base_dir, 'data', self.timestamp)
            # self.simulation = Network(def_info, self.sim_dir)
        elif self.scale == "filament":
            self.sim_dir = os.path.join(self.base_dir, 'data/filament', self.timestamp)
            print(f"Simulation directory: {self.sim_dir}")
            self.simulation = Filament(def_info, self.sim_dir)
        # print(f"------------\nSimulation directory: {self.sim_dir}\n-------------")
        if not os.path.exists(self.sim_dir):
            os.makedirs(self.sim_dir)
        self.mat_props_manager = Sampler(self.sim_dir, 
                                        base_mat_props, 
                                        study_props_info,
                                        def_info,
                                        joint_dist,
                                        self.sampling_method,
                                        self.n_samples,
                                        )

    def run_study(self):
        """Run the parametric study by looping through material configurations."""
        # Generate samples
        # samples = self.mat_props_manager.generate_samples()
        self.mat_props_manager.write_mat_props()
        if self.scale == "element":
            for c, config in enumerate(self.mat_props_manager.mat_configs):
                print(f"Running simulation for config {c}")
                mat_props = self.mat_props_manager.get_material_props(config)
                os.chdir(self.abaqus_dir)
                # Create input files and run simulation
                self.simulation.create_material_file(mat_props)
                self.simulation.create_rnd_inp_file(mat_props)
                self.simulation.create_section_file()
                self.simulation.delete_old_files()
                self.simulation.run_simulation()
                self.simulation.extract_results()
                self.simulation.move_results(c)
                os.chdir(self.base_dir)
                print(f"Simulation for config {c} completed.")
        elif self.scale == "GP":
            # Create empty dataframe for results
            df_results = pd.DataFrame()
            dict_results = {}
            for iSample in range(self.n_samples):
                print(f"Running simulation for sample {iSample}")
                mat_props = self.mat_props_manager.get_material_props(iSample)
                # os.chdir(self.sim_dir)
                # Run the UMAT
                iResults = self.simulation.run_umat(mat_props)
                # Add results dict to the empty dataframe
                dict_results[iSample] = iResults
                print(f"Simulation for sample {iSample} completed.")
            # Save results to csv
            df_results = pd.DataFrame.from_dict(dict_results, orient='index')
            self.simulation.save_results(df_results, self.mat_props_manager.study_props_info)
        elif self.scale == "filament":
            # Create empty dataframe for results
            df_results = pd.DataFrame()
            dict_results = {}
            for iSample in range(self.n_samples):
                print(f"Running simulation for sample {iSample}")
                mat_props = self.mat_props_manager.get_material_props(iSample)
                # os.chdir(self.sim_dir)
                # Run filament relations
                # iResults = self.simulation.inextensible(mat_props)
                iResults = self.simulation.fil_force(mat_props, dw=True, ddw=True)
                dict_results[iSample] = iResults
                print(f"Simulation for sample {iSample} completed.")
            # Save results to dataframe
            df_results = pd.DataFrame.from_dict(dict_results, orient='index')
            # Add columns for each random input
            df_results = pd.concat([df_results, pd.DataFrame(self.mat_props_manager.samples)], axis=1)
            # Explode dataframe
            df_results = self.mat_props_manager.explode_df(df_results)
            df_results = self.mat_props_manager.split_samples(df_results,
                                                              train_ratio=self.train_ratio,
                                                              stratify=self.stratify)
            self.simulation.save_results(df_results)