import argparse
import chaospy as cp
import numpy as np

from stochastic_biopolymers.simulation_manager import SimulationManager



if __name__ == "__main__":

    # Script arguments
    parser = argparse.ArgumentParser(description='Run a parametric study of a material model in Abaqus.')
    parser.add_argument('--compile', type=bool, default=False, help='Whether to compile the UMAT.')
    parser.add_argument('--scale', type=str, default='filament', help='The scale of the study (element, GP or filament).')
    parser.add_argument('--n_samples', type=int, default=100, help='The number of samples to generate.')
    parser.add_argument('--seed', type=int, default=42, help='The seed value for random number generation.')
    args = parser.parse_args()

    N_samples = args.n_samples
    seed = args.seed
    compile = args.compile
    scale  = args.scale
    N_samples_list = list(range(100, 2501, 100))

    # Base Material Properties
    base_mat_props = {
        # K: Penalty parameter
        'K': 1000.0,
        # C10, C01: Isotropic matrix constants
        'C10': 1.0, 'C01': 1.0,
        # PHI: Filament volume fraction 
        'PHI': 1.0,
        # L: Contour length | CACTIN: Actin concentration
        'L': 1.96, # CACTIN: 9.5
        # ROF: Filament end-to-end distance | R: Ratio between actin and crosslinker concentration
        'R0F': 1.63, # 'R': 0.1'
        # ROC: Crosslinker end-to-end distance
        'R0C': 0.014,
        # ETA: Crosslinker relative stiffness
        'ETA': 0.5, 
        # MU0: Stretch modulus [pN] (inextensible filaments: MU0 -> inf.)
        'MU0': 38600.0, 
        # BETA: Exponent parameter
        'BETA': 0.5,
        # LP: Filament persistence length| B0: Filmanent bending stiffness (T*Lp*kb) [pN * microm**2]
        'Lp': 16.0, # 'B0': 294.0 * 16.0 * 1.38065e-5, 
        # LAMBDA0: Network pre-stretch
        'LAMBDA0': 1.0,
        # NA: Isotropic filaments per unit volume [microm**-3] | A: Ratio between filament contour length and end-to-end distance L/ROF
        'NA': 7.66, # 'A': 1.2,
        # BDISP: Filament dispersion parameter
        'BDISP': 0.001
    }

    np.random.seed(seed)

    distributions = {
        'normal': cp.Normal,
        'uniform': cp.Uniform, 
        'lognormal': cp.LogNormal,
    }

    methods_list = ['random', 'sobol', 'korobov','latin_hypercube'] # 'additive_recursion', 'sobol', 'halton', 'hammersley', 'korobov', 'latin_hypercube']
    sampling_method = 'random'  # 'random', 'additive_recursion', 'sobol', 'halton', 'hammersley', 'korobov'    

    ''' Scale ["filament", "GP", "element"]''' 
    # filament: study filament relations (filament level - 1D)
    # GP: study umat response (Network / Gauss Point level)
    # element: study response of a unit cube in abaqus (element level)

    ''' Approach ["deterministic", "random"] '''
    # "deterministic": same properties for all filament directions in a GP
    # "random": properties of each direction determined by a PDF
    approach = "deterministic"
    
    ''' Random Parameters '''
    study_props_info = {
        # 'C10': {'distribution': 'normal',
        #         'mean': 1.0,
        #         'std': 0.2},
        # 'C01': {'distribution': 'uniform',
        #         'low': 0.7,
        #         'high': 1.3},
        # 'CACTIN': {'distribution': 'normal',
        #         'mean': 9.5,
        #         'std': 1.5},
        # 'R': {'distribution': 'lognormal',
        #         'mean': 0.1,
        #         'std': 0.01},
        'R0C': {'distribution': 'normal',
                'mean': 0.014,
                'std': 0.0014},
        'ETA': {'distribution': 'normal',
                'mean': 0.5,
                'std': 0.05},
        # 'BDISP': {'distribution': 'uniform',
        #         'low': 0.8,
        #         'high': 1.2},
        'MU0': {'distribution': 'normal',
                'mean': 38.6e3,
                'std': 3.86e3},
        # 'L': {'distribution': 'uniform',
        #         'low': 1.95,
        #         'high': 2.05},
        'L': {'distribution': 'uniform',
                'low': 1.9,
                'high': 2.1},
        'R0F': {'distribution': 'uniform',
                'low': 1.55,
                'high': 1.7},
        # 'Lp': {'distribution': 'normal',
        #         'mean': 16.0,
        #         'std': 1.6},
        'Lp': {'distribution': 'normal',
               'mean': 16.0,
               'std': 1.6},
        # 'BETA': {'distribution': 'normal',
        #         'mean': 0.5,
        #         'std': 0.1},
        'STRETCH': {'distribution': 'uniform',  
                    'low': 1.00, # not working from stretch higher than 1.175
                    'high': 1.05}, 
    }

    # Create a joint distribution for the random parameters
    joint_dist = cp.J(
                *[distributions[value['distribution']](*list(value.values())[1:]) 
                for value in study_props_info.values()]
                     )

    # Dictionary for random material properties: {'eta': [eta_1, ..., eta_nsamples]}
    study_props = {key: distributions[value['distribution']](*list(value.values())[1:])
                for key, value in study_props_info.items()}

    '''Deformation info'''
    def_info = {
        'def_mode': 'U', # 'U', 'SSx', 'SSy'
        'def_initial': 1.0,
        'def_max': 1.2,
        'increments': 1,
        'time_initial': 0.0,
        'time_final': 1.0,
        'ramp_time': 1.0,
        'sampling_method': sampling_method,
    }

    # Add dict with all info for train/test split (stratify, train_ratio) then pass it to simulation manager
    stratify = True  # Stratified sampling

    # if approach == "deterministic":
    #     study_props = study_props_deterministic
    # elif approach == "random":
    #     study_props = study_props_random
    # else:
    #     raise ValueError("Approach not valid. Please choose 'deterministic' or 'random'.")

    # Initialize and run the simulation manager
    #     for sampling_method in methods_list:
#     print(f"Using sampling method: {sampling_method}")
#     for N_samples in N_samples_list:
    print(f"Running simulation with N_samples: {N_samples}")
    sim_manager = SimulationManager(base_mat_props,
                                study_props_info,
                                scale, 
                                def_info,
                                N_samples,
                                joint_dist,
                                sampling_method,
                                stratify)
    sim_manager.run_study()