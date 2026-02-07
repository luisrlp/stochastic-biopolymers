from stochastic_biopolymers.pce_model import PCEConstructor
import os
import pickle
import chaospy as cp
import json
from sklearn import linear_model
import numpy as np
import pandas as pd

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run PCE model construction and evaluation.')
    parser.add_argument('--train', action='store_true', help='Whether to train the PCE model.')
    parser.add_argument('--test', action='store_true', help='Whether to test the PCE model.')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    parser.add_argument('--scale', type=str, default='filament', help='Scale of the study (element, GP, or filament)')
    parser.add_argument('--degree', type=int, default=2, help='Degree of the polynomial basis for PCE')
    parser.add_argument('--data_timestamp', type=str, default=None, help='Timestamp for the data file')
    parser.add_argument('--model_timestamp', type=str, default=None, help='Timestamp for the PCE model file')
    parser.add_argument('--truncation', type=float, default=1.0, help='Truncation value for polynomial basis')
    parser.add_argument('--output', type=list, default=['force'], help='Output quantity for PCE results')
    args = parser.parse_args()

    seed = args.seed
    scale = args.scale
    degree = args.degree
    truncation = args.truncation
    output_columns = args.output

    # Distributions
    distributions = {
        'normal': cp.Normal,
        'uniform': cp.Uniform,
        'lognormal': cp.LogNormal,
    }

    def_min = None
    def_max = None

    regression_model_list = [#None, 
                             #linear_model.OrthogonalMatchingPursuit(fit_intercept=False), 
                             linear_model.Lars(fit_intercept=False), 
                             #linear_model.LassoLars(fit_intercept=False),]
    ]

    regression_label_list = ['LARS_NON_UNI', 
                             #'LASSO-LARS2'] # ['OLS', 'OMP', 'LARS', 'LASSO-LARS']
    ]
    # regression_model = None # linear_model.OrthogonalMatchingPursuit(fit_intercept=False, n_nonzero_coefs=20) # linear_model.Ridge(fit_intercept=False)
    # regression_model.set_params(n_nonzero_coefs=20)

    output_columns = ['force', 'dw', 'ddw'] # 'lambda_f' ]

    # Get folder names in directory data/filament

    data_timestamp = '' if args.data_timestamp is None else args.data_timestamp
    model_timestamp = '' if args.model_timestamp is None else args.model_timestamp

    base_dir = os.getcwd()
    # data_timestamp_list = os.listdir(os.path.join(base_dir, 'data', scale))[2:3] #, '4k_samples')) #[:1]
    # data_timestamp_list = ['3k_110_115']
    data_timestamp_uniform = '6k_full'
    data_timestamp = '6k_non_uniform' # '6k_full'
    data_file_list = os.listdir(os.path.join(base_dir, 'data', scale, data_timestamp))
    # step 100 from 100 to 3500
    n_train_list = list(range(100, 3001, 100))


    # Data loading from 1 file
    data_path = os.path.join(base_dir, 'data', scale, data_timestamp_uniform)
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data file not found at {data_path}")
    print(f"Loading data from {data_path}")
    # with open(f"{data_path}/results.pkl", 'rb') as file:
    #     df = pickle.load(file)
    with open(f"{data_path}/study_props_info.csv", 'rb') as file:
        rnd_inputs = json.load(file)

    # Data loading from list of files
    df = pd.DataFrame()
    for data_file_id in data_file_list:
        data_path = os.path.join(base_dir, 'data', scale, data_timestamp, data_file_id)
        print(f"Loading data from {data_path}")
    ## Load test results
        with open(f'{data_path}/results.pkl', 'rb') as file:
            df_aux = pickle.load(file)

        # with open(f'{data_path}/study_props_info.csv', 'r') as file:
        #     rnd_inputs = json.load(file)

        df = pd.concat([df, df_aux], ignore_index=True)

    if (def_min and def_max) is not None:
        df = df[(df['deformation'] >= def_min) & (df['deformation'] <= def_max)].reset_index(drop=True)
    
    joint_dist = cp.J(
            *[distributions[value['distribution']](*list(value.values())[1:]) 
            for value in rnd_inputs.values()]
                    )
    
    # Instantiate a class variable using a given dataframe
    pce_constructor = PCEConstructor(df=df, 
                                    rnd_inputs=rnd_inputs, 
                                    joint_dist=joint_dist,
                                    y_columns=output_columns,
                                    seed=seed)
    
    # for data_timestamp in data_timestamp_list:
    # for degree in range(2, 5):
    for iRegression, regression_model in enumerate(regression_model_list):
        regression_label = regression_label_list[iRegression]
        print(f"Using regression model: {regression_label}")
        
        for degree in range(4,6):
            print(f"Using polynomial degree: {degree}")
            if degree > 4 and regression_model is None:
                print(f"Skipping OLS for degree {degree} due to computational cost.")
                continue

            poly_basis, norms = pce_constructor.poly_basis(degree=degree, 
                                        truncation=truncation)
            
            if regression_label == 'OMP' or regression_label == 'LARS_NON_UNI':
                coefs_list = np.linspace(5, 200, 20, dtype=int)
            elif regression_label == 'LASSO-LARS2':
                coefs_list = np.logspace(-5, -2, 10)
            else:
                coefs_list = [0.0]

            for n_train in list(range(250, len(df[df['split'] == 'train']) + 1, 250)):
                print(f"Running PCE model with {n_train} training samples.")

                for coef in coefs_list:

                    # Load dataframe and information about random inputs
                    if args.train:
                        print("Training PCE model...")

                        # Set the n_nonzero_coefs parameter if using OrthogonalMatchingPursuit
                        if regression_label == 'OMP' or regression_label == 'LARS_NON_UNI':
                            print(f"Setting n_nonzero_coefs to {coef}")
                            n_nonzero_coefs = coef
                            regression_model.set_params(n_nonzero_coefs=n_nonzero_coefs)
                            alpha = None
                        elif regression_label == 'LASSO-LARS2':
                            print(f"Setting alpha to {coef}")
                            n_nonzero_coefs = None
                            alpha = coef
                            regression_model.set_params(alpha=alpha)
                        else:
                            n_nonzero_coefs = None
                            alpha = None

                        # Train the PCE model
                        pce = pce_constructor.train_pce(poly_basis=poly_basis,
                                                        regression_model=regression_model,
                                                        n_train=n_train) # "None" to use all training data

                        # Save the PCE model
                        model_timestamp = pce_constructor.save_pce(pce=pce,
                                                            poly_basis=poly_basis,
                                                            data_timestamp=data_timestamp,
                                                            regression_model=regression_model,
                                                            dir_info=regression_label,
                                                            n_train=n_train,
                                                            n_nonzero_coeffs=n_nonzero_coefs,
                                                            alpha=alpha)

                    if args.test:
                        print("Testing PCE model...")
                        # Load a PCE model
                        pce = pce_constructor.load_pce(timestamp=model_timestamp)
                        # Load a polynomial basis
                        poly_basis = pce_constructor.load_poly_basis(timestamp=model_timestamp)
                        if pce is None or poly_basis is None:
                            raise ValueError("PCE model or polynomial basis could not be loaded.")
                        
                        # Test
                        test_dict = pce_constructor.test_pce(pce_dict=pce,
                                                            deformation=None,
                                                            split='test')
                        
                        # Save the test results
                        pce_constructor.save_test_results(test_results=test_dict,
                                                        timestamp=model_timestamp)
