from stochastic_biopolymers.pce_model import PCEConstructor
import os
import pickle
import chaospy as cp
import json
from sklearn import linear_model

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

    regression_model_list = [None, linear_model.Ridge(fit_intercept=False)]
    regression_model = None # linear_model.Ridge(fit_intercept=False)

    output_columns = ['lambda_f', 'force', 'dw', 'ddw']

    # Get folder names in directory data/filament

    data_timestamp = '' if args.data_timestamp is None else args.data_timestamp
    model_timestamp = '' if args.model_timestamp is None else args.model_timestamp

    base_dir = os.getcwd()
    data_timestamp_list = os.listdir(os.path.join(base_dir, 'data', scale))[:10]

    for data_timestamp in data_timestamp_list:
    # for regression_model in regression_model_list:
    # for degree in range(2, 5):
        print(f"Running PCE model with degree {degree}, regression model {regression_model} on data timestamp {data_timestamp}")
        data_path = os.path.join(base_dir, 'data', scale, data_timestamp)
        print(f"Loading data from {data_path}")
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"Data file not found at {data_path}")
        with open(f"{data_path}/results.pkl", 'rb') as file:
            df = pickle.load(file)
        with open(f"{data_path}/study_props_info.csv", 'rb') as file:
            rnd_inputs = json.load(file)

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
        
        # Load dataframe and information about random inputs
        if args.train:
            print("Training PCE model...")
            # Create the polynomial basis
            poly_basis, norms = pce_constructor.poly_basis(degree=degree, 
                                                        truncation=truncation)
            
            # Train the PCE model
            pce = pce_constructor.train_pce(poly_basis=poly_basis,
                                            regression_model=regression_model)

            # Save the PCE model
            model_timestamp = pce_constructor.save_pce(pce=pce,
                                                poly_basis=poly_basis,
                                                data_timestamp=data_timestamp,
                                                regression_model=regression_model)

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
