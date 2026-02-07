import pickle
import json
import numpy as np
import pandas as pd
import chaospy as cp
import numpoly
from datetime import datetime
import os


class PCEConstructor:
    def __init__(self, df: pd.DataFrame, 
                 rnd_inputs: dict, 
                 joint_dist: cp.J = None,
                 # y_column: str = 'force',
                 y_columns: list = ['force'],
                 seed: int = 42):
        self.df = df
        self.rnd_inputs = rnd_inputs
        self.joint_dist = joint_dist
        self.seed = seed
        # self.y_column = y_column
        self.y_columns = y_columns

    def x(self, split: str = 'train', deformation: float = None) -> np.ndarray:
        if split not in self.df['split'].unique():
            raise ValueError(f"Split '{split}' not found in DataFrame. Available splits: {self.df['split'].unique()}")
        x_columns = list(self.rnd_inputs.keys())
        split_df = self.df[self.df['split'] == split]
        if deformation is not None:
            split_df = split_df[split_df['deformation'] == deformation]
        return np.array(split_df[x_columns]).T

    def y(self, split: str = 'train', deformation: float = None) -> np.ndarray:
        if split not in self.df['split'].unique():
            raise ValueError(f"Split '{split}' not found in DataFrame. Available splits: {self.df['split'].unique()}")
        split_df = self.df[self.df['split'] == split]
        if deformation is not None:
            split_df = split_df[split_df['deformation'] == deformation]
        # Normalize output -> ln(1 + y)
        # return np.log1p(np.array(split_df[self.y_columns].values.tolist()))
        return np.array(split_df[self.y_columns].values.tolist())
    
    def standardize_x(self, x: np.ndarray) -> np.ndarray:
        ''' Standardizes the input data according to the distribution of each variable (in rnd_inputs).'''
        standardized = np.empty_like(x)
        for i, col in enumerate(self.rnd_inputs.keys()):
            dist_type = self.rnd_inputs[col]['distribution']
            if dist_type == 'normal':
                mean = self.rnd_inputs[col]['mean']
                std = self.rnd_inputs[col]['std']
                standardized[i] = (x[i] - mean) / std
            elif dist_type == 'uniform':
                lower = self.rnd_inputs[col]['low']
                upper = self.rnd_inputs[col]['high']
                # Standardize to [-1, 1]
                standardized[i] = 2 * (x[i] - lower) / (upper - lower) - 1
            else:
                raise ValueError(f"Unknown distribution type '{dist_type}' for variable '{col}'")
        return standardized

    def poly_basis(self, degree: int,
                   rule: str = 'three_terms_recurrence',
                   normed: bool = False,
                   truncation: float = 1.0):
        
        '''Generates polynomial basis for PCE expansion.'''
        poly_basis, norms = cp.generate_expansion(order=degree,
                                                  dist=self.joint_dist,
                                                  rule=rule,
                                                  normed=True, # normed,
                                                  retall=True,
                                                  cross_truncation=truncation)
        return poly_basis, norms

    def train_pce(self, 
                  poly_basis: numpoly.baseclass.ndpoly,
                  regression_model = None,
                  n_train: int = None) -> dict:
        np.random.seed(self.seed)
        pce_dict = {}
        if n_train is None:
            n_train = len(self.df[self.df['split'] == 'train'])
        # Get unique values of the deformation column
        if 'STRETCH' not in self.df.columns:
            deformation_values = self.df['deformation'].unique()
            pce_dict = {
                deformation: cp.fit_regression(
                poly_basis,
                self.x(split='train', deformation=deformation)[:, :n_train],
                # self.standardize_x(self.x(split='train', deformation=deformation)[:, :n_train]),
                self.y(split='train', deformation=deformation)[:n_train, :],
                model=regression_model
                )
                for deformation in deformation_values
            }
        else:
            # Dummy value for pipeline compatibility
            dummy_deformation = 1.0
            pce_dict[dummy_deformation] = cp.fit_regression(
                poly_basis,
                # self.standardize_x(self.x(split='train')[:, :n_train]),
                self.x(split='train')[:, :n_train],
                self.y(split='train')[:n_train, :],
                model=regression_model,
            )
        return pce_dict
    
    def test_pce(self, pce_dict: dict, 
                 deformation: float = None, 
                 split: str = 'test') -> dict:
        '''Get the mean error, mean squared error, and r2 score for the PCE predictions.'''
        results = {}
        n_test = len(self.df[self.df['split'] == split])
        # print(f"Number of test samples: {n_test}")
        if deformation is not None or len(pce_dict.keys()) == 1:
            # PCE with stretch as input
            if deformation is None:
                deformation = list(pce_dict.keys())[0]
                # x = self.standardize_x(self.x(split=split))
                x = self.x(split=split)
                y = self.y(split=split)
            else:            
                # x = self.standardize_x(self.x(split=split, deformation=deformation))
                x = self.x(split=split, deformation=deformation)
                y = self.y(split=split, deformation=deformation)
            pce = pce_dict[deformation]
            pce_output = pce(*x)
            pce_output = np.array(pce_output).T   #.reshape(y.shape)
            # Unnormalize output before comparison
            # pce_output = np.expm1(pce_output)
            # y = np.expm1(y)
            MAE = np.mean(abs(pce_output - y), axis=0)
            MAPE = np.mean(abs((pce_output - y) / y), axis=0)
            top20_MAPE = [sorted(abs((pce_output - y) / y)[:,i], reverse=True)[:20] for i in range(y.shape[1])]
            MSE = np.mean(((pce_output - y)) ** 2, axis=0)
            RMSE = np.sqrt(MSE)
            R2 = 1 - (np.sum((pce_output - y) ** 2, axis=0) / np.sum((y - np.mean(y, axis=0)) ** 2, axis=0))
            results[deformation] = {
            'real_output': y,
            'predictions': pce_output,
            'MAE': MAE,
            'MAPE': MAPE,
            'RMSE': RMSE,
            'R2': R2,
            'HIGH_MAPE': top20_MAPE,
            }
        else:
            pce = pce_dict
            for deformation_value, pce_model in pce.items():
                x_val = self.x(split=split, deformation=deformation_value)
                y_val = self.y(split=split, deformation=deformation_value)
                pce_output = pce_model(*x_val)
                mean_error = np.mean((pce_output - y_val) / y_val)
                r2_score = 1 - (np.sum((pce_output - y_val) ** 2) / np.sum((y_val - np.mean(y_val)) ** 2))
                results[deformation_value] = {
                    'real_output': y_val,
                    'predictions': pce_output,
                    'mean_error': mean_error,
                    'r2_score': r2_score,
                }
        return results

    def sobol_indices(self, pce_dict: dict,
                      deformation: float = None) -> dict:
        """        Computes Sobol indices for the PCE model.
        If deformation is None, it computes Sobol indices for all deformations in pce_dict.
        """
        sobol_indices = {}
        if deformation is not None or len(pce_dict.keys()) == 1:
            if deformation is None:
                deformation = list(pce_dict.keys())[0]
            pce = pce_dict[deformation]
            sobol_indices[deformation] = cp.Sens_m(pce, self.joint_dist)
        else:
            for deformation_value, pce_model in pce_dict.items():
                sobol_indices[deformation_value] = cp.Sens_m(pce_model, self.joint_dist)
        return sobol_indices
    
    def save_pce(self, pce: dict, 
                 poly_basis: numpoly.baseclass.ndpoly,
                 data_timestamp: str = None,
                 regression_model = None,
                 dir_info: str = None,
                 n_train: int = None,
                 n_nonzero_coeffs: int = None,
                 alpha: float = None
                 ) -> str:
        """
        Saves the PCE model and the polynomials basis as pickle files in the 'models' folder with a timestamped subfolder.
        """
        now = datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S_%f")
        if dir_info is not None:
            timestamp = f"{dir_info}/{timestamp}"
        save_dir = f'models/{timestamp}'
        os.makedirs(save_dir, exist_ok=True)
        with open(f'{save_dir}/pce_model.pkl', 'wb') as file:
            print(f"Saving PCE model to {save_dir}/pce_model.pkl")
            pickle.dump(pce, file)
        with open(f'{save_dir}/poly_basis.pkl', 'wb') as file:
            pickle.dump(poly_basis, file)
        # Save training information
        degree = max(sum(exp) for exp in poly_basis.exponents)
        training_info = {
            'data_timestamp': data_timestamp,
            'n_train_samples': n_train,
            'degree': int(degree),
            'regression_model': str(regression_model),
            'rnd_inputs': self.rnd_inputs,
            'y_columns': self.y_columns,
            'seed': self.seed,
            'n_nonzero_coeffs': int(n_nonzero_coeffs) if n_nonzero_coeffs is not None else None,
            'alpha': float(alpha) if alpha is not None else None
        }
        with open(f'{save_dir}/training_data.csv', 'w') as file:
            json.dump(training_info, file)
        return timestamp
    
    def save_test_results(self, test_results: dict,
                          timestamp: str):
        """
        Saves the test results as a pickle file.
        """
        save_dir = f'models/{timestamp}'
        os.makedirs(save_dir, exist_ok=True)
        with open(f'{save_dir}/test_results.pkl', 'wb') as file:
            pickle.dump(test_results, file)

    def load_pce(self, timestamp: str):
        """
        Loads the PCE model and the polynomials basis from the 'models' folder using a timestamped subfolder.
        """
        save_dir = f'models/{timestamp}'
        with open(f'{save_dir}/pce_model.pkl', 'rb') as file:
            pce = pickle.load(file)
        return pce
    
    def load_poly_basis(self, timestamp: str):
        """
        Loads the polynomial basis from the 'models' folder using a timestamped subfolder.
        """
        save_dir = f'models/{timestamp}'
        with open(f'{save_dir}/poly_basis.pkl', 'rb') as file:
            poly_basis = pickle.load(file)
        return poly_basis

    # def save_pce(self, pce: dict):
    #     '''Saves pce_dict {int: cp.fit_regression(...)}, joint_dist and info about training data
    #     in folder "models"
        
    #     Can be loaded with numpy.load()
    #     '''
    #     now = datetime.now()
    #     timestamp = now.strftime("%Y%m%d_%H%M%S")
    #     filename = f"pce_{timestamp}.npz"
    #     filepath = f"models/{filename}"