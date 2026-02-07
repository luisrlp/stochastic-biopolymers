import os
import json
import pandas as pd
import numpy as np
import chaospy as cp
from sklearn.model_selection import train_test_split

class Sampler:
    """
    - Manage material properties and their different configurations.
    - Generate material-related input files
    """

    def __init__(self, data_dir, base_mat_props: dict, 
                 study_props_info: dict = None,
                 def_info: dict = None, 
                 joint_dist: cp.J = None,
                 sampling_method: str = 'random',
                 n_samples: int = 1):
        
        self.data_dir = data_dir
        self.base_mat_props = base_mat_props
        self.study_props_info = study_props_info
        self.def_info = def_info
        self.joint_dist = joint_dist
        self.sampling_method = sampling_method
        self.n_samples = n_samples
        self.samples = self._generate_samples()

    def _generate_samples(self) -> dict:
        """Generate samples of material properties based on the joint distribution and sampling method."""
        valid_methods = ['random', 'additive_recursion', 'sobol', 'halton', 'hammersley', 'korobov', 'latin_hypercube']
        if self.sampling_method not in valid_methods:
            raise ValueError(f"Invalid sampling method. Choose from: {valid_methods}")
    
        if self.joint_dist is None:
            raise ValueError("joint_dist must be defined to generate samples.")
        
        if self.study_props_info is None:
            raise ValueError("study_props_info must be provided to map sample dimensions.")
        
        if self.joint_dist.sample(1).shape[0] != len(self.study_props_info):
            raise ValueError("Mismatch between joint_dist dimensions and study_props_info keys.")
        
        samples_list = self.joint_dist.sample(self.n_samples, rule=self.sampling_method)
        # If there is only one random property, reshape (100,) to (1,100)
        if samples_list.ndim == 1:
            samples_list = samples_list.reshape(1, -1)
        samples_dict = {key: samples_list[i, :] for i, key in enumerate(self.study_props_info.keys())}
        return samples_dict

    def get_material_props(self, iSample: int) -> dict:
        """Material properties of a given configuration"""
        mat_props = self.base_mat_props.copy()
        for key, value in self.samples.items():
            mat_props[key] = value[iSample]
        return mat_props
    
    def write_mat_props(self):
        file_names = ['base_mat_props.csv', 
                      'study_props_info.csv', 
                      'def_info.csv']
        dicts = [self.base_mat_props, self.study_props_info, self.def_info]
        for i_file, file in enumerate(file_names):
            with open(os.path.join(self.data_dir, file), 'w') as json_file:
                json.dump(dicts[i_file], json_file, indent=4)

    def explode_df(self, df: pd.DataFrame,
                   columns: list = None) -> pd.DataFrame:
        """
        Explode the DataFrame by expanding the 'samples' column into multiple rows.
        
        Parameters:
        - df: DataFrame containing the results.
        - columns: List of columns to keep in the exploded DataFrame.
        
        Returns:
        - Exploded DataFrame with each sample in a separate row.
        """
        if columns is None:
            explode_columns = [col for col in df.columns 
                               if isinstance(df[col].iloc[0], (list, np.ndarray))]
        else:
            explode_columns = df.columns.tolist()
        
        exploded_df = df.explode(explode_columns)
        exploded_df.reset_index(drop=True, inplace=True)
        return exploded_df

    # def split_samples(self, 
    #                   df: pd.DataFrame, 
    #                   train_ratio: float = 0.8,
    #                   output: str = 'force',
    #                   stratify: bool = False,
    #                   n_bins: int = 10,
    #                   seed: int = 42) -> dict:
    #     """
    #     Split the samples into training and testing sets based on the provided DataFrame.
        
    #     Parameters:
    #     - df: DataFrame containing the results.
    #     - train_ratio: Proportion of samples to include in the training set.
    #     - output: Output column.
    #     - stratify: Whether to stratify the split based on the output column.
    #     - n_bins: Number of bins for quantile-based stratification.
        
    #     Returns:
    #     - Concatenated DataFrame with 'train' and 'test' splits, with 'split' column.
    #     """
    #     print(f"output column: {output}")
    #     if output not in df.columns:
    #         raise ValueError(f"Output column '{output}' not found in DataFrame.")
    #     df[f'{output}_binned'] = pd.qcut(df[output], q=n_bins, duplicates='drop')

    #     if stratify:
    #         train_df, test_df = train_test_split(
    #             df, 
    #             train_size=train_ratio, 
    #             stratify=df[f'{output}_binned'], 
    #             random_state=seed
    #         )
    #     else:
    #         train_df, test_df = train_test_split(
    #             df, 
    #             train_size=train_ratio, 
    #             random_state=seed
    #         )
        
    #     # Concatenate the train and test DataFrames
    #     train_df['split'] = 'train'
    #     test_df['split'] = 'test'
    #     df_split = pd.concat([train_df, test_df], ignore_index=True)
    #     df_split.drop(columns=[f'{output}_binned'], inplace=True)
    #     return df_split
    
    def split_samples(self, 
                    df: pd.DataFrame, 
                    train_ratio: float = 0.8,
                    output: str = 'force',
                    stratify: bool = False,
                    n_bins: int = 10,
                    seed: int = 42) -> pd.DataFrame:
        """
        Split the samples into training and testing sets for each unique deformation value.
        Allows for stratified sampling based on quantiles of 1 output variable.

        Parameters:
        - df: DataFrame containing the results.
        - train_ratio: Proportion of samples to include in the training set.
        - output: Output column.
        - stratify: Whether to stratify the split based on the output column.
        - n_bins: Number of bins for quantile-based stratification.
        - seed: Random seed for reproducibility.

        Returns:
        - Concatenated DataFrame with 'train' and 'test' splits, with 'split' column.
        """
        if output not in df.columns:
            raise ValueError(f"Output column '{output}' not found in DataFrame.")

        split_dfs = []
        # ADD: IGNORE FIRST IF , IN CASE STRETCH IS A COLUMN IN THE DATAFRAME
        if 'STRETCH' not in df.columns:
            # When stretch is not treated as a random input
            for deformation_value in df['deformation'].unique():
                df_sub = df[df['deformation'] == deformation_value].copy()
                if len(df_sub) < 2:
                    # Not enough samples to split, assign all to train
                    df_sub['split'] = 'train'
                    split_dfs.append(df_sub)
                    continue

                df_sub[f'{output}_binned'] = pd.qcut(df_sub[output], q=min(n_bins, len(df_sub)), duplicates='drop')

                if stratify and len(df_sub[f'{output}_binned'].unique()) > 1:
                    train_df, test_df = train_test_split(
                        df_sub, 
                        train_size=train_ratio, 
                        stratify=df_sub[f'{output}_binned'], 
                        random_state=seed
                    )
                else:
                    train_df, test_df = train_test_split(
                        df_sub, 
                        train_size=train_ratio, 
                        random_state=seed
                    )

                train_df['split'] = 'train'
                test_df['split'] = 'test'
                split_df = pd.concat([train_df, test_df], ignore_index=True)
                split_df.drop(columns=[f'{output}_binned'], inplace=True)
                split_dfs.append(split_df)

            df_split = pd.concat(split_dfs, ignore_index=True)
        else:
            # When stretch is treated as a random input
            df_sub = df.copy()

            df_sub[f'{output}_binned'] = pd.qcut(df_sub[output], q=min(n_bins, len(df_sub)), duplicates='drop')

            if stratify and len(df_sub[f'{output}_binned'].unique()) > 1:
                train_df, test_df = train_test_split(
                    df_sub, 
                    train_size=train_ratio, 
                    stratify=df_sub[f'{output}_binned'], 
                    random_state=seed
                )
            else:
                train_df, test_df = train_test_split(
                    df_sub, 
                    train_size=train_ratio, 
                    random_state=seed
                )

            train_df['split'] = 'train'
            test_df['split'] = 'test'
            df_split = pd.concat([train_df, test_df], ignore_index=True)
            df_split.drop(columns=[f'{output}_binned'], inplace=True)

        return df_split

