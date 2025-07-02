import os
import numpy as np
import pandas as pd
import umatpy 

class Network:
    def __init__(self, def_info : dict, sim_dir : str):
        self.sim_dir = sim_dir
        self.def_mode = def_info['def_mode']
        self.ramp_time = def_info['ramp_time']
        self.def_initial = def_info['def_initial']
        self.def_max = def_info['def_max']
        self.increments = def_info['increments']
        self.time_initial = def_info['time_initial']
        self.time_final = def_info['time_final']

    def deformation_gradient(self, deformation : float):
        """Computes the deformation gradient F for simple deformation modes
        INPUTS:
        - def_mode: deformation mode (str) - 'U', 'SSx', 'SSy' 
        - deformation: deformation [stretch, shear] value (float)
        OUTPUTS:
        - F: deformation gradient (np.array, 3x3)
        """
        if self.def_mode not in ['U', 'SSx', 'SSy']:
            raise ValueError("Invalid deformation mode. Choose 'U', 'SSx', or 'SSy'.")
        if not isinstance(deformation, (int, float)):
            raise TypeError("Stretch/shear value must be a number.")
        F = np.eye(3)
        if self.def_mode == 'U':
            F[0, 0] = deformation
            F[1, 1] = 1/np.sqrt(deformation)
            F[2, 2] = 1/np.sqrt(deformation)
        elif self.def_mode == 'SSx':
            F[0, 1] = deformation
        elif self.def_mode == 'SSy':
            F[1, 0] = deformation
        return F
    
    def ramp(self, time : float, ramp_time : float, def_initial : float, def_max : float):
        """Ramp function for deformation"""
        if time < ramp_time:
            return def_initial + (def_max - def_initial) * time / ramp_time
        else:
            return def_max

    def run_umat(self, mat_props: dict):
        """Run UMAT with given deformation and ramp time."""
        time_array = np.linspace(self.time_initial, self.time_final, self.increments)
        dtime = time_array[1] - time_array[0]
        # s_xx = np.zeros((self.increments))
        # s_yy = np.zeros((self.increments))
        # s_zz = np.zeros((self.increments))
        # s_xy = np.zeros((self.increments))
        # s_shear = np.zeros((self.increments))
        # sef_array = np.zeros((self.increments))
        # def_array = np.zeros((self.increments))
        stress_array = np.zeros((self.increments, 6))
        sef_array = np.zeros(self.increments)
        def_array = np.zeros(self.increments)
        mat_props_list = [mat_props[key] for key in mat_props]
        for i, t in enumerate(time_array):
            deformation = self.ramp(t, self.ramp_time, self.def_initial, self.def_max)
            def_array[i] = deformation
            F = self.deformation_gradient(deformation)
            # Convert mat_props to list
            stress, sef = umatpy.run_umat_py(mat_props_list, 
                                             F, 
                                             [t, t + dtime], dtime, 
                                             i+1)
            print(f'F: {F}')
            print(f"Stress: {stress}")
            print(f"SEF: {sef}")
            sef_array[i] = sef
            if self.def_mode == 'SSx':
                stress[1] = stress[1] - stress[2]
                stress[2] = 0 # ????
            stress_array[i,:] = stress
            # if t == self.time_final:
            #     print(f"Stretch/shear: {deformation}")
            #     print(f"Final stress: {stress}")
            #     print(f"Final SEF: {sef}")
        return {'stress': np.array(stress_array),
                'sef': np.array(sef_array),
                'deformation': np.array(def_array),
                'time': np.array(time_array),
                }
    
    def save_results(self, results: pd.DataFrame, study_props: dict):
        """Save results to csv file"""
        # Add a column for each study property        
        for key, _ in study_props.items():
            results[key] = study_props[key]
        # Save results to a dataframe
        results = pd.DataFrame(results)
        # Save dataframe to pkl
        results.to_pickle(os.path.join(self.sim_dir, 'results.pkl'))