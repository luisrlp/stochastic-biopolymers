import os
import glob
import subprocess
import numpy as np
import numpy as np
import shutil
import globpy

class AbaqusSimulation:
    """
    - Manage Abaqus simulations, from input creation to output extraction
    """

    def __init__(self, job_name='cube_anl', material_file='material_param.inp', 
                 base_dir=None, abaqus_dir=None, sim_dir=None):
        self.job_name = job_name
        self.material_file = material_file
        self.base_dir = base_dir
        self.abaqus_dir = abaqus_dir
        self.sim_dir = sim_dir
        [self.nsdv, self.nelem, self.ndir, self.ngp] = globpy.global_read_py()
        self.abq_extensions = ['.com', '.dat', '.msg', '.prt', '.sim', '.sta', '.odb', '.par']
    
    def delete_old_files(self):
        """Delete temp files and old odbs."""    
        for ext in self.abq_extensions:
            files = glob.glob(f'{self.job_name}*{ext}')
            for file in files:
                try:
                    os.remove(file)
                except OSError as e:
                    print(f"Error deleting file {file}: {e}")
    
    def create_material_file(self, mat_props):
        """Create the material parameter inp file from given material properties."""
        with open(self.material_file, 'w') as file:
            file.write("*parameter\n")
            for key, value in mat_props.items():
                file.write(f"{key} = {value}\n")

    def create_rnd_inp_file(self, mat_props):
        eta = mat_props['ETA']
        with open('etadir.inp', 'w') as inp_file:
            for el in range(1, self.nelem + 1):
                for gp in range(1, self.ngp + 1):
                    rnd_array = np.random.normal(eta, eta*0.1, self.ndir)
                    rnd_array = np.clip(rnd_array * (eta*self.ndir) / np.sum(rnd_array), 0, 1)
                    line = f"{el}, {gp}, " + ', '.join(f"{prop}" for prop in rnd_array) + '\n'
                    inp_file.write(line)

    def create_section_file(self):
        # Create section input file
        filename = 'sec_anl.inp'
        with open(filename, 'w') as file:
            file.write('*Solid Section, elset=main_element, material=UD\n')
            # Write the Material Definition
            file.write('*Material, name=UD\n')
            file.write('*User Material, constants=14\n')
            file.write('<K>,<C10>,<C01>,<PHI>,<L>,<R0F>,<R0C>,<ETA>\n')
            file.write('<MU0>,<BETA>,<B0>,<LAMBDA0>,<NA>,<BDISP>\n')
            file.write('*DEPVAR\n')
            file.write(f'{self.nsdv},\n')
            file.write('1, DET, "DET"\n')
            if self.nsdv > 1:
                for i in range(2, self.nsdv+1):
                    file.write(f'{i}, ETAC{i-1}, "ETAC{i-1}"\n')
    
    def run_simulation(self):
        """Run Abaqus simulation."""
        try:
            # Add gpu support
            subprocess.run(['abaqus', 'job=' + self.job_name + '.inp', 'user=umat_anl_ai.f', '-interactive'], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running simulation: {e}")

    def extract_results(self):
        """Extract results with getoutput.py"""
        try:
            subprocess.run(['abaqus', 'cae', '-noGUI', 'getoutput.py'], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error extracting results: {e}")

    def move_results(self, index):
        """ Move abaqus (including random input) files and numpy outputs to sim directory """
        dest_path = os.path.join(self.sim_dir, str(index))
        files = glob.glob(f'*npy')
        os.makedirs(dest_path)
        for ext in self.abq_extensions:
            ext_files = glob.glob(f'{self.job_name}*{ext}')
            files.extend(ext_files)
        files.append('etadir.inp')
        for file in files:
            filepath = os.path.join(self.abaqus_dir, file)
            shutil.move(filepath, dest_path)