import subprocess
import os
import argparse

'''Compiles material routines in main directory and /test_in_abaqus, for deterministic or random approach.'''

def build_umat(approach: str = 'deterministic', compile_umat: bool = False, compile_modules: bool = False):

    if approach == 'deterministic':
        ignore_umat = '_umat_rnd_.f90'
        ignore_affcl = 'affclnetfic_discrete_rnd.f90'
        ignore_uexternal = 'uexternaldb_rnd.f90'
    elif approach == 'random':
        ignore_umat = '_umat_det_.f90'
        ignore_affcl = 'affclnetfic_discrete_det.f90'
        ignore_uexternal = 'uexternaldb_det.f90'

    if compile_umat:
        # 1. Concatenate all .f90 files from *src/* to umat_anl_ai.f90
        # Global module must be at the beginning
        
        subprocess.run(
            "find . -type f -path '*src/umat/*' -name '_global.f90' -exec cat {} + > modules/umat_anl_ai.f90", 
            shell=True, check=True
        )

        subprocess.run(
        f"find . -type f -path '*src/umat/*' -name '*.f90' "
        f"-not \\( -name '_global.f90' -o -name '{ignore_umat}' -o -name '{ignore_affcl}' -o -name '{ignore_uexternal}' \\) "
        f"-exec cat {{}} + >> modules/umat_anl_ai.f90",
        shell=True, check=True
        )

        # 2. Concatenate all .f90 files from *src/* to test_in_abaqus folder
        subprocess.run(
            "find . -type f -path '*src/umat/*' -name '_global.*' -exec cat {} + > test_in_abaqus/umat_anl_ai.f", 
            shell=True, check=True
        )

        subprocess.run(
        f"find . -type f -path '*src/umat/*' -name '*.f90' "
        f"-not \\( -name '_global.f90' -o -name '{ignore_umat}' -o -name '{ignore_affcl}' -o -name '{ignore_uexternal}' \\) "
        f"-exec cat {{}} + >> test_in_abaqus/umat_anl_ai.f",
        shell=True, check=True
        )

        subprocess.run(
            "find . -type f -path '*src/umat/*' -name '_global.*' -exec cat {} + > modules/glb_py.f90", 
            shell=True, check=True
        )

        subprocess.run(
            "find . -type f -name 'global_read_py.f90' -exec cat {} + >> modules/glb_py.f90", 
            shell=True, check=True
        )


        subprocess.run(
            "find . -type f -path '*src/umat/*' -name '_global.*' -exec cat {} + > test_in_abaqus/glb_py.f90", 
            shell=True, check=True
        )

        subprocess.run(
            "find . -type f -name 'global_read_py.f90' -exec cat {} + >> test_in_abaqus/glb_py.f90", 
            shell=True, check=True
        )

        print("UMAT compiled successfully.")

    if compile_modules:
        print("Compiling f2py modules...")
        modules_dir = os.path.abspath('modules')
        os.makedirs(modules_dir, exist_ok=True)
        command1 = ['python', '-m', 'numpy.f2py', 
                    '-c', '-m', 'globpy', 'glb_py.f90']
        command2 = ['python', '-m', 'numpy.f2py', 
                '-c', '-m', 'umatpy', 'py_anl_ai.f90', 'umat_anl_ai.f90']
        try:
            subprocess.run(command1, check=True, cwd=modules_dir) 
            subprocess.run(command2, check=True, cwd=modules_dir)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred: {e}")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build UMAT")
    parser.add_argument("--approach", choices=["deterministic", "random"], default="deterministic", help="Approach to use")
    parser.add_argument("--compile_umat", action="store_true", help="Compile the UMAT")
    parser.add_argument("--compile_modules", action="store_true", help="Compile f2py interface modules")
    args = parser.parse_args()

    build_umat(approach=args.approach, compile_umat=args.compile_umat, compile_modules=args.compile_modules)