"""
PyMOL conda environment management module.

This module handles:
1. Detecting existing PyMOL conda environments
2. Creating a new pymol-env if needed
3. Getting the Python interpreter path for PyMOL execution
"""

import os
import subprocess
import shutil
from typing import Optional, Tuple


# Default environment name for PyMOL
PYMOL_ENV_NAME = "pymol-env"


def find_conda_executable() -> Optional[str]:
    """
    Find conda or mamba executable.

    Returns:
        Path to conda/mamba executable or None
    """
    # Prefer mamba over conda
    for exe in ['mamba', 'conda']:
        path = shutil.which(exe)
        if path:
            return path

    # Check common locations
    common_paths = [
        os.path.expanduser("~/miniforge3/bin/mamba"),
        os.path.expanduser("~/miniforge3/bin/conda"),
        os.path.expanduser("~/mambaforge/bin/mamba"),
        os.path.expanduser("~/mambaforge/bin/conda"),
        os.path.expanduser("~/miniconda3/bin/conda"),
        os.path.expanduser("~/anaconda3/bin/conda"),
        "/opt/conda/bin/conda",
    ]

    for path in common_paths:
        if os.path.exists(path):
            return path

    return None


def get_conda_envs(conda_exe: str) -> dict:
    """
    Get list of conda environments.

    Args:
        conda_exe: Path to conda/mamba executable

    Returns:
        Dict mapping env name to path
    """
    try:
        result = subprocess.run(
            [conda_exe, 'env', 'list', '--json'],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            import json
            data = json.loads(result.stdout)
            envs = {}
            for env_path in data.get('envs', []):
                env_name = os.path.basename(env_path)
                envs[env_name] = env_path
            return envs
    except Exception as e:
        print(f"Warning: Could not list conda environments: {e}")

    return {}


def find_pymol_env(conda_exe: Optional[str] = None) -> Optional[str]:
    """
    Find an existing PyMOL conda environment.

    Args:
        conda_exe: Optional path to conda/mamba executable

    Returns:
        Path to PyMOL environment or None
    """
    if conda_exe is None:
        conda_exe = find_conda_executable()

    if conda_exe is None:
        return None

    envs = get_conda_envs(conda_exe)

    # Look for pymol-env or similar
    pymol_env_names = ['pymol-env', 'pymol', 'pymol-open-source']

    for name in pymol_env_names:
        if name in envs:
            return envs[name]

    # Check if any env has pymol installed
    for env_name, env_path in envs.items():
        python_path = os.path.join(env_path, 'bin', 'python')
        if os.path.exists(python_path):
            try:
                result = subprocess.run(
                    [python_path, '-c', 'import pymol; print("ok")'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if result.returncode == 0 and 'ok' in result.stdout:
                    return env_path
            except Exception:
                pass

    return None


def create_pymol_env(conda_exe: Optional[str] = None, env_name: str = PYMOL_ENV_NAME) -> Tuple[bool, str]:
    """
    Create a new conda environment with PyMOL installed.

    Args:
        conda_exe: Optional path to conda/mamba executable
        env_name: Name for the new environment

    Returns:
        Tuple of (success, message or path)
    """
    if conda_exe is None:
        conda_exe = find_conda_executable()

    if conda_exe is None:
        return False, "No conda/mamba executable found"

    print(f"Creating conda environment '{env_name}' with PyMOL...")

    try:
        # Create environment with pymol-open-source from conda-forge
        result = subprocess.run(
            [conda_exe, 'create', '-n', env_name, '-c', 'conda-forge',
             'pymol-open-source', 'python=3.11', '-y'],
            capture_output=True,
            text=True,
            timeout=600  # 10 minutes timeout
        )

        if result.returncode != 0:
            return False, f"Failed to create environment: {result.stderr}"

        # Get the new environment path
        envs = get_conda_envs(conda_exe)
        if env_name in envs:
            return True, envs[env_name]

        return False, "Environment created but not found in list"

    except subprocess.TimeoutExpired:
        return False, "Environment creation timed out"
    except Exception as e:
        return False, f"Error creating environment: {e}"


def get_pymol_python(env_path: Optional[str] = None) -> Optional[str]:
    """
    Get the Python interpreter path for a PyMOL environment.

    Args:
        env_path: Path to conda environment (if None, will search)

    Returns:
        Path to Python interpreter or None
    """
    if env_path is None:
        env_path = find_pymol_env()

    if env_path is None:
        return None

    # Check for Python executable
    python_paths = [
        os.path.join(env_path, 'bin', 'python'),  # Linux/macOS
        os.path.join(env_path, 'python.exe'),      # Windows
        os.path.join(env_path, 'Scripts', 'python.exe'),  # Windows alt
    ]

    for python_path in python_paths:
        if os.path.exists(python_path):
            return python_path

    return None


def setup_pymol_env(auto_create: bool = True) -> Tuple[bool, Optional[str], str]:
    """
    Setup PyMOL environment - find existing or create new.

    Args:
        auto_create: If True, create environment if not found

    Returns:
        Tuple of (success, python_path, message)
    """
    # First check for macOS PyMOL.app
    macos_pymol_python = "/Applications/PyMOL.app/Contents/bin/python"
    if os.path.exists(macos_pymol_python):
        return True, macos_pymol_python, "Using PyMOL.app"

    # Look for conda executable
    conda_exe = find_conda_executable()
    if conda_exe is None:
        return False, None, "No conda/mamba found. Install mamba/conda or use --no-pymol flag."

    # Look for existing PyMOL environment
    env_path = find_pymol_env(conda_exe)

    if env_path:
        python_path = get_pymol_python(env_path)
        if python_path:
            return True, python_path, f"Found existing PyMOL environment at {env_path}"

    # Create new environment if requested
    if auto_create:
        print("No PyMOL environment found. Creating one...")
        success, result = create_pymol_env(conda_exe)

        if success:
            python_path = get_pymol_python(result)
            if python_path:
                return True, python_path, f"Created PyMOL environment at {result}"
            return False, None, f"Environment created at {result} but Python not found"
        else:
            return False, None, result

    return False, None, "No PyMOL environment found and auto_create is False"


def verify_pymol(python_path: str) -> Tuple[bool, str]:
    """
    Verify that PyMOL works with the given Python interpreter.

    Args:
        python_path: Path to Python interpreter

    Returns:
        Tuple of (success, version or error message)
    """
    try:
        result = subprocess.run(
            [python_path, '-c',
             'import pymol; print(f"PyMOL {pymol.cmd.get_version()[0]}")'],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            version = result.stdout.strip()
            return True, version
        else:
            return False, result.stderr

    except Exception as e:
        return False, str(e)


if __name__ == "__main__":
    print("PyMOL Environment Setup")
    print("=" * 50)

    # Check conda
    conda_exe = find_conda_executable()
    print(f"Conda/Mamba: {conda_exe or 'Not found'}")

    if conda_exe:
        envs = get_conda_envs(conda_exe)
        print(f"\nFound {len(envs)} conda environments")

        pymol_env = find_pymol_env(conda_exe)
        print(f"PyMOL environment: {pymol_env or 'Not found'}")

    # Try setup
    print("\n" + "=" * 50)
    print("Running setup_pymol_env()...")
    success, python_path, message = setup_pymol_env(auto_create=False)

    print(f"Success: {success}")
    print(f"Python: {python_path}")
    print(f"Message: {message}")

    if success and python_path:
        print("\nVerifying PyMOL...")
        ok, version = verify_pymol(python_path)
        print(f"Verification: {version if ok else 'Failed - ' + version}")
