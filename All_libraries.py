import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# List of libraries to install
libraries = [
    "math",  # Standard library, no need to install
    "csv",  # Standard library, no need to install
    "typing",  # Standard library, no need to install
    "rdkit",  # Chemistry toolkit
]

# Installing the non-standard libraries
for lib in libraries:
    if lib not in ["math", "csv", "typing"]:
        install(lib)

print("Installation complete.")
