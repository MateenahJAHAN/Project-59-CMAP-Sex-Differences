"""Run the complete Phase 6 screen. See fellowship/README.md for inputs and setup."""
from pathlib import Path
import subprocess
import sys

if __name__ == '__main__':
    entry = Path(__file__).resolve().parents[2] / 'fellowship/run.py'
    subprocess.run([sys.executable, str(entry), *sys.argv[1:], '--phases', '6'], check=True)
