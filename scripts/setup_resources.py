import os
import urllib.request
import sys
from pathlib import Path

# Define paths relative to this script
# Script is in scripts/, so parent is root, then ppr_analysis/data
ROOT_DIR = Path(__file__).parent.parent
DATA_DIR = ROOT_DIR / "ppr_analysis" / "data"
BIN_DIR = DATA_DIR / "bin"
PROFILES_DIR = DATA_DIR / "profiles"

def download_file(url, dest):
    if dest.exists():
        print(f"  [Skip] Already exists: {dest.name}")
        return
    
    print(f"  [Download] {url} -> {dest.name}...")
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:
        print(f"  [Error] Failed to download {url}: {e}")
        sys.exit(1)

def main():
    print("Setting up resources for PPR Motif Analyzer...")
    
    # Create directories
    os.makedirs(BIN_DIR, exist_ok=True)
    os.makedirs(PROFILES_DIR, exist_ok=True)
    
    # URLs
    PPRFINDER_JAR_URL = 'https://raw.githubusercontent.com/YaoYinYing/PPRfinder/master/out/production/PPRfinder/PPRfinder.jar'
    HMM_PROFILE_BASE_URL = 'https://raw.githubusercontent.com/YaoYinYing/PPRfinder/master/data/'
    PROSITE_URL = 'https://ftp.expasy.org/databases/prosite/'

    print(f"Target Directory: {DATA_DIR}")

    # Download PPRFinder.jar (No longer needed, using Python implementation)
    # download_file(PPRFINDER_JAR_URL, BIN_DIR / "PPRfinder.jar")

    # Download HMM Profiles
    for p in ['morf.hmm', 'dyw_kp.hmm', 'dyw.hmm', 'all_PPR.hmm']:
        download_file(f"{HMM_PROFILE_BASE_URL}{p}", PROFILES_DIR / p)

    # Download Prosite Data
    for p in ['prosite.dat', 'prorule.dat']:
        download_file(f"{PROSITE_URL}{p}", PROFILES_DIR / p)

    print("Resource setup complete.")

if __name__ == "__main__":
    main()
