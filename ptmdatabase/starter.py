import os
import subprocess
import argparse

def starter():
    parser = argparse.ArgumentParser(description="Start the QCMSPyCloud with a given result data folder of MSPyCloud.")
    parser.add_argument("data_folder", type=str, help="Path to the result data folder.")
    
    # Removing --help since argparse already provides help functionality
    args = parser.parse_args()

    path = os.path.dirname(os.path.abspath(__file__))
    
    # Checking if the provided data_folder is valid
    if os.path.isdir(args.data_folder):
        subprocess.call(['streamlit', 'run', os.path.join(path, 'Home_page.py'), "--server.enableXsrfProtection", "false", args.data_folder])
    else:
        print(f"Provided data_folder path does not exist: {args.data_folder}")

if __name__ == "__main__":
    starter()
