import subprocess
import sys
import os
import platform
import shutil

def main():
    # Check for PyInstaller
    if not shutil.which("pyinstaller"):
        print("PyInstaller not found. Please install it first:")
        print("uv pip install pyinstaller")
        sys.exit(1)

    # Create temporary entry point
    entry_script = "run_t_level.py"
    print(f"Creating temporary entry point: {entry_script}")
    with open(entry_script, "w") as f:
        f.write("from t_level import main\n\nif __name__ == '__main__':\n    main()\n")

    # Determine separator for --add-data
    sep = ";" if platform.system() == "Windows" else ":"
    
    # Path to database
    db_source = os.path.join("t_level", "ecmprobs.db")
    db_dest = "t_level"
    add_data_arg = f"{db_source}{sep}{db_dest}"

    common_args = [
        "pyinstaller",
        "--noconfirm",
        "--clean",
        "--add-data", add_data_arg,
        entry_script
    ]

    # Build Directory Version (Standard)
    # Output: dist/dir/t-level/t-level
    print("\nBuilding Directory Version...")
    subprocess.run(common_args + ["--onedir", "--name", "t-level", "--distpath", os.path.join("dist", "dir")], check=True)

    # Build Single File Version (Portable)
    # Output: dist/portable/t-level
    print("\nBuilding Single File Version (Portable)...")
    subprocess.run(common_args + ["--onefile", "--name", "t-level", "--distpath", os.path.join("dist", "portable")], check=True)

    # Cleanup
    print("\nCleaning up...")
    if os.path.exists(entry_script):
        os.remove(entry_script)
    
    for spec in ["t-level.spec"]:
        if os.path.exists(spec):
            os.remove(spec)

    print("\n" + "="*50)
    print("Build Complete!")
    print("="*50)
    
    dist_dir = os.path.join("dist")
    dir_exe = os.path.join(dist_dir, "dir", "t-level", "t-level")
    portable_exe = os.path.join(dist_dir, "portable", "t-level")
    
    if platform.system() == "Windows":
        dir_exe += ".exe"
        portable_exe += ".exe"

    print(f"Directory version:   {dir_exe}")
    print(f"Portable version:    {portable_exe}")
    print("="*50)

if __name__ == "__main__":
    main()