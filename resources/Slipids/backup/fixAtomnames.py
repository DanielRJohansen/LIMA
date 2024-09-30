import os

def processAtomLines(filePath):
    with open(filePath, 'r') as file:
        lines = file.readlines()

    with open(filePath, 'w') as file:
        for line in lines:
            if line.startswith("ATOM") and len(line) > 16 and line[12].isdigit():
                print(line)
                line = line[:11] + " " + line[13:16] + line[12] + line[16:]
                print(line)
                print("\n")
            file.write(line)

def find_and_process_pdb_files():
    # Get the current working directory
    current_dir = os.getcwd()

    # Loop through all files in the directory
    for file_name in os.listdir(current_dir):
        if file_name.endswith(".pdb"):
            file_path = os.path.join(current_dir, file_name)
            processAtomLines(file_path)

if __name__ == "__main__":
    find_and_process_pdb_files()
