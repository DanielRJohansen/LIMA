import os
import shutil

def ListFilenames(folder_path):
    # List all files in the folder
    filenames = os.listdir(folder_path)

    # Remove file extensions and return
    #return [os.path.splitext(filename)[0] for filename in filenames if            os.path.isfile(os.path.join(folder_path, filename))]
    return list(set(os.path.splitext(filename)[0] for filename in filenames if os.path.isfile(os.path.join(folder_path, filename))))


def copyFiles(sourceDir, destDir, names, extension):
    for name in names:
        srcFile = os.path.join(sourceDir, name + extension)
        destPath = os.path.join(destDir, name)
        if not os.path.exists(destPath):
            os.makedirs(destPath)
        shutil.copy2(srcFile, os.path.join(destPath, name + extension))

def findPdbFiles(directory):
    pdb_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.pdb'):
                pdb_files.append(os.path.join(root, file))
    return pdb_files

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

def MakeLipids():
    filenames = ListFilenames(r"F:\LIMA\SLipids_2020\itp_files")
    print(filenames)
    copyFiles(r"F:\LIMA\SLipids_2020\itp_files", r"C:\Users\Daniel\git_repo\LIMA\resources\Lipids", filenames, ".itp")
    copyFiles(r"F:\LIMA\SLipids_2020\itp_files", r"C:\Users\Daniel\git_repo\LIMA\resources\Lipids", filenames, ".pdb")

    pdbfiles = findPdbFiles(r"C:\Users\Daniel\git_repo\LIMA\resources\Lipids")
    for file in pdbfiles:
        processAtomLines(file)






