import os

# Get the current folder path
folder_path = os.getcwd()

# Loop over all .itp files in the current folder
for filename in os.listdir(folder_path):
    if filename.endswith(".itp"):
        file_path = os.path.join(folder_path, filename)

        # Read the file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Replace "lipid_section" with ";lipid_section"
        updated_lines = [line.replace('lipid_section', ';lipid_section') for line in lines]

        # Write the changes back to the file
        with open(file_path, 'w') as file:
            file.writelines(updated_lines)

print("Replacement complete.")