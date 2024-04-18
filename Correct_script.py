import os

# specify the directory containing the folders
root_dir = '.'

# iterate through all subdirectories in root_dir
for subdir, dirs, files in os.walk(root_dir):
    for file in files:
        # check if file has the .dat extension
        if file.endswith('.dat'):
            print("Rewriting file: " + file)
            file_path = subdir + os.sep + file
            with open(file_path, 'r') as f:
                # read the file and store the content in a list
                lines = f.readlines()
                # remove the first 10 lines
                lines = lines[10:]
                for i in range(len(lines)):
                    # divide the third column by 1000
                    columns = lines[i].split()
                    columns[2] = str(float(columns[2])/1000)
                    event_name = columns[0]
                    element = columns[1]
                    energy = columns[2]
                    x = columns[3]
                    y = columns[4]
                    z = columns[5]
                    cosx = columns[6]
                    cosy = columns[7]
                    cosz = columns[8]
                    # Determine the number of spaces to insert
                    if len(event_name) == 1:
                        spaces = "      "
                    elif len(event_name) == 2:
                        spaces = "     "
                    elif len(event_name) == 3:
                        spaces = "    "
                    elif len(event_name) == 4:
                        spaces = "   "
                    elif len(event_name) == 5:
                        spaces = "  "
                    # rebuild the line with the original formatting
                    lines[i] = event_name + spaces + element + " " + energy + " " + x + " " + y + " " + z + " " + cosx + " " + cosy + " " + cosz + "\r\n"
            with open(file_path, 'w') as f:
                # write the modified content back to the file
                f.writelines(lines)