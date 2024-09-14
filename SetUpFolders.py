import os
import shutil

datapath = "/Volumes/LAB DRIVE/240829NTA with folder structure"
dat_extension = ".dat"
xml_extension = ".xml"
remove = ["ConstantBinsTable_", "Videos_"]
for filename in os.listdir(datapath):
    if filename.endswith(dat_extension) == False or filename.startswith("."): continue
    filename = filename[:-len(dat_extension)]
    for remove_str in remove:
        if filename.startswith(remove_str) == False: continue
        foldername = filename[len(remove_str):]
    filepath = f"{datapath}/{filename}"
    folderpath = f"{datapath}/{foldername}/Results"
    os.makedirs(folderpath, exist_ok = True)
    shutil.move(f"{filepath}{dat_extension}", f"{folderpath}/{filename}{dat_extension}")
foldernames = [item for item in os.listdir(datapath) if os.path.isfile(f"{datapath}/{item}") == False]
for filename in os.listdir(datapath):
    if filename.endswith(xml_extension) == False or filename.startswith("."): continue
    filename = filename[:-len(xml_extension)]
    filepath = f"{datapath}/{filename}"
    for foldername in foldernames:
        print(filename, foldername)
        if filename.startswith(foldername):
            folderpath = f"{datapath}/{foldername}"
            shutil.move(f"{filepath}{xml_extension}", f"{folderpath}/{filename}{xml_extension}")
            break
for foldername in foldernames:
    with open(f"{datapath}/{foldername}/info.md", 'w') as file:
        pass
print("Done")