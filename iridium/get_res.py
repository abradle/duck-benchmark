import os,glob,yaml
import numpy as np
import os,sys

import  matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt

  
def conv_to_dict(out_data):
  new_data = {}
  for key in out_data:                     
    this_data = out_data[key]
    new_data[key] = {}                                    
    for data in this_data:                    
      new_data[key][data[0]] = data[1] 
  return new_data

def get_wqb_simple(file_duck_dat):
    f = open(file_duck_dat,'r')
    data = []
    for line in f:
        a = line.split()
        data.append([float(a[1]), float(a[3]), float(a[5]), float(a[8])])
    f.close()
    data = np.array(data[1:])
    Work = data[:,3]
    Wqb_max = max(Work[400:])
    Wqb_min = min(Work[:400])
    Wqb_value = Wqb_max - Wqb_min
    return(Wqb_value, data, Wqb_min)

def get_Wqb_value_all(input_dir,data_dir):
    file_list = []
    for fil in os.listdir(input_dir):
        if fil[-3:] == 'dat':
            file_list.append(os.path.join(input_dir,fil))

    Wqb_values = []
    plt.figure(figsize = (7,7))
    for fil in file_list:
        Wqb_data = get_wqb_simple(fil)
        Wqb_values.append(Wqb_data[0])
        plt.plot(1*Wqb_data[1][:,0], Wqb_data[1][:,3]-Wqb_data[2])

    plt.xlabel('HB Distance (A)')
    plt.ylabel('Work (kcal/mol)')
    plt.savefig(data_dir+'.png')
    Wqb = min(Wqb_values)
    return(Wqb)


def convert_name(pdb_id,res_name):
   return pdb_id + "_" + res_name 


run_files = glob.glob("*/*/run.yaml")
checkpoints = {"ligand":"ANTECHAMBER_AM1BCC.AC","solvated": "complex_solvated.pdb", "merged": "complex.pdb",
     "equil": "equil.chk", "heating": "heating.csv","density": "density.csv","md_$i": "md_$i.chk"}
check_list = ["ligand","solvated", "merged", "equil", "heating","density","md_$i"]

base_data = {}
for x in open("baseline.csv").readlines():
  base_data[convert_name(x.split()[1],x.split()[2].strip())] = float(x.split()[0])
print(base_data)
out_d = {}
for yaml_file in run_files:
    base_dir = os.path.dirname(yaml_file)
    run_name = os.path.split(base_dir)[-1]
    run_dict = yaml.load(open(yaml_file))
    checkpoints_list = []
    for check in check_list:
        file_paths= []
        if "$i" in check:
            for i in range(int(run_dict["num_smd_cycles"])):
                check_p = checkpoints[check].replace("$i", str(i))
                file_p = os.path.join(base_dir,check_p)
                checkpoints_list.append((file_p,check.replace("$i", str(i))))
        else:
            check_p = checkpoints[check]
            file_p = os.path.join(base_dir,check_p)
            checkpoints_list.append((file_p, check))
    results = []
    for data in checkpoints_list: 
        file_p = data[0]
        check_p = data[1]
        if os.path.isfile(file_p):
            results.append((check_p,True))
        else:
            results.append((check_p,False))
    out_d[run_name] = results
    if run_name in base_data:
        results.append(("ref_wqb",base_data[run_name]))
    try:
        wqb = get_Wqb_value_all(base_dir,run_name)
        results.append(("wqb",wqb))
    except Exception as e:
        print(e)
        if os.path.isfile(run_name+'.png'): os.remove(run_name+'.png')
        results.append(("wqb",None))
new_d = conv_to_dict(out_d)

for key in new_d:
    header = list(new_d[key].keys())
header.append("title")
out_f = open("data.csv","w")
out_f.write(",".join(header)+"\n")
for key in new_d:
    w_list = []
    for head in header:
        if head =="title":continue
        w_list.append(str(new_d[key][head]))
    w_list.append(key)
    out_f.write(",".join(w_list)+"\n")


