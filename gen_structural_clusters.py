import os
import sys 
import numpy as np  
from scipy.spatial.distance import cdist
import pandas as pd
import math
from pathlib import Path
import mdtraj as md
import Bio.PDB as bpdb
from sklearn.cluster import AgglomerativeClustering
import getopt
from glob import glob

def validate_numerical(parameter, val):

        try:
                float(val)
        except ValueError:
                print("ERROR: input value %s is not numerical [%s]" % (val, parameter))
                print("Please ensure this value is numerical and try again")
                sys.exit()


def get_fd_info(pdb_file, backbone_only, start_residue, end_residue):

    coords_list = []
    structure = bpdb.PDBParser(QUIET=True).get_structure('temp', pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_num = residue.get_full_id()[3][1]
                residue_name = residue.get_resname()
                if residue_name in ['GLU', 'ASP']:
                    charge = -1
                elif residue_name in ['LYS', 'ARG']:
                    charge = 1 
                else:
                    charge = 0 
                if residue_num >= start_residue and residue_num <= end_residue:
                    for atom in residue:
                        atom_type = atom.get_name()
                        atom_num = atom.get_serial_number()
                        coords = list(residue[atom_type].get_vector())
                        coords.append(residue_num)
                        coords.append(residue_name)
                        coords.append(atom_type)
                        coords.append(charge)
                        if backbone_only:
                            if atom_type == 'CA':
                                coords_list.append(coords)
                        else:
                            coords_list.append(coords)

    pdb_df = pd.DataFrame(data=coords_list, columns=["x", "y", "z", "residue_num", "residue_name", "atom", "charge"])
    pdb_df = pdb_df.round({'x': 2, 'y': 2, 'z': 2})

    pdb_df.to_csv('./pdb_info.csv', index=False)



def get_atom_idx(pdb_file, backbone_only, start_residue, end_residue):

    atom_num_list = []
    structure = bpdb.PDBParser(QUIET=True).get_structure('temp', pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_num = residue.get_full_id()[3][1]
                if residue_num >= start_residue and residue_num <= end_residue:
                    for atom in residue:
                        atom_type = atom.get_name()
                        atom_num = atom.get_serial_number()
                        if backbone_only:
                            if atom_type == 'CA':
                                atom_num_list.append(atom_num)
                        else:
                            atom_num_list.append(atom_num)

    return(list(np.array(atom_num_list)-1))

def get_pairwise_rmsd(traj_file_list, pdb_file, idr_start_residue, idr_end_residue, fd_start_residue, fd_end_residue, rmsd_ca_resolution, num_clusters, frame_stride):

    ca_atom_idx_idr = get_atom_idx(pdb_file, True, idr_start_residue, idr_end_residue)
    ca_atom_idx_fd = get_atom_idx(pdb_file, True, fd_start_residue, fd_end_residue)
    ca_atom_idx_all = get_atom_idx(pdb_file, True, min(idr_start_residue, fd_start_residue), max(idr_end_residue, fd_end_residue))

    atom_idx_idr = get_atom_idx(pdb_file, False, idr_start_residue, idr_end_residue)
    atom_idx_fd = get_atom_idx(pdb_file, False, fd_start_residue, fd_end_residue)
    atom_idx_all = get_atom_idx(pdb_file, False, 0, max(idr_end_residue, fd_end_residue))

    ca_atom_idx_subset = np.linspace(0,len(ca_atom_idx_idr)-1,rmsd_ca_resolution,dtype=int) #for calculating rmsd

    traj_data_list = [] 
    for i in range(0,len(traj_file_list)):
        print('loading traj %i' % (i+1))
        traj_data_list.append(md.load(traj_file_list[i], top = pdb_file))

    if len(traj_data_list) > 1:
        traj = md.join(traj_data_list, check_topology=False, discard_overlapping_frames=False)
    else:
        traj = traj_data_list[0]

    traj.restrict_atoms(atom_idx_all) #exclude hetatm

    traj.superpose(traj, atom_indices = atom_idx_fd)

    num_frames = traj.n_frames

    if frame_stride > num_frames:
        print('Frame Stride > Total # of Frames')
        sys.exit()

    dist_matrix = np.zeros((num_frames, num_frames))

    #print(ca_atom_idx_idr)
    #print(ca_atom_idx_subset)
    print("Generating Distance Matrix")
    for i in range(0, num_frames, frame_stride):
        #if i % 100 == 0:
            #print('Frame %i' % i)
        rmsd = md.rmsd(traj, traj, i, atom_indices = np.array(ca_atom_idx_idr)[ca_atom_idx_subset]) #when atom_indices != None, precentered is ignored
        dist_matrix[i,:] = rmsd


    clustering = AgglomerativeClustering(affinity='precomputed', n_clusters=num_clusters, linkage='complete').fit(dist_matrix)


    traj_xyz_fd = traj.xyz[0,atom_idx_fd,:]

    if (fd_start_residue != idr_start_residue) and (fd_end_residue != idr_end_residue):
        print("Saving PDB Coordinate Representation")
        pdb_df = pd.DataFrame(data=traj_xyz_fd, columns=["x", "y", "z"])
        pdb_df = pdb_df.round({'x': 2, 'y': 2, 'z': 2})
        pdb_df.to_csv('./pdb_coordinate_representation.csv', index=False)
    else:
        print("Saving (EMPTY) PDB Coordinate Representation")
        pdb_df = pd.DataFrame()
        pdb_df.to_csv('./pdb_coordinate_representation.csv', index=False)


    relevant_frames = np.arange(0,num_frames,frame_stride)
    traj_xyz_ca_idr = traj.xyz[relevant_frames[:,None],np.array(ca_atom_idx_idr)[None,:],:]
    traj_xyz_ca_idr_df = traj_xyz_ca_idr.reshape(traj_xyz_ca_idr.shape[0]*traj_xyz_ca_idr.shape[1],traj_xyz_ca_idr.shape[2])

    aa_num = np.tile(range(1,traj_xyz_ca_idr.shape[1]+1), traj_xyz_ca_idr.shape[0])
    aa_num = aa_num.reshape(aa_num.shape[0],1)

    cluster_labels = np.zeros((traj_xyz_ca_idr.shape[0]*traj_xyz_ca_idr.shape[1],1))
    frame_num = np.zeros((traj_xyz_ca_idr.shape[0]*traj_xyz_ca_idr.shape[1],1))
    for i in range(0,traj_xyz_ca_idr.shape[0]):
        start_idx = traj_xyz_ca_idr.shape[1]*i
        end_idx = start_idx + traj_xyz_ca_idr.shape[1]
        frame_num[start_idx:end_idx] = i+1
        cluster_labels[start_idx:end_idx] = clustering.labels_[i]


    print("Saving Cluster Output")
    output = np.hstack((traj_xyz_ca_idr_df,aa_num, frame_num, cluster_labels))
    output_df = pd.DataFrame(data=output, columns=["x", "y", "z", "aa_num", "frame_num", "cluster"])
    output_df = output_df.round({'x': 2, 'y': 2, 'z': 2})

    output_df.to_csv('./cluster_output.csv', index=False)


def main(argv):

        idr_start_residue = '' 
        idr_end_residue = ''
        fd_start_residue = '' 
        fd_end_residue = ''
        pdb_file = ''
        traj_file = '' 
        rmsd_ca_resolution = ''
        num_clusters = ''
        frame_stride = ''

        try:
            opts, args = getopt.getopt(argv,"hs:e:t:f:p:j:r:n:d:",["idr_start_residue=", "idr_end_residue=", "fd_start_residue=", "fd_end_residue=", "pdb_file=", "traj_file=", "rmsd_ca_resolution=", "num_clusters=", "frame_stride="])
        except getopt.GetoptError:
                print('Usage: gen_structural_clusters.py -s <idr_start_residue> -e <idr_end_residue> -t <fd_start_residue> -f <fd_end_residue> -p <pdb_file> -j <traj_file> -r <rmsd_ca_resolution> -n <num_clusters> -d <frame_stride>')
                sys.exit(2)
        for opt, arg in opts:
                if opt == '-h':
                        print('Usage: gen_structural_clusters.py -s <idr_start_residue> -e <idr_end_residue> -t <fd_start_residue> -f <fd_end_residue> -p <pdb_file> -j <traj_file> -r <rmsd_ca_resolution> -n <num_clusters> -d <frame_stride>')
                        sys.exit()
                elif opt in ("-s", "--idr_start_residue"):
                        idr_start_residue = arg 
                elif opt in ("-e", "--idr_end_residue"):
                        idr_end_residue = arg 
                elif opt in ("-t", "--fd_start_residue"):
                        fd_start_residue = arg                
                elif opt in ("-f", "--fd_end_residue"):
                        fd_end_residue = arg 
                elif opt in ("-p", "--pdb_file"):
                        pdb_file = arg 
                elif opt in ("-j", "--traj_file"):
                        traj_file = arg 
                elif opt in ("-r", "--rmsd_ca_resolution"):
                        rmsd_ca_resolution = arg 
                elif opt in ("-n", "--num_clusters"):
                        num_clusters = arg 
                elif opt in ("-d", "--frame_stride"):
                        frame_stride = arg 

        
        #idr_start_residue = 2
        #idr_end_residue = 31
        #fd_start_residue = 32 
        #fd_end_residue = 247
        #fd_start_residue = 'NA'
        #fd_end_residue = 'NA'

        #rmsd_ca_resolution = 10
        #num_clusters = 3
        #frame_stride = 10

        #pdb_file = '/Users/ishan/Holehouse/IDR_FD_wo_NCterm/tail_len30/GSn_GFPw15/equilibriation/coil_start/1/__START.pdb'
        #pdb_file = '/Users/ishan/Holehouse/IDR_only/GSnr15_1/equilibriation/coil_start/1/__START.pdb'
        #traj_file = '/Users/ishan/Holehouse/IDR_FD_wo_NCterm/tail_len30/GSn_GFPw15/equilibriation/coil_start/1/__traj.xtc'
        #traj_file = '/Users/ishan/Holehouse/IDR_FD_wo_NCterm/tail_len30/GSn_GFPw15/equilibriation/coil_start'
        #traj_file = '/Users/ishan/Holehouse/IDR_only/GSnr15_1/equilibriation/coil_start'

        if fd_start_residue == 'NA' and fd_end_residue != 'NA':
            print('FD Start set to NA but FD End was not')
            sys.exit()
        if fd_start_residue != 'NA' and fd_end_residue == 'NA':
            print('FD End set to NA but FD Start was not')
            sys.exit()

        trajectory_type = ''
        if fd_start_residue == 'NA':
            fd_start_residue = idr_start_residue
        if fd_end_residue == 'NA':
            fd_end_residue = idr_end_residue
            trajectory_type = 'idr_only'
            print("IDR only trajectories")
        else:
            trajectory_type = 'idr_fd'
            print("IDR + FD trajectories")

        validate_numerical('idr_start_residue', idr_start_residue)
        validate_numerical('idr_end_residue', idr_end_residue)
        validate_numerical('fd_start_residue', fd_start_residue)
        validate_numerical('fd_end_residue', fd_end_residue)
        validate_numerical('rmsd_ca_resolution', rmsd_ca_resolution)
        validate_numerical('num_clusters', num_clusters)
        validate_numerical('frame_stride', frame_stride)

        idr_start_residue = int(idr_start_residue)
        idr_end_residue = int(idr_end_residue)
        fd_start_residue = int(fd_start_residue)
        fd_end_residue = int(fd_end_residue)
        rmsd_ca_resolution = int(rmsd_ca_resolution)
        num_clusters = int(num_clusters)
        frame_stride = int(frame_stride)

        if rmsd_ca_resolution == 0:
            print('RSMD CA resolution cannot equal 0')
            sys.exit()


        if trajectory_type == 'idr_fd':
            residue_check_cterm = fd_start_residue < fd_end_residue and fd_end_residue < idr_start_residue and idr_start_residue < idr_end_residue 
            residue_check_nterm = fd_start_residue < fd_end_residue and fd_start_residue > idr_end_residue  and idr_start_residue < idr_end_residue 
            residue_check_linker = fd_start_residue < fd_end_residue and idr_start_residue > fd_start_residue  and idr_end_residue < fd_end_residue and idr_start_residue < idr_end_residue
        else:
            residue_check_cterm = idr_start_residue < idr_end_residue 
            residue_check_nterm = idr_start_residue < idr_end_residue
            residue_check_linker = idr_start_residue < idr_end_residue

        if residue_check_cterm == False and residue_check_nterm == False and residue_check_linker == False:
            print('Invalid Residue Ranges')
            sys.exit()

        if os.path.exists(pdb_file) == False:
            print('Invalid PDB File Path')
            sys.exit()

        if os.path.exists(traj_file) == False:
            print('Invalid Traj File Path')
            sys.exit()

        #filepath inputted
        if 'xtc' not in traj_file and 'dcd' not in traj_file:
            xtc_files = [y for x in os.walk(traj_file) for y in glob(os.path.join(x[0], '*.xtc'))]
            dcd_files = [y for x in os.walk(traj_file) for y in glob(os.path.join(x[0], '*.dcd'))]
            traj_file_list = xtc_files + dcd_files
        else:
            traj_file_list = [traj_file]

        get_fd_info(pdb_file, False, fd_start_residue, fd_end_residue)
        get_pairwise_rmsd(traj_file_list, pdb_file, idr_start_residue, idr_end_residue, fd_start_residue, fd_end_residue, rmsd_ca_resolution, num_clusters, frame_stride)



if __name__ == "__main__":
   main(sys.argv[1:])



