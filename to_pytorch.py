"""
Converts raw ProteinNet records to Pytorch-saved dictionaries.
"""

from glob import glob
import torch
import os


def read_protein_from_file(file_pointer, include_tertiary=False):
    """ Modified from github.com/OpenProtein/openprotein:preprocessing.py on June 20, 2019."""
    dict_ = {}
    _dssp_dict = {'L': 0, 'H': 1, 'B': 2, 'E': 3, 'G': 4, 'I': 5, 'T': 6, 'S': 7}
    _mask_dict = {'-': 0, '+': 1}

    while True:
        next_line = file_pointer.readline()
        if next_line == '[ID]\n':
            id_ = file_pointer.readline()[:-1]
            dict_.update({'id': id_})
        elif next_line == '[PRIMARY]\n':
            primary = file_pointer.readline()[:-1]
            dict_.update({'primary': primary})
        elif next_line == '[EVOLUTIONARY]\n':
            evolutionary = []
            for residue in range(21): evolutionary.append(
                [float(step) for step in file_pointer.readline().split()])
            dict_.update({'evolutionary': evolutionary})
        elif next_line == '[SECONDARY]\n':
            secondary = list([_dssp_dict[dssp] for dssp in file_pointer.readline()[:-1]])
            dict_.update({'secondary': secondary})
        elif next_line == '[TERTIARY]\n' and include_tertiary:
            tertiary = []
            # 3 dimension
            for axis in range(3): tertiary.append(
                [float(coord) for coord in file_pointer.readline().split()])
            dict_.update({'tertiary': tertiary})
        elif next_line == '[MASK]\n':
            mask = list([_mask_dict[aa] for aa in file_pointer.readline()[:-1]])
            dict_.update({'mask': mask})
        elif next_line == '\n':
            return dict_
        elif next_line == '':
            return None


def main():
    """ Converts ProteinNet data in data/raw into Pytorch-saved dictionaries that contain non-tertiary data. """
    if not os.path.exists("data/torch/"):
        os.mkdir("data/torch/")
    input_files = glob("data/raw/*")
    for input_filename in input_files:
        print(input_filename)
        input_file = open(input_filename, "r")
        meta_dict = {}
        while True:
            next_protein = read_protein_from_file(input_file)
            if next_protein is None:
                break
            id_ = next_protein["id"]
            del next_protein["id"]
            meta_dict.update({id_: next_protein})
        torch.save(meta_dict, "data/torch/{}".format(os.path.basename(input_filename) + ".pt"))
        input_file.close()


if __name__ == "__main__":
    main()
