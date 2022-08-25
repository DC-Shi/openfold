import argparse
import logging
import os

from openfold.data import mmcif_parsing
from openfold.np import protein, residue_constants

from tqdm import tqdm
from multiprocessing.pool import ThreadPool
from multiprocessing import Pool

# Retrieve a single page and report the URL and contents
def parseFile(names_path):
    basename, ext, fpath, _ = names_path
    ret = []
    if(ext == ".cif"):
        with open(fpath, 'r') as fp:
            mmcif_str = fp.read()
        
        mmcif = mmcif_parsing.parse(
            file_id=basename, mmcif_string=mmcif_str
        )
        if(mmcif.mmcif_object is None):
            logging.warning(f'Failed to parse {basename}.{ext}...')
            if(args.raise_errors):
                raise list(mmcif.errors.values())[0]
            else:
                # continue
                return []

        mmcif = mmcif.mmcif_object
        for chain, seq in mmcif.chain_to_seqres.items():
            chain_id = '_'.join([basename, chain])
            ret.append(f">{chain_id}")
            ret.append(seq)
    elif(ext == ".core"):
        with open(fpath, 'r') as fp:
            core_str = fp.read()

        core_protein = protein.from_proteinnet_string(core_str)
        aatype = core_protein.aatype
        seq = ''.join([
            residue_constants.restypes_with_x[aatype[i]] 
            for i in range(len(aatype))
        ])
        ret.append(f">{basename}")
        ret.append(seq)

    return ret
            

def main(args):
    fasta = []
    # Construct names and paths first
    names_paths = []
    count = 0
    for fname in os.listdir(args.data_dir):
        basename, ext = os.path.splitext(fname)
        basename = basename.upper()
        fpath = os.path.join(args.data_dir, fname)
        names_paths.append((basename, ext, fpath, count))
        #fasta.append(count)
        count = count + 1
    print(f"Searched all mmcif, total {count} files.")
        
    with Pool(processes=40) as pool:         # start 40 worker processes
        threadsIter = pool.imap(parseFile, names_paths, chunksize=10)     # Loop over each element in array
        with tqdm(total = len(names_paths)+1) as pbar:
            for result in threadsIter:
                fasta.extend(result)
                pbar.update(1)

    with open(args.output_path, "w") as fp:
        fp.write('\n'.join(fasta))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_dir", type=str,
        help="Path to a directory containing mmCIF or .core files"
    )
    parser.add_argument(
        "output_path", type=str,
        help="Path to output FASTA file"
    )
    parser.add_argument(
        "--raise_errors", type=bool, default=False,
        help="Whether to crash on parsing errors"
    )

    args = parser.parse_args()

    main(args)



## Running commands:
# python3 scripts/data_dir_to_fasta.py /data/pdb_mmcif/mmcif_files /data/example_DB_fasta
