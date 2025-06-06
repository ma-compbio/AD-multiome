import os
import torch
import numpy as np
import pickle
import pandas as pd

from torch.utils.data import Dataset, DataLoader
from utils.utils import read_DNAseq_tsv, read_Expre_tsv, read_Expre_mtx, read_1D_HiC, read_pbulk_exp, hic_h5_coo

from scipy.sparse import coo_matrix

# /** Training the Model from Scratch **/
class bulk_mBC(Dataset):
    def __init__(self, seq, exp):
        self.seq = seq
        self.exp = exp

    def __getitem__(self, index):
        seq = torch.tensor(self.seq[index])
        exp = torch.tensor(self.exp[index])

        return seq.long(), exp.float()
    
    def __len__(self):
        return self.seq.shape[0]

class bulk_mBC_hic1d(Dataset):
    def __init__(self, seq, exp, hic_1d):
        self.seq = seq
        self.exp = exp
        self.hic_1d = hic_1d

    def __getitem__(self, index):
        seq = torch.tensor(self.seq[index])
        exp = torch.tensor(self.exp[index])
        hic_1d = torch.tensor(self.hic_1d[index])

        return seq.long(), exp.float(), hic_1d.float()
    
    def __len__(self):
        return self.seq.shape[0]

# /** For Using Enformer's CNN **/   
class bulk_mBC_pretrain(Dataset):
    def __init__(self, seq, exp, hic_1d, hic_2d):
        self.seq = seq
        self.exp = exp
        self.hic_1d = hic_1d
        self.hic_2d = hic_2d

    def __getitem__(self, index):
        seq = self.seq[index]
        exp = torch.tensor(self.exp[index])
        if self.hic_1d != None and self.hic_2d != None:
            hic_1d = self.hic_1d[index]
            hic_2d = self.hic_2d[index]
            return seq.float(), exp.float(), hic_1d.float()/255, hic_2d
        elif self.hic_1d != None:
            hic_1d = self.hic_1d[index]
            return seq.float(), exp.float(), hic_1d.float()/255
        elif self.hic_2d != None:
            hic_2d = self.hic_2d[index]
            return seq.float(), exp.float(), hic_2d
        else:
            return seq.float(), exp.float()
    
    def __len__(self):
        return self.seq.shape[0]

class sc_mBC(Dataset):
    def __init__(self, seq, seq_indice, exp, hic_1d):
        self.seq = seq
        self.seq_indice = seq_indice
        self.exp = exp
        self.hic_1d = hic_1d

    def __getitem__(self, index):
        seq_indice = self.seq_indice[index]
        seq = self.seq[seq_indice]
        exp = torch.tensor(self.exp[index])
        hic_1d = self.hic_1d[index]

        return seq.float(), exp.float(), hic_1d.float()/255

    def __len__(self):
        return self.exp.shape[0]

class pbulk_mBC(Dataset):
    def __init__(self, seq, seq_indice, exp, hic_1d, hic_2d):
        self.seq = seq
        self.seq_indice = seq_indice
        self.exp = exp
        self.hic_1d = hic_1d
        self.hic_2d = hic_2d
    def __getitem__(self, index):
        seq_indice = self.seq_indice[index]
        seq = self.seq[seq_indice]
        exp = torch.tensor(self.exp[index])
        if self.hic_1d != None and self.hic_2d != None:
            hic_1d = self.hic_1d[index]
            hic_2d = self.hic_2d[index]
            return seq.float(), exp.float(), hic_1d.float()/255, hic_2d
        elif self.hic_1d != None:
            hic_1d = self.hic_1d[index]
            return seq.float(), exp.float(), hic_1d.float()/255
        elif self.hic_2d != None:
            hic_2d = self.hic_2d[index]
            return seq.float(), exp.float(), hic_2d
        else:
            return seq.float(), exp.float()
    def __len__(self):
        return self.exp.shape[0]

def collate_fn_pbulk(has_hic_1d, has_hic_2d):
    def collate_fn(batch):
        if has_hic_1d and has_hic_2d:
            seq, exp, hic_1d, hic_2d = map(list, zip(*batch))
            seq = torch.stack(seq, dim=0)
            exp = torch.stack(exp, dim=0)
            hic_1d = torch.stack(hic_1d, dim=0)
            return seq, exp, hic_1d, hic_2d
        elif has_hic_1d:
            seq, exp, hic_1d = map(list, zip(*batch))
            seq = torch.stack(seq, dim=0)
            exp = torch.stack(exp, dim=0)
            hic_1d = torch.stack(hic_1d, dim=0)
            return seq, exp, hic_1d
        elif has_hic_2d:
            seq, exp, hic_2d = map(list, zip(*batch))
            seq = torch.stack(seq, dim=0)
            exp = torch.stack(exp, dim=0)
            return seq, exp, hic_2d
        else:
            seq, exp = map(list, zip(*batch))
            seq = torch.stack(seq, dim=0)
            exp = torch.stack(exp, dim=0)
            return seq, exp
    return collate_fn

#/** Using Enformer's Pretrain CNN **/
def load_data_pbulk(path, seed, batch_size, num_workers, target_len, hic_1d, hic_2d, gt_mode='raw'):
    total_sequences = torch.load(os.path.join(path ,'sequence_vector.pt'))
    total_expressions = read_pbulk_exp(os.path.join(path, 'expression_cov_1024_200_celltypebulk.pkl'))
    gene_df = pd.read_csv(os.path.join(path, 'genes.tsv'), sep="\t")
    gene_num = gene_df.shape[0]
    if hic_1d:
        total_ab_score = read_1D_HiC(os.path.join(path, '1d-score-celltypebulk-10kb-ab_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_ins_score_25 = read_1D_HiC(os.path.join(path, '1d-score-celltypebulk-10kb-is-hw5_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_ins_score_50 = read_1D_HiC(os.path.join(path, '1d-score-celltypebulk-10kb-is-hw10_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_ins_score_100 = read_1D_HiC(os.path.join(path, '1d-score-celltypebulk-10kb-is-hw25_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_genebody = read_1D_HiC(os.path.join(path, '1d-score-celltypebulk-10kb-genebody-hw0_1024_200_uint8.pkl')).reshape(-1, 400, 1)

        total_ab_score = torch.from_numpy(total_ab_score)
        total_ins_score_25 = torch.from_numpy(total_ins_score_25)
        total_ins_score_50 = torch.from_numpy(total_ins_score_50)
        total_ins_score_100 = torch.from_numpy(total_ins_score_100)
        total_genebody = torch.from_numpy(total_genebody)

        total_1D_HiC = torch.concat((total_ab_score, total_ins_score_25, total_ins_score_50, total_ins_score_100, total_genebody), axis=2)
    if hic_2d:
        if os.path.exists(os.path.join(path, 'contact_1024_200_celltypebulk.pkl')):
            print('Read the contact map pickle')
            with open(os.path.join(path, 'contact_1024_200_celltypebulk.pkl'), 'rb') as f:
                total_2D_HiC = pickle.load(f)
        else:
            print('Generate the contact map pickle the first time')
            total_2D_HiC = hic_h5_coo(os.path.join(path, 'contact_1024_200_celltypebulk.h5'))
            with open(os.path.join(path, 'contact_1024_200_celltypebulk.pkl'), 'wb') as f:
                pickle.dump(total_2D_HiC, f)

    if gt_mode == 'raw': 
        # Normalize the Expression data
        row_min = np.min(total_expressions, axis=1, keepdims=True)
        row_max = np.max(total_expressions, axis=1, keepdims=True)
        total_expressions = (total_expressions - row_min) / (row_max - row_min)
        total_expressions = np.log1p(total_expressions * 1e4)
        # modify here
        # total_expressions = total_expressions[:, :240]
        # total_expressions = np.mean(total_expressions, axis=1, keepdims=True)
        total_expressions = total_expressions.reshape(-1, 400)
    elif gt_mode == 'binary':
        total_expressions = (total_expressions != 0).astype(int)
        total_expressions = total_expressions.reshape(-1, 400)

    # # crop the DNA-sequence from two sides
    trim = (target_len - total_expressions.shape[1]) // 2
    total_expressions = total_expressions[:, -trim:trim]

    # split the dataset
    # sample the train cell type
    torch.manual_seed(seed)
    perm = torch.randperm(26)
    # cell_indice_train = perm[:int(26*0.8)]
    # cell_indice_valid = perm[int(26*0.8):]
    cell_indice_train = [0,1,2,4,5,6,7,8,9,11,13,14,15,16,17,18,19,22,23,24]
    cell_indice_valid = [3,10,12,20,21,25]
    # cell_indice_train = [1,5,7,9,11,13,17,19,23,25]
    # cell_indice_valid = [3,15,21]
    # cell_indice_train = []
    # cell_indice_test = []
    # cell_indice_train = [(item - 1) for item in cell_indice_train]
    # cell_indice_valid = [(item - 1) for item in cell_indice_valid]
    cell_indice_train = torch.tensor(cell_indice_train)
    cell_indice_valid = torch.tensor(cell_indice_valid)
    cell_indice_train = cell_indice_train[torch.randperm(len(cell_indice_train))]
    cell_indice_valid = cell_indice_valid[torch.randperm(len(cell_indice_valid))]
    
    # shuffled_indices = torch.randperm(gene_num)
    # train_size = int(gene_num * 0.8)
    # seq_indice_train = shuffled_indices[:train_size]
    # seq_indice_valid = shuffled_indices[train_size:]
    seq_indice_train = torch.arange(int(gene_num*0.8))
    seq_indice_valid = torch.arange(int(gene_num*0.8), gene_num)

    train_indice = torch.cartesian_prod(cell_indice_train, seq_indice_train)
    train_indice = train_indice[:, 0] * gene_num + train_indice[:, 1]

    valid_indice_1 = torch.cartesian_prod(cell_indice_train, seq_indice_valid)
    valid_indice_1 = valid_indice_1[:, 0] * gene_num + valid_indice_1[:, 1]
    valid_indice_2 = torch.cartesian_prod(cell_indice_valid, seq_indice_train)
    valid_indice_2 = valid_indice_2[:, 0] * gene_num + valid_indice_2[:, 1]
    valid_indice = torch.concat((valid_indice_1, valid_indice_2), dim=0)

    test_indice = torch.cartesian_prod(cell_indice_valid, seq_indice_valid)
    test_indice = test_indice[:, 0] * gene_num + test_indice[:, 1]

    train_seq_indice = train_indice % gene_num
    valid_seq_indice = valid_indice % gene_num
    test_seq_indice = test_indice % gene_num

    # split the data
    train_exp = total_expressions[train_indice]
    valid_exp = total_expressions[valid_indice]
    test_exp = total_expressions[test_indice]

    train_1d_hic = None
    valid_1d_hic = None
    test_1d_hic = None 
    train_2d_hic = None
    valid_2d_hic = None
    test_2d_hic = None

    if hic_1d:
        train_1d_hic = total_1D_HiC[train_indice]
        valid_1d_hic = total_1D_HiC[valid_indice]
        test_1d_hic = total_1D_HiC[test_indice]
    if hic_2d:
        train_2d_hic = [total_2D_HiC[i] for i in train_indice]
        valid_2d_hic = [total_2D_HiC[i] for i in valid_indice]
        test_2d_hic = [total_2D_HiC[i] for i in test_indice]

    train_dataset = pbulk_mBC(total_sequences, train_seq_indice, train_exp, train_1d_hic, train_2d_hic)
    valid_dataset = pbulk_mBC(total_sequences, valid_seq_indice, valid_exp, valid_1d_hic, valid_2d_hic)
    test_dataset = pbulk_mBC(total_sequences, test_seq_indice, test_exp, test_1d_hic, test_2d_hic)

    train_loader = DataLoader(
        dataset = train_dataset, batch_size=batch_size, shuffle=True, num_workers=num_workers, collate_fn=collate_fn_pbulk(hic_1d, hic_2d))
    valid_loader = DataLoader(
        dataset = valid_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, collate_fn=collate_fn_pbulk(hic_1d, hic_2d))
    test_loader = DataLoader(
        dataset = test_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, collate_fn=collate_fn_pbulk(hic_1d, hic_2d))

    return train_loader, valid_loader, test_loader
        

def load_data_sc(path, seed, batch_size, num_workers, target_len, split=3):
    total_sequences = torch.load(os.path.join(path ,'sequence_vector.pt'))
    total_expressions = read_Expre_mtx(os.path.join(path, 'expression_cov_1024_200.mtx')).X.toarray()
    total_ab_score = read_1D_HiC(os.path.join(path, '1d-score-10kb-ab_1024_200_uint8.pkl')).reshape(-1, 400, 1)
    total_ins_score_25 = read_1D_HiC(os.path.join(path, '1d-score-10kb-is-hw25_1024_200_uint8.pkl')).reshape(-1, 400, 1)
    total_ins_score_50 = read_1D_HiC(os.path.join(path, '1d-score-10kb-is-hw50_1024_200_uint8.pkl')).reshape(-1, 400, 1)
    total_ins_score_100 = read_1D_HiC(os.path.join(path, '1d-score-10kb-is-hw100_1024_200_uint8.pkl')).reshape(-1, 400, 1)
    total_genebody = read_1D_HiC(os.path.join(path, '1d-score-10kb-genebody_1024_200_uint8.pkl')).reshape(-1, 400, 1)

    total_ab_score = torch.from_numpy(total_ab_score)
    total_ins_score_25 = torch.from_numpy(total_ins_score_25)
    total_ins_score_50 = torch.from_numpy(total_ins_score_50)
    total_ins_score_100 = torch.from_numpy(total_ins_score_100)
    total_genebody = torch.from_numpy(total_genebody)

    total_1D_HiC = torch.concat((total_ab_score, total_ins_score_25, total_ins_score_50, total_ins_score_100, total_genebody), axis=2)
    # total_1D_HiC = torch.concat((total_ab_score, total_ins_score_100, total_genebody), dim=2)

    # Normalize the Expression data
    row_min = np.min(total_expressions, axis=1, keepdims=True)
    row_max = np.max(total_expressions, axis=1, keepdims=True)
    total_expressions = (total_expressions - row_min) / (row_max - row_min)
    total_expressions = np.log1p(total_expressions * 1e4)
    total_expressions = total_expressions.reshape(-1, 400)

    # crop the DNA-sequence from two sides
    trim = (target_len - total_expressions.shape[1]) // 2
    total_expressions = total_expressions[:, -trim:trim]

    # transform the 1D HiC data
    # total_1D_HiC = np.concatenate((total_ab_score, total_ins_score_25, total_ins_score_50, total_ins_score_100, total_genebody), axis=2)

    # split the dataset
    cell_indice_train = torch.arange(int(3105/10))
    cell_indice_valid = torch.arange(int(3105/10), 3105)
    seq_indice_train = torch.arange(int(3740/4))
    seq_indice_valid = torch.arange(int(3740/4), 3740)

    train_indice = torch.cartesian_prod(cell_indice_train, seq_indice_train)
    train_indice = train_indice[:, 0] * 3740 + train_indice[:, 1]

    valid_indice_1 = torch.cartesian_prod(cell_indice_train, seq_indice_valid)
    valid_indice_1 = valid_indice_1[:, 0] * 3740 + valid_indice_1[:, 1]
    valid_indice_2 = torch.cartesian_prod(cell_indice_valid, seq_indice_train)
    valid_indice_2 = valid_indice_2[:, 0] * 3740 + valid_indice_2[:, 1]
    valid_indice = torch.concat((valid_indice_1, valid_indice_2), dim=0)

    # sample validation sample
    torch.manual_seed(seed)
    perm = torch.randperm(valid_indice.shape[0])
    sample_indice = perm[:int(valid_indice.shape[0]/100)]
    valid_indice = valid_indice[sample_indice]

    test_indice = torch.cartesian_prod(cell_indice_valid, seq_indice_valid)
    test_indice = test_indice[:, 0] * 3740 + test_indice[:, 1]

    # sample test sample
    torch.manual_seed(seed)
    perm = torch.randperm(test_indice.shape[0])
    sample_indice = perm[:int(test_indice.shape[0]/200)]
    test_indice = test_indice[sample_indice]

    train_seq_indice = train_indice % 3740
    valid_seq_indice = valid_indice % 3740
    test_seq_indice = test_indice % 3740

    # split the data
    train_exp = total_expressions[train_indice]
    train_1d_hic = total_1D_HiC[train_indice]

    valid_exp = total_expressions[valid_indice]
    valid_1d_hic = total_1D_HiC[valid_indice]

    test_exp = total_expressions[test_indice]
    test_1d_hic = total_1D_HiC[test_indice]

    train_dataset = sc_mBC(total_sequences, train_seq_indice, train_exp, train_1d_hic)
    valid_dataset = sc_mBC(total_sequences, valid_seq_indice, valid_exp, valid_1d_hic)
    test_dataset = sc_mBC(total_sequences, test_seq_indice, test_exp, test_1d_hic)

    train_loader = DataLoader(
        dataset = train_dataset, batch_size=batch_size, shuffle=True, num_workers=num_workers)
    valid_loader = DataLoader(
        dataset = valid_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers)
    test_loader = DataLoader(
        dataset = test_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers)

    return train_loader, valid_loader, test_loader

def load_data_bulk_pretrain(path, seed, batch_size, num_workers, target_len, hic_1d, hic_2d):
    total_sequences = torch.load(os.path.join(path,'sequence_vector.pt'))
    total_expressions = read_Expre_tsv(os.path.join(path, 'expression_cov_1024_200_bulk.tsv'))
    if hic_1d:
        total_ab_score = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-ab_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_ins_score_25 = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-is-hw5_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_ins_score_50 = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-is-hw10_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_ins_score_100 = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-is-hw25_1024_200_uint8.pkl')).reshape(-1, 400, 1)
        total_genebody = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-genebody_1024_200_uint8.pkl')).reshape(-1, 400, 1)

        total_ab_score = torch.from_numpy(total_ab_score)
        total_ins_score_25 = torch.from_numpy(total_ins_score_25)
        total_ins_score_50 = torch.from_numpy(total_ins_score_50)
        total_ins_score_100 = torch.from_numpy(total_ins_score_100)
        total_genebody = torch.from_numpy(total_genebody)

        total_1D_HiC = torch.concat((total_ab_score, total_ins_score_25, total_ins_score_50, total_ins_score_100, total_genebody), axis=2)
    if hic_2d:
        if os.path.exists(os.path.join(path, 'contact_1024_200_celltypebulk.pkl')):
            print('Read the contact map pickle')
            with open(os.path.join(path, 'contact_1024_200_celltypebulk.pkl'), 'rb') as f:
                total_2D_HiC = pickle.load(f)
        else:
            print('Generate the contact map pickle the first time')
            total_2D_HiC = hic_h5_coo(os.path.join(path, 'contact_1024_200_celltypebulk.h5'))
            with open(os.path.join(path, 'contact_1024_200_celltypebulk.pkl'), 'wb') as f:
                pickle.dump(total_2D_HiC, f)

    # Normalize the Expression data
    row_min = np.min(total_expressions, axis=1, keepdims=True)
    row_max = np.max(total_expressions, axis=1, keepdims=True)
    total_expressions = (total_expressions - row_min) / (row_max - row_min)
    total_expressions = np.log1p(total_expressions * 1e4)

    # crop the DNA-sequence from two sides
    trim = (target_len - total_expressions.shape[1]) // 2
    total_expressions = total_expressions[:, -trim:trim]

    # generate random indice
    k = int(total_sequences.shape[0]/20)
    indice = np.zeros((total_sequences.shape[0], 1))
    indice[:k] = 1
    indice[k:2*k] = 2
    np.random.seed(seed)
    indice = np.random.permutation(indice)
    train_indice = np.where(indice==0)[0]
    valid_indice = np.where(indice==2)[0]
    test_indice = np.where(indice==1)[0]

    # split the data
    train_seq = total_sequences[train_indice]
    train_exp = total_expressions[train_indice]

    valid_seq = total_sequences[valid_indice]
    valid_exp = total_expressions[valid_indice]

    test_seq = total_sequences[test_indice]
    test_exp = total_expressions[test_indice]

    train_1d_hic = None
    valid_1d_hic = None
    test_1d_hic = None 
    train_2d_hic = None
    valid_2d_hic = None
    test_2d_hic = None

    if hic_1d:
        train_1d_hic = total_1D_HiC[train_indice]
        valid_1d_hic = total_1D_HiC[valid_indice]
        test_1d_hic = total_1D_HiC[test_indice]
    if hic_2d:
        train_2d_hic = [total_2D_HiC[i] for i in train_indice]
        valid_2d_hic = [total_2D_HiC[i] for i in valid_indice]
        test_2d_hic = [total_2D_HiC[i] for i in test_indice]

    train_dataset = bulk_mBC_pretrain(train_seq, train_exp, train_1d_hic, train_2d_hic)
    valid_dataset = bulk_mBC_pretrain(valid_seq, valid_exp, valid_1d_hic, valid_2d_hic)
    test_dataset = bulk_mBC_pretrain(test_seq, test_exp, test_1d_hic, test_2d_hic)

    train_loader = DataLoader(
        dataset = train_dataset, batch_size=batch_size, shuffle=True, num_workers=num_workers, collate_fn=collate_fn_pbulk(hic_1d, hic_2d))
    valid_loader = DataLoader(
        dataset = valid_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, collate_fn=collate_fn_pbulk(hic_1d, hic_2d))
    test_loader = DataLoader(
        dataset = test_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, collate_fn=collate_fn_pbulk(hic_1d, hic_2d))

    return train_loader, valid_loader, test_loader

# /** Train the Whole Model From Scratch **/
def load_data_bulk_enf(path, seed, batch_size, num_workers, target_len):
    total_sequences = read_DNAseq_tsv(os.path.join(path, 'sequence_1024_200.tsv'))
    total_expressions = read_Expre_tsv(os.path.join(path, 'expression_cov_1024_200_bulk.tsv'))

    # # manually convert 1024-resolution to 128-resolution
    # total_expressions = np.repeat(total_expressions, 8).reshape(total_expressions.shape[0], -1)

    # crop the DNA-sequence from two sides
    trim = (target_len - total_expressions.shape[1]) // 2
    total_expressions = total_expressions[:, -trim:trim]

    # Normalize the Expression data
    row_min = np.min(total_expressions, axis=1, keepdims=True)
    row_max = np.max(total_expressions, axis=1, keepdims=True)
    total_expressions = (total_expressions - row_min) / (row_max - row_min)
    total_expressions = np.log1p(total_expressions * 1e4)

    # generate random indice
    k = int(total_sequences.shape[0]/20)
    indice = np.zeros((total_sequences.shape[0], 1))
    indice[:k] = 1
    indice[k:2*k] = 2
    np.random.seed(seed)
    indice = np.random.permutation(indice)
    train_indice = np.where(indice==0)[0]
    valid_indice = np.where(indice==2)[0]
    test_indice = np.where(indice==1)[0]

    # split the data
    train_seq = total_sequences[train_indice]
    train_exp = total_expressions[train_indice]

    valid_seq = total_sequences[valid_indice]
    valid_exp = total_expressions[valid_indice]

    test_seq = total_sequences[test_indice]
    test_exp = total_expressions[test_indice]

    train_dataset = bulk_mBC(train_seq, train_exp)
    valid_dataset = bulk_mBC(valid_seq, valid_exp)
    test_dataset = bulk_mBC(test_seq, test_exp)

    train_loader = DataLoader(
        dataset = train_dataset, batch_size=batch_size, shuffle=True, num_workers=num_workers)
    valid_loader = DataLoader(
        dataset = valid_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers)
    test_loader = DataLoader(
        dataset = test_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers)

    return train_loader, valid_loader, test_loader

def load_data_bulk_hcf(path, seed, batch_size, num_workers, target_len):
    
    total_sequences = read_DNAseq_tsv(os.path.join(path, 'sequence_1024_200.tsv'))
    total_expressions = read_Expre_tsv(os.path.join(path, 'expression_cov_1024_200_bulk.tsv'))
    total_ab_score = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-ab_1024_200.pkl')).reshape(-1, 400, 1)
    total_ins_score_25 = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-is-hw25_1024_200.pkl')).reshape(-1, 400, 1)
    total_ins_score_50 = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-is-hw50_1024_200.pkl')).reshape(-1, 400, 1)
    total_ins_score_100 = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-is-hw100_1024_200.pkl')).reshape(-1, 400, 1)
    total_genebody = read_1D_HiC(os.path.join(path, '1d-score-bulk-10kb-genebody_1024_200.pkl')).reshape(-1, 400, 1)

    # # manually convert 1024-resolution to 128-resolution
    # total_expressions = np.repeat(total_expressions, 8).reshape(total_expressions.shape[0], -1)

    # crop the DNA-sequence from two sides
    trim = (target_len - total_expressions.shape[1]) // 2
    total_expressions = total_expressions[:, -trim:trim]

    # Normalize the Expression data
    row_min = np.min(total_expressions, axis=1, keepdims=True)
    row_max = np.max(total_expressions, axis=1, keepdims=True)
    total_expressions = (total_expressions - row_min) / (row_max - row_min)
    total_expressions = np.log1p(total_expressions * 1e4)

    # transform the 1D HiC data
    total_1D_HiC = np.concatenate((total_ab_score, total_ins_score_25, total_ins_score_50, total_ins_score_100, total_genebody), axis=2)

    # generate random indice
    k = int(total_sequences.shape[0]/20)
    indice = np.zeros((total_sequences.shape[0], 1))
    indice[:k] = 1
    indice[k:2*k] = 2
    np.random.seed(seed)
    indice = np.random.permutation(indice)
    train_indice = np.where(indice==0)[0]
    valid_indice = np.where(indice==2)[0]
    test_indice = np.where(indice==1)[0]

    # split the data
    train_seq = total_sequences[train_indice]
    train_exp = total_expressions[train_indice]
    train_1d_hic = total_1D_HiC[train_indice]

    valid_seq = total_sequences[valid_indice]
    valid_exp = total_expressions[valid_indice]
    valid_1d_hic = total_1D_HiC[valid_indice]

    test_seq = total_sequences[test_indice]
    test_exp = total_expressions[test_indice]
    test_1d_hic = total_1D_HiC[test_indice]

    train_dataset = bulk_mBC_hic1d(train_seq, train_exp, train_1d_hic)
    valid_dataset = bulk_mBC_hic1d(valid_seq, valid_exp, valid_1d_hic)
    test_dataset = bulk_mBC_hic1d(test_seq, test_exp, test_1d_hic)

    train_loader = DataLoader(
        dataset = train_dataset, batch_size=batch_size, shuffle=True, num_workers=num_workers)
    valid_loader = DataLoader(
        dataset = valid_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers)
    test_loader = DataLoader(
        dataset = test_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers)

    return train_loader, valid_loader, test_loader