from collections import defaultdict, Counter
import os, time, pickle, itertools, copy, re, h5py, json, logging
from pathlib import Path
from tqdm.auto import tqdm
from util import config_logger
from collections.abc import Iterable
from typing import Dict, Tuple, Union, List, Optional

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix
from scipy.io import mmwrite
from umap import UMAP

import scanpy as sc
from anndata import AnnData

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns


logger = config_logger(logging.getLogger(__name__))


def load_chrom_size(genome='hg38_mm10'):
	with open(f'/work/magroup/xinyuelu/genome/{genome}/{genome}.chrom.sizes', 'r') as f:
		d = dict(line.strip().split() for line in f)
	d = {k: int(v) for k, v in d.items()}
	return d


def contact2sparse(m: np.ndarray, num_cells: int, num_chroms: int, l=None):
	if l is None:
		l = max(m['pos1'].max(), m['pos2'].max()) + 1
	row = np.ravel_multi_index((m['cell_id'], m['pos1']), dims=(num_cells, l))
	col = np.ravel_multi_index((m['chrom1'], m['pos2']), dims=(num_chroms, l))
	m = coo_matrix((m['count'].astype(int), (row, col)), shape=(num_cells * l, num_chroms * l))
	m = m.tocsr()
	return m, l


def sparse2contact(m: coo_matrix, l: int, num_cells: int, num_chroms: int, dtype: np.dtype):
	m = m.tocoo()
	cell_id, pos1 = np.unravel_index(m.row, shape=(num_cells, l))
	chrom, pos2 = np.unravel_index(m.col, shape=(num_chroms, l))
	assert np.all(pos1 <= pos2)
	c = np.empty(len(cell_id), dtype=dtype)
	assert pos1.max() <= np.iinfo(dtype.fields['pos1'][0]).max, (pos1.max(), dtype['pos1'])
	assert pos2.max() <= np.iinfo(dtype.fields['pos2'][0]).max, (pos2.max(), dtype['pos2'])
	assert m.data.max() <= np.iinfo(dtype.fields['count'][0]).max, (m.data.max(), dtype['count'])
	c['chrom1'] = c['chrom2'] = chrom
	c['pos1'] = pos1
	c['pos2'] = pos2
	c['cell_id'] = cell_id
	c['count'] = m.data
	return c


def binning(m: np.ndarray, r: float, num_cells: int, num_chroms: int):
	m = m.copy()
	dtype = m.dtype
	r = int(r)
	m = m[m['chrom1'] == m['chrom2']]
	assert np.all(m['pos1'] <= m['pos2'])
	m['pos1'] //= r
	m['pos2'] //= r
	m, l = contact2sparse(m, num_cells, num_chroms)
	m = m.tocsr()
	m.sum_duplicates()
	m = m.tocoo()
	m = sparse2contact(m, l, num_cells, num_chroms, dtype)
	m['pos1'] *= r
	m['pos2'] *= r
	return m


class Library:
	def __init__(self, lib_DNA: str, lib_RNA: str, lib_name: str, path2root: Path, wellid_list: Optional = None):
		self.lib_DNA = lib_DNA
		self.lib_RNA = lib_RNA
		self.lib_name = lib_name
		self.path2root = path2root

		if wellid_list is None:
			wellid_list = path2root / 'results' / 'wellid_list_all.tsv'
		if isinstance(wellid_list, str) or isinstance(wellid_list, Path):
			wellid_list = str(wellid_list)
			if wellid_list.endswith('pkl'):
				self.load_meta(wellid_list)
				# wellid_list = self.meta['well id'].values
				wellid_list = self.meta.index.values
			else:
				with open(wellid_list, 'r') as f:
					wellid_list = f.read().strip().split()
		else:
			wellid_list = wellid_list
		self.cell_i2a = np.array(wellid_list)
		self.cell_a2i = dict(_[::-1] for _ in enumerate(wellid_list))
		self.num_cells = len(self.cell_i2a)
		logger.info(f'# of cells: {self.num_cells}')

		self.gene_i2a = self.gene_a2i = self.gene_i2c = self.num_genes = self.gene_expr = self.gene_expr_norm = None
		self.chrom_i2a = self.chrom_a2i = self.num_chroms = self.contact = None
		self.meta = pd.DataFrame(index=pd.Series(self.cell_i2a, name='well id'))
		self.mask_gene = slice(None)

	def save_meta(self, path2file: Path):
		logger.info(f'Saving to {path2file}')
		if path2file.suffix == '.csv':
			self.meta.to_csv(path2file)
		elif path2file.suffix == '.pkl':
			with open(str(path2file), 'wb') as f:
				pickle.dump(self.meta, f)
		else:
			raise NotImplementedError()
		logger.info('Done')

	def load_meta(self, path2file: Path):
		logger.info(f'Loading {path2file}')
		path2file = Path(path2file)
		if path2file.suffix == '.csv':
			self.meta = pd.read_csv(path2file, index_col=0)
		elif path2file.suffix == '.pkl':
			with open(path2file, 'rb') as f:
				self.meta = pd.DataFrame(pickle.load(f))
		else:
			raise NotImplementedError()

	def load_RNA(self, path2file: Optional[Path] = None):
		logger.info(f'Loading {path2file}')
		with open(path2file, 'rb') as f:
			self.gene_expr, self.gene_i2a, self.gene_i2c = pickle.load(f)
		self.num_genes = len(self.gene_i2a)
		self.gene_a2i = dict(_[::-1] for _ in enumerate(self.gene_i2a))
		logger.info(
			f'Loaded {self.num_genes} genes, '
			f'{self.gene_expr.nnz} cell-gene, '
			f'and {self.gene_expr.sum()} mRNAs'
		)

	def save_RNA(self, path2file: Path):
		logger.info(f'Saving to {path2file}')
		with open(str(path2file), 'wb') as f:
			pickle.dump((self.gene_expr, self.gene_i2a, self.gene_i2c), f)
		logger.info('Done')

	def load_DNA(self, path2file: Optional[Path] = None):
		logger.info(f'Loading {path2file}')
		with open(path2file, 'rb') as f:
			self.contact, self.chrom_i2a = pickle.load(f)
		self.num_chroms = len(self.chrom_i2a)
		self.chrom_a2i = dict(_[::-1] for _ in enumerate(self.chrom_i2a))
		logger.info(f'Loaded {self.num_chroms} chroms and {len(self.contact)} contacts')

	def save_DNA(self, path2file: str):
		logger.info(f'Saving to {path2file}')
		with open(path2file, 'wb') as f:
			pickle.dump((self.contact, self.chrom_i2a), f)
		logger.info('Done')

	def calc_DNA_coverage(self):
		if self.contact is None: return
		logger.info('Calculating DNA coverage')
		m = self.contact
		shape = (self.num_cells, self.num_chroms)
		x1 = np.ravel_multi_index((m['cell_id'], m['chrom1']), dims=shape)
		x2 = np.ravel_multi_index((m['cell_id'], m['chrom2']), dims=shape)
		def fn(idx):
			c = np.bincount(x1.ravel()[idx], weights=m['count'][idx], minlength=np.prod(shape)) + \
				np.bincount(x2.ravel()[idx], weights=m['count'][idx], minlength=np.prod(shape))
			return c.reshape(*shape).astype(float) / 2
		cnt = {
			'all': fn(slice(None)),
			'inter': fn(self.contact['chrom1'] != self.contact['chrom2']),
		}
		cnt['intra'] = cnt['all'] - cnt['inter']
		for k, c in cnt.items():
			self.meta[f'DNA coverage {k}'] = c.sum(1)
			self.meta[[f'DNA coverage {k} {chrom}' for chrom in self.chrom_i2a]] = c
		logger.info('Done')

	def calc_DNA_coverage_nnz(self, resolution=int(5e5), resolution_text='500kb'):
		logger.info('Calculating DNA coverage nnz')
		shape = (self.num_cells, self.num_chroms)
		m = binning(self.contact, resolution, *shape)
		m = m[m['pos1'] != m['pos2']]
		data = np.ravel_multi_index((m['cell_id'], m['chrom1']), dims=shape)
		data = np.bincount(data, minlength=np.prod(shape)).reshape(*shape)
		self.meta[f'DNA coverage nnz{resolution_text}'] = data.sum(1)
		self.meta[[f'DNA coverage nnz{resolution_text} {chrom}' for chrom in self.chrom_i2a]] = data
		logger.info('Done')

	def calc_DNA_longshort(self, thr=1e7):
		if self.contact is None: return
		logger.info('Calculating DNA long/short interactions')
		m = self.contact
		m = m[m['chrom1'] == m['chrom2']]
		mask_short = np.abs(m['pos1'] - m['pos2']) <= thr
		mask_long = ~mask_short
		shape = (self.num_cells, self.num_chroms)
		t = np.ravel_multi_index((m['cell_id'], m['chrom1']), dims=shape)
		n_short = np.bincount(t[mask_short], weights=m['count'][mask_short], minlength=np.prod(shape)).reshape(*shape)
		n_long  = np.bincount(t[mask_long ], weights=m['count'][mask_long ], minlength=np.prod(shape)).reshape(*shape)
		self.meta['DNA long'] = n_long.sum(1)
		self.meta['DNA short'] = n_short.sum(1)
		self.meta['DNA long:short'] = n_long.sum(1) / np.maximum(n_short.sum(1) + n_long.sum(1), 1)
		self.meta[[f'DNA long {chrom}' for chrom in self.chrom_i2a]] = n_long
		self.meta[[f'DNA short {chrom}' for chrom in self.chrom_i2a]] = n_short
		self.meta[[f'DNA long:short {chrom}' for chrom in self.chrom_i2a]] = n_long / np.maximum(n_short + n_long, 1)
		logger.info('Done')

	def calc_DNA_cistrans(self):
		if self.contact is None: return
		logger.info('Calculating DNA cis/trans ratio')
		m = self.contact
		n_all = np.bincount(m['cell_id'], weights=m['count'], minlength=self.num_cells)
		mask = m['chrom1'] != m['chrom2']
		n_trans = np.bincount(m['cell_id'][mask], weights=m['count'][mask],  minlength=self.num_cells)
		n_cis = n_all - n_trans
		self.meta['DNA cis'] = n_cis
		self.meta['DNA trans'] = n_trans
		self.meta['DNA cis:trans'] = n_cis / np.maximum(n_trans, 1)
		logger.info('Done')

	def calc_DNA_mitosis(self):
		if self.contact is None: return
		logger.info('Calculating DNA mitosis band:short')
		m = self.contact
		shape = (self.num_cells, self.num_chroms)
		m = m[m['chrom1'] == m['chrom2']]
		x = np.ravel_multi_index((m['cell_id'], m['chrom1']), dims=shape)
		mask = np.abs(m['pos2'] - m['pos1']) <= 2e6
		n_short = np.bincount(x[mask], weights=m['count'][mask], minlength=np.prod(shape)).reshape(*shape)
		mask = ~mask & (np.abs(m['pos2'] - m['pos1']) <= 1.2e7)
		n_mitosis = np.bincount(x[mask], weights=m['count'][mask], minlength=np.prod(shape)).reshape(*shape)
		self.meta['DNA mitosis short'] = n_short.sum(1)
		self.meta['DNA mitosis band'] = n_mitosis.sum(1)
		self.meta['DNA mitosis short %'] = n_short.sum(1) / np.maximum(self.meta['DNA coverage all'], 1)
		self.meta['DNA mitosis band %'] = n_mitosis.sum(1) / np.maximum(self.meta['DNA coverage all'], 1)
		self.meta['DNA mitosis band:short'] = n_mitosis.sum(1) / np.maximum(n_short.sum(1), 1)
		logger.info('Done')

	def calc_RNA_coverage(self):
		if self.gene_expr is None: return
		logger.info('Calculating RNA coverage')
		m = self.gene_expr
		self.meta['RNA coverage gene'] = np.bincount(m.row, weights=None, minlength=self.num_cells)
		self.meta['RNA coverage mRNA'] = np.bincount(m.row, weights=m.data, minlength=self.num_cells)
		mask = np.array([chrom.endswith('M') for chrom in self.gene_i2c])
		idx = mask[m.col]
		self.meta['RNA coverage mito'] = np.bincount(m.row[idx], weights=m.data[idx], minlength=self.num_cells)
		self.meta['RNA coverage %mito'] = self.meta['RNA coverage mito'] / np.maximum(self.meta['RNA coverage mRNA'], 1)
		logger.info('Done')

	def calc_DNA_coverage_species(self, species_list=('hg38', 'mm10')):
		if self.contact is None: return
		logger.info('Calculating DNA coverage per species')
		m = self.contact
		for species in species_list:
			mask = np.array([_.startswith(species) for _ in self.chrom_i2a])
			idx = mask[m['chrom1']] & mask[m['chrom2']]
			self.meta[f'DNA coverage all {species}'] = \
				np.bincount(m['cell_id'][idx], weights=m['count'][idx], minlength=self.num_cells)
		logger.info('Done')

	def calc_RNA_coverage_species(self, species_list=(None, 'hg38', 'mm10')):
		if self.gene_expr is None: return
		logger.info('Calculating RNA coverage per species')
		m = self.gene_expr
		for species in species_list:
			if species is None:
				mask = np.ones(len(self.gene_i2c), dtype=bool)
				suffix = ''
			else:
				mask = np.array([chrom.startswith(species) for chrom in self.gene_i2c])
				suffix = f' {species}'
			idx = mask[m.col]
			self.meta[f'RNA coverage gene{suffix}'] = \
				np.bincount(m.row[idx], weights=None, minlength=self.num_cells)
			self.meta[f'RNA coverage mRNA{suffix}'] = \
				np.bincount(m.row[idx], weights=m.data[idx], minlength=self.num_cells)
			mask &= np.array([chrom.endswith('M') for chrom in self.gene_i2c])
			idx = mask[m.col]
			self.meta[f'RNA coverage mito{suffix}'] = \
				np.bincount(m.row[idx], weights=m.data[idx], minlength=self.num_cells)
			self.meta[f'RNA coverage %mito{suffix}'] = \
				self.meta[f'RNA coverage mito{suffix}'] / np.maximum(self.meta[f'RNA coverage mRNA{suffix}'], 1)
		logger.info('Done')    

	def calc_all_in_one(self):
		self.calc_DNA_coverage()
		self.calc_DNA_coverage_nnz()
		self.calc_DNA_longshort()
		self.calc_DNA_cistrans()
		self.calc_DNA_mitosis()
		# self.calc_DNA_coverage_species()
		# self.calc_RNA_coverage_species()
		self.calc_RNA_coverage()
		self.get_cell_label()

	def filter_cell(self, chrom_size, thrs: Dict[str, Union[int, float, tuple]], species_list=('mm10', 'hg38')):
		species_flags = {
			k: np.ones(self.num_cells, dtype=bool) for k in itertools.product(['DNA', 'RNA'], species_list)}
		doublet_flags = {k: np.ones(self.num_cells, dtype=bool) for k in ['DNA', 'RNA']}
		species_list = np.array(species_list)

		def fn(f: np.ndarray, values: np.ndarray, key: str):
			if key not in thrs:
				return
			thr = thrs[key]
			assert np.isfinite(values).all(), (np.isfinite(values).mean(), values)
			if isinstance(thr, tuple):
				f &= (values >= thr[0]) & (values <= thr[1])
			else:
				f &= values >= thr

		for s in species_list:
			if self.gene_expr is not None:
				logger.info(f'Filtering for RNA {s}')
				f = species_flags[('RNA', s)]
				g = doublet_flags['RNA']
				for key in [
					'RNA coverage gene','RNA coverage mRNA', 'RNA coverage %mito'
				]:
					if key not in self.meta:
						continue
					fn(f, self.meta[key], key)
					fn(g, self.meta[key], key)

			if self.contact is not None:
				logger.info(f'Filtering for DNA {s}')
				f = species_flags[('DNA', s)]
				g = doublet_flags['DNA']
				for key in ['DNA coverage all']:
					if key not in self.meta:
						continue
					fn(f, self.meta[key], key)
					fn(g, self.meta[key], key)
				f &= self.meta['DNA cis'] > self.meta['DNA trans']
				g &= self.meta['DNA cis'] > self.meta['DNA trans']
				for chrom, l in chrom_size.items():
					if not chrom.startswith(s):
						continue
					if f'DNA coverage nnz500kb %l' in thrs:
						fn(f, self.meta[f'DNA coverage nnz500kb {chrom}'] / l, f'DNA coverage nnz500kb %l')
						fn(g, self.meta[f'DNA coverage nnz500kb {chrom}'] / l, f'DNA coverage nnz500kb %l')
					if f'DNA coverage nnz500kb %l {s}' in thrs:
						fn(f, self.meta[f'DNA coverage nnz500kb {chrom}'] / l, f'DNA coverage nnz500kb %l {s}')
						fn(g, self.meta[f'DNA coverage nnz500kb {chrom}'] / l, f'DNA coverage nnz500kb %l {s}')
					if f'DNA coverage all %l' in thrs:
						fn(f, self.meta[f'DNA coverage all {chrom}'] / l, f'DNA coverage all %l')
						fn(g, self.meta[f'DNA coverage all {chrom}'] / l, f'DNA coverage all %l')
					if f'DNA coverage all %l {s}' in thrs:
						fn(f, self.meta[f'DNA coverage all {chrom}'] / l, f'DNA coverage all %l {s}')
						fn(g, self.meta[f'DNA coverage all {chrom}'] / l, f'DNA coverage all %l {s}')

		for lib in ['RNA', 'DNA']:
			f = np.stack([species_flags[(lib, s)] for s in species_list])
			g = doublet_flags[lib]
			self.meta[f'{lib} species'] = np.where(
				~f.any(0),
				np.where(g, 'doublet', 'empty'),
				np.where(f.sum(0) == 1, species_list[np.argmax(f, 0)], 'doublet')
			)

		if self.lib_RNA is None:
			self.meta['species'] = self.meta['DNA species']
		elif self.lib_DNA is None:
			self.meta['species'] = self.meta['RNA species']
		else:
			self.meta['species'] = np.where(
				self.meta['RNA species'] == self.meta['DNA species'], self.meta['RNA species'],
				np.where(
					(self.meta['RNA species'] == 'doublet') | (self.meta['DNA species'] == 'doublet'),
					'doublet', 'empty'),
			)

		logger.info(self.meta.groupby(['species', 'DNA species', 'RNA species']).size())

	def subset(
			self, cell_mask=None, chrom_mask=None, species=None,
			columns2keep: Union[str, Tuple[str]] = ('DNA species', 'RNA species', 'species'),
	):
		if cell_mask is None:
			cell_mask = np.ones(self.num_cells, dtype=bool)
		cell_mapping = np.cumsum(cell_mask)-1
		self.num_cells = cell_mask.sum()
		self.cell_i2a = self.cell_i2a[cell_mask]
		self.cell_a2i = dict(_[::-1] for _ in enumerate(self.cell_i2a))

		if self.gene_expr is not None:
			m = self.gene_expr
			i = cell_mask[m.row]
			if species is not None:
				gene_mask = np.array([_.split('_')[0] in species for _ in self.gene_i2c])
				i &= gene_mask[m.col]
			gene_mask = np.bincount(m.col[i], minlength=self.num_genes) > 0
			gene_mapping = np.cumsum(gene_mask)-1
			self.num_genes = gene_mask.sum()
			self.gene_i2a = self.gene_i2a[gene_mask]
			self.gene_a2i = dict(_[::-1] for _ in enumerate(self.gene_i2a))
			self.gene_i2c = self.gene_i2c[gene_mask]
			self.gene_expr = coo_matrix(
				(m.data[i], (cell_mapping[m.row[i]], gene_mapping[m.col[i]])), shape=(self.num_cells, self.num_genes))
			del m, i

		if self.contact is not None:
			if chrom_mask is None: chrom_mask = np.ones(self.num_chroms, dtype=bool)
			chrom_mapping = np.cumsum(chrom_mask)-1
			self.chrom_i2a = self.chrom_i2a[chrom_mask]
			self.chrom_a2i = dict(_[::-1] for _ in enumerate(self.chrom_i2a))
			self.num_chroms = chrom_mask.sum()
			m = self.contact
			self.contact = m[cell_mask[m['cell_id']] & chrom_mask[m['chrom1']] & chrom_mask[m['chrom2']]]
			self.contact['cell_id'] = cell_mapping[self.contact['cell_id']]
			assert np.all(self.contact['cell_id'] >= 0)
			self.contact['chrom1'] = chrom_mapping[self.contact['chrom1']]
			self.contact['chrom2'] = chrom_mapping[self.contact['chrom2']]
			del m

		if self.meta is not None:
			self.meta = self.meta.iloc[cell_mask]
			if columns2keep == 'all':
				logger.info('Keeping all columns')
				pass
			elif columns2keep == 'none':
				logger.info('Removing all columns')
				self.meta = pd.DataFrame(index=self.meta.index)
			elif isinstance(columns2keep, tuple):
				cols = np.array(columns2keep)
				cols = cols[np.isin(cols, self.meta.columns)]
				logger.info(f'Keeping columns `{cols}`')
				self.meta = self.meta[cols]
			for col in self.meta.columns:
				if self.meta[col].dtype.name == 'category':
					self.meta[col] = self.meta[col].cat.remove_unused_categories()

def collect_singlecell_stats(objs: Dict[str, Library], lib_name: str, batch: str, path2root: Path, species_list=('hg38', 'mm10')):
	metas = {k: o.meta for k, o in objs.items()}
	obj = objs['complete']
	meta = metas['complete']
	df = pd.DataFrame(index=meta.index)
	path2dir = path2root / 'alignment' / batch

	# -- DNA --
	if obj.lib_DNA is not None:
		df['# of raw HiC read pairs'] = pd.read_csv(
			path2dir / f'statSingleCell_nHiCRawReadPair_{obj.lib_DNA}_all_hg38_mm10.txt',
			header=None, index_col=0, sep='\t',
		)
		df['# of uniquely mapped HiC read pairs'] = pd.read_csv(
			path2dir / f'statSingleCell_nHiCUniquelyMappedAllReadPair_{obj.lib_DNA}_all_hg38_mm10.txt',
			header=None, index_col=0, sep='\t',
		)
		df['# of dedup HiC read pairs'] = pd.read_csv(
			path2dir / f'statSingleCell_nHiCDedupAllReadPair_{obj.lib_DNA}_all_hg38_mm10.txt',
			header=None, index_col=0, sep='\t',
		)
		df['# of HiC cis-contacts'] = meta['DNA cis']
		df['# of HiC trans-contacts'] = meta['DNA trans']
		df['# of HiC cis:trans'] = meta['DNA cis:trans']
		df['# of HiC walk contacts'] = meta['DNA coverage all'] - metas['wowalk']['DNA coverage all']
	else:
		df['# of raw HiC read pairs'] = 0
		df['# of uniquely mapped nonwalk HiC read pairs'] = 0
		df['# of uniquely mapped walk HiC read pairs'] = 0
		df['# of uniquely mapped HiC read pairs'] = 0
		for s in species_list:
			df[f'# of HiC contacts {s}'] = 0
		df['# of HiC cis-contacts'] = 0
		df['# of HiC trans-contacts'] = 0
		df['# of HiC cis:trans'] = 0
		df['# of HiC walk contacts'] = 0
  
	# -- RNA --
	if obj.lib_RNA is not None:
		df['# of raw RNA reads'] = pd.read_csv(
			path2dir / f'statSingleCell_nRNARawRead_{obj.lib_RNA}_all_hg38_mm10.txt',
			header=None, index_col=0, sep='\t',
		)
		df['# of uniquely mapped RNA reads'] = pd.read_csv(
			path2dir / f'statSingleCell_nRNAUniquelyMappedRead_{obj.lib_RNA}_all_hg38_mm10.txt',
			header=None, index_col=0, sep='\t',
		)
		df[f'# of unique RNA molecules'] = meta[f'RNA coverage mRNA']
		df[f'# of genes'] = meta[f'RNA coverage gene']
		df[f'% of mitochondrial mRNA molecules'] = meta[f'RNA coverage %mito']
	else:
		df['# of raw RNA reads'] = 0
		df['# of uniquely mapped RNA reads'] = 0
		for s in species_list:
			df[f'# of unique RNA molecules {s}'] = 0
		for s in species_list:
			df[f'# of genes {s}'] = 0
		for s in species_list:
			df[f'% of mitochondrial mRNA molecules {s}'] = 0

	for col in df.columns:
		if np.all(np.rint(df[col]) == df[col]):
			df[col] = df[col].astype(int)

	# -- species --
	df['species'] = meta['species']

	df = df.fillna(0)
	path2file = path2root / 'stats' / f'statSingleCell_{lib_name}.csv'
	logger.info(f'Saving to {path2file}')
	df.to_csv(path2file)
	return df

def merge_libraries(objs: Dict[str, Library], lib_name: str, path2root: Path):
	wellid_list = [f'{batch}_{wellid}' for batch, obj in objs.items() for wellid in obj.meta.index.values]
	self = Library('none', 'none', lib_name, path2root, wellid_list=wellid_list)
	cell_offset = np.cumsum([0] + [obj.num_cells for obj in objs.values()])
	self.num_cells = cell_offset[-1]
	self.cell_i2a = np.array(wellid_list)
	self.cell_a2i = dict(_[::-1] for _ in enumerate(self.cell_i2a))

	self.gene_i2a = np.unique(np.concatenate([obj.gene_i2a for obj in objs.values()]))
	self.gene_a2i = dict(_[::-1] for _ in enumerate(self.gene_i2a))
	self.num_genes = len(self.gene_i2a)
	self.gene_i2c = np.array([
		next(
			obj.gene_i2c[obj.gene_a2i[g]]
			for obj in objs.values() if g in obj.gene_a2i
		)
		for g in self.gene_i2a
	])
	data, row, col = [], [], []
	for obj, co in zip(objs.values(), cell_offset):
		m = obj.gene_expr
		data.append(m.data)
		row.append(m.row + co)
		gene_mapping = np.searchsorted(self.gene_i2a, obj.gene_i2a)
		assert np.all(self.gene_i2a[gene_mapping] == obj.gene_i2a)
		col.append(gene_mapping[m.col])
	data, row, col = list(map(np.concatenate, [data, row, col]))
	self.gene_expr = coo_matrix((data, (row, col)), shape=(self.num_cells, self.num_genes))

	self.chrom_i2a = np.unique(np.concatenate([obj.chrom_i2a for obj in objs.values()]))
	self.chrom_a2i = dict(_[::-1] for _ in enumerate(self.chrom_i2a))
	self.num_chroms = len(self.chrom_i2a)
	num_contacts = np.cumsum([0] + [len(obj.contact) for obj in objs.values()])
	self.contact = np.concatenate([obj.contact for obj in objs.values()])
	for obj, co, l, r in zip(objs.values(), cell_offset[:-1], num_contacts[:-1], num_contacts[1:]):
		c = self.contact[l: r]
		c['cell_id'] += co
		chrom_mapping = np.searchsorted(self.chrom_i2a, obj.chrom_i2a)
		assert np.all(self.chrom_i2a[chrom_mapping] == obj.chrom_i2a)
		c['chrom1'] = chrom_mapping[c['chrom1']]
		c['chrom2'] = chrom_mapping[c['chrom2']]

	self.meta = []
	for batch, obj in objs.items():
		df = obj.meta.copy()
		df['batch'] = batch
		df.index = pd.Series([f'{batch}_{wellid}' for wellid in df.index], name='well id')
		self.meta.append(df)
	self.meta = pd.concat(self.meta)
	assert all(self.meta.index.values == wellid_list)
	self.meta['batch'] = pd.Categorical(self.meta['batch'], categories=list(objs.keys()))
	return self











