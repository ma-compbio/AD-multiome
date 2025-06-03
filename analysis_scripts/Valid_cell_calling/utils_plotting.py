from collections import defaultdict, Counter
import os, pickle, itertools, copy, re, logging
from pathlib import Path
from tqdm.auto import tqdm
from util import config_logger
from collections.abc import Iterable
from typing import Dict, List, Tuple, Union, Optional

import numpy as np
import pandas as pd

import scanpy as sc
from anndata import AnnData

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns

from backend import subset_expr_by_gene, calc_log_bins


logger = config_logger(logging.getLogger(__name__))

def plot_mapping_statistics(self, thrs, spec, batch, show=True):
	palette = {spec: 'C1', 'empty': 'darkgray'}
	kwargs_scatter = dict(
		s=10, markers=dict((_, 'x') for _ in [spec, 'empty']),
		hue_order=[spec, 'empty'], palette=palette, data=self.meta,
	)

	def draw_thr_line(ax, mode, key, label, **kwargs):
		if key not in thrs: return
		thr = thrs[key]
		kwargs = dict(**kwargs, **dict(linestyle=':'))
		for thr in thr if isinstance(thr, Iterable) else [thr]:
			if mode in ['h', 'v']: getattr(ax, f'ax{mode}line')(thr, label=label.format(thr), **kwargs)
			elif mode in ['x', 'y']:
				ax.axline(
					(thr/(1-thr), 1) if mode == 'x' else (1, thr/(1-thr)),
					slope=1, label=label.format(thr), **kwargs,
				)
			else: raise NotImplementedError(mode)

	def plot_histogram(ax, lib, key):
		bins = np.geomspace(self.meta[key].min()+1, self.meta[key].max(), 51)
		ax.hist(self.meta[key], bins=bins)
		ax.set(title=f'{lib} - Count distribution', xlabel='# of contacts', ylabel='Cell counts', xscale='log')
		draw_thr_line(ax, 'v', f'{key}', c=palette['hg38'], label='#hg38={}')

	def plot_scatter_DNARNA(ax, hue):
		sns.scatterplot(x='DNA coverage all', y='RNA coverage mRNA', ax=ax, hue=hue, **kwargs_scatter)
		ax.set(title=f'Colored by {hue}', xscale='log', yscale='log', xlabel='# of Hi-C reads', ylabel='# of RNA reads')
		draw_thr_line(ax, 'v', 'DNA coverage all', '#hg38 DNA={}', c=palette['hg38'])
		# draw_thr_line(ax, 'v', 'DNA coverage all mm10', '#mm10 DNA={}', c=palette['mm10'])
		draw_thr_line(ax, 'h', 'RNA coverage mRNA', '#hg38 RNA={}', c=palette['hg38'])
		# draw_thr_line(ax, 'h', 'RNA coverage mRNA mm10', '#mm10 RNA={}', c=palette['mm10'])
		ax.legend(loc='lower left', bbox_to_anchor=(.01, .01))

	fig, axes = plt.subplots(3, 3, figsize=(18, 15))
	axes_iter = iter(axes.ravel())

	if self.lib_DNA is not None:
		ax = next(axes_iter)
		plot_histogram(ax, 'DNA', 'DNA coverage all')

		ax = next(axes_iter)
		sns.scatterplot(ax=ax, data=self.meta, x='DNA cis', y='DNA trans', s=10, hue='DNA species', palette=palette)
		ax.axline([1, 1], slope=1, c='grey', label='cis:trans=1', linestyle=':')
		ax.set(
			title='DNA - cis:trans', xlabel='# of cis contacts', ylabel='# of trans contacts', xscale='log', yscale='log')

		ax = next(axes_iter)
		sns.scatterplot(
			ax=ax, data=self.meta[self.meta['DNA species'] != 'empty'],
			x=self.meta['DNA mitosis band %'], y=self.meta['DNA mitosis short %'],
			s=10, hue='DNA species', palette=palette)
		ax.set(title='DNA - mitosis', xlabel='% of <=2Mb contacts', ylabel='% of 2-12Mb contacts')

	if self.lib_RNA is not None:
		plot_histogram(next(axes_iter), 'RNA', 'RNA coverage mRNA')

		ax = next(axes_iter)
		sns.scatterplot(
			ax=ax, data=self.meta, x='RNA coverage mRNA', y='RNA coverage gene',
			s=10, hue='RNA species', palette=palette)
		ax.set(title='RNA - gene:mRNA', xlabel='# of mRNA', ylabel='# of gene')
		ax.set(xscale='log', yscale='log')

		ax = next(axes_iter)
		sns.scatterplot(
			ax=ax, data=self.meta, x='RNA coverage gene', y='RNA coverage %mito',
			s=10, hue='RNA species', palette=palette)
		ax.set(title='RNA - %mito:gene', xlabel='# of genes', ylabel='% of mito mRNAs')
		ax.set(xscale='log', yscale='linear')

	if self.lib_DNA is not None and self.lib_RNA is not None:
		plot_scatter_DNARNA(next(axes_iter), 'DNA species')
		plot_scatter_DNARNA(next(axes_iter), 'RNA species')
		plot_scatter_DNARNA(next(axes_iter), 'species')

	i = 1
	for ax in axes_iter: ax.remove(); i += 1
	# assert i > 1
	for ax in axes.flat: ax.grid(True)
	for ax in axes.ravel()[:-i]:
		try: ax.get_legend().remove()
		except: pass
	ax = axes.ravel()[-i]
	ax.legend(loc='center left', bbox_to_anchor=(1.01, .5), fontsize=20)

	fig.subplots_adjust(top=.85, bottom=.15, left=.04, right=.98, wspace=.35, hspace=.35)
	ax = fig.add_axes([.05, .05, .9, .85], frameon=False)
	ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	ax.set_title(self.lib_name, fontsize=20)
	path2figdir = self.path2root / 'figures-stats-overall' / batch
	os.makedirs(path2figdir, exist_ok=True)
	# for suffix in ['png', 'pdf']:
	for suffix in ['png']:
		path2file = path2figdir / f'stats_{self.lib_name}.{suffix}'
		logger.info(f'Saving to {path2file}')
		plt.savefig(path2file, dpi=600)
	if show:
		plt.show()
	plt.close()


def plot_RNA_statistics(self, spec, batch, bins_cell=20, bins_gene=20, show=True):
	fig, axes = plt.subplots(2, 3, figsize=(16, 8))
	iter_ax = iter(axes.ravel())
	m = self.gene_expr

	for key in ['gene', 'mRNA']:
		ax = next(iter_ax)
		x = self.meta[f'RNA coverage {key}']
		sns.histplot(ax=ax, x=x, bins=calc_log_bins(x, bins=bins_cell))
		ax.set(title=f'RNA coverage {key} per cell', xlabel=f'# of {key}s', ylabel='# of cells', xscale='log', yscale='log')

	ax = next(iter_ax)
	sns.scatterplot(
		ax=ax, data=self.meta, x='RNA coverage mRNA', y='RNA coverage gene', s=10, hue='RNA species',
		palette={spec: 'C1', 'empty': 'darkgray'},
		hue_order=[spec, 'empty'],
	)
	ax.set(title='RNA coverage mRNA:gene per cell', xlabel='# of mRNAs', ylabel='# of genes')

	for w, text in [(None, 'cell'), (m.data, 'mRNA')]:
		ax = next(iter_ax)
		x = np.bincount(m.col, weights=w, minlength=self.num_genes)
		sns.histplot(ax=ax, x=x, bins=calc_log_bins(x, bins=bins_gene))
		ax.set(title=f'# of {text}:genes', xlabel=f'# of {text}s', ylabel='# of genes', xscale='log', yscale='log')

	ax = next(iter_ax)
	sns.scatterplot(
		ax=ax, data=self.meta, x='RNA coverage gene', y='RNA coverage %mito', s=10, hue='RNA species',
		palette={spec: 'C1', 'empty': 'darkgray'},
		hue_order=[spec, 'empty'],
	)
	ax.set(title='RNA coverage gene:%mito per cell', xlabel=f'# of genes', ylabel='% of mito mRNAs')

	for ax in axes.flat: ax.grid()

	fig.subplots_adjust(hspace=.4, wspace=.4, top=.9)
	fig.suptitle(self.lib_name, fontsize=20)
	path2figdir = self.path2root / 'figures-stats-rna' / batch
	os.makedirs(path2figdir, exist_ok=True)
	for suffix in ['png', 'pdf']:
		path2file = path2figdir / f'stats_{self.lib_name}.{suffix}'
		logger.info(f'Saving to {path2file}')
		plt.savefig(path2file, dpi=300)
	if show:
		plt.show()
	plt.close()


def plot_DNA_coverage(self, chrom_size, spec, batch, bins_cell=20, show=True):
	ncol = 6
	nrow = (self.num_chroms*2 + ncol - 1) // ncol
	fig, axes = plt.subplots(nrow, ncol, figsize=(16, 16/ncol*nrow), squeeze=False)
	iter_ax = iter(axes.ravel())
	for chrom, l in chrom_size.items():
		if not chrom.startswith(spec):
			continue
		chrom = chrom.split('_')[-1]
		if not chrom in self.chrom_i2a:
			continue
		for key, text in [('all', 'contacts'), ('nnz500kb', 'loci pairs')]:
			ax = next(iter_ax)
			x = self.meta[f'DNA coverage {key} {chrom}'] / l
			if np.all(x == 0):
				continue
			sns.histplot(ax=ax, x=x, bins=calc_log_bins(x, bins=bins_cell, method='log'))
			ax.set(title=f'{chrom} contacts {key}', xlabel=f'# of {text}', ylabel='# of cells', xscale='log', yscale='log')

	for ax in axes.flat: ax.grid()

	fig.subplots_adjust(hspace=.4, wspace=.4, top=.96)
	fig.suptitle(self.lib_name, y=.98)
	path2figdir = self.path2root / 'figures-stats-dna' / batch
	os.makedirs(path2figdir, exist_ok=True)
	# for suffix in ['png', 'pdf']:
	for suffix in ['png']:
		path2file = path2figdir / f'stats_{self.lib_name}.{suffix}'
		logger.info(f'Saving to {path2file}')
		plt.savefig(path2file, dpi=300)
	if show:
		plt.show()
	plt.close()
