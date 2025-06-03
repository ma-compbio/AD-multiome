from collections import defaultdict, Counter
import gzip, os, sys, time, pickle, itertools, copy, datetime, gc, importlib, re, h5py, json
from multiprocessing import Pool
from pathlib import Path
from tqdm.auto import tqdm
from collections.abc import Iterable
from typing import Dict, Union, List, Tuple

import numpy as np
import pandas as pd

from utils import merge_libraries, Library


def merge(obj_configs: List[Tuple[str, str, str, str]], out_lib_name: str, species: str):
	path2root = Path('/work/magroup/xinyuelu/AD-proj')
	print(f'path to root = {path2root}')

	for walk_policy in [
		'complete',
		# 'linear',
		# 'wowalk',
		# 'walkonly',
	]:
		print(f'walk policy = {walk_policy}')
		objs_single = {}
		for lib_name, hic_name, rna_name, s in obj_configs:
			path2dir = path2root / f'results'
			suffix = f'{lib_name}_{walk_policy}_{s}_filtered'
			obj = Library(hic_name, rna_name, lib_name, path2dir, wellid_list=path2dir / 'results' / f'meta_{suffix}.pkl')
			obj.load_RNA(path2dir / 'results' / f'RNA_{suffix}.pkl')
			obj.load_DNA(path2dir / 'results' / f'contact_{suffix}.pkl')
			obj.load_meta(path2dir / 'results' / f'meta_{suffix}.pkl')
			objs_single[lib_name] = obj
			del obj
		print('merging')
		obj = merge_libraries(objs_single, out_lib_name.format(walk_policy=walk_policy), path2root)
		print('writing raw data')
		path2dir = Path('/work/magroup/xinyuelu/AD-proj/results')
		obj.save_meta(path2dir / 'results' / f'meta_{obj.lib_name}.pkl')
		obj.save_meta(path2dir / 'results' / f'meta_{obj.lib_name}.csv')
		obj.save_RNA(path2dir / 'results' / f'RNA_{obj.lib_name}.pkl')
		obj.save_DNA(path2dir / 'results' / f'contact_{obj.lib_name}.pkl')


if __name__ == '__main__':
 
	merge(
		[
			('hPFC-lib-1', 'hPFC-scHiC-lib-1', 'hPFC-scRNA-lib-1', 'hg38'),
			('hPFC-lib-2', 'hPFC-scHiC-lib-2', 'hPFC-scRNA-lib-2', 'hg38'),
			('hPFC-lib-3', 'hPFC-scHiC-lib-3', 'hPFC-scRNA-lib-3', 'hg38'),
			('hPFC-lib-4', 'hPFC-scHiC-lib-4', 'hPFC-scRNA-lib-4', 'hg38'),
			('hPFC-lib-5', 'hPFC-scHiC-lib-5', 'hPFC-scRNA-lib-5', 'hg38'),
			('hPFC-lib-6', 'hPFC-scHiC-lib-6', 'hPFC-scRNA-lib-6', 'hg38'),
			('hPFC-lib-7', 'hPFC-scHiC-lib-7', 'hPFC-scRNA-lib-7', 'hg38'),
			('hPFC-lib-8', 'hPFC-scHiC-lib-8', 'hPFC-scRNA-lib-8', 'hg38')
		],
		'hPFC-first-8_complete_hg38_filtered_202409',
		'hg38',
	)

	# merge(
	# 	[
	# 		('hPFC-lib-9', 'hPFC-scHiC-lib-9', 'hPFC-scRNA-lib-9', 'hg38'),
	# 		('hPFC-lib-10', 'hPFC-scHiC-lib-10', 'hPFC-scRNA-lib-10', 'hg38'),
	# 		('hPFC-lib-11', 'hPFC-scHiC-lib-11', 'hPFC-scRNA-lib-11', 'hg38'),
	# 		('hPFC-lib-12', 'hPFC-scHiC-lib-12', 'hPFC-scRNA-lib-12', 'hg38'),
	# 		('hPFC-lib-13', 'hPFC-scHiC-lib-13', 'hPFC-scRNA-lib-13', 'hg38'),
	# 		('hPFC-lib-14', 'hPFC-scHiC-lib-14', 'hPFC-scRNA-lib-14', 'hg38')
	# 	],
	# 	'hPFC-last-12_complete_hg38_filtered_202409',
	# 	'hg38',
	# )
