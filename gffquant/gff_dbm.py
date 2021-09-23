import sys
import gzip
from functools import lru_cache

from intervaltree import IntervalTree

class EmapperFormat:
	def __init__(self, query_field, fields, categories):
		self.query_field = query_field
		self.categories = dict(zip(fields, categories))
	def get_category(self, index):
		return self.categories.get(index)


EMAPPER_FORMATS = {
	"v1": EmapperFormat(
		0,
		(5, 6, 7, 9, 11),
		("Gene_Ontology_terms", "KEGG_ko", "BiGG_Reaction", "eggNOG_OGs", "COG_Functional_Category")
	),
	"v2": EmapperFormat(
		0,
		(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20),
		("Gene_Ontology_terms", "EC_number", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "eggNOG_OGs", "COG_Functional_Category")
	)
}


class GffDatabaseManager:

	def iterate(self):
		with self.db as db_stream:
			for line in db_stream:
				if not line.startswith("#"):
					line = line.strip().split("\t")
					categories = list()
					for i, col in enumerate(line):
						category = self.emapper_format.get_category(i)
						if category is not None:
							features = tuple(feature for feature in col.strip().split(",") if feature)
							if features:
								categories.append((category, features))
					if categories:
						yield line[self.emapper_format.query_field], (("strand", None), ) + tuple(categories)

	def iterate_old(self):
		header = None
		with self.db as db_stream:
			for line in db_stream:
				line = line.strip()
				if line.startswith("#"):
					header = line.strip("#").split("\t")
				else:
					line = line.split("\t")
					features = list()
					for feat_cat, subfeatures in zip(header[6:], line[6:]):
						subfeatures = tuple(sf for sf in subfeatures.strip().split(",") if sf)
						if subfeatures:
							features.append((feat_cat, subfeatures))
					yield line[0], (("strand", None),) + tuple(features)

	def _read_index(self, f):
		db_index = dict()
		for line in open(f, "rt"):
			line = line.strip().split("\t")
			db_index.setdefault(line[0], list()).append(list(map(int, line[1:3])))
		return db_index

	def __init__(self, db, reference_type, db_index=None, emapper_version="v2"):
		gz_magic = b"\x1f\x8b\x08"
		gzipped = open(db, "rb").read(3).startswith(gz_magic)
		_open = gzip.open if gzipped else open
		self.reference_type = reference_type
		if db_index:
			if gzipped:
				raise ValueError(f"Database {db} is gzipped. This doesn't work together with an index. Please unzip and re-index.")
			self.db_index = self._read_index(db_index)
			self.db = _open(db, "rt")
		elif self.reference_type == "domain": # bed
			self.db = dict()
			for line in _open(db, "rt"):
				line = line.strip().split("\t")
				self.db.setdefault(line[0], dict()).setdefault((int(line[1]), int(line[2])), list()).append(line[3])
		else:
			_open = gzip.open if gzipped else open
			self.db = _open(db, "rt")
			self.db_index = None

		self.emapper_format = EMAPPER_FORMATS.get(emapper_version)
		if not self.emapper_format:
			raise ValueError(f"Cannot find emapper parse instructions for version {emapper_version}.")

	@lru_cache(maxsize=4096)
	def _read_data(self, ref_id, include_payload=False):
		gff_annotation = dict()
		for offset, size in self.db_index.get(ref_id, list()):
			self.db.seek(offset)
			for line in self.db.read(size).strip("\n").split("\n"):
				if not line.startswith("#"):
					line = line.strip().split("\t")
					features = dict()
					if include_payload:
						features = (("strand", line[6]),)
						features += tuple((item.split("=")[0], tuple(sorted(item.split("=")[1].split(",")))) for item in line[8].strip().split(";"))
					key = (line[0], int(line[3]), int(line[4]) + 1)
					gff_annotation[key] = features
		if not gff_annotation and not include_payload:
			print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id), file=sys.stderr, flush=True)
		return gff_annotation

	@lru_cache(maxsize=4096)
	def _get_tree(self, ref, cache_data=False):
		if self.reference_type == "domain":
			return IntervalTree.from_tuples(sorted((start, end) for start, end in self.db.get(ref, dict())))

		return IntervalTree.from_tuples(sorted([key[1:] for key in self._read_data(ref, include_payload=cache_data)]))

	def get_data(self, ref, start, end):
		if self.reference_type == "domain":
			dom_features = self.db.get(ref, dict()).get((start, end), list())
			features = (("strand", None), ("ID", ref))
			features += (("domtype", tuple(dom_features)),)
			return features if dom_features else dict()

		return self._read_data(ref, include_payload=True).get((ref, start, end), dict())

	def get_overlaps(self, ref, start, end, cache_data=False):
		def calc_covered_fraction(start, end, interval):
			if interval.begin <= start <= end <= interval.end:
				return start, end

			if start < interval.begin:
				return interval.begin, min(end, interval.end)
			elif interval.end < end:
				return max(start, interval.begin), interval.end
			return start, end
		overlaps = self._get_tree(ref, cache_data=cache_data)[start:end]
		covered = [calc_covered_fraction(start, end, interval) for interval in overlaps]
		# print("OVL_DEBUG", ref, start, end, overlaps, covered)

		return overlaps, covered

	def clear_caches(self):
		print(self._read_data.cache_info(), flush=True)
		self._read_data.cache_clear()
		print(self._get_tree.cache_info(), flush=True)
		self._get_tree.cache_clear()
