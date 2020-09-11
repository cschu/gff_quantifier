import gzip
import sys
import argparse


class GffIndexer:
	def __init__(self, gff, gzipped=True):
		self._refs = dict()
		f_open = gzip.open if gzipped else open
		offset = 0
		cur_ref = None
		with f_open(gff, "rt") as f:
			for line in f:
				if not line.startswith("#"):
					ref = line.split("\t")[0]
					if ref != cur_ref:
						#print(self._refs)
						if cur_ref is not None:
							self._refs[cur_ref][-1][1] = offset - self._refs[cur_ref][-1][0]
							print(cur_ref, *self._refs[cur_ref][-1], sep="\t", flush=True)
						cur_ref = ref
						self._refs.setdefault(cur_ref, list()).append([offset, 0])
				offset += len(line)

			self._refs[cur_ref][-1][1] = offset - self._refs[cur_ref][-1][0]
			print(cur_ref, *self._refs[cur_ref][-1], sep="\t", flush=True)

	def get_index(self, out=sys.stdout):
		for r, p in self._refs.items():
			for pp in p:
				print(r, *pp, file=out, sep="\t", flush=True)

def index_gff(gff, gzipped=True, chunksize=1024):
	f_open = gzip.open if gzipped else open
	refs = dict()
	cur_ref = None
	offset = 0
	with f_open(gff, "rb") as f:
		while True:
			try:
				chunk = f.read(chunksize)
			except:
				break

			chunk = chunk.decode()

			#print(chunk)
			p = chunk.rfind("\n")
			if not chunk.endswith("\n") and p != -1:
				chunk, part = chunk[:p], chunk[p + 1:]
				if part.strip():
					f.seek(-len(part), 1)
			#print(chunk)
			#print(part)
			if not chunk.strip():
				break

			for line in chunk.split("\n"):
				print(line, len(line) + 1)
				if not line.startswith("#"):
					ref = line.split("\t")[0]
					if ref != cur_ref:
						print(refs)
						if cur_ref is not None:
							refs[cur_ref][1] = offset - refs[cur_ref][0]
						cur_ref = ref
						refs[cur_ref] = [offset, 0]
				offset += len(line) + 1
			#break
		refs[cur_ref][1] = offset - refs[cur_ref][0]

	for r, p in refs.items():
		print(r, *p)


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("gff_file", type=str)
	args = ap.parse_args()

	GffIndexer(args.gff_file)#.get_index()

if __name__ == "__main__":
	main()
