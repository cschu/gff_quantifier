from dataclasses import dataclass, asdict


@dataclass(slots=True)
class ReferenceHit:
    rid: int = None
    start: int = None
    end: int = None
    rev_strand: bool = None
    cov_start: int = None
    cov_end: int = None
    has_annotation: bool = None
    n_aln: int = None
    is_ambiguous: bool = None
    library_mod: int = None
    mate_id: int = None

    def __hash__(self):
        return hash(tuple(asdict(self).values()))

    def __eq__(self, other):
        return all(
            item[0][1] == item[1][1]
            for item in zip(
                sorted(asdict(self).items()),
                sorted(asdict(other).items())
            )
        )

    def __str__(self):
        return "\t".join(map(str, asdict(self).values()))

    def __repr__(self):
        return str(self)
