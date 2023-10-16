# pylint: disable=C0103,R0902,R0913,W2301,W1203

""" module docstring """

import logging

from .database_importer import GqDatabaseImporter
from ..models import db


logger = logging.getLogger(__name__)

SPIRE_ANNOTATION_COLUMNS = ",".join(
    (
        "COG_category",
        "GOs",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "PFAMs",
    )
)

class SmallGenomeDatabaseImporter(GqDatabaseImporter):
    def __init__(
        self,
        db_path=None,
        db_session=None,
        columns=SPIRE_ANNOTATION_COLUMNS,
    ):
        self.columns = set(columns.split(","))

        super().__init__(db_path=db_path, db_session=db_session)

    def parse_annotations(self, input_data, input_data2=None):

        genes_in = (
            line.split(";")[0]
            for line in input_data.read().decode().strip().split("\n")
            if line[0] != "#"
        )

        seq_features = {}            
        for self.nseqs, line in enumerate(genes_in, start=1):
            # k141_36810	Prodigal_v2.6.3	CDS	211	366	5.7	-	0	ID=1_1
            contig, _, _, start, end, _, strand, _, seq_id = line.split("\t")
            gene_id = "_".join((contig, seq_id[seq_id.rfind("_") + 1:]))
            seq_features[gene_id] = db.AnnotatedSequence(                
                seqid=contig,
                featureid=gene_id,
                start=int(start),
                end=int(end),
                strand=int(strand == "+") if strand is not None else None,
            )
        
        annotations_in = (
            line 
            for line in input_data2.read().decode().strip().split("\n")
            if line[:2] != "##"
        )

        header_line = next(annotations_in, None)
        if header_line is None or header_line[0] != "#":
            raise ValueError("Missing header line.")
        
        header_line = header_line[1:].split("\t")        

        for self.nseqs, line in enumerate(annotations_in, start=1):
            if line and line[0] != "#":
                line = line.strip().split("\t")

                seq_feature = seq_features.get(line[0])
                if seq_feature is None:
                    raise ValueError(f"Found entry without matching gene {line}.")

                columns = {
                    colname: coldata.strip()
                    for colname, coldata in zip(header_line, line)
                    if colname in self.columns
                }

                annotation = [
                    (category, tuple(features.split(",")))
                    for category, features in columns.items()
                    if features and features != "-"
                ]

                if annotation:                    
                    yield seq_feature, annotation
