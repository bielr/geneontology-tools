import csv
from gzip import GzipFile
import networkx as nx
from os import path

from geneontology import get_curated_evidence_codes


class GAF(object):
    def __init__(self, stream):
        self._stream = stream

    @property
    def stream(self):
        return _stream

    def stream_simplified(self):
        for gaf_row in self._stream:
            db_obj_id = gaf_row[1]
            qualifiers = gaf_row[3].split('|')
            go_id = gaf_row[4]
            evidence_code = gaf_row[6]
            species_ncbi_ids = [int(taxon.replace('taxon:', '')) for taxon in gaf_row[12].split('|')]

            yield (species_ncbi_ids, db_obj_id, qualifiers, go_id, evidence_code)

    def uniprot_proteins_only(self):
        def stream():
            for gaf_row in self._stream:
                db = gaf_row[0]
                obj_type = gaf_row[11]

                if db.startswith('UniProt') and obj_type == 'protein':
                    yield gaf_row

        return GAF(stream())

    def positively_qualified_only(self):
        def stream():
            for gaf_row in self._stream:
                qualifiers = gaf_row[3].split('|')

                if 'NOT' not in qualifiers:
                    yield gaf_row

        return GAF(stream())

    def for_evidence_codes_only(self, evidence_codes):
        def stream():
            for gaf_row in self._stream:
                evidence_code = gaf_row[6]

                if evidence_code in evidence_codes:
                    yield gaf_row

        return GAF(stream())

    def for_object_ids_only(self, db_obj_ids):
        db_obj_ids = frozenset(db_obj_ids)

        def stream():
            for gaf_row in self._stream:
                db_obj_id = gaf_row[1]

                if db_obj_id in db_obj_ids:
                    yield gaf_row

        return GAF(stream())

    def for_species_only(self, species_ncbi):
        req_taxons = [f'taxon:{species_id}' for species_id in species_ncbi]

        def stream():
            for gaf_row in self._stream:
                taxons = gaf_row[12].split('|')

                if any(req_taxon in taxons for req_taxon in req_taxons):
                    yield gaf_row

        return GAF(stream())

    def for_transitively_annotated(self, rel_g, go_ids):
        go_ids = frozenset(go_ids)

        def stream():
            for gaf_row in self._stream:
                go_id = gaf_row[4]

                if go_id in go_ids:
                    yield gaf_row

                elif go_id in rel_g and (go_ids & nx.descendants(rel_g, go_id)):
                    yield gaf_row

        return GAF(stream())


def get_uniprot_gaf_path():
    script_dir = path.dirname(path.realpath(__file__))
    geneontology_dir = path.dirname(script_dir)

    return path.join(geneontology_dir, 'extra', 'goa_uniprot_all.gaf.gz')


def open_gaf(local=True):
    if local:
        return open(get_uniprot_gaf_path(), 'rb')
    else:
        import urllib
        return urllib.urlopen('ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz', 'rb')

def load_gaf_entries(local=True):
    with open_gaf(local=local) as zf:
        with GzipFile(fileobj=zf) as f:
            f_without_comments = (line.decode('utf-8') for line in f if not line.startswith(b'!'))
            yield from csv.reader(f_without_comments, delimiter='\t')


def get_gaf_protein_annotations(gaf_rows, rel_g, uniprot_ids, species_ncbi=None, evidence_codes=get_curated_evidence_codes()):
    gaf = GAF(gaf_rows)                          \
        .for_object_ids_only(uniprot_ids)        \
        .for_evidence_codes_only(evidence_codes) \
        .uniprot_proteins_only()                 \
        .positively_qualified_only()

    if species_ncbi is not None:
        gaf = gaf.for_species_only(species_ncbi)

    return gaf

def get_gaf_transitively_annotatated_uniprot_proteins(gaf_rows, rel_g, go_ids, evidence_codes=get_curated_evidence_codes()):
    gaf = GAF(gaf_rows)                             \
        .for_evidence_codes_only(evidence_codes)    \
        .uniprot_proteins_only()                    \
        .positively_qualified_only()                \
        .for_transitively_annotated(rel_g, go_ids)

    return gaf

