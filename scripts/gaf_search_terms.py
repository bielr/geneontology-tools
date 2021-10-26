import csv
import sys

import geneontology as godb
import gaf_util as gaf


if __name__ == '__main__':
    uniprot_acs = sys.argv[1:]
    assert len(uniprot_acs) > 0

    go_onto = godb.load_go_obo()
    go_is_a_g = godb.onto_rel_graph(go_onto)
    gaf_rows = gaf.load_gaf_entries()

    selected = gaf.get_gaf_protein_annotations(gaf_rows, go_is_a_g, uniprot_acs)

    out_csv = csv.writer(sys.stdout, delimiter='\t')

    for species_ncbi_ids, db_obj_id, qualifiers, go_id, evidence_code in selected.stream_simplified():
        out_csv.writerow(['|'.join(map(str, species_ncbi_ids)), db_obj_id, go_id, evidence_code])
