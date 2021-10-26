import networkx as nx
import pronto
from os import path
import sys


def get_obo_path():
    script_dir = path.dirname(path.realpath(__file__))
    geneontology_dir = path.dirname(script_dir)

    # '/data/pinaweb/geneontology/extra/go.obo'
    return path.join(geneontology_dir, 'extra', 'go.obo')

def get_curated_evidence_codes():
    return ('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'IC')

def load_onto(obo_path):
    return pronto.Ontology(obo_path)

def load_go_obo():
    return load_onto(get_obo_path())

def is_obsolete(term, onto=None):
    if isinstance(term, str):
        if term not in onto.terms:
            return True
        term = onto.terms[term]

    return 'is_obsolete' in term.other and 'true' in term.other['is_obsolete']

def is_disconnected(term, onto):
    if isinstance(term, str):
        if term not in onto.terms:
            return True
        term = onto[term]

    if term.relations:
        return False

    return not any(term in term2_rels for term2      in onto
                                      for term2_rels in term2.relations.values())

def onto_rel_graph(onto, rel_type = pronto.Relationship('is_a')):
    obo_is_a = nx.DiGraph()

    nodes = {u for u,term in onto.terms.items() if not is_obsolete(term)}

    obo_is_a.add_nodes_from(nodes)

    rels = [(u,vs.id) for u in obo_is_a.nodes()
                      for rel, vs in onto.terms[u].relations.items()
                      if rel == rel_type]

    edges = ((u,v) for u, vs in rels for v in vs if v in nodes)

    obo_is_a.add_edges_from(edges)

    return obo_is_a

def term_alt_ids(term):
    return term.other['alt_id'] if 'alt_id' in term.other else []

def onto_alt_id_graph(onto, onto_rel_graph):
    g = nx.Graph()

    g.add_nodes_from(onto.terms.keys())
    g.add_edges_from((term.id, alt_id) for term in onto for alt_id in term_alt_ids(term))

    return g

def find_alternatives(go, alt_id_g):
    return list(alt_id for alt_id in nx.node_connected_component(alt_id_g, go))

def find_valid_alternatives(go, alt_id_g, rel_g):
    r = [alt_id for alt_id in find_alternatives(go, alt_id_g)
                if alt_id in rel_g and not nx.is_isolate(rel_g, alt_id)]

    if len(r) > 1:
        print(f'found {len(r)} valid alternatives of {go}', file=sys.stderr)

    return r


