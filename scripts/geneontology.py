import pronto
import networkx as nx


def get_obo_path():
    return '/data/pinaweb/geneontology/go.obo'

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


def connect_to_docker():
    import MySQLdb
    import docker

    docker_client = docker.from_env()

    containers = docker_client.containers.list(
        filters = {'label': [
            'com.docker.compose.project=geneontology',
            'com.docker.compose.service=geneontology'
        ]})

    assert len(containers) == 1

    container = containers[0]
    network_settings = container.attrs['NetworkSettings']

    host = network_settings['Networks']['geneontology-net']['IPAddress']
    env = dict(e.split('=', 1) for e in container.attrs['Config']['Env'])

    return MySQLdb.connect(host=host, port=3306, user=env['MYSQL_USER'], password=env['MYSQL_PASSWORD'], database=env['MYSQL_DATABASE'])

def get_uniprot_gene_products(cursor, uniprot_ids):
    cursor.execute("""
        select
          gene_product.id     product_id,
          dbxref.xref_dbname  xref_dbname,
          dbxref.xref_key     xref_key
        from
          gene_product
          inner join dbxref on (dbxref.id = gene_product.dbxref_id)
        where
          dbxref.xref_dbname like '%%uniprot%%'
          and
          dbxref.xref_key in %(uniprot_ids)s;
        """,
        {'uniprot_ids': tuple(uniprot_ids)})

    return cursor.fetchall()

def get_explicit_uniprot_annotations(cursor, uniprot_ids, species_ncbi, evidence_codes=get_curated_evidence_codes()):
    cursor.execute("""
        select distinct
          dbxref.xref_key     product_key,
          association.is_not  is_not,
          term.acc            term_acc,
          term.term_type      term_type,
          evidence.code       evidence_code
        from
          association
          inner join gene_product on (gene_product.id         = association.gene_product_id)
          inner join dbxref       on (dbxref.id               = gene_product.dbxref_id)
          inner join evidence     on (evidence.association_id = association.id)
          inner join term         on (term.id                 = association.term_id)
        where
          gene_product.species_id = (select id from species where ncbi_taxa_id = %(species_ncbi)s)
          and
          dbxref.xref_dbname like '%%uniprot%%'
          and
          dbxref.xref_key in %(uniprot_ids)s
          and
          evidence.code in %(evidence_codes)s
          and
          (not term.is_obsolete);
        """,
        {'species_ncbi': species_ncbi, 'uniprot_ids': tuple(uniprot_ids), 'evidence_codes': tuple(evidence_codes)})

    return cursor.fetchall()

def get_transitive_uniprot_annotations(cursor, uniprot_ids, species_ncbi, evidence_codes=get_curated_evidence_codes()):
    cursor.execute("""
        select distinct
          dbxref.xref_key     product_key,
          association.is_not  is_not,
          term2.acc           term_acc,
          term2.term_type     term_type,
          evidence.code       evidence_code
        from
          association
          inner join gene_product   on (gene_product.id         = association.gene_product_id)
          inner join dbxref         on (dbxref.id               = gene_product.dbxref_id)
          inner join evidence       on (evidence.association_id = association.id)
          left join graph_path      on (graph_path.term1_id     = association.term_id)
          inner join term term2     on (term2.id                = graph_path.term2_id)
        where
          gene_product.species_id = (select id from species where ncbi_taxa_id = %(species_ncbi)s)
          and
          graph_path.relationship_type_id = (select id from term where acc = 'is_a')
          and
          dbxref.xref_dbname like '%%uniprot%%'
          and
          dbxref.xref_key in %(uniprot_ids)s
          and
          evidence.code in %(evidence_codes)s
          and
          (not term2.is_obsolete);
        """,
        {'species_ncbi': species_ncbi, 'uniprot_ids': tuple(uniprot_ids), 'evidence_codes': tuple(evidence_codes)})

    return cursor.fetchall()

def count_protein_annotations(cursor, term_acc, evidence_codes=get_curated_evidence_codes()):
    cursor.execute("""
        select
          count(distinct gene_product.id)
        from
          association
          inner join graph_path   on (graph_path.term1_id     = association.term_id)
          inner join gene_product on (gene_product.id         = association.gene_product_id)
          inner join evidence     on (evidence.association_id = association.id)
        where
          gene_product.type_id = (select id from term where acc = 'protein')
          and
          term2_id = (select id from term where acc = %(term2_acc)s)
          and
          evidence.code in %(evidence_codes)s;
        """,
        {'term2_acc': term_acc, 'evidence_codes': evidence_codes})

    #inner join dbxref       on (dbxref.id               = gene_product.dbxref_id)
    #inner join term term2   on (term2.id                = graph_path.term2_id)
    #dbxref.xref_dbname like '%%uniprot%%' and

    cnt, = cursor.fetchall()
    cnt, = cnt
    return cnt
