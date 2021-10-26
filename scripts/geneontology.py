def connect_through_docker_network():
    import MySQLdb
    return MySQLdb.connect(host='geneontology', port=3306, user='geneontology', password='geneontology', database='geneontology')


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


def get_explicit_uniprot_annotations(cursor, uniprot_ids, species_ncbi=None, evidence_codes=get_curated_evidence_codes()):
    sql = """
        select distinct
          species.ncbi_taxa_id ncbi_taxa_id,
          dbxref.xref_key      product_key,
          association.is_not   is_not,
          term.acc             term_acc,
          term.term_type       term_type,
          evidence.code        evidence_code
        from
          association
          inner join gene_product on (gene_product.id         = association.gene_product_id)
          inner join dbxref       on (dbxref.id               = gene_product.dbxref_id)
          inner join evidence     on (evidence.association_id = association.id)
          inner join term         on (term.id                 = association.term_id)
          inner join species     on (species.id              = gene_product.species_id)
        where
          dbxref.xref_dbname like '%%uniprot%%'
          and
          dbxref.xref_key in %(uniprot_ids)s
          and
          evidence.code in %(evidence_codes)s
          and
          (not term.is_obsolete)
      """

    if species_ncbi is None:
        sql += ";"
    else:
        sql += """
          and
          species.ncbi_taxa_id = %(species_ncbi)s;
          """
    cursor.execute(sql,
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


