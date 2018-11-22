!#/usr/bin/env python3

from intermine.webservice import Service

cycle_genes_file = snakemake.output['cycle_genes']

service = Service("https://apps.araport.org:443/thalemine/service")
query = service.new_query("Gene")
query.add_view("primaryIdentifier")
query.add_constraint("curatorSummary", "CONTAINS", "cell cycle", code="A")

cell_cycle_genes = [x["primaryIdentifier"] for x in query.rows()]

with open(cycle_genes_file, 'wt') as f:
    f.write('\n'.join(cell_cycle_genes))
