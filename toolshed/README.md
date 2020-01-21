Personal orthology db based on data from OrthoDB v9 (https://www.orthodb.org/).
See https://www.orthodb.org/v9/download/README.txt for file and table info

# Create DB

```sh
# odb9_levels.tab:
sqlite3 ortho.db <<EOF
CREATE TABLE levels(
  NCBI_tax_id INTEGER PRIMARY KEY,
  scientific_name TEXT,
  total_non_redundant_count_of_genes_in_all_underneath_clustered_species INTEGER,
  total_count_of_OGs_built_on_it INTEGER,
  total_non_redundant_count_of_species_underneath INTEGER
);
EOF

# odb9_species.tab
sqlite3 ortho.db <<EOF
CREATE TABLE species(
  NCBI_tax_id INTEGER,
  scientific_name TEXT,
  total_count_of_clustered_genes_in_this_species INTEGER,
  total_count_of_the_OGs_it_participates INTEGER,
  mapping_type_Clustered_Mapped TEXT
);
EOF

# odb9_genes.tab
sqlite3 ortho.db <<EOF
CREATE TABLE genes(
  Ortho_DB_gene_id TEXT PRIMARY KEY,
  organism_tax_id INTEGER,
  protein_original_sequence_id TEXT,
  uniprot_id TEXT,
  ensembl_id TEXT,
  ncbi_id TEXT,
  description TEXT
);
CREATE INDEX genes__uniprot_id ON genes (uniprot_id);
CREATE INDEX genes__ensembl_id ON genes (ensembl_id);
CREATE INDEX genes__ncbi_id ON genes (ncbi_id);
EOF

# odb9_OGs.tab
sqlite3 ortho.db <<EOF
CREATE TABLE OGs(
  OG_unique_id TEXT PRIMARY KEY,
  level_tax_id_on_which_the_group_was_built INTEGER,
  OG_name TEXT
);
CREATE INDEX OGs__level_tax_id_on_which_the_group_was_built ON OGs (level_tax_id_on_which_the_group_was_built);
EOF

# odb9_OG2genes.tab
sqlite3 ortho.db <<EOF
CREATE TABLE OG2genes(
  OG_unique_id TEXT,
  Ortho_DB_gene_id TEXT,
  PRIMARY KEY (OG_unique_id, Ortho_DB_gene_id)
);
CREATE INDEX OG2genes__Ortho_DB_gene_id ON OG2genes (Ortho_DB_gene_id);
EOF

# odb9_fasta_fungi.tgz
sqlite3 ortho.db <<EOF
CREATE TABLE fasta_fungi(
  Ortho_DB_gene_id_with_NCBI_tax_id TEXT PRIMARY KEY,
  uniprot_id TEXT,
  seq TEXT
);
CREATE INDEX fasta_fungi__uniprot_id ON fasta_fungi (uniprot_id);
EOF
```

# Download data and insert on the fly to the db

Get it from [OrthoDB v9](https://www.orthodb.org/v9/download/)

```sh
cat <<EOF > sqlite.init
.mode tabs
.separator "\t"
EOF
wget -qO- https://www.orthodb.org/v9/download/odb9_levels.tab.gz | gzip -cd | sqlite3 -init sqlite.init ortho.db ".import /dev/stdin levels"
wget -qO- https://www.orthodb.org/v9/download/odb9_species.tab.gz | gzip -cd | sqlite3 -init sqlite.init ortho.db ".import /dev/stdin species"
wget -qO- https://www.orthodb.org/v9/download/odb9_genes.tab.gz | gzip -cd | sqlite3 -init sqlite.init ortho.db ".import /dev/stdin genes"
wget -qO- https://www.orthodb.org/v9/download/odb9_OGs.tab.gz | gzip -cd | sqlite3 -init sqlite.init ortho.db ".import /dev/stdin OGs"
wget -qO- https://www.orthodb.org/v9/download/odb9_OG2genes.tab.gz | gzip -cd | sqlite3 -init sqlite.init ortho.db ".import /dev/stdin OG2genes"
wget -qO- https://www.orthodb.org/v9/download/odb9_fasta_fungi.tgz | tar xz --to-stdout | paste -s -d '\t\n' | sed 's/^>//' | sed 's/\s/\t/' | sqlite3 -init sqlite.init ortho.db ".import /dev/stdin fasta_fungi"
rm sqlite.init
```

# An example query

Get protein sequences of the orthotologs of 2 proteins in some species of Saccharomyces

```sql
-- NOTES:
--   * Taxonomy level 4893: Saccharomycetaceae
--   * Species ids:
--     - 1064592: Naumovozyma castellii
--     - 1160507: Saccharomyces arboricola
--     - 226230: Saccharomyces kudriavzevii
--     - 559292: Saccharomyces cerevisiae S288c

-- Get sequences of orthologous genes for the 4 species up here
--EXPLAIN QUERY PLAN  -- check that all steps scan through indexes
SELECT a.organism_tax_id, a.uniprot_id, a.description, b.OG_unique_id, c.seq
  FROM genes AS a
  JOIN OG2genes AS b ON a.Ortho_DB_gene_id=b.Ortho_DB_gene_id
  JOIN fasta_fungi AS c ON a.uniprot_id=c.uniprot_id
 WHERE b.OG_unique_id IN (
   -- Get orthologous groups of our gene of interest, calculated at the Saccharomycetaceae level
   SELECT b.OG_unique_id
     FROM genes AS a
     JOIN OG2genes AS b ON a.Ortho_DB_gene_id=b.Ortho_DB_gene_id
     JOIN OGs AS c ON b.OG_unique_id=c.OG_unique_id
    WHERE a.uniprot_id IN (
      "Q07807",  --PUF3
      "P04147"   --PAB1
    )                       
      AND c.level_tax_id_on_which_the_group_was_built=4893
 )
   AND a.organism_tax_id IN (1064592, 226230, 559292, 1160507)
 ORDER BY b.OG_unique_id, a.organism_tax_id, a.uniprot_id;
```
