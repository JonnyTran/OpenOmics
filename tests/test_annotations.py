from openomics.database import GENCODE, RNAcentral, MirBase, GTEx, GeneOntology
from .test_multiomics import *


@pytest.fixture
def generate_GENCODE_ftp():
    return GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                   file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                   "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                   "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                   agg_sequences="shortest")


def test_import_gencode_db(generate_GENCODE_ftp):
    assert generate_GENCODE_ftp.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/'


def test_annotate_gencode(generate_TCGA_LUAD, generate_GENCODE_ftp):
    generate_TCGA_LUAD.LncRNA.annotate_genomics(generate_GENCODE_ftp, index="gene_id",
                                                columns=['feature', 'start', 'end', 'strand', 'tag', 'havana_gene'])
    assert {'feature', 'start', 'end', 'strand', 'tag', 'havana_gene'}.issubset(
        generate_TCGA_LUAD.LncRNA.get_annotations().columns)


@pytest.fixture
def generate_RNACentral_ftp():
    return RNAcentral(path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/",
                      file_resources={
                          "rnacentral_rfam_annotations.tsv": "go_annotations/rnacentral_rfam_annotations.tsv.gz",
                          "gencode.tsv": "id_mapping/database_mappings/gencode.tsv"},
                      )


def test_import_rnacentral_db(generate_RNACentral_ftp):
    assert generate_RNACentral_ftp.data_path == 'ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/'


def test_rnacentral_annotate(generate_TCGA_LUAD, generate_RNACentral_ftp):
    generate_TCGA_LUAD.MessengerRNA.annotate_genomics(database=generate_RNACentral_ftp, index="gene_name",
                                                      columns=['gene_name', 'transcript_id', 'RNA type', 'go_id',
                                                               'Rfams'])
    assert {'transcript_id', 'RNA type', 'go_id', 'Rfams'}.issubset(
        generate_TCGA_LUAD.MessengerRNA.get_annotations().columns)


@pytest.fixture
def generate_MirBase_ftp():
    return MirBase(path="ftp://mirbase.org/pub/mirbase/CURRENT/")


def test_import_mirbase_db(generate_MirBase_ftp):
    assert generate_MirBase_ftp.data_path == "ftp://mirbase.org/pub/mirbase/CURRENT/"


@pytest.fixture
def generate_GTEx_expressions():
    return GTEx(path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/", )


def test_import_GTEx(generate_GTEx_expressions):
    assert generate_GTEx_expressions.data_path == "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/"
    assert not generate_GTEx_expressions.get_expressions(index="gene_name").empty
    assert not generate_GTEx_expressions.get_expressions(index="gene_id").empty


def test_GTEx_annotate(generate_TCGA_LUAD, generate_GTEx_expressions):
    generate_TCGA_LUAD.MessengerRNA.annotate_expressions(database=generate_GTEx_expressions, index="gene_name")
    generate_TCGA_LUAD.LncRNA.annotate_expressions(database=generate_GTEx_expressions, index="gene_id")
    generate_TCGA_LUAD.MicroRNA.annotate_expressions(database=generate_GTEx_expressions, index="gene_name")
    assert not generate_TCGA_LUAD.MessengerRNA.get_annotation_expressions().empty
    assert not generate_TCGA_LUAD.LncRNA.get_annotation_expressions().empty
    assert not generate_TCGA_LUAD.MicroRNA.get_annotation_expressions().empty


@pytest.fixture
def generate_GeneOntology():
    return GeneOntology(path="http://geneontology.org/gene-associations/",
                        file_resources={"goa_human.gaf": "goa_human.gaf.gz",
                                        "goa_human_rna.gaf": "goa_human_rna.gaf.gz"}
                        )


def test_import_GeneOntology(generate_GeneOntology):
    assert generate_GeneOntology.data_path == "http://geneontology.org/gene-associations/"


def test_annotate_GeneOntology(generate_TCGA_LUAD, generate_GeneOntology):
    generate_TCGA_LUAD.MessengerRNA.annotate_genomics(database=generate_GeneOntology, index="gene_name",
                                                      columns=['go_id'])
    generate_TCGA_LUAD.LncRNA.annotate_genomics(database=generate_GeneOntology, index="gene_id",
                                                columns=['go_id'])
    generate_TCGA_LUAD.MicroRNA.annotate_genomics(database=generate_GeneOntology, index="gene_name",
                                                  columns=['go_id'])
    assert {'go_id'}.issubset(generate_TCGA_LUAD.MessengerRNA.get_annotations().columns)
    assert {'go_id'}.issubset(generate_TCGA_LUAD.LncRNA.get_annotations().columns)
    assert {'go_id'}.issubset(generate_TCGA_LUAD.MicroRNA.get_annotations().columns)
    assert not generate_TCGA_LUAD.LncRNA.annotations["go_id"].empty
    assert not generate_TCGA_LUAD.MessengerRNA.annotations["go_id"].empty
    assert not generate_TCGA_LUAD.MicroRNA.annotations["go_id"].empty
