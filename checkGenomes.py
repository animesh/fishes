import pandas as pd
import requests

# -------------------------------
# Load previous UniProt audit
# -------------------------------

proteomes = pd.read_csv(
    "all_proteomes.tsv",
    sep="\t"
)

species_list = (
    proteomes.loc[
        proteomes["status"] == "NO_PROTEOME_FOUND",
        "input_species"
    ]
    .dropna()
    .unique()
)

rows = []

for species in species_list:

    print(species)

    try:

        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/"
            f"{species}/dataset_report"
        )

        r = requests.get(
            url,
            timeout=120
        )

        r.raise_for_status()

        data = r.json()

        reports = data.get("reports", [])

        if not reports:

            rows.append(
                {
                    "species": species,
                    "assembly_found": False,
                }
            )
            continue

        # Prefer RefSeq assembly
        best = None

        for rep in reports:

            acc = rep.get("accession", "")

            if acc.startswith("GCF_"):
                best = rep
                break

        if best is None:
            best = reports[0]

        assembly = best.get("accession")

        annotation = (
            best.get("annotation_info")
            or {}
        )

        gene_counts = (
            annotation.get("stats", {})
            .get("gene_counts", {})
        )

        protein_coding = gene_counts.get(
            "protein_coding"
        )

        annotated = bool(annotation)

        # NCBI accession page

        assembly_page = (
            "https://www.ncbi.nlm.nih.gov/datasets/genome/"
            f"{assembly}/"
            if assembly else None
        )

        #
        # Datasets package URLs
        #
        # NCBI packaging changes occasionally, but
        # these URLs are the easiest to test.
        #

        genome_fasta_url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{assembly}/download?include_annotation_type=GENOME_FASTA"
            if assembly else None
        )

        protein_fasta_url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{assembly}/download?include_annotation_type=PROTEIN_FASTA"
            if assembly else None
        )

        gff_url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{assembly}/download?include_annotation_type=GFF3"
            if assembly else None
        )

        gbff_url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{assembly}/download?include_annotation_type=GENOME_GBFF"
            if assembly else None
        )

        rows.append(
            {
                "species": species,
                "assembly_found": True,
                "assembly": assembly,
                "annotated": annotated,
                "protein_coding_genes": protein_coding,
                "assembly_page": assembly_page,
                "genome_fasta_url": genome_fasta_url,
                "protein_fasta_url": protein_fasta_url,
                "gff3_url": gff_url,
                "gbff_url": gbff_url,
            }
        )

    except Exception as e:

        rows.append(
            {
                "species": species,
                "assembly_found": False,
                "error": str(e),
            }
        )

report = pd.DataFrame(rows)

report.to_excel(
    "fish_genome_proteome_links.xlsx",
    index=False
)

report.to_csv(
    "fish_genome_proteome_links.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("  fish_genome_proteome_links.xlsx")
print("  fish_genome_proteome_links.tsv")