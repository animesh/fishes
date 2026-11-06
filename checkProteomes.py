import requests
import pandas as pd
import time
import json
from pathlib import Path

SPECIES = [
    "Gadus morhua",
    "Merluccius merluccius",
    "Pollachius virens",
    "Melanogrammus aeglefinus",
    "Pollachius pollachius",
    "Molva molva",
    "Trisopterus minutus",
    "Trisopterus esmarkii",
    "Micromesistius poutassou",
    "Gadiculus argenteus",
    "Clupea harengus",
    "Scomber scombrus",
    "Trachurus trachurus",
    "Microstomus kitt",
    "Pleuronectes platessa",
    "Limanda limanda",
    "Lepidorhombus whiffiagonis",
    "Scophthalmus maximus",
    "Cyclopterus lumpus",
    "Maurolicus muelleri",
    "Ammodytes marinus",
    "Anarhichas lupus",
    "Argentina sphyraena",
    "Squalus acanthias",
    "Amblyraja radiata",
    "Callionymus lyra",
    "Argentina silus",
    "Merlangius merlangus",
    "Eutrigla gurnardus",
]

session = requests.Session()

Path("uniprot_json").mkdir(exist_ok=True)

def get_taxid(species):

    url = (
        "https://rest.uniprot.org/taxonomy/search"
        f'?query="{species}"&format=json&size=50'
    )

    r = session.get(url, timeout=120)
    r.raise_for_status()

    js = r.json()

    results = js.get("results", [])

    if not results:
        return None

    species_lower = species.lower()

    for x in results:

        sci = x.get("scientificName", "")

        if sci.lower() == species_lower:
            return x["taxonId"]

    return results[0]["taxonId"]


def find_species_taxon(record):

    for node in record.get("taxonLineage", []):

        if node.get("rank") == "species":
            return node.get("taxonId")

    return None


def rank_proteome(record):

    score = 0

    ptype = str(record.get("proteomeType", "")).lower()

    if "reference" in ptype:
        score += 100000

    if "representative" in ptype:
        score += 50000

    score += int(record.get("proteinCount", 0))
    score += int(record.get("annotationScore", 0)) * 1000

    return score


all_rows = []

for species in SPECIES:

    print(f"\n{species}")

    taxid = get_taxid(species)

    if taxid is None:

        all_rows.append(
            {
                "input_species": species,
                "status": "TAXON_NOT_FOUND",
            }
        )

        continue

    url = (
        "https://rest.uniprot.org/proteomes/search"
        f"?query=taxonomy_id:{taxid}&format=json&size=500"
    )

    r = session.get(url, timeout=120)
    r.raise_for_status()

    results = r.json().get("results", [])

    with open(
        Path("uniprot_json")
        / f"{species.replace(' ','_')}.json",
        "w",
        encoding="utf-8"
    ) as fh:
        json.dump(results, fh, indent=2)

    if not results:

        all_rows.append(
            {
                "input_species": species,
                "query_taxid": taxid,
                "status": "NO_PROTEOME_FOUND",
            }
        )

        continue

    for p in results:

        taxonomy = p.get("taxonomy", {})

        assembly = p.get("genomeAssembly", {})

        pid = p.get("id")

        all_rows.append(
            {
                "input_species": species,
                "query_taxid": taxid,

                "proteome_id": pid,

                "taxonomy_name":
                    taxonomy.get("scientificName"),

                "taxonomy_taxid":
                    taxonomy.get("taxonId"),

                "species_taxid":
                    find_species_taxon(p),

                "proteome_type":
                    p.get("proteomeType"),

                "protein_count":
                    p.get("proteinCount"),

                "gene_count":
                    p.get("geneCount"),

                "annotation_score":
                    p.get("annotationScore"),

                "assembly_id":
                    assembly.get("assemblyId"),

                "assembly_source":
                    assembly.get("source"),

                "assembly_level":
                    assembly.get("level"),

                "assembly_url":
                    assembly.get("genomeAssemblyUrl"),

                "reviewed_proteins":
                    p.get(
                        "proteomeStatistics",
                        {}
                    ).get(
                        "reviewedProteinCount"
                    ),

                "unreviewed_proteins":
                    p.get(
                        "proteomeStatistics",
                        {}
                    ).get(
                        "unreviewedProteinCount"
                    ),

                "rank":
                    rank_proteome(p),

                "fasta_url":
                    f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=false&query=proteome:{pid}"
                    if pid else None,

                "proteome_url":
                    f"https://www.uniprot.org/proteomes/{pid}"
                    if pid else None,

                "status":
                    "FOUND",
            }
        )

    time.sleep(0.25)

all_df = pd.DataFrame(all_rows)

all_df.to_csv(
    "all_proteomes.tsv",
    sep="\t",
    index=False,
)

all_df.to_excel(
    "all_proteomes.xlsx",
    index=False,
)

best_rows = []

for species in SPECIES:

    tmp = all_df[
        (all_df["input_species"] == species)
        &
        (all_df["status"] == "FOUND")
    ]

    if len(tmp) == 0:
        continue

    tmp = tmp.sort_values(
        "rank",
        ascending=False
    )

    best_rows.append(
        tmp.iloc[0].to_dict()
    )

best_df = pd.DataFrame(best_rows)

best_df.to_csv(
    "best_proteomes.tsv",
    sep="\t",
    index=False,
)

best_df.to_excel(
    "best_proteomes.xlsx",
    index=False,
)

summary = (
    all_df
    .groupby("input_species")
    .size()
    .reset_index(name="n_proteomes")
)

summary.to_csv(
    "summary.tsv",
    sep="\t",
    index=False,
)

print("\nSaved:")
print("  all_proteomes.tsv")
print("  all_proteomes.xlsx")
print("  best_proteomes.tsv")
print("  best_proteomes.xlsx")
print("  summary.tsv")
print("  uniprot_json/")