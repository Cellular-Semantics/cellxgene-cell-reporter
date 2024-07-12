import os
import cellxgene_census

import pandas as pd


TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../templates/cellxgene_subset.tsv")
README_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../README.md")


def report_cl_coverage(species):
    """
    Queries cellxgene census to report mentioned CL terms.
    :param species: species to report. CxG census only covers human and mouse for now.
    :return: unique list of mentioned CL terms.
    """
    with cellxgene_census.open_soma() as census:
        # Reads SOMADataFrame as a slice
        cell_metadata = census["census_data"][species].obs.read(
            column_names=["cell_type_ontology_term_id"]
        )

        # Concatenates results to pyarrow.Table
        cell_metadata = cell_metadata.concat()
        cell_metadata = cell_metadata.to_pandas()
        # print(cell_metadata)

        mentioned_terms = cell_metadata.cell_type_ontology_term_id.unique().tolist()
        if "unknown" in mentioned_terms:
            mentioned_terms.remove("unknown")
        print("{} unique CL terms mentioned at {} datasets.".format(len(mentioned_terms), species))
        # print(mentioned_terms[:10])

        return mentioned_terms


def generate_robot_template(mentioned_terms):
    robot_template_seed = {'ID': 'ID',
                           'subset': 'AI oboInOwl:inSubset',
                           }
    dl = [robot_template_seed]
    for cl_term in mentioned_terms:
        d = dict()
        d["ID"] = cl_term
        d["subset"] = "http://purl.obolibrary.org/obo/cl#cellxgene_subset"
        dl.append(d)
    robot_template = pd.DataFrame.from_records(dl)
    robot_template.to_csv(TEMPLATE_PATH, sep="\t", index=False)


def update_read_me(total_count, counts_dict):
    readme = "# CellxGene Cell Ontology Coverage \n"
    readme += "There are **{}** unique Cell Ontology terms referenced within CellxGene.\n\n".format(total_count)
    readme += "| Species | CL terms |\n"
    readme += "|---------|----------|\n"
    for species in counts_dict:
        readme += "| {} | {} |\n".format(species, counts_dict[species])
    readme += "\n"

    with open(README_PATH, "w") as readme_file:
        readme_file.write(readme)


if __name__ == "__main__":
    hs_cl_terms = set(report_cl_coverage("homo_sapiens"))
    mm_cl_terms = set(report_cl_coverage("mus_musculus"))
    all_cl_terms = hs_cl_terms.union(mm_cl_terms)

    print("Total {} unique CL terms mentioned in CxG census.".format(len(all_cl_terms)))
    generate_robot_template(all_cl_terms)
    update_read_me(len(all_cl_terms), {"human": len(hs_cl_terms), "mouse": len(mm_cl_terms)})
