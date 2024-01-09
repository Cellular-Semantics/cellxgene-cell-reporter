import os
import cellxgene_census

import pandas as pd


TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../templates/cellxgene_subset.tsv")


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


if __name__ == "__main__":
    cl_terms = set(report_cl_coverage("homo_sapiens"))
    cl_terms.update(set(report_cl_coverage("mus_musculus")))
    print("Total {} unique CL terms mentioned in CxG census.".format(len(cl_terms)))

    generate_robot_template(cl_terms)
