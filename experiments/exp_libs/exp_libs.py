from Bio.Restriction import BsaI
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
import re
from os import path
import os
import copy
from QUEEN.qobj import QUEEN
from QUEEN.queen import cropdna
from sangerseq_viewer import recursive_alignment, abi_to_dict, generate_consensusseq
from typing import List, Optional, Tuple, Dict
import json


DNA_STUFFER = 'A' * 5


class SequencingResSetting:
    def __init__(self, primer: str, reference: str):
        self.primer = primer
        self.reference = reference
        self.ANALYSE_REGION_DICT = {
            ('ALL482', 'AzR_T7_assembly'): {'analyse_region': (3450, 4300), 'coding_region': (3568, 4248)},
            ('ALL486', 'AzR_T7_assembly'): {'analyse_region': (3450, 4300), 'coding_region': (3568, 4248)},
        }
        self.start = self.ANALYSE_REGION_DICT[(self.primer, self.reference)]['analyse_region'][0]
        self.end = self.ANALYSE_REGION_DICT[(self.primer, self.reference)]['analyse_region'][1]
        self.coding_start = self.ANALYSE_REGION_DICT[(self.primer, self.reference)]['coding_region'][0]
        self.coding_end = self.ANALYSE_REGION_DICT[(self.primer, self.reference)]['coding_region'][1]


def extract_expected_variant_dict(expected_variants_path: str) -> dict:
    expected_variants = pd.read_excel(
        expected_variants_path, 'source_plt_all_target_plt_all'
    )[['Target Well', 'Mutations']].drop_duplicates()
    wells = [f"{s[0]}{int(s[1:]):02d}" for s in expected_variants['Target Well']]
    expected_variants_dict = dict(zip(wells, expected_variants['Mutations']))

    return expected_variants_dict


def align_abi_to_ref(gbkpath: str, abipath_list: List[str], start: Optional[int] = None, end: Optional[int] = None):
    template = QUEEN(record=gbkpath)

    if start is None and end is None:
        pass
    else:
        if start is None:
            start = 0
        if end is None:
            end = len(template.seq)

    abidata_list = []
    query_list = []
    for abifile_path in abipath_list:
        print(abifile_path)
        if ".ab1" == abifile_path[-4:]:
            abidata = abi_to_dict(abifile_path)
            query = generate_consensusseq(abidata)
            abidata_list.append(abidata)
            query_list.append(query)
    template_aligned, query_aligned_list, nongap, region_start, region_end = recursive_alignment(template, query_list)
    if region_start > start:
        nstart = 0
    else:
        nstart = start - region_start
    if region_end < end:
        nend = region_end - region_start
    else:
        nend = end - region_start
    m = 0
    n_list = []
    for n, char in enumerate(template_aligned.seq):
        if char == "-":
            pass
        else:
            if nstart <= m <= nend:
                n_list.append(n)
            else:
                pass
            m += 1
    new_abidata_list = []
    new_query_aligned_list = []
    for (aquery, strand), aabidata in zip(query_aligned_list, abidata_list):
        m = 0
        for n, char in enumerate(aquery):
            if n == min(n_list):
                ms = m
            if n == max(n_list):
                me = m
            if char == "-":
                pass
            else:
                m += 1
        new_query_aligned_list.append((aquery[min(n_list):max(n_list)], strand))
        new_abidata = copy.deepcopy(aabidata)
        if strand == 1:
            nstart = ms
            nend = me
        else:
            nstart, nend = (len(aabidata["conf"][0]) - me, len(aabidata["conf"][0]) - ms - 1)
        new_abidata["conf"][0] = new_abidata["conf"][0][nstart:nend]
        new_abidata["conf"][1] = new_abidata["conf"][1][nstart:nend]
        new_abidata["channel"]["A"][0] = new_abidata["channel"]["A"][0][nstart:nend]
        new_abidata["channel"]["A"][1] = new_abidata["channel"]["A"][1][nstart:nend]
        new_abidata["channel"]["T"][0] = new_abidata["channel"]["T"][0][nstart:nend]
        new_abidata["channel"]["T"][1] = new_abidata["channel"]["T"][1][nstart:nend]
        new_abidata["channel"]["G"][0] = new_abidata["channel"]["G"][0][nstart:nend]
        new_abidata["channel"]["G"][1] = new_abidata["channel"]["G"][1][nstart:nend]
        new_abidata["channel"]["C"][0] = new_abidata["channel"]["C"][0][nstart:nend]
        new_abidata["channel"]["C"][1] = new_abidata["channel"]["C"][1][nstart:nend]
        for nucl in "ATGC":
            positions = [p - aabidata["channel"][nucl][0][nstart] for p in aabidata["_channel"][nucl][0] if
                         aabidata["channel"][nucl][0][nstart] <= p <= aabidata["channel"][nucl][0][nend]]
            values = [v for p, v in zip(aabidata["_channel"][nucl][0], aabidata["_channel"][nucl][1]) if
                      aabidata["channel"][nucl][0][nstart] <= p <= aabidata["channel"][nucl][0][nend]]
            new_abidata["_channel"][nucl][0] = positions
            new_abidata["_channel"][nucl][1] = values
        new_abidata_list.append(new_abidata)
    template_aligned = cropdna(template_aligned, min(n_list), max(n_list), quinable=0)
    query_aligned_list = new_query_aligned_list
    query_aligned_list = list(list(zip(*query_aligned_list))[0])
    for i in range(len(query_aligned_list)):
        aquery = query_aligned_list[i]
        for j, char in enumerate(aquery):
            if char == "-":
                pass
            else:
                break
        for k, char in enumerate(aquery[::-1]):
            if char == "-":
                pass
            else:
                break
        if k == 0:
            newquery = j * "/" + aquery[j:]
        else:
            newquery = j * "/" + aquery[j:-1 * k] + k * "/"
        query_aligned_list[i] = newquery

    return template_aligned, query_aligned_list


def extract_region_of_circular_dna(start: int, end: int, dna_seq: str) -> str:
    """0-indexed."""
    if start < end:
        return dna_seq[start:end]
    else:
        return dna_seq[start:] + dna_seq[:end]


def modify_position_on_circular_dna(position: int, dna_length: int) -> int:
    return position - dna_length if position > dna_length else position


def extract_moclo_entities(plasmids: List[SeqRecord]) -> List[dict]:
    overhang_len = len(re.sub("^[ATCGN]+\\^|_[ATCGN]+$", "", BsaI.elucidate()))
    sel_fragments = []
    for rec in plasmids:
        dna_length = len(rec.seq)
        sites = BsaI.search(rec.seq, linear=False)
        # The position returned by the method search is the first base of the downstream segment produced by a
        # restriction (i.e. the first base after the position where the enzyme will cut).
        if len(sites) != 2:
            raise ValueError(
                f'BsaI needs to cut exactly two times on the target sequence, '
                f'but {len(sites)} are found on the target sequence at {sites}.'
            )
        for site_id, site in enumerate(sorted(sites)):

            if site_id == len(sites) - 1:  # last site
                fragment = extract_region_of_circular_dna(start=site - 1, end=sites[0] - 1, dna_seq=rec.seq)
            else:
                fragment = extract_region_of_circular_dna(start=site - 1, end=sites[site_id + 1] - 1, dna_seq=rec.seq)
            if len(BsaI.search(DNA_STUFFER + fragment + DNA_STUFFER)) != 0:
                # select fragment that does not contain BsaI sites
                continue
            # BsaI cuts at the position, leaving a specific overhang
            overhang5_start = site - 1  # Adjust for 0-based indexing
            overhang5_end = site + overhang_len - 1  # BsaI leaves a 4 bp overhang
            overhang3_start = site + len(fragment) - 1  # Adjust for 0-based indexing
            overhang3_end = site + len(fragment) + overhang_len - 1  # BsaI leaves a 4 bp overhang
            # Extract the overhang
            overhang5 = extract_region_of_circular_dna(
                modify_position_on_circular_dna(overhang5_start, dna_length),
                modify_position_on_circular_dna(overhang5_end, dna_length),
                rec.seq
            )
            overhang3 = extract_region_of_circular_dna(
                modify_position_on_circular_dna(overhang3_start, dna_length),
                modify_position_on_circular_dna(overhang3_end, dna_length),
                rec.seq
            )
            fragment_dict = {"Fragment": fragment, "Name": rec.description, "Overhang5": overhang5,
                             "Overhang3": overhang3}
            sel_fragments.append(fragment_dict)
    return sel_fragments


def find_overlap(frag1, frag2):
    return frag1["Overhang3"] == frag2["Overhang5"]


def order_moclo_fragments(moclo_fragments: List[dict]) -> List[dict]:
    """Orders DNA fragments based on overhangs."""
    ordered_fragments = []
    used = set()  # To keep track of used fragments

    # Start with the backbone fragment
    current_fragment = [frag for frag in moclo_fragments if re.search('backbone', frag['Name'], re.IGNORECASE)][0]
    ordered_fragments.append(current_fragment)
    used.add(current_fragment['Name'])

    while len(ordered_fragments) < len(moclo_fragments):
        for next_fragment in moclo_fragments:
            if next_fragment['Name'] not in used and find_overlap(current_fragment, next_fragment):
                ordered_fragments.append(next_fragment)
                used.add(next_fragment['Name'])
                current_fragment = next_fragment
                break
        else:
            # If no more overlaps are found, we can break the loop
            break

    return ordered_fragments


def generate_expected_variant_seqrecords(expected_variants_path: str, toolbox_plasmids: str) -> str:
    expected_variants = pd.read_excel(
        expected_variants_path, 'source_plt_all_target_plt_all'
    )[['Source Well', 'Target Well', 'Mutations']].groupby(['Target Well'])

    module_tbl = pd.read_excel(
        expected_variants_path, 'source module table'
    )[['Sequence Name', 'Well Position', 'Sequence']]

    toolbox_plasmids_records = list(SeqIO.parse(toolbox_plasmids, "gb"))
    toolbox_fragments = extract_moclo_entities(toolbox_plasmids_records)

    expected_variant_records = []
    for well, grouped in expected_variants:
        merged = pd.merge(left=grouped, right=module_tbl, how="left", left_on="Source Well", right_on="Well Position")
        mutation = grouped['Mutations'].drop_duplicates().values[0]
        sub_protein_fragments = [SeqRecord(seq=Seq(frag), description=seq_name) for frag, seq_name in
                                 zip(merged['Sequence'], merged['Sequence Name'])]
        sub_protein_fragments = extract_moclo_entities(sub_protein_fragments)
        ordered_fragments = order_moclo_fragments(moclo_fragments=toolbox_fragments + sub_protein_fragments)
        # Concat all digested fragments
        assembeled_dna = Seq('')
        start_position = 0
        features = []
        for frag in ordered_fragments:
            assembeled_dna = Seq('').join([assembeled_dna, frag["Fragment"]])
            end_position = start_position + len(frag["Fragment"])
            feature = SeqFeature(FeatureLocation(start_position, end_position), type="module", id=frag['Name'])
            features.append(feature)
            start_position = end_position
        expected_variant_record = SeqRecord(
            seq=assembeled_dna,
            id=f"{well[0][0]}{int(well[0][1:]):02d}",
            name=mutation,
            description=mutation,
            features=features,
            annotations={"molecule_type": "dna", "topology": "circular"}
        )
        expected_variant_record_shifted = expected_variant_record[500:] + expected_variant_record[:500]
        expected_variant_records.append(expected_variant_record_shifted)

        outfile = path.join(
            path.dirname(expected_variants_path),
            f"{expected_variant_record.id}_expected_variant_record.gb"
        )
        with open(outfile, "w") as handle:
            SeqIO.write(expected_variant_record_shifted, handle, "gb")

    return path.dirname(expected_variants_path)


def analyse_sanger_seq_res(seq_res_dir: path, ref_dir: str, out_dir: str, primer: str, ref_dna_name: str) -> str:
    seq_res_dir = path.join(seq_res_dir, primer)
    sequencing_res = SequencingResSetting(primer, ref_dna_name)
    res_dict_by_wel_position = aggregate_seq_res_and_ref_by_well_position(seq_res_dir, ref_dir)
    res_all = {}
    for well, ref_n_res_path_dict in res_dict_by_wel_position.items():
        if len(ref_n_res_path_dict["res_path"]) == 0:
            continue
        template_aligned, query_aligned_list = align_abi_to_ref(
            gbkpath=ref_n_res_path_dict["ref_path"],
            abipath_list=ref_n_res_path_dict["res_path"],
            start=sequencing_res.start,
            end=sequencing_res.end
        )
        res_per_well = analyse_seq_res_per_well(
            template_aligned=template_aligned,
            query_aligned_list=query_aligned_list,
            start=sequencing_res.start,
            abipath_list=ref_n_res_path_dict["res_path"]
        )
        res_all[well] = {'template_aligned': template_aligned.seq, 'positions': res_per_well}

    # Convert and write JSON object to file
    outpath = path.join(out_dir, f"{primer}_sanger_sequencing_analysis2.json")
    with open(outpath, "w") as outfile:
        json.dump(res_all, outfile)
    return outpath


def validate_clones(
        seq_analysis_file_primer_1: str,
        seq_analysis_file_primer_2: str,
        primer_1: str,
        primer_2: str,
        seq_res_dir: str,
        ref_dir: str
) -> Tuple[Dict[str, set], Dict[str, set]]:
    res_dict_by_wel_position = aggregate_seq_res_and_ref_by_well_position(seq_res_dir, ref_dir)

    with open(seq_analysis_file_primer_1) as json_file:
        seq_analysis_primer_1 = json.load(json_file)
    with open(seq_analysis_file_primer_2) as json_file:
        seq_analysis_primer_2 = json.load(json_file)
    validated_by_well: Dict[str, set] = {}
    not_validated_by_well: Dict[str, set] = {}
    for well in seq_analysis_primer_1.keys():
        all_clones = [path.basename(file_name) for file_name in res_dict_by_wel_position[well]["res_path"]]
        all_clones = [re.sub(f"-\\d_{primer_1}|-\\d_{primer_2}", "", file_name) for file_name in all_clones]
        intersecting_positions = (
                seq_analysis_primer_1[well]['positions'].keys() &
                seq_analysis_primer_2[well]['positions'].keys()
        )
        if len(intersecting_positions) == 0:  # no mutations are found.
            validated_by_well[well] = set(all_clones)
            continue
        validated_by_position: Dict[str, set] = {}
        for position in intersecting_positions:
            genotypes_primer1 = seq_analysis_primer_1[well]['positions'][position]
            genotypes_primer2 = seq_analysis_primer_2[well]['positions'][position]
            validated_at_this_position = \
                genotypes_primer1['query_genotypes'].get(genotypes_primer1['ref_genotype'], []) + \
                genotypes_primer2['query_genotypes'].get(genotypes_primer2['ref_genotype'], [])
            validated_at_this_position = [
                re.sub(f"-\\d_{primer_1}|-\\d_{primer_2}", "", file_name) for file_name in validated_at_this_position
            ]
            validated_by_position[position] = set(validated_at_this_position)
        validated_by_well[well] = set.intersection(*validated_by_position.values())
        not_validated_by_well[well] = set(all_clones) - validated_by_well[well]

    return validated_by_well, not_validated_by_well


def aggregate_seq_res_and_ref_by_well_position(seq_res_dir: str, ref_dir: str) -> dict:
    res_dict_by_wel_position = {}
    for ref_file in os.listdir(ref_dir):
        if not re.search("[.]gb$", ref_file): # skip non-Genbank file
            continue
        well_position = ref_file[:3]
        if well_position in res_dict_by_wel_position:
            raise ValueError(f"Multiple reference files were found for well: {well_position}!")    
        res_dict_by_wel_position[well_position] = {"ref_path": os.path.join(ref_dir, ref_file),
                                                   "res_path": []}

    for res_file in os.listdir(seq_res_dir):
        well_position = res_file[:3]
        if well_position not in res_dict_by_wel_position:
            print(f"No reference sequence file was found for sequencing result: {well_position}.")
        res_dict_by_wel_position[well_position]["res_path"].append(os.path.join(seq_res_dir, res_file))

    return res_dict_by_wel_position


def analyse_seq_res_per_well(template_aligned, query_aligned_list: list, start: int, abipath_list: list) -> dict:
    abifile_list = [os.path.basename(full_path) for full_path in abipath_list]
    positions = {}
    for idx in range(len(template_aligned.seq)):
        ref_nucl_at_current_idx = template_aligned.seq[idx]
        nucl_at_current_idx = [q_alned[idx] for q_alned in query_aligned_list if q_alned[idx] != '/']
        nucl_at_current_idx_set = set(nucl_at_current_idx)
        # if all nucleotide at current position are the same as reference, skip to next nucleotide
        if len(nucl_at_current_idx_set) == 1 and list(nucl_at_current_idx_set)[0] == ref_nucl_at_current_idx:
            continue
        query_genotypes = {}
        for qry_idx, nucl in enumerate(nucl_at_current_idx):
            if nucl not in query_genotypes:
                query_genotypes[nucl] = [abifile_list[qry_idx]]
            else:
                query_genotypes[nucl].append(abifile_list[qry_idx])

        positions[idx + start + 1] = {
            'ref_genotype': ref_nucl_at_current_idx,
            'query_genotypes': query_genotypes
        }
    return positions
