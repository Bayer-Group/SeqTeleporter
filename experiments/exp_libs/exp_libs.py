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

# class MoCloFragment:
#
#     def __init__(self, seqrecord: SeqRecord):
#         self.parent_seqrecord = seqrecord
#         moclo_fragment_dict = self.extract_moclo_entity(seqrecord)
#         self.sequence = moclo_fragment_dict["Fragment"]
#         self.name = moclo_fragment_dict["Name"]
#         self.fusion_site5 = moclo_fragment_dict["Overhang5"]
#         self.fusion_site3 = moclo_fragment_dict["Overhang3"]
#
#     @staticmethod
#     def extract_region_of_circular_dna(start: int, end: int, dna_seq: str) -> str:
#         """0-indexed."""
#         if start < end:
#             return dna_seq[start:end]
#         else:
#             return dna_seq[start:] + dna_seq[:end]
#
#     @staticmethod
#     def modify_position_on_circular_dna(position: int, dna_length: int) -> int:
#         return position - dna_length if position > dna_length else position
#
#     def extract_moclo_entity(self, rec: SeqRecord) -> dict:
#         overhang_len = len(re.sub("^[ATCGN]+\\^|_[ATCGN]+$", "", BsaI.elucidate()))
#         dna_length = len(rec.seq)
#         sites = BsaI.search(rec.seq, linear=False)
#         # The position returned by the method search is the first base of the downstream segment produced by a
#         # restriction (i.e. the first base after the position where the enzyme will cut).
#         if len(sites) != 2:
#             raise ValueError(
#                 f'Error in record {rec}!'
#                 f'BsaI needs to cut exactly two times on the target sequence, '
#                 f'but {len(sites)} cut sites are found on the target sequence at {sites}.'
#             )
#         for site_id, site in enumerate(sorted(sites)):
#             if site_id == len(sites) - 1:  # last site
#                 fragment = self.extract_region_of_circular_dna(start=site - 1, end=sorted(sites)[0] - 1,
#                                                                dna_seq=rec.seq)
#             else:
#                 fragment = self.extract_region_of_circular_dna(start=site - 1, end=sorted(sites)[site_id + 1] - 1,
#                                                                dna_seq=rec.seq)
#             if len(BsaI.search(DNA_STUFFER + fragment + DNA_STUFFER)) != 0:
#                 # select fragment that does not contain BsaI sites
#                 continue
#             # BsaI cuts at the position, leaving a specific overhang
#             overhang5_start = site - 1  # Adjust for 0-based indexing
#             overhang5_end = site + overhang_len - 1  # BsaI leaves a 4 bp overhang
#             overhang3_start = site + len(fragment) - 1  # Adjust for 0-based indexing
#             overhang3_end = site + len(fragment) + overhang_len - 1  # BsaI leaves a 4 bp overhang
#             # Extract the overhang
#             overhang5 = self.extract_region_of_circular_dna(
#                 self.modify_position_on_circular_dna(overhang5_start, dna_length),
#                 self.modify_position_on_circular_dna(overhang5_end, dna_length),
#                 rec.seq
#             )
#             overhang3 = self.extract_region_of_circular_dna(
#                 self.modify_position_on_circular_dna(overhang3_start, dna_length),
#                 self.modify_position_on_circular_dna(overhang3_end, dna_length),
#                 rec.seq
#             )
#             fragment_dict = {"Fragment": fragment, "Name": rec.description, "Overhang5": overhang5,
#                              "Overhang3": overhang3}
#         return fragment_dict


class MoCloFragment:
    """
    A class to represent a MoClo (Modular Cloning) fragment derived from a DNA sequence.

    Attributes:
        parent_seqrecord (SeqRecord): The original sequence record from which the fragment is extracted.
        sequence (str): The DNA sequence of the MoClo fragment.
        name (str): The name/description of the MoClo fragment.
        fusion_site5 (str): The 5' overhang sequence of the MoClo fragment.
        fusion_site3 (str): The 3' overhang sequence of the MoClo fragment.
    """

    def __init__(self, seqrecord: SeqRecord):
        """
        Initializes a MoCloFragment object with the given sequence record.

        Args:
            seqrecord (SeqRecord): The sequence record containing the DNA sequence from which the fragment is extracted.
        """
        self.parent_seqrecord = seqrecord
        moclo_fragment_dict = self.extract_moclo_entity(seqrecord)
        self.sequence = moclo_fragment_dict["Fragment"]
        self.name = moclo_fragment_dict["Name"]
        self.fusion_site5 = moclo_fragment_dict["Overhang5"]
        self.fusion_site3 = moclo_fragment_dict["Overhang3"]

    @staticmethod
    def extract_region_of_circular_dna(start: int, end: int, dna_seq: str) -> str:
        """Extracts a region from a circular DNA sequence.

        Args:
            start (int): The starting index (0-indexed) of the region to extract.
            end (int): The ending index (0-indexed) of the region to extract.
            dna_seq (str): The circular DNA sequence.

        Returns:
            str: The extracted region of the DNA sequence.
        """
        if start < end:
            return dna_seq[start:end]
        else:
            return dna_seq[start:] + dna_seq[:end]

    @staticmethod
    def modify_position_on_circular_dna(position: int, dna_length: int) -> int:
        """Adjusts a position on a circular DNA sequence to ensure it falls within the valid range.

        Args:
            position (int): The original position to adjust.
            dna_length (int): The length of the DNA sequence.

        Returns:
            int: The adjusted position on the circular DNA sequence.
        """
        return position - dna_length if position > dna_length else position

    def extract_moclo_entity(self, rec: SeqRecord) -> dict:
        """Extracts the MoClo entity from the given sequence record.

        Args:
            rec (SeqRecord): The sequence record from which to extract the MoClo entity.

        Returns:
            dict: A dictionary containing the MoClo fragment, its name, and the overhang sequences.

        Raises:
            ValueError: If the number of BsaI cut sites found is not equal to two.
        """
        overhang_len = len(re.sub("^[ATCGN]+\\^|_[ATCGN]+$", "", BsaI.elucidate()))
        dna_length = len(rec.seq)
        sites = BsaI.search(rec.seq, linear=False)

        if len(sites) != 2:
            raise ValueError(
                f'Error in record {rec}!'
                f'BsaI needs to cut exactly two times on the target sequence, '
                f'but {len(sites)} cut sites are found on the target sequence at {sites}.'
            )

        for site_id, site in enumerate(sorted(sites)):
            if site_id == len(sites) - 1:  # last site
                fragment = self.extract_region_of_circular_dna(start=site - 1, end=sorted(sites)[0] - 1,
                                                               dna_seq=rec.seq)
            else:
                fragment = self.extract_region_of_circular_dna(start=site - 1, end=sorted(sites)[site_id + 1] - 1,
                                                               dna_seq=rec.seq)
            if len(BsaI.search(DNA_STUFFER + fragment + DNA_STUFFER)) != 0:
                # select fragment that does not contain BsaI sites
                continue

            # BsaI cuts at the position, leaving a specific overhang
            overhang5_start = site - 1  # Adjust for 0-based indexing
            overhang5_end = site + overhang_len - 1  # BsaI leaves a 4 bp overhang
            overhang3_start = site + len(fragment) - 1  # Adjust for 0-based indexing
            overhang3_end = site + len(fragment) + overhang_len - 1  # BsaI leaves a 4 bp overhang

            # Extract the overhang
            overhang5 = self.extract_region_of_circular_dna(
                self.modify_position_on_circular_dna(overhang5_start, dna_length),
                self.modify_position_on_circular_dna(overhang5_end, dna_length),
                rec.seq
            )
            overhang3 = self.extract_region_of_circular_dna(
                self.modify_position_on_circular_dna(overhang3_start, dna_length),
                self.modify_position_on_circular_dna(overhang3_end, dna_length),
                rec.seq
            )
            fragment_dict = {
                "Fragment": fragment,
                "Name": rec.description,
                "Overhang5": overhang5,
                "Overhang3": overhang3
            }
        return fragment_dict

# class PlasmidAssembler:
#
#     @staticmethod
#     def export_assembled_plasmids(records: list, output_file_path: str) -> None:
#         with open(output_file_path, "w") as handle:
#             SeqIO.write(records, handle, "gb")
#
#     @staticmethod
#     def find_overlap(frag1: MoCloFragment, frag2: MoCloFragment) -> bool:
#         return frag1.fusion_site3 == frag2.fusion_site5
#
#     @staticmethod
#     def validate_ordered_fragments_group(ordered_frags_group: list[MoCloFragment]) -> bool:
#         # validate that all fusion sites fit
#         # validate that the lenth of ordered_fragments_group is the same as the length of the frag name group
#         is_fitting = True
#         for idx, f in enumerate(ordered_frags_group):
#             if idx == len(ordered_frags_group) - 1:
#                 if f.fusion_site3 != ordered_frags_group[0].fusion_site5:
#                     is_fitting = False
#                     print(f)
#                 continue
#             if f.fusion_site3 != ordered_frags_group[idx + 1].fusion_site5:
#                 is_fitting = False
#                 print(f)
#         return is_fitting
#
#     @staticmethod
#     def group_fragments(grouped_frag_names: list[str], all_frags: list[MoCloFragment]) -> list[MoCloFragment]:
#         frags_group = []
#         for frag_name in grouped_frag_names:
#             for f in all_frags:
#                 if frag_name == f.name:
#                     frags_group.append(f)
#         if len(frags_group) != len(grouped_frag_names):
#             raise ValueError(f'Inconsistent fragment number in {grouped_frag_names} and \n {frags_group}')
#         return frags_group
#
#     def order_moclo_fragments(self, moclo_fragments: List[MoCloFragment]) -> List[MoCloFragment]:
#         """Orders DNA fragments based on overhangs."""
#         ordered_fragments = []
#         used = set()  # To keep track of used fragments
#
#         # Start with the backbone fragment
#         current_fragment = [frag for frag in moclo_fragments if re.search('backbone', frag.name, re.IGNORECASE)][0]
#         ordered_fragments.append(current_fragment)
#         used.add(current_fragment.name)
#
#         while len(ordered_fragments) < len(moclo_fragments):
#             for next_fragment in moclo_fragments:
#                 if next_fragment.name not in used and self.find_overlap(current_fragment, next_fragment):
#                     ordered_fragments.append(next_fragment)
#                     used.add(next_fragment.name)
#                     current_fragment = next_fragment
#                     break
#             else:
#                 # If no more overlaps are found, we can break the loop
#                 break
#
#         return ordered_fragments
#
#     # Concat all digested fragments
#     def assemble_digested_frags(self, ordered_fragments: list[MoCloFragment]) -> SeqRecord:
#         assembeled_dna = Seq('')
#         start_position = 0
#         features = []
#         name_ls = []
#         for frag in ordered_fragments:
#             assembeled_dna = Seq('').join([assembeled_dna, frag.sequence])
#             end_position = start_position + len(frag.sequence)
#             feature = SeqFeature(FeatureLocation(start_position, end_position), type="module", id=frag.name)
#             features.append(feature)
#             start_position = end_position
#             if not re.search('backbone', frag.name):
#                 name_ls.append(frag.name)
#         name = '_'.join(name_ls)
#         expected_variant_record = SeqRecord(
#             seq=assembeled_dna,
#             id=name,
#             name=name,
#             description=name,
#             features=features,
#             annotations={"molecule_type": "dna", "topology": "circular"}
#         )
#         expected_variant_record_shifted = expected_variant_record[500:] + expected_variant_record[:500]
#         return expected_variant_record_shifted
#
#     # group fragments based on grouped modules and then assemble
#     def assemble_plasmids(
#             self,
#             frag_name_groups: list[list[str]],
#             all_frags: list[MoCloFragment],
#             bb_fragments: list[MoCloFragment]
#     ) -> list[SeqRecord]:
#         assembled_plasmids = []
#         for group in frag_name_groups:
#             frags_group = self.group_fragments(group, all_frags)
#             ordered_frags_group = self.order_moclo_fragments(moclo_fragments=frags_group + bb_fragments)
#             if not self.validate_ordered_fragments_group(ordered_frags_group):
#                 raise ValueError(f'Incompatible fusion sites in {group} \n {frags_group + bb_fragments}')
#             assembled = self.assemble_digested_frags(ordered_frags_group)
#             assembled_plasmids.append(assembled)
#
#         return assembled_plasmids
#
#     def generate_expected_variant_seqrecords(self, expected_variants_path: str, toolbox_plasmids: str) -> str:
#         expected_variants = pd.read_excel(
#             expected_variants_path, 'source_plt_all_target_plt_all'
#         )[['Source Well', 'Target Well', 'Mutations']].groupby(['Target Well'])
#
#         module_tbl = pd.read_excel(
#             expected_variants_path, 'source module table'
#         )[['Sequence Name', 'Well Position', 'Sequence']]
#
#         toolbox_plasmids_records = list(SeqIO.parse(toolbox_plasmids, "gb"))
#         toolbox_fragments = extract_moclo_entities(toolbox_plasmids_records)
#
#         expected_variant_records = []
#         for well, grouped in expected_variants:
#             merged = pd.merge(left=grouped, right=module_tbl, how="left", left_on="Source Well",
#                               right_on="Well Position")
#             mutation = grouped['Mutations'].drop_duplicates().values[0]
#             sub_protein_fragments_seqrecords = [SeqRecord(seq=Seq(frag), description=seq_name) for frag, seq_name in
#                                                 zip(merged['Sequence'], merged['Sequence Name'])]
#             sub_protein_fragments = extract_moclo_entities(sub_protein_fragments_seqrecords)
#             ordered_fragments = self.order_moclo_fragments(
#                 moclo_fragments=toolbox_fragments + sub_protein_fragments)
#             # Concat all digested fragments
#             assembeled_dna = Seq('')
#             start_position = 0
#             features = []
#             for frag in ordered_fragments:
#                 assembeled_dna = Seq('').join([assembeled_dna, frag.sequence])
#                 end_position = start_position + len(frag.sequence)
#                 feature = SeqFeature(FeatureLocation(start_position, end_position), type="module", id=frag.name)
#                 features.append(feature)
#                 start_position = end_position
#             expected_variant_record = SeqRecord(
#                 seq=assembeled_dna,
#                 id=f"{well[0][0]}{int(well[0][1:]):02d}",
#                 name=mutation,
#                 description=mutation,
#                 features=features,
#                 annotations={"molecule_type": "dna", "topology": "circular"}
#             )
#             expected_variant_record_shifted = expected_variant_record[500:] + expected_variant_record[:500]
#             expected_variant_records.append(expected_variant_record_shifted)
#
#             outfile = path.join(
#                 path.dirname(expected_variants_path),
#                 f"{expected_variant_record.id}_expected_variant_record.gb"
#             )
#             with open(outfile, "w") as handle:
#                 SeqIO.write(expected_variant_record_shifted, handle, "gb")
#
#         return path.dirname(expected_variants_path)
class PlasmidAssembler:
    """
    A class to assemble plasmids from modular cloning fragments.

    This class provides methods for assembling DNA fragments into plasmids, validating the order of fragments,
    and exporting assembled plasmids. It handles the grouping of fragments based on their names and ensures
    that the fusion sites are compatible for successful assembly.

    Methods:
        export_assembled_plasmids(records: list, output_file_path: str) -> None:
            Exports assembled plasmids to a specified file in GenBank format.

        find_overlap(frag1: MoCloFragment, frag2: MoCloFragment) -> bool:
            Checks if two MoClo fragments overlap based on their fusion sites.

        validate_ordered_fragments_group(ordered_frags_group: list[MoCloFragment]) -> bool:
            Validates that all fragments in the ordered group fit together based on their fusion sites.

        group_fragments(grouped_frag_names: list[str], all_frags: list[MoCloFragment]) -> list[MoCloFragment]:
            Groups fragments based on a list of fragment names.

        order_moclo_fragments(moclo_fragments: List[MoCloFragment]) -> List[MoCloFragment]:
            Orders MoClo fragments based on their overhangs.

        assemble_digested_frags(ordered_fragments: list[MoCloFragment]) -> SeqRecord:
            Assembles the DNA sequence from ordered fragments and creates a SeqRecord.

        assemble_plasmids(frag_name_groups: list[list[str]], all_frags: list[MoCloFragment],
                          bb_fragments: list[MoCloFragment]) -> list[SeqRecord]:
            Assembles plasmids from groups of fragments and backbone fragments.

        generate_expected_variant_seqrecords(expected_variants_path: str, toolbox_plasmids: str) -> str:
            Generates expected variant SeqRecords based on provided input files.
    """

    @staticmethod
    def export_assembled_plasmids(records: list, output_file_path: str) -> None:
        """Exports assembled plasmids to a specified file in GenBank format.

        Args:
            records (list): A list of SeqRecord objects representing the assembled plasmids.
            output_file_path (str): The file path where the assembled plasmids will be saved.
        """
        with open(output_file_path, "w") as handle:
            SeqIO.write(records, handle, "gb")

    @staticmethod
    def find_overlap(frag1: MoCloFragment, frag2: MoCloFragment) -> bool:
        """Checks if two MoClo fragments overlap based on their fusion sites.

        Args:
            frag1 (MoCloFragment): The first MoClo fragment.
            frag2 (MoCloFragment): The second MoClo fragment.

        Returns:
            bool: True if the fragments overlap, False otherwise.
        """
        return frag1.fusion_site3 == frag2.fusion_site5

    @staticmethod
    def validate_ordered_fragments_group(ordered_frags_group: list[MoCloFragment]) -> bool:
        """Validates that all fragments in the ordered group fit together based on their fusion sites.

        Args:
            ordered_frags_group (list[MoCloFragment]): A list of ordered MoClo fragments.

        Returns:
            bool: True if all fragments fit together, False otherwise.
        """
        is_fitting = True
        for idx, f in enumerate(ordered_frags_group):
            if idx == len(ordered_frags_group) - 1:
                if f.fusion_site3 != ordered_frags_group[0].fusion_site5:
                    is_fitting = False
                    print(f)
                continue
            if f.fusion_site3 != ordered_frags_group[idx + 1].fusion_site5:
                is_fitting = False
                print(f)
        return is_fitting

    @staticmethod
    def group_fragments(grouped_frag_names: list[str], all_frags: list[MoCloFragment]) -> list[MoCloFragment]:
        """Groups fragments based on a list of fragment names.

        Args:
            grouped_frag_names (list[str]): A list of fragment names to group.
            all_frags (list[MoCloFragment]): A list of all available MoClo fragments.

        Returns:
            list[MoCloFragment]: A list of MoClo fragments that match the given names.

        Raises:
            ValueError: If the number of grouped fragments does not match the names provided.
        """
        frags_group = []
        for frag_name in grouped_frag_names:
            for f in all_frags:
                if frag_name == f.name:
                    frags_group.append(f)
        if len(frags_group) != len(grouped_frag_names):
            raise ValueError(f'Inconsistent fragment number in {grouped_frag_names} and \n {frags_group}')
        return frags_group

    def order_moclo_fragments(self, moclo_fragments: List[MoCloFragment]) -> List[MoCloFragment]:
        """Orders MoClo fragments based on overhangs.

        Args:
            moclo_fragments (List[MoCloFragment]): A list of MoClo fragments to be ordered.

        Returns:
            List[MoCloFragment]: A list of ordered MoClo fragments.
        """
        ordered_fragments = []
        used = set()  # To keep track of used fragments

        # Start with the backbone fragment
        current_fragment = [frag for frag in moclo_fragments if re.search('backbone', frag.name, re.IGNORECASE)][0]
        ordered_fragments.append(current_fragment)
        used.add(current_fragment.name)

        while len(ordered_fragments) < len(moclo_fragments):
            for next_fragment in moclo_fragments:
                if next_fragment.name not in used and self.find_overlap(current_fragment, next_fragment):
                    ordered_fragments.append(next_fragment)
                    used.add(next_fragment.name)
                    current_fragment = next_fragment
                    break
            else:
                # If no more overlaps are found, we can break the loop
                break

        return ordered_fragments

    def assemble_digested_frags(self, ordered_fragments: list[MoCloFragment]) -> SeqRecord:
        """Assembles the DNA sequence from ordered fragments and creates a SeqRecord.

        Args:
            ordered_fragments (list[MoCloFragment]): A list of ordered MoClo fragments.

        Returns:
            SeqRecord: A SeqRecord object representing the assembled DNA sequence.
        """
        assembled_dna = Seq('')
        start_position = 0
        features = []
        name_ls = []
        for frag in ordered_fragments:
            assembled_dna = Seq('').join([assembled_dna, frag.sequence])
            end_position = start_position + len(frag.sequence)
            feature = SeqFeature(FeatureLocation(start_position, end_position), type="module", id=frag.name)
            features.append(feature)
            start_position = end_position
            if not re.search('backbone', frag.name):
                name_ls.append(frag.name)
        name = '_'.join(name_ls)
        expected_variant_record = SeqRecord(
            seq=assembled_dna,
            id=name,
            name=name,
            description=name,
            features=features,
            annotations={"molecule_type": "dna", "topology": "circular"}
        )
        expected_variant_record_shifted = expected_variant_record[500:] + expected_variant_record[:500]
        return expected_variant_record_shifted

    def assemble_plasmids(
            self,
            frag_name_groups: list[list[str]],
            all_frags: list[MoCloFragment],
            bb_fragments: list[MoCloFragment]
    ) -> list[SeqRecord]:
        """Assembles plasmids from groups of fragments and backbone fragments.

        Args:
            frag_name_groups (list[list[str]]): A list of lists, where each inner list contains names of fragments to assemble.
            all_frags (list[MoCloFragment]): A list of all available MoClo fragments.
            bb_fragments (list[MoCloFragment]): A list of backbone fragments.

        Returns:
            list[SeqRecord]: A list of SeqRecord objects representing the assembled plasmids.

        Raises:
            ValueError: If the fusion sites in the ordered fragments group are incompatible.
        """
        assembled_plasmids = []
        for group in frag_name_groups:
            frags_group = self.group_fragments(group, all_frags)
            ordered_frags_group = self.order_moclo_fragments(moclo_fragments=frags_group + bb_fragments)
            if not self.validate_ordered_fragments_group(ordered_frags_group):
                raise ValueError(f'Incompatible fusion sites in {group} \n {frags_group + bb_fragments}')
            assembled = self.assemble_digested_frags(ordered_frags_group)
            assembled_plasmids.append(assembled)

        return assembled_plasmids

    def generate_expected_variant_seqrecords(self, expected_variants_path: str, toolbox_plasmids: str) -> str:
        """Generates expected variant SeqRecords based on provided input files.

        Args:
            expected_variants_path (str): The file path to the Excel file containing expected variants.
            toolbox_plasmids (str): The file path to the toolbox plasmids in GenBank format.

        Returns:
            str: The directory path where expected variant records are saved.
        """
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
            merged = pd.merge(left=grouped, right=module_tbl, how="left", left_on="Source Well",
                              right_on="Well Position")
            mutation = grouped['Mutations'].drop_duplicates().values[0]
            sub_protein_fragments_seqrecords = [SeqRecord(seq=Seq(frag), description=seq_name) for frag, seq_name in
                                                zip(merged['Sequence'], merged['Sequence Name'])]
            sub_protein_fragments = extract_moclo_entities(sub_protein_fragments_seqrecords)
            ordered_fragments = self.order_moclo_fragments(
                moclo_fragments=toolbox_fragments + sub_protein_fragments)
            # Concat all digested fragments
            assembled_dna = Seq('')
            start_position = 0
            features = []
            for frag in ordered_fragments:
                assembled_dna = Seq('').join([assembled_dna, frag.sequence])
                end_position = start_position + len(frag.sequence)
                feature = SeqFeature(FeatureLocation(start_position, end_position), type="module", id=frag.name)
                features.append(feature)
                start_position = end_position
            expected_variant_record = SeqRecord(
                seq=assembled_dna,
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
    template_aligned, query_aligned_list, nongap, region_start, region_end, cannot_align = recursive_alignment(
        template,
        query_list
    )
    # Removing elements by index in reverse order
    for qry_id in sorted(cannot_align, reverse=True):
        abidata_list.pop(qry_id)
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
    if len(n_list) == 0:
        cannot_align = [i for i, data in enumerate(abidata_list)]
        return "", [], cannot_align
    
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

    return template_aligned.seq, query_aligned_list, cannot_align


def extract_moclo_entities(plasmids: List[SeqRecord]) -> List[MoCloFragment]:
    sel_fragments = []
    for rec in plasmids:
        moclo_fragment = MoCloFragment(rec)
        sel_fragments.append(moclo_fragment)
    return sel_fragments


def analyse_sanger_seq_res(
        seq_res_dir: str,
        filter_seq_res_by: str,
        ref_dir: str,
        out_dir: str,
        primer: str,
        ref_dna_name: str
) -> str:
    seq_res_dir = path.join(seq_res_dir, primer)
    sequencing_res = SequencingResSetting(primer, ref_dna_name)
    res_dict_by_wel_position = aggregate_seq_res_and_ref_by_well_position(seq_res_dir, filter_seq_res_by, ref_dir)
    res_all = {}
    for well, ref_n_res_path_dict in res_dict_by_wel_position.items():
        abipath_list = ref_n_res_path_dict["res_path"]
        if len(ref_n_res_path_dict["res_path"]) == 0:
            continue
        template_aligned_seq, query_aligned_list, cannot_align = align_abi_to_ref(
            gbkpath=ref_n_res_path_dict["ref_path"],
            abipath_list=abipath_list,
            start=sequencing_res.start,
            end=sequencing_res.end
        )
        if len(template_aligned_seq) == 0:
            res_per_well = {}
        else:
            res_per_well = analyse_seq_res_per_well(
                template_aligned_seq=template_aligned_seq,
                query_aligned_list=query_aligned_list,
                start=sequencing_res.start,
                abipath_list=abipath_list,
                cannot_align=cannot_align
            )
        cannot_align_files = [
            os.path.basename(full_path) for id, full_path in enumerate(abipath_list) if id in cannot_align
        ]
        res_all[well] = {
            'template_aligned': template_aligned_seq,
            'positions': res_per_well,
            'cannot_align': cannot_align_files
        }

    # Convert and write JSON object to file
    filter_seq_res_by_ = re.sub("[*.:$|^(?!)]", "_", filter_seq_res_by)
    outpath = path.join(out_dir, f"{primer}_sanger_sequencing_analysis_{filter_seq_res_by_}.json")
    with open(outpath, "w") as outfile:
        json.dump(res_all, outfile)
    return outpath


def validate_clones(
        seq_analysis_file_primer_1: str,
        seq_analysis_file_primer_2: str,
        primer_1: str,
        primer_2: str,
        seq_res_dir: str,
        filter_seq_res_by: str,
        ref_dir: str
) -> Tuple[Dict[str, set], Dict[str, set]]:
    res_dict_by_wel_position = aggregate_seq_res_and_ref_by_well_position(seq_res_dir, filter_seq_res_by, ref_dir)

    with open(seq_analysis_file_primer_1) as json_file:
        seq_analysis_primer_1 = json.load(json_file)
    with open(seq_analysis_file_primer_2) as json_file:
        seq_analysis_primer_2 = json.load(json_file)
    validated_by_well: Dict[str, set] = {}
    not_validated_by_well: Dict[str, set] = {}
    for well in seq_analysis_primer_1.keys():
        all_clones = [path.basename(file_name) for file_name in res_dict_by_wel_position[well]["res_path"]]
        all_clones = [re.sub(f"-\\d_{primer_1}|-\\d_{primer_2}", "", file_name) for file_name in all_clones]
        if (seq_analysis_primer_1[well]['template_aligned'] == "" and
                seq_analysis_primer_2[well]['template_aligned'] == ""):
            validated_by_well[well] = set()
            continue
        if (seq_analysis_primer_1[well]['template_aligned'] == "" and
                seq_analysis_primer_2[well]['template_aligned'] != ""):
            intersecting_positions = seq_analysis_primer_2[well]['positions'].keys()
        if (seq_analysis_primer_1[well]['template_aligned'] != "" and
                seq_analysis_primer_2[well]['template_aligned'] == ""):
            intersecting_positions = seq_analysis_primer_1[well]['positions'].keys()
        else:
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


def aggregate_seq_res_and_ref_by_well_position(seq_res_dir: str, filter_seq_res_by: str, ref_dir: str) -> dict:
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
        if well_position[-1] == "-":
            well_position = well_position[0] + "0" + well_position[1]

        if well_position not in res_dict_by_wel_position:
            print(f"No reference sequence file was found for sequencing result: {well_position}.")
        if re.search(filter_seq_res_by, res_file):
            res_dict_by_wel_position[well_position]["res_path"].append(os.path.join(seq_res_dir, res_file))

    return res_dict_by_wel_position


def analyse_seq_res_per_well(
        template_aligned_seq,
        query_aligned_list: list,
        start: int,
        abipath_list: list,
        cannot_align: list
) -> dict:
    abifile_list = [os.path.basename(full_path) for id, full_path in enumerate(abipath_list) if id not in cannot_align]
    positions = {}
    for idx in range(len(template_aligned_seq)):
        ref_nucl_at_current_idx = template_aligned_seq[idx]
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


def make_plate_layout(validated_by_well: dict) -> pd.DataFrame:
    data = {}
    seen = []
    seen_idx = []
    index = []
    for well, clones in validated_by_well.items():
        if well[:1] not in seen_idx:
            index.append(well[:1])
        seen_idx.append(well[:1])
        col = well[1:]
        if col not in seen:
            data[col] = [clones]
        else:
            data[col].append(clones)
        seen.append(col)
    return pd.DataFrame(data, index=index)

