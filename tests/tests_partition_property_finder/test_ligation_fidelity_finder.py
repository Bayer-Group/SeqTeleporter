from unittest import TestCase
import pandas as pd

from proseqteleporter.partition_property_finder.ligation_fidelity_finder import compute_ligation_fidelity


class TestComputeLigationFidelity(TestCase):
    def setUp(self) -> None:

        self.fidelity_data = pd.DataFrame(
            {
                'GCCC': [259, 0, 0, 0, 0],
                'AGTA': [0, 350, 0, 0, 0],
                'GGGC': [0, 0, 259, 0, 0],
                'TACT': [0, 0, 0, 350, 0],
                'AGCT': [0, 0, 0, 0, 299]
            }, index=['GGGC', 'TACT', 'GCCC', 'AGTA', 'AGCT']
        )
        self.fusion_site_indices = {site: idx for idx, site in enumerate(self.fidelity_data.index)}
        self.fusion_site_cols = {site: idx for idx, site in enumerate(self.fidelity_data.columns)}
        self.correct_lig_freq_dict = {site: self.fidelity_data.iloc[i, i] for site, i in self.fusion_site_indices.items()}


    def test_compute_ligation_fidelity_simple(self):
        all_fusion_sites_of_a_partition = ('GGGC', 'TACT')
        expected_output = 1
        output = compute_ligation_fidelity(all_fusion_sites_of_a_partition=all_fusion_sites_of_a_partition,
                                           fidelity_data=self.fidelity_data.values,
                                           fusion_site_indices=self.fusion_site_indices,
                                           fusion_site_cols=self.fusion_site_cols,
                                           correct_lig_freq_dict=self.correct_lig_freq_dict)
        # output = compute_ligation_fidelity(all_fusion_sites_of_a_partition, self.fidelity_data)
        self.assertEqual(expected_output, output)

    def test_compute_ligation_fidelity_reverse_compliment(self):
        all_fusion_sites_of_a_partition = ('GGGC', 'TACT', 'GCCC', 'AGTA')
        expected_output = 1
        output = compute_ligation_fidelity(all_fusion_sites_of_a_partition=all_fusion_sites_of_a_partition,
                                           fidelity_data=self.fidelity_data.values,
                                           fusion_site_indices=self.fusion_site_indices,
                                           fusion_site_cols=self.fusion_site_cols,
                                           correct_lig_freq_dict=self.correct_lig_freq_dict)
        # output = compute_ligation_fidelity(
        #     all_fusion_sites_of_a_partition, self.fidelity_data
        # )
        self.assertEqual(expected_output, output)

    def test_compute_ligation_fidelity_own_reverse_complement(self):
        all_fusion_sites_of_a_partition = ('TACT', 'GCCC', 'AGCT')
        expected_output = 0.5
        output = compute_ligation_fidelity(all_fusion_sites_of_a_partition=all_fusion_sites_of_a_partition,
                                           fidelity_data=self.fidelity_data.values,
                                           fusion_site_indices=self.fusion_site_indices,
                                           fusion_site_cols=self.fusion_site_cols,
                                           correct_lig_freq_dict=self.correct_lig_freq_dict)
        # output = compute_ligation_fidelity(all_fusion_sites_of_a_partition, self.fidelity_data)
        self.assertEqual(expected_output, output)

