import random
import tempfile
import unittest
from pathlib import Path

from scripts.host_network import (add_identical_profiles, bootstrap, category,
                                  normalize_host, permutation, score_and_select)


class HostNetworkTests(unittest.TestCase):
    def test_excludes_nonbiological_hosts(self):
        self.assertEqual(normalize_host("Dermacentor silvarum")[1], "tick_environment_host")
        self.assertEqual(normalize_host("Snownnecut virus")[0], "")
        self.assertEqual(normalize_host("Coniella diplodiella")[0], "Coniella")
        self.assertEqual(normalize_host("Vitis vinifera")[0], "Vitis")
        self.assertEqual(normalize_host("Yellow silver pine")[0], "Yellow silver pine")
        self.assertEqual(category("Yellow silver pine"), "plant-associated")

    def test_score_and_statistics_are_reproducible(self):
        rows = []
        for query, family_name in (("CV1", "Chrysoviridae"), ("BV1", "bunya-like lineage")):
            for index, genus in enumerate(("Fusarium", "Coniella")):
                rows.append({"query": query, "family": family_name, "accession": f"A{query}{index}",
                             "host_genus": genus, "host_category": "fungal-associated",
                             "threshold_status": True, "bitscore": 100 - index,
                             "identity": 0.5, "query_coverage": 0.8})
        accepted, selected = score_and_select(rows, 10)
        self.assertEqual(len(accepted), 4)
        self.assertEqual(len(selected), 4)
        support_a = bootstrap(selected, 200, random.Random(4))
        support_b = bootstrap(selected, 200, random.Random(4))
        self.assertEqual(support_a, support_b)
        enrichment = permutation(selected, 200, random.Random(4))
        self.assertTrue(enrichment)

    def test_cv2_is_explicitly_derived_from_cv1(self):
        rows = [{"query": "CV1", "manuscript_virus": "TsCV1",
                 "candidate_labels": "TsCV1", "family": "Chrysoviridae",
                 "accession": "A1", "profile_provenance": "observed_input"}]
        derived = add_identical_profiles(rows)
        self.assertEqual([row["query"] for row in derived], ["CV1", "CV2"])
        self.assertEqual(derived[-1]["profile_provenance"], "copied_from_CV1_identical_RdRp")


if __name__ == "__main__":
    unittest.main()
