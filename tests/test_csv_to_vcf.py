import importlib.util
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_TO_VCF_PATH = REPO_ROOT / "workflow" / "scripts" / "csv_to_vcf.py"


def load_csv_to_vcf_module():
    spec = importlib.util.spec_from_file_location("csv_to_vcf", CSV_TO_VCF_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


csv_to_vcf = load_csv_to_vcf_module()


class CsvToVcfTests(unittest.TestCase):
    def test_build_mapping_detects_advp_aliases(self):
        df = pd.DataFrame(
            columns=[
                "#dbSNP_hg38_chr",
                "dbSNP_hg38_position",
                "nonref_allele",
                "Top SNP",
                "P-value",
            ]
        )
        mapping = csv_to_vcf.build_mapping(df)
        self.assertEqual(mapping["CHROM"], "#dbSNP_hg38_chr")
        self.assertEqual(mapping["POS"], "dbSNP_hg38_position")
        self.assertEqual(mapping["ALT"], "nonref_allele")
        self.assertEqual(mapping["ID"], "Top SNP")
        self.assertEqual(mapping["P"], "P-value")
        self.assertNotIn("REF", mapping)

    def test_derive_missing_ref_fills_ref_and_normalizes_contig(self):
        df_norm = pd.DataFrame(
            [
                {"CHROM": "19", "POS": "44906745", "REF": "", "ALT": "A"},
                {"CHROM": "chr1", "POS": "100", "REF": "C", "ALT": "T"},
            ]
        )
        report_rows = []

        with mock.patch.object(csv_to_vcf, "load_fasta_contigs", return_value={"chr19", "chr1"}), mock.patch.object(
            csv_to_vcf,
            "fetch_reference_bases",
            return_value={"chr19:44906745-44906745": "G"},
        ):
            unresolved = csv_to_vcf.derive_missing_ref(df_norm, "fake.fa", report_rows)

        self.assertEqual(unresolved, [])
        self.assertEqual(df_norm.loc[0, "CHROM"], "chr19")
        self.assertEqual(df_norm.loc[0, "REF"], "G")
        self.assertTrue(any(item == "ref_derivation" and status == "OK" for item, status, _ in report_rows))

    def test_derive_missing_ref_reports_non_snp_alt_as_unresolved(self):
        df_norm = pd.DataFrame([{"CHROM": "1", "POS": "10", "REF": "", "ALT": "NR"}])
        report_rows = []

        with mock.patch.object(csv_to_vcf, "load_fasta_contigs", return_value={"1"}):
            unresolved = csv_to_vcf.derive_missing_ref(df_norm, "fake.fa", report_rows)

        self.assertEqual(len(unresolved), 1)
        self.assertIn("cannot derive REF", unresolved[0][1])
        self.assertEqual(csv_to_vcf.clean_text(df_norm.loc[0, "REF"]), "")

    def test_normalize_rows_for_vcf_drops_invalid_rows(self):
        df_norm = pd.DataFrame(
            [
                {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "T"},
                {"CHROM": "1", "POS": "bad", "REF": "A", "ALT": "T"},
                {"CHROM": "1", "POS": "101", "REF": "G", "ALT": "G"},
                {"CHROM": "", "POS": "102", "REF": "C", "ALT": "A"},
            ]
        )
        report_rows = []

        df_vcf, row_errors = csv_to_vcf.normalize_rows_for_vcf(df_norm, report_rows)

        self.assertEqual(len(df_vcf), 1)
        self.assertEqual(len(row_errors), 3)
        self.assertEqual(df_vcf.iloc[0]["POS"], "100")
        self.assertTrue(any(item == "row_validation" and status == "WARNING" for item, status, _ in report_rows))

    def test_build_normalized_table_marks_missing_ref_as_warning(self):
        df = pd.DataFrame(
            [
                {
                    "#dbSNP_hg38_chr": "19",
                    "dbSNP_hg38_position": "44906745",
                    "nonref_allele": "A",
                }
            ]
        )
        mapping = csv_to_vcf.build_mapping(df)
        report_rows = []
        df_norm, missing_required = csv_to_vcf.build_normalized_table(df, mapping, report_rows)

        self.assertIn("CHROM", df_norm.columns)
        self.assertIn("POS", df_norm.columns)
        self.assertIn("ALT", df_norm.columns)
        self.assertEqual(missing_required, ["REF"])
        self.assertTrue(
            any(
                item == "required_fields"
                and status == "WARNING"
                and "Missing required fields" in detail
                for item, status, detail in report_rows
            )
        )


if __name__ == "__main__":
    unittest.main()
