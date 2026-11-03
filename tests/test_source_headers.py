"""Source files change layout between releases: a renamed column must fail the build, a moved one must be read right."""

import openpyxl
import pytest

from kraken.harmonizers.loinc import LOINC_FILES, LoincHarmonizer
from kraken.harmonizers.pgs_catalog import PGSCatalogHarmonizer


def test_loinc_fails_when_a_column_is_renamed(tmp_path):
    # SCALE_TYP is read only inside the category callable, so it is the easiest one to leave unchecked
    loinc_csv = LOINC_FILES[0]
    path = tmp_path / "Loinc.csv"
    path.write_text(",".join(sorted(loinc_csv.columns - {"SCALE_TYP"})) + "\n")
    with pytest.raises(ValueError, match="SCALE_TYP"):
        next(LoincHarmonizer._read_csv(path, loinc_csv.columns))


def test_pgs_reads_columns_by_name():
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "Performance Metrics"
    sheet.append(["PGS Sample Set (PSS)", "Evaluated Score"])  # the reverse of the real sheet's order
    sheet.append(["PSS000001", "PGS000001"])
    sheet.append(["PSS000002", "PGS000001"])
    harmonizer = object.__new__(PGSCatalogHarmonizer)  # skips __init__, which loads Biolink and biomapper2
    assert harmonizer._load_evaluations(workbook) == {"PGS000001": {"PSS000001", "PSS000002"}}

    sheet["B1"] = "Score"
    with pytest.raises(KeyError, match="Evaluated Score"):
        harmonizer._load_evaluations(workbook)
