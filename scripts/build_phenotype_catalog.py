"""Render all retained phenotype rows without refitting or selecting outcomes."""

import json
from pathlib import Path
from string import Template

from data_io import read_csv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/context/phenotypes"
COORDINATE_COUNTS = (3, 6, 10, 20)
DOMAIN_NAMES = {
    "blood_chemistry": "Blood chemistry",
    "blood_pressure": "Blood pressure",
    "body_size": "Body measurements",
    "cognition": "Cognition",
    "personality_affect": "Personality and affect",
}


def index_label_rows(rows: list[dict[str, str]]) -> dict:
    """Group by domain and original label, requiring all four coordinate counts."""
    grouped = {}
    for row in rows:
        key = (row["domain"], row["label"])
        coordinate_count = int(row["k_identity_dims"])
        by_dimension = grouped.setdefault(key, {})
        assert coordinate_count not in by_dimension, (key, coordinate_count)
        by_dimension[coordinate_count] = row
    assert len(rows) == 1164 and len(grouped) == 291
    assert all(
        set(by_dimension) == set(COORDINATE_COUNTS) for by_dimension in grouped.values()
    )
    return grouped


def format_domain_table(domain_rows: list[dict[str, str]]) -> str:
    lines = [
        "| Domain | Fields | K = 3 | K = 6 | K = 10 | K = 20 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for domain, name in DOMAIN_NAMES.items():
        by_dimension = {
            int(row["k_identity_dims"]): row
            for row in domain_rows
            if row["domain"] == domain
        }
        assert set(by_dimension) == set(COORDINATE_COUNTS)
        values = [name, by_dimension[3]["n_labels"]]
        values.extend(
            f"{float(by_dimension[k]['loocv_r2']):+.4f}" for k in COORDINATE_COUNTS
        )
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def format_label_tables(grouped: dict) -> str:
    sections = []
    for domain, name in DOMAIN_NAMES.items():
        keys = sorted(key for key in grouped if key[0] == domain)
        lines = [
            f"### {name} · {len(keys)} fields",
            "",
            "| Source group / field | Observed n | K = 3 | K = 6 | K = 10 | K = 20 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for key in keys:
            by_dimension = grouped[key]
            assert len({row["n"] for row in by_dimension.values()}) == 1
            label = " / ".join(key[1].split("__")[-2:]).replace("|", "\\|")
            values = [label, by_dimension[3]["n"]]
            values.extend(
                f"{float(by_dimension[k]['cv_corr']):+.4f}" for k in COORDINATE_COUNTS
            )
            lines.append("| " + " | ".join(values) + " |")
        sections.append("\n".join(lines))
    return "\n\n".join(sections)


def coverage_receipt(
    label_rows: list[dict], domain_rows: list[dict], grouped: dict
) -> dict:
    return {
        "operation": "Retained-table catalog generation, not statistical recomputation",
        "source_label_rows": len(label_rows),
        "unique_fields": len(grouped),
        "dimension_settings": list(COORDINATE_COUNTS),
        "domain_summary_rows": len(domain_rows),
        "complete_label_dimension_coverage": True,
        "counts_by_domain": {
            domain: sum(key[0] == domain for key in grouped) for domain in DOMAIN_NAMES
        },
    }


def main() -> None:
    label_rows = read_csv(DATA / "domain_label_reconstruction.csv")
    domain_rows = read_csv(DATA / "domain_recovery_summary.csv")
    grouped = index_label_rows(label_rows)
    template = Template(
        (Path(__file__).parent / "templates/phenotype_atlas.md").read_text()
    )
    catalog = template.substitute(
        domain_table=format_domain_table(domain_rows),
        label_tables=format_label_tables(grouped),
    )
    (ROOT / "PHENOTYPE_ATLAS.md").write_text(catalog)
    receipt = coverage_receipt(label_rows, domain_rows, grouped)
    (DATA / "catalog_coverage.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
