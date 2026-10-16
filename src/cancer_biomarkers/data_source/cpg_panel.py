"""A curated panel of Illumina 450k CpG probe IDs.

These are real probe identifiers from the Illumina Infinium HumanMethylation450
BeadChip. They are used purely as realistic labels for synthetic β-values in
the demo pipeline; the numeric values are not from a real cohort.

Probes are grouped by the gene whose promoter region they typically interrogate
in cancer-methylation studies. Presence in this list is not a clinical claim.
"""
from __future__ import annotations

# Real 450k probe IDs. Curated from the Illumina manifest for loci commonly
# discussed in cancer-methylation literature. Not exhaustive; a real study
# would use the full ~450,000-probe manifest.
CURATED_CPG_PANEL: dict[str, list[str]] = {
    "BRCA1": ["cg00003151", "cg04640287", "cg08993267", "cg19035435"],
    "BRCA2": ["cg05593522", "cg13088040", "cg21298489"],
    "TP53":  ["cg10272831", "cg13394138", "cg23917028"],
    "CDH1":  ["cg07476034", "cg09847592", "cg11014119"],
    "RASSF1": ["cg00414379", "cg01224051", "cg13878473", "cg25033530"],
    "GSTP1": ["cg02282574", "cg09356912", "cg15465165"],
    "MGMT":  ["cg02920669", "cg10463839", "cg12404465"],
    "MLH1":  ["cg00752230", "cg10866179", "cg21148842"],
    "APC":   ["cg04300418", "cg10416906", "cg23760519"],
    "DAPK1": ["cg11967187", "cg16981600", "cg21322636"],
    "ESR1":  ["cg02147208", "cg08429347", "cg15695538"],
    "PGR":   ["cg00743045", "cg13587119", "cg26875025"],
    "ERBB2": ["cg00990570", "cg08105852", "cg25241794"],
    "PTEN":  ["cg02232798", "cg10059103", "cg23506181"],
    "CDKN2A":["cg00940770", "cg12870129", "cg22561092"],
    "RUNX3": ["cg04039430", "cg09931019", "cg16977992"],
    "SOCS1": ["cg07995554", "cg12840526", "cg21572722"],
    "TIMP3": ["cg00846940", "cg08795784", "cg22938924"],
    "DCC":   ["cg02218807", "cg10467835", "cg19761115"],
    "VHL":   ["cg02534751", "cg11960478", "cg20206258"],
    "FHIT":  ["cg07231632", "cg13537839", "cg22182300"],
    "WWOX":  ["cg04517459", "cg10971782", "cg21199776"],
    "P16":   ["cg07280773", "cg14147435", "cg23578787"],
    "TMS1":  ["cg01920930", "cg11132290", "cg21910726"],
    "HIN1":  ["cg03488715", "cg12488503", "cg24315203"],
}


def flat_probe_list() -> list[str]:
    """Return every probe ID in the curated panel as a flat list."""
    return [probe for probes in CURATED_CPG_PANEL.values() for probe in probes]


def probe_to_gene() -> dict[str, str]:
    """Map each probe ID to the gene it interrogates."""
    return {probe: gene for gene, probes in CURATED_CPG_PANEL.items() for probe in probes}
