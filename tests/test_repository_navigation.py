"""Entry-page compatibility and catalogue completeness, not scientific evidence."""
from __future__ import annotations

import re
import unittest
from pathlib import Path

from tests.test_reviewer_narrative import MARKDOWN_LINK, local_target


ROOT = Path(__file__).resolve().parents[1]
# Public navigation anchors are an interface, including aliases left by older
# titles. These plain-language entry headings intentionally avoid the complex
# mathematical slug conventions of individual proof chapters.
COMPATIBILITY_ANCHORS = {
    "README.md": [
        "additional-resource-points-and-remaining-depth-gap",
        "endpoint-count-claim-and-remaining-depth-gap",
        "evidence-and-reproduction",
        "exact-and-fault-tolerant-compilation-of-hopf-differential-frames",
        "four-principal-results",
        "linear-count-endpoint-and-remaining-depth-gap",
        "one-clean-count-frontier-and-remaining-depth-gap",
        "operational-consequence",
        "target-and-resource-contract",
        "two-clean-count-theorem-and-remaining-depth-gap",
        "two-open-resource-gaps"
    ],
    "REVIEW.md": [
        "1-why-the-complete-operator-is-the-target",
        "2-result-a-exact-schedules-cover-every-clean-budget",
        "3-result-b-pay-fine-precision-once-with-sufficient-clean-work",
        "4-result-c-operator-precision-and-conditional-suffix-work",
        "4-result-c-the-optimal-one-clean-count-frontier",
        "5-result-d-count-and-depth-on-the-same-circuit",
        "6-operational-consequence-endpoint-status-and-evidence",
        "6-operational-consequence-open-gaps-and-evidence",
        "proof-roadmap-exact-and-fault-tolerant-hopf-frame-compilation"
    ],
    "WORKSPACE.md": [
        "continuing-toward-the-general-answer",
        "continuing-workspace",
        "current-checkpoint-native-source-scope-and-complete-readout",
        "current-depth-frontier",
        "destination-and-strongest-established-statements",
        "main-unfinished-route-and-exact-missing-operation",
        "mandate",
        "read-by-claim",
        "refurnishing-policy",
        "research-decision-and-stopping-rules",
        "restore-and-verify",
        "stopped-routes-reopening-conditions-and-superseded-instructions",
        "working-contract-and-reading-map"
    ],
    "docs/README.md": [
        "a-exact-compilation-at-every-clean-workspace-budget",
        "b-matching-t-count-with-sufficient-clean-workspace",
        "c-one-clean-compilation-and-its-corollaries",
        "completed-changed-decoder-supplement",
        "d-simultaneous-frame-t-count-and-t-depth",
        "d-simultaneous-real-frame-t-count-and-t-depth",
        "documentation-map",
        "evidence-and-sources",
        "fixed-decoder-qbp-consequence",
        "proof-chapters",
        "retained-two-clean-count-corollaries",
        "selected-endpoint-count-theorem-and-complex-magnitude-corollary",
        "selected-growing-source-count-theorems",
        "selected-precision-uniform-two-clean-count-theorem",
        "separate-endpoint-theorem",
        "separate-supplements-and-optional-research",
        "shared-target-and-error-contract"
    ],
    "docs/OPEN_PROBLEM.md": [
        "a-regular-fixed-tree-cayley-family",
        "bounded-diagonal-factorization",
        "constant-clean-high-precision-count",
        "constant-depth-frontier",
        "current-status",
        "earlier-endpoint-constructions",
        "endpoint-theorem-and-remaining-depth-question",
        "full-circuit-requirement",
        "large-width-t-depth",
        "middle-crossing-rank-and-the-remaining-gap",
        "one-clean-count-frontier-and-remaining-depth-question",
        "one-clean-count-frontier-and-sublinear-large-width-depth",
        "one-clean-inner-spread",
        "one-clean-midpoint-restrictions",
        "optimal-banked-count",
        "outer-axis-overlap",
        "precision-uniform-one-clean-count",
        "precision-uniform-two-clean-count",
        "proved-linear-count-families",
        "resource-frontiers-and-the-remaining-depth-question",
        "result-c-the-one-clean-count-frontier",
        "source-program-resynthesis",
        "two-clean-count-theorems-and-remaining-depth-question",
        "two-resource-questions-for-complete-hopf-frames",
        "zero-clean-error-floors-and-fixed-clean-thresholds"
    ],
    "docs/SOURCE_MAP.md": [
        "1-imported-state-preparation-toolkit",
        "2-inherited-hopf-operator-interface",
        "3-inherited-qbp-interface",
        "4-results-established-in-this-repository",
        "5-fault-tolerant-sources-and-contribution-boundaries",
        "6-evidence-classification",
        "7-provenance-records",
        "research-only-group-algebra-premise",
        "source-and-dependency-map"
    ],
    "docs/VERIFICATION.md": [
        "1-evidence-levels",
        "10-evidence-boundary",
        "11-fault-tolerant-evidence-and-approximate-qbp",
        "2-what-is-represented-locally",
        "3-geometry-and-compiler-contract",
        "4-strict-zero-construction",
        "5-binaryone-hot-prefix-decoder",
        "6-coherent-routeoperateunroute",
        "7-resource-theorem",
        "8-gradient-records",
        "9-reproduce-the-checks",
        "claim-to-evidence-index",
        "detailed-evidence-records",
        "endpoint-count-evidence",
        "one-clean-count-evidence",
        "proof-and-executable-correspondence",
        "sublinear-depth-evidence",
        "subtree-projector-decoder-evidence",
        "two-clean-count-evidence"
    ],
    "docs/reference/VERIFICATION_CATALOGUE.md": [
        "catalogue-index",
        "complete-frame-depth-checks",
        "count-source-checks",
        "detailed-finite-checks",
        "endpoint-source-and-representation-checks",
        "reading-index",
        "rendering-checks",
        "shared-source-and-recovery-checks",
        "state-based-qbp-coverage",
        "verification-catalogue"
    ],
    "research/README.md": [
        "attribution-and-reproducibility",
        "constant-clean-high-precision-count",
        "diagonal-and-reflection-architectures",
        "earlier-full-frame-reductions-and-their-native-limits",
        "endpoint-results",
        "full-frame-reductions-and-remaining-native-steps",
        "large-workspace-t-depth",
        "native-primitives-and-promised-families",
        "optimal-banked-count-frontier",
        "precision-uniform-one-clean-count",
        "precision-uniform-two-clean-count",
        "research-results-for-complete-frame-compilation",
        "retained-primitives-and-promised-families",
        "source-programs-and-encodings"
    ]
}

EVIDENCE_CATALOGUES = (
    "tests/README.md",
    "docs/reference/VERIFICATION_CATALOGUE.md",
    "docs/reference/DEPTH_FIXTURES.md",
    "docs/reference/SHARED_SOURCE_FIXTURES.md",
    "docs/reference/SOURCE_AND_OPERATOR_FIXTURES.md",
)

PAGE_LINE_BUDGETS = {
    "WORKSPACE.md": 220,
    "docs/README.md": 200,
    "docs/SOURCE_MAP.md": 300,
    "docs/OPEN_PROBLEM.md": 420,
    "docs/VERIFICATION.md": 480,
    "research/README.md": 650,
    "llms.txt": 80,
}


def entry_page_anchors(text: str) -> set[str]:
    """Collect explicit IDs and ordinary GitHub-style entry-heading slugs.

    This is deliberately not a whole-repository Markdown renderer. Explicit
    aliases protect any unusual historical title; fenced examples are ignored.
    """
    anchors: set[str] = set()
    headings: dict[str, int] = {}
    fence = None
    for line in text.splitlines():
        marker = re.match(r"^\s*(`{3,}|~{3,})", line)
        if marker:
            token = marker.group(1)
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            continue
        if fence is not None:
            continue
        anchors.update(re.findall(r"""<a\s+id=["']([^"']+)["']""", line))
        heading = re.match(r"^#{1,6}\s+(.+?)(?:\s+#+)?$", line)
        if heading is None:
            continue
        title = re.sub(r"<[^>]+>", "", heading.group(1))
        title = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", title)
        slug = re.sub(r"[^\w -]", "", title.lower()).replace(" ", "-")
        occurrence = headings.get(slug, 0)
        headings[slug] = occurrence + 1
        anchors.add(slug + (f"-{occurrence}" if occurrence else ""))
    return anchors


class RepositoryNavigationTests(unittest.TestCase):
    def test_entry_anchor_parser_ignores_examples_and_keeps_aliases(self):
        source = (
            "# Reading map\n"
            '<a id="earlier-reading-map"></a>\n'
            "## Counts and `T`-depth\n"
            "## Counts and `T`-depth\n"
            "```markdown\n"
            "# Not a heading\n"
            '<a id="not-an-anchor"></a>\n'
            "```\n"
        )
        self.assertEqual(entry_page_anchors(source), {
            "reading-map", "earlier-reading-map",
            "counts-and-t-depth", "counts-and-t-depth-1",
        })

    def test_public_entry_anchors_remain_compatible(self):
        for relative, expected in COMPATIBILITY_ANCHORS.items():
            with self.subTest(page=relative):
                actual = entry_page_anchors((ROOT / relative).read_text())
                self.assertFalse(
                    set(expected) - actual,
                    msg=f"{relative}: lost anchors {sorted(set(expected) - actual)}",
                )

    def test_research_catalogue_links_every_retained_document(self):
        catalogue = ROOT / "research/README.md"
        targets = {
            target
            for raw in MARKDOWN_LINK.findall(catalogue.read_text())
            if (target := local_target(catalogue, raw)) is not None
        }
        expected = {
            path.resolve() for path in (ROOT / "research").rglob("*.md")
            if path != catalogue
            and not any(part.startswith(".")
                        for part in path.relative_to(ROOT).parts)
        }
        missing = sorted(str(path.relative_to(ROOT)) for path in expected - targets)
        self.assertEqual(missing, [], msg="Uncatalogued research documents")

    def test_evidence_catalogues_link_every_test_module(self):
        targets = {
            target
            for relative in EVIDENCE_CATALOGUES
            for catalogue in [ROOT / relative]
            for raw in MARKDOWN_LINK.findall(catalogue.read_text())
            if (target := local_target(catalogue, raw)) is not None
        }
        expected = {
            path.resolve() for path in (ROOT / "tests").rglob("test_*.py")
            if not any(part.startswith(".")
                       for part in path.relative_to(ROOT).parts)
        }
        missing = sorted(str(path.relative_to(ROOT)) for path in expected - targets)
        self.assertEqual(missing, [], msg="Uncatalogued test modules")

    def test_entry_pages_route_instead_of_accumulating_history(self):
        for relative, limit in PAGE_LINE_BUDGETS.items():
            with self.subTest(page=relative):
                self.assertLessEqual(
                    len((ROOT / relative).read_text().splitlines()), limit,
                    msg="Move detail to its proof home or route register",
                )
        # Long single-line theorem summaries must not bypass the routing budget.
        self.assertLessEqual(len((ROOT / "llms.txt").read_text()), 8_000)


if __name__ == "__main__":
    unittest.main()
