from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree


class CafeParser:
    """CAFE5 Output File Parser Class"""

    def __init__(self, target_dir: str | Path):
        """
        Parameters
        ----------
        target_dir : str | Path
            Parse target CAFE5 result directory
        """
        target_dir = Path(target_dir)
        self._asr_tree_file = next(target_dir.glob("*_asr.tre"))
        self._branch_prob_file = next(target_dir.glob("*_branch_probabilities.tab"))
        self._clade_result_file = next(target_dir.glob("*_clade_results.txt"))
        self._family_change_file = next(target_dir.glob("*_change.tab"))
        self._family_count_file = next(target_dir.glob("*_count.tab"))
        self._family_result_file = next(target_dir.glob("*_family_results.txt"))

        # Check duplicate family id
        self.check_family_id_duplication()

        # Parse CAFE5 result files
        _ = self.famid2asr_tree
        _ = self.famid2branch_prob
        _ = self.famid2family_change
        _ = self.famid2family_count
        _ = self.famid2family_result

    ############################################################
    # CAFE Ouput File Path
    ############################################################

    @property
    def asr_tree_file(self) -> Path:
        """Ancestral state reconstruction tree file

        `{model}_asr.tre`
        """
        return self._asr_tree_file

    @property
    def branch_prob_file(self) -> Path:
        """Branch probability file

        `{model}_branch_probabilities.tab`
        """
        return self._branch_prob_file

    @property
    def clade_result_file(self) -> Path:
        """All families clade increase/decrease results file

        `{model}_clade_results.txt`
        """
        return self._clade_result_file

    @property
    def family_change_file(self) -> Path:
        """Gene family change file

        `{model}_change.tab`
        """
        return self._family_change_file

    @property
    def family_count_file(self) -> Path:
        """Gene family count file

        `{model}_count.tab`
        """
        return self._family_count_file

    @property
    def family_result_file(self) -> Path:
        """Gene family significance result file

        `{model}_family_results.txt`
        """
        return self._family_result_file

    ############################################################
    # Parse Dictionaly Results
    ############################################################

    @cached_property
    def famid2asr_tree(self) -> dict[str, Tree]:
        """Gene family id & tree dict"""
        return {str(t.name): t for t in self.asr_trees}

    @cached_property
    def famid2branch_prob(self) -> dict[str, BranchProb]:
        """Gene family id & branch probability dict"""
        return {bp.famid: bp for bp in self.branch_probs}

    @cached_property
    def famid2family_change(self) -> dict[str, FamilyChange]:
        """Gene family id & change dict"""
        return {fc.famid: fc for fc in self.family_changes}

    @cached_property
    def famid2family_count(self) -> dict[str, FamilyCount]:
        """Gene family id & count dict"""
        return {fc.famid: fc for fc in self.family_counts}

    @cached_property
    def famid2family_result(self) -> dict[str, FamilyResult]:
        """Gene family id & significance result dict"""
        return {fr.famid: fr for fr in self.family_results}

    ############################################################
    # Parse List Results
    ############################################################

    @cached_property
    def asr_trees(self) -> list[Tree]:
        """Parse ancestral state reconstruction tree file"""
        return CafeTree(self.asr_tree_file).trees

    @cached_property
    def branch_probs(self) -> list[BranchProb]:
        """Parse branch probability file"""
        return BranchProb.parse(self.branch_prob_file)

    @cached_property
    def clade_results(self) -> list[CladeResult]:
        """Parse all families clade increase/decrease results file"""
        return CladeResult.parse(self.clade_result_file)

    @cached_property
    def family_changes(self) -> list[FamilyChange]:
        """Parse gene family change file"""
        return FamilyChange.parse(self.family_change_file)

    @cached_property
    def family_counts(self) -> list[FamilyCount]:
        """Parse gene family count file"""
        return FamilyCount.parse(self.family_count_file)

    @cached_property
    def family_results(self) -> list[FamilyResult]:
        """Parse gene family significance result file"""
        return FamilyResult.parse(self.family_result_file)

    ############################################################
    # Summary Results
    ############################################################

    @cached_property
    def all_famid_list(self) -> list[str]:
        """All gene family id"""
        return [fr.famid for fr in self.family_results]

    @cached_property
    def signif_famid_list(self) -> list[str]:
        """Significant gene family id"""
        return [fr.famid for fr in self.family_results if fr.is_signif]

    @cached_property
    def taxonid_list(self) -> list[str]:
        """Taxonid list"""
        return [cr.taxonid for cr in self.clade_results]

    def write_result_summary(self, outfile: str | Path) -> None:
        """Write result summary

        Table of gene counts, changes, and expansion/contraction p-values
        for each family and taxon.

        Parameters
        ----------
        outfile : str | Path
            Result summary file
        """
        result_summary_content = "FamilyID\tTaxonID\tCount\tChange\tPvalue\n"
        for famid in self.signif_famid_list:
            taxonid2count = self.famid2family_count[famid].taxonid2count
            taxonid2change = self.famid2family_change[famid].taxonid2change
            taxonid2prob = self.famid2branch_prob[famid].taxonid2prob

            for taxonid in self.taxonid_list:
                row_values = []
                row_values.append(famid)
                row_values.append(taxonid)
                row_values.append(str(taxonid2count[taxonid]))
                row_values.append(str(taxonid2change[taxonid]))
                row_values.append(str(taxonid2prob[taxonid]))

                result_summary_content += "\t".join(row_values) + "\n"

        with open(outfile, "w") as f:
            f.write(result_summary_content)

    def get_taxonid2total_increase(self, only_signif: bool = False) -> dict[str, int]:
        """Get taxonid & total gene family increase dict

        This result is same as `{model}_clade_result.txt` when only_signif=False.
        I implemented this method to check and understand CAFE5 result for details.

        Parameters
        ----------
        only_signif : bool, optional
            If True, count only significant increase gene family

        Returns
        -------
        taxonid2total_increase : dict[str, int]
            taxonid & total gene family increase dict
        """
        taxonid2total_increase = defaultdict(int)
        famid_list = self.signif_famid_list if only_signif else self.all_famid_list
        for famid in famid_list:
            family_change = self.famid2family_change[famid]
            for taxonid, change in family_change.taxonid2change.items():
                if change > 0:
                    taxonid2total_increase[taxonid] += 1
        return defaultdict(int, sorted(taxonid2total_increase.items()))

    def get_taxonid2total_decrease(self, only_signif: bool = False) -> dict[str, int]:
        """Get taxonid & total gene family decrease dict

        This result is same as `{model}_clade_result.txt` when only_signif=False.
        I implemented this method to check and understand CAFE5 result for details.

        Parameters
        ----------
        only_signif : bool, optional
            If True, count only significant decrease gene family

        Returns
        -------
        taxonid2total_decrease : dict[str, int]
            taxonid & total gene family decrease dict
        """
        taxonid2total_decrease = defaultdict(int)
        famid_list = self.signif_famid_list if only_signif else self.all_famid_list
        for famid in famid_list:
            family_change = self.famid2family_change[famid]
            for taxonid, change in family_change.taxonid2change.items():
                if change < 0:
                    taxonid2total_decrease[taxonid] += -1
        return defaultdict(int, sorted(taxonid2total_decrease.items()))

    def check_family_id_duplication(self) -> None:
        """Check family id duplication

        Raises
        ------
        ValueError
            If Family IDs are duplicated, raise ValueError
        """
        dup_famid_list = []
        counter = Counter(self.all_famid_list)
        for famid, count in counter.items():
            if count >= 2:
                dup_famid_list.append(famid)
        if len(dup_famid_list) > 0:
            raise ValueError(f"Following Family IDs are duplicated!!\n{dup_famid_list}")


class ParserBase:
    """Parser Base Class"""

    def rename_taxonid(self, taxonid: str) -> str:
        """Rename taxonid to simple

        Parameters
        ----------
        taxonid : str
            Taxon id (e.g. `human<0>`, `<21>`)

        Returns
        -------
        rename_taxonid : str
            Renamed taxonid (e.g. `human`, `21`)
        """
        if re.search(r"^<\d+>$", taxonid):
            # Rename "<2>" to "2"
            return taxonid.replace("<", "").replace(">", "")
        else:
            # Rename "human<0>" to "human"
            return re.sub(r"<\d+>", "", taxonid)


class CafeTree(ParserBase):
    def __init__(self, asr_tree_file: str | Path):
        trees: list[Tree] = list(Phylo.parse(asr_tree_file, "nexus"))
        for tree in trees:
            for node in tree.find_clades():
                node_name = str(node.name).replace("*", "")
                taxonid = "_".join(node_name.split("_")[0:-1])
                gene_count = node_name.split("_")[-1]
                node.name = self.rename_taxonid(taxonid)
                node.confidence = gene_count
        self._trees = trees

    @property
    def trees(self) -> list[Tree]:
        """CAFE asr trees"""
        return self._trees


@dataclass
class BranchProb(ParserBase):
    """Branch Probability Dataclass"""

    famid: str
    taxonid2prob: dict[str, float | None]

    def __post_init__(self):
        # Rename taxonid
        rename_taxonid2prob = {}
        for taxonid, prob in self.taxonid2prob.items():
            rename_taxonid2prob[self.rename_taxonid(taxonid)] = prob
        self.taxonid2prob = rename_taxonid2prob

    @staticmethod
    def parse(branch_prob_file: str | Path) -> list[BranchProb]:
        """Parse branch probabilities file

        Parameters
        ----------
        branch_prob_file : str | Path
            Branch probabilities file (`{model}_branch_probabilities.tab`)

        Returns
        -------
        branch_probs : list[BranchProb]
            Family Branch probabilities
        """
        with open(branch_prob_file) as f:
            reader = csv.reader(f, delimiter="\t")
            taxonid_list = next(reader)[1:]
            branch_probs = []
            for row in reader:
                famid, probs = row[0], row[1:]
                probs = [float(p) if p != "N/A" else None for p in probs]
                taxonid2prob = dict(zip(taxonid_list, probs))
                branch_probs.append(BranchProb(famid, taxonid2prob))
            return branch_probs


@dataclass
class CladeResult(ParserBase):
    """All families Clade Increase/Decrease Result Dataclass"""

    taxonid: str
    increase: int
    decrease: int

    def __post_init__(self):
        # Rename taxonid
        self.taxonid = self.rename_taxonid(self.taxonid)

    @property
    def increase_label(self) -> str:
        """Increase label (e.g. `+100`)"""
        return f"+{self.increase}"

    @property
    def decrease_label(self) -> str:
        """Decrease label (e.g. `-100`)"""
        return f"-{self.decrease}"

    @staticmethod
    def parse(clade_result_file: str | Path) -> list[CladeResult]:
        """Parse all families clade increase/decrease results file

        Parameters
        ----------
        clade_result_file : str | Path
            Clade results file (`{model}_clade_results.txt`)

        Returns
        -------
        clade_results : list[CladeResult]
            Clade results
        """
        clade_results = []
        with open(clade_result_file) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # Skip header
            for row in reader:
                taxonid, increase, decrease = str(row[0]), int(row[1]), int(row[2])
                clade_result = CladeResult(taxonid, increase, decrease)
                clade_results.append(clade_result)
        return clade_results


@dataclass
class FamilyChange(ParserBase):
    """Gene Family Change Dataclass"""

    famid: str
    taxonid2change: dict[str, int]

    def __post_init__(self):
        # Rename taxonid
        rename_taxonid2change = {}
        for taxonid, change in self.taxonid2change.items():
            rename_taxonid2change[self.rename_taxonid(taxonid)] = change
        self.taxonid2change = rename_taxonid2change

    @staticmethod
    def parse(change_file: str | Path) -> list[FamilyChange]:
        """Parse gene family change file

        Parameters
        ----------
        change_file : str | Path
            Gene family change file (`{model}_change.tab`)

        Returns
        -------
        family_changes : list[FamilyChange]
            Family changes
        """
        family_changes = []
        with open(change_file) as f:
            reader = csv.reader(f, delimiter="\t")
            taxonid_list = next(reader)[1:]
            for row in reader:
                famid, change = row[0], list(map(int, row[1:]))
                taxonid2change = dict(zip(taxonid_list, change))
                family_changes.append(FamilyChange(famid, taxonid2change))
            return family_changes


@dataclass
class FamilyCount(ParserBase):
    """Gene Family count Dataclass"""

    famid: str
    taxonid2count: dict[str, int]

    def __post_init__(self):
        # Rename taxonid
        rename_taxonid2count = {}
        for taxonid, count in self.taxonid2count.items():
            rename_taxonid2count[self.rename_taxonid(taxonid)] = count
        self.taxonid2count = rename_taxonid2count

    @staticmethod
    def parse(count_file: str | Path) -> list[FamilyCount]:
        """Parse gene family count file

        Parameters
        ----------
        count_file : str | Path
            Gene family count file (`{model}_count.tab`)

        Returns
        -------
        family_counts : list[FamilyCount]
            Family counts
        """
        family_counts = []
        with open(count_file) as f:
            reader = csv.reader(f, delimiter="\t")
            taxonid_list = next(reader)[1:]
            for row in reader:
                famid, count = row[0], list(map(int, row[1:]))
                taxonid2count = dict(zip(taxonid_list, count))
                family_counts.append(FamilyCount(famid, taxonid2count))
            return family_counts


@dataclass
class FamilyResult:
    """Gene Family Significance Result DataClass"""

    famid: str
    pvalue: float
    is_signif: bool

    @property
    def display_pvalue(self) -> str:
        """pvalue string for display"""
        if self.pvalue == 0:
            return "0"
        elif self.pvalue >= 0.001:
            return f"{self.pvalue:.3f}"
        else:
            return f"{self.pvalue:.1e}"

    @staticmethod
    def parse(family_result_file: str | Path) -> list[FamilyResult]:
        """Parse gene family significance result file

        Parameters
        ----------
        family_result_file : str | Path
            Gene family significance result file (`{model}_family_results.tab`)

        Returns
        -------
        family_results : list[FamilyCount]
            Family results
        """
        family_results = []
        with open(family_result_file) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)[1:]  # Skip header
            for row in reader:
                famid, pvalue, signif = row[0], float(row[1]), row[2]
                is_signif = True if signif == "y" else False
                family_results.append(FamilyResult(famid, pvalue, is_signif))
            return family_results
