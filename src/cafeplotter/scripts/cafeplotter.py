from __future__ import annotations

import argparse
from pathlib import Path

import cafeplotter
from cafeplotter import CafeParser, TreePlotter

# P-value thresholds
HIGH_P, MID_P, LOW_P = 0.05, 0.01, 0.001


def main():
    """Main function called from CLI"""
    args = get_args()
    run(**args.__dict__)


def run(
    indir: Path,
    outdir: Path,
    format: str = "png",
    fig_height: float = 0.5,
    fig_width: float = 8.0,
    leaf_label_size: int = 12,
    count_label_size: int = 8,
    innode_label_size: int = 0,
    p_label_size: int = 0,
    ignore_branch_length: bool = False,
    expansion_color: str = "red",
    contraction_color: str = "blue",
    dpi: int = 300,
) -> None:
    """Plot CAFE5 gene family expansion/contraction result"""
    # Parse CAFE5 output result
    print("\nParsing CAFE5 result files...\n")
    cp = CafeParser(indir)
    outdir.mkdir(exist_ok=True)

    # Output CAFE5 gene expansion/contraction result summary for each family and taxon.
    cafe_result_summary_file = outdir / "result_summary.tsv"
    print("Write expansion/contraction result summary for each family and taxon:")
    print(f" => `{cafe_result_summary_file}`\n")
    cp.write_result_summary(cafe_result_summary_file)

    # Set general gene family tree plot properties
    tree_plot_props: dict = dict(
        height=fig_height,
        width=fig_width,
        leaf_label_size=leaf_label_size,
        innode_label_size=innode_label_size,
        ignore_branch_length=ignore_branch_length,
    )

    # Plot all gene family expansion/contraction tree
    all_gene_family_plot_file = outdir / f"summary_all_gene_family.{format}"
    print("Plot Summary of All Expansion/Contraction Gene Family Tree:")
    print(f" => `{all_gene_family_plot_file}`\n")
    tp = TreePlotter(cp.asr_trees[0], **tree_plot_props)
    for clade_result in cp.clade_results:
        label = f"{clade_result.increase_label} / {clade_result.decrease_label}"
        tp.text_on_node(clade_result.taxonid, label, size=count_label_size)
    tp.set_title("Summary of All Expansion/Contraction Gene Family")
    tp.savefig(all_gene_family_plot_file, dpi=dpi)

    # Plot each gene family expansion/contraction tree
    gene_family_outdir = outdir / "gene_family"
    gene_family_outdir.mkdir(exist_ok=True)
    for famid, asr_tree in cp.famid2asr_tree.items():
        # Ignore not significant increase/decrease gene family
        if famid not in cp.signif_famid_list:
            continue
        # Print console
        gene_family_plot_file = gene_family_outdir / f"{famid}_gene_family.{format}"
        print(f"[FamilyID={famid}] ", end="")
        print("Plot Significant Expansion/Contraction Gene Family Tree:")
        print(f" => `{gene_family_plot_file}`")

        # Plot gene family expansion/contraction tree
        tp = TreePlotter(asr_tree, **tree_plot_props)
        taxonid2count = cp.famid2family_count[famid].taxonid2count
        taxonid2change = cp.famid2family_change[famid].taxonid2change
        taxonid2prob = cp.famid2branch_prob[famid].taxonid2prob
        for taxonid in taxonid2count.keys():
            count = taxonid2count[taxonid]
            change = taxonid2change[taxonid]
            prob = taxonid2prob[taxonid]
            # Plot gene count label
            gene_count_label = make_gene_count_label(count, prob)
            gene_count_color = "black"
            if gene_count_label.startswith("*"):
                gene_count_color = expansion_color if change >= 0 else contraction_color
            tp.text_on_node(
                taxonid, gene_count_label, size=count_label_size, color=gene_count_color
            )
            # Plot branch probability
            if p_label_size > 0 and prob is not None:
                prob_label = f"{prob:.3f}" if prob > 0.001 else f"{prob:.1e}"
                tp.text_on_node(taxonid, prob_label, size=p_label_size, va="top")

        # Plot title & text
        pvalue = cp.famid2family_result[famid].display_pvalue
        title = "Significant Expansion/Contraction Gene Family\n"
        title += f"(FamilyID={famid}, p-value={pvalue})"
        tp.set_title(title)
        signif_label = f"(*) p < {HIGH_P}  (**) p < {MID_P}  (***) p < {LOW_P}"
        tp.text(
            x=max(tp.xlim), y=min(tp.ylim), s=signif_label, size=7, va="top", ha="right"
        )
        tp.savefig(gene_family_plot_file, dpi=dpi)


def make_gene_count_label(gene_count: int, branch_prob: float | None) -> str:
    """Make gene count label with branch probability

    - branch_prob < 0.001: `***{count}`
    - branch_prob < 0.01: `**{count}`
    - branch_prob < 0.05: `*{count}`
    - branch_prob >= 0.05: `{count}`

    Parameters
    ----------
    gene_count : int
        Gene count
    branch_prob : float | None
        Branch probability

    Returns
    -------
    label : str
        Gene count label
    """
    label = str(gene_count)
    if branch_prob is None:
        return label
    if branch_prob < LOW_P:
        return f"***{label}"
    elif branch_prob < MID_P:
        return f"**{label}"
    elif branch_prob < HIGH_P:
        return f"*{label}"
    else:
        return label


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """

    class CustomHelpFormatter(argparse.RawTextHelpFormatter):
        def __init__(self, prog, indent_increment=2, max_help_position=40, width=None):
            super().__init__(prog, indent_increment, max_help_position, width)

    desc = "A tool for plotting CAFE5 gene family expansion/contraction result"
    parser = argparse.ArgumentParser(
        description=desc, add_help=False, formatter_class=CustomHelpFormatter
    )

    #######################################################
    # General options
    #######################################################
    general_opts = parser.add_argument_group("General Options")
    general_opts.add_argument(
        "-i",
        "--indir",
        type=Path,
        help="CAFE5 result directory as input",
        required=True,
        metavar="IN",
    )
    general_opts.add_argument(
        "-o",
        "--outdir",
        type=Path,
        help="Output directory for plotting CAFE5 result",
        required=True,
        metavar="OUT",
    )
    default_format = "png"
    format_list = ["png", "jpg", "pdf", "svg"]
    general_opts.add_argument(
        "--format",
        type=str,
        help="Output image format ('png'[default]|'jpg'|'svg'|'pdf')",
        default=default_format,
        choices=format_list,
        metavar="",
    )
    general_opts.add_argument(
        "-v",
        "--version",
        version=f"v{cafeplotter.__version__}",
        help="Print version information",
        action="version",
    )
    general_opts.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )

    #######################################################
    # Figure appearence options
    #######################################################
    fig_opts = parser.add_argument_group("Figure Appearence Options")
    default_height = 0.5
    fig_opts.add_argument(
        "--fig_height",
        type=float,
        help=f"Figure height per leaf node of tree (Default: {default_height})",
        default=default_height,
        metavar="",
    )
    default_width = 8.0
    fig_opts.add_argument(
        "--fig_width",
        type=float,
        help=f"Figure width (Default: {default_width})",
        default=default_width,
        metavar="",
    )
    default_leaf_label_size = 12
    fig_opts.add_argument(
        "--leaf_label_size",
        type=int,
        help=f"Leaf label size (Default: {default_leaf_label_size})",
        default=default_leaf_label_size,
        metavar="",
    )
    default_count_label_size = 8
    fig_opts.add_argument(
        "--count_label_size",
        type=int,
        help=f"Gene count label size (Default: {default_count_label_size})",
        default=default_count_label_size,
        metavar="",
    )
    default_innode_label_size = 0
    fig_opts.add_argument(
        "--innode_label_size",
        type=int,
        help=f"Internal node label size (Default: {default_innode_label_size})",
        default=default_innode_label_size,
        metavar="",
    )
    default_p_label_size = 0
    fig_opts.add_argument(
        "--p_label_size",
        type=int,
        help=f"Branch p-value label size (Default: {default_p_label_size})",
        default=default_p_label_size,
        metavar="",
    )
    fig_opts.add_argument(
        "--ignore_branch_length",
        help="Ignore branch length for plotting tree (Default: OFF)",
        action="store_true",
    )
    default_expansion_color = "red"
    fig_opts.add_argument(
        "--expansion_color",
        type=str,
        help="Plot color of gene family expansion "
        f"(Default: '{default_expansion_color}')",
        default=default_expansion_color,
        metavar="",
    )
    default_contraction_color = "blue"
    fig_opts.add_argument(
        "--contraction_color",
        type=str,
        help="Plot color of gene family contraction "
        f"(Default: '{default_contraction_color}')",
        default=default_contraction_color,
        metavar="",
    )
    default_dpi = 300
    fig_opts.add_argument(
        "--dpi",
        type=int,
        help=f"Figure DPI (Default: {default_dpi})",
        default=default_dpi,
        metavar="",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
