from __future__ import annotations

from copy import deepcopy
from functools import cached_property
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Rectangle


class TreePlotter:
    """Phylogenetic Tree Plot Class"""

    def __init__(
        self,
        tree: str | Path | Tree,  # type: ignore
        *,
        format: str = "newick",
        height: float = 0.5,
        width: float = 8,
        ignore_branch_length: bool = False,
        leaf_label_size: float = 12,
        innode_label_size: float = 0,
    ):
        """
        Parameters
        ----------
        tree : str | Path | Tree
            Tree file or Tree object
        format : str, optional
            Tree format (e.g. `newick`, `phylip`, `phyloxml`, ...)
        height : float, optional
            Figure height per leaf node of tree
        width : float, optional
            Figure width
        ignore_branch_length : bool, optional
            If True, Ignore branch length for plotting tree.
        leaf_label_size : float, optional
            Leaf label size
        innode_label_size : float, optional
            Internal node label size
        """
        if isinstance(tree, (str, Path)):
            tree: Tree = Phylo.read(tree, format=format)

        # Set unique node name and branch length if not exists
        tree = self._set_unique_node_name(tree)
        max_tree_depth = max(tree.depths().values())
        if ignore_branch_length or max_tree_depth == 0:
            tree = self._to_ultrametric_tree(tree)
        self._tree = tree

        # Plot objects
        self._plot_patches: list[Patch] = []
        self._plot_funcs: list[Callable[[Axes], None]] = []

        # Plot parameters
        self._figsize = (width, height * self.tree.count_terminals())
        self._leaf_label_size = leaf_label_size
        self._innode_label_size = innode_label_size

        self._tree_line_kws: dict[str, Any] = {}
        self.update_plot_props(
            tree_line_kws=dict(color="black", lw=1, clip_on=False),
        )

    ############################################################
    # Properties
    ############################################################

    @property
    def tree(self) -> Tree:
        """BioPython's Tree Object"""
        return self._tree

    @property
    def figsize(self) -> tuple[float, float]:
        """Figure size"""
        return self._figsize

    @property
    def max_tree_depth(self) -> float:
        """Max tree depth (root-leaf max branch length)"""
        return max(self.tree.depths().values())

    @property
    def xlim(self) -> tuple[float, float]:
        """Axes xlim"""
        return (0, self.max_tree_depth)

    @property
    def ylim(self) -> tuple[float, float]:
        """Axes ylim"""
        return (0, self.tree.count_terminals() + 1)

    @property
    def name2confidence(self) -> dict[str, str | None]:
        """Tree node name & confidence dict"""
        name2confidence: dict[str, str | None] = {}
        node: Clade
        for node in self.tree.find_clades():
            if node.confidence is None:
                name2confidence[str(node.name)] = None
            else:
                name2confidence[str(node.name)] = str(node.confidence)
        return name2confidence

    @cached_property
    def name2xy(self) -> dict[str, tuple[float, float]]:
        """Tree node name & xy coordinates dict"""
        name2xy: dict[str, tuple[float, float]] = {}
        node: Clade
        for idx, node in enumerate(self.tree.get_terminals(), 1):
            # Leaf node xy coordinates
            name2xy[str(node.name)] = (self.tree.distance(node.name), idx)
        for node in self.tree.get_nonterminals("postorder"):
            # Internal node xy coordinates
            y = sum([name2xy[n.name][1] for n in node.clades]) / len(node.clades)
            name2xy[str(node.name)] = (self.tree.distance(node.name), y)
        return name2xy

    @cached_property
    def name2rect(self) -> dict[str, Rectangle]:
        """Tree node name & rectangle dict"""
        name2rect = {}
        for name, xy in self.name2xy.items():
            # Get parent node
            node: Clade = next(self.tree.find_clades(name))
            if node == self.tree.root:
                parent_node = node
            else:
                tree_path: list[Clade] = self.tree.get_path(node.name)  # type: ignore
                tree_path = [self.tree.root] + tree_path
                parent_node: Clade = tree_path[-2]  # type: ignore

            # Get child node xy coordinates
            child_node_names = [str(n.name) for n in node.find_clades()]
            x_list: list[float] = []
            y_list: list[float] = []
            for child_node_name in child_node_names:
                x, y = self.name2xy[child_node_name]
                x_list.append(x)
                y_list.append(y)

            # Calculate rectangle min-max xy coordinate
            min_x = (xy[0] + self.name2xy[str(parent_node.name)][0]) / 2
            max_x = max(x_list)
            min_y = min(y_list) - 0.5
            max_y = max(y_list) + 0.5

            # Set rectangle
            rect = Rectangle(
                xy=(min_x, min_y),
                width=max_x - min_x,
                height=max_y - min_y,
            )
            name2rect[name] = rect

        return name2rect

    ############################################################
    # Public Method
    ############################################################

    def update_plot_props(
        self,
        *,
        tree_line_kws: dict[str, Any] | None = None,
    ) -> None:
        """Update plot properties

        Parameters
        ----------
        tree_line_kws : dict[str, Any], optional
            Axes.plot properties (e.g. `dict(color="red", lw=0.5, ...)`)
        """
        tree_line_kws = {} if tree_line_kws is None else tree_line_kws
        self._tree_line_kws.update(tree_line_kws)

    def highlight(self, node_names: list[str], color: str, **kwargs) -> None:
        """Plot highlight for target nodes

        Parameters
        ----------
        node_names : list[str]
            Node names for highlight. If multiple nodes are set, MRCA node is set.
        color : str
            Highlight color
        **kwargs : dict, optional
            Rectangle properties
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Rectangle.html>
        """
        mrca_node_name = str(self.tree.common_ancestor(node_names).name)
        rect = self.name2rect[mrca_node_name]
        rect.set_color(color)
        rect.set(**kwargs)
        self._plot_patches.append(rect)

    def text(self, x: float, y: float, s: str, **kwargs) -> None:
        """Plot text

        Parameters
        ----------
        x : float
            X position
        y : float
            Y position
        s : str
            Text label
        **kwargs : dict, optional
            Text properties (e.g. `size=12, color="red"`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """

        def plot_text(ax: Axes) -> None:
            ax.text(x, y, s, **kwargs)

        self._plot_funcs.append(plot_text)

    def text_on_node(
        self,
        node_name: str | list[str] | tuple[str],
        text: str,
        *,
        size: int = 8,
        va: str = "bottom",
        ha: str = "right",
        xmargin_ratio: float = 0.01,
        ymargin_ratio: float = 0.05,
        **kwargs,
    ) -> None:
        """Plot text on target node"""
        # Get text plot target node & xy coordinate
        if isinstance(node_name, (list, tuple)):
            target_node_name = self.tree.common_ancestor(*node_name)
        else:
            target_node_name = node_name
        x, y = self.name2xy[target_node_name]
        node: Clade = next(self.tree.find_clades(target_node_name))

        def plot_text(ax: Axes) -> None:
            # Calculate x,y margin
            if node.is_terminal():
                xmargin = 0
            else:
                xmargin = max(ax.get_xlim()) * xmargin_ratio
            ymargin = max(ax.get_ylim()) / self.tree.count_terminals() * ymargin_ratio
            if va == "center":
                ymargin = 0
            if va == "top":
                ymargin *= -1

            # Plot text
            xy = (x - xmargin, y + ymargin)
            ax.text(*xy, s=text, size=size, va=va, ha=ha, **kwargs)

        self._plot_funcs.append(plot_text)

    def set_title(
        self,
        label: str,
        **kwargs,
    ) -> None:
        """Set title

        Parameters
        ----------
        label : str
            Title text label
        **kwargs : dict, optional
            Axes.set_title properties (e.g. `size=12, color="red", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_title.html>
        """

        def set_title(ax: Axes):
            ax.set_title(label, **kwargs)

        self._plot_funcs.append(set_title)

    def plotfig(
        self,
        *,
        dpi=100,
        ax: Axes | None = None,
    ) -> Figure:
        """Plot figure

        Parameters
        ----------
        dpi : int, optional
            Figure DPI
        ax : Axes | None, optional
            Matplotlib axes for plotting. If None, figure & axes are newly created.

        Returns
        -------
        figure : Figure
            Matplotlib figure
        """
        # Initialize axes
        if ax is None:
            # Create matplotlib Figure & Axes
            fig, ax = self._init_figure(self.figsize, dpi=dpi)
            self._init_axes(ax)
        else:
            # Get matplotlib Figure & Axes
            self._init_axes(ax)
            fig: Figure = ax.get_figure()  # type: ignore

        # Plot tree line
        node: Clade
        for node in self.tree.get_nonterminals():
            parent_x, parent_y = self.name2xy[str(node.name)]
            for child_node in node.clades:
                child_x, child_y = self.name2xy[child_node.name]
                v_line_points = (parent_x, parent_x), (parent_y, child_y)
                ax.plot(*v_line_points, **self._tree_line_kws)
                h_line_points = (parent_x, child_x), (child_y, child_y)
                ax.plot(*h_line_points, **self._tree_line_kws)

        # Plot node label
        for node in self.tree.find_clades():
            x, y = self.name2xy[str(node.name)]
            margin = self.max_tree_depth / 100
            if node.is_terminal():
                label_size = self._leaf_label_size
            else:
                label_size = self._innode_label_size
            if label_size > 0:
                text_kws = dict(size=label_size, ha="left", va="center")
                ax.text(x + margin, y, s=str(node.name), **text_kws)  # type: ignore

        # Plot all patches
        for patch in self._get_plot_patches():
            ax.add_patch(patch)

        # Execute all plot functions
        for plot_func in self._get_plot_funcs():
            plot_func(ax)

        return fig

    def savefig(
        self,
        savefile: str | Path,
        *,
        dpi: int = 100,
        pad_inches: float = 0.1,
    ) -> None:
        """Save figure to file

        Parameters
        ----------
        savefile : str | Path
            Save file
        dpi : int, optional
            DPI
        pad_inches : float, optional
            Padding inches
        """
        fig = self.plotfig(dpi=dpi)
        fig.savefig(
            fname=str(savefile),
            dpi=dpi,
            pad_inches=pad_inches,
            bbox_inches="tight",
        )
        # Clear & close figure to suppress memory leak
        fig.clear()
        plt.close(fig)

    ############################################################
    # Private Method
    ############################################################

    def _init_figure(
        self,
        figsize: tuple[float, float],
        dpi: int = 100,
    ) -> tuple[Figure, Axes]:
        """Initialize matplotlib figure

        Parameters
        ----------
        figsize : tuple[float, float]
            Figure size
        dpi : int, optional
            Figure DPI

        Returns
        -------
        figure, axes : tuple[Figure, Axes]
            Matplotlib Figure & Axes
        """
        fig = plt.figure(figsize=figsize, dpi=dpi, tight_layout=True)
        ax: Axes = fig.add_subplot()
        return fig, ax  # type: ignore

    def _init_axes(self, ax: Axes) -> None:
        """Initialize matplotlib axes

        - xlim = (0, `root-leaf max branch length`)
        - ylim = (0, `total tree leaf count + 1`)

        Parameters
        ----------
        ax : Axes
            Matplotlib axes
        """
        ax.set_xlim(*self.xlim)
        ax.set_ylim(*self.ylim)
        ax.axis("off")

    def _get_plot_patches(self) -> list[Patch]:
        """Plot patches"""
        return deepcopy(self._plot_patches)

    def _get_plot_funcs(self) -> list[Callable[[Axes], None]]:
        """Plot functions"""
        return self._plot_funcs

    def _set_unique_node_name(self, tree: Tree) -> Tree:
        """Set unique node name (N_1, N_2, ..., N_XXX)

        Parameters
        ----------
        tree : Tree
            Tree object

        Returns
        -------
        tree: Tree
            Unique node name set tree object
        """
        tree = deepcopy(tree)
        for idx, node in enumerate(tree.get_nonterminals(), 1):
            node.name = f"N_{idx}" if node.name is None else node.name
        return tree

    def _to_ultrametric_tree(self, tree: Tree) -> Tree:
        """Convert to ultrametric tree

        Parameters
        ----------
        tree : Tree
            Tree

        Returns
        -------
        tree : Tree
            Ultrametric tree
        """
        tree = deepcopy(tree)
        # Get unit branch depth info
        name2depth = {n.name: d for n, d in tree.depths(True).items()}
        name2depth = dict(sorted(name2depth.items(), key=lambda t: t[1], reverse=True))
        max_tree_depth = max(name2depth.values())
        # Reset node branch length
        for node in tree.find_clades():
            node.branch_length = None
        tree.root.branch_length = 0
        # Calculate appropriate ultrametric tree branch length
        for name, depth in name2depth.items():
            node = next(tree.find_clades(name))
            if not node.is_terminal():
                continue
            path: list[Clade] | None = tree.get_path(node)
            if path is None:
                raise ValueError(f"{name=} node not exists?")
            if depth == max_tree_depth:
                for path_node in path:
                    path_node.branch_length = 1
            else:
                # Collect nodes info which has branch length
                bl_sum, bl_exist_node_count = 0, 0
                for path_node in path:
                    if path_node.branch_length is not None:
                        bl_sum += path_node.branch_length
                        bl_exist_node_count += 1
                # Set branch length to no branch length nodes
                other_bl = (max_tree_depth - bl_sum) / (len(path) - bl_exist_node_count)
                for path_node in path:
                    if path_node.branch_length is None:
                        path_node.branch_length = other_bl
        return tree
