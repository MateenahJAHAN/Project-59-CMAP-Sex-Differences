"""Study design."""


import p59style as S
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

PANELS = (
    (
        "A",
        "Audit the recorded labels against measured expression",
        (
            ("Vehicle cells", "Tahoe-100M, EmeraldBay\nand sci-Plex"),
            ("Score Y genes and XIST", "Summarize by cell line\nand classify as Y only, XIST only,\nboth or neither"),
            ("Compare the patterns", "Recorded donor sex, external\nbulk RNA-seq and STR AMELY\nstatus for the Y call"),
        ),
    ),
    (
        "B",
        "Compare LINCS2020 drug responses by expression pattern",
        (
            ("Two inputs", "Bulk RNA-seq for marker calls\nLINCS2020 L1000 for response"),
            ("Compare within lineage", "Y only and XIST only each\ncompared with neither marker\nwith adjustment for TP53"),
            ("Estimate differences", "Panel and per-compound tests\nPooled vehicle-referenced\nresponse-size estimates"),
        ),
    ),
)


def build():
    width, height = S.FULL, 2.44
    fig = S.plt.figure(figsize=(width, height))
    coord = fig.dpi_scale_trans
    margin, gap, box_height = 0.08, 0.22, 0.80
    box_width = (width - 2 * margin - 2 * gap) / 3
    cards = []
    for top, (letter, heading, boxes) in zip((2.36, 1.17), PANELS):
        fig.text(margin, top, letter, transform=coord, va="top", fontsize=10, weight="bold", color=S.TEXT)
        fig.text(margin + 0.22, top, heading, transform=coord, va="top",
                 fontsize=8.5, weight="bold", color=S.TEXT)
        y = top - 0.26 - box_height
        for col, (title, detail) in enumerate(boxes):
            x = margin + col * (box_width + gap)
            patch = FancyBboxPatch((x, y), box_width, box_height, transform=coord,
                                   boxstyle="round,pad=0,rounding_size=0.06",
                                   facecolor=S.BAND, edgecolor="none")
            fig.add_artist(patch)
            title_artist = fig.text(x + 0.10, y + box_height - 0.11, title, transform=coord, va="top",
                                    fontsize=8, weight="bold", color=S.TEXT)
            detail_artist = fig.text(x + 0.10, y + box_height - 0.31, detail, transform=coord, va="top",
                                     fontsize=7.6, linespacing=1.23, color=S.TEXT2)
            cards.append((patch, title_artist, detail_artist))
            if col < len(boxes) - 1:
                fig.add_artist(FancyArrowPatch(
                    (x + box_width + 0.035, y + box_height / 2),
                    (x + box_width + gap - 0.035, y + box_height / 2),
                    transform=coord, arrowstyle="-|>", mutation_scale=7,
                    linewidth=0.75, color=S.TEXT2, shrinkA=0, shrinkB=0,
                ))
    fig.canvas.draw()
    return fig, cards


def save(fig, stem):
    # Preserve the explicit full-width frame with an inset around all artists.
    with S.plt.rc_context({"savefig.bbox": None}):
        S.save(fig, str(stem))
