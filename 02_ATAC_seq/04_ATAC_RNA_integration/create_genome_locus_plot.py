# import required packages
from pathlib import Path
import sys
import subprocess
import shutil

import numpy as np
import pandas as pd
import pyBigWig
import matplotlib.pyplot as plt

# project directory
proj_dir = Path("/Volumes/ritd-ag-project-rd0227-nhhan20/Michelle/ATAC_seq/PRJ26157")

# input directories
bw_dir = proj_dir / "bigwig"
ferris_dir = proj_dir / "ferris_2025"

# output directories
plot_dir = proj_dir / "plots"
locus_dir = proj_dir / "genome_tracks"
chip_dir = ferris_dir / "H3k27ac_ChIP"

plot_dir.mkdir(parents=True, exist_ok=True)
locus_dir.mkdir(parents=True, exist_ok=True)

# verify input directories
for name, directory in {
    "Project": proj_dir,
    "ATAC bigWig": bw_dir,
    "Ferris reference": ferris_dir,
    "H3K27ac ChIP-seq": chip_dir
}.items():
    if not directory.is_dir():
        raise FileNotFoundError(f"{name} directory not found: {directory}")

# find ATAC-seq bigWig files
atac_bw_files = sorted(bw_dir.glob("*.bw"))

if not atac_bw_files:
    raise FileNotFoundError(f"No ATAC-seq bigWig files found in {bw_dir}")

# create sample metadata
atac_metadata = pd.DataFrame({
    "file": atac_bw_files,
    "sample": [f.name.split("-")[0] for f in atac_bw_files],
    "sex": ["Male" if "-m-" in f.name else "Female" for f in atac_bw_files],
    "condition": ["60 min" if "-60_" in f.name else "Basal" for f in atac_bw_files]
})

# exclude samples that failed QC
failed_samples = ["2095", "2124", "2127", "2131"]

atac_metadata = atac_metadata[~atac_metadata["sample"].isin(failed_samples)].copy()

# check sample metadata
print(atac_metadata.to_string(index=False))
print("\nSamples per group:")
print(atac_metadata.groupby(["sex", "condition"]).size())

# define paths for genomic annotation and reference files
gtf_file = Path("/Volumes/ritd-ag-project-rd0227-nhhan20/Michelle/genomes/mm39/gencode.vM36.basic.annotation.gtf")
rna_atac_plac_file = proj_dir / "csv" / "rna_atac_plac_all.csv"
plac_file = ferris_dir / "PLAC" / "plac_interactions_mm39.bedpe"
plac_bins_file = ferris_dir / "PLAC" / "plac_bin_mapping_mm39.csv"
h3k27ac_file = chip_dir / "h3k27ac_reference_mm39.bed"

# check that all files exist
annotation_files = {
    "Gene annotation": gtf_file,
    "DAR-gene associations": rna_atac_plac_file,
    "PLAC interactions": plac_file,
    "PLAC distal bins": plac_bins_file,
    "H3K27ac peaks": h3k27ac_file
}

for name, path in annotation_files.items():
    if not path.is_file():
        raise FileNotFoundError(f"{name} not found: {path}")

print("All annotation and reference files found.")

# import genomic annotation files
rna_atac_plac = pd.read_csv(rna_atac_plac_file)

plac_interactions = pd.read_csv(plac_file,
                                sep="\t",
                                header=None,
                                names=["chr1", "start1", "end1", "chr2", "start2", "end2"])

plac_bins = pd.read_csv(plac_bins_file)

h3k27ac_peaks = pd.read_csv(h3k27ac_file,
                            sep="\t",
                            header=None,
                            names=["chr", "start", "end"])

# check imported datasets
print("Peak-gene associations:", len(rna_atac_plac))
print("PLAC distal bins:", len(plac_bins))
print("PLAC interactions:", len(plac_interactions))
print("H3K27ac peaks:", len(h3k27ac_peaks))

# define a function to retrieve gene-associated DARs
def get_gene_peaks(gene, significant_only=False):
    gene_peaks = rna_atac_plac.loc[rna_atac_plac["symbol"] == gene].copy()

    if significant_only:
        gene_peaks = gene_peaks.loc[gene_peaks["atac_padj"] < 0.05]

    return gene_peaks

# import gene-level annotations from GENCODE GTF
genes = pd.read_csv(gtf_file,
                    sep="\t",
                    comment="#",
                    header=None,
                    names=["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])

# retain gene-level annotations
genes = genes.loc[genes["feature"] == "gene"].copy()

# extract gene symbols and Ensembl IDs
genes["symbol"] = genes["attributes"].str.extract(r'gene_name "([^"]+)"')
genes["ensembl"] = genes["attributes"].str.extract(r'gene_id "([^"]+)"')

# remove Ensembl version suffixes
genes["ensembl"] = genes["ensembl"].str.replace(r"\.\d+$", "", regex=True)

# determine transcription start site
genes["tss"] = np.where(genes["strand"] == "+",
                        genes["start"],
                        genes["end"])

print("Number of annotated genes:", len(genes))

# define a function to retrieve gene annotation
def get_gene_annotation(gene):
    return genes.loc[genes["symbol"] == gene].copy()

# define a function to determine the plotting window
def get_gene_window(gene, padding=20000, mode="all"):

    # retrieve gene annotation
    annotation = get_gene_annotation(gene)

    if annotation.empty:
        raise ValueError(f"Gene {gene} not found in GENCODE annotation")

    # retrieve linked peaks according to selected mode
    if mode not in ("all", "significant"):
        raise ValueError("mode must be 'all' or 'significant'")

    peaks = get_gene_peaks(gene,
                           significant_only=(mode == "significant"))

    # start with gene boundaries
    chrom = annotation["chr"].iloc[0]
    start = int(annotation["start"].min())
    end = int(annotation["end"].max())

    # expand window to include linked peaks
    peaks = peaks.loc[peaks["chr"] == chrom]

    if not peaks.empty:
        start = min(start, int(peaks["start"].min()))
        end = max(end, int(peaks["end"].max()))

    # add flanking sequence
    start = max(0, start - padding)
    end = end + padding

    return chrom, start, end

# define a function to extract mean bigWig signal across genomic bins
def get_bw_signal(bw_file, chrom, start, end, n_bins=2000):

    with pyBigWig.open(str(bw_file)) as bw:
        signal = bw.stats(chrom,
                          start,
                          end,
                          type="mean",
                          nBins=n_bins,
                          exact=True)

    return np.array([0 if value is None else value for value in signal])

# define a function to calculate mean ATAC-seq coverage across biological replicates
def get_group_signal(sex, condition, chrom, start, end, n_bins=2000):

    # select samples belonging to the group
    samples = atac_metadata.loc[(atac_metadata["sex"] == sex) & (atac_metadata["condition"] == condition)]

    if samples.empty:
        raise ValueError(f"No ATAC-seq samples found for {sex} {condition}")

    # extract coverage from each biological replicate
    signals = [get_bw_signal(file, chrom, start, end, n_bins) for file in samples["file"]]

    # calculate mean coverage across samples
    mean_signal = np.mean(signals, axis=0)

    return mean_signal

# define a function to export group-mean ATAC-seq coverage as bedGraph
def export_group_bedgraph(gene, sex, condition, bin_size=25, mode="significant"):

    # retrieve genomic window
    chrom, start, end = get_gene_window(gene, mode=mode)

    # calculate number of bins
    n_bins = int(np.ceil((end - start) / bin_size))

    # calculate mean coverage across biological replicates
    signal = get_group_signal(sex, condition, chrom, start, end, n_bins=n_bins)

    # create output directory
    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)

    # define output filename
    condition_name = condition.replace(" ", "")
    output_file = output_dir / f"{gene}_{sex}_{condition_name}.bedgraph"

    # match genomic intervals to the bins returned by pyBigWig
    edges = np.linspace(start, end, n_bins + 1)

    with open(output_file, "w") as f:
        for i, value in enumerate(signal):
            bin_start = int(np.floor(edges[i]))
            bin_end = int(np.floor(edges[i + 1]))

            if bin_end > bin_start:
                f.write(f"{chrom}\t{bin_start}\t{bin_end}\t{value:.6f}\n")

    return output_file

# define a function to export a gene-specific GTF
def export_region_gtf(gene):

    annotation = get_gene_annotation(gene)

    if annotation.empty:
        raise ValueError(f"Gene {gene} not found in GENCODE annotation")

    gene_id = annotation["ensembl"].iloc[0]

    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{gene}_region.gtf"

    with open(gtf_file, "r") as source, open(output_file, "w") as output:
        for line in source:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            attributes = fields[8]

            if (f'gene_id "{gene_id}.' in attributes
                or f'gene_id "{gene_id}"' in attributes):
                output.write(line)
        
    # verify that the exported GTF contains annotations
    if output_file.stat().st_size == 0:
        raise ValueError(f"No GTF annotations exported for {gene}")

    return output_file

# define a function to export significant DARs for highlighting
def export_dar_highlights(gene, padding=1000):

    # retrieve significant DARs across both sexes
    peaks = (get_gene_peaks(gene, significant_only=True)
             .drop_duplicates(subset=["chr", "start", "end"])
             .sort_values(["chr", "start"])
             .copy())

    # add visual padding around each DAR
    peaks["start"] = (peaks["start"] - padding).clip(lower=0)
    peaks["end"] = peaks["end"] + padding

    # define output file
    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{gene}_DAR_highlights.bed"

    peaks[["chr", "start", "end"]].to_csv(output_file,
                                          sep="\t",
                                          header=False,
                                          index=False)

    return output_file

# verify that UCSC liftOver is available
if shutil.which("liftOver") is None:
    raise FileNotFoundError("UCSC liftOver executable not found in PATH")

# define genome assembly conversion files
forward_chain = chip_dir / "mm10ToMm39.over.chain"
reverse_chain = chip_dir / "mm39ToMm10.over.chain"

# check that chain files exist
for chain in [forward_chain, reverse_chain]:
    if not chain.is_file():
        raise FileNotFoundError(f"Chain file not found: {chain}")

# define a function to convert an mm39 plotting window to mm10
def get_mm10_window(gene, mode="significant"):

    # retrieve the gene's mm39 plotting window
    chrom, start, end = get_gene_window(gene, mode=mode)

    # define input and output files
    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)

    window_file = output_dir / f"{gene}_mm39_window.bed"
    mapped_file = output_dir / f"{gene}_mm10_window.bed"
    unmapped_file = output_dir / f"{gene}_unmapped_window.bed"

    # write the mm39 genomic interval
    with open(window_file, "w") as f:
        f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")

    # convert mm39 coordinates to mm10
    result = subprocess.run(["liftOver",
                             str(window_file),
                             str(reverse_chain),
                             str(mapped_file),
                             str(unmapped_file)],
                             capture_output=True,
                             text=True)

    if result.returncode != 0:
        raise RuntimeError(f"liftOver failed for {gene}: {result.stderr}")

    # read mapped genomic intervals
    mapped = pd.read_csv(mapped_file,
                         sep="\t",
                         header=None,
                         names=["chr", "start", "end", "gene"])

    if len(mapped) != 1:
        raise ValueError(f"Expected one mm10 interval for {gene}, found {len(mapped)}")

    return (mapped["chr"].iloc[0],
            int(mapped["start"].iloc[0]),
            int(mapped["end"].iloc[0]))

# identify H3K27ac ChIP-seq bigWig replicates
chip_files = sorted(chip_dir.glob("GSM890771[345]_*.bw"))

if len(chip_files) != 3:
    raise ValueError(f"Expected 3 H3K27ac replicates, found {len(chip_files)}")

print("H3K27ac bigWig files:")
for file in chip_files:
    print(file.name)

# define a function to lift H3K27ac signal intervals from mm10 to mm39
def liftover_h3k27ac(file, gene, chrom10, start10, end10):

   # extract H3K27ac signal from the mm10 bigWig
   with pyBigWig.open(str(file)) as bw:
    bw_chrom = chrom10.removeprefix("chr")
    intervals = bw.intervals(bw_chrom, start10, end10)

    if not intervals:
        raise ValueError(f"No H3K27ac signal found for {gene} in {file.name}")

    signal = pd.DataFrame(intervals,
                          columns=["start", "end", "signal"])

    signal.insert(0, "chr", chrom10)

    # clip intervals to the selected mm10 window
    signal["start"] = signal["start"].clip(lower=start10)
    signal["end"] = signal["end"].clip(upper=end10)
    signal = signal.loc[signal["end"] > signal["start"]].copy()

    # define output files
    sample = file.name.split("_")[0]
    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)

    input_file = output_dir / f"{gene}_{sample}_mm10.bed"
    output_file = output_dir / f"{gene}_{sample}_mm39.bed"
    unmapped_file = output_dir / f"{gene}_{sample}_unmapped.bed"

    signal.to_csv(input_file, sep="\t", header=False, index=False)

    # lift signal intervals to mm39
    result = subprocess.run(["liftOver",
                             str(input_file),
                             str(forward_chain),
                             str(output_file),
                             str(unmapped_file)],
                             capture_output=True,
                             text=True)

    if result.returncode != 0:
        raise RuntimeError(f"H3K27ac liftOver failed for {gene}, {sample}: {result.stderr}")

    # read mapped signal intervals
    if output_file.stat().st_size == 0:
        raise ValueError(f"No H3K27ac intervals mapped for {gene}, {sample}")

    mapped = pd.read_csv(output_file, sep="\t", header=None, names=["chr", "start", "end", "signal"])

    mapped = mapped.sort_values(["chr", "start"]).reset_index(drop=True)

    print(f"{gene} {sample}: {len(mapped)}/{len(signal)} intervals mapped")

    return mapped

# define a function to export mean H3K27ac coverage
def export_h3k27ac_mean(gene, bin_size=25, mode="significant"):

    # retrieve mm39 and mm10 plotting windows
    chrom, start, end = get_gene_window(gene, mode=mode)
    chrom10, start10, end10 = get_mm10_window(gene, mode=mode)

    # lift H3K27ac signal from all three replicates to mm39
    replicates = [liftover_h3k27ac(file, gene, chrom10, start10, end10) for file in chip_files]

    # define common mm39 genomic bins
    n_bins = int(np.ceil((end - start) / bin_size))
    edges = np.floor(np.linspace(start, end, n_bins + 1)).astype(int)

    # calculate coverage for each replicate using interval overlaps
    binned = []

    for replicate in replicates:
        signal = np.zeros(n_bins, dtype=float)
        covered = np.zeros(n_bins, dtype=float)

        intervals = replicate.loc[(replicate["chr"] == chrom) &
                                  (replicate["end"] > start) &
                                  (replicate["start"] < end)]

        for row in intervals.itertuples(index=False):
            interval_start = max(int(row.start), start)
            interval_end = min(int(row.end), end)

            first_bin = max(0, np.searchsorted(edges, interval_start, side="right") - 1)
            last_bin = min(n_bins - 1, np.searchsorted(edges, interval_end, side="left"))

            for i in range(first_bin, last_bin + 1):
                overlap = max(0, min(interval_end, edges[i + 1]) - max(interval_start, edges[i]))

                if overlap > 0:
                    signal[i] += overlap * float(row.signal)
                    covered[i] += overlap

        # calculate mean signal within each genomic bin
        binned.append(np.divide(signal,
                                covered,
                                out=np.full(n_bins, np.nan),
                                where=covered > 0))

    # average the three biological replicates
    with np.errstate(invalid="ignore"):
        h3k27ac_mean = np.nanmean(np.vstack(binned), axis=0)

    # replace uncovered bins with zero for plotting
    h3k27ac_mean = np.nan_to_num(h3k27ac_mean, nan=0.0)

    # export mean coverage as bedGraph
    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{gene}_H3K27ac_mean.bedgraph"

    with open(output_file, "w") as f:
        for i, value in enumerate(h3k27ac_mean):
            if edges[i + 1] > edges[i]:
                f.write(f"{chrom}\t{edges[i]}\t{edges[i + 1]}\t{value:.6f}\n")

    return output_file

# define a function to export promoter-DAR PLAC-seq interactions
def export_plac_interactions(gene, promoter_window=2000, significant_only=True):

    # retrieve gene annotation and transcription start site
    annotation = get_gene_annotation(gene)

    if annotation.empty:
        raise ValueError(f"Gene {gene} not found in GENCODE annotation")

    chrom = annotation["chr"].iloc[0]
    tss = int(annotation["tss"].iloc[0])

    # identify interactions overlapping the promoter region
    plac_gene = plac_interactions.loc[(plac_interactions["chr1"] == chrom) &
                                      (plac_interactions["chr2"] == chrom) &
                                      (((plac_interactions["start1"] < tss + promoter_window) &
                                        (plac_interactions["end1"] > tss - promoter_window)) |
                                        ((plac_interactions["start2"] < tss + promoter_window) &
                                         (plac_interactions["end2"] > tss - promoter_window)))].copy()

    # retrieve associated atac peaks using the selected significance filter
    dar_gene = (get_gene_peaks(gene, significant_only=significant_only)
                [["chr", "start", "end"]].drop_duplicates())

    # retain promoter interactions with an anchor overlapping a DAR
    plac_dar = plac_gene.loc[plac_gene.apply(
            lambda contact: (((dar_gene["chr"] == contact["chr1"]) &
                              (dar_gene["start"] < contact["end1"]) &
                              (dar_gene["end"] > contact["start1"])) |
                              ((dar_gene["chr"] == contact["chr2"]) &
                               (dar_gene["start"] < contact["end2"]) &
                               (dar_gene["end"] > contact["start2"]))).any(),
                               axis=1)].drop_duplicates()

    # export interactions as BEDPE
    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{gene}_PLAC_links.bedpe"

    plac_dar[["chr1", "start1", "end1", "chr2", "start2", "end2"]].to_csv(output_file,
                                                                          sep="\t",
                                                                          header=False,
                                                                          index=False)

    return output_file

# define a function to generate pyGenomeTracks configuration files
def create_atac_config(gene, tracks):

    output_dir = locus_dir / gene
    output_dir.mkdir(parents=True, exist_ok=True)
    config_file = output_dir / f"{gene}_atac.ini"

    groups = [("Male", "Basal", "ATAC-seq", "#78A8D1"),
              ("Male", "60 min", "ATAC-seq", "#47749B"),
              ("Female", "Basal", "ATAC-seq", "#D78FA9"),
              ("Female", "60 min", "ATAC-seq", "#B65679")]

    # shared ATAC-seq y-axis maximum
    ymax = int(np.ceil(max(pd.read_csv(file, sep="\t", header=None, usecols=[3]).iloc[:, 0].max()
                           for file in tracks.values()) * 1.1))

    with open(config_file, "w") as f:

        # PLAC-seq interactions
        plac_file = output_dir / f"{gene}_PLAC_links.bedpe"

        if plac_file.exists() and plac_file.stat().st_size > 0:
            f.write(f"\n[PLAC_links]\n"
                    f"file = {plac_file}\n"
                    f"title = PLAC-seq Ferris et al.)\n"
                    f"file_type = links\n"
                    f"height = 1.5\n"
                    f"color = #5FA9A2\n"
                    f"line_width = 1\n"
                    f"links_type = arcs\n")

        # genomic coordinate axis
        f.write("\n[x-axis]\n"
                "fontsize = 6\n")

        # gene annotation
        f.write(f"\n[gene_annotation]\n"
                f"file = {output_dir / f'{gene}_region.gtf'}\n"
                f"title = {gene}\n"
                f"height = 0.3\n"
                f"prefered_name = gene_name\n"
                f"merge_transcripts = true\n"
                f"labels = false\n"
                f"style = flybase\n"
                f"display = collapsed\n"
                f"color = #9F9F9F\n"
                f"border_color = #9F9F9F\n"
                f"color_backbone = #9F9F9F\n"
                f"line_width = 0.5\n")

        # ATAC-seq coverage tracks
        for sex, condition, assay, color in groups:

            track_file = tracks[(sex, condition)]

            f.write(f"\n[{sex}_{condition.replace(' ', '')}]\n"
                    f"file = {track_file}\n"
                    f"title = {sex} {condition} ({assay})\n"
                    f"height = 1.5\n"
                    f"color = {color}\n"
                    f"min_value = 0\n"
                    f"max_value = {ymax}\n"
                    f"type = fill\n")

            # space between coverage tracks
            f.write("\n[spacer]\nheight = 0.2\n")

        # H3K27ac coverage track
        h3k27ac_file = output_dir / f"{gene}_H3K27ac_mean.bedgraph"

        if h3k27ac_file.exists():

            h3k27ac_max = pd.read_csv(h3k27ac_file,
                                      sep="\t",
                                      header=None,
                                      usecols=[3]).iloc[:, 0].max()

            h3k27ac_ymax = int(np.ceil(h3k27ac_max * 1.05 / 25) * 25)

            f.write(f"\n[H3K27ac_coverage]\n"
                    f"file = {h3k27ac_file}\n"
                    f"title = H3K27ac (Ferris et al.)\n"
                    f"height = 1.5\n"
                    f"color = #A995C9\n"
                    f"min_value = 0\n"
                    f"max_value = {h3k27ac_ymax}\n"
                    f"type = fill\n")

        # significant DAR highlights, with existing BED padding
        dar_highlights = output_dir / f"{gene}_DAR_highlights.bed"

        if dar_highlights.exists() and dar_highlights.stat().st_size > 0:
            f.write(f"\n[vhighlight]\n"
                    f"file = {dar_highlights}\n"
                    f"type = vhighlight\n"
                    f"color = #D3D3D3\n"
                    f"alpha = 0.5\n")

    return config_file

# create a pyGenomeTracks wrapper for figure styling
wrapper_file = locus_dir / "run_pygenometracks.py"

wrapper_file.write_text("""import numpy as np
np.float = float

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pygenometracks.tracks import GenomeTrack

# typography
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['font.weight'] = 'normal'
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['axes.titleweight'] = 'normal'
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['pdf.fonttype'] = 42

# thin coverage scale lines to 0.5 pt
_original_plot_y_axis = GenomeTrack.plot_y_axis

def _thin_y_axis(self, ax, plot_axis, *args, **kwargs):
    before = len(ax.lines)
    _original_plot_y_axis(self, ax, plot_axis, *args, **kwargs)

    for line in ax.lines[before:]:
        line.set_linewidth(0.5)

GenomeTrack.plot_y_axis = _thin_y_axis

# put Ferris et al. on a new line in track labels
_original_plot_label = GenomeTrack.plot_label

def _wrapped_plot_label(self, ax, *args, **kwargs):
    title = self.properties.get('title', '')

    if '(Ferris et al.)' in title:
        self.properties['title'] = title.replace(' (Ferris et al.)', '\\n(Ferris et al.)')

    return _original_plot_label(self, ax, *args, **kwargs)

GenomeTrack.plot_label = _wrapped_plot_label

# final figure styling
_original_savefig = Figure.savefig

def _styled_savefig(self, *args, **kwargs):

    for ax in self.axes:

        # thin standard axis spines
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

        # remove x-axis tick marks but retain coordinate labels
        ax.tick_params(
            axis='x',
            which='both',
            length=0,
            width=0
        )

        for tick in ax.xaxis.get_major_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)

        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_visible(False)
            tick.tick2line.set_visible(False)

    return _original_savefig(self, *args, **kwargs)

Figure.savefig = _styled_savefig

from pygenometracks.plotTracks import main
main()
""")

# verify that the generated wrapper has valid Python syntax
import ast
ast.parse(wrapper_file.read_text())
print("Wrapper syntax OK")

# define a function to generate a locus plot for each gene
def plot_locus(gene, mode="significant", significant_only=True):

    print(f"\nGenerating locus plot for {gene}...")

    # define genomic region
    chrom, start, end = get_gene_window(gene, mode=mode)

    # export mean ATAC-seq coverage for each group
    tracks = {}

    for sex in ["Male", "Female"]:
        for condition in ["Basal", "60 min"]:
            tracks[(sex, condition)] = export_group_bedgraph(gene, sex, condition, mode=mode)

    # export gene annotation
    export_region_gtf(gene)

    # export significant DAR highlights
    export_dar_highlights(gene)

    # export promoter-DAR PLAC-seq interactions
    export_plac_interactions(gene, significant_only=significant_only)

    # export mean H3K27ac coverage
    export_h3k27ac_mean(gene, mode=mode)

    # generate pyGenomeTracks configuration
    config_file = create_atac_config(gene, tracks)

    # render locus plot as PDF
    output_file = locus_dir / gene / f"{gene}_atac.pdf"

    result = subprocess.run([sys.executable,
                             str(wrapper_file),
                             "--tracks", str(config_file),
                             "--region", f"{chrom}:{start}-{end}",
                             "--fontSize", "5.5",
                             "--outFileName", str(output_file),
                             "--width", "15.92",
                             "--dpi", "300"],
                             capture_output=True,
                             text=True)

    if result.returncode != 0:
        raise RuntimeError(f"pyGenomeTracks failed for {gene}:\n{result.stderr}")

    print(f"Saved: {output_file}")

    return output_file

# generate locus plots for selected genes
if __name__ == "__main__":

    # define plotting settings for candidate genes
    gene_settings = {
        "Sgk1":   {"mode": "significant", "significant_only": True},
        "Dio2":   {"mode": "significant", "significant_only": True},
        "Zfp36":  {"mode": "significant", "significant_only": True},
        "Ctla2a": {"mode": "all",         "significant_only": False},
        "Ska1":   {"mode": "significant", "significant_only": True},
        "Htr3a":  {"mode": "significant", "significant_only": True},
        "Htr3b":  {"mode": "significant", "significant_only": True},
        "Fkbp5":  {"mode": "significant", "significant_only": True}
    }

    # generate locus plots using gene-specific settings
    for gene, settings in gene_settings.items():
        plot_locus(gene, **settings)
