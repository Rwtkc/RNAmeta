% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAmeta_gene_type_plot.R
\name{RNAmeta_gene_type_plot}
\alias{RNAmeta_gene_type_plot}
\title{Gene Type Plot}
\usage{
RNAmeta_gene_type_plot(preliminary_analysis_data = NULL, set_group_name = NULL)
}
\arguments{
\item{preliminary_analysis_data}{A list containing the preliminary analysis data. This typically includes genomic ranges, GFF file, and other results.
\itemize{
\item The second element should be genomic ranges (\code{peak_gr}).
\item The third element should be a GFF file (\code{gff_file}).
}}

\item{set_group_name}{A vector of user-defined group names to map to the unique \code{Sample} identifiers in the \code{all_sites} data.
The length of \code{set_group_name} must match the number of unique \code{Sample} values in the data.}
}
\value{
A \code{ggplot} object representing a bar plot of gene types in the dataset. The plot shows the count of each gene type across samples,
with bars colored based on the sample groups. The labels indicate the count values, and the position of the legend depends on the number
of unique samples (if more than one group is present, the legend is shown).
}
\description{
Gene Type Plot
}
