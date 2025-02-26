% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAmeta_get_transcription_plot.R
\name{RNAmeta_get_transcription_plot}
\alias{RNAmeta_get_transcription_plot}
\title{RNAmeta_get_transcription_plot}
\usage{
RNAmeta_get_transcription_plot(
  transcription_data = NULL,
  set_group_name = NULL
)
}
\arguments{
\item{transcription_data}{A list containing two elements:
\itemize{
\item The first element is a data frame or data table containing information on transcription coverage (\code{cvg_list}).
\item The second element is a list with data related to peak heatmap data for the transcriptional regions (i.e., \code{peak_heatmap_data}).
}}

\item{set_group_name}{An optional vector of custom group names to replace the default sample names in the data. The length of this vector must match the number of unique sample groups in the data.
If \code{NULL}, the function will use the original sample names.}
}
\value{
A list of ggplot objects. The list contains the following plots:
\enumerate{
\item A density plot of transcription coverage around the transcription start site (TSS).
\item A density plot of transcription coverage around the transcription end site (TES).
\item A heatmap of transcription signal for the TSS region.
\item A heatmap of transcription signal for the TES region.
}
}
\description{
RNAmeta_get_transcription_plot
}
