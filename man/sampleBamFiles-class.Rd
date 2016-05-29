\name{sampleBamFiles-class}
\Rdversion{1.1}
\docType{class}
%
\alias{sampleBamFiles-class}
\alias{initialize, sampleBamFiles}
%
\alias{[,sampleBamFiles,ANY,ANY,ANY-method}
%
\alias{show,sampleBamFiles-method}
\alias{length,sampleBamFiles-method}
%
\alias{sampleBamFiles}
\alias{sampleBamFiles-methods}
\alias{sampleBamFiles,integer-method}
\alias{sampleBamFiles,numeric-method}
\alias{sampleBamFiles,character-method}
\alias{sampleBamFiles,factor-method}
%
\alias{bamFiles}
\alias{bamFiles-methods}
\alias{bamFiles,sampleBamFiles-method}
%
\alias{bamFiles<-}
\alias{bamFiles<--methods}
\alias{bamFiles<-,sampleBamFiles,character-method}
\alias{bamFiles<-,sampleBamFiles,character,logical-method}
%
\alias{bamIdxFiles}
\alias{bamIdxFiles-methods}
\alias{bamIdxFiles,sampleBamFiles-method}
%
\alias{bamIdxFiles<-}
\alias{bamIdxFiles<--methods}
\alias{bamIdxFiles<-,sampleBamFiles,character-method}
%
\alias{createIndex,sampleBamFiles,character-method}
\alias{createIndex,sampleBamFiles,missing-method}
%
\alias{nAligns,sampleBamFiles-method}
%
\alias{nAligns<-}
\alias{nAligns<--methods}
\alias{nAligns<-,sampleBamFiles,numeric-method}
\alias{nAligns<-,sampleBamFiles,integer-method}
\alias{nAligns<-,sampleBamFiles,character,numeric-method}
%
\alias{sampleLabels}
\alias{sampleLabels-methods}
\alias{sampleLabels,sampleBamFiles-method}
%
\alias{sampleLabels<-}
\alias{sampleLabels<--methods}
\alias{sampleLabels<-,sampleBamFiles,character-method}
\alias{sampleLabels<-,sampleBamFiles,factor-method}
\alias{sampleLabels<-,sampleBamFiles,integer-method}
%
\alias{sampleGroups}
\alias{sampleGroups-methods}
\alias{sampleGroups,sampleBamFiles-method}
%
\alias{sampleGroups<-}
\alias{sampleGroups<--methods}
\alias{sampleGroups<-,sampleBamFiles,character-method}
\alias{sampleGroups<-,sampleBamFiles,factor-method}
%
\alias{groupTable}
\alias{groupTable-methods}
\alias{groupTable,sampleBamFiles-method}
%
\alias{groupTable<-}
\alias{groupTable<--methods}
\alias{groupTable<-,sampleBamFiles,character-method}
\alias{groupTable<-,sampleBamFiles,data.frame-method}
%
\alias{checkBamFiles}
\alias{checkBamFiles-methods}
\alias{checkBamFiles,sampleBamFiles-method}
%
\alias{bamCountAll,sampleBamFiles,logical-method}
\alias{bamCountAll,sampleBamFiles,missing-method}
%
\title{Class \code{"sampleBamFiles"}: Data on multiple BAM files from
an experimental setting.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
The object contains data on location of BAM files, experimental
group assignment and number of aligns per BAM file.
The align-numbers are used for normalisation when alignment
depth is plotted for single genes.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Creation of objects
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Objects from the Class}{Objects can be created by calls of the form
    \code{bs<-sampleBamFiles(object)} where object is a numeric vector of
    length1.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Slots
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Slots}{
    \describe{
        \item{\code{bamFiles}:}{\code{"character"}:
            Location of BAM files}
        %
        \item{\code{bamIdxFiles}:}{\code{"character"}:
            Location of BAM index files}
        \item{\code{nAligns}:}{\code{"numeric"}:
            Total number of alignments in each BAM file.}
        %
        \item{\code{group}:}{\code{"factor"}:
            Group assignment}
        %
        \item{label}{\code{"character"}:
            Short textual identifier for each sample.}
        %
        \item{length}{\code{"integer"}:
            Vector length for bamFiles, bamIdxFiles, nAligns,
            group and label.}
        %
        \item{ev}{\code{"environment"}:
            Contains additional data  (e.g. group table).}
  }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{show}{\code{signature(object="sampleBamFiles")}:
            Prints a short message with some summarizing data.}
        %
        \item{bamFiles}{\code{signature(object="sampleBamFiles"}:
            Returns names of BAM files.}
        %
        \item{bamIdxFiles}{\code{signature(object="sampleBamFiles")}:
            Returns names of BAM index files.}
        %
        \item{length}{\code{signature(object="sampleBamFiles")}:
            Returns number of BAM files.}
        %
        \item{nAligns}{\code{signature(object="sampleBamFiles")}:
            Returns total number of aligns in BAM files.}
        %
        \item{sampleLabels}{\code{signature(object="sampleBamFiles")}:
            Returns sample labels for BAM files.}
        %
        \item{sampleGroups}{\code{signature(object="sampleBamFiles")}:
            Returns group assignment for BAM files.}
        %
        \item{groupTable}{\code{signature(object="sampleBamFiles")}:
            Returns group table (if present).}
        %
}}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam<-system.file("extdata", "accepted_hits.bam", package="rbamtools")
bs <- sampleBamFiles(bam)
sampleLabels(bs) <- "s1"
sampleGroups(bs) <- "g1"
checkBamFiles(bs)
nAligns(bs) <- bamCountAll(bs)
bs
}
\keyword{sampleBamFiles}
