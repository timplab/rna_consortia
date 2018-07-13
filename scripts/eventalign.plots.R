library(tidyverse)
library(stringr)
library(cowplot)
library(ggridges)
library(ggplot2)

##load in eventalign files
oligo = read.delim("/dilithium/Data/Nanopore/Analysis/180323_m6aoligo.train/180612_evalign.m6aoligo.sub.txt", sep = "\t")
colnames(oligo) = c("contig", "position", "reference_kmer", "read_index", "strand", "event_index", "event_level_mean","event_stdv","event_length","model_kmer", "model_mean", "model_stdv", "standardized_level")

dRNA = read.delim("dRNA.EEF2.GGACU.evalign.txt", sep = "\t") %>%
    filter(position == 3976325) %>%
    mutate(type = "nvRNA") %>%
    mutate(modified = "TRUE")

IVT = read.delim("IVT.EEF2.GGACU.evalign.txt", sep = "\t") %>%
    filter(position == 3976325) %>%
    mutate(type = "nvRNA") %>%
    mutate(modified = "FALSE")

oligo.T =
    oligo %>%
    mutate(modified=str_detect(reference_kmer, "GGACT")) %>%
    filter(position == 1754) %>%
    mutate(type = "oligo")

oligo.F =
    oligo %>%
    mutate(modified=str_detect(reference_kmer, "GGCCT"))%>%
    filter(position == 1773) %>%
    mutate(type = "oligo")

joy = rbind(oligo.T, oligo.F, dRNA, IVT)

##plot

pdf("EEFcomp.hg38.pdf")

ggplot(joy, aes(x=event_level_mean, y=type, fill=modified)) +
    geom_density_ridges2(alpha = 0.9, panel_scaling=FALSE, size = 1) +
    scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    theme_bw()

ggplot(joy, aes(x=event_level_mean, colour=modified)) +
    stat_ecdf(alpha = 0.9, size = 1) +
    facet_grid(. ~ type, scales = "free_y") +
    #scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed", size = 1) +
    theme_classic()

ggplot(joy, aes(x=event_level_mean, colour=modified)) +
    geom_density(alpha = 0.9, size = 1) +
    facet_grid(. ~ type, scales = "free_y") +
    #scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed", size = 1) +
    theme_classic()

ggplot(joy, aes(x=event_level_mean, colour=modified)) +
    geom_freqpoly(alpha = 0.9, size = 1) +
    facet_grid(. ~ type, scales = "free_y") +
    #scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed", size = 1) +
    theme_classic()

dev.off()

#####for EEF1A1 isoforms

readmatch = read_tsv("/dilithium/Data/Nanopore/rna/isoforms/EEF1A1.reads.txt") %>%
    `colnames<-`(c("read_name", "isoform", "gene"))

pdf("EEF1A1.GGACUbyiso.pdf")

######First, plots with annotated transcripts pooled to their txID and discounting alternative stops and starts

dRNA = read_tsv("/dilithium/Data/Nanopore/rna/isoforms/dRNA.EEF1A1.mq5.GGACU.evalign.txt") %>%
    inner_join(readmatch, by = "read_name") %>%
    filter(model_kmer != "NNNNN") %>%
    filter(reference_kmer == "AGUCC") %>%
    filter(grepl("ENST", isoform)) %>%
    separate(., col=isoform, into=c("isoform", "gene"), sep = "_") %>%
    unite(coordinates, contig, position, sep = ":") %>%
    filter(coordinates == "chr6:73518769")

ggplot(dRNA, aes(x=event_level_mean, y=isoform, colour=isoform)) +
    geom_density_ridges2(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="right") +
    labs(title = "Byannot.tx")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_density(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    facet_grid(. ~isoform) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="right") +
    labs(title = "Byannot.tx")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_freqpoly(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    #theme(legend.position="bottom")+
    labs(title = "Byannot.tx")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    stat_ecdf(alpha = 0.9, size = 1) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    #theme(legend.position="bottom")+
    labs(title = "Byannot.tx")

####Second, throw in the alternative stops and starts

dRNA = read_tsv("/dilithium/Data/Nanopore/rna/isoforms/dRNA.EEF1A1.mq5.GGACU.evalign.txt") %>%
    inner_join(readmatch, by = "read_name") %>%
    filter(model_kmer != "NNNNN") %>%
    filter(reference_kmer == "AGUCC") %>%
    filter(grepl("ENST", isoform)) %>%
    unite(coordinates, contig, position, sep = ":") %>%
    filter(coordinates == "chr6:73518769")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_density(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    facet_grid(. ~coordinates) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_freqpoly(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    stat_ecdf(alpha = 0.9, size = 1) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE")

###ok, so let's get rid of isoforms with low counts in this context

dRNA = read_tsv("/dilithium/Data/Nanopore/rna/isoforms/dRNA.EEF1A1.mq5.GGACU.evalign.txt") %>%
    inner_join(readmatch, by = "read_name") %>%
    filter(model_kmer != "NNNNN") %>%
    filter(reference_kmer == "AGUCC") %>%
    filter(grepl("ENST", isoform)) %>%
    unite(coordinates, contig, position, sep = ":") %>%
    filter(coordinates == "chr6:73518769") %>%
    group_by(isoform) %>%
    filter(n() > 5)

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_density(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    facet_grid(. ~coordinates) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE.5X")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_freqpoly(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE.5X")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    stat_ecdf(alpha = 0.9, size = 1) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE.5X")

dRNA = read_tsv("/dilithium/Data/Nanopore/rna/isoforms/dRNA.EEF1A1.mq5.GGACU.evalign.txt") %>%
    inner_join(readmatch, by = "read_name") %>%
    filter(model_kmer != "NNNNN") %>%
    filter(reference_kmer == "AGUCC") %>%
    filter(grepl("ENST", isoform)) %>%
    unite(coordinates, contig, position, sep = ":") %>%
    filter(coordinates == "chr6:73518769") %>%
    group_by(isoform) %>%
    filter(n() > 10)

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_density(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    facet_grid(. ~coordinates) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE.10X")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    geom_freqpoly(alpha = 0.9, panel_scaling=FALSE, size = 1, fill = NA) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE.10X")

ggplot(dRNA, aes(x=event_level_mean, colour=isoform)) +
    stat_ecdf(alpha = 0.9, size = 1) +
    #scale_y_discrete(expand=c(0.01, 0)) +
    facet_grid(. ~coordinates) +
    scale_x_continuous(limits=c(100,150))+
    geom_vline(xintercept = 123.8, linetype = "dashed") +
    scale_fill_brewer(palette="Set3") +
    theme(legend.position="none")+
    labs(title = "Byannot.tx.plusTSS.TSE.10X")

dev.off()
