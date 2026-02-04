calcDensityGFF <- function (gff_file, karyotype_info, feature = "gene", window = 1e+05) 
{
    gff <- read.table(gff_file, header = FALSE, comment.char = "#")
    if (is.character(karyotype_info)) {
        karyotype <- read.table(karyotype_info, header = FALSE)
    }
    else {
        karyotype <- karyotype_info
    }
    all_chrs <- karyotype[[1]]
    gff <- gff[gff[[1]] %in% all_chrs & gff[[3]] == feature, 
        ]
    density <- data.frame()
    for (chr in all_chrs) {
        chr_end <- karyotype[karyotype[[1]] == chr, 3]
        tmp <- data.frame(table(cut(gff[gff[[1]] == chr, 4], 
            breaks = unique(c(seq(0, chr_end, window), chr_end)))))
        tmp <- tidyr::separate(tmp, 1, into = c("Start", "End"), 
            sep = ",")
        tmp$Start <- as.numeric(gsub("\\(", "", tmp$Start)) + 
            1
        tmp$End <- as.numeric(gsub("\\]", "", tmp$End))
        tmp$Chr <- chr
        tmp <- tmp[c(4, 1:3)]
        colnames(tmp) <- c("Chr", "Start", "End", "Count")
        tmp[nrow(tmp), "End"] <- chr_end
        density <- rbind(density, tmp)
    }
    return(density)
}