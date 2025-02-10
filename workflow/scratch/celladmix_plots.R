#' @export
de_volcano_plot <- function(de_df, tumor_marker_genes, genes_to_label) {
  ## making volcano plot for original DE
  de_df$padj <- 2*pnorm(abs(de_df$Z), mean = 0, sd = 1, lower.tail = FALSE) # convert z back to padj
  de_df$logpv <- -log10(de_df$padj)
  ndx_lab <- which(de_df$Gene %in% genes_to_label)
  de_df$top_mark_other <- NA
  de_df[ndx_lab,'top_mark_other'] <- de_df$Gene[ndx_lab]
  de_df$logpv[de_df$logpv>30] <- 30
  ymax <- max(de_df$logpv,na.rm=TRUE) + 1
  de_df$mark_other <- de_df$Gene %in% tumor_marker_genes
  de_df$mark_other <- factor(de_df$mark_other,levels=c(TRUE,FALSE))
  levels(de_df$mark_other) <- c('Malignant marker', 'other')

  myColors <- c('red4','grey70')

  # Helper variables
  limits <- range(de_df$logpv,na.rm = TRUE,finite = TRUE)
  step   <- diff(limits) * 0.1
  size   <- 0.45 * step

  de_df$mark_other <- factor(de_df$mark_other,levels=c('Malignant marker','other'))
  bottom_lim <- min(limits[1] - as.numeric(de_df$mark_other) * step - size) - .5

  p <- ggplot(de_df,aes(x=M,y=logpv,color=mark_other,label=top_mark_other)) +
    geom_segment(aes(y=ymax, yend=0, x =  log2(1.5), xend=log2(1.5)),color='red',linetype='dashed') +
    geom_segment(aes(y=ymax, yend=0, x = -log2(1.5), xend=-log2(1.5)),color='red',linetype='dashed') +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=-log10(.05),color='red',linetype='dashed') +
    geom_segment(
      aes(
        color = mark_other,
        xend = M,
        y    = limits[1] - as.numeric(mark_other) * step + size,
        yend = limits[1] - as.numeric(mark_other) * step - size
      )
    ) +
    ggrepel::geom_text_repel(size=4, show.legend = FALSE) +
    geom_vline(xintercept=0) +
    geom_point(data=de_df[de_df$mark_other=='other',],alpha=.5,size=1.5) +
    geom_point(data=de_df[de_df$mark_other=='Malignant marker',],alpha=.5,size=1.5) +
    scale_y_continuous(limits=c(bottom_lim,ymax),expand = c(0, 0)) +
    xlab('logFC') +
    ylab('-log10(Padj)') +
    scale_colour_manual(breaks = c('Malignant marker','other'), values = myColors) +
    theme_classic(base_line_size = .5) +
    theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank())

  return(p)
}

#' @export
plot_nmf_loading <- function(loading, factor.name=NULL) {
    Hi <- loading %>% {data.frame(frac=.)} %>% tibble::rownames_to_column('gene')

    gene_order <- Hi %>% arrange(-frac) %>% pull(gene) %>% head(20) %>% rev

    if (is.null(factor.name)) {
        factor.name <- 'Factor'
    }

    gg <- Hi %>%
      filter(gene %in% gene_order) %>%
      mutate(gene = factor(gene, gene_order)) %>%
      ggplot(aes(x = frac, y = gene)) +
      xlab('Loadings fraction') +
      ylab(factor.name) +
      geom_col() +
      theme_classic(base_line_size = .5) +
      scale_x_continuous(
        breaks=round(seq(0, max(Hi$frac), max(Hi$frac)/2),digits = 2),
        position="top",
        expand=c(0, 0, 0.025, 0)
      ) +
      theme(
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8,face="bold")
      )

    return(gg)
}

#' @export
plot_expression_comparison <- function(
    cm_orig, cm_clean, markers_native, markers_admix, admixture_type, native_type,
    combine=TRUE
  ) {
  dat_orig_norm <- sweep(cm_orig, 2, colSums(cm_orig), '/')
  dat_cln_norm <- sweep(cm_clean, 2, colSums(cm_clean), '/')

  admix_av_expr_before <- rowMeans(dat_orig_norm[markers_admix,])
  native_av_expr_before <- rowMeans(dat_orig_norm[markers_native,])

  admix_av_expr_after <- rowMeans(dat_cln_norm[markers_admix,])
  native_av_expr_after <- rowMeans(dat_cln_norm[markers_native,])

  tmp1a <- cbind.data.frame(markers_admix, admix_av_expr_before,'original')
  colnames(tmp1a) <- c('gene','av_expr','type')
  tmp1b <- cbind.data.frame(markers_admix,admix_av_expr_after,'cleaned')
  colnames(tmp1b) <- c('gene','av_expr','type')
  tmp1 <- rbind.data.frame(tmp1a,tmp1b)
  tmp1$type <- factor(tmp1$type,levels=c('original','cleaned'))

  p1 <- ggplot(tmp1,aes(x=gene,y=av_expr,fill=type)) +
      geom_bar(stat="identity", position=position_dodge()) +
      xlab(paste(admixture_type, 'marker')) +
      ylab('Average expression') +
      ggtitle(paste(admixture_type, 'marker expression\n in', native_type)) +
      theme(plot.title = element_text(hjust = 0.5))

  tmp1a <- cbind.data.frame(markers_native,native_av_expr_before,'original')
  colnames(tmp1a) <- c('gene','av_expr','type')
  tmp1b <- cbind.data.frame(markers_native,native_av_expr_after,'cleaned')
  colnames(tmp1b) <- c('gene','av_expr','type')
  tmp1 <- rbind.data.frame(tmp1a,tmp1b)
  tmp1$type <- factor(tmp1$type,levels=c('original','cleaned'))

  p2 <- ggplot(tmp1,aes(x=gene,y=av_expr,fill=type)) +
      geom_bar(stat="identity", position=position_dodge()) +
      xlab(paste(native_type, 'marker')) +
      ylab('Average expression') +
      ggtitle(paste(native_type, 'marker expression\n in', native_type)) +
      theme(plot.title = element_text(hjust = 0.5))

  if (!combine) {
    return(list(p1,p2))
  }

  return(cowplot::plot_grid(p1, p2, nrow=1))
}

#' @export
plot_de_comparison <- function(de_orig, de_clean, genes_to_label, tumor_marker_genes) {
  z_thresh <- qnorm(.01/2,lower.tail = FALSE)
  tmp <- de_orig %>%
    {data.frame(orig=.[,'Z'], clean=de_clean[rownames(.),'Z'], row.names=rownames(.))}

  tmp %<>% as.matrix() %>% pmin(10) %>% pmax(-10) %>% as.data.frame()

  ## trying to color points by tumor marker set belonging
  tmp$markers <- ifelse(rownames(tmp) %in% tumor_marker_genes, 'Malignant marker', 'other') %>%
    factor(levels=c('other', 'Malignant marker'))

  ndx_lab <- match(genes_to_label, rownames(tmp))
  genes_to_label <- genes_to_label[!is.na(ndx_lab)]
  ndx_lab <- ndx_lab[!is.na(ndx_lab)]
  tmp$glab <- NA
  tmp[ndx_lab,'glab'] <- genes_to_label

  p <- ggplot(tmp,aes(x=orig,y=clean,color=markers,label=glab)) +
    geom_point(data=tmp[tmp$markers=='other',],alpha=.5,size=1.5) +
    geom_point(data=tmp[tmp$markers=='Malignant marker',],alpha=.5,size=1.5) +
    ggrepel::geom_text_repel(size=4, show.legend = FALSE) +
    xlab('Z-score (original)') +
    ylab('Z-score (cleaned)') +
    geom_hline(yintercept = z_thresh,color='red',linetype='dashed') +
    geom_hline(yintercept = -1*z_thresh,color='red',linetype='dashed') +
    geom_hline(yintercept = 0) +
    geom_abline(slope = 1,intercept=0) +
    scale_colour_manual(breaks = c('Malignant marker','other'), values=c('red4','grey70')) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

  return(p)
}