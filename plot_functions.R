library(ggplot2)

#######################################################################################
# DNAm AGE VS CHRONOLOGICAL AGE
#######################################################################################
DNAmAge.ChronoAge.plot = function(df.ttype, legname, colors, symbol.shapes = "", labels, x.label, 
                                  y.label, title, x.min = 0, x.max = 100, y.min = 0, 
                                  y.max = 100, leg.x = 0.8, leg.y = 0.87, text.size = 18){
  
  if(symbol.shapes == ""){ symbol.shapes <- c( rep(16, length(colors)) ) }
  
  p <- ggplot(df.ttype, aes(x = Chrono.age, y = DNAm.age, color = ttype, shape = ttype)) +
    geom_point() + 
    scale_color_manual(
      name = legname, 
      values = colors, 
      breaks = labels
    ) +
    scale_shape_manual( 
      name = legname,
      values = symbol.shapes
    ) +
    scale_x_continuous(
      expand=c(0, 0),
      limits = c(x.min, x.max)
    ) +
    scale_y_continuous(
      expand=c(0, 0),
      limits = c(y.min, y.max)
    ) + 
    labs(x = x.label) + 
    labs(y = y.label) + 
    labs(title = title) + 
    theme_bw(base_size = text.size) + 
    theme(legend.key.width = unit(3, "line"), legend.position=c(leg.x, leg.y)) +
    theme(
      axis.ticks.length=unit(-0.25, "cm"), 
      axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black", size = text.size), 
      axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black", size = text.size), 
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    ) 
  
  p
  
}

#######################################################################################
# AGE ACCELERATION HISOGRAM
#######################################################################################
accel.hist.plot = function(df.ttype, bw = 20, legname, linetypes, colors, labels, x.label, y.label, title,
                           annot.x = 90, annot.y = 0.03, annot.sep = 0.002, x.min = -40, x.max = 125, 
                           y.min = 0, y.max = 0.04, leg.x = 0.8, leg.y = 0.87, text.size = 18){

  
  types <- unique(df.ttype$ttype)
  median.list <- list()
  median.type <- list()
  i = 1
  for(t in types){
    res <- as.numeric(as.vector(df.ttype[df.ttype$ttype == t, "res"]))
    median.list[[i]] <- median(res)
    median.type[[i]] <- t
    i = i+1
  }
  
  p <- ggplot(df.ttype, aes(x = res, stat(density), color = ttype, linetype = ttype)) +
    geom_freqpoly(size = 1.2, binwidth = bw) +
    scale_linetype_manual(
      name = legname, 
      values = linetypes,
      labels = labels
    ) +
  scale_color_manual(
    name = legname, 
    values = colors, 
    labels = labels
  ) +
  scale_x_continuous(
    expand=c(0, 0),
    limits = c(x.min, x.max)
  ) +
  scale_y_continuous(
    expand=c(0, 0),
    limits = c(y.min, y.max)
  ) + 
  labs(x = x.label) + 
  labs(y = y.label) + 
  labs(title = title)
  
  i = 1
  for(med in median.list){
    p <- p +   
      annotate(
        "text", x = annot.x, y = (annot.y - (i-1)*annot.sep), 
        label = paste("Median Accel. ", "(", median.type[[i]], ") = ", round(med,  digits = 1), sep = ""), 
        color = colors[i]
      )
    i = i+1
  }
  p <- p + 
    geom_vline(xintercept = 0, color = "black", size = 1, linetype = "dotted") + 
    theme_bw(base_size = text.size) + 
    theme(legend.key.width = unit(3, "line"), legend.position=c(leg.x, leg.y)) +
    theme(
      axis.ticks.length=unit(-0.25, "cm"), 
      axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black", size = text.size), 
      axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black", size = text.size), 
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    ) 
  
  p

}


#######################################################################################
# AGE ACCELERATION BOX PLOT
#######################################################################################
accel.box.plot = function(df.ttype, residuals, width = 0.75, pos.col = "blue", neg.col = "red", 
                          x.label = "", y.label = "", title = "", show.leg = TRUE, leg.x = 0.25, 
                          leg.y = 0.87, text.size = 18){

  w = width
  df.list <- list()

  # Set up dataframes 
  i <- 1
  for(r in residuals){
    r.sorted <- sort(r)
    lower <- i - w/2
    upper <- i + w/2
    step  <- w/(length(r.sorted)-1)
    x <- seq(lower, upper, step)
    df <- data.frame(x, x, r.sorted)
    colnames(df) <- c("x0", "x1", "y1")
    df$y0 <- 0
    df$clr <- "Positive Acceleration"
    df[df$y1 < 0, "clr"] <- "Negative Acceleration"
    df <- df[ c("x0", "x1", "y0", "y1", "clr") ]
    df$clr <- factor( df$clr, levels = rev(unique(as.character(df$clr))) )
    df.list[[i]] <- df
    i <- i + 1
  }

  p <- ggplot(df.ttype, aes(x=ttype, y=res)) +
    geom_boxplot(aes(x=ttype, y=res), data = df.ttype, width = 1.05*width) +
    geom_hline(yintercept = 0, color = "gray", size = 1, linetype = "dashed")
  for( i in df.list ){
    p <- p + geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1, color = clr), data = i)
  }
  p <- p +   
  scale_linetype_manual(
    name = "", 
    values = c("solid", "solid"),
    labels = c(expression("EAA">=0), expression("EAA"<0))
  ) + 
  scale_color_manual(
    name = "", 
    values = c(pos.col, neg.col), 
    labels = c(expression("EAA">=0), expression("EAA"<0))
  ) +
  geom_boxplot(aes(x=ttype, y=res), data = df.ttype, width = 1.05*width, alpha = 0.2) +
  labs(x = x.label) + 
  labs(y = y.label) + 
  labs(title = title) + 
  theme_bw(base_size = text.size) +
  theme(
    legend.background = element_rect(color = "transparent", fill = "transparent"),
    legend.key.width = unit(1, "line"), 
    legend.position=c(leg.x, leg.y),
    axis.ticks.length=unit(-0.25, "cm"), 
    axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black", size = text.size), 
    axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"), color = "black", size = text.size), 
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) 
  
  if(!show.leg){ p <- p + theme(legend.position="none") }
  
  p
}