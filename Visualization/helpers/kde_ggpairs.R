kdeggpairs <- function(x, n_1d=20, n_2d=25, labels, density=TRUE,
                             contour=TRUE, ...) {
  require(MASS)
  require(sp)
  require(GGally)
  require(ggplot2)
  require(dplyr)
  require(tidyverse)
  
  par(cex.axis=0.95) # magnify the axes font
  
  
  fun.lower <- function(x1, x2, ...) {
    if (is.factor(x1)) x1 <- as.integer(x1)
    if (is.factor(x2)) x1 <- as.integer(x2)
    OK <- length(unique(x1)) > 2 && length(unique(x2)) > 2
    
    if (!density && !contour) n_2d <- 0
    
    if (n_2d > 0 && OK) {
      if (density || contour) {
        d <- tryCatch(kde2d(x1, x2, n=n_2d), error=function(e) {
          if (e$message == "bandwidths must be strictly positive") {
            cat("Warning: resetting the bandwidths")
            hx <- 0.1*(max(x1)-min(x1)) #had 0.1*
            hy <- 0.1*(max(x2)-min(x2))
            d <- kde2d(x1, x2, n=n_2d, h=c(hx, hy))
          } else cat("Error caught! Handling not implemented!")
        })
      }
      
      if (density) {
        # library("colorspace"); vv<-c("white", heat_hcl(12))
        # library(RColorBrewer); vv<-c("white", rev(brewer.pal(9, "YlOrRd")))
        #vv<-c(rev(bpy.colors(n=10, cutoff.tails=0.3, alpha=0.7))) #"white", 
        vv<-c( bpy.colors(n=100, cutoff.tails=0.3, alpha=0.8)) # 
        
        image(d, col=vv, add=TRUE)
      }
      
      if (contour) graphics:::contour(d, add=TRUE, nlevels=8) #nlevels=8
    } else points(x1, x2)
  }
  
  fun.upper <- function(x1, x2, ...) {
    if (is.factor(x1)) x1 <- as.integer(x1)
    if (is.factor(x2)) x1 <- as.integer(x2)
    
    vv<-bpy.colors(n=100, cutoff.tails=0.3, alpha=0.2)
    
    # This adds a column of color values based on the last column values
    #datc <- vv[as.numeric(cut(as.numeric(loglh), breaks=100))]
    datc <- vv[cut(loglh, breaks=100)]
    points(x1, x2, col=datc, pch=16, cex=0.5)
  }
  

    nd <- dim(x)[2]
    #loglh <- x[, nd]
    #pairs(x[, 1:(nd-2)], labels=labels, lower.panel=fun.lower,
    #      upper.panel=fun.upper, diag.panel=fun.diag.kde)
    # was x[,1:(nd-1)]
    
    data<-data.frame(x)

    ci_eti <- ci(data[,1:(nd-2)],ci=0.95, method = "HDI",verbose = TRUE)
    
    
    plotList <- list()
    for(i in 1:(nd-2))
    {
      for(j in 1:(nd-2))
      {
        plotId = i + (nd-2) * (j-1)
        
        if(i==j)
        {
          plotList[[plotId]] <- data[,i] %>% 
            estimate_density(extend=TRUE) %>%
            ggplot(aes(x = x, y = y)) +
            geom_area(fill = "orange") +
            theme_bw() +
            # Quantile in red
            geom_vline(xintercept = ci_eti$CI_low[i], color = "red", size = 1) +
            geom_vline(xintercept = ci_eti$CI_high[i], color = "red", size = 1)+
            theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            annotate("text", x=-Inf, y = Inf, label = l[i], parse = TRUE,  hjust = -0.2, vjust = 1.2, size=10)
            
        }
        else if(i>j)
        {
          d<-data.frame(data[,i],data[,j])
          x<-colnames(d)[1]
          y<-colnames(d)[2]
          plotList[[plotId]] <- ggplot(d, 
                                       aes_(x=as.name(x),y=as.name(y))) +
                            geom_point(aes(colour = -data[,nd-1],alpha = 0.01))+
                            geom_smooth(method = "lm")+
                            theme_bw() +
                            scale_x_continuous(position = "top") + 
                            scale_y_continuous(position = "right")

        }
        else if(i<j)
        {
          d<-data.frame(data[,i],data[,j])
          x<-colnames(d)[1]
          y<-colnames(d)[2]
          plotList[[plotId]] <- ggplot(d, 
                                       aes_(x=as.name(x),y=as.name(y))) +
                                       geom_density_2d_filled() + theme_bw() +
                                       scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
        }
      }
    }
    
    pm <- ggmatrix(plotList, nd-2, nd-2) + theme(axis.text=element_text(size=14-(nd-2)), #change font size of axis text
                                                 axis.title=element_text(size=16),axis.text.x = element_text(angle = 90)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3))# +
      #scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
    
    return(pm)
    
}
