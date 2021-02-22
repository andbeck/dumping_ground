library(igraph)
library(NetIndices)
library(patchwork)

source("new.adbm.r")
source("Plot.matrix.r")

input <- 
  data.frame(
    # sampling body sizes
    ref_seed = 1,
    # species richness
    ref_num_S = 20,
    # mean Body Mass
    ref_mean_BM = 10,
    # SD body mass (log)
    ref_sd_log_BM = 5,
    # log 10 attack rate mass scaling CONSTANT
    ref_a = -2,
    # Prey Mass attach rate EXPONENT
    ref_ai = 0.5,
    # Predator Mass attach rate EXPONENT
    ref_aj = 0.5,
    # Handling Time CONSTANT
    ref_r.a = 1,
    # Critical Predator Prey Mass Ratio (optimal in terms of profit of prey)
    # value of 2 makes most profitable item same size as predator
    ref_r.b = 2,
    # Prey Energy Content CONSTANT
    ref_e = 1,
    # Prey Energy Content EXPONENT
    ref_ei = 1,
    # Prey Abundance Mass Scaling CONSTANT
    ref_n = 1,
    # Prey Abundance Mass Scaling EXPONENT
    ref_ni = -0.75)


makeADBM <- function(input) {
  set.seed(input$ref_seed)
  M <- sort(rlnorm(input$ref_num_S, input$ref_mean_BM, input$ref_sd_log_BM))

  EHL <- Ratio.allometric.EHL(M=M,
                                e=input$ref_e,
                                r.a=input$ref_r.a, r.b=input$ref_r.b,
                                a=10^input$ref_a, ai=input$ref_ai, aj=input$ref_aj,
                                n=input$ref_n, ni=input$ref_ni)
  webout <- Get.web(EHL)
  connectance <- round(sum(webout/input$ref_num_S^2),2)
  output <- list(web = webout, connectance = connectance, Masses = M)
  return(output)
}

# use fun
ADBM_out <- makeADBM(input)


## Calculate trophic level... funky when there is no strict hierarchy
makeIgraph <- function(web){
  TL <- TrophInd(ADBM_out$web)$TL-1
    oo <- order(TL) ## save order, so we can later reorder the food web
    ## Chunk to make x_pos with width according to num spp per TL
    d_TL <- round(TL)
    spp_per_TL <- table(d_TL)
    x_pos <- rep(NA, length=input$ref_num_S)
    st <- cumsum(spp_per_TL)
    x_pos[1:st[1]] <- 1:st[1] - mean(1:st[1])
    for(i in 2:length(spp_per_TL))
      x_pos[(st[i-1]:(st[i]-1))+1] <- st[i-1]:(st[i]-1)-mean(st[i-1]:(st[i]-1))
    ## make and fill the layout matrix
    lay<-matrix(nrow=input$ref_num_S, ncol=2)
    #lay[,1] <- 1:input$num_S
    lay[,1] <- x_pos
    lay[,2] <- (TrophInd(web)$TL-1)
    ## reorder the web to correspond to layout
    oo_web <- web[oo,oo]
    ## convert web to igraph object
    gg <- graph_from_adjacency_matrix(web)
    plot.igraph(gg,layout=lay,
                vertex.label=NA,vertex.size=40,
                edge.arrow.size=.5,edge.width=.5,
                rescale=FALSE,
                ylim=c(0,20), xlim=c(-10,10),
                frame=TRUE)
    return(d_TL)
}

makeIgraph(ADBM_out$web)
Plot.matrix(ADBM_out$web)
ADBM_out$connectance
ADBM_out$Masses


library(tidyverse)
PPMR_df <- data.frame(round(TrophInd(ADBM_out$web)), Masses = ADBM_out$Masses)

PPMR_df %>% filter(TL == 1|TL == 2) %>% 
