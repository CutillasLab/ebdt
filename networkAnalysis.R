

# Author: Pedro R. Cutillas 2018


library(tidyverse)
library(tidygraph)
library(igraph)

# 1. get enrichbmet values of edges from folds using Kinase.Substrate.Enrichment

# 2. the output from Kinase.Substrate.Enrichment is input to Network.Analysis.From.Edge.Enrichment
# which returns centrality values

getNetwork <- function(df){
  k1 <- df %>%
    distinct(k1) %>%
    rename(label = k1)
  # Get distinct destination names
  k2 <- df %>%
    distinct(k2) %>%
    rename(label = k2)
  # Join the two data to create node
  # Add unique ID for each country

  nodes <- full_join(k1, k2, by = "label")
  nodes <- nodes %>%
    mutate(id = 1:nrow(nodes)) %>%
    select(id, everything())
  ##############
  # (a) Join nodes id for source column
  edges <- df %>%
    left_join(nodes, by = c("k1" = "label")) %>%
    rename(from = id)
  # (b) Join nodes id for destination column
  edges <- edges %>%
    left_join(nodes, by = c("k2" = "label")) %>%
    rename(to = id)
  # (c) Select/keep only the columns from and to
  edges <- select(edges, from, to, weight)
  net.igraph <- graph_from_data_frame(
    d = edges, vertices = nodes,
    directed = FALSE)
  return(net.igraph)
}


network.plots <- function(edges, mynetwork, mytitle){
  kk <- mytitle
  pp1 <- ggplot() +
    geom_segment(data=edges,aes(x=X1, y=Y1, xend = X2, yend = Y2),
                 size = 0.5, colour="grey") +
    geom_point(data=mynetwork, aes(x=X1, y=X2, color=betweenness,
                                   size=strength) )+
    ggrepel::geom_text_repel(data=mynetwork,
                             aes(x=X1, y=X2, label=kinase, size=strength))+
    theme_void()+
    scale_color_gradient2(low = "grey",mid = "orange",high = "red")+
    ggtitle(kk)+
    theme(plot.title = element_text(hjust = 0.5, color="blue"))

  pp2 <- ggplot(data=mynetwork,aes(x=strength,y=log10(betweenness)))
  pp2 <- pp2+geom_point(aes(color=degree, size=strength))+
    ggrepel::geom_text_repel(aes(label=kinase, size=strength))+
    theme_light()+
    scale_color_gradient2(low = "grey",mid = "orange",high = "red")+
    ggtitle(kk)+
    theme(plot.title = element_text(hjust = 0.5, color="blue"))

  return(list(pp1,pp2))
}


Network.Analysis.From.Edge.Enrichment <- function (df.edges, fcutoff, drawPlots=0){

  # df.edges is resutls of KSEA using nodes.edges database
  # fcutoff is recomended at 0.2
  # drawPlots, if == 1 , it will draw networks and will plot node strength comparisons
  ## returns: df.deg, #1
  #           df.bet,#2
  #           df.strength,#3
  #           df.index, #4
  #           plotList1, #5
  #           plotList2, #6
  #           plotList3, #7
  #           plotList4 #8
  ##
  df.edges <- scale(df.edges)

  tt1 <- read.table(text = rownames(df.edges), sep = ".",
                    colClasses = "character")
  colnames(tt1) <- c("k1","k2")
  df.edges <- cbind(tt1,df.edges)
  nc <- ncol(df.edges)
  df.kinases <- data.frame(table(rbind(df.edges$k1, df.edges$k2)))
  colnames(df.kinases) <- c("kinase","freq")
  df.deg <- data.frame(kinase=character(nrow(df.kinases)),
                       matrix(nrow = nrow(df.kinases), ncol=nc),stringsAsFactors = F)
  df.bet <- data.frame(kinase=character(nrow(df.kinases)),
                       matrix(nrow = nrow(df.kinases), ncol=nc),stringsAsFactors = F)
  df.strength <- data.frame(kinase=character(nrow(df.kinases)),
                            matrix(nrow = nrow(df.kinases), ncol=nc),stringsAsFactors = F)

  df.index <- data.frame(kinase=character(nrow(df.kinases)),
                         matrix(nrow = nrow(df.kinases), ncol=nc),stringsAsFactors = F)

  df.deg$kinase <- df.kinases[,1]
  df.bet$kinase <- df.kinases[,1]
  df.index$kinase <- df.kinases[,1]
  df.strength$kinase <- df.kinases[,1]

  rownames(df.deg) <- df.kinases[,1]
  rownames(df.bet) <- df.kinases[,1]
  rownames(df.index) <- df.kinases[,1]
  rownames(df.strength) <- df.kinases[,1]

  kinases <- df.kinases[,1]
  i=1
  x=7
  plotList1 <- list()
  plotList2 <- list()
  plotList3 <- list()
  plotList4 <- list()
  for (x in 3:nc){
    kk <- colnames(df.edges)[x]
    sh.s <- df.edges[,c(1,2,x)]
    colnames(sh.s) <- c("k1","k2","weight")
    colnames(df.bet)[x-1] <- kk
    colnames(df.deg)[x-1] <- kk
    colnames(df.index)[x-1] <- kk
    colnames(df.strength)[x-1] <- kk
    #sh.s$weight <- as.numeric(sh.s$weight) *(-1)
    if (sum(sh.s$weight, na.rm = T)!=0){
      #fcutoff <- 0.2#median(sh.s$weight, na.rm = T)+0.1
      sh.s <- subset(sh.s,weight>fcutoff & sh.s$weight !="NA")
      net <- getNetwork(na.omit(sh.s))
      d=degree(net) # report degree
      #V(net)$size= d#10
      ww <- 1/as.numeric(E(net)$weight)#*2
      E(net)$width <- ww#(max(sh.s$weight)/2)
      dd <- betweenness(net, directed = F ,weights = 1/ E(net)$width ) # betwennness for all nodes
      mystr <- igraph::strength(net,weights = 1/E(net)$width) # strength for all nodes
      mynetwork <- data.frame(V(net)$label,d,dd,mystr)
      colnames(mynetwork) <- c("kinase","degree","betweenness","strength")
      mynetwork$cen.index <- d*(dd)
      r=1
      for (r in 1:nrow(mynetwork)){
        k <- as.character(mynetwork$kinase[r])
        df.deg[k,x-1] <- mynetwork$degree[r]
        df.bet[k,x-1] <- mynetwork$betweenness[r]
        #df.bet[k,1] <- k
        df.index[k,x-1] <- mynetwork$cen.index[r]
        df.strength[k,x-1] <- mynetwork$strength[r]
      }
      if (drawPlots==1){
        plotcord <- data.frame(layout.fruchterman.reingold(net))
        edgelist <- get.edgelist(net)
        node.labels <- V(net)$label
        mynetwork1 <- cbind(mynetwork,plotcord)
        edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
        colnames(edges) <- c("X1","Y1","X2","Y2")

        myplots1 <- network.plots(edges,mynetwork1,kk)
        plotList1[i] <- plot(myplots1[[1]])
        plotList2[i] <-plot(myplots1[[2]])

        #append(plotList1, myplots2[[1]])
        #append(plotList2, myplots2[[2]])
        ##############################################################
        ## Reverse network
        sh.s2 <- df.edges[,c(1,2,x)]
        colnames(sh.s2) <- c("k1","k2","weight")
        sh.s2$weight <- sh.s2$weight * -1
        fcutoff =0.2#<- median(sh.s2$weight, na.rm = T)+0.1
        sh.s2 <- subset(sh.s2,weight>fcutoff & sh.s2$weight !="NA")
        net <- getNetwork(na.omit(sh.s2))
        d2=degree(net2)
        ww2 <- 1/as.numeric(E(net2)$weight)#*2
        E(net2)$width <- ww2
        dd2 <- betweenness(net2, directed = F ,weights = 1/ E(net2)$width )
        mystr2 <- igraph::strength(net2,weights = 1/E(net2)$width)
        mynetwork2 <- data.frame(V(net2)$label,d2,dd2,mystr2)
        colnames(mynetwork2) <- c("kinase","degree","betweenness","strength")
        mynetwork2$cen.index <- d2*log10(dd2)


        plotcord2 <- data.frame(layout.fruchterman.reingold(net2))
        edgelist2 <- get.edgelist(net2)
        node.labels2 <- V(net2)$label
        mynetwork22 <- cbind(mynetwork2,plotcord2)
        edges2 <- data.frame(plotcord2[edgelist2[,1],], plotcord2[edgelist2[,2],])
        colnames(edges2) <- c("X1","Y1","X2","Y2")

        myplots2 <- network.plots(edges2,mynetwork22,paste(kk, "Control") )
        plotList3[i] <-  plot(myplots2[[1]])
        #append(plotList3, myplots2[[1]])
        #mynetwork1$type <- "Case"
        # mynetwork22$type <- "Control"
        net.comp.stren <- data.frame(kinase =kinases ,
                                     case.s=numeric(nrow(df.kinases)),
                                     control.s=numeric(nrow(df.kinases)),
                                     case.b=numeric(nrow(df.kinases)),
                                     control.b=numeric(nrow(df.kinases)),
                                     stringsAsFactors = F)
        #k <- "ERN1"
        ii=1
        for (k in kinases){
          case.b <- mynetwork1[mynetwork1$kinase==k,"betweenness"]
          control.b <- mynetwork1[mynetwork22$kinase==k,"betweenness"]
          case.s <- mynetwork1[mynetwork1$kinase==k,"strength"]
          control.s <- mynetwork1[mynetwork22$kinase==k,"strength"]
          net.comp.stren[ii,1] <- k
          tryCatch({net.comp.stren[ii,2] <- case.s}, error=function(e){})
          tryCatch({net.comp.stren[ii,3] <- control.s}, error=function(e){})
          tryCatch({net.comp.stren[ii,4] <- case.b}, error=function(e){})
          tryCatch({net.comp.stren[ii,5] <- control.b}, error=function(e){})
          ii=ii+1
        }

        net.comp.stren[is.na(net.comp.stren)] <- 0
        net.comp.stren$delta.stren.case.control <- net.comp.stren$case.s-net.comp.stren$control.s
        net.comp.stren$delta.bet.case.control <- net.comp.stren$case.b-net.comp.stren$control.b

        xx <- 0
        if (xx==1){
          pp5 <- ggplot(net.comp.stren,aes(y=case.s+control.s,
                                           x=delta.stren.case.control))
          pp5+geom_text(aes(label=kinase))+ggtitle(kk)

          pp6 <- ggplot(net.comp.stren,aes(y=case.b+control.b,
                                           x=delta.bet.case.control))
          pp6+geom_text(aes(label=kinase))+ggtitle(kk)
          pp7 <- ggplot(net.comp.stren,aes(x=delta.bet.case.control,y=delta.stren.case.control, label=kinase))
          pp7+geom_point()+ ggrepel::geom_text_repel()+ggtitle(kk)

        }
        df.ss <- subset(net.comp.stren, delta.stren.case.control>1|
                          delta.stren.case.control<(-1))
        pp8 <- ggplot(df.ss, aes(x=reorder(kinase,
                                           delta.stren.case.control),
                                 y=delta.stren.case.control))
        pp8 <- pp8+geom_point()+coord_flip()+
          theme_bw()+
          geom_hline(yintercept = 0)+
          ggtitle(paste(kk, "differences in strength case - control"))
        plotList4[i] <- plot(pp8)
      }
      i=i+1
    }
  }

  return(list(df.deg, #1
              df.bet,#2
              df.strength,#3
              df.index, #4
              plotList1, #5
              plotList2, #6
              plotList3, #7
              plotList4))

}


get.protein.set.data <- function(dataset){

  ## data sets
    edges <- "https://www.dropbox.com/s/ttmzd40mnjgh1iu/edges.csv?dl=1"
    pdts <- "https://www.dropbox.com/s/86jfnayv0qa1n2q/pdts.csv?dl=1"
    psite <- "https://www.dropbox.com/s/eb1qoofz793f4tq/psite.csv?dl=1"
    signor <- "https://www.dropbox.com/s/alpbq880emz1z2t/signor.csv?dl=1"
    reactome <- "https://www.dropbox.com/s/jdcc1355cz73mmi/reactome.csv?dl=1"
    process <- "https://www.dropbox.com/s/z8zef96q45vi0je/process.csv?dl=1"
    myfunction <- "https://www.dropbox.com/s/ev50zpd0g41pds3/function.csv?dl=1"
    location <- "https://www.dropbox.com/s/js8vwyamlhbnqqf/location.csv?dl=1"
    nci <- "https://www.dropbox.com/s/fe8t4nyhbljsn5y/nci.csv?dl=1"
    human.cell.markers.full <- "https://www.dropbox.com/s/a2tvvuh8kgi340v/human_cell_markers_full_lineage.csv?dl=1"
    human.cell.markers.short <- "https://www.dropbox.com/s/4njofwfs5uu3wya/human_cell_markers_short_lineage.csv?dl=1"
    mouse.cell.markers <- "https://www.dropbox.com/s/223bpbc5p24j919/mouse_cell_markers_full_lineage.csv?dl=1"
    blood.cell.markers.human1 <- "https://www.dropbox.com/s/kl2apjfnbqdw7gj/blood_cell_markers.csv?dl=1"
    blood.cell.markers.human2 <- "https://www.dropbox.com/s/o28o4fbmmuwtgs3/blood_cell_markers2.csv?dl=1"
    bone.marrow.cell.markers <- "https://www.dropbox.com/s/gl0vim65cz9804a/bone_marrow_cell_markers.csv?dl=1"
    tf.targets <- "https://www.dropbox.com/s/34w7e8ute0w21n5/tf_targets.csv?dl=1"

    chromatin <- "https://www.dropbox.com/s/y2gkjn45oxy72nw/chromatin.csv?dl=1"

    selected <- "https://www.dropbox.com/s/nrraoj8hlvzap87/selected.csv?dl=1"

    dataset.names <- c("edges","pdts","psite","signor","reactome","process","function","location",
                       "nci","human.cell.markers.full",
                       "human.cell.markers.short","mouse.cell.markers",
                       "blood.cell.markers.human1",
                       "blood.cell.markers.human2",
                       "bone.marrow.cell.markers",
                       "tf.targets", "chromatin",
                       "selected")
    datasets <- c(edges,pdts,psite,signor,reactome,process,myfunction,location,nci,human.cell.markers.full,
                       human.cell.markers.short,mouse.cell.markers,
                  blood.cell.markers.human1,
                  blood.cell.markers.human2,
                  bone.marrow.cell.markers,
                  tf.targets, chromatin,selected)

    df.datasets <- data.frame(dataset.names,datasets)

    myfile <- as.character(df.datasets[df.datasets$dataset.names==dataset,2])
    if (length(myfile)==0){
      print (paste(dataset, "not found. Check spelling"))
    }
    df.out <- read.csv(myfile)
    return(df.out)
}

Kinase.Substrate.Enrichment <- function(df.fold, ks_db){



  # df.fold == dataset of fold changes
  # ks_db == database of kinase-substrate relationships
  #       possibilities are "edges", "ctams" "pdts", "pSite"
  #
  # returns 4 data frames:
  #         results.zscores,
  #         results.pvalues,
  #         results.m,
  #         results.q
  #         results.sites

  #mydir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  #mydir <- "C:/Users/cutill01/Dropbox/01pedro work/03_MACROS/KSEAR"
  ##mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
  #myfile <- paste(mydir,"/",ks_db,".csv",sep = "")

  nc <- ncol(df.fold)



  df.ks <- get.protein.set.data(ks_db)



  nr <- nrow(df.ks)
  #print(nr)
  results.pvalues <-data.frame(matrix(nrow = nr, ncol = nc-1))
  rownames(results.pvalues) <- make.names(df.ks[,1], unique = T)
  colnames(results.pvalues) <- colnames(df.fold)[2:nc]

  results.zscores <-data.frame(matrix(nrow = nr, ncol = nc-1))
  rownames(results.zscores) <- make.names(df.ks[,1],unique = T)
  colnames(results.zscores) <- colnames(df.fold)[2:nc]


  results.q <-data.frame(matrix(nrow = nr, ncol = nc-1))
  rownames(results.q) <- make.names(df.ks[,1], unique=T)
  colnames(results.q) <- colnames(df.fold)[2:nc]


  results.m <-data.frame(matrix(nrow = nr, ncol = nc-1))
  rownames(results.m) <- make.names(df.ks[,1],unique = T)
  colnames(results.m) <- colnames(df.fold)[2:nc]

  results.sites <-data.frame(matrix(nrow = nr, ncol = 2))
  rownames(results.sites) <- make.names(df.ks[,1],unique = T)


  df.selected.folds <- data.frame(matrix(nrow=1,ncol=nc+1))
  # df.selected.pvalues <- data.frame(matrix(nrow=1,ncol=nc+1))
  # df.selected.fdr <- data.frame(matrix(nrow=1,ncol=nc+1))

  colnames(df.selected.folds) <- c(colnames(df.fold),
                                   "prot.group")

  df.selected.folds$pathway.pvalue <- 1
  df.selected.folds$pathway.zscore <- 0
  df.selected.folds$pathway.m <- 0
  df.selected.folds$pathway.q <- 0
  df.selected.folds$pathway.median <- 0
  df.selected.folds$all.median <- 0
  df.selected.folds$all.sd <- 0

  #  colnames(df.selected.pvalues) <- c(colnames(df.pvalues),
  #             "prot.group")
  # colnames(df.selected.fdr) <- c(colnames(df.fdr),
  #         "prot.group")

  all.sites <- data.frame(matrix(nrow=1,ncol=2))
  colnames(all.sites) <- c("site","kinase")

  r=1
  for (r in 1:nr) {
    mym <- df.ks[r,2]
    kinase <- df.ks[r,1]
    if (is.na(mym) == F){
      if(mym>2){
        substrates <- as.character(df.ks[r,3])
        ss <- c(unlist(strsplit(substrates,";")))
        start.time <- Sys.time()
        df.xx <- subset(df.fold,df.fold[,1] %in% paste(ss,";",sep=""))

        # df.xxp<- subset(df.pvalue,df.pvalue[,1] %in% paste(ss,";",sep=""))
        # df.xxf<- subset(df.fdr,df.fdr[,1] %in% paste(ss,";",sep=""))


        sites <-paste(unlist(df.xx[,1]),collapse=";")

        end.time <- Sys.time()

        #print (paste(kinase, round(end.time-start.time, 2),"sec"))
        if (nrow(df.xx)>2){
          df.xx$prot.group <- kinase
          #df.xxp$prot.group <- kinase
          #df.xxf$prot.group <- kinase

          sites.x <- data.frame(site=df.xx[,1])
          sites.x$kinase <- kinase
          all.sites <- rbind(all.sites,sites.x)



          c=3
          for (c in 2:nc){
            values.all <- as.numeric(subset(df.fold[,c], df.fold[,c]!=0))
            myvalues <- as.numeric(subset(df.xx[,c],df.xx[,c]!=0))
            #hist(myvalues)
            #hist(av.all)
            pval <- 1
            tryCatch(
              myks <- ks.test(values.all,myvalues)
              , error=function(e){}
            )
            pval <- myks$p.value
            m <- 0
            q <- 0
            m <- nrow(df.xx)
            mysd <- sd(values.all)
            mymedian <- median(myvalues)
            mymedian.all <- median(values.all)
            sd.all <- sd(values.all)
            zscore <- ((mymedian-mymedian.all)*((sqrt(m)))/mysd)
            q <- sum(myvalues>(mymedian.all+1))+sum(myvalues<(mymedian.all-1))

            results.pvalues[r,c-1] <- pval
            results.zscores[r,c-1] <- zscore
            results.m[r,c-1] <- m
            results.q[r,c-1] <- q
            results.sites[r,1] <- as.character(sites)

            ###############################
            df.xx$pathway.pvalue <- pval
            df.xx$pathway.zscore <- zscore
            df.xx$pathway.m <- m
            df.xx$pathway.q <- q
            df.xx$pathway.median <- mymedian
            df.xx$all.median <- mymedian.all
            df.xx$all.sd <- mysd



            #df.selected.folds <- rbind(df.selected.folds,df.xx)
            #df.selected.pvalues <- rbind(df.selected.pvalues,df.xxp)
            #df.selected.fdr <- rbind(df.selected.fdr,df.xxf)

          }
        }
      }
    }
  }

  return(list(results.zscores,results.pvalues, results.m, results.q,
              results.sites,
              all.sites))

  ##################################
}
