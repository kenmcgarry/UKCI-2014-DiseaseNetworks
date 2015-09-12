## UKCI2014_hdn.r     2/06/14
## load in files for disease and gene lists
## The deadline is 8th June 2014
library(igraph)
library(linkcomm)
library(bipartite)
library(NetIndices)
library(NCBI2R)
#library(networksis)
#library(network)

#----------------- GENE CARDS: diabetes and diseases related to diabetes ----------------------------
related <- read.csv('C:\\R-files\\disease\\GCARDS_diabetes_related_diseases.csv', header = TRUE)
#related <- as.matrix(related)     
#genesDia <- read.csv('C:\\R-files\\disease\\GCARDS_diabetes_related_genes.csv', header = TRUE)
#genesAlz <- read.csv('C:\\R-files\\disease\\GCARDS_alzheimers_related_genes.csv', header = TRUE)
#----------------------------------------------------------------------------------------------------

# Create the data structure to build a bipartite graph of all known diseases related to diabetes,
# will build a three column dataset [disease name;associated genes; score]
# THIS CODE NOW DEBUGGED.....Ken 31/05/2014
disease=10;score=100;genes=200
BPG<-data.frame(disease=disease,genes=genes,score=score) # must create the first disease-gene record
pindex<- 1;
for (i in 1:nrow(related)){
    genes<-unlist(strsplit(as.character(related[i,2]),";"))
    dStr<-as.character(related[i,1])
    scoreStr<-as.character(related[i,3])
    if(length(genes)>0){
      for (j in 1:length(genes)){
      geneStr  <-as.character(genes[j])
      BPG[pindex,]<-c(dStr,geneStr,scoreStr)
      pindex<-pindex+1
      }
    }  
}

# ------------------------------------------------------------------
# plot histo of gene scores 31/05/2014
score<-as.numeric(BPG[,3]) # convert from string to numeric
# histograms in R are rubbish without setting the axis - basically columns are bigger than axis!
hist(score,xlab="Protein scores",main=NULL,xlim = range(pretty(c(0,score))),ylim=c(0,3000),col = "lightgray")

# ------------------------------------------------------------------

# finally write results to a CSV file
write.csv(BPG,file="C:\\R-files\\disease\\bipartite_related.csv")

# ok, now get it back for conversion from edge-list to matrix
edge.list <- read.csv("C:\\R-files\\disease\\bipartite_related.csv",header=TRUE)
edge.list <- edge.list[,2:4] # we just need disease,genes , score columns
g <- graph.data.frame(edge.list, directed=FALSE)
#get.adjacency(g, type="both", attr="weight") 
am<-get.adjacency(g,type="both") 

ig <- graph.data.frame(edge.list)
V(ig)$type <- V(ig)$name %in% edge.list[,1] # uses matrix mode i.e. %in%
ig

subg <- subgraph.edges(ig, 1:25, 1:50)  # Nodes; Links
tkplot(subg,
  layout = layout.fruchterman.reingold,
  vertex.label = V(subg)$name,
  vertex.label.color= "black",
  ##vertex.size=nodesize,
  #vertex.colors=nodecolor,
  #vertex.shape=
  edge.arrow.size=0,
  edge.curved=FALSE
)

#------ search for other diseases with proteins in the diabetes type 2 group ---------
am<-get.adjacency(ig,sparse=F)
properties<-GenInd(am) # takes a while but calculates some metrics such as...
                                # No Nodes, No links, link density, 'connectance'

deg.distr<-degree.distribution(ig,cumulative=T,mode="all")
  #Degree distribution is the cumulative frequency of nodes with a given degree
    # this, like degree() can be specified as "in", "out", or "all"

plot(betweenness(am))
degree(am)
#------------- USING BIPARTITE PACKAGE TO PRINT OUT BIPARTITE GRAPHS --------
# plotweb(Safariland, abuns.type='independent',arrow="down.center")
# using the bipartite package function 'plotweb' we can actually plot the bipartite graph.
plotweb(am[1:20,],col.high=c("white","white"),text.rot=90, low.lablength=20, high.lablength=5)

#-----------
plotweb(am[1:20,],col.high=c("red","green"),text.rot = 45)
plotweb(am[21:40,],col.high=c("red","green"),text.rot = 45)

plotweb(am[1:30,],col.low=c("orange","green"),col.high=c("white","grey","purple"),
        text.high.col=c("blue","red"), col.interaction=c("red",rep("green",26),rep("brown",242)),
        bor.col.interaction=c(rep("green",26),rep("brown",242)),method="normal",
        text.rot=90, low.lablength=10, high.lablength=5)

# we also use visweb to plot a confusion-matrix like grid.
#visweb(am[21:40,],square="b",box.col="green",box.border="black")


# closeness from bipartite
plot(closeness(am))

# Network notation V = vertex or node, E = edge or connections
V(ig)[V(ig)$name == 'LRP5' ]
E(ig) [ from(1) ]

### Start to search the network for: 
# 1. node that contains diabetes
# 2. the genes linked with diabetes 
# 3. diseases that have these genes in common.
y<-V(ig)$name == 'type 2 diabetes mellitus'
n<-as.numeric(V(ig)[V(ig)$name == 'type 2 diabetes mellitus' ]) # Gives us the node (vertex) number 
                                                                # of diabetes

V(ig)[n]  # diabetes is node 1, just to prove it #V(ig)[2]  # carcinoma #V(ig)[3]  # pancreatitis
                                    #V(ig)[2]  # obesity

theGenes<- E(ig) [ from(n) ] # Get the genes associated with diabetes.
length(theGenes)           # How many genes?

for (i in 1:length(theGenes)){
  #print(theGenes[i])
  #print(V(ig)$name[i])
  
}


f<-V(ig)[V(ig)$name == 'PIK3R1' ]
E(ig) [ from(j) ]

V(ig)$name[905]  # This gets the vertex name as a string which is 'PIK3R1'
E(ig)[ V(ig)]   # Appears to list every edge
el<-get.edgelist(ig)


#------- using linkcomm package for further analysis of connectivity patterns ------
lc <- getLinkCommunities(edge.list[,1:2],hcmethod="average")
print(lc) # plot a summary for bullet list in paper
plot(lc,type="members")  # plots the community membership 'blockmap'

cc <- getCommunityCentrality(lc)
plot(cc)  # index v CC scatterplot plot

hc <- getClusterRelatedness(lc) # the nice looking cluster dendrogram

cm <- getCommunityConnectedness(lc,conn="modularity") # barchart of community
plot(lc,type="commsumm",summary="modularity")

#res <- computeModules(am) # takes several hours! Plots a grid style diagram. 
plotModuleWeb(res) # NB. plotmoduleweb() is a better function.

#----- NOW GET THE CLUSTERS -----------------------------------     

ad<-get.adjacency(g, names=TRUE, sparse=FALSE)

res <- computeModules(el) 

resA<-computeModules(am[1:550,], deep = FALSE, deleteOriginalFiles = TRUE,
                      steps = 500000, tolerance = 1e-10, experimental = FALSE)
     
resB<-computeModules(am[550:1229,], deep = FALSE, deleteOriginalFiles = TRUE,
                       steps = 500000, tolerance = 1e-10, experimental = FALSE)
 
#plotModuleWeb(resA)
#plotModuleWeb(resbig3)
     
# ---- compute CZ values and do the plot -----
cz <- czvalues(res)
plot(cz[[1]],cz[[2]],pch=16,xlab="protein coefficient",ylab="z-score",cex=0.8,xlim=c(0,1),las=1)
abline(v=0.62) # threshold of Olesen et al. 2007
abline(h=2.5) # dito
text(cz[[1]], cz[[2]], names(cz[[2]]), pos=4, cex=0.7)
     
#--------- start to delve deeper into pathways 
# however, we can only use entrez ids
paths1<-GetPathways(c(5294,55634,5286,4041))# PIK3CG 5294;ZNF673 55634;PIK3C2A 5286;LRP5 4041
paths2<-GetPathways(c(5294,5286,5972,7422,7450,56729)) # 
     
# d<-GetIDs("REN[sym] human")     
# GetGOs(5294)     
     
#-------------------------------------------------------------------------
# rm(list = ls())  # kill all variables in memory

subg <- subgraph.edges(g, 60:65, 1:50)  # Nodes; Links

## plot histogram, useful to see  distribution of nodes per neighbourhood
ns <- neighborhood.size(subg, 1)
hist(ns, xlab="Neighbourhood size of subgraph", main=NULL)

ns <- neighborhood.size(g, 1)
hist(ns, xlab="Neighbourhood size of full graph", main=NULL)

#-----------------------------------
# the adjacency matrix was 1229 x 1229 when in fact a matrix of 325 x 905 should be the case.
disease<-unique(edge.list[,1])
proteins<-unique(edge.list[,2])

#creates a matrix of 0s with the node IDs as rows and columns
mat=matrix(0,nrow=length(disease),ncol=length(proteins),dimnames=list(disease,proteins)) 


#for each row in the edgelist, find the appropriate cell in the empty matrix and assign it a 1
for (i in 1:nrow(edge.list)) {
  for (j in 1:nrow(diseases))
    mat[,j]<-1
  #mat[edge.list[i,1],edge.list[i,2]]<-1 #mat[edge.list[i,1],edge.list[i,2]]
}

# simple but correct?
z <- table(edge.list[,1:2]) # just the diaeses and proteins

y<-as.data.frame.matrix(z) 

date()
res <- computeModules(y[1:900,]) # takes a loooong time, in fact 900 did not complete.
date()

# do the nice cz plots, recall the cut-off points (0.62 & 2.5) are arbitary 
cz <- czvalues(res)
plot(cz[[1]],cz[[2]],pch=16,xlab="protein coefficient",ylab="z-score",cex=0.8,xlim=c(0,1),las=1)
abline(v=0.62) # threshold of Olesen et al. 2007
abline(h=2.5) # dito
text(cz[[1]], cz[[2]], names(cz[[2]]), pos=4, cex=0.7)

#------------ NOT USED , BUT CONSIDER USING FOR NEXT AROUND --------------
# Random bipartite graph from igraph example website
inc <- matrix(sample(0:1, 50, replace=TRUE, prob=c(2,1)), 10, 5)
g <- graph.incidence(inc)
plot(g, layout=layout.bipartite,vertex.color=c("green","cyan")[V(g)$type+1])
# Two columns
lay <- layout.bipartite(g)
plot(g, layout=lay[,2:1])
#---------------------------------------------------

# hist(score,xlab="Protein scores",main=NULL,xlim = range(pretty(c(0,score))),ylim=c(0,3000),col = "lightgray")
#----------------------------------------------------------------------------
# Plot the simple example bipartite graph based on the adjaceny matrix, Fig 3. in the paper.

edgelist <- read.table(text="Gene    Disease
                         G1  D2
                         G1  D5
                         G2  D1
                        G2  D3
                        G2  D4
                        G3  D2
                        G3  D4
                        G4  D1
                        G4  D3
                        G4  D4
                        G5  D1
                        G5  D4", 
                       header=TRUE)
ig <- graph.data.frame(edgelist)
V(ig)$type <- V(ig)$name %in% edgelist[,1]

nodecolor=character(nrow(edgelist))  # create a character for every column in adjaceny matrix,
nodecolor<-c(rep("lightblue",5),rep("pink",5))

tkplot(ig,
       layout = layout.fruchterman.reingold,
       vertex.label = V(ig)$name,
       vertex.label.color= "black",
       ##vertex.size=nodesize,
       vertex.color=nodecolor,
       #vertex.shape=
       edge.arrow.size=0,
       edge.curved=FALSE
)

# scrap.r   try out crazy ideas here!
disease<-unique(edge.list[,1])
proteins<-unique(edge.list[,2])

#creates a matrix of 0s with the node IDs as rows and columns
mat=matrix(0,nrow=length(disease),ncol=length(proteins),dimnames=list(disease,proteins)) 


#for each row in the edgelist, find the appropriate cell in the empty matrix and assign it a 1
for (i in 1:nrow(edge.list)) {
  for (j in 1:nrow(diseases))
    mat[,j]<-1
  #mat[edge.list[i,1],edge.list[i,2]]<-1 #mat[edge.list[i,1],edge.list[i,2]]
}

# simple but correct?
z <- table(edge.list[,1:2]) # just the diaeses and proteins

y<-as.data.frame.matrix(z) 

date()
res <- computeModules(y[1:900,]) 
date()


cz <- czvalues(res)
plot(cz[[1]],cz[[2]],pch=16,xlab="protein coefficient",ylab="z-score",cex=0.8,xlim=c(0,1),las=1)
abline(v=0.62) # threshold of Olesen et al. 2007
abline(h=2.5) # dito
text(cz[[1]], cz[[2]], names(cz[[2]]), pos=4, cex=0.7)


     