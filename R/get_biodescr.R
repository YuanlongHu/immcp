#' Extract Biological descriptor
#'
#' @rdname extr_biodescr
#' @title Extract Biological descriptor
#' @param BasicData BasicData object.
#' @param geneset Charactor vector, one of "kegg"(KEGG), "mkegg"(KEGG Module), "go"(GO-BP), and "wp"(WikiPathways); a data frame and list.
#' @param arguments A list of the arguments of `clusterProfiler`, including `minGSSize`, `maxGSSize`, `pvalue`, and `qvalue`.
#' @param ref_type Charactor vector, one of "drug", "herb", "compound" or "target", defaults to "drug".
#' @param ref Charactor vector, reference drug, herb, compound or target, defaults to `NULL`.
#' @param to_ENTREZID Logical, whether to translate to ENTREZID from SYMBOL, defaults to TRUE.
#' @return A BioDescr object.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @exportMethod extr_biodescr
#' @examples
#' \dontrun{
#' data(drugdemo)
#' drug_herb <- PrepareData(drugdemo$drug_herb, from = "drug", to="herb")
#' herb_compound <- PrepareData(drugdemo$herb_compound, from = "herb", to="compound")
#' compound_target <- PrepareData(drugdemo$compound_target, from = "compound", to="target")
#' disease <- PrepareData(drugdemo$disease, diseaseID = "disease",from = "target", to="target")
#' BasicData <- CreateBasicData(drug_herb, herb_compound, compound_target, diseasenet = disease)
#' biodescr <- extr_biodescr(BasicData, geneset= "kegg")
#' }

setMethod("extr_biodescr", signature(BasicData = "BasicData"),
          function(BasicData, geneset= c("kegg", "mkegg","go","wp"),
                   arguments= list(minGSSize=5,maxGSSize=500,
                                   pvalue=0.05,qvalue=0.1),
                   ref_type="drug", ref=NULL, to_ENTREZID = TRUE) {

            message("<<<<< Prepare >>>>>")

            if (!is.null(ref)) {
              ref <- BasicData@vertices %>% filter(.data$type==ref_type & .data$name %in% ref)
            }else{
              ref <- BasicData@vertices %>% filter(.data$type==ref_type)
            }

            Druglist <- lapply(as.list(ref$name), function(x){
              subgraph <- subset_network(BasicData, from=x)
              node <- subgraph@vertices
              node <- node$name[node$type=="target"]
              return(node)
            })
            names(Druglist) <- ref$name

            res <- extr_biodescr_internal(Druglist = Druglist,
                                          DisBiomarker = BasicData@biomarker,
                                          geneset= geneset[1],
                                          arguments= arguments,
                                          ref_type=ref_type, ref=ref,
                                          to_ENTREZID = to_ENTREZID)
            return(res)
          })


#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom magrittr %<>%
#' @importFrom igraph union
#' @importFrom igraph as_data_frame
#' @importFrom pbapply pblapply
#' @importFrom rlang .data


extr_biodescr_internal <- function(Druglist, DisBiomarker,
                                  geneset= c("kegg", "mkegg","go", "wp"),
                                  arguments= list(minGSSize = 5, maxGSSize = 500,
                                                    pvalue = 0.05, qvalue = 0.1),
                                  ref_type = "drug", ref = NULL, to_ENTREZID = TRUE){

  message("<<<<< Extract >>>>>")
  message("> Extract biological descriptor of Disease")

  # "DisBiomarker" is a list of disease-related gene.
  pathway_disease <- pblapply(DisBiomarker, function(x){
    res <- enrich_f(x, geneset=geneset,
                    arguments=arguments,
                    to_ENTREZID=to_ENTREZID)
    res <- res[,c("ID", "geneID", "Description")]
    return(res)
    })
  message("> Extract biological descriptor of Drug")
  pathway_drug <- pblapply(Druglist, function(x){
    res <- enrich_f(x, geneset=geneset,
                    arguments=arguments,
                    to_ENTREZID=to_ENTREZID)
    res <- res[,c("ID", "geneID", "Description")]
    return(res)
  })

  message("<<<<< Output >>>>>")
  # a igraph object
  # from drug to geneset ID
  message("> Output Disease Graph")
  pathway_disease1 <- lapply(pathway_disease, function(x) x$ID)
  pathway_drug1 <- lapply(pathway_drug, function(x) x$ID)
  message("> Output Drug Graph")
  pathway_disease1 <- to_df(pathway_disease1)
  names(pathway_disease1) <- c("from", "to")
  pathway_drug1 <- to_df(pathway_drug1)
  names(pathway_drug1) <- c("from", "to")

  v1 <- rbind(
    data.frame(name=unique(pathway_disease1$from), type="disease"),
    data.frame(name=unique(pathway_disease1$to), type="pathway"),
    data.frame(name=unique(pathway_drug1$from), type="drug"),
    data.frame(name=unique(pathway_drug1$to), type="pathway")
  ) %>% distinct()

  pathway <- graph.data.frame(rbind(pathway_disease1, pathway_drug1), directed = TRUE, vertices = v1)

  # a list of igraph object
  # from geneset ID to gene, for each disease
  pathway_disease2 <- lapply(pathway_disease, function(x){
    res <- to_list(x, input = "basket", sep="/") %>%
      to_df()
    names(res) <- c("from", "to")
    res <- graph.data.frame(res, directed = TRUE)
    return(res)
  })
  # a list of igraph object
  # from geneset ID to gene, for each drug
  pathway_drug2 <- lapply(pathway_drug, function(x){
    res <- to_list(x, input = "basket", sep="/") %>%
      to_df()
    names(res) <- c("from", "to")
    res <- graph.data.frame(res, directed = TRUE)
    return(res)
  })
  IDtoDescription <- lapply(c(pathway_disease,pathway_drug), function(x) x[,c("ID", "Description")])
  IDtoDescription <- Reduce(rbind, IDtoDescription) %>% distinct()

  res <- new("BioDescr",
             drug_geneset = pathway,
             geneset_gene = list(disease = pathway_disease2,
                                 drug = pathway_drug2),
             anno = IDtoDescription
  )
  return(res)
}


#' Enrich Analysis
#'
#'
#' @title enrich_f
#' @param target_character Charactor vector of gene.
#' @param geneset Charactor vector, one of "kegg"(KEGG), "mkegg"(KEGG Module), "go"(GO-BP), and "wp"(WikiPathways); a data frame and list.
#' @param arguments A list of the arguments of `clusterProfiler`, including `minGSSize`, `maxGSSize`, `pvalue`, and `qvalue`.
#' @param out_dataframe Logical, whether to output data frame,defaults to `FALSE`.
#' @param to_ENTREZID Logical, whether to translate to ENTREZID from SYMBOL, defaults to `TRUE`.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom clusterProfiler enricher
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler enrichMKEGG
#' @importFrom clusterProfiler enrichGO
#' @importFrom clusterProfiler enrichWP
#' @importFrom clusterProfiler bitr
#' @importFrom utils data
#' @importFrom DOSE setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom rlang .data
#' @return data frame
#' @export
#' @author Yuanlong Hu

enrich_f <- function(target_character,
                     geneset = c("kegg", "mkegg","go", "wp"),
                     arguments=list(minGSSize=5,maxGSSize = 500,
                                    pvalue=0.05,qvalue=0.1),
                     out_dataframe=TRUE, to_ENTREZID = TRUE){

  if(to_ENTREZID){
    target_character <- bitr(target_character,
                             fromType = "SYMBOL", toType = "ENTREZID",
                             OrgDb = org.Hs.eg.db)
    target_character <- target_character$ENTREZID
  }

if (is.data.frame(geneset)) {
    res <- enricher(target_character,
                    TERM2GENE = geneset, # a data.frame of 2 column with term and gene
                    minGSSize = arguments$minGSSize,
                    maxGSSize = arguments$maxGSSize,
                    pvalueCutoff = 1,qvalueCutoff = 1)

  }else if(is.list(geneset)){

    geneset <- to_df(geneset)
    res <- enricher(target_character,
                    TERM2GENE = geneset, # a data.frame of 2 column with term and gene
                    minGSSize = arguments$minGSSize,
                    maxGSSize = arguments$maxGSSize,
                    pvalueCutoff = 1,qvalueCutoff = 1)

  }else if(geneset=="kegg"){
    res <- enrichKEGG(
      target_character,
      organism = "hsa",keyType = "kegg",
      pvalueCutoff = 1,
      minGSSize = arguments$minGSSize,
      maxGSSize = arguments$maxGSSize,
      qvalueCutoff = 1,
      use_internal_data = FALSE
    ) %>%
      setReadable(org.Hs.eg.db, keyType="ENTREZID")
  }else if(geneset=="mkegg"){
    res <- enrichMKEGG(
      target_character,
      organism = "hsa",keyType = "kegg",
      pvalueCutoff = 1,minGSSize = 5,
      maxGSSize = 500,qvalueCutoff = 1
    ) %>%
      setReadable(org.Hs.eg.db, keyType="ENTREZID")
  }else if(geneset=="go"){
    res <- enrichGO(
      target_character,
      OrgDb='org.Hs.eg.db',
      keyType = "ENTREZID", ont = "BP",
      pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 10, maxGSSize = 500,
      readable = TRUE
    )
  }else if(geneset=="wp"){
    # WikiPathways
    res <- enrichWP(target_character, organism = "Homo sapiens",
                    minGSSize = arguments$minGSSize,
                    maxGSSize = arguments$maxGSSize,
                    pvalueCutoff = 1,qvalueCutoff = 1)
  }
if (out_dataframe) {
    res <- res@result %>%
       filter(.data$pvalue < arguments$pvalue & .data$qvalue < arguments$qvalue)
  }
  return(res)
}
