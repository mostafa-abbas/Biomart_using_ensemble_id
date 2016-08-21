library("biomaRt")
ensembl37 <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL') #select the ensembl version that contains hg19 (GRch37) 
ensembl <- useDataset("hsapiens_gene_ensembl",mart= ensembl37) #select human genes
attributes <- c("ensembl_gene_id", "ensembl_transcript_id","chromosome_name","uniprot_genename","hgnc_symbol","hgnc_transcript_name","external_gene_name") #select the attributes
filter <- "ensembl_gene_id" #select the filter
gene_id <- read.table("ensembl_gene_id.txt") #upload the values of our filter
gene_id <- as.matrix(gene_id)
geneInfo = data.frame() #creat an empty dataframe
#the followinf while loop to handel the gene_id arry 500 by 500 (500 based on the advise on ensmble website beacuse whene we call the getBM function by the whole gene_id array we will miss some values)
finished = FALSE
limit = 500
k=0
while(!finished)
{
	if((length(gene_id)-(k*limit))>limit)
	{
		gene_id_tmp = gene_id[c((limit*k+1):(limit*(k+1)))]
		k=k+1
	}
	else
	{
		gene_id_tmp = gene_id[c((limit*k+1):length(gene_id))]
		finished = TRUE
	}
	geneInfo <- rbind(geneInfo,getBM(attributes = attributes, filters = filter, values = gene_id_tmp, mart = ensembl))
}
row.na <- apply(geneInfo, 1, function(x){is.na(x[4])}) #select the rows that contains NA in gc content column(also trancript_length)
geneInfo.f <- geneInfo[!row.na,] 
write.csv(geneInfo.f, file = "geneInfo_bimomart_from_ensemble.csv")
