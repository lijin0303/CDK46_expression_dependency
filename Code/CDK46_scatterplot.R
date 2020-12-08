require(taigr)
require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)
#===== 1) Data Loading =====
CCLE.expression <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='CCLE_expression')
Achilles.gene.effect <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='Achilles_gene_effect')
proteomics <- load.from.taiga(data.name='total-proteome--5c50',data.version=1, data.file='protein_quant_current_normalized')
CCLE.mutations <- load.from.taiga(data.name='public-20q1-c3b6', data.version=9, data.file='CCLE_mutations')
gene.effect <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=19, data.file='gene_effect')
#===== 2) Data Wrangling =====
mut_annot <- CCLE.mutations%>%
  select(Hugo_Symbol,isCOSMIChotspot,DepMap_ID)%>%
  filter(Hugo_Symbol %in% c("RB1") & isCOSMIChotspot=="True")%>%
  select(Hugo_Symbol,DepMap_ID)%>%distinct()

DepD <- Achilles.gene.effect[,c("CDK4 (1019)","CDK6 (1021)")]%>%
  as.data.frame()%>%
  magrittr::set_colnames(c("CDK4","CDK6"))%>%
  tibble::rownames_to_column("DepMap_ID")%>%
  gather(Dgene,Dependency,-DepMap_ID)%>%
  mutate(source="CERES")%>%
  rbind(
    gene.effect[,c("CDK4 (1019)","CDK6 (1021)")]%>%
      as.data.frame()%>%
      magrittr::set_colnames(c("CDK4","CDK6"))%>%
      tibble::rownames_to_column("DepMap_ID")%>%
      gather(Dgene,Dependency,-DepMap_ID)%>%
      mutate(source="DEMETER"))

CCLE <- CCLE.expression[,c("CDK4 (1019)","CDK6 (1021)")]%>%
  as.data.frame()%>%
  magrittr::set_colnames(c("CDK4","CDK6"))%>%
  tibble::rownames_to_column("DepMap_ID")%>%
  gather(Egene,expr,-DepMap_ID)%>%
  mutate(esource="CCLE")

prote <- proteomics%>%
  arrange(Gene_Symbol)%>%
  filter(Gene_Symbol %in% c("CDK4","CDK6"))%>%
  t()%>%
  as.data.frame(stringsAsFactors = FALSE)%>%
  magrittr::set_colnames(c("CDK4","CDK6"))%>%
  tibble::rownames_to_column("cl")%>%
  filter(grepl("_TenP",cl))%>%
  mutate(cclename = gsub("_TenP.*","",cl))%>%
  mutate(DepMap_ID = purrr::pmap_chr(list(cclename),celllinemapr::ccle.to.arxspan))%>%
  select(DepMap_ID,CDK4,CDK6)%>%
  gather(Egene,expr,-DepMap_ID)%>%
  mutate(esource="proteomics")

expr <- CCLE%>%rbind(prote)
  
combinedD <- expr%>%
  inner_join(DepD,by="DepMap_ID")%>%
  left_join(mut_annot,by="DepMap_ID")%>%
  mutate(Situ = purrr::pmap_chr(list(Egene,Dgene),~paste0("x:",c(...)[2],", ","y:",c(...)[1])))%>%
  mutate_at(vars(expr),as.numeric)%>%
  mutate(Mutation =factor(case_when(is.na(Hugo_Symbol)~"WT",
                                    TRUE~Hugo_Symbol),
                          levels=c("RB1","WT")))

#===== 3) Scatter Plotting =====
sp <- combinedD%>%
  filter(!(esource=="proteomics"&source=="DEMETER"))%>%
  group_by(esource,source,Situ)%>%
  do(plot = if(unique(.$esource)=="CCLE"){
    ggplot(data=.,aes(Dependency,expr,colour = Mutation)) +
      scale_color_manual(values=c("firebrick4","grey10"),name = "Mutations")+
      geom_point(data=.%>%filter(Mutation=="WT"),alpha=0.7,size=1)+
      geom_point(data=.%>%filter(Mutation!="WT"),alpha=0.7,size=1)+
      geom_smooth(method=lm, se=FALSE,color="steelblue")+
      theme(plot.title = element_text(color='black', hjust = 0.5),
            plot.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", size = 0.25, fill = NA),
            text = element_text(color='black',family="serif"),
            panel.grid = element_blank(),
            axis.text = element_text(color='black'),
            strip.background = element_rect(colour="black", fill="gray81", 
                                            size=1, linetype="solid"),
            legend.key = element_rect(fill = NA))+
      guides(colour = guide_legend(override.aes = list(size=5)))+
      geom_vline(xintercept = 0, linetype="dashed")+
      labs(x=paste0(unique(.$Dgene)," ",unique(.$source)),y=paste0(unique(.$Egene)," expression"))
    
    }else{
        ggplot(data=.,aes(Dependency,expr,colour = Mutation)) +
          scale_color_manual(values=c("firebrick4","grey10"),name = "Mutations")+
          geom_point(data=.%>%filter(Mutation=="WT"),alpha=0.7,size=1)+
          geom_point(data=.%>%filter(Mutation!="WT"),alpha=0.7,size=1)+
          geom_smooth(method=lm, se=FALSE,color="steelblue")+
          theme(plot.title = element_text(color='black', hjust = 0.5),
                plot.background = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color = "black", size = 0.25, fill = NA),
                text = element_text(color='black',family="serif"),
                panel.grid = element_blank(),
                axis.text = element_text(color='black'),
                strip.background = element_rect(colour="black", fill="gray81", 
                                                size=1, linetype="solid"),
                legend.key = element_rect(fill = NA))+
          guides(colour = guide_legend(override.aes = list(size=5)))+
          geom_vline(xintercept = 0, linetype="dashed")+
          labs(x=paste0(unique(.$Dgene)," ",unique(.$source)),
               y=paste0(unique(.$Egene)," proteomics"))+ylim(c(-4,4))})

ggarrange(plotlist=sp$plot, nrow = 3,ncol=4,common.legend = T,align="hv")+
  ggsave(width = 12,height = 9,dpi = 2000,
         filename = "~/Desktop/CDK46_scatter.pdf")

combinedD%>%
  filter((esource=="proteomics"&source=="CERES")|(esource=="CCLE"&source=="DEMETER"))%>%
  select(DepMap_ID,
         Expression_gene_symbol = Egene,
         Expression_value  = expr,
         Expression_data_source = esource,
         Dependency_gene_symbol = Dgene,
         Dependency_value = Dependency,
         Dependency_data_source = source,
         Mutation)%>%
  openxlsx::write.xlsx(file = "~/Desktop/Figure_3E.xlsx")

combinedD%>%
  filter((esource=="proteomics"&source=="CERES")|(esource=="CCLE"&source=="DEMETER"))%>%
  select(DepMap_ID,
         Expression_gene_symbol = Egene,
         Expression_value  = expr,
         Expression_data_source = esource,
         Dependency_gene_symbol = Dgene,
         Dependency_value = Dependency,
         Dependency_data_source = source,
         Mutation)%>%
  group_by(Expression_gene_symbol,Expression_data_source,Dependency_gene_symbol,Dependency_data_source)%>%
  summarise(n=n(),p.cor = cor(Expression_value,Dependency_value,method="pearson"),
            pval = summary(lm(Expression_value~Dependency_value))$coefficients[2,4])%>%
  arrange(Expression_data_source,Dependency_gene_symbol)

combinedD%>%
  filter((esource=="CCLE"&source=="CERES"))%>%
  select(DepMap_ID,
         Expression_gene_symbol = Egene,
         Expression_value  = expr,
         Expression_data_source = esource,
         Dependency_gene_symbol = Dgene,
         Dependency_value = Dependency,
         Dependency_data_source = source,
         Mutation)%>%
  openxlsx::write.xlsx(file = "~/Desktop/Figure_S3B.xlsx")

combinedD%>%
  filter((esource=="CCLE"&source=="CERES"))%>%
  select(DepMap_ID,
         Expression_gene_symbol = Egene,
         Expression_value  = expr,
         Expression_data_source = esource,
         Dependency_gene_symbol = Dgene,
         Dependency_value = Dependency,
         Dependency_data_source = source,
         Mutation)%>%
  group_by(Expression_gene_symbol,Expression_data_source,Dependency_gene_symbol,Dependency_data_source)%>%
  summarise(n=n(),p.cor = cor(Expression_value,Dependency_value,method="pearson"),
            pval = summary(lm(Expression_value~Dependency_value))$coefficients[2,4])%>%
  arrange(Expression_data_source,Dependency_gene_symbol)

