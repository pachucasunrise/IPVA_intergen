################################################################################
## Project: Intergen of IPVA
## Script purpose: Take outputs and put into nice tables. And plot graphs for Tables 3 and 4
## Date: 04.01.24
## Author: Annie Herbert
## Email: annie.herbert@bristol.ac.uk
################################################################################

cd C:/Program Files/R/R-4.2.2/bin
R

#Packages
packages<-c('rlang','installr','tidyverse','readxl','ggplot2','gridExtra','ggpubr','dplyr')
	#'readstata13','data.table','tidyr','formattable','tidyverse','dplyr','gdata','foreign',
	#'matrixStats','tableone','Rcmdr','mice','magrittr','varhandle','zoo','mice','backports')
#source("http://bioconductor.org/biocLite.R")
for(pkg in packages){
  if(!require(pkg,character.only=TRUE)){
    BiocManager::install(pkg,suppressUpdates=TRUE, lib="C:/Users/is19262/R/win-library/4.1")
    library(pkg,character.only=TRUE, lib.loc="C:/Users/is19262/R/win-library/4.1")
    }
  }
rm(pkg,packages)


#Read in data
loc <- "results/"
results <- read.table(paste0(loc,'tidy_tables1to4_R.csv'),sep=",",header=TRUE)
results <- results %>% distinct()

#Empty table
sex <- c(2,1)
dva <- c("dv_mrep","phys_mrep","psych_mrep","cc_mrep")
noyes <- c("no","yes")
num_row <- length(sex)*length(dva)*length(noyes)
outcome <- c("vic_1821","per_1821")

for (o in outcome){
	table <- data.frame(matrix(NA, nrow = num_row, ncol = 10))
	colnames(table) <- c("sex","dva_type","level","perc_w_outcome","crude_rr","95ci","adj_rr_A","95ci_A","adj_rr_B","95ci_B")
	table[,1]<-c(rep(2,length(dva)*length(noyes)),rep(1,length(dva)*length(noyes)))
	table[,2]<-rep(c("dv_mrep","dv_mrep","phys_mrep","phys_mrep","psych_mrep","psych_mrep","cc_mrep","cc_mrep"),2)
	table[,3]<-rep(noyes,length(dva))

	#Loop over and fill in
	table[table$level=="no",c("crude_rr","95ci","adj_rr_A","95ci_A","adj_rr_B","95ci_B")]<-" "
	for (s in sex){
		for (d in dva){
			subset <- results[results$Sex==s & results$Exposure==d & results$Adjusted=="No" & results$Outcome==o,]
			table[table$sex==s & table$dva_type==d & table$level=="no","perc_w_outcome"]<-round(subset$no_perc,digits=1)
			table[table$sex==s & table$dva_type==d & table$level=="yes","perc_w_outcome"]<-round(subset$yes_perc,digits=1)
			table[table$sex==s & table$dva_type==d & table$level=="yes","crude_rr"]<-round(subset$RR,digits=2)
			table[table$sex==s & table$dva_type==d & table$level=="yes","95ci"]<-paste0("(",round(subset$LCI,digits=2)," to ",round(subset$UCI,digits=2),")")
			
			subset <- results[results$Sex==s & results$Exposure==d & results$Adjusted=="Yes" & results$Outcome==o,]			
			table[table$sex==s & table$dva_type==d & table$level=="yes","adj_rr_A"]<-round(subset$RR,digits=2)
			table[table$sex==s & table$dva_type==d & table$level=="yes","95ci_A"]<-paste0("(",round(subset$LCI,digits=2)," to ",round(subset$UCI,digits=2),")")

			subset <- results[results$Sex==s & results$Exposure==d & results$Adjusted=="Yes + ACEs" & results$Outcome==o,]			
			table[table$sex==s & table$dva_type==d & table$level=="yes","adj_rr_B"]<-round(subset$RR,digits=2)
			table[table$sex==s & table$dva_type==d & table$level=="yes","95ci_B"]<-paste0("(",round(subset$LCI,digits=2)," to ",round(subset$UCI,digits=2),")")
			}
		}
	#Save out table
	assign(paste0("neat_table_rrs_",o),table)
	write.csv(table, paste0(loc,"neat_table_rrs_",o,".csv"), row.names = TRUE)
	}


#PAF/Joint PAF tables
results_PAFs <- read.table(paste0(loc,'tidy_tables1to4_PAFs_R.csv'),sep=",",header=TRUE)
results_PAFs <- results_PAFs %>% distinct()

#Temporary
results_PAFs$PAF <- results_PAFs$PAF/100
results_PAFs$LPAF <- results_PAFs$LPAF/100
results_PAFs$UPAF <- results_PAFs$UPAF/100

results_PAFs$RD <- results_PAFs$RD*100
results_PAFs$LRD <- results_PAFs$LRD*100
results_PAFs$URD <- results_PAFs$URD*100
results_PAFs$Risk.ref <- results_PAFs$Risk.ref*100
results_PAFs$pc <- results_PAFs$pc
head(results_PAFs)

exp2 <- c("maltreat","aces3","restoftrio")
exp2_shorter <- c("maltreat","aces3") #As restoftrio causes issues
num_row <- length(dva)*(1+length(exp2))

firstcol <- NA
i<-0
for (d in dva){
	i<-i+1
	firstcol[i] <- d
	for (e in exp2){
		i<-i+1
		firstcol[i] <- paste0(d," + ",e)
		}
	}

#Combinations to loop over
list <- list(outcome = NA, dva = NA, exp = NA)
for (o in c("vic_1821")){
	for (d in dva){
		for (e in exp2_shorter){
			list$outcome <- c(list$outcome,o)
			list$dva <- c(list$dva,d)
			list$exp <- c(list$exp,e)
		}
	}
	for (d in dva){
		for (e in c("restoftrio")){
			list$outcome <- c(list$outcome,o)
			list$dva <- c(list$dva,d)
			list$exp <- c(list$exp,e)
		}
	}
}
for (o in c("per_1821")){
	for (d in dva){
		for (e in exp2_shorter){
			list$outcome <- c(list$outcome,o)
			list$dva <- c(list$dva,d)
			list$exp <- c(list$exp,e)
		}
	}
	for (d in c("dv_mrep","phys_mrep")){
		for (e in c("restoftrio")){
			list$outcome <- c(list$outcome,o)
			list$dva <- c(list$dva,d)
			list$exp <- c(list$exp,e)
		}
	}
}
list$outcome <- list$outcome[2:length(list$outcome)]
list$dva <- list$dva[2:length(list$dva)]
list$exp <- list$exp[2:length(list$exp)]

for (out in c("vic_1821","per_1821")){
	table <- data.frame(matrix(NA, nrow = num_row, ncol = 19))
	colnames(table) <- c("factor","prev_exp_w","paf_w","95ci_w","risk_ref_w","rd_w","rd_95ci_w","prev_exp_m","paf_m","95ci_m","risk_ref_m","rd_m","rd_95ci_m","lci_w","uci_w","lci_m","uci_m","exp1","exp2")
	table[,1]<-firstcol
	#Loop over and fill in
	for (i in 1:length(list$outcome)){
		o <- list$outcome[i]
		d <- list$dva[i]
		e <- list$exp[i]
		if(o == out){		
			#Single PAFs come from 'results'
			#Females
			subset <- results[results$Sex==2 & results[,"Exposure"]==d & results$Adjusted=="Yes" & results$Outcome==o,] #(reduced set)
			table[table$factor==d,"prev_exp_w"]<-round(subset$prev,digits=1)
			table[table$factor==d,"paf_w"]<-round(subset$PAF,digits=2)
			table[table$factor==d,"95ci_w"]<-paste0("(",round(subset$LPAF,digits=2)," to ",round(subset$UPAF,digits=2),")")
			table[table$factor==d,"risk_ref_w"]<-round(subset$Risk.ref,digits=2)
			table[table$factor==d,"rd_w"]<-round(subset$RD,digits=1)
			table[table$factor==d,"rd_95ci_w"]<-paste0("(",round(subset$LRD,digits=1)," to ",round(subset$URD,digits=1),")")

			#Want individual raw ci limits for plot, later	
			table[table$factor==d,"lci_w"]<-subset$LPAF
			table[table$factor==d,"uci_w"]<-subset$UPAF

			#Males
			subset <- results[results$Sex==1 & results[,"Exposure"]==d & results$Adjusted=="Yes" & results$Outcome==o,] #(reduced set)
			table[table$factor==d,"prev_exp_m"]<-round(subset$prev,digits=1)
			table[table$factor==d,"paf_m"]<-round(subset$PAF,digits=1)
			table[table$factor==d,"95ci_m"]<-paste0("(",round(subset$LPAF,digits=1)," to ",round(subset$UPAF,digits=1),")")
			table[table$factor==d,"risk_ref_m"]<-round(subset$Risk.ref,digits=2)
			table[table$factor==d,"rd_m"]<-round(subset$RD,digits=1)
			table[table$factor==d,"rd_95ci_m"]<-paste0("(",round(subset$LRD,digits=1)," to ",round(subset$URD,digits=1),")")

			table[table$factor==d,"lci_m"]<-subset$LPAF
			table[table$factor==d,"uci_m"]<-subset$UPAF

			subset <- results_PAFs[results_PAFs$Sex==2 & results_PAFs[,"Joint.exposure.level"]==paste0(d,"=1 ",e,"=1") & results_PAFs$Adjusted=="Yes" & results_PAFs$Outcome==o,]
			table[table$factor==paste0(d," + ",e),"prev_exp_w"]<-round(subset$prev,digits=1)
			table[table$factor==paste0(d," + ",e),"paf_w"]<-round(subset$PAF,digits=1)
			table[table$factor==paste0(d," + ",e),"95ci_w"]<-paste0("(",round(subset$LPAF,digits=1)," to ",round(subset$UPAF,digits=1),")")
			table[table$factor==paste0(d," + ",e),"risk_ref_w"]<-round(subset$Risk.ref,digits=2)
			table[table$factor==paste0(d," + ",e),"rd_w"]<-round(subset$RD,digits=1)
			table[table$factor==paste0(d," + ",e),"rd_95ci_w"]<-paste0("(",round(subset$LRD,digits=1)," to ",round(subset$URD,digits=1),")")
			table[table$factor==paste0(d," + ",e),"lci_w"]<-subset$LPAF
			table[table$factor==paste0(d," + ",e),"uci_w"]<-subset$UPAF

			subset <- results_PAFs[results_PAFs$Sex==1 & results_PAFs[,"Joint.exposure.level"]==paste0(d,"=1 ",e,"=1") & results_PAFs$Adjusted=="Yes" & results_PAFs$Outcome==o,]
			table[table$factor==paste0(d," + ",e),"prev_exp_m"]<-round(subset$prev,digits=1)
			table[table$factor==paste0(d," + ",e),"paf_m"]<-round(subset$PAF,digits=1)
			table[table$factor==paste0(d," + ",e),"95ci_m"]<-paste0("(",round(subset$LPAF,digits=1)," to ",round(subset$UPAF,digits=1),")")
			table[table$factor==paste0(d," + ",e),"risk_ref_m"]<-round(subset$Risk.ref,digits=2)
			table[table$factor==paste0(d," + ",e),"rd_m"]<-round(subset$RD,digits=1)
			table[table$factor==paste0(d," + ",e),"rd_95ci_m"]<-paste0("(",round(subset$LRD,digits=1)," to ",round(subset$URD,digits=1),")")
			table[table$factor==paste0(d," + ",e),"lci_m"]<-subset$LPAF
			table[table$factor==paste0(d," + ",e),"uci_m"]<-subset$UPAF

			table[table$factor==d,"exp1"] <- d
			table[table$factor==d,"exp2"] <- "any"
			table[table$factor==paste0(d," + ",e),"exp1"] <- d
			table[table$factor==paste0(d," + ",e),"exp2"] <- e
			}
		}
	assign(paste0("neat_table_PAFs_",out),table)
	write.csv(table, paste0(loc,"neat_table_PAFs_",out,".csv"), row.names = TRUE)
  }

# Figures
library(ggplot2)
library(ggpubr)

for (o in c("vic_1821","per_1821")){
	#Get data
	table <- get(paste0("neat_table_PAFs_",o))
	table <- table[is.na(table$prev_exp_w)==FALSE,]

	#Rename variables
	table[table$exp1=="dv_mrep","exp1"] <- "IPVA"
	table[table$exp1=="phys_mrep","exp1"] <- "Phys"
	table[table$exp1=="psych_mrep","exp1"] <- "Psych"
	table[table$exp1=="cc_mrep","exp1"] <- "Control"

	table[table$exp2=="any","exp2"] <- "Any"
	table[table$exp2=="maltreat","exp2"] <- "+ maltreatment"
	table[table$exp2=="restoftrio","exp2"] <- "+ MH + SA"
	table[table$exp2=="aces3","exp2"] <- "+ 3+ other ACEs"

	#Order levels
	table$factor <- factor(table$factor, levels=firstcol)
	table$exp1 <- factor(table$exp1, levels=c("IPVA","Phys","Psych","Control"))
	table$exp2 <- factor(table$exp2, levels=c("Any","+ maltreatment","+ 3+ other ACEs","+ MH + SA"))

		#And prev_exp_w
	ymin <- round(min(c(table$lci_w,table$lci_m#,table$prev_exp_w,table$prev_exp_m
	), na.rm=TRUE),digits=1)-0.1
	ymax <- round(max(c(table$uci_w,table$uci_m#,table$prev_exp_w,table$prev_exp_m
	), na.rm=TRUE),digits=1)+0.1

	ymin_prev <- round(min(c(table$prev_exp_w,table$prev_exp_m), na.rm=TRUE),digits=1)-0.1
	ymax_prev <- round(max(c(table$prev_exp_w,table$prev_exp_m), na.rm=TRUE),digits=1)+0.1

	# plots <- lapply(c("vic_1821","per_1821"), 
	#   function(x) 
	#   assign(paste0("plot_w_",x), 

	assign(paste0("plot_w_",o), 
		ggplot(table,
	  aes(x = exp1, y = paf_w, ymin = lci_w, ymax = uci_w)) + 
	  geom_pointrange(aes(col=exp2), 
	    position=position_dodge(width=0.30)) +
	  scale_color_brewer(palette="Dark2", name="Combination") +
	  geom_hline(yintercept=0, linetype="dashed", color="black") +
	  ylim(ymin,ymax) + 
	  xlab("") +
	  ylab("PAF (95% CI)") +
	  theme(axis.text.x = element_blank()) +
	  ggtitle("Women")
	  )

	assign(paste0("plot_m_",o), 
		ggplot(table,
	  aes(x = exp1, y = paf_m, ymin = lci_m, ymax = uci_m)) + 
	  geom_pointrange(aes(col=exp2), 
	    position=position_dodge(width=0.30)) +
	  scale_color_brewer(palette="Dark2") +
	  geom_hline(yintercept=0, linetype="dashed", color="black") +
	  ylim(ymin,ymax) + 
	  xlab("") +
	  ylab("") +
		theme(axis.text.x = element_blank()) +
	  ggtitle("Men")
		)

	assign(paste0("plot_prev_w_",o), 
		ggplot(table,
	  aes(x = exp1, y = prev_exp_w)) + 
	  geom_point(aes(col=exp2), 
	    position=position_dodge(width=0.30)) +
	  scale_color_brewer(palette="Dark2") +
	  ylim(ymin_prev,ymax_prev) + 
	  xlab("") +
	  ylab("Prevalence")
	  )

	assign(paste0("plot_prev_m_",o), 
		ggplot(table,
	  aes(x = exp1, y = prev_exp_m)) + 
	  geom_point(aes(col=exp2), 
	    position=position_dodge(width=0.30)) +
	  scale_color_brewer(palette="Dark2") +
	  ylim(ymin_prev,ymax_prev) + 
	  xlab("") +
	  ylab("")
	  )

	g_legend<-function(a.gplot){
	  tmp <- ggplot_gtable(ggplot_build(a.gplot))
	  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	  legend <- tmp$grobs[[leg]]
	  return(legend)}

	mylegend<-g_legend(get(paste0("plot_w_",o)))

	#Combine plots
	tiff(paste0(loc,'fig_',o,'.tif'), width = 800, height = 500, unit = "px")
	grid.arrange(arrangeGrob(get(paste0("plot_w_",o)) + theme(legend.position="none"),
	                         get(paste0("plot_m_",o)) + theme(legend.position="none"),
	                         get(paste0("plot_prev_w_",o)) + theme(legend.position="none"),
	                         get(paste0("plot_prev_m_",o)) + theme(legend.position="none"),
	                         nrow=2, ncol=2, widths=c(5,5), heights=c(10,5)),
	                bottom = text_grob("IPVA subtype"),
	                #left = text_grob("PAF (95% CI)",rot=90), 
	                mylegend, 
	                nrow=1, widths=c(8.5,1.5)
	         )
	dev.off()
	}


#Descriptives
for (s in sex){
	results <- read.table(paste0(loc,'tidy_descriptives_',s,'.csv'),sep=",",header=TRUE)
	results <- results %>% distinct()
	
	list <- list(orig = c("COHORT","X0","n","p","X2901","ALSPAC","X","X.1","X.2","X.3"), new = c("cohort_name","cohort_level","cohort_n","cohort_p","cohort_tot","alspac_name","alspac_level","alspac_n","alspac_p","alspac_tot"))

	for (i in 1:length(list$orig)){
		names(results)[i]<-list$new[i]
	  }
	
	results<-results[-1,]
	
	list <- list(vars = c("marstat","smokpreg","matage","parity","mated","hhsocclas",
												"dv_mrep","phys_mrep","psych_mrep","cc_mrep",
												"kz021","ethnicity","birthwgt",
												"vic_1821","vic_phys_1821","vic_emo_1821","vic_emo_co_1821","vic_sex_1821","per_1821","per_phys_1821","per_emo_1821","per_emo_co_1821","per_sex_1821",
											  "maltreat","emotional_ab","emotional_ne","physical_abu","sexual_abu","bullying_0_1","mentl_hlth_p","substance_ho","parent_convi","parental_sep",
											  "aces3","aces4","dv_3plus","aces4_dv_mal","aces4_trio","dv_mal","dv_trio"), 
							class = c("bin","bin","cont","4","bin","5",
												rep("bin",4),
												"bin","bin","cont",
												rep("bin",27)))

	n <- 0
	names <- NA
	levels <- NA
	for (i in 1:length(list$vars)) {
		if (list$class[i]=="bin" | list$class[i]=="cont") {
			n <- n+1 
			names <- c(names,list$vars[i])
			levels <- c(levels,"yes")
			}
		else if (list$class[i]!="bin" & list$class[i]!="cont") {
			n <- n+as.numeric(list$class[i])  
			names <- c(names,rep(list$vars[i],as.numeric(list$class[i])))
			levels <- c(levels,as.numeric(results[results$cohort_name==list$var[i],"cohort_level"]))
			}
		}
	names <- names[2:length(names)]
	levels <- levels[2:length(levels)]

	table <- data.frame(matrix(NA, nrow = n, ncol = 6))
	colnames(table) <- c("var","level","cohort_n","cohort_p","alspac_n","alspac_p")
	table[,1]<-names
	table[,2]<-levels

	for (i in 1:length(list$vars)) {
		if (list$class[i]=="cont") {
			table[table$var==list$var[i],"cohort_n"] <- round(as.numeric(results[results$cohort_name==list$var[i],"cohort_n"]))
			table[table$var==list$var[i],"alspac_n"] <- round(as.numeric(results[results$cohort_name==list$var[i],"alspac_n"]))

			table[table$var==list$var[i],"cohort_p"] <- paste0(#"(",
				round(as.numeric(results[results$cohort_name==list$var[i],"cohort_p"]),digits=1)," to ",round(as.numeric(results[results$cohort_name==list$var[i],"cohort_tot"]),digits=1)#,")"
			)
			table[table$var==list$var[i],"alspac_p"] <- paste0(#"(",
				round(as.numeric(results[results$cohort_name==list$var[i],"alspac_p"]),digits=1)," to ",round(as.numeric(results[results$cohort_name==list$var[i],"alspac_tot"]),digits=1)#,")"
			)
			}
		else if (list$class[i]=="bin") {
			table[table$var==list$var[i],"cohort_n"] <- round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==1,"cohort_n"]))
			table[table$var==list$var[i],"alspac_n"] <- round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==1,"alspac_n"]))
			
			table[table$var==list$var[i],"cohort_p"] <- paste0(#"(",
				round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==1,"cohort_n"])/3243*100,digits=1)#,")"
			)
			table[table$var==list$var[i],"alspac_p"] <- paste0(#"(",
				round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==1,"alspac_n"])/14837*100,digits=1)#,")"
			)
			}
		else {
			for (j in as.numeric(results[results$cohort_name==list$var[i],"cohort_level"])) {
				table[table$var==list$var[i] & table$level==j,"cohort_n"] <- round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==j,"cohort_n"]))
				table[table$var==list$var[i] & table$level==j,"alspac_n"] <- round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==j,"alspac_n"]))
				table[table$var==list$var[i] & table$level==j,"cohort_p"] <- paste0(#"(",
					round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==j,"cohort_n"])/3243*100,digits=1)#,")"
				)
				table[table$var==list$var[i] & table$level==j,"alspac_p"] <- paste0(#"(",
					round(as.numeric(results[results$cohort_name==list$var[i] & results$cohort_level==j,"alspac_n"])/14837*100,digits=1)#,")"
				)
				}
			}
		}

	write.csv(table, paste0(loc,"descriptives_",s,".csv"), row.names = TRUE)
	}
